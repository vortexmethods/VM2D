/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: World2D.cpp                                                      |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Файл кода с описанием класса World2D
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#include "World2D.h"

#include "Airfoil2DRigid.h"
#include "Airfoil2DDeformable.h"

#include "Boundary2DConstLayerAver.h"
#include "Boundary2DLinLayerAver.h"
#include "Boundary2DVortexCollocN.h"

#include "MeasureVP2D.h"

#include "Mechanics2DRigidImmovable.h"
#include "Mechanics2DRigidGivenLaw.h"
#include "Mechanics2DRigidOscillPart.h"
#include "Mechanics2DRigidRotatePart.h"
#include "Mechanics2DDeformable.h"

#include "Passport2D.h"
#include "StreamParser.h"

#include "Velocity2DBiotSavart.h"

#include "Wake2D.h"

#include "intersection.cuh"
#include "GMRES.h"
#include "wrapper.h"


#include <type_traits>

using namespace VM2D;


int multipoleOrderGMRES;


//Конструктор
World2D::World2D(const VMlib::PassportGen& passport_) :
	WorldGen(passport_),
	passport(dynamic_cast<const Passport&>(passport_)),
	cuda(Gpu(*this))
{
	std::stringstream ss;
	ss << "#" << passport.problemNumber << " (" << passport.problemName << ")";		
	info.assignStream(defaults::defaultWorld2DLogStream, ss.str());	
		
	passport.physicalProperties.setCurrTime(passport.timeDiscretizationProperties.timeStart);
	currentStep = 0;
		
	timestat.reset(new Times(*this)),
	
	wake.reset(new Wake(*this));
	// загрузка пелены из файла
	if (passport.wakeDiscretizationProperties.fileWake != "")
		wake->ReadFromFile(passport.wakesDir, passport.wakeDiscretizationProperties.fileWake); //Считываем из каталога с пеленой
	
	source.reset(new WakeDataBase(*this));
	// загрузка положений источников из файла
	if (passport.wakeDiscretizationProperties.fileSource != "")
		source->ReadFromFile(passport.dir, passport.wakeDiscretizationProperties.fileSource); //Считываем из текущего каталога

	//считываем массив точек для подсчета и вывода поля скоростей и давлений
	measureVP.reset(new MeasureVP(*this));
	if (passport.timeDiscretizationProperties.saveVPstep != 0)
		measureVP->ReadPointsFromFile(passport.dir);

	switch (passport.numericalSchemes.velocityComputation.second)
	{
	case 0:
		velocity.reset(new VelocityBiotSavart(*this));
		break;
	//case 1:
	//	velocity.reset(new VelocityBarnesHut(*this));
	//	break;
	}

	velocity->virtualVortexesParams.resize(passport.airfoilParams.size());

	auto CreateBoundary = [this](size_t i) {
		switch (passport.numericalSchemes.boundaryCondition.second)
		{
		case 0:
			this->boundary.emplace_back(new BoundaryVortexCollocN(*this, i));
			//info('e') << "BoundaryMDV is not implemented now! " << std::endl;
			//exit(1);
			break;

		case 1:
			this->boundary.emplace_back(new BoundaryConstLayerAver(*this, i));
			break;

		case 2:
			this->boundary.emplace_back(new BoundaryLinLayerAver(*this, i));
			break;

		default:
			info('e') << "Unknown scheme!" << std::endl;
			exit(1);
		}
	};


	for (size_t i = 0; i < passport.airfoilParams.size(); ++i)
	{
		
		switch (passport.airfoilParams[i].mechanicalSystemType)
		{
		case 0:
			airfoil.emplace_back(new AirfoilRigid(*this, i));
			airfoil[i]->ReadFromFile(passport.airfoilsDir);	
			CreateBoundary(i);
			mechanics.emplace_back(new MechanicsRigidImmovable(*this, i));
			break;

		case 1:
			airfoil.emplace_back(new AirfoilRigid(*this, i));
			airfoil[i]->ReadFromFile(passport.airfoilsDir);	
			CreateBoundary(i);
			mechanics.emplace_back(new MechanicsRigidGivenLaw(*this, i));
			break;

		case 2:
			airfoil.emplace_back(new AirfoilRigid(*this, i));
			airfoil[i]->ReadFromFile(passport.airfoilsDir);	
			CreateBoundary(i);
			mechanics.emplace_back(new MechanicsRigidOscillPart(*this, i));
			break;

		case 3:
			airfoil.emplace_back(new AirfoilRigid(*this, i));
			airfoil[i]->ReadFromFile(passport.airfoilsDir);	
			CreateBoundary(i);
			mechanics.emplace_back(new MechanicsRigidRotatePart(*this, i));
			break;
			
		case 4:
			airfoil.emplace_back(new AirfoilDeformable(*this, i));
			airfoil[i]->ReadFromFile(passport.airfoilsDir);
			CreateBoundary(i);
			mechanics.emplace_back(new MechanicsDeformable(*this, i));
			break;
		}	
		
		
	}

	IQ.resize(passport.airfoilParams.size());
	for (size_t i = 0; i < passport.airfoilParams.size(); ++i)
		IQ[i].resize(passport.airfoilParams.size());


	info.endl();
	info('i') << "Start solving problem " << passport.problemName << std::endl;
	info.endl();

	VMlib::PrintLogoToStream(info('_') << std::endl);
}//World2D(...)



///// Функция выполнения предварительного шага
//void World2D::ZeroStep()
//{
//	getTimestat().ToZero();
//	CalcPanelsVeloAndAttachedSheets();
//
//	getNonConstMeasureVP().Initialization();
//
//	CalcAndSolveLinearSystem();
//}//ZeroStep()



//Основная функция выполнения одного шага по времени
void World2D::Step() // ЮИ
{
	try {
		//Очистка статистики
		getTimestat().ToZero();

#ifdef USE_CUDA
		cuSetCurrentStep((int)currentStep, 1);
#endif // USE_CUDA


		//Засечка времени в начале шага
		getTimestat().timeWholeStep.first += omp_get_wtime();


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		size_t countStrongCoupling = 0;
		for (size_t m = 0; m < mechanics.size(); ++m)
		{
			MechanicsRigidOscillPart* mechVar = dynamic_cast<MechanicsRigidOscillPart*>(mechanics[m].get());
			if (mechVar && getPassport().airfoilParams[m].addedMass.length2() > 0)
			{
				mechVar->getStrongCoupling() = true;
				++countStrongCoupling;
			}
		}

		bool semiImplicitStrategy = (countStrongCoupling == mechanics.size());
		if ((currentStep == 0) && (semiImplicitStrategy))
			info('i') << "Strong (semi-implicit) coupling strategy" << std::endl;
	
//НЕПОДВИЖНЫЕ ТЕЛА
//*
		if (!semiImplicitStrategy)
		{
			CalcPanelsVeloAndAttachedSheets();
			measureVP->Initialization();
			CalcAndSolveLinearSystem();


			//Вычисление скоростей вихрей: для тех, которые в следе, и виртуальных, а также в точках wakeVP
			CalcVortexVelo(false);

			//Расчет и сохранение поля давления
			if (ifDivisible(passport.timeDiscretizationProperties.saveVPstep))
			{
				measureVP->CalcPressure();
				measureVP->SaveVP();
			}

			//Вычисление сил, действующих на профиль и сохранение в файл	
			for (auto& mech : mechanics)
			{
				mech->GetHydroDynamForce();
				mech->GenerateForcesString();
				mech->GeneratePositionString();
			}

			//Движение вихрей (сброс вихрей + передвижение пелены)
			WakeAndAirfoilsMotion(true);
			wake->Restruct();
			wake->SaveKadrVtk();
		}
//*/

		
//СОПРЯЖЕННАЯ ПОСТАНОВКА
		if (semiImplicitStrategy)
		{
			CalcPanelsVeloAndAttachedSheets();
			measureVP->Initialization();

			CalcAndSolveLinearSystem();

			//Вычисление скоростей вихрей: для тех, которые в следе, и виртуальных, а также в точках wakeVP
			CalcVortexVelo(false);

			//for (size_t s = 0; s < mechanics.size(); ++s) {// проверка: вычисление полной силы по циркуляциям новых вихрей
			//	mechanics[s]->GetHydroDynamForce();
			//	mechanics[s]->GenerateForcesString();
			//}

			//Расчет и сохранение поля давления
			if (ifDivisible(passport.timeDiscretizationProperties.saveVPstep))
			{
				measureVP->CalcPressure();
				measureVP->SaveVP();
			}

			std::vector<MechanicsRigidOscillPart*> mechOscilPart(mechanics.size(), nullptr);
			if (mechanics.size() > 0)
			{
				for (size_t s = 0; s < mechanics.size(); ++s)
					mechOscilPart[s] = dynamic_cast<MechanicsRigidOscillPart*>(mechanics[0].get());

				if ((mechOscilPart[0]) && (getPassport().airfoilParams[0].addedMass.length2() > 0))
				{
					WakeAndAirfoilsMotion(false); //профиль - кинематически, здесь Vold := V;  
					wake->Restruct();

					CalcPanelsVeloAndAttachedSheets();
					CalcAndSolveLinearSystem();
				}

			}
			//Вычисление сил, действующих на профиль и сохранение в файл	

			for (size_t m = 0; m < mechanics.size(); ++m)
			{
				mechanics[m]->GetHydroDynamForce();

				MechanicsRigidOscillPart* mechVar = dynamic_cast<MechanicsRigidOscillPart*>(mechanics[m].get());
				if (mechVar)
				{
					mechVar->MoveOnlyVelo();
					Point2D accel = (mechVar->getV() - mechVar->getVOld()) * (1.0 / getPassport().timeDiscretizationProperties.dt);
					mechVar->hydroDynamForce -=
						Point2D{ accel[0] * getPassport().airfoilParams[m].addedMass[0],
								 accel[1] * getPassport().airfoilParams[m].addedMass[1] };
					mechVar->GenerateForcesString();
					mechVar->GeneratePositionString();
				}
				else
					exit(3333);
			}
		}
		
		
		wake->SaveKadrVtk();		
//*/
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//Засечка времени в конце шага
		getTimestat().timeWholeStep.second += omp_get_wtime();

		info('i') << "Step = " << currentStep \
			<< " PhysTime = " << passport.physicalProperties.getCurrTime() \
			<< " StepTime = " << Times::dT(getTimestat().timeWholeStep) << std::endl;
		getTimestat().GenerateStatString();

		passport.physicalProperties.addCurrTime(passport.timeDiscretizationProperties.dt);
		++currentStep;
	}
	catch (...)
	{
		info('e') << "!!! Exception from unknown source !!!" << std::endl;
		exit(-1);
	}
}//Step()


//Проверка проникновения вихрей внутрь профиля
void World2D::CheckInside(std::vector<Point2D>& newPos, const std::vector<std::unique_ptr<Airfoil>>& oldAirfoil)
{
	//getTimestat().timeCheckInside.first += omp_get_wtime();
		

	bool movable = false;
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
		movable = movable || mechanics[afl]->isMoves;

#if (!defined(USE_CUDA))	
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
		wake->Inside(newPos, *airfoil[afl], mechanics[afl]->isMoves, *oldAirfoil[afl]);
#else	
//	//////////////////////// CUDA ///////////////////////
	if (true || movable)
	{
		getTimestat().timeCheckInside.first += omp_get_wtime();
		for (size_t afl = 0; afl < airfoil.size(); ++afl)
			wake->Inside(newPos, *airfoil[afl], mechanics[afl]->isMoves, *oldAirfoil[afl]);
		getTimestat().timeCheckInside.second += omp_get_wtime();
	}
	else
	if (newPos.size() > 0)
	{
		double* devNewpos_ptr;
		cuReserveDevMem((void*&)devNewpos_ptr, newPos.size() * sizeof(double) * 2, 0);
		cuCopyFixedArray(devNewpos_ptr, newPos.data(), newPos.size() * sizeof(double) * 2, 119);

		for (size_t afl = 0; afl < airfoil.size(); ++afl)
		{
			std::vector<double> gamma(airfoil[afl]->getNumberOfPanels(), 0.0);
			
			//Через поиск ближайшей панели и псевдонормали
			getTimestat().timeCheckInside.first += omp_get_wtime();
			//std::vector<int> hit = lbvh_check_inside((int)currentStep, (int)newPos.size(), devNewpos_ptr, (int)airfoil[afl]->getNumberOfPanels(), airfoil[afl]->devRPtr, movable);			
			std::vector<int> hit = lbvh_check_inside_ray((int)currentStep, (int)newPos.size(), this->getWake().devVtxPtr, devNewpos_ptr, (int)airfoil[afl]->getNumberOfPanels(), airfoil[afl]->devRPtr, movable);
			getTimestat().timeCheckInside.second += omp_get_wtime();

			/*for (size_t q = 0; q < hit.size(); ++q)
			{
				if (hit[q] != hit2[q])
					std::cout << "q = " << q << " " << hit[q] << " / " << hit2[q] << std::endl;
			}
			*/
			
			//std::stringstream sss;
			//sss << "hit_" << currentStep;
			//std::ofstream of(getPassport().dir + "dbg/" + sss.str());
			//for (size_t i = 0; i < gamma.size(); ++i)
			//	of << hit[i] << std::endl;
			//of.close();

			for (int i = 0; i < newPos.size(); ++i)
			{
				if (hit[i] != (unsigned int)(-1))
				{
					gamma[hit[i]] += wake->vtx[i].g();
					wake->vtx[i].g() = 0.0;
				}
			}
			airfoil[afl]->gammaThrough = gamma;
		}//for afl

		cuDeleteFromDev(devNewpos_ptr);
	}
#endif

	//getTimestat().timeCheckInside.second += omp_get_wtime();
}




#ifdef USE_CUDA
void World2D::GMRES(
	std::vector<std::vector<double>>& X,
	std::vector<double>& R,
	const std::vector<std::vector<double>>& rhs,
	const std::vector<double> rhsReg,
	int& niter,
	bool linScheme)
{
	size_t nAfl = airfoil.size();

	size_t totalVsize = 0;
	for (size_t i = 0; i < airfoil.size(); ++i)
		totalVsize += airfoil[i]->getNumberOfPanels();

	//std::vector<double> currentSol(totalVsize, 0.0);

	if (linScheme)
		totalVsize *= 2;

	std::vector<double> w(totalVsize + nAfl), g(totalVsize + nAfl + 1);
	size_t m;
	const size_t iterSize = 50;

	std::vector<std::vector<double>> V;
	std::vector<std::vector<double>> H;

	V.reserve(iterSize);
	H.reserve(iterSize + 1);

	V.resize(1);
	H.resize(1);

	double beta;

	std::vector<std::vector <double>> diag(nAfl);

	//#ifndef OLD_OMP
	//#pragma omp simd
	//#endif
	for (int p = 0; p < (int)nAfl; ++p)
	{
		if (!linScheme)
			diag[p].resize(airfoil[p]->getNumberOfPanels());
		else
			diag[p].resize(2 * airfoil[p]->getNumberOfPanels());

		for (int i = 0; i < airfoil[p]->getNumberOfPanels(); ++i)
		{
			diag[p][i] = 0.5 / airfoil[p]->len[i];

			if (linScheme)
				diag[p][i + airfoil[p]->getNumberOfPanels()] = (1.0 / 24.0) / airfoil[p]->len[i];
		}
	}

	//Проверка, что предобуславливатель работает корректно (влияния соседних панелей рассчитаны по прямой схеме)
	/*
	for (int i = 0; i < (int)BH.pointsCopyPan.size(); ++i)
	{
		if ((BH.pointsCopyPan[0][i].a.length2() == 0.0) || (BH.pointsCopyPan[0][i].c.length2() == 0.0))
		{
			std::cout << "eps is too small!" << std::endl;
			exit(-1);
		}
	}
	*/

	std::vector<std::vector<double>> residual(nAfl);
	for (size_t p = 0; p < nAfl; ++p)
		residual[p] = rhs[p];

	for (size_t i = 0; i < nAfl; ++i)
		V[0].insert(V[0].end(), residual[i].begin(), residual[i].end());

	for (size_t i = 0; i < nAfl; ++i)
		V[0].push_back(rhsReg[i]); /// суммарная гамма

	// PRECONDITIONER
	//std::vector<PointsCopy> buf1;
	std::vector<double> vbuf1;
	size_t np = 0;

	for (size_t p = 0; p < nAfl; ++p)
	{
		if (!linScheme)
			vbuf1.resize(airfoil[p]->getNumberOfPanels() + 1);
		else
			vbuf1.resize(2 * airfoil[p]->getNumberOfPanels() + 1);

		if (!linScheme)
			np += ((p == 0) ? 0 : airfoil[p-1]->getNumberOfPanels());
		else
			np += ((p == 0) ? 0 : 2 * airfoil[p-1]->getNumberOfPanels());

		for (size_t i = 0; i < airfoil[p]->getNumberOfPanels(); ++i)
		{
			vbuf1[i] = V[0][np + i];
			if (linScheme)
				vbuf1[i + airfoil[p]->getNumberOfPanels()] = V[0][np + airfoil[p]->getNumberOfPanels() + i];
		}

		if (!linScheme)
			vbuf1[airfoil[p]->getNumberOfPanels()] = V[0][V[0].size() - (nAfl - p)];
		else
			vbuf1[2 * airfoil[p]->getNumberOfPanels()] = V[0][V[0].size() - (nAfl - p)];

		//SolM(vbuf1, vbuf1, /*buf1*/ BH, p, n[p]);

		for (size_t j = np; j < np + airfoil[p]->getNumberOfPanels(); ++j)
		{
			V[0][j] = vbuf1[j - np];
			if (linScheme)
				V[0][j + airfoil[p]->getNumberOfPanels()] = vbuf1[j - np + airfoil[p]->getNumberOfPanels()];
		}

		V[0][V[0].size() - (nAfl - p)] = vbuf1[vbuf1.size() - 1];
	}

	beta = norm(V[0]);
	if (beta > 0)
		V[0] = (1.0 / beta) * V[0];
  	else
	 exit(-200);      

	double gs = beta;
	std::vector<double> c, s;
	c.reserve(iterSize);
	s.reserve(iterSize);

	size_t nTotPan = 0;
	for (size_t s = 0; s < getNumberOfAirfoil(); ++s)
		nTotPan += getAirfoil(s).getNumberOfPanels();
	std::vector<double> bufnewSol((linScheme ? 2 : 1) * nTotPan);

	std::vector<double> bufcurrentSol;
	if (linScheme)
		bufcurrentSol.resize(totalVsize, 0.0);

	this->cuda.AllocateSolution(cuda.dev_sol, nTotPan);
	if (linScheme)	
		this->cuda.AllocateSolution(cuda.dev_solLin, nTotPan);	
	else
		cuda.dev_solLin = nullptr;

	for (int j = 0; j < totalVsize - 1; ++j) //+ n.size()
	{
		double t1 = omp_get_wtime();

		//BH.UpdateGams(V[j]);
		//currentSol = V[j];

		//if (afl.numberInPassport == 0)
		double*& dev_ptr_pt = airfoil[0]->devRPtr;

		double*& dev_ptr_rhs = airfoil[0]->devRhsPtr;
		double*& dev_ptr_rhsLin = airfoil[0]->devRhsLinPtr;

		double timingsToRHS[7];

		double* linPtr = (double*)(!linScheme ? nullptr : dev_ptr_rhsLin);

		double t2 = omp_get_wtime();

		if (!linScheme)
		{
			this->cuda.SetSolution(V[j].data(), cuda.dev_sol, nTotPan);
			
			//std::ofstream bufnewfile(passport.dir + "/dbg/mul" + std::to_string(currentStep) + "_" + std::to_string(j) + "_new.txt");
			//for (auto& q : V[j])
			//	bufnewfile << q << '\n';
			//bufnewfile.close();
			
		}
		
		if (linScheme) 
		{
			size_t npred = 0;
			for (size_t i = 0; i < getNumberOfAirfoil(); ++i)
			{
				for (size_t p = 0; p < getAirfoil(i).getNumberOfPanels(); ++p) 
				{
					bufcurrentSol[npred + p] = V[j][2*npred + p];
					bufcurrentSol[nTotPan + npred + p] = V[j][2*npred + getAirfoil(i).getNumberOfPanels() + p];
				}
				npred += getAirfoil(i).getNumberOfPanels();
			}
			this->cuda.SetSolution(bufcurrentSol.data(), cuda.dev_sol, nTotPan);
			this->cuda.SetSolution(bufcurrentSol.data() + nTotPan, cuda.dev_solLin, nTotPan);

			
			//std::ofstream bufnewfile(passport.dir + "/dbg/mul" + std::to_string(currentStep) + "_" + std::to_string(j) + "_new.txt");
			//for (auto& q : bufcurrentSol)
			//	bufnewfile << q << '\n';
			//bufnewfile.close();
			
		}
		double t2a = omp_get_wtime();



		BHcu::wrapperMatrixToVector wrapper(
			(double*)dev_ptr_pt,    //начала и концы панелей
			(double*)dev_ptr_rhs,   //куда сохранить результат 
			linPtr,
			getNonConstCuda().CUDAptrsAirfoilVrt[0],  //указатели на дерево вихрей
			true,                   //признак перестроения дерева вихрей				
			(int)nTotPan,           //общее число панелей на всех профилях
			timingsToRHS,           //засечки времени				
			1.2,//multipoleTheta,                    //theta
			multipoleOrderGMRES,//multipoleOrder,                    //order
			getPassport().numericalSchemes.boundaryCondition.second
		);

		wrapper.calculate(
			(double*)cuda.dev_sol,
			(double*)cuda.dev_solLin
		);

		//exit(-100);


		double t2b = omp_get_wtime();

		int mul = (linScheme ? 2 : 1);
		
		//w.resize(0);
		//w.resize(totalVsize + nAfl, 0.0);

		if (!linScheme)
			getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_rhs, w.data(), 22);

		if (linScheme)
		{
			getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_rhs, bufnewSol.data(), 22);
			getCuda().CopyMemFromDev<double, 1>(nTotPan, linPtr, bufnewSol.data() + nTotPan, 22);
		}

		double t3 = omp_get_wtime();

		if (linScheme)
		{
			size_t npred = 0;
			for (size_t i = 0; i < getNumberOfAirfoil(); ++i)
			{
				for (size_t j = 0; j < getAirfoil(i).getNumberOfPanels(); ++j)
				{
					w[2 * npred + j] = bufnewSol[npred + j];
					w[2 * npred + getAirfoil(i).getNumberOfPanels() + j] = bufnewSol[nTotPan + npred + j];
				}
				npred += getAirfoil(i).getNumberOfPanels();
			}
		}

		//for (int s = 0; s < totalVsize; ++s)
		//	w[s] = newSol[s];
		size_t cntr = 0;
		for (size_t pi = 0; pi < nAfl; ++pi)
		{
#pragma omp parallel for
				for (int i = 0; i < (int)airfoil[pi]->getNumberOfPanels(); ++i)
				{	
					w[cntr + i] += V[j][totalVsize + pi] - V[j][cntr + i] * diag[pi][i] * airfoil[pi]->len[i];
					if (linScheme)
						w[cntr + i + airfoil[pi]->getNumberOfPanels()] += -V[j][cntr + i + airfoil[pi]->getNumberOfPanels()] * diag[pi][i + airfoil[pi]->getNumberOfPanels()] * airfoil[pi]->len[i];
				}

			cntr += airfoil[pi]->getNumberOfPanels();
			if (linScheme)
				cntr += airfoil[pi]->getNumberOfPanels();
		}	

		//{
		//	std::ofstream bufnewfile(passport.dir + "/dbg/s" + std::to_string(currentStep) + "_" + std::to_string(j) + "_new.txt");
		//	for (auto& q : w)
		//		bufnewfile << q << '\n';
		//	bufnewfile.close();
		//}



		for (size_t i = 0; i < airfoil.size(); ++i)
			w[totalVsize + i] = 0.0;

		cntr = 0;
		for (size_t i = 0; i < airfoil.size(); ++i)
		{
			for (size_t k = 0; k < airfoil[i]->getNumberOfPanels(); ++k)
				w[totalVsize + i] += V[j][cntr + k] * airfoil[i]->len[k];
			cntr += airfoil[i]->getNumberOfPanels();
			if (linScheme)
				cntr += airfoil[i]->getNumberOfPanels();
		}


		// PRECONDITIONER
		std::vector<double> vbuf2;
		size_t np = 0;

		for (size_t p = 0; p < nAfl; ++p)
		{
			if (!linScheme)
				vbuf2.resize(airfoil[p]->getNumberOfPanels() + 1);
			else
				vbuf2.resize(2 * airfoil[p]->getNumberOfPanels() + 1);

			if (!linScheme)
				np += ((p == 0) ? 0 : airfoil[p-1]->getNumberOfPanels());
			else
				np += ((p == 0) ? 0 : 2 * airfoil[p-1]->getNumberOfPanels());

			for (size_t i = 0; i < airfoil[p]->getNumberOfPanels(); ++i)
			{				
				vbuf2[i] = w[np + i];
				if (linScheme)
					vbuf2[i + airfoil[p]->getNumberOfPanels()] = w[np + airfoil[p]->getNumberOfPanels() + i];
			}

			if (!linScheme)
				vbuf2[airfoil[p]->getNumberOfPanels()] = w[w.size() - (nAfl - p)];
			else
				vbuf2[2 * airfoil[p]->getNumberOfPanels()] = w[w.size() - (nAfl - p)];

			//SolM(vbuf2, vbuf2, /*buf2*/ BH, p, n[p]);

			for (size_t j = np; j < np + airfoil[p]->getNumberOfPanels(); ++j)
			{
				w[j] = vbuf2[j - np];
				if (linScheme)
					w[j + airfoil[p]->getNumberOfPanels()] = vbuf2[j - np + airfoil[p]->getNumberOfPanels()];
			}

			w[w.size() - (nAfl - p)] = vbuf2[vbuf2.size() - 1];
		}


		H.resize(j + 2);
		for (int i = 0; i < j + 2; ++i)
			H[i].resize(j + 1);

		for (int i = 0; i <= j; ++i)
		{
			double scal = 0.0;
#pragma omp parallel 
			{
#pragma omp for reduction(+:scal)
				for (int q = 0; q < (int)w.size(); ++q)
					scal += w[q] * V[i][q];

				H[i][j] = scal;
#pragma omp for				
				for (int q = 0; q < (int)w.size(); ++q)
					w[q] -= scal * V[i][q];
			}
		}

		H[j + 1][j] = norm(w);
		V.push_back((1 / H[j + 1][j]) * w);
		m = j + 1;

		if (IterRot(H, rhs[0], gs, c, s, j + 1, (int)totalVsize, 1e-8, (int)m, true))
			break;

		double t4 = omp_get_wtime();

		printf("beforeGPU = %f, copy = %f, GPU = %f, copy = %f, afterGPU = %f\n", t2 - t1, t2a-t2, t2b - t2a, t3-t2b, t4 - t3);

	}

	this->cuda.ReleaseSolution(cuda.dev_sol);
	if (linScheme)
		this->cuda.ReleaseSolution(cuda.dev_solLin);

	//printf("------------------------------------------------------------------------\n");
	niter = (int)m;

	g[0] = beta;
	for (int i = 1; i < m + 1; i++)
		g[i] = 0.;

	//GivensRotations
	double oldValue;
	for (int i = 0; i < m; i++)
	{
		oldValue = g[i];
		g[i] = c[i] * oldValue;
		g[i + 1] = -s[i] * oldValue;
	}
	//end of GivensRotations

	std::vector<double> Y(m);
	double sum;

	// Solve HY=g
	Y[m - 1] = g[m - 1] / H[m - 1][m - 1];

	for (int k = (int)m - 2; k >= 0; --k)
	{
		sum = 0;
		for (int s = k + 1; s < m; ++s)
			sum += H[k][s] * Y[s];
		Y[k] = (g[k] - sum) / H[k][k];
	}
	// end of Solve HY=g

	size_t cntr = 0;
	for (size_t p = 0; p < airfoil.size(); p++) {
		if (!linScheme)
			for (size_t i = 0; i < airfoil[p]->getNumberOfPanels(); i++)
			{
				sum = 0.0;
				for (size_t j = 0; j < m; j++)
					sum += V[j][i + cntr] * Y[j];
				X[p][i] += sum;// / airfoil[p]->len[i];
			}
		else
			for (size_t i = 0; i < 2 * airfoil[p]->getNumberOfPanels(); i++)
			{
				sum = 0.0;
				for (size_t j = 0; j < m; j++)
					sum += V[j][i + cntr] * Y[j];
				X[p][i] += sum;// / airfoil[p]->len[i % airfoil[p]->getNumberOfPanels()];
			}
		cntr += airfoil[p]->getNumberOfPanels();

		if (linScheme)
			cntr += airfoil[p]->getNumberOfPanels();

	}
	sum = 0.0;
	for (size_t p = 0; p < airfoil.size(); p++)
	{
		for (size_t j = 0; j < m; j++)
			sum += V[j][totalVsize + p] * Y[j];
		R[p] += sum;
	}
}

#endif






//Решение системы линейных алгебраических уравнений
void World2D::SolveLinearSystem()
{
	getTimestat().timeSolveLinearSystem.first += omp_get_wtime();

/*
	if (currentStep == 0)
	{
		invMatr = matr;
		std::ifstream fileMatrix("matrix.txt");
		int nx, ny;
		fileMatrix >> nx;
		fileMatrix >> ny;
		for (int i = 0; i < matr.rows(); ++i)
		{
			for (int j = 0; j < matr.cols(); ++j)
			{
				fileMatrix >> matr(i, j);
			}
		}
		fileMatrix.close();
	}
*/

/*
	if (currentStep == 0)
	{
		std::ofstream fileMatrix(passport.dir + "/dbg/matrix.txt");
		//int nx, ny;
		//fileMatrix >> nx;
		//fileMatrix >> ny;
		for (int i = 0; i < matr.rows(); ++i)
		{
			for (int j = 0; j < matr.cols(); ++j)
			{
				fileMatrix << matr(i, j) << " ";
			}
			fileMatrix << std::endl;
		}
		fileMatrix.close();
	}
	*/

	//if (currentStep == 0)
	//{
	//	std::ofstream fileMatrix(passport.dir + "/dbg/matrixReord.txt");
	//	//int nx, ny;
	//	//fileMatrix >> nx;
	//	//fileMatrix >> ny;
	//	fileMatrix.precision(16);
	//	for (int i = 0; i < matrReord.rows(); ++i)
	//	{
	//		for (int j = 0; j < matrReord.cols(); ++j)
	//		{
	//			fileMatrix << matrReord(i, j) << " ";
	//		}
	//		fileMatrix << std::endl;
	//	}
	//	fileMatrix.close();
	//}

	//{
	//	std::ofstream fileMatrix(passport.dir + "/dbg/RhsReord"+std::to_string(currentStep)+".txt");
	//	//int nx, ny;
	//	//fileMatrix >> nx;
	//	//fileMatrix >> ny;
	//	fileMatrix.precision(16);
	//	for (int i = 0; i < rhsReord.size(); ++i)
	//	{			
	//		fileMatrix << rhsReord(i) << std::endl;
	//	}
	//	fileMatrix.close();
	//}

	
//////

	bool GAUSSIAN = true;
	multipoleOrderGMRES = 12;

	if (GAUSSIAN)
	{

		if (useInverseMatrix && (currentStep == 0))
		{
			info('t') << "Inverting matrix... ";

#if (defined(USE_CUDA))	
			invMatr.resize(matrReord.rows(), matrReord.cols());
			for (int i = 0; i < (int)matrReord.rows(); ++i)
				for (int j = 0; j < (int)matrReord.cols(); ++j)
					invMatr(i, j) = (i == j) ? 1.0 : 0.0;
			cuInverseMatrix((int)matrReord.rows(), matrReord.data(), invMatr.data());
#else
			invMatr = matrReord.inverse();
#endif
			info('t') << "done" << std::endl;
		}

		if (currentStep == 0)
			info('t') << "Solving system at step #0... ";
		
		if (useInverseMatrix)
			sol = invMatr * rhsReord;
		else
			sol = matrReord.partialPivLu().solve(rhsReord);
		
		if (currentStep == 0)
			info('t') << "done" << std::endl;
		


		//std::ofstream solFile(getPassport().dir + "/dbg/sol" + std::to_string(currentStep) + "-old.txt");
		//solFile.precision(16);
		//for (int i = 0; i < sol.size(); ++i)
		//	solFile << sol(i) << std::endl;
		//solFile.close();

	}

//////


/*
	if (currentStep == 0)
	{
		invMatr = matr;
		std::ifstream fileMatrix("invMatrix.txt");
		int nx, ny;
		fileMatrix >> nx;
		fileMatrix >> ny;
		for (int i = 0; i < invMatr.rows(); ++i)
		{
			for (int j = 0; j < invMatr.cols(); ++j)
			{
				fileMatrix >> invMatr(i, j);
			}
		}
		fileMatrix.close();
	}
*/

/*
	if (currentStep == 0)
	{
		Eigen::MatrixXd mmul = matr * invMatr;
		for (int i = 0; i < invMatr.rows(); ++i)
			mmul(i,i) -= 1.0;
			
		double maxElem = 0.0;
		for (int i = 0; i < mmul.rows(); ++i)
			for (int j = 0; j < mmul.cols(); ++j)
				if (fabs(mmul(i,j)) > maxElem )
					maxElem = mmul(i,j);

		std::cout << "A*A^-1: MAX ELEMENT = " << maxElem << std::endl;
	}
*/
	
/*
	if (currentStep == 0)
	{
		std::ofstream fileMatrix("invMatrix4x3200.txt");
		fileMatrix << invMatr.rows() << " " << invMatr.cols() << std::endl;
		fileMatrix.precision(18);
		for (int i = 0; i < invMatr.rows(); ++i)
		{
			for (int j = 0; j < invMatr.cols(); ++j)
			{
				fileMatrix << invMatr(i, j) << " ";
			}
			fileMatrix << std::endl;
		}
		fileMatrix.close();
	}
*/


#ifdef USE_CUDA
	if (!GAUSSIAN)
	{
		int nFullVars = (int)getNumberOfBoundary();
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			nFullVars += (int)(boundary[i]->GetUnknownsSize());


		std::vector<double> Grhs(rhsReord.size());
		for (int i = 0; i < rhsReord.size(); ++i)
			Grhs[i] = rhsReord(i);

		std::vector<int> Gpos(getNumberOfBoundary(), 0);
		for (int i = 1; i < getNumberOfBoundary(); ++i)
			Gpos[i] = Gpos[i - 1] + (int)(boundary[i - 1]->GetUnknownsSize());

		std::vector<int> Gvsize(getNumberOfBoundary());
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			Gvsize[i] = (int)(boundary[i]->GetUnknownsSize());

		std::vector<std::vector<double>> Ggam(getNumberOfBoundary());
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			Ggam[i].resize(boundary[i]->GetUnknownsSize());

		std::vector<double> GR(getNumberOfBoundary());


		// for direct GMRES
		//std::vector<double> Gmatr(nFullVars * nFullVars);
		//for (size_t i = 0; i < nFullVars; ++i)
		//	for (size_t j = 0; j < nFullVars; ++j)
		//		Gmatr[i * nFullVars + j] = matrReord(i, j);
		//
		//std::vector<int> Gn(getNumberOfBoundary());
		//for (int i = 0; i < getNumberOfBoundary(); ++i)
		//	Gn[i] = getAirfoil(i).getNumberOfPanels();
		//
		//GMRES_Direct(nFullVars, getNumberOfBoundary(), Gmatr, Grhs, Gn, Gpos, Gvsize, Ggam, GR, currentStep);

		// for fast GMRES
		std::vector<std::vector<double>> GGrhs(getNumberOfBoundary());
		int cntrRhs = 0;
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			GGrhs[i].resize(boundary[i]->GetUnknownsSize());
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			for (int j = 0; j < boundary[i]->GetUnknownsSize(); ++j)
				GGrhs[i][j] = rhsReord(cntrRhs++);

		std::vector<double> GrhsReg(getNumberOfBoundary());
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			GrhsReg[i] = rhsReord[nFullVars - (getNumberOfBoundary() - i)];


		double time_GMRES = -omp_get_wtime();

		int niter;
		bool linScheme = (passport.numericalSchemes.boundaryCondition.second == 2);
		GMRES(Ggam, GR, GGrhs, GrhsReg, niter, linScheme);

		time_GMRES += omp_get_wtime();

		std::cout << "Time_GMRES = " << time_GMRES << std::endl;


		sol.resize(nFullVars);
		int cntr = 0;
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			for (int j = 0; j < boundary[i]->GetUnknownsSize(); ++j)
				sol(cntr++) = Ggam[i][j];
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			sol(cntr++) = GR[i];

		//std::ofstream solFile(getPassport().dir + "/dbg/sol-new" + std::to_string(currentStep) + ".txt");
		//solFile.precision(16);
		//for (int i = 0; i < sol.size(); ++i)
		//	solFile << sol(i) << std::endl;
		//solFile.close();
	}
#endif

	getTimestat().timeSolveLinearSystem.second += omp_get_wtime();
}//SolveLinearSystem()

//Заполнение матрицы системы для всех профилей
void World2D::FillMatrixAndRhs()
{
	getTimestat().timeFillMatrixAndRhs.first += omp_get_wtime();

	Eigen::MatrixXd locMatr;
	Eigen::MatrixXd otherMatr;
	Eigen::VectorXd locLastLine, locLastCol;
	
	std::vector<std::vector<Point2D>> locIQ;
	std::vector<std::vector<Point2D>> locOtherIQ;

	//обнуляем матрицу на первом шаге расчета
	if (currentStep == 0)
	{
		for (int i = 0; i < matrReord.rows(); ++i)
		for (int j = 0; j < matrReord.cols(); ++j)
			matrReord(i, j) = 0.0;
	}

	size_t currentRow = 0;
	size_t currentSkosRow = 0;

	size_t nAllVars = 0;
	for (size_t bou = 0; bou < boundary.size(); ++bou)
		nAllVars += boundary[bou]->GetUnknownsSize();
	
	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
		size_t nVars = boundary[bou]->GetUnknownsSize();
		if (currentStep == 0 || mechanics[bou]->isDeform)
		{
			locMatr.resize(nVars, nVars);
			locLastLine.resize(nVars);
			locLastCol.resize(nVars);
		}

		if (currentStep == 0 || mechanics[bou]->isDeform)
		{
			boundary[bou]->FillMatrixSelf(locMatr, locLastLine, locLastCol);
		}

		//размазываем матрицу
		for (size_t i = 0; i < nVars; ++i)
		{
			if (currentStep == 0 || mechanics[bou]->isDeform)
			{
				for (size_t j = 0; j < nVars; ++j)
				{
					//matr(i + currentRow, j + currentRow) = locMatr(i, j);
					//matrSkos(i + currentSkosRow, j + currentSkosRow) = locMatr(i, j);
					matrReord(i + currentSkosRow, j + currentSkosRow) = locMatr(i, j);
				}
				
				//matr(currentRow + nVars, i + currentRow) = locLastLine(i);
				//matr(i + currentRow, currentRow + nVars) = locLastCol(i);

				matrReord(nAllVars + bou, i + currentSkosRow) = locLastLine(i);
				matrReord(i + currentSkosRow, nAllVars + bou) = locLastCol(i);
			}
		}
				
		if ( (currentStep == 0) || (!useInverseMatrix) )
		{
			size_t currentCol = 0;
			size_t currentSkosCol = 0;
			for (size_t oth = 0; oth < boundary.size(); ++oth)
			{
				size_t nVarsOther = boundary[oth]->GetUnknownsSize();

				if (bou != oth)
				{
					otherMatr.resize(nVars, nVarsOther);

					boundary[bou]->FillMatrixFromOther(*boundary[oth], otherMatr);

					//размазываем матрицу
					for (size_t i = 0; i < nVars; ++i)
					{
						for (size_t j = 0; j < nVarsOther; ++j)
						{
							//matr(i + currentRow, j + currentCol) = otherMatr(i, j);
							//matrSkos(i + currentSkosRow, j + currentSkosCol) = otherMatr(i, j);
							matrReord(i + currentSkosRow, j + currentSkosCol) = otherMatr(i, j);
						}
					}
				}// if (bou != oth)
				currentCol += nVarsOther + 1;
				currentSkosCol += nVarsOther;
			}// for oth
		}// if (currentStep == 0 || mechanics[oth]->isMoves)

		currentRow += nVars + 1;
		currentSkosRow += nVars;
	}// for bou

	velocity->FillRhs(/*rhs, rhsSkos,*/ rhsReord);
	
	getTimestat().timeFillMatrixAndRhs.second += omp_get_wtime();
}//FillMatrixAndRhs()

//вычисляем размер матрицы и резервируем память под нее и под правую часть
void World2D::ReserveMemoryForMatrixAndRhs()
{
	getTimestat().timeReserveMemoryForMatrixAndRhs.first += omp_get_wtime();

	
	if (currentStep == 0)
	{
		dispBoundaryInSystem.resize(boundary.size());
		dispBoundaryInSystem[0] = 0;
		
		for (size_t i = 1; i < boundary.size(); ++i)
		{
			dispBoundaryInSystem[i] = dispBoundaryInSystem[i - 1] + boundary[i - 1]->GetUnknownsSize() + 1;
		}

		size_t matrSize = boundary.size();
		size_t matrSkosSize = 0;

		for (auto it = boundary.begin(); it != boundary.end(); ++it)
		{
			matrSize += (*it)->GetUnknownsSize();
			matrSkosSize += (*it)->GetUnknownsSize();
		}

		//matr.resize(matrSize, matrSize);
		matrReord.resize(matrSize, matrSize);
		//matrSkos.resize(matrSkosSize, matrSkosSize);

		//matr.setZero();
		matrReord.setZero();
		//matrSkos.setZero();


		for (size_t i = 0; i < getNumberOfAirfoil(); ++i)
		{
			size_t nVari = boundary[i]->GetUnknownsSize();
			for (size_t j = 0; j < getNumberOfAirfoil(); ++j)
			{
				size_t nVarj = boundary[j]->GetUnknownsSize();
				IQ[i][j].first.resize(nVari, nVarj);
				IQ[i][j].second.resize(nVari, nVarj);
			}
		}

		//rhs.resize(matrSize);
		rhsReord.resize(matrSize);
		//rhsSkos.resize(matrSkosSize);
	}
	//rhs.setZero();
	rhsReord.setZero();
	//rhsSkos.setZero();

	getTimestat().timeReserveMemoryForMatrixAndRhs.second += omp_get_wtime();
}//ReserveMemoryForMatrixAndRhs()



// Вычисление скоростей (и конвективных, и диффузионных) вихрей (в пелене и виртуальных), а также в точках вычисления VP 
void World2D::CalcVortexVelo(bool shiftTime)
{
	getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
	
	velocity->ResizeAndZero();
	
	//Конвективные скорости всех вихрей (и в пелене, и виртуальных), индуцируемые вихрями в пелене
	//Подготовка CUDA
#if defined(__CUDACC__) || defined(USE_CUDA)
	cuda.RefreshWake(2);	

	cuda.RefreshAfls(2);	
	cuda.RefreshVirtualWakes(2);
#endif
	getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();

	//Вычисление конвективных скоростей вихрей (и в пелене, и виртуальных), а также в точках wakeVP
	if (shiftTime)
		--currentStep;

	velocity->CalcConvVelo();
	
	if (shiftTime)
		++currentStep;

	//Расчет средних значений eps для каждой панели и их передача на видеокарту
	getTimestat().timeCalcVortexDiffVelo.first += omp_get_wtime();
	for (size_t bou = 0; bou < getNumberOfBoundary(); ++bou)
		getNonConstAirfoil(bou).calcMeanEpsOverPanel();

#if defined(__CUDACC__) || defined(USE_CUDA)
	for (size_t i = 0; i < airfoil.size(); ++i)
		cuda.CopyMemToDev<double, 1>(airfoil[i]->getNumberOfPanels(), airfoil[i]->meanEpsOverPanel.data(), airfoil[i]->devMeanEpsOverPanelPtr);
#endif
	getTimestat().timeCalcVortexDiffVelo.second += omp_get_wtime();
	
	//Вычисление диффузионных скоростей вихрей (и в пелене, и виртуальных)
	velocity->CalcDiffVelo();		
	
	getTimestat().timeSaveKadr.first += omp_get_wtime();
	/*
	//Сохранение всех параметров для вихрей в пелене
	if(!(currentStep % 1))
	{
		VMlib::CreateDirectory(passport.dir, "dbg");
		std::ostringstream sss;
		sss << "prmWake";
		sss << currentStep;
		std::ofstream prmtFile(passport.dir + "dbg/" + sss.str());
		prmtFile << "i x y g epsast convVeloX convVeloY diffVeloX diffVeloY I0 I1 I2X I2Y I3X I3Y" << std::endl;
		for (size_t i = 0; i < wake->vtx.size(); ++i)
			prmtFile << i << " " \
			<< wake->vtx[i].r()[0] << " " << wake->vtx[i].r()[1] << " " \
			<< wake->vtx[i].g() << " " \
			<< velocity->wakeVortexesParams.epsastWake[i] << " " \
			<< velocity->wakeVortexesParams.convVelo[i][0] << " " << velocity->wakeVortexesParams.convVelo[i][1] << " "\
			<< velocity->wakeVortexesParams.diffVelo[i][0] << " " << velocity->wakeVortexesParams.diffVelo[i][1] << " "\
			//<< std::endl;
			<< velocity->wakeVortexesParams.I0[i] << " " \
			<< velocity->wakeVortexesParams.I1[i] << " " \
			<< velocity->wakeVortexesParams.I2[i][0] << " " << velocity->wakeVortexesParams.I2[i][1] << " " \
			<< velocity->wakeVortexesParams.I3[i][0] << " " << velocity->wakeVortexesParams.I3[i][1] << " " \
			<< std::endl;

		prmtFile.close();
	}
//*/

/*
	//Сохранение всех параметров для виртуальных вихрей
	if (!(currentStep % 1))
	{
		for (size_t b = 0; b < boundary.size(); ++b)
		{
			std::ostringstream sss;
			sss << "prmVirtual_";
			sss << b << "-";
			sss << currentStep;
			std::ofstream prmFileVirt(passport.dir + "dbg/" + sss.str());
			//prmFileVirt.precision(16);
			prmFileVirt << "i x y g epsast convVeloX convVeloY diffVeloX diffVeloY I0 I1 I2X I2Y I3X I3Y" << std::endl;
			for (size_t i = 0; i < boundary[b]->virtualWake.vtx.size(); ++i)
				prmFileVirt << i << " " \
				<< boundary[b]->virtualWake.vtx[i].r()[0] << " " << boundary[b]->virtualWake.vtx[i].r()[1] << " " \
				<< boundary[b]->virtualWake.vtx[i].g() << " " \
				<< velocity->virtualVortexesParams[b].epsastWake[i] << " " \
				<< velocity->virtualVortexesParams[b].convVelo[i][0] << " " << velocity->virtualVortexesParams[b].convVelo[i][1] << " "\
				<< velocity->virtualVortexesParams[b].diffVelo[i][0] << " " << velocity->virtualVortexesParams[b].diffVelo[i][1] << " "\
				//<< std::endl;				
				<< velocity->virtualVortexesParams[b].I0[i] << " " \
				<< velocity->virtualVortexesParams[b].I1[i] << " " \
				<< velocity->virtualVortexesParams[b].I2[i][0] << " " << velocity->virtualVortexesParams[b].I2[i][1] << " " \
				<< velocity->virtualVortexesParams[b].I3[i][0] << " " << velocity->virtualVortexesParams[b].I3[i][1] << " " \
				<< std::endl;
			prmFileVirt.close();
		}
		if (currentStep==3) exit(-123);
	}
//*/

	
	//printf("conv[140] = {%f, %f}\n", velocity->wakeVortexesParams.convVelo[140][0], velocity->wakeVortexesParams.convVelo[140][1]);
	//printf("diff[140] = {%f, %f}\n", velocity->wakeVortexesParams.diffVelo[140][0], velocity->wakeVortexesParams.diffVelo[140][1]);
	
	getTimestat().timeSaveKadr.second += omp_get_wtime();
}//CalcVortexVelo()



// Вычисление скоростей панелей и интенсивностей присоединенных слоев вихрей и источников
void World2D::CalcPanelsVeloAndAttachedSheets()
{
	getTimestat().timeOther.first += omp_get_wtime();

	//вычисляем скорости панелей
	for (size_t i = 0; i < airfoil.size(); ++i)
		mechanics[i]->VeloOfAirfoilPanels(currentStep * passport.timeDiscretizationProperties.dt);

	//вычисляем интенсивности присоединенных слоев
	for (size_t i = 0; i < airfoil.size(); ++i)
		boundary[i]->ComputeAttachedSheetsIntensity();

	getTimestat().timeOther.second += omp_get_wtime();
}//CalcPanelsVeloAndAttachedSheets(...)


//Вычисляем новые положения вихрей (в пелене и виртуальных)
void World2D::MoveVortexes(std::vector<Point2D>& newPos)
{
	getTimestat().timeMoveVortexes.first += omp_get_wtime();

	size_t nvt = wake->vtx.size();
	size_t nVirtVortex = 0;
	for (size_t i = 0; i < getNumberOfBoundary(); ++i)
		nVirtVortex += boundary[i]->virtualWake.vtx.size();

	//newPos.clear();
	newPos.resize(nvt);

#pragma omp parallel for
	for (int i = 0; i < (int)wake->vtx.size(); ++i)
	{
		newPos[i] = wake->vtx[i].r() \
			+ (velocity->wakeVortexesParams.convVelo[i] +
				velocity->wakeVortexesParams.diffVelo[i] +
				passport.physicalProperties.V0()) * passport.timeDiscretizationProperties.dt;
	}

	size_t counter = wake->vtx.size() - nVirtVortex;
	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
		for (int i = 0; i < (int)boundary[bou]->virtualWake.vtx.size(); ++i)
		{
			wake->vtx[counter].g() = boundary[bou]->virtualWake.vtx[i].g();
			++counter;
		}
	}

	getTimestat().timeMoveVortexes.second += omp_get_wtime();
}//MoveVortexes(...)

//Метод - обертка для вызова метода генерации заголовка файла нагрузок и заголовка файла положения(последнее --- если профиль движется)
void World2D::GenerateMechanicsHeader(size_t mechanicsNumber)
{
	mechanics[mechanicsNumber]->GenerateForcesHeader();
	mechanics[mechanicsNumber]->GeneratePositionHeader();
}//GenerateMechanicsHeader(...)

// Заполнение матрицы, состоящей из интегралов от (r-xi) / |r-xi|^2
void World2D::FillIQ()
{
	getTimestat().timeFillMatrixAndRhs.first += omp_get_wtime();
	
	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
		if (currentStep == 0 || mechanics[bou]->isDeform)
		{
			boundary[bou]->FillIQSelf(IQ[bou][bou]);
		}

		for (size_t oth = 0; oth < boundary.size(); ++oth)
		{			

#ifdef INITIAL
			if (currentStep == 0 || !useInverseMatrix)
			{
				//size_t nVarsOther = boundary[oth]->GetUnknownsSize();
				if (bou != oth)
				{
					//std::cout << "matr!" << std::endl;
					boundary[bou]->FillIQFromOther(*boundary[oth], IQ[bou][oth]);
				}
			}// if (currentStep == 0 || !useInverseMatrix)
#endif

#ifdef BRIDGE
			if (currentStep == 0)
			{
				//size_t nVarsOther = boundary[oth]->GetUnknownsSize();
				if (bou != oth)
					boundary[bou]->FillIQFromOther(*boundary[oth], IQ[bou][oth]);
			}// if (currentStep == 0)
#endif

		}// for oth

	}// for bou
	   
	getTimestat().timeFillMatrixAndRhs.second += omp_get_wtime();
}//FillIQ()

void World2D::CalcAndSolveLinearSystem()
{
	if (airfoil.size() > 0)
	{
		ReserveMemoryForMatrixAndRhs();

#if defined(__CUDACC__) || defined(USE_CUDA)
		cuda.setAccelCoeff(passport.physicalProperties.accelCft());
		cuda.setMaxGamma(passport.wakeDiscretizationProperties.maxGamma);

		int sch = passport.numericalSchemes.boundaryCondition.second;

		if ((sch == 1) || (sch == 2) || (sch == 0))
			cuda.setSchemeSwitcher(sch);
		else
		{
			info('e') << "schemeSwitcher is not 0, or 1, or 2! " << std::endl;
			exit(1);
		}

		cuda.RefreshWake(1);
		cuda.RefreshAfls(1);
		cuda.RefreshVirtualWakes(1);

#endif
		if (currentStep == 0)
		{
#ifdef BRIDGE
			useInverseMatrix = true;
#endif //BRIDGE

#ifdef INITIAL
			useInverseMatrix = (
				((mechanics.size() == 1) && (!mechanics.front()->isDeform))
				||
				(mechanics.size() > 1 && !std::any_of(mechanics.begin(), mechanics.end(), [](const std::unique_ptr<Mechanics>& m) { return m->isMoves; }))
				);
#endif //INITIAL
		}

		FillIQ();
		FillMatrixAndRhs();

		//{
		//	std::stringstream ss;
		//	ss << "IQ1-" << currentStep;
		//	std::ofstream of(passport.dir + "dbg/" + ss.str());
		//	for (size_t i = 0; i < (IQ[0][0]).first.rows(); ++i)
		//	{
		//		for (size_t j = 0; j < (IQ[0][0]).first.cols(); ++j)
		//			of << (IQ[0][0]).first(i, j) << " ";
		//		of << std::endl;
		//	}
		//	of.close();
		//}

		/*
		{
			std::stringstream ss;
			ss << "matr-" << currentStep;
			std::ofstream of(passport.dir + "dbg/" + ss.str());
			of.precision(16);
			for (size_t i = 0; i < matr.rows(); ++i)
			{
				for (size_t j = 0; j < matr.cols(); ++j)
					of << matr(i, j) << " ";
				of << std::endl;
			}
			of.close();
		}

		{
			std::stringstream ss;
			ss << "rhs-" << currentStep;
			std::ofstream of(passport.dir + "dbg/" + ss.str());
			of.precision(16);
			for (size_t i = 0; i < rhs.rows(); ++i)
			{
				of << rhs(i) << std::endl;
			}
			of.close();
		}

		//*/

		SolveLinearSystem();

		getTimestat().timeOther.first += omp_get_wtime();

		size_t currentRow = 0;
		for (size_t bou = 0; bou < boundary.size(); ++bou)
		{
			size_t nVars = boundary[bou]->GetUnknownsSize();
			Eigen::VectorXd locSol;
			locSol.resize(nVars);
			for (size_t i = 0; i < nVars; ++i)
				locSol(i) = sol(currentRow + i);

			boundary[bou]->SolutionToFreeVortexSheetAndVirtualVortex(locSol);
			currentRow += nVars;// +1;
		}

		//14-05-2024
		size_t nVirtVortices = 0;
		for (size_t bou = 0; bou < boundary.size(); ++bou)
			nVirtVortices += boundary[bou]->virtualWake.vtx.size();

		wake->vtx.reserve(wake->vtx.size() + nVirtVortices);
		for (size_t bou = 0; bou < boundary.size(); ++bou)
			for (size_t v = 0; v < boundary[bou]->virtualWake.vtx.size(); ++v)
				wake->vtx.push_back(Vortex2D{ boundary[bou]->virtualWake.vtx[v].r(), 0.0 });

		/*
		if (currentStep == 0)
		{
			Point2D addMass = { 0.0, 0.0 };
			for (size_t q = 0; q < airfoil[0]->getNumberOfPanels(); ++q)
			{
				addMass += (sol(q) + boundary[0]->sheets.attachedVortexSheet(q,0)) * 0.5 * (airfoil[0]->getR(q) + airfoil[0]->getR(q + 1)).kcross() * airfoil[0]->len[q];			
			}
			addMass *= passport.physicalProperties.rho;

			std::cout << "AddMass = " << addMass << std::endl;
			//exit(-42);
		}
		*/

		getTimestat().timeOther.second += omp_get_wtime();
	}
}//CalcAndSolveLinearSystem()

void World2D::WakeAndAirfoilsMotion(bool dynamics)
{
	std::vector<Point2D> newPos;

	MoveVortexes(newPos);

//#ifdef BRIDGE
//	double totalForce = 0;
//	for (size_t afl = 0; afl < airfoil.size(); ++afl)
//	{
//		totalForce += mechanics[afl]->hydroDynamForce[1];
//	}
//	totalForce *= 1.2;
//
//	mechanics[0]->hydroDynamForce[1] = totalForce;
//#endif

	getTimestat().timeOther.first += omp_get_wtime();

	oldAirfoil.resize(0);
	for (auto& afl : airfoil)
	{
		if (dynamic_cast<AirfoilRigid*>(afl.get()))
			oldAirfoil.emplace_back(new AirfoilRigid(*afl));

		if (dynamic_cast<AirfoilDeformable*>(afl.get()))
			oldAirfoil.emplace_back(new AirfoilDeformable(*afl));
		

#ifdef BRIDGE
		if (dynamics)
		{
			if (afl->numberInPassport == 0)
			{
				mechanics[afl->numberInPassport]->Move();
			}
			else
			{
				MechanicsRigidOscillPart* mechTest = dynamic_cast<MechanicsRigidOscillPart*>(mechanics[afl->numberInPassport].get());
				if (mechTest != nullptr)
				{
					Mechanics& mechGen = *mechanics[0];
					MechanicsRigidOscillPart& mech0 = dynamic_cast<MechanicsRigidOscillPart&>(mechGen);

					Point2D dr = mech0.getR() - mech0.getROld();
					Point2D dv = mech0.getV() - mech0.getVOld();

					double dphi = mech0.getPhi() - mech0.getPhiOld();
					double dw = mech0.getW() - mech0.getWOld();

					Mechanics& mechGenI = *mechanics[afl->numberInPassport];
					MechanicsRigidOscillPart& mechI = dynamic_cast<MechanicsRigidOscillPart&>(mechGenI);

					mechI.getROld() = mechI.getR();
					mechI.getVOld() = mechI.getV();
					mechI.getPhiOld() = mechI.getPhi();
					mechI.getWOld() = mechI.getW();

					//std::cout << "afl = " << afl << ", dy = " << dy << std::endl;

					airfoil[afl->numberInPassport]->Move(dr);
					airfoil[afl->numberInPassport]->Rotate(dphi);

					mechI.getR() += dr;
					mechI.getV() += dv;

					mechI.getPhi() += dphi;
					mechI.getW() += dw;
				}
			}
		}
#endif

#ifdef INITIAL
		if (dynamics)
			mechanics[afl->numberInPassport]->Move();
		else
		{
			MechanicsRigidOscillPart* mech = dynamic_cast<MechanicsRigidOscillPart*>(mechanics[afl->numberInPassport].get());
			if (mech)
				mech->MoveKinematic();
			else
				exit(2222);
		}
#endif
	}//for

	for (auto& bou : boundary)
		bou->virtualWake.vtx.clear();

	getTimestat().timeOther.second += omp_get_wtime();

	CheckInside(newPos, oldAirfoil);
	
	getTimestat().timeOther.first += omp_get_wtime();	

	//передача новых положений вихрей в пелену	
	for (size_t i = 0; i < wake->vtx.size(); ++i)
		wake->vtx[i].r() = newPos[i];

//	getWake().SaveKadrVtk();
	getTimestat().timeOther.second += omp_get_wtime();
}//WakeAndAirfoilsMotion()

