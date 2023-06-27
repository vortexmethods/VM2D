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

#include "Airfoil2DRect.h"

#include "Boundary2DConstLayerAver.h"
#include "Boundary2DLinLayerAver.h"
#include "Boundary2DVortexCollocN.h"

#include "MeasureVP2D.h"

#include "Mechanics2DRigidImmovable.h"
#include "Mechanics2DRigidGivenLaw.h"
#include "Mechanics2DRigidOscillPart.h"
#include "Mechanics2DRigidRotatePart.h"

#include "Passport2D.h"
#include "StreamParser.h"

#include "Velocity2DBiotSavart.h"

#include "Wake2D.h"

#include "intersection.cuh"

using namespace VM2D;

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

	for (size_t i = 0; i < passport.airfoilParams.size(); ++i)
	{
		switch (passport.numericalSchemes.panelsType.second)
		{
		case 0:			
			airfoil.emplace_back(new AirfoilRect(*this, i));
			break;
		//case 1:
		//	airfoil.emplace_back(new AirfoilCurv(*this, i));
		//	break;
		}

		airfoil[i]->ReadFromFile(passport.airfoilsDir);	//Считываем из каталога с коллекцией профилей
		
		switch (passport.numericalSchemes.boundaryCondition.second)
		{
		case 0:
			boundary.emplace_back(new BoundaryConstLayerAver(*this, i));
			break;

		case 1:
			boundary.emplace_back(new BoundaryLinLayerAver(*this, i));
			break;
		
		case 10:			
			boundary.emplace_back(new BoundaryVortexCollocN(*this, i));
			//info('e') << "BoundaryMDV is not implemented now! " << std::endl;
			//exit(1);
			break;

		default:
			info('e') << "Unknown scheme!" << std::endl;
			exit(1);
		}


		switch (passport.airfoilParams[i].mechanicalSystemType)
		{
		case 0:
			mechanics.emplace_back(new MechanicsRigidImmovable(*this, i));
			break;
		case 1:
			mechanics.emplace_back(new MechanicsRigidGivenLaw(*this, i));
			break;
		case 2:
			mechanics.emplace_back(new MechanicsRigidOscillPart(*this, i));
			break;

		case 3:
			mechanics.emplace_back(new MechanicsRigidRotatePart(*this, i));
			break;
		/*
		case 4:
			mechanics.emplace_back(new MechanicsRigidOscillMon(*this, i));
			break;
		*/
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

		//Засечка времени в начале шага
		getTimestat().timeWholeStep.first += omp_get_wtime();

		CalcPanelsVeloAndAttachedSheets();		
		measureVP->Initialization();
		
		//for (int i = 0; i < wake->vtx.size(); ++i)
		//{
		//	if (std::isnan(wake->vtx[i].r()[0]) || std::isnan(wake->vtx[i].r()[1]) || std::isnan(wake->vtx[i].g()))
		//	{
		//		std::cout << wake->vtx[i].r() << " " << wake->vtx[i].g() << std::endl;
		//		exit(-4);
		//	}
		//}
		
		
		CalcAndSolveLinearSystem();

		//Вычисление скоростей вихрей: для тех, которые в следе, и виртуальных, а также в точках wakeVP
		CalcVortexVelo(false);

		//Расчет и сохранение поля давления
		if (ifDivisible(passport.timeDiscretizationProperties.saveVPstep))
		{
			measureVP->CalcPressure();
			measureVP->SaveVP();
		}

		std::vector<MechanicsRigidOscillPart*> mechOscilPart;
		if (mechanics.size() > 0)
		{
			mechOscilPart.resize(mechanics.size(), nullptr);

			for (size_t s = 0; s < mechanics.size(); ++s)
				mechOscilPart[s] = dynamic_cast<MechanicsRigidOscillPart*>(mechanics[0].get());

			if (getPassport().airfoilParams[0].addedMass.length2() > 0)
			{
				for (size_t s = 0; s < mechanics.size(); ++s)
					mechanics[s]->hydroDynamForce = { 0.0, 0.0 }; //12-01

				WakeAndAirfoilsMotion();                      //12-01
				CalcPanelsVeloAndAttachedSheets(); //12-01
				CalcAndSolveLinearSystem();        //12-01
			}

		}
		//Вычисление сил, действующих на профиль и сохранение в файл	

		for (auto& mech : mechanics)
		{
			mech->GetHydroDynamForce();			
			mech->GenerateForcesString();
			mech->GeneratePositionString();
		}
		
		for (size_t s = 0; s < mechanics.size(); ++s)
			if (mechOscilPart[s])
			{
				if (getPassport().airfoilParams[0].addedMass.length2() > 0)
					mechOscilPart[s]->UpdateU();     //12-01
			}

		//Движение вихрей (сброс вихрей + передвижение пелены)
		WakeAndAirfoilsMotion();  //12-01	

		wake->Restruct();

		wake->SaveKadrVtk();

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
	getTimestat().timeCheckInside.first += omp_get_wtime();
		

#if (!defined(USE_CUDA))	
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
		wake->Inside(newPos, *airfoil[afl], mechanics[afl]->isMoves, *oldAirfoil[afl]);
#else	
	//////////////////////// CUDA ///////////////////////
	if (newPos.size() > 0)
	{
		double* devNewpos_ptr;
		cuReserveDevMem((void*&)devNewpos_ptr, newPos.size() * sizeof(double) * 2, 0);
		cuCopyFixedArray(devNewpos_ptr, newPos.data(), newPos.size() * sizeof(double) * 2, 0);

		for (size_t afl = 0; afl < airfoil.size(); ++afl)
		{
			std::vector<double> gamma(airfoil[afl]->getNumberOfPanels(), 0.0);
			std::vector<unsigned int> hit = lbvh_check_inside((int)currentStep, (int)newPos.size(), devNewpos_ptr, (int)airfoil[afl]->getNumberOfPanels(), airfoil[afl]->devRPtr);
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

	getTimestat().timeCheckInside.second += omp_get_wtime();
}

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


	if (useInverseMatrix && (currentStep == 0))
	{
		info('t') << "Inverting matrix" << std::endl;

#if (defined(USE_CUDA))	
		invMatr.resize(matr.rows(), matr.cols());
		for (int i = 0; i < (int)matr.rows(); ++i)
			for (int j = 0; j < (int)matr.cols(); ++j)
				invMatr(i, j) = (i == j) ? 1.0 : 0.0;
		cuInverseMatrix((int)matr.rows(), matr.data(), invMatr.data());
#else
		invMatr = matr.inverse();
#endif
	}



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


	if (useInverseMatrix)
		sol = invMatr * rhs;
	else
		sol = matr.partialPivLu().solve(rhs);

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
		for (int i = 0; i < matr.rows(); ++i)
		for (int j = 0; j < matr.cols(); ++j)
			matr(i, j) = 0.0;
	}

	size_t currentRow = 0;

	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
		size_t nVars = boundary[bou]->GetUnknownsSize();
		if (currentStep == 0)
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
					matr(i + currentRow, j + currentRow) = locMatr(i, j);
				matr(currentRow + nVars, i + currentRow) = locLastLine(i);
				matr(i + currentRow, currentRow + nVars) = locLastCol(i);
			}
		}
				
		if ( (currentStep == 0) || (!useInverseMatrix) )
		{
			size_t currentCol = 0;
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
							matr(i + currentRow, j + currentCol) = otherMatr(i, j);
					}
				}// if (bou != oth)
				currentCol += nVarsOther + 1;
			}// for oth
		}// if (currentStep == 0 || mechanics[oth]->isMoves)

		currentRow += nVars + 1;
	}// for bou

	velocity->FillRhs(rhs);
	
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
		for (auto it = boundary.begin(); it != boundary.end(); ++it)
			matrSize += (*it)->GetUnknownsSize();

		matr.resize(matrSize, matrSize);
		matr.setZero();

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

		rhs.resize(matrSize);
	}
	rhs.setZero();

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
	/*//Сохранение всех параметров для вихрей в пелене

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
	{
		for (size_t b = 0; b < boundary.size(); ++b)
		{
			std::ostringstream sss;
			sss << "prmVirtual_";
			sss << b << "-";
			sss << currentStep;
			std::ofstream prmFileVirt(passport.dir + "dbg/" + sss.str());
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
		//if (currentStep==2) exit(-123);
	}
*/
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
	for (size_t i = 0; i < getNumberOfBoundary(); ++i)
		nvt += boundary[i]->virtualWake.vtx.size();

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

	int curCounterNewPos = (int)wake->vtx.size();

	wake->vtx.resize(nvt);

	
	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
#pragma omp parallel for
		for (int i = 0; i < (int)boundary[bou]->virtualWake.vtx.size(); ++i)
		{
			wake->vtx[curCounterNewPos + i] = (boundary[bou]->virtualWake.vtx[i]);
			newPos[curCounterNewPos + i] = (boundary[bou]->virtualWake.vtx[i].r() \
				+ (velocity->virtualVortexesParams[bou].convVelo[i] +
					velocity->virtualVortexesParams[bou].diffVelo[i] +
					passport.physicalProperties.V0()) * passport.timeDiscretizationProperties.dt);
		}
		curCounterNewPos += (int)boundary[bou]->virtualWake.vtx.size();
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
				size_t nVarsOther = boundary[oth]->GetUnknownsSize();
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

		if ((sch == 0) || (sch == 1) || (sch == 10))
			cuda.setSchemeSwitcher(sch + 1);
		else
		{
			info('e') << "schemeSwitcher is not 0, or 1, or 10! " << std::endl;
			exit(1);
		}

		cuda.RefreshWake(1);
		cuda.RefreshAfls(1);
		cuda.RefreshVirtualWakes(1);
#endif

		if (currentStep == 0)
		{
			useInverseMatrix = (
				(mechanics.size() == 1)
				||
				(mechanics.size() > 1 && !std::any_of(mechanics.begin(), mechanics.end(), [](const std::unique_ptr<Mechanics>& m) { return m->isMoves; }))
				);
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
			for (size_t i = 0; i < matr.rows(); ++i)
			{
				for (size_t j = 0; j < matr.cols(); ++j)
					of << matr(i, j) << " ";
				of << std::endl;
			}
			of.close();
		}
		*/

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
			currentRow += nVars + 1;
		}

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

		getTimestat().timeOther.second += omp_get_wtime();
	}
}//CalcAndSolveLinearSystem()

void World2D::WakeAndAirfoilsMotion()
{
	std::vector<Point2D> newPos;

	MoveVortexes(newPos);

#ifdef BRIDGE
	double totalForce = 0;
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
	{
		totalForce += mechanics[afl]->hydroDynamForce[1];
	}
	totalForce *= 1.2;

	mechanics[0]->hydroDynamForce[1] = totalForce;
#endif

	getTimestat().timeOther.first += omp_get_wtime();

	oldAirfoil.resize(0);
	for (auto& afl : airfoil)
	{
		switch (passport.numericalSchemes.panelsType.second)
		{
		case 0:
			oldAirfoil.emplace_back(new AirfoilRect(*afl));
			break;
		//case 1:
		//	oldAirfoil.emplace_back(new AirfoilCurv(*afl));
		//	break;
		}


#ifdef BRIDGE
		if (afl->numberInPassport == 0)
		{
			mechanics[afl->numberInPassport]->Move();
		}
		else
		{
			Mechanics& mechGen = *mechanics[0];
			MechanicsRigidOscillPart& mech0 = dynamic_cast<MechanicsRigidOscillPart&>(mechGen);
			double dy = mech0.getY() - mech0.getYOld();
			double du = mech0.getU() - mech0.getUOld();
			
			Mechanics& mechGenI = *mechanics[afl->numberInPassport];
			MechanicsRigidOscillPart& mechI = dynamic_cast<MechanicsRigidOscillPart&>(mechGenI);
			
			mechI.getUOld() = mechI.getU();
			mechI.getYOld() = mechI.getY();
			
			//std::cout << "afl = " << afl << ", dy = " << dy << std::endl;
			
			airfoil[afl->numberInPassport]->Move({ 0.0, dy });
			mechI.Vcm[1] += du;
			
			mechI.getY() += dy;
			mechI.getU() += du;
		}
#endif

#ifdef INITIAL
		mechanics[afl->numberInPassport]->Move();		
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

