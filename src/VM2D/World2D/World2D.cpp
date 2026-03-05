/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
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
\Version 1.14
\date 6 марта 2026 г.
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

#include "StreamParser.h"

#include "Velocity2DBiotSavart.h"
#include "Velocity2DBarnesHut.h"

#include "Wake2D.h"

#include "Gmres2D.h"
#include "treeKernels.cuh"
#include <type_traits>
#include <numeric>
#include "nummatrix.h"

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
		
	currentTime = passport.timeDiscretizationProperties.timeStart;
	currentStep = 0;

	std::vector<std::string> timerLabels = { "Step", "MatRhs", "Solve", "ConvVel", "Nut", /*"Nbody",*/ "DiffVel", "Force", "VelPres", "Inside", "Restr", "Save"};
	timers = std::make_unique<VMlib::TimersGen>(*this, timerLabels);
		
	wake.reset(new Wake(*this));
	// загрузка пелены из файла
	if (passport.wakeDiscretizationProperties.fileWake != "") 
		wake->ReadFromFile(passport.wakesDir, passport.wakeDiscretizationProperties.fileWake); //Считываем из каталога с пеленой
	
	source.reset(new WakeDataBase(*this));
	// загрузка положений источников из файла
	if (passport.wakeDiscretizationProperties.fileSource != "")
		source->ReadFromFile(passport.dir, passport.wakeDiscretizationProperties.fileSource); //Считываем из текущего каталога


	switch (passport.numericalSchemes.velocityComputation.second)
	{
	case 0:
		velocity.reset(new VelocityBiotSavart(*this));
		break;
	case 1:
		velocity.reset(new VelocityBarnesHut(*this));
		inflTreeWake.reset(new CpuTreeInfo(tree_T::vortex, object_T::point4, scheme_T::noScheme));
		cntrTreeWake.reset(new CpuTreeInfo(tree_T::contr, object_T::point4, scheme_T::noScheme));
		cntrTreeVP.reset(new CpuTreeInfo(tree_T::contr, object_T::point2, scheme_T::noScheme));
		break;
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

	if (getPassport().wakeDiscretizationProperties.sigma0 == 0)
	{
		double sumLength = 0.0;
		size_t totNumPan = 0;
		for (size_t bou = 0; bou < getNumberOfAirfoil(); ++bou)
		{
			for (size_t pnl = 0; pnl < getAirfoil(bou).getNumberOfPanels(); ++pnl)
				sumLength += getAirfoil(bou).len[pnl];
			totNumPan += getAirfoil(bou).getNumberOfPanels();
		}
		double eps = 0.5 * sumLength / totNumPan;
		getNonConstPassport().wakeDiscretizationProperties.sigma0 = eps;
		info('i') << "sigma0 = " << getPassport().wakeDiscretizationProperties.sigma0 << " is calculated automatically" << std::endl;
	}

	if (getPassport().wakeDiscretizationProperties.epscol == 0)
	{
		getNonConstPassport().wakeDiscretizationProperties.epscol = (2.0 / 3.0) * getPassport().wakeDiscretizationProperties.sigma0;

		info('i') << "epscol = " << getPassport().wakeDiscretizationProperties.epscol << " is calculated automatically" << std::endl;
	}

	for (size_t afl = 0; afl < passport.airfoilParams.size(); ++afl)
		if (passport.airfoilParams[afl].chord == 0)
		{
			auto& prm = getNonConstPassport().airfoilParams[afl];
			prm.chord = (prm.initialGab.second[0] - prm.initialGab.first[0]) * prm.scale[0];
			info('i') << "airfoil #" << afl << " chord = " << prm.chord << " is calculated automatically" << std::endl;
		}


#ifdef USE_CUDA
	if (getNumberOfAirfoil() > 0)
		if (passport.numericalSchemes.linearSystemSolver.second == 1 || passport.numericalSchemes.linearSystemSolver.second == 2)
			cuda.Gmres.reset(new GmresSolver(*this));
#endif		

	//2026-03-28
	//считываем массив точек для подсчета и вывода поля скоростей и давлений
	measureVP.reset(new MeasureVP(*this));

	if (passport.timeDiscretizationProperties.saveVPstep != 0)
		measureVP->ReadPointsFromFile(passport.dir);

#ifdef OPTIMIZER
//	std::vector<Point2D> VPPoints, VPHistory;
//	for (size_t a = 0; a < getNumberOfAirfoil(); ++a)
//	{
//		const auto& afl = *airfoil[a];
//		for (size_t p = 0; p < afl.getNumberOfPanels(); ++p)
//			VPHistory.push_back(0.5 * (afl.getR(p) + afl.getR(p + 1)) + afl.nrm[p] * passport.airfoilParams[a].chord * 0.01);
//	}
//	if (passport.timeDiscretizationProperties.saveVPstep != 0)
//		measureVP->SetPoints(VPPoints, VPHistory);
#endif

	IQ.resize(passport.airfoilParams.size());
	for (size_t i = 0; i < passport.airfoilParams.size(); ++i)
		IQ[i].resize(passport.airfoilParams.size());


	info.endl();
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
	//try {
		//Очистка статистики
		getTimers().resetAll();

#ifdef USE_CUDA
		cuSetCurrentStep((int)currentStep, 1);
#endif // USE_CUDA


		//Засечка времени в начале шага
		getTimers().start("Step");

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

		bool semiImplicitStrategy = ((countStrongCoupling == mechanics.size()) && (mechanics.size() > 0));
		if ((currentStep == 0) && (semiImplicitStrategy))
			info('i') << "Strong (semi-implicit) coupling strategy" << std::endl;
	
//НЕПОДВИЖНЫЕ ТЕЛА
//*
		int nTotPan = 0;
		for (size_t s = 0; s < getNumberOfAirfoil(); ++s)
			nTotPan += (int)getAirfoil(s).getNumberOfPanels();


		if (!semiImplicitStrategy)
		{
			CalcPanelsVeloAndAttachedSheets();			

			timerInitialBuild.reset();
			timerInitialBuild.start();
			if (getPassport().numericalSchemes.velocityComputation.second == 0 && getPassport().numericalSchemes.linearSystemSolver.second == 2)
			{
#ifdef USE_CUDA
				cuda.RefreshAfls(3);
				if (nTotPan > 0)
				{
					auto& afl = getAirfoil(0);
					auto& treePnlVrt = *getCuda().inflTreePnlVortex;
					if (getCurrentStep() == 0)
						treePnlVrt.MemoryAllocate((int)getCuda().n_CUDA_pnls);
					
					if (getCurrentStep() == 0 || isAnyMovableOrDeformable())
					{
						//Построение дерева для влияющих панелей (вихревых слоев)
						treePnlVrt.UpdatePanelGeometry(nTotPan, (double4*)afl.devRPtr);
						treePnlVrt.Build();
					}

					treePnlVrt.UpdatePanelAttachedVortexIntensity(afl.devAttachedVortexSheetPtr, afl.devAttachedVortexSheetLinPtr);
					treePnlVrt.UpwardTraversal(multipoleOrder);
				}
#endif
			}

			if (getPassport().numericalSchemes.velocityComputation.second == 1)
			{
#ifdef USE_CUDA
				cuda.RefreshWake(3);
				cuda.RefreshAfls(3);

				//Построение дерева для вихрей
				auto& treeWake = *getCuda().inflTreeWake;
				treeWake.MemoryAllocate((int)getCuda().n_CUDA_wake);
				treeWake.Update((int)getWake().vtx.size(), getWake().devVtxPtr);
				treeWake.Build();
				treeWake.UpwardTraversal(getPassport().numericalSchemes.nbodyMultipoleOrder);

				if (nTotPan > 0)
				{
					auto& afl = getAirfoil(0);
					auto& treePnl = *getCuda().cntrTreePnl;
					auto& treePnlVrt = *getCuda().inflTreePnlVortex;
					auto& treePnlSrc = *getCuda().inflTreePnlSource;
					auto& treePnlAux = *getCuda().auxTreePnl;					

					if (getCurrentStep() == 0)
					{
						treePnl.MemoryAllocate((int)getCuda().n_CUDA_pnls);
						treePnlAux.MemoryAllocate((int)getCuda().n_CUDA_pnls);
						treePnlVrt.MemoryAllocate((int)getCuda().n_CUDA_pnls);
						if (isAnyMovableOrDeformable())
							treePnlSrc.MemoryAllocate((int)getCuda().n_CUDA_pnls);
					}

					if (getCurrentStep() == 0 || isAnyMovableOrDeformable())
					{					
						//Построение дерева для влияющих панелей (вихревых слоев)
						treePnlVrt.UpdatePanelGeometry(nTotPan, (double4*)afl.devRPtr);
						treePnlVrt.Build();

						//Построение контрольного дерева панелей
						treePnl.UpdatePanelGeometry((int)nTotPan, (double4*)afl.devRPtr);
						treePnl.Build();

						//Построение дерева для влияющих панелей (источников)
						if (isAnyMovableOrDeformable())
						{
							treePnlSrc.UpdatePanelGeometry(nTotPan, (double4*)afl.devRPtr);
							treePnlSrc.UpdatePanelAttachedSourceIntensity(afl.devAttachedSourceSheetPtr, afl.devAttachedSourceSheetLinPtr);
							treePnlSrc.Build();
							treePnlSrc.UpwardTraversal(multipoleOrder);
						}
					}

					treePnlVrt.UpdatePanelAttachedVortexIntensity(afl.devAttachedVortexSheetPtr, afl.devAttachedVortexSheetLinPtr);
					treePnlVrt.UpwardTraversal(multipoleOrder);
				}
#else
				//inflTreeWake->Update(getWake().vtx, getPassport().wakeDiscretizationProperties.eps);
				//inflTreeWake->Build();
				//inflTreeWake->UpwardTraversal(getPassport().numericalSchemes.nbodyMultipoleOrder);
#endif
			}
			timerInitialBuild.stop();

			measureVP->Initialization();
			CalcAndSolveLinearSystem();

			//added masses
			if (getPassport().physicalProperties.typeAccel.second == 3)
			{
				std::vector<Point2D> lambdaAdd(getNumberOfAirfoil(), {0.0, 0.0});
				std::vector<double> muAdd(getNumberOfAirfoil());
				for (size_t bou = 0; bou < getNumberOfAirfoil(); ++bou)
				{
					if (dynamic_cast<MechanicsRigidGivenLaw*>(&getNonConstMechanics(bou)))
					{
						for (size_t pnl = 0; pnl < getAirfoil(bou).getNumberOfPanels(); ++pnl)
						{
							double sumGam = getBoundary(bou).sheets.freeVortexSheet(pnl, 0) + getBoundary(bou).sheets.attachedVortexSheet(pnl, 0);
							Point2D rpnl = 0.5 * (getAirfoil(bou).getR(pnl) + getAirfoil(bou).getR(pnl + 1));
							lambdaAdd[bou] += rpnl.kcross() * sumGam * getAirfoil(bou).len[pnl];
							muAdd[bou] += (rpnl - getAirfoil(bou).rcm).length2() * sumGam * getAirfoil(bou).len[pnl];
						}
						lambdaAdd[bou] *= getPassport().physicalProperties.rho;
						muAdd[bou] *= 0.5 * getPassport().physicalProperties.rho;
						info('i') << "Added Masses for airfoil #" << bou << " = { " << lambdaAdd[bou][0] << ", " << lambdaAdd[bou][1] << ", " << muAdd[bou] << " }" << std::endl;

						char direction;
						switch ((int)(getPassport().physicalProperties.timeAccel))
						{
						case 0:
							direction = 'x';
							break;
						case 1:
							direction = 'y';
							break;
						case 2:
							direction = 'w';
							break;
						default:
							direction = '?';
							info('e') << "Wrong dirfection is specified!" << std::endl;
							exit(1);
						}

						std::string addMassFileName = getPassport().dir + "addMass-" + std::to_string(getAirfoil(bou).numberInPassport);
						std::ofstream addMassFile;
						if (!VMlib::fileExistTest(addMassFileName, info) || getPassport().physicalProperties.timeAccel == 0)
						{
							addMassFile.open(addMassFileName);
							addMassFile << "Added Masses for airfoil:" << std::endl;
						}
						else
							addMassFile.open(addMassFileName, std::ios_base::app);

						addMassFile << direction << "-direction: " << lambdaAdd[bou][0] << " " << lambdaAdd[bou][1] << " " << muAdd[bou] << std::endl;

						addMassFile.close();
					}
				}
				currentTime = getPassport().timeDiscretizationProperties.timeStop;
			}

			timerInitialBuild.start();
			if (nTotPan > 0 && getPassport().numericalSchemes.velocityComputation.second == 1)
			{
				auto& afl = getAirfoil(0);
#ifdef USE_CUDA
				getCuda().inflTreePnlVortex->UpdatePanelFreeAndAttachedVortexIntensity(afl.devFreeVortexSheetPtr, afl.devFreeVortexSheetLinPtr, afl.devAttachedVortexSheetPtr, afl.devAttachedVortexSheetLinPtr);
				getCuda().inflTreePnlVortex->UpwardTraversal(multipoleOrder);
#endif
			}
			timerInitialBuild.stop();

			//Вычисление скоростей вихрей: для тех, которые в следе, и виртуальных, а также в точках wakeVP
			CalcVortexVelo();


#ifdef TURB
			getTimers().start("Nut");
			std::vector<double> nut;
			CalcVeloDifference(nut);		
			getTimers().stop("Nut");

			for (size_t q = 0; q < getWake().vtx.size(); ++q)			
				getNonConstWake().vtx[q].sigma() = 4.48364 * sqrt(passport.timeDiscretizationProperties.dt * (passport.physicalProperties.nu + nut[q]));
#endif

			 
//#include "gammaCirc.h"						

			//Расчет и сохранение поля давления
			if (ifDivisible(passport.timeDiscretizationProperties.saveVPstep))
			
#ifdef OPTIMIZER			
				if ((getCurrentStep() >= OPTIMIZER_START_STEP && getCurrentStep() <= OPTIMIZER_STOP_STEP) || getCurrentStep()==0)  //2026-03-28
#endif
			{
				const double& vRef = getPassport().physicalProperties.vRef;
				//scaleV = 1.0 / vRef;
				double scaleP = 1.0 / (0.5 * getPassport().physicalProperties.rho * sqr(vRef));

#ifdef USE_CUDA			
				measureVP->GPUCalcPressure();
#else
				measureVP->CalcPressure();
#endif
				
				{
					//std::cout << "presForces are calculated\n";
					//MechanicsDeformable* ptr = dynamic_cast<MechanicsDeformable*>(mechanics[0].get());
					Mechanics* ptr = mechanics[0].get();
					
					if (measureVP->getTotalNumberOfRealPoints() > 0)
						measureVP->SaveVP();

					//if (ptr && !ptr->beam->fsi)

						//сохранение главного вектора силы, вычисленного как интеграл от давления
						//*
						if (measureVP->elasticPoints.size() > 0)
						{
							std::ofstream presForcesFile;
							if (currentStep == 0)
							{
								presForcesFile.open(getPassport().dir + "presForcesFile.csv");
								presForcesFile << "time,Fx,Fy,Py" << std::endl;
							}
							else
								presForcesFile.open(getPassport().dir + "presForcesFile.csv", std::ios_base::app);

							auto r = measureVP->GetVPinElasticPoints();
							Point2D presForce = { 0.0, 0.0 };
							double yPower = 0.0;

							for (int q = 0; q < r.size(); ++q)
							{
								presForce -= r[q].second * getAirfoil(0).len[q] * getAirfoil(0).nrm[q];
								yPower -= (r[q].second * getAirfoil(0).len[q] * getAirfoil(0).nrm[q])[1] * (0.5 * (getAirfoil(0).getV(q) + getAirfoil(0).getV(q + 1))[1]);
							}

							presForcesFile << currentTime << "," << scaleP * presForce[0] << "," << scaleP * presForce[1] << "," << scaleP * yPower << std::endl;
							presForcesFile.close();
						}
				}
				//*/

				//Коэффициенты разложения силы по балочным функциям				
				MechanicsDeformable* ptr = dynamic_cast<MechanicsDeformable*>(mechanics[0].get());
				if (ptr && ptr->beam->fsi)
				{
					auto r = measureVP->GetVPinElasticPoints();

					/*
					std::ofstream elastPresFile;
					elastPresFile.open(getPassport().dir + "elastPresFile-" + std::to_string(currentStep) + ".txt");
					//if (currentStep == 0)
					//	elastPresFile.open(getPassport().dir + "elastPresFile.txt");
					//else
					//	elastPresFile.open(getPassport().dir + "elastPresFile.txt", std::ios_base::app);					
					//elastPresFile << currentStep << " ";				
					
					for (size_t q = 0; q < r.size(); ++q)
						elastPresFile << measureVP->elasticPoints[q][0] << " " << measureVP->elasticPoints[q][1] << " " << r[q].second << std::endl;
					//elastPresFile << std::endl;
					elastPresFile.close();
					*/

					//Сохраняем давление на последних nLastSteps шагах по дисциплине очереди
					std::vector<double> currentPres(ptr->chord.size());
					for (size_t j = 0; j < ptr->chord.size(); ++j)
					{
						//currentPres[j] = -(ptr->beam->rho * 2 * ptr->beam->F); // gravity force for g = 2.0;
						currentPres[j] = -(r[2 * j + 0].second - r[2 * j + 1].second);
					}

					if (ptr->beam->presLastSteps.size() < ptr->beam->nLastSteps)
						ptr->beam->presLastSteps.push_back(currentPres);
					else
					{
						for (int w = 1; w < ptr->beam->nLastSteps; ++w)
							ptr->beam->presLastSteps[w - 1] = std::move(ptr->beam->presLastSteps[w]);
						ptr->beam->presLastSteps.back() = currentPres;
					}				

					for (int q = 0; q < ptr->beam->R; ++q)
					{
						ptr->beam->qCoeff[q] = 0;

						if (ptr->beam->presLastSteps.size() == ptr->beam->nLastSteps)
						{
							for (size_t j = 0; j < ptr->chord.size(); ++j)
							{
								double averpres = 0.0;

								//осредняем давление
								for (int i = 0; i < ptr->beam->nLastSteps; ++i)
									averpres += ptr->beam->presLastSteps[i][j];
								averpres /= ptr->beam->nLastSteps;

								ptr->beam->qCoeff[q] += -averpres * (ptr->initialChord[j].beg - ptr->initialChord[j].end).length() * ptr->beam->shape(q, 0.5 * (ptr->initialChord[j].beg + ptr->initialChord[j].end)[0]);
							}
							ptr->beam->qCoeff[q] /= (ptr->beam->intSqUnitShape * ptr->beam->L);
						}

						//std::cout << "q[" << q << "] = " << ptr->beam->qCoeff[q] << std::endl;
					}

					std::ofstream phiFile;
					if (currentStep == 0)
					{
						phiFile.open(getPassport().dir + "phiFile.csv");
						phiFile << "time";
						for (int p = 0; p < ptr->beam->R; ++p)
							phiFile << ",phi-" << std::to_string(p + 1);
						for (int p = 0; p < ptr->beam->R; ++p)
							phiFile << ",q-" << std::to_string(p + 1);
						phiFile << std::endl;
					}
					else
						phiFile.open(getPassport().dir + "phiFile.csv", std::ios_base::app);					
					
					phiFile << currentTime;
					for (int p = 0; p < ptr->beam->R; ++p)
						phiFile << "," << ptr->beam->phi(p, 0);
					for (int p = 0; p < ptr->beam->R; ++p)
						phiFile << "," << ptr->beam->qCoeff[p];
					
					phiFile << std::endl;
										
					phiFile.close();
				}
//*/
				//getInfo('e') << "Deformable airfoils are not supported now!";
			}

			//Вычисление сил, действующих на профиль и сохранение в файл	
			for (auto& mech : mechanics)
			{
				mech->GetHydroDynamForce();
				mech->GenerateForcesString();
				mech->GeneratePositionString();
			}

			//Движение вихрей (сброс вихрей + передвижение пелены)
#ifdef TURB
			WakeAndAirfoilsMotion(true, &nut);
			

#else
			WakeAndAirfoilsMotion(true);
#endif
			wake->Restruct();
			wake->SaveKadrVtk();

			//std::string fname = getPassport().dir + "/airfoil-" + std::to_string(getAirfoil(0).numberInPassport) + "-" + std::to_string(currentStep) + ".pfl";
			////if (!fileExistTest(fname, W.getInfo()))
			//{
			//	std::ofstream of(fname);
			//	for (size_t i = 0; i < getAirfoil(0).getNumberOfPanels(); ++i)
			//		of << getAirfoil(0).getR(i)[0] << " " << getAirfoil(0).getR(i)[1] << std::endl;
			//	of.close();
			//}

		}//if (!semiImplicitStrategy)

		
//СОПРЯЖЕННАЯ ПОСТАНОВКА
		if (semiImplicitStrategy)
		{
			CalcPanelsVeloAndAttachedSheets();
			measureVP->Initialization();

			CalcAndSolveLinearSystem();

			//Вычисление скоростей вихрей: для тех, которые в следе, и виртуальных, а также в точках wakeVP
			CalcVortexVelo();

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

			WakeAndAirfoilsMotion(false); //профиль - кинематически, здесь Vold := V;  
			wake->Restruct();

			CalcPanelsVeloAndAttachedSheets();
			CalcAndSolveLinearSystem();

			//Вычисление сил, действующих на профиль и сохранение в файл	

			for (size_t m = 0; m < mechanics.size(); ++m)
			{
				mechanics[m]->GetHydroDynamForce();

				MechanicsRigidOscillPart* mechVar = dynamic_cast<MechanicsRigidOscillPart*>(mechanics[m].get());
				if (mechVar)
				{
					if (getPassport().airfoilParams[m].addedMass.length2() == 0)
					{
						getInfo('e') << "Added mass of the airfoil should be non-zero!" << std::endl;
						exit(1);
					}
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
		getTimers().stop("Step");
		/*
		std::cout << "nvt = " << nVtxBeforeMerging << std::endl;
		std::cout << "tIBT = " << timerInitialBuild.duration() << std::endl;
		std::cout << "tRHS = " << timerRhs.duration() << std::endl;
		std::cout << "tFIL = " << timerFillMatrix.duration() << std::endl;
		std::cout << "tLIN = " << timerSlaeSolve.duration() << std::endl;
		std::cout << "tVEL = " << timerConvVelo.duration() << std::endl;
		std::cout << "tINS = " << timerInside.duration() << std::endl;
		std::cout << "tMRG = " << timerMerging.duration() << std::endl;
		*/

		info('i') << "Step = " << getCurrentStep() \
			<< " PhysTime = " << getCurrentTime() \
			<< std::setprecision(3) \
			<< " StepTime = " << getTimers().durationStep() \
			<< std::setprecision(6) \
			<< std::endl;
		
		getTimers().GenerateStatString(getCurrentStep(), getCurrentTime(), getWake().vtx.size());

		currentTime += passport.timeDiscretizationProperties.dt;
		currentStep += 1;
	//}
	//catch (...)
	//{
	//	info('e') << "!!! Exception from unknown source !!!" << std::endl;
	//	exit(-1);
	//}
}//Step()


//Проверка проникновения вихрей внутрь профиля
void World2D::CheckInside(std::vector<Point2D>& newPos, const std::vector<std::unique_ptr<AirfoilGeometry>>& oldAirfoil)
{
	getTimers().start("Inside");

	gabb = 0;
	check01 = 0;
	check02 = 0;
	checkPan = 0;

#if (!defined(USE_CUDA))	
	nVtxBeforeMerging = newPos.size();
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
		wake->Inside(newPos, *airfoil[afl], mechanics[afl]->isMoves, *oldAirfoil[afl]);


	//Тренируемся работать с древом панелей
	/*
	int nTotPan = 0;	
	for (size_t s = 0; s < getNumberOfAirfoil(); ++s)
		nTotPan += (int)getAirfoil(s).getNumberOfPanels();

	auxTreePnl.reset(new CpuTreeInfo(tree_T::aux, object_T::panel, scheme_T::noScheme));

	std::vector<std::pair<Point2D, Point2D>> pnls(nTotPan);
	int counter = 0;
	for (size_t s = 0; s < getNumberOfAirfoil(); ++s)
		for (size_t p = 0; p < getAirfoil(s).getNumberOfPanels(); ++p)
			pnls[counter++] = { getAirfoil(s).getR(p), getAirfoil(s).getR(p + 1) };

	auxTreePnl->UpdatePanelGeometry(pnls, 0);
	auxTreePnl->Build();
	auxTreePnl->UpwardTraversal(0);

	std::unique_ptr<CpuTreeInfo> cntrTreePnt;
	cntrTreePnt.reset(new CpuTreeInfo(tree_T::contr, object_T::point2, scheme_T::noScheme));

	std::vector<Vortex2D> newPosVtx(newPos.size());
	for (size_t i = 0; i < newPos.size(); ++i)
		newPosVtx[i].r() = newPos[i];

	cntrTreePnt->Update(newPosVtx);
	float tBuild = cntrTreePnt->Build();
	float tUpward = cntrTreePnt->UpwardTraversal(0);

	std::vector<std::pair<int,double>> closeIndexDist(cntrTreePnt->object.size());

	
	VMlib::vmTimer timerA;
	timerA.start();
	double tDownward = auxTreePnl->DownwardTraversalClosestPanelToPoints(*cntrTreePnt, closeIndexDist, false, nullptr);
	timerA.stop();

	std::ofstream treeTimeFile;
	if (getCurrentStep() == 0)
	{
		treeTimeFile.open(passport.dir + "/dbg/treeTime.csv");
		treeTimeFile << "step,time,N,tBUI,tUPW,tDNW\n";
	}
	else
		treeTimeFile.open(passport.dir + "/dbg/treeTime.csv", std::ios::app);

	treeTimeFile << getCurrentStep() << ',' << getCurrentTime() << ',' << newPos.size() << ',' << tBuild << ',' << tUpward << ',' << tDownward << '\n';

	treeTimeFile.close();

	//std::cout << "timeFastNeib = " << timerA.duration() << ", N = " << newPos.size() << "\n";
	/*
	std::ofstream closeFile(passport.dir + "closeFile.txt");
	for (size_t i = 0; i < newPos.size(); ++i)
		closeFile << newPos[i][0] << " " << newPos[i][1] << " " << closeIndexDist[i].first << "\n";
	closeFile.close();
	*/

	//*/ 

	//std::cout << "gab: " << gabb << " check1: " << check01 << " check2: " << check02 << " checkPan: " << checkPan << std::endl;
#else	
//	//////////////////////// CUDA ///////////////////////
	timerInside.reset();
	timerInside.start();
	
	if ((newPos.size() > 0) && (getNumberOfAirfoil() > 0))
	{
		int nTotPanels = 0;
		for (size_t afl = 0; afl < airfoil.size(); ++afl)
			nTotPanels += (int)airfoil[afl]->getNumberOfPanels();

		nVtxBeforeMerging = newPos.size();

		auto& auxTree = *getNonConstCuda().auxTreePnl;

		if ((currentStep == 0) || isAnyMovableOrDeformable())
		{
			auxTree.MemoryAllocate((int)getCuda().n_CUDA_pnls);
			auxTree.UpdatePanelGeometry(nTotPanels, (double4*)airfoil[0]->devRPtr);
			auxTree.Build();
			auxTree.UpwardTraversal(0);
		}

		std::vector<int> hit(newPos.size());
		
		if (isAnyMovableOrDeformable()) //Если профили подвижны - метод псевдонормалей
		{
			double* devNewpos_ptr;
			cudaMalloc(&devNewpos_ptr, newPos.size() * sizeof(double) * 2);
			cudaMemcpy(devNewpos_ptr, newPos.data(), newPos.size() * sizeof(double) * 2, cudaMemcpyHostToDevice);

			auto& cntrTreePnt = *getCuda().cntrTreePoint;
			cntrTreePnt.MemoryAllocate((int)getCuda().n_CUDA_wake);
			cntrTreePnt.Update((int)newPos.size(), devNewpos_ptr);
			cntrTreePnt.Build();

			BHcu::treeClosestPanelToPointsCalculationWrapper(auxTree, cntrTreePnt, getWake().devNearestPanelPtr, true, getAirfoil(0).devPsnPtr);

			cudaMemcpy(hit.data(), getWake().devNearestPanelPtr, newPos.size() * sizeof(int), cudaMemcpyDeviceToHost);
			cudaFree(devNewpos_ptr);


			//std::ofstream vortexFile(getPassport().dir + "vortexFile.txt");
			//for (size_t i = 0; i < getWake().vtx.size(); ++i)
			//	vortexFile << getWake().vtx[i].r()[0] << " " << getWake().vtx[i].r()[1] << " " << hit[i] << '\n';
			//vortexFile.close();

			//std::ofstream panelFile(getPassport().dir + "panelFile.txt");
			//for (size_t i = 0; i < getAirfoil(0).getNumberOfPanels(); ++i)
			//	panelFile << getAirfoil(0).getR(i)[0] << " " << getAirfoil(0).getR(i)[1] << " " << \
			//		getAirfoil(0).getR(i + 1)[0] << " " << getAirfoil(0).getR(i + 1)[1] << " " << '\n';
			//panelFile.close();


		}
		else //иначе (для неподвижного профиля - метод трассировки лучей
		{
			std::vector<std::pair<Point2D, Point2D>> segments(newPos.size());
			for (size_t i = 0; i < newPos.size(); ++i)
			{
				segments[i].first = wake->vtx[i].r();
				segments[i].second = newPos[i];
			}
			double* devSegments_ptr;

			cudaMalloc(&devSegments_ptr, newPos.size() * sizeof(double) * 4);
			cudaMemcpy(devSegments_ptr, segments.data(), newPos.size() * sizeof(double) * 4, cudaMemcpyHostToDevice);

			auto& cntrTreeSeg = *getCuda().cntrTreeSegment;
			cntrTreeSeg.MemoryAllocate((int)getCuda().n_CUDA_wake);
			cntrTreeSeg.UpdatePanelGeometry((int)newPos.size(), (double4*)devSegments_ptr);
			cntrTreeSeg.Build();

			BHcu::treePanelsSegmentsIntersectionCalculationWrapper(auxTree, cntrTreeSeg, getWake().devNearestPanelPtr);
			cudaMemcpy(hit.data(), getWake().devNearestPanelPtr, newPos.size() * sizeof(int), cudaMemcpyDeviceToHost);
			cudaFree(devSegments_ptr);
		}

		size_t panCounter = 0;
		for (size_t afl = 0; afl < airfoil.size(); ++afl)
		{
			std::vector<double> gamma(airfoil[afl]->getNumberOfPanels(), 0.0);
			for (int i = 0; i < newPos.size(); ++i)
			{
				if (hit[i] != -1)
				{
					if ((hit[i] >= panCounter) && (hit[i] < panCounter + airfoil[afl]->getNumberOfPanels()))
					{
						if (fabs(wake->vtx[i].g()) > 1.0)
							std::cout << "Too large gamma is through: i = " << i << ", " << hit[i] << ", " << "gamma[hit[i] - panCounter] += " << wake->vtx[i].g() << std::endl;

						gamma[hit[i] - panCounter] += wake->vtx[i].g();
						wake->vtx[i].g() = 0.0;
					}
				}
			}
			airfoil[afl]->gammaThrough = gamma;
			panCounter += airfoil[afl]->getNumberOfPanels();
		}//for afl
	}
	timerInside.stop();

#endif

	getTimers().stop("Inside");
}





//Решение системы линейных алгебраических уравнений
void World2D::SolveLinearSystem()
{
	getTimers().start("Solve");
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

//*
	//if (currentStep == 0)
	//{
	//	std::ofstream fileMatrix(passport.dir + "/dbg/matrix.txt");
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
	//	std::ofstream fileMatrix(passport.dir + "/dbg/rhs"+std::to_string(currentStep)+".txt");
	//	//int nx, ny;
	//	//fileMatrix >> nx;
	//	//fileMatrix >> ny;
	//	fileMatrix.precision(16);
	//	for (int i = 0; i < rhsReord.size(); ++i)
	//	{			
	//		fileMatrix << rhsReord(i) << std::endl;
	//	}
	//	fileMatrix.close();
	//}//*/

	//exit(-100500);
	
//////
	const int linSystemScheme = passport.numericalSchemes.linearSystemSolver.second;

	if (linSystemScheme == 0) // Gaussian elimination
	{
		double t1 = -omp_get_wtime();
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
		{
			/*
			std::cout << "Solution via invMatrix" << std::endl;
			for (int i = 0; i < invMatr.rows(); ++i)
				for (int j = 0; j < invMatr.cols(); ++j)
					if (fabs(invMatr(i, j)) > 1e+5)
						std::cout << "invMatrixError: " << "invA(" << i << ", " << j << ") = " << invMatr(i, j) << std::endl;

			for (int i = 0; i < rhsReord.size(); ++i)
				if (fabs(rhsReord(i)) > 1e+5)
					std::cout << "rhsReordError: " << "rhsReord(" << i << ") = " << rhsReord(i) << std::endl;
			*/
			timerSlaeSolve.reset();
			timerSlaeSolve.start();
			sol = invMatr * rhsReord;
			timerSlaeSolve.stop();
		}
		else
		{
			timerSlaeSolve.reset();
			timerSlaeSolve.start();
			sol = matrReord.partialPivLu().solve(rhsReord);
			timerSlaeSolve.stop();
		}
		
		if (currentStep == 0)
			info('t') << "done" << std::endl;
		
		t1 += omp_get_wtime();

		//std::cout << " Time in Gauss = " << t1 << std::endl;
/*
		std::ofstream solFile(getPassport().dir + "/sol" + std::to_string(currentStep) + "-Gauss.txt");
		solFile.precision(16);
		for (int i = 0; i < sol.size(); ++i)
			solFile << sol(i) << std::endl;
		solFile.close();	
		exit(10101);
//*/
	}//Gauss

/*
	if (currentStep == 0)
	{  N 
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

	

	if ((linSystemScheme == 1 || linSystemScheme == 2)) // GMRES
	{
		int nFullVars = (int)getNumberOfBoundary();
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			nFullVars += (int)(boundary[i]->GetUnknownsSize());

		std::vector<std::vector<double>> Ggam(getNumberOfBoundary());
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			Ggam[i].resize(boundary[i]->GetUnknownsSize());

		std::vector<double> GR(getNumberOfBoundary());

#ifndef USE_CUDA
		std::vector<double> Grhs(rhsReord.size());
		for (int i = 0; i < rhsReord.size(); ++i)
			Grhs[i] = rhsReord(i);

		std::vector<int> Gpos(getNumberOfBoundary(), 0);
		for (int i = 1; i < getNumberOfBoundary(); ++i)
			Gpos[i] = Gpos[i - 1] + (int)(boundary[i - 1]->GetUnknownsSize());

		std::vector<int> Gvsize(getNumberOfBoundary());
		for (int i = 0; i < getNumberOfBoundary(); ++i)
			Gvsize[i] = (int)(boundary[i]->GetUnknownsSize());

		if (passport.numericalSchemes.linearSystemSolver.second == 1)
		{
			//for direct GMRES
			std::vector<double> Gmatr(nFullVars * nFullVars);
			for (size_t i = 0; i < nFullVars; ++i)
				for (size_t j = 0; j < nFullVars; ++j)
					Gmatr[i * nFullVars + j] = matrReord(i, j);

			for (int i = 0; i < getNumberOfBoundary(); ++i)
				for (int j = 0; j < boundary[i]->GetUnknownsSize(); ++j)
					auto y_test = airfoil[i]->len[j % airfoil[i]->getNumberOfPanels()];

			//Gmres->GMRES_Direct(*this, nFullVars, (int)getNumberOfBoundary(), Gmatr, Grhs, Gpos, Gvsize, Ggam, GR);
			std::cout << "GMRES_Direct should be modified!" << std::endl;
			exit(-100);


			sol.resize(nFullVars);
			int cntr = 0;
			for (int i = 0; i < getNumberOfBoundary(); ++i)
				for (int j = 0; j < boundary[i]->GetUnknownsSize(); ++j)
					sol(cntr++) = Ggam[i][j] /*/ airfoil[i]->len[j % airfoil[i]->getNumberOfPanels()]*/;
			for (int i = 0; i < getNumberOfBoundary(); ++i)
				sol(cntr++) = GR[i];

			//std::ofstream solFile(getPassport().dir + "/dbg/sol-direct-gmres" + std::to_string(currentStep) + ".txt");
			//solFile.precision(16);
			//for (int i = 0; i < sol.size(); ++i)
			//	solFile << sol(i) << std::endl;
			//solFile.close();
		}
		if (passport.numericalSchemes.linearSystemSolver.second == 2)
		{
			getInfo('e') << "Fast GMRES without CUDA is not implemented" << std::endl;
			exit(-2);
		}
#else		
		if (passport.numericalSchemes.linearSystemSolver.second == 1)
		{
			getInfo('e') << "Direct GMRES for CUDA is not implemented" << std::endl;
			exit(-2);
		}

		if ( (passport.numericalSchemes.linearSystemSolver.second == 2))
		{
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


			int niter;
			bool linScheme = (passport.numericalSchemes.boundaryCondition.second == 2);

			double time_GMRES = -omp_get_wtime();
			getNonConstCuda().Gmres->GMRES(Ggam, GR, GGrhs, GrhsReg, niter);//, linScheme);
			time_GMRES += omp_get_wtime();
			//std::cout << "Time_GMRES = " << time_GMRES << std::endl;
 
			sol.resize(nFullVars);
			int cntr = 0;
			for (int i = 0; i < getNumberOfBoundary(); ++i)
				for (int j = 0; j < boundary[i]->GetUnknownsSize(); ++j)
					sol(cntr++) = Ggam[i][j] / airfoil[i]->len[j % airfoil[i]->getNumberOfPanels()];
			for (int i = 0; i < getNumberOfBoundary(); ++i)
				sol(cntr++) = GR[i];
			
			/*
			std::ofstream solFile(getPassport().dir + "/sol-fast-new-gmres" + std::to_string(currentStep) + ".txt");
			solFile.precision(16);
			for (int i = 0; i < sol.size(); ++i)
				solFile << sol(i) << std::endl;
			solFile.close();
			exit(1010110);
			//*/
		}		
#endif
	}

	getTimers().stop("Solve");	
}//SolveLinearSystem()

//Заполнение матрицы системы для всех профилей
void World2D::FillMatrixAndRhs()
{
	getTimers().start("MatRhs");	

	const int linSystemScheme = passport.numericalSchemes.linearSystemSolver.second;

	if (linSystemScheme == 0 || linSystemScheme == 1)
	{
		timerFillMatrix.reset();
		timerFillMatrix.start();

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
				boundary[bou]->FillMatrixSelf(locMatr, locLastLine, locLastCol);
			

			//размазываем матрицу
			for (size_t i = 0; i < nVars; ++i)
			{
				if (currentStep == 0 || mechanics[bou]->isDeform)
				{
					for (size_t j = 0; j < nVars; ++j)
						matrReord(i + currentSkosRow, j + currentSkosRow) = locMatr(i, j);	

					matrReord(nAllVars + bou, i + currentSkosRow) = locLastLine(i);
					matrReord(i + currentSkosRow, nAllVars + bou) = locLastCol(i);
				}
			}

			if ((currentStep == 0) || (!useInverseMatrix))
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
								matrReord(i + currentSkosRow, j + currentSkosCol) = otherMatr(i, j);							
						}
					}// if (bou != oth)
					currentCol += nVarsOther + 1;
					currentSkosCol += nVarsOther;
				}// for oth
			}// if (currentStep == 0 || mechanics[oth]->isMoves)

			currentRow += nVars + 1;
			currentSkosRow += nVars;
		}// for bou
		timerFillMatrix.stop();
	}

	velocity->FillRhs(rhsReord);
	
	getTimers().stop("MatRhs");
}//FillMatrixAndRhs()

//вычисляем размер матрицы и резервируем память под нее и под правую часть
void World2D::ReserveMemoryForMatrixAndRhs()
{
	//getTimers().start("Mem");
		
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

		if ((getPassport().numericalSchemes.linearSystemSolver.second == 0 || getPassport().numericalSchemes.linearSystemSolver.second == 1 || (getPassport().numericalSchemes.linearSystemSolver.second == 2 && isAnyMovableOrDeformable())))
			{
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
		}

		//rhs.resize(matrSize);
		rhsReord.resize(matrSize);
		//rhsSkos.resize(matrSkosSize);
	}
	//rhs.setZero();
	rhsReord.setZero();
	//rhsSkos.setZero();

	//getTimers().stop("Mem");	
}//ReserveMemoryForMatrixAndRhs()



// Вычисление скоростей (и конвективных, и диффузионных) вихрей (в пелене и виртуальных), а также в точках вычисления VP 
void World2D::CalcVortexVelo()
{
	velocity->ResizeAndZero();
	
	getTimers().start("ConvVel");

	//Конвективные скорости всех вихрей (и в пелене, и виртуальных), индуцируемые вихрями в пелене
	//Подготовка CUDA
#if defined(__CUDACC__) || defined(USE_CUDA)
	cuda.RefreshWake(2);	

	cuda.RefreshAfls(2);	
	cuda.RefreshVirtualWakes(2);
#endif
	getTimers().stop("ConvVel");

	velocity->CalcConvVelo();

	//Расчет средних значений eps для каждой панели и их передача на видеокарту
	
	getTimers().start("DiffVel");
	
	for (size_t bou = 0; bou < getNumberOfBoundary(); ++bou)
		getNonConstAirfoil(bou).calcMeanEpsOverPanel();

#if defined(__CUDACC__) || defined(USE_CUDA)
	for (size_t i = 0; i < airfoil.size(); ++i)
		cuda.CopyMemToDev<double, 1>(airfoil[i]->getNumberOfPanels(), airfoil[i]->meanEpsOverPanel.data(), airfoil[i]->devMeanEpsOverPanelPtr);
#endif
	
	getTimers().stop("DiffVel");
	
	//Вычисление диффузионных скоростей вихрей (и в пелене, и виртуальных)
	velocity->CalcDiffVelo();	
 
 	//Обнуление вязких напряжений, если они не были вычислены
	for (auto& afl : airfoil)
	{
		if (afl->viscousStress.size() == 0)
			afl->viscousStress.resize(afl->getNumberOfPanels(), 0.0);
	}
	
	getTimers().start("Save");
	
	/*
	//Сохранение всех параметров для вихрей в пелене
	if (!(currentStep % 1))
	{
		VMlib::CreateDirectory(passport.dir, "dbg");
		std::ostringstream sss, sss2;
		sss << "prmWake";
		sss << currentStep;
		sss2 << "testWake";
		sss2 << currentStep;

		std::ofstream prmtFile(passport.dir + "dbg/" + sss.str());
		prmtFile.precision(11);
		prmtFile << "i x y g sigma epsast convVeloX convVeloY diffVeloX diffVeloY I0 I1 I2X I2Y I3X I3Y" << std::endl;
		for (size_t i = 0; i < wake->vtx.size(); ++i)
			prmtFile << i << " " \
			<< wake->vtx[i].r()[0] << " " << wake->vtx[i].r()[1] << " " \
			<< wake->vtx[i].g() << " " \
			<< wake->vtx[i].sigma() << " " \
			<< velocity->wakeVortexesParams.epsastWake[i] << " " \
			<< velocity->wakeVortexesParams.convVelo[i][0] << " " << velocity->wakeVortexesParams.convVelo[i][1] << " "\
			<< velocity->wakeVortexesParams.diffVelo[i][0] << " " << velocity->wakeVortexesParams.diffVelo[i][1] << " "\
			//<< std::endl;
			<< velocity->wakeVortexesParams.I0[i] << " " \
			<< velocity->wakeVortexesParams.I1[i] << " " \
			<< velocity->wakeVortexesParams.I2[i][0] << " " << velocity->wakeVortexesParams.I2[i][1] << " " \
			<< velocity->wakeVortexesParams.I3[i][0] << " " << velocity->wakeVortexesParams.I3[i][1] << " " \
			<< "\n";

		prmtFile.close();

		std::ofstream prmtFile2(passport.dir + "dbg/" + sss2.str());
		prmtFile2 << wake->vtx.size() << "\n";
		prmtFile2.precision(16);
		for (size_t i = 0; i < wake->vtx.size(); ++i)
			prmtFile2 << wake->vtx[i].r()[0] << " " << wake->vtx[i].r()[1] << " " << wake->vtx[i].g() << "\n";	
		prmtFile2.close();
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
			prmFileVirt << "i x y g sigma epsast convVeloX convVeloY diffVeloX diffVeloY I0 I1 I2X I2Y I3X I3Y" << std::endl;
			for (size_t i = 0; i < boundary[b]->virtualWake.vtx.size(); ++i)
				prmFileVirt << i << " " \
				<< boundary[b]->virtualWake.vtx[i].r()[0] << " " << boundary[b]->virtualWake.vtx[i].r()[1] << " " \
				<< boundary[b]->virtualWake.vtx[i].g() << " " \
				<< boundary[b]->virtualWake.vtx[i].sigma() << " " \
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
		//if (currentStep==3) exit(-123);
	}
//*/

	getTimers().stop("Save");
	
}//CalcVortexVelo()


#ifdef TURB
void World2D::CalcVeloDifference(std::vector<double>& nut) const
{
	nut.resize(getWake().vtx.size());

	//std::vector<Point2D> basePoints = { {0.,0.}, {0.7960681614061415,0.28974511519885826},{0.4235790213720913,0.7336603860367653},{-0.14710745031841074,0.8342878085211184},{-0.6489607110877307,0.5445426933222601},{-0.8471580427441825,0.},{-0.6489607110877307, -0.5445426933222601},{-0.14710745031841083, -0.8342878085211185},{0.42357902137209114, -0.7336603860367653},{0.7960681614061413, -0.28974511519885826},{1.5783500941642636,0.3354886691271114},{1.3054390240228053,0.9484569686703798},{0.8068056870865843,1.3974284418694807},{0.16866831725352804,1.6047718422305033},{-0.4986333369362209,1.534635612175369},{-1.0797167572280424,1.1991469430482575},{-1.474107341276333,0.6563148735601236},{-1.6136113741731684, 0.0},{-1.4741073412763332, -0.6563148735601233},{-1.0797167572280428, -1.1991469430482573},{-0.49863333693622125, -1.5346356121753688},{0.1686683172535277, -1.6047718422305033},{0.8068056870865841, -1.397428441869481},{1.3054390240228053, -0.94845696867038},{1.5783500941642634, -0.3354886691271116},{2.374139932721701,0.3578440178114266},{2.1631871809578898,1.041736042079128},{1.7600257551753844,1.6330652474247764},{1.2004783172997366,2.0792894389479355},{0.5342631127196844,2.340759639009162},{-0.1794237140005648,2.394243072889942},{-0.8771679563438638,2.234987502781785},{-1.496971976377837,1.8771434849703577},{-1.9837634669573248,1.3525070308108138},{-2.2942888678950686,0.7076943915843854},{-2.4009566345994733,0.},{-2.2942888678950686, -0.7076943915843852},{-1.983763466957325, -1.3525070308108136},{-1.4969719763778373, -1.8771434849703577},{-0.877167956343864, -2.2349875027817845},{-0.1794237140005649, -2.394243072889942},{0.5342631127196842, -2.340759639009162},{1.2004783172997364, -2.0792894389479355},{1.7600257551753842, -1.6330652474247767},{2.1631871809578898, -1.041736042079128},{2.3741399327217017, -0.3578440178114266},{3.174096170885181,0.35763474709931864},{3.0149337350056395,1.054970941725745},{2.7045899323860105,1.6994064881590272},{2.258626697172506,2.2586266971725055},{1.6994064881590272,2.7045899323860105},{1.054970941725745,3.0149337350056395},{0.35763474709931853,3.174096170885181},{-0.35763474709931864,3.174096170885181},{-1.054970941725745,3.0149337350056395},{-1.6994064881590272,2.7045899323860105},{-2.2586266971725055,2.258626697172506},{-2.7045899323860105,1.6994064881590272},{-3.0149337350056395,1.054970941725745},{-3.174096170885181,0.35763474709931853},{-3.174096170885181, -0.35763474709931864},{-3.0149337350056395, -1.054970941725745},{-2.7045899323860105, -1.6994064881590272},{-2.258626697172506, -2.2586266971725055},{-1.6994064881590272, -2.7045899323860105},{-1.054970941725745, -3.0149337350056395},{-0.35763474709931853, -3.174096170885181},{0.35763474709931864, -3.174096170885181},{1.054970941725745, -3.0149337350056395},{1.6994064881590272, -2.7045899323860105},{2.2586266971725055, -2.258626697172506},{2.7045899323860105, -1.6994064881590272},{3.0149337350056395, -1.054970941725745},{3.174096170885181, -0.35763474709931853},{3.9723533560766993,0.3680927461358311},{3.8370796268482272,1.0917432814164258},{3.5711387529972023,1.7782158086501256},{3.1835870320217006,2.404133329255034},{2.687622072962485,2.948180972772273},{2.100133367655915,3.39183184921509},{1.4411271393461673,3.719977960086802},{0.7330450547474776,3.921444683149553},{0.0, 3.9893713107821425},{-0.7330450547474779,3.921444683149553},{-1.4411271393461675,3.719977960086802},{-2.1001333676559155,3.39183184921509},{-2.687622072962485,2.948180972772273},{-3.1835870320217,2.4041333292550338},{-3.5711387529972023,1.7782158086501254},{-3.8370796268482272,1.0917432814164256},{-3.9723533560766993,0.368092746135831},{-3.9723533560766993, -0.3680927461358311},{-3.8370796268482272, -1.0917432814164258},{-3.5711387529972023, -1.7782158086501256},{-3.1835870320217006, -2.404133329255034},{-2.687622072962485, -2.948180972772273},{-2.100133367655915, -3.39183184921509},{-1.4411271393461673, -3.719977960086802},{-0.7330450547474776, -3.921444683149553},{0.0, -3.9893713107821425},{0.7330450547474779, -3.921444683149553},{1.4411271393461675, -3.719977960086802},{2.1001333676559155, -3.39183184921509},{2.687622072962485, -2.948180972772273},{3.1835870320217, -2.4041333292550338},{3.5711387529972023, -1.7782158086501254},{3.8370796268482272, -1.0917432814164256},{3.9723533560766993, -0.368092746135831} };

	//std::vector<Point2D> basePoints = { {0.,0.}, { -2.47356, -2.00346 }, { -1.79985,1.13964 }, { -2.95238,-0.054986 }, { -1.70333,-0.910983 }, { -2.71457,1.15701 }, { -2.01605,-3.05827 }, { -1.31647,1.96372 }, { -2.99071,1.97152 }, { -0.0334834,0.767391 }, { -3.2381,0.623172 }, { -1.9359,2.42903 }, { -3.52381,-0.678158 }, { 1.05002,-2.98592 }, { 2.01605,-3.05827 }, { 1.90533,-1.62888 }, { 0.666259,0.130419 }, { 3.52381,-0.678157 }, { 1.33293,0.258719 }, { -0.678158,3.52381 }, { 1.61226,1.64718 }, { -0.699742,0.318031 }, { 0.623172,3.2381 }, { -0.376369,-0.567049 }, { 1.13573,2.72822 }, { 3.2381,0.623172 }, { 2.47356,-2.00346 }, { 2.99071,1.97152 }, { 0.351358,-1.84626 }, { -3.13329,-1.62097 }, { -1.09647,-2.98592 }, { -1.36641,0.446331 }, { -0.054986,2.95238 }, { -0.518096,-3.44188 }, { -1.90533,-1.62888 }, { 3.13329,-1.62097 }, { 0.92973,-2.31605 }, { 1.95024,3.00436 }, { 2.71457,1.15701 }, { 2.95238,-0.0549855 }, { 1.99542,-0.214456 }, { -1.62097,3.13329 }, { -0.378359,2.09655 }, { -1.99542,-0.214456 }, { -0.750952,-1.13528 }, { -0.397807,-1.84626 }, { 0.750952,-1.13528 }, { -0.976178,-2.31605 }, { 0.471648,-3.44188 }, { 0.376369,-0.567049 }, { 1.70333,-0.910983 }, { 1.76636,0.952026 }, { 0.0948167,1.43406 }, { 0.930748,1.78001 }, { 2.26517,2.3001 } };

	std::vector<Point2D> basePoints = { {0.,0.}, { -2.10684, -1.84373 }, { 0.494872,-3.43051 }, { -1.8307,-0.772619 }, { -2.97035,-1.77037 }, { -0.869454,-2.66541 }, { -1.81226,-2.74913 }, { -0.93585,0.525783 }, { -2.84757,1.56598 }, { -0.140883,1.06415 }, { -3.2381,0.623172 }, { -1.70686,2.57839 }, { -3.36087,-0.827564 }, { 0.869454,-2.66541 }, { 1.81226,-2.74913 }, { 1.8307,-0.772618 }, { 0.794967,0.794967 }, { 2.78945,-0.204392 }, { 1.04125,-0.439927 }, { -0.366572,3.2381 }, { 1.39868,1.39868 }, { -1.04125,-0.439927 }, { 0.623172,3.2381 }, { 0.0, -1.13645 }, { 1.22688,2.50847 }, { 3.2381,0.623172 }, { 2.97035,-1.77037 }, { 2.78461,2.04139 }, { 0.0, -1.99554 }, { -2.78945,-0.204392 }, { -0.494871,-3.43051 }, { 2.10684,-1.84373 }, { 3.36087,-0.827563 }, { -0.635755,2.30225 }, { -1.87866,1.46859 }, { 2.50847,1.22688 }, { 2.04139,2.78461 } };


	//std::vector<Point2D> basePoints(101);
	//basePoints[0] = { 0.0, 0.0 };
	//for (size_t i = 0; i < 100; ++i)
	//	basePoints[i+1] = { cos(DPI * i / 100), sin(DPI * i / 100) };


	WakeDataBase controlPoints(*this);
	controlPoints.vtx.resize(getWake().vtx.size() * basePoints.size());

	for (size_t q = 0; q < getWake().vtx.size(); ++q)
		for (size_t i = 0; i < basePoints.size(); ++i)
			controlPoints.vtx[q * basePoints.size() + i].r() = getWake().vtx[q].r() + basePoints[i] * getWake().vtx[q].sigma();

	cuReserveDevMem((void*&)controlPoints.devVtxPtr, controlPoints.vtx.size() * sizeof(Vortex2D), 333);
	cuCopyWakeToDev(controlPoints.vtx.size(), controlPoints.vtx.data(), controlPoints.devVtxPtr, 444);
	cuReserveDevMem((void*&)controlPoints.devVelPtr, controlPoints.vtx.size() * 2 * sizeof(double), 555);

	auto& treeNut = getNonConstCuda().cntrTreeNut;
	treeNut->MemoryAllocate((int)getCuda().n_CUDA_wake * (int)basePoints.size());
	treeNut->Update((int)controlPoints.vtx.size(), controlPoints.devVtxPtr);
	treeNut->Build();

	std::vector<Point2D> veloVar(getWake().vtx.size() * basePoints.size(), Point2D{0.0, 0.0});
	std::vector<double> epsVar;
	velocity->GPUCalcConvVeloToSetOfPointsFromWake(treeNut, controlPoints, veloVar, epsVar, true, false);
	if (getNumberOfAirfoil() > 0)
		velocity->GPUCalcConvVelocityToSetOfPointsFromSheets(treeNut, controlPoints, veloVar);

	cuDeleteFromDev((void*&)controlPoints.devVtxPtr, 777);
	cuDeleteFromDev((void*&)controlPoints.devVelPtr, 888);

	for (size_t q = 0; q < getWake().vtx.size(); ++q)
	{
		double dv2 = 0.0;
		for (size_t i = 1; i < basePoints.size(); ++i)
			dv2 += (veloVar[q * basePoints.size()] - veloVar[q * basePoints.size() + i]).length2() * pow(getWake().vtx[q].sigma() / (controlPoints.vtx[q * basePoints.size() + i].r() - getWake().vtx[q].r()).length(), 0.6666667);
		
		dv2 *= 1.0 / (basePoints.size() - 1);
		nut[q] = 0.105 * 0.604 * wake->vtx[q].sigma() * sqrt(dv2);
	}
}
#endif 



// Вычисление скоростей панелей и интенсивностей присоединенных слоев вихрей и источников
void World2D::CalcPanelsVeloAndAttachedSheets()
{
	//вычисляем скорости панелей
	for (size_t i = 0; i < airfoil.size(); ++i)
		mechanics[i]->VeloOfAirfoilPanels(currentStep * passport.timeDiscretizationProperties.dt);

	//вычисляем интенсивности присоединенных слоев
	for (size_t i = 0; i < airfoil.size(); ++i)
		boundary[i]->ComputeAttachedSheetsIntensity();

}//CalcPanelsVeloAndAttachedSheets(...)




#include <random>

//Вычисляем новые положения вихрей (в пелене и виртуальных)
void World2D::MoveVortexes(std::vector<Point2D>& newPos, std::vector<double>* nutPtr)
{
	//getTimers().start("Move");	

	size_t nvt = wake->vtx.size();
	size_t nVirtVortex = 0;
	for (size_t i = 0; i < getNumberOfBoundary(); ++i)
		nVirtVortex += boundary[i]->virtualWake.vtx.size();




	size_t counter = wake->vtx.size() - nVirtVortex;
	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
		for (int i = 0; i < (int)boundary[bou]->virtualWake.vtx.size(); ++i)
		{
			wake->vtx[counter].g() = boundary[bou]->virtualWake.vtx[i].g();
			wake->vtx[counter].sigma() = getPassport().wakeDiscretizationProperties.sigma0;
			++counter;
		}
	}




	/////////////////////////////////////////
#if defined (SPH)

	const double c = 4.5;

	const auto& eps = velocity->wakeVortexesParams.epsastWake;

	double minEpsast = *std::min_element(eps.begin(), eps.end());
	double meanEpsast = std::accumulate(eps.begin(), eps.end(), 0.0) / eps.size();
	
	double uniformH = minEpsast + 0.2 * (meanEpsast - minEpsast);
	
	std::vector<std::vector<int>> neib(wake->vtx.size());
	for (auto& nb : neib)
		nb.reserve(100);
	
#pragma omp parallel for schedule(dynamic, 100)
	for (int i = 0; i < (int)wake->vtx.size(); ++i)
	{
		const double hati = c * eps[i];
		const double hati2 = sqr(hati);
		const int dist = std::ceil(hati / uniformH);
		
		Point2D vtx = wake->vtx[i].r();
		int cx = (int)(vtx[0] / uniformH);
		int cy = (int)(vtx[1] / uniformH);

		for (int j = 0; j < (int)wake->vtx.size(); ++j)
		{
			Point2D test = wake->vtx[j].r();
			int tx = (int)(test[0] / uniformH);
			int ty = (int)(test[1] / uniformH);

			if (abs(cx - tx) <= dist && abs(cy - ty) <= dist)
				if ((vtx - test).length2() < hati2)
					neib[i].push_back(j);
		}
	}

	//Ядро иранца
	//auto W = [](double xi, double R) {double t = xi / R; return (t > 1 ? 0.0 : (8.0 / ( PI * R * R)) * (0.75 * t * t - 1.5 * t + 0.75)); };
	//auto gW = [](double xi, double R) {double t = xi / R; return (t > 1 ? 0.0 : (12.0 / (PI * sqr(R * R))) * (1.0 - 1.0/t)); };
	
	//M4
	Point2D zeroVec{ 0.0, 0.0 };
	auto W = [](double xi, double R) {double t = xi / R; return (t > 1 ? 0.0 : 40.0 / (7.0 * PI * R * R) * (t > 0.5 ? 2.0 * cubPower(1.0 - t) : 1.0 - 6.0 * sqr(t) * (1.0 - t))); };	
	auto WW = [zeroVec](const Point2D& v, double R) {double t = v.length() / R; return (t > 1 ? zeroVec : v * (240.0 / (7.0 * PI * sqr(sqr(R))) * (t > 0.5 ? -sqr(1.0 - t) / t : -2.0 + 3.0 * t))); };

	//Poly6		
	//auto W = [](double xi, double R) {double t = xi / R; return (t > 1 ? 0.0 : 4.0 / (PI * R * R) * cubPower(1.0 - t * t)); };
	//auto gW = [](double xi, double R) {double t = xi / R; return (t > 1 ? 0.0 : -24.0 / (PI * sqr(R * R)) * sqr(1.0 - t * t)); };

	//Poly6	+ grad Spicky	
	//auto W = [](double xi, double R) {double t = xi / R; return (t > 1 ? 0.0 : 4.0 / (PI * R * R) * cubPower(1.0 - t * t)); };	
	//auto gW = [](double xi, double R) {double t = xi / R; return (t > 1 ? 0.0 : -30.0 / (PI * sqr(R * R)) * sqr(1.0 - t) / t); };

	//WendlandC6
	/*
	Point2D zeroVec{ 0.0, 0.0 };
	auto W = [](double xi, double R) {double t = xi / R; 
		return (t > 1 ? 0.0 : 78.0 / (7.0 * PI * sqr(R)) * sqr(sqr(sqr((1.0 - t)))) * (32.0 * cubPower(t)+ 25.0 * sqr(t) + 8.0*t + 1.0)); };
	auto WW = [zeroVec](const Point2D& v, double R) -> Point2D {double t = v.length() / R; 
		return (t > 1 ? zeroVec : v * (-1716.0 / (7.0 * PI * sqr(sqr(R))) * sqr(cubPower(1.0 - t)) * (1.0 - t) * (16.0 * sqr(t) + 7.0 * t + 1.0))); };
	*/

	//fp5
	/*
	auto W = [](double xi, double R) {
		double t = xi / R;
		if (t >= 1) return 0.0;
		if (t >= 2.0 / 3.0)
			return 1701.0 / (478.0 * PI) * sqr(sqr(1.0 - t)) * (1.0 - t);
		if (t >= 1.0 / 3.0)
			return 1701.0 / (478.0 * PI) * (sqr(sqr(1.0 - t)) * (1.0 - t) - 6.0 * sqr(sqr((2.0 / 3.0 - t))) * (2.0 / 3.0 - t));

		return 1701.0 / (478.0 * PI) * (sqr(sqr(1.0 - t)) * (1.0 - t) - 6.0 * sqr(sqr((2.0 / 3.0 - t))) * (2.0 / 3.0 - t) + 15.0 * sqr(sqr((1.0 / 3.0 - t))) * (1.0 / 3.0 - t));
		};

	auto gW = [](double xi, double R) {
		double t = xi / R;
		if (t >= 1) return 0.0;
		if (t >= 2.0 / 3.0)
			return -8505.0 / (478.0 * PI * R * R) * sqr(sqr(1.0 - t)) / t;
		if (t >= 1.0 / 3.0)
			return -8505.0 / (478.0 * PI * R * R) * (sqr(sqr(1.0 - t)) - 6.0 * sqr(sqr((2.0 / 3.0 - t)))) / t;

		return -8505.0 / (478.0 * PI * R * R) * (sqr(sqr(1.0 - t)) - 6.0 * sqr(sqr((2.0 / 3.0 - t))) + 15.0 * sqr(sqr((1.0 / 3.0 - t)))) / t;
		};
	*/


	std::vector<double> OIP(wake->vtx.size(), 0.0);
	std::vector<double> VIP(wake->vtx.size(), 0.0);
	std::vector<Point2D> GIP(wake->vtx.size(), { 0.0, 0.0 });
	std::vector<Point2D> DIFFV(wake->vtx.size(), { 0.0, 0.0 });
	
	
	const auto& vtx = wake->vtx;
	const auto& vel = velocity->wakeVortexesParams.convVelo;


//#ifdef TURB
//	std::vector<Point2D> grNu(wake->vtx.size(), { 0.0, 0.0 });
//	
//	//std::vector<Point2D> grVx(wake->vtx.size(), { 0.0, 0.0 });
//	//std::vector<Point2D> grVy(wake->vtx.size(), { 0.0, 0.0 });
//	std::vector<nummatrix<double,2,2>> grV(wake->vtx.size());
//	std::vector<Point2D> TP(wake->vtx.size(), { 0.0, 0.0 });
//	std::vector<double> RIP(wake->vtx.size(), 0.0);
//	std::vector<double> DELTAG(wake->vtx.size(), 0.0);
//#endif //TURB

	std::vector<nummatrix<double, 2, 2>> mtrL(wake->vtx.size(), { {0.0, 0.0}, {0.0, 0.0} });
	std::vector<std::vector<double>> weight(wake->vtx.size());
	
	
	std::vector<Point2D> sumDIFFV(wake->vtx.size(), {0.0, 0.0});
	std::vector<double> sumWeight(wake->vtx.size(), 0.0);

#pragma omp parallel for schedule(dynamic, 100)
	for (int i = 0; i < (int)vtx.size(); ++i)
	{
		const double hati = c * eps[i];

		double resultOIP = 0.0;
		//double resultIVIP = 0.0;

		const Point2D myPos = vtx[i].r();

		double I1 = 0.0, I0 = 1.0;


		weight[i].reserve(neib[i].size());

		for (auto& nbr : neib[i])
		{
			double dist2 = (myPos - vtx[nbr].r()).length2();
			double dist = sqrt(dist2);
			double w = W(dist, hati);
			I1 += vtx[nbr].g() * w;
			weight[i].push_back(w);
		}

		OIP[i] = I1 / I0;
	}

#pragma omp parallel for schedule(dynamic, 100)
	for (int i = 0; i < (int)vtx.size(); ++i)
	{
		const Point2D myPos = vtx[i].r();
		const double hati = c * eps[i];

		nummatrix<double, 2, 2> M = { {0.0, 0.0}, {0.0, 0.0} };
		GIP[i].toZero();

		for (auto& nbr : neib[i])
		{
			Point2D dr = myPos - vtx[nbr].r();
			Point2D grweight = WW(dr, hati);
			double vol = (vtx[nbr].g() / OIP[nbr]);
			if (vol != vol) 
				vol = 0.0;

			M += (dr | grweight) * vol;
			GIP[i] -= grweight * ((OIP[nbr] - OIP[i]) * vol);
		}

		double detM = M[0][0] * M[1][1] - M[1][0] * M[0][1];
		nummatrix<double, 2, 2> L{ {0.0, 0.0}, {0.0, 0.0} };
		if (fabs(detM) > 1e-2)
			L = nummatrix<double, 2, 2>{ { M[1][1], -M[0][1] }, { -M[1][0], M[0][0] } } *(1.0 / detM);
		
		GIP[i] = (L & GIP[i]);
		DIFFV[i] = GIP[i] * (1.0 / OIP[i]);
		//GIP[i] = -GIP[i];
	}

	//осредняем
	//*
//#pragma omp parallel for schedule(dynamic, 100)
	for (int i = 0; i < (int)vtx.size(); ++i)
		for (int j = 0; j < neib[i].size(); ++j)
			if ((vtx[i].r() - vtx[neib[i][j]].r()).length() <= 1.0 * eps[i])
			{
				sumDIFFV[i] += weight[i][j] * DIFFV[neib[i][j]];
				sumWeight[i] += weight[i][j];
			}

	//#pragma omp parallel for schedule(dynamic, 100)
	for (int i = 0; i < (int)vtx.size(); ++i)
	{		
		DIFFV[i] = sumDIFFV[i] * (1.0 / sumWeight[i]);
	}
	//*/

//#ifdef TURB
//#pragma omp parallel for schedule(dynamic, 100)
//	for (int i = 0; i < (int)vtx.size(); ++i)
//	{
//		Point2D myPos = vtx[i].r();
//		Point2D resultGradNu = { 0.0, 0.0 };
//		Point2D resultGradVx = { 0.0, 0.0 };
//		Point2D resultGradVy = { 0.0, 0.0 };
//
//		//if (i == 0)
//		//{
//		//	std::cout << "h = " << h << std::endl;
//		//	std::cout << "VIP = " << VIP[0] << std::endl;
//		//	std::cout << "my_coord = " << vtx[i].r() << std::endl;
//		//	std::cout << "mtrC = " << mtrC[i] << std::endl;
//		//}
//
//		for (auto& nbr : neib[i])
//		{
//			//if (i == 0)
//			//	std::cout << "nb_" << nbr << "_coord = " << vtx[nbr].r() << std::endl;
//
//			Point2D dr = myPos - vtx[nbr].r();
//			double dist = dr.length();			
//			double gweight = gW(dist, h);
//
//			if (nbr != i)
//			{
//				double nutNbr = (*nutPtr)[nbr];
//				double nutI = (*nutPtr)[i];
//				//double nutNbr = vtx[nbr].r() & Point2D { 1.0, 1.0 };
//				//double nutSlf = vtx[i].r() & Point2D { 1.0, 1.0 };
//
//				resultGradNu += (nutNbr - nutI) * (mtrL[i] & dr) * (gweight * /*vtx[nbr].g() / OIP[nbr]*/ VIP[nbr]);
//				resultGradVx += (vel[nbr][0] - vel[i][0]) * (mtrL[i] & dr) * (gweight * /*vtx[nbr].g() / OIP[nbr]*/ VIP[nbr]);
//				resultGradVy += (vel[nbr][1] - vel[i][1]) * (mtrL[i] & dr) * (gweight * /*vtx[nbr].g() / OIP[nbr]*/ VIP[nbr]);
//			}
//		}
//
//		grNu[i] = resultGradNu;
//		grV[i][0] = resultGradVx;
//		grV[i][1] = resultGradVy;
//
//		nummatrix<double, 2, 2> regGrV = grV[i];
//		regGrV[0][0] -= grV[i][1][1];
//		regGrV[0][1] += grV[i][1][0];
//		regGrV[1][0] += grV[i][0][1];
//		regGrV[1][1] -= grV[i][0][0];
//
//		regGrV *= 0.5;
//
//		regGrV[0][1] += 0.5 * OIP[i];
//		regGrV[1][0] -= 0.5 * OIP[i];
//
//		TP[i] = (regGrV & grNu[i]);
//	}//for i
//#endif //TURB


//#ifdef TURB
//
//#pragma omp parallel for schedule(dynamic, 100)
//	for (int i = 0; i < (int)vtx.size(); ++i)
//	{
//		Point2D myPos = vtx[i].r();
//
//		Point2D resultGradTPx = { 0.0, 0.0 };
//		Point2D resultGradTPy = { 0.0, 0.0 };
//
//		for (auto& nbr : neib[i])
//		{
//			Point2D dr = myPos - vtx[nbr].r();
//			double dist = dr.length();
//			double gweight = gW(dist, h);
//
//			if (nbr != i)
//			{			
//				resultGradTPx += (TP[nbr][0] - TP[i][0]) * (mtrL[i] & dr) * (gweight * /*vtx[nbr].g() / OIP[nbr]*/ VIP[nbr]);
//				resultGradTPy += (TP[nbr][1] - TP[i][1]) * (mtrL[i] & dr) * (gweight * /*vtx[nbr].g() / OIP[nbr]*/ VIP[nbr]);
//			}
//		}
//
//		RIP[i] = resultGradTPy[0] - resultGradTPx[1];
////*
//		double dG = RIP[i] * VIP[i]  * getPassport().timeDiscretizationProperties.dt;
//		//if (i == 0)
//		//	std::cout << "nVirtVortex = " << nVirtVortex << std::endl;
//		DELTAG[i] = 0.0;//dG;// (i < (wake->vtx.size() - nVirtVortex) ? dG : 0.0);
//			 
//		//if (fabs(dG / getWake().vtx[i].g()) <= 0.25)
//			getNonConstWake().vtx[i].g() += DELTAG[i];
//		//else
//		//{
//		//	if (dG * getWake().vtx[i].g() > 0)
//		//		getNonConstWake().vtx[i].g() *= 1.25;
//		//	else
//		//		getNonConstWake().vtx[i].g() /= 1.25;
//		//}
//
//		getNonConstWake().vtx[i].g() = std::clamp(getWake().vtx[i].g(), -1.5 * getPassport().wakeDiscretizationProperties.maxGamma, 1.5 * getPassport().wakeDiscretizationProperties.maxGamma);
//	}//for i
//	
//	wake->SaveScalarFields("SPH", VIP, OIP, *nutPtr, grNu, RIP, DELTAG);
//#endif //TURB

#endif //SPH

	

	//newPos.clear();
	newPos.resize(nvt);


	// 1. Инициализация генератора случайных чисел (ПСЧГ)
	//std::random_device rd; // Источник энтропии для инициализации
	//std::mt19937 gen(rd()); // Mersenne Twister, инициализированный случайным значением

	// 2. Определение параметров нормального распределения
	// μ (среднее) = 0.0, σ (стандартное отклонение) = 1.0
	//std::normal_distribution<> d(0.0, sqrt(2.0 * getPassport().physicalProperties.nu * getPassport().timeDiscretizationProperties.dt));


	//std::ofstream checkVel   (passport.dir + "/dbg/checkVel"    + std::to_string(getCurrentStep()) + ".txt");
	//std::ofstream checkVelDyn(passport.dir + "/dbg/checkVelDyn" + std::to_string(getCurrentStep()) + ".txt");
	
	//checkVel << "i x y epsast Vx Vy Wx Wy corrWx corrWy Om gOmx gOmy\n";
	//checkVelDyn << "i x y epsast Vx Vy Wx Wy corrWx corrWy Om gOmx gOmy\n";
	
//#pragma omp parallel for
	for (int i = 0; i < (int)wake->vtx.size(); ++i)
	{
#ifdef SPH
		//Point2D W = -(passport.physicalProperties.nu / OIP[i]) * GIP[i];
		Point2D W = -passport.physicalProperties.nu * DIFFV[i];
		
		checkVel << i << " " << wake->vtx[i].r()[0] << " " << wake->vtx[i].r()[1] << " " << velocity->wakeVortexesParams.epsastWake[i] << " " \
			<< velocity->wakeVortexesParams.convVelo[i][0] << " " << velocity->wakeVortexesParams.convVelo[i][1] << " " \
			<< W[0] << " " << W[1];
		
		if (W.length() > 1.5 * passport.physicalProperties.vRef)
		{
			std::cout << "i = " << i << ", W was " << W << std::endl;
			W.toZero();
			//W.normalize(1.5 * passport.physicalProperties.vRef);
		}

		checkVel << " " << W[0] << " " << W[1] << " " << OIP[i] << " " << GIP[i][0] << " " << GIP[i][1] << "\n";

#else
		Point2D W = velocity->wakeVortexesParams.diffVelo[i] * (nutPtr ? (1.0 + (*nutPtr)[i] / passport.physicalProperties.nu) : 1.0);
		
		const auto& diff = velocity->wakeVortexesParams;
		Point2D grad = diff.I2[i] * (1.0 / (diff.epsastWake[i] * diff.I0[i])) - diff.I3[i] * (diff.I1[i] / sqr(diff.I0[i]));
		
		//checkVelDyn << i << " " << wake->vtx[i].r()[0] << " " << wake->vtx[i].r()[1] << " " << velocity->wakeVortexesParams.epsastWake[i] << " " \
		//	<< velocity->wakeVortexesParams.convVelo[i][0] << " " << velocity->wakeVortexesParams.convVelo[i][1] << " " \
		//	<< W[0] << " " << W[1] << " " << W[0] << " " << W[1] << " " << diff.I1[i]/diff.I0[i] << " "  << grad[0] << " " << grad[1] << "\n";
#endif

		newPos[i] = wake->vtx[i].r() + (velocity->wakeVortexesParams.convVelo[i] + W + getV0()) * passport.timeDiscretizationProperties.dt;

		//Chorin Random Walk method
		//newPos[i] = wake->vtx[i].r() + (velocity->wakeVortexesParams.convVelo[i] + getV0()) * passport.timeDiscretizationProperties.dt + Point2D{ d(gen), d(gen) };
	}
	//checkVel.close();
	//checkVelDyn.close();


	//for (int i = 0; i < (int)wake->vtx.size(); ++i)
	//{
	//	std::cout << i << " " << wake->vtx[i].r() << " " << \
	//		velocity->wakeVortexesParams.convVelo[i] << " " << \
	//		velocity->wakeVortexesParams.diffVelo[i] << " " << \
	//		getV0() << "\n";
	//}




	//getTimers().stop("Move");	
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
	getTimers().start("MatRhs");
		
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
	   
	getTimers().stop("MatRhs");	
}//FillIQ()

void World2D::CalcAndSolveLinearSystem()
{
	if (airfoil.size() > 0)
	{
		ReserveMemoryForMatrixAndRhs();

#if defined(__CUDACC__) || defined(USE_CUDA)
		cuda.setAccelCoeff(passport.physicalProperties.accelCft(getCurrentTime()));
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
		const int linSystemScheme = passport.numericalSchemes.linearSystemSolver.second;		

		if (linSystemScheme == 0) //Gauss
		{
			if (currentStep == 0)
			{
#ifdef BRIDGE
				useInverseMatrix = true;
#endif //BRIDGE

#ifdef INITIAL
				useInverseMatrix = (
					((mechanics.size() == 1) && (!mechanics.front()->isDeform))
					||
					(mechanics.size() > 1 && !isAnyMovableOrDeformable())
					);
#endif //INITIAL
			}
		}//if Gauss

		if (linSystemScheme == 0 || linSystemScheme == 1 || (linSystemScheme == 2 && isAnyMovableOrDeformable()))
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
			for (size_t i = 0; i < matrReord.rows(); ++i)
			{
				for (size_t j = 0; j < matrReord.cols(); ++j)
					of << matrReord(i, j) << " ";
				of << std::endl;
			}
			of.close();
		}

		{
			std::stringstream ss;
			ss << "rhs-" << currentStep;
			std::ofstream of(passport.dir + "dbg/" + ss.str());
			of.precision(16);
			for (size_t i = 0; i < rhsReord.rows(); ++i)
			{
				of << rhsReord(i) << std::endl;
			}
			of.close();
		}

		//*/

		SolveLinearSystem();
		
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
				wake->vtx.push_back(Vortex2D{ boundary[bou]->virtualWake.vtx[v].r(), 0.0, getPassport().wakeDiscretizationProperties.sigma0});

		//for (size_t v = 0; v < wake->vtx.size(); ++v)
		//{
		//	if (fabs(wake->vtx[v].g()) > 1.0)
		//		std::cout << "Error gam! " << "v = " << v << ", gam[v] = " << wake->vtx[v].g() << std::endl;
		//}

		//for (size_t bou = 0; bou < boundary.size(); ++bou)
		//	for (size_t v = 0; v < airfoil[bou]->getNumberOfPanels(); ++v)
		//	{
		//		if (fabs(boundary[bou]->sheets.freeVortexSheet(v, 0)) > 1000.0)
		//			std::cout << "Error gamma sheet! " << "bou = " << bou << ", v = " << v << ", gam[v] = " << boundary[bou]->sheets.freeVortexSheet(v, 0) << std::endl;
		//	}

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
	}
}//CalcAndSolveLinearSystem()

void World2D::WakeAndAirfoilsMotion(bool dynamics, std::vector<double>* nutPtr)
{
	std::vector<Point2D> newPos;

	MoveVortexes(newPos, nutPtr);

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

	oldAirfoil.resize(0);
	for (auto& afl : airfoil)
	{
		if (dynamic_cast<AirfoilRigid*>(afl.get()))
			oldAirfoil.emplace_back(new AirfoilGeometry(*afl));

		if (dynamic_cast<AirfoilDeformable*>(afl.get()))
			oldAirfoil.emplace_back(new AirfoilGeometry(*afl));
		

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

#if defined(__CUDACC__) || defined(USE_CUDA)
	cuda.RefreshAfls(2);
#endif

	for (auto& bou : boundary)
		bou->virtualWake.vtx.clear();
	

	CheckInside(newPos, oldAirfoil);
		

	//передача новых положений вихрей в пелену	
	for (size_t i = 0; i < wake->vtx.size(); ++i)
	{
		wake->vtx[i].r() = newPos[i];
		//wake->vtx[i].sigma() = getPassport().wakeDiscretizationProperties.sigma0;
	}

//	getWake().SaveKadrVtk();
	
}//WakeAndAirfoilsMotion()

// Возврат признака того, что хотя бы один из профилей подвижный
bool World2D::isAnyMovable() const
{
	return std::any_of(getMechanicsVector().begin(), getMechanicsVector().end(),
		[](const std::unique_ptr<Mechanics>& m) {return (m->isMoves); });
}

// Возврат признака того, что хотя бы один из профилей подвижный или деформируемый
bool World2D::isAnyMovableOrDeformable() const
{
	return std::any_of(getMechanicsVector().begin(), getMechanicsVector().end(),
		[](const std::unique_ptr<Mechanics>& m) {return (m->isMoves || m->isDeform); });
}