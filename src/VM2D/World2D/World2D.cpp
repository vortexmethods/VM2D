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

	std::vector<std::string> timerLabels = { "Step", "MatRhs", "Solve", "ConvVel", "DiffVel", "Force", "VelPres", "Inside", "Restr", "Save"};
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

	if (getPassport().wakeDiscretizationProperties.eps == 0)
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
		getNonConstPassport().wakeDiscretizationProperties.eps = eps;
		getNonConstPassport().wakeDiscretizationProperties.eps2 = eps * eps;
		info('i') << "eps = " << getPassport().wakeDiscretizationProperties.eps << " is calculated automatically" << std::endl;
	}

	if (getPassport().wakeDiscretizationProperties.epscol == 0)
	{
		getNonConstPassport().wakeDiscretizationProperties.epscol = (2.0 / 3.0) * getPassport().wakeDiscretizationProperties.eps;

		info('i') << "epscol = " << getPassport().wakeDiscretizationProperties.epscol << " is calculated automatically" << std::endl;
	}

	for (size_t afl = 0; afl < passport.airfoilParams.size(); ++afl)
		if (passport.airfoilParams[afl].chord == 0)
		{
			auto& prm = getNonConstPassport().airfoilParams[afl];
			prm.chord = (prm.initialGab.second[0] - prm.initialGab.first[0]) * prm.scale[0];
			info('i') << "airfoil #" << afl << " chord = " << prm.chord << " is calculated automatically" << std::endl;
		}

	//2026-03-28
	//считываем массив точек для подсчета и вывода поля скоростей и давлений
	measureVP.reset(new MeasureVP(*this));

	if (passport.timeDiscretizationProperties.saveVPstep != 0)
		measureVP->ReadPointsFromFile(passport.dir);

        //optimizer
	//2026-03-28
	//std::vector<Point2D> VPPoints, VPHistory;
	//for (size_t a = 0; a < getNumberOfAirfoil(); ++a)
	//{
	//	const auto& afl = *airfoil[a];
	//	for (size_t p = 0; p < afl.getNumberOfPanels(); ++p)
	//		VPHistory.push_back(0.5 * (afl.getR(p) + afl.getR(p + 1)) + afl.nrm[p] * passport.airfoilParams[a].chord * 0.01);
	//}
	//if (passport.timeDiscretizationProperties.saveVPstep != 0)
	//	measureVP->SetPoints(VPPoints, VPHistory);


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
				treeWake.Update((int)getWake().vtx.size(), getWake().devVtxPtr, getPassport().wakeDiscretizationProperties.eps);
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
				//info('i') << "Goodbye!" << std::endl;
				//exit(0);
				currentTime = getPassport().timeDiscretizationProperties.timeStop;
			}

			timerInitialBuild.start();
			if (nTotPan > 0 && getPassport().numericalSchemes.velocityComputation.second == 1)
			{
				auto& afl = getAirfoil(0);
#if USE_CUDA
				getCuda().inflTreePnlVortex->UpdatePanelFreeAndAttachedVortexIntensity(afl.devFreeVortexSheetPtr, afl.devFreeVortexSheetLinPtr, afl.devAttachedVortexSheetPtr, afl.devAttachedVortexSheetLinPtr);
				getCuda().inflTreePnlVortex->UpwardTraversal(multipoleOrder);
#endif
			}
			timerInitialBuild.stop();

			//Вычисление скоростей вихрей: для тех, которые в следе, и виртуальных, а также в точках wakeVP
			CalcVortexVelo();

//#include "gammaCirc.h"						

			//Расчет и сохранение поля давления
			if (ifDivisible(passport.timeDiscretizationProperties.saveVPstep))
			//optimizer
				//if ((getCurrentStep() >= 0 && getCurrentStep() <= 12000) || getCurrentStep()==0)  //2026-03-28
			{
				const double& vRef = getPassport().physicalProperties.vRef;
				//scaleV = 1.0 / vRef;
				double scaleP = 1.0 / (0.5 * getPassport().physicalProperties.rho * sqr(vRef));

				measureVP->CalcPressure();
				
				{
					MechanicsDeformable* ptr = dynamic_cast<MechanicsDeformable*>(mechanics[0].get());
					if (measureVP->getTotalNumberOfRealPoints() > 0)
						measureVP->SaveVP();

					if (ptr && !ptr->beam->fsi)

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
						phiFile << std::endl;
					}
					else
						phiFile.open(getPassport().dir + "phiFile.csv", std::ios_base::app);					
					
					phiFile << currentTime;
					for (int p = 0; p < ptr->beam->R; ++p)
						phiFile << "," << ptr->beam->phi(p, 0);
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
			WakeAndAirfoilsMotion(true);
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
			cntrTreePnt.Update((int)newPos.size(), devNewpos_ptr, 0.0);
			cntrTreePnt.Build();

			BHcu::treeClosestPanelToPointsCalculationWrapper(auxTree, cntrTreePnt, getWake().devNearestPanelPtr, true, getAirfoil(0).devPsnPtr);

			cudaMemcpy(hit.data(), getWake().devNearestPanelPtr, newPos.size() * sizeof(int), cudaMemcpyDeviceToHost);
			cudaFree(devNewpos_ptr);
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

/*
	if (currentStep == 0)
	{
		std::ofstream fileMatrix(passport.dir + "/dbg/matrix.txt");
		//int nx, ny;
		//fileMatrix >> nx;
		//fileMatrix >> ny;
		fileMatrix.precision(16);
		for (int i = 0; i < matrReord.rows(); ++i)
		{
			for (int j = 0; j < matrReord.cols(); ++j)
			{
				fileMatrix << matrReord(i, j) << " ";
			}
			fileMatrix << std::endl;
		}
		fileMatrix.close();
	}
	
		

	{
		std::ofstream fileMatrix(passport.dir + "/dbg/rhs"+std::to_string(currentStep)+".txt");
		//int nx, ny;
		//fileMatrix >> nx;
		//fileMatrix >> ny;
		fileMatrix.precision(16);
		for (int i = 0; i < rhsReord.size(); ++i)
		{			
			fileMatrix << rhsReord(i) << std::endl;
		}
		fileMatrix.close();
	}//*/

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

			GMRES_Direct(*this, nFullVars, (int)getNumberOfBoundary(), Gmatr, Grhs, Gpos, Gvsize, Ggam, GR);

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
			GMRES(*this, Ggam, GR, GGrhs, GrhsReg, niter, linScheme);
			time_GMRES += omp_get_wtime();
			std::cout << "Time_GMRES = " << time_GMRES << std::endl;
 
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
	if(!(currentStep % 1))
	{
		VMlib::CreateDirectory(passport.dir, "dbg");
		std::ostringstream sss;
		sss << "prmWake";
		sss << currentStep;
		std::ofstream prmtFile(passport.dir + "dbg/" + sss.str());
		prmtFile.precision(11);
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
			<< "\n";

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
		//if (currentStep==3) exit(-123);
	}
//*/

	getTimers().stop("Save");
	
}//CalcVortexVelo()



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
void World2D::MoveVortexes(std::vector<Point2D>& newPos)
{
	//getTimers().start("Move");	

	size_t nvt = wake->vtx.size();
	size_t nVirtVortex = 0;
	for (size_t i = 0; i < getNumberOfBoundary(); ++i)
		nVirtVortex += boundary[i]->virtualWake.vtx.size();

	//newPos.clear();
	newPos.resize(nvt);


	// 1. Инициализация генератора случайных чисел (ПСЧГ)
	//std::random_device rd; // Источник энтропии для инициализации
	//std::mt19937 gen(rd()); // Mersenne Twister, инициализированный случайным значением

	// 2. Определение параметров нормального распределения
	// μ (среднее) = 0.0, σ (стандартное отклонение) = 1.0
	//std::normal_distribution<> d(0.0, sqrt(2.0 * getPassport().physicalProperties.nu * getPassport().timeDiscretizationProperties.dt));


#pragma omp parallel for
	for (int i = 0; i < (int)wake->vtx.size(); ++i)
	{
		newPos[i] = wake->vtx[i].r() 
			+ (velocity->wakeVortexesParams.convVelo[i] +
				velocity->wakeVortexesParams.diffVelo[i] +
				getV0()) * passport.timeDiscretizationProperties.dt;
		
		//newPos[i] = wake->vtx[i].r() + (velocity->wakeVortexesParams.convVelo[i] + getV0()) * passport.timeDiscretizationProperties.dt + Point2D{ d(gen), d(gen) };
	}

	//for (int i = 0; i < (int)wake->vtx.size(); ++i)
	//{
	//	std::cout << i << " " << wake->vtx[i].r() << " " << \
	//		velocity->wakeVortexesParams.convVelo[i] << " " << \
	//		velocity->wakeVortexesParams.diffVelo[i] << " " << \
	//		getV0() << "\n";
	//}


	size_t counter = wake->vtx.size() - nVirtVortex;
	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
		for (int i = 0; i < (int)boundary[bou]->virtualWake.vtx.size(); ++i)
		{
			wake->vtx[counter].g() = boundary[bou]->virtualWake.vtx[i].g();
			++counter;
		}
	}

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
				wake->vtx.push_back(Vortex2D{ boundary[bou]->virtualWake.vtx[v].r(), 0.0 });

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
		wake->vtx[i].r() = newPos[i];

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