/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
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
\version 1.10
\date 17 мая 2021 г.
*/

#include "World2D.h"

#include "Airfoil2DRect.h"
#include "Airfoil2DCurv.h"

#include "Boundary2DConstLayerAver.h"
#include "Boundary2DLinLayerAver.h"
#include "Boundary2DVortexCollocN.h"

#include "MeasureVP2D.h"

#include "Mechanics2DRigidImmovable.h"
#include "Mechanics2DRigidGivenLaw.h"
#include "Mechanics2DRigidOscillPart.h"
#include "Mechanics2DRigidRotatePart.h"

#include "Parallel.h"

#include "Passport2D.h"

#include "StreamParser.h"


#include "Velocity2DBiotSavart.h"
#include "Velocity2DBarnesHut.h"

#include "Wake2D.h"


using namespace VM2D;

//Конструктор
World2D::World2D(const VMlib::PassportGen& passport_, const VMlib::Parallel& parallel_) :
	WorldGen(passport_, parallel_),
	passport(dynamic_cast<const Passport&>(passport_)),
	cuda(Gpu(*this))
{
	if (parallel.myidWork == 0)
	{
		std::stringstream ss;
		ss << "#" << passport.problemNumber << " (" << passport.problemName << ")";		
		info.assignStream(defaults::defaultWorld2DLogStream, ss.str());
	}
		
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
	if (passport.timeDiscretizationProperties.saveVP != 0)
		measureVP->ReadPointsFromFile(passport.dir);

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

	for (size_t i = 0; i < passport.airfoilParams.size(); ++i)
	{
		switch (passport.numericalSchemes.panelsType.second)
		{
		case 0:			
			airfoil.emplace_back(new AirfoilRect(*this, i));
			break;
		case 1:
			airfoil.emplace_back(new AirfoilCurv(*this, i));
			break;
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
			
		case 11:			
			//boundary.emplace_back(new BoundaryVortColl(*this, i));			
			info('e') << "BoundaryVortColl is not implemented now! " << std::endl;
			exit(1);
			break;

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


//Основная функция выполнения одного шага по времени
void World2D::Step()
{	
	try {
		//double tt1, tt2;
		
		//Очистка статистики
		getTimestat().ToZero();

		//Засечка времени в начале шага
		getTimestat().timeWholeStep.first += omp_get_wtime();

		CalcPanelsVeloAndAttachedSheets();

		measureVP->Initialization();
		BuildAllTrees(); //Строго после measureVP->Initialization
				
		CalcAndSolveLinearSystem();

		wake->WakeSynchronize();

		for (auto& bou : boundary)
			bou->virtualWake.WakeSynchronize();
		
		//tt1 = omp_get_wtime();
		//Вычисление скоростей вихрей: для тех, которые в следе, и виртуальных, а также в точках wakeVP
		CalcVortexVelo();
		//tt2 = omp_get_wtime();
		//std::cout << tt2 - tt1 << std::endl;

		
//#include "gammaCirc.h"


		//Расчет и сохранение поля давления
		if (ifDivisible(passport.timeDiscretizationProperties.saveVP))
		{
			measureVP->CalcPressure();
			measureVP->SaveVP();
		}		
		

		//Вычисление сил, действующих на профиль и сохранение в файл	
		if (parallel.myidWork == 0)
			for (auto& mech : mechanics)
			{
				mech->GetHydroDynamForce();
				mech->GenerateForcesString();
				mech->GeneratePositionString();
			}

//		boundary[0]->virtualWake.SaveKadrVtk("VirtWake");
		
		//Движение вихрей (сброс вихрей + передвижение пелены)
		WakeAndAirfoilsMotion();		
		

		wake->Restruct();
		wake->SaveKadrVtk();
		

		//Засечка времени в конце шага
		getTimestat().timeWholeStep.second += omp_get_wtime();


		if (parallel.myidWork == 0)
		{
			info('i') << "Step = " << currentStep \
				<< " PhysTime = " << passport.physicalProperties.getCurrTime() \
				<< " StepTime = " << Times::dT(getTimestat().timeWholeStep) << std::endl;				

			getTimestat().GenerateStatString();
		}

		passport.physicalProperties.addCurrTime(passport.timeDiscretizationProperties.dt);
		currentStep++;


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
		
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
		wake->Inside(newPos, *airfoil[afl], mechanics[afl]->isMoves, *oldAirfoil[afl]);

	getTimestat().timeCheckInside.second += omp_get_wtime();
}

//Решение системы линейных алгебраических уравнений
void World2D::SolveLinearSystem()
{
	getTimestat().timeSolveLinearSystem.first += omp_get_wtime();
	
	/// \todo Нет возможности считать взаимно движущиеся профили
	//sol = matr.partialPivLu().solve(rhs);

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

	if (currentStep == 0)
		invMatr = matr.inverse();


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

	sol = invMatr*rhs;

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
		size_t nVars;

		if (parallel.myidWork == 0)
		{			
			nVars = boundary[bou]->GetUnknownsSize();

			if (currentStep == 0)
			{
				locMatr.resize(nVars, nVars);
				locLastLine.resize(nVars);
				locLastCol.resize(nVars);
			}
		}
		
		if (currentStep == 0 || mechanics[bou]->isDeform)
		{
			if (parallel.myidWork == 0)
			{
				boundary[bou]->FillMatrixSelf(locMatr, locLastLine, locLastCol);
			}
		}

		//размазываем матрицу
		if (parallel.myidWork == 0)
		{
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

			/// \todo пока матрица считается только на первом шаге
			if (currentStep == 0)
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

		}// if (parallel.myidWork == 0)
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

// Построение дерева для точек типа type
void World2D::BuildTree(const PointType type) const
{
	/*

	//const int nLevelWake = 14;
	//const int nLevelVP = 14;
	//const int nLevelSorces = 5;
	//const int nLevelSheets = 8;

	const int nLevelWake = 14;
	const int nLevelVP = 1;
	const int nLevelSources = 1;
	const int nLevelSheets = 10;



	switch (type)
	{
	case PointType::wake:
	{
		if (wake->vtx.size() > 0)
			treeWake.reset(new Tree(*this, *wake, PointType::wake, nLevelWake));
		else
			treeWake = nullptr;
		break;
	}	
	case PointType::wakeVP:
	{
		if (measureVP->getWakeVP().vtx.size() > 0)
			treeVP.reset(new Tree(*this, measureVP->getWakeVP(), PointType::wakeVP, nLevelVP));
		else
			treeVP = nullptr;
		break;
	}
	case PointType::sourceWake:
	{
		if (source->vtx.size() > 0)
			treeSourcesWake.reset(new Tree(*this, *source, PointType::sourceWake, nLevelSources));
		else
			treeSourcesWake = nullptr;
		break;
	}

	case PointType::sheetGam:
	{
		if (boundary.size() > 0)
			treeSheetsGam.reset(new Tree(*this, PointType::sheetGam, nLevelSheets));
		else
			treeSheetsGam = nullptr;
		break;
	}
	case PointType::source:
	{
		if (boundary.size() > 0)
			treeSheetsSource.reset(new Tree(*this, PointType::source, nLevelSheets));
		else
			treeSheetsSource = nullptr;
		break;
	}

	default:
		break;
	}
	*/
}//BuildTree(...)


void World2D::BuildAllTrees()
{
	/*
	if (passport.numericalSchemes.velocityComputation.second == 1)
	{
		getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();

		BuildTree(PointType::wake);
		if (treeWake)
			treeWake->rootCell.CalculateCellsParams();

		BuildTree(PointType::sourceWake);
		if (treeSourcesWake)
			treeSourcesWake->rootCell.CalculateCellsParams();

		if (getNumberOfAirfoil() > 0)
		{
			BuildTree(PointType::sheetGam);
			BuildTree(PointType::source);

			treeSheetsSource->rootCell.CalculateCellsParams();
		}
		getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();

		getTimestat().timeVP.first += omp_get_wtime();
		if ((getPassport().timeDiscretizationProperties.saveVP > 0) && (!(getCurrentStep() % getPassport().timeDiscretizationProperties.saveVP)))
			BuildTree(PointType::wakeVP);
		getTimestat().timeVP.second += omp_get_wtime();
	}
	*/

}//BuildAllTrees()

// Вычисление скоростей (и конвективных, и диффузионных) вихрей (в пелене и виртуальных), а также в точках вычисления VP 
void World2D::CalcVortexVelo()
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
	velocity->CalcConvVelo();

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
/*	//Сохранение всех параметров для вихрей в пелене
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
*/

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

	newPos.clear();
	newPos.reserve(nvt);	

	if (parallel.myidWork == 0)
		for (size_t i = 0; i < wake->vtx.size(); ++i)
		{
			newPos.push_back(wake->vtx[i].r() \
				+ (velocity->wakeVortexesParams.convVelo[i] + 
				   velocity->wakeVortexesParams.diffVelo[i] + 
				   passport.physicalProperties.V0())*passport.timeDiscretizationProperties.dt);			
		}

	if (parallel.myidWork == 0)
		for (size_t bou = 0; bou < boundary.size(); ++bou)
			for (size_t i = 0; i < boundary[bou]->virtualWake.vtx.size(); ++i)
			{
				wake->vtx.push_back(boundary[bou]->virtualWake.vtx[i]);
				newPos.push_back(boundary[bou]->virtualWake.vtx[i].r() \
					+ (velocity->virtualVortexesParams[bou].convVelo[i] +
					   velocity->virtualVortexesParams[bou].diffVelo[i] +
					   passport.physicalProperties.V0())*passport.timeDiscretizationProperties.dt);
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
			if (parallel.myidWork == 0)
			{
				boundary[bou]->FillIQSelf(IQ[bou][bou]);
			}
		}

		if (parallel.myidWork == 0)
		{
			for (size_t oth = 0; oth < boundary.size(); ++oth)
			{
				//\todo пока матрица считается только на первом шаге

#ifdef INITIAL
				if (currentStep == 0 || mechanics[bou]->isMoves || mechanics[oth]->isMoves)				
				{
					//size_t nVarsOther = boundary[oth]->GetUnknownsSize();
					if (bou != oth)
						boundary[bou]->FillIQFromOther(*boundary[oth], IQ[bou][oth]);
				}// if (currentStep == 0 || mechanics[oth]->isMoves)
#endif
				
#ifdef BRIDGE
				if (currentStep == 0)
				{
					size_t nVarsOther = boundary[oth]->GetUnknownsSize();
					if (bou != oth)
						boundary[bou]->FillIQFromOther(*boundary[oth], IQ[bou][oth]);
				}// if (currentStep == 0 || mechanics[oth]->isMoves)
#endif

			}// for oth
		}// if (parallel.myidWork == 0)
	}// for bou
	   
	getTimestat().timeFillMatrixAndRhs.second += omp_get_wtime();
}//FillIQ()

void World2D::CalcAndSolveLinearSystem()
{
	if (airfoil.size() > 0)
	{
		if (parallel.myidWork == 0)
			ReserveMemoryForMatrixAndRhs();

#if defined(__CUDACC__) || defined(USE_CUDA)
		cuda.setAccelCoeff(passport.physicalProperties.accelCft());
		cuda.setMaxGamma(passport.wakeDiscretizationProperties.maxGamma);

		int sch = passport.numericalSchemes.boundaryCondition.second;

		if( (sch == 0) || (sch == 1) || (sch == 10) )
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
		
		{
			std::stringstream ss;
			ss << "rhs-" << currentStep;
			std::ofstream of(passport.dir + "dbg/" + ss.str());
			for (size_t i = 0; i < rhs.size(); ++i)
				of << rhs(i) << std::endl;
			of.close();
		}
		exit(10);
		*/

		if (parallel.myidWork == 0)
			SolveLinearSystem();

	/*	{
			std::stringstream ss;
			ss << "sol-" << currentStep;
			std::ofstream of(passport.dir + "dbg/" + ss.str());
			for (size_t i = 0; i < sol.size(); ++i)
				of << sol(i) << std::endl;
			of.close();
		}
	*/

		getTimestat().timeOther.first += omp_get_wtime();
		if (parallel.myidWork == 0)
		{
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
				for (size_t q = 0; q < boundary[0]->virtualWake.vtx.size(); ++q)
				{
					addMass += boundary[0]->virtualWake.vtx[q].g() * boundary[0]->virtualWake.vtx[q].r().kcross();
				}
				std::cout << "AddMass = " << addMass << std::endl;
				//exit(-42);
			}
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

	
	for (auto& afl : airfoil)
	{
		oldAirfoil.resize(0);

		switch (passport.numericalSchemes.panelsType.second)
		{
		case 0:
			oldAirfoil.emplace_back(new AirfoilRect(*afl));
			break;
		case 1:
			oldAirfoil.emplace_back(new AirfoilCurv(*afl));
			break;
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
	if (parallel.myidWork == 0)
		for (size_t i = 0; i < wake->vtx.size(); ++i)
			wake->vtx[i].r() = newPos[i];
//	getWake().SaveKadrVtk();
	getTimestat().timeOther.second += omp_get_wtime();
}//WakeAndAirfoilsMotion()

