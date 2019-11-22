/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
| VM is distributed in the hope that it will be useful, but WITHOUT           |
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
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7   
\date 22 ноября 2019 г.
*/

#include "World2D.h"

#include "Airfoil2DRect.h"
#include "Airfoil2DCurv.h"

#include "Boundary2DConstLayerAver.h"

#include "Mechanics2DRigidImmovable.h"
#include "Mechanics2DRigidGivenLaw.h"
#include "Mechanics2DRigidOscillPart.h"

#include "Velocity2DBiotSavart.h"
#include "Velocity2DBarnesHut.h"

#include "MeasureVP2D.h"

#include "Parallel.h"

#include "Passport2D.h"

#include "StreamParser.h"

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

	//airfoil.resize(GetPassport().airfoilParams.size());
	//boundary.resize(GetPassport().airfoilParams.size());
	//mechanics.resize(GetPassport().airfoilParams.size());


	//считываем массив точек для подсчета и вывода поля скоростей и давлений
	if (passport.timeDiscretizationProperties.saveVP == 0)
		measureVP.reset(nullptr);
	else
	{
		measureVP.reset(new MeasureVP(*this));
		measureVP->ReadPointsFromFile(passport.dir);
	}

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

//			info('e') << "AirfoilCurv is not implemented now! " << std::endl;
//			exit(1);
			break;
		}

		airfoil[i]->ReadFromFile(passport.airfoilsDir);	//Считываем из каталога с коллекцией профилей
		
		switch (passport.numericalSchemes.boundaryCondition.second)
		{
		case 0:
			boundary.emplace_back(new BoundaryConstLayerAver(*this, i));
			break;
		
		case 10:			
			//boundary.emplace_back(new BoundaryMDV(*this, i));			
			info('e') << "BoundaryMDV is not implemented now! " << std::endl;
			exit(1);
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
		/*
		case 3:
			mechanics.emplace_back(new MechanicsRigidOscillMon(*this, i));
			break;
		case 4:
			mechanics.emplace_back(new MechanicsRigidRotateMon(*this, i));
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

	timePeriod& time = getTimestat().timeWholeStep;

	//Очистка статистики
	getTimestat().ToZero();

	//Засечка времени в начале шага
	time.first = omp_get_wtime();

	const double& dt = passport.timeDiscretizationProperties.dt;


	getTimestat().timeOther.first += omp_get_wtime();
	//вычисляем скорости панелей
	for (size_t i = 0; i < airfoil.size(); ++i)
		mechanics[i]->VeloOfAirfoilPanels(currentStep * dt);
		
	//вычисляем интенсивности присоединенных слоев
	for (size_t i = 0; i < airfoil.size(); ++i)
		boundary[i]->ComputeAttachedSheetsIntensity();

	getTimestat().timeOther.second += omp_get_wtime();

	// если быстрый метод
	if (passport.numericalSchemes.velocityComputation.second == 1)
	{
		auto velBH = dynamic_cast<VelocityBarnesHut*>(velocity.get());
		// строим все деревья
		velBH->BuildTrees(PointType::wake);
		velBH->BuildTrees(PointType::sourceWake);
		velBH->BuildTrees(PointType::wakeVP);

		if (airfoil.size() > 0)
			velBH->BuildTrees(PointType::sheetGam);
	}
	if (airfoil.size() > 0)
	{
		if (parallel.myidWork == 0)
			ReserveMemoryForMatrixAndRhs(getTimestat().timeReserveMemoryForMatrixAndRhs);

#if defined(__CUDACC__) || defined(USE_CUDA)
		cuda.setAccelCoeff(passport.physicalProperties.accelCft());
		cuda.setMaxGamma(passport.wakeDiscretizationProperties.maxGamma);
		cuda.setMinEpsAst2(sqr(2.0*passport.wakeDiscretizationProperties.epscol));

		cuda.RefreshWake();
		cuda.RefreshAfls();	
		cuda.RefreshVirtualWakes();
#endif

//		double timeStartFillIQ = omp_get_wtime();
		if (passport.numericalSchemes.panelsType.second == 0)
			FillIQ();
//		double timeEndFillIQ = omp_get_wtime();
//		info('t') << "FillIQ_time = " << timeEndFillIQ - timeStartFillIQ << std::endl;

//		double timeStartFillMatrixAndRhs = omp_get_wtime();
		FillMatrixAndRhs(getTimestat().timeFillMatrixAndRhs);
//		double timeEndFillMatrixAndRhs = omp_get_wtime();
//		info('t') << "FillMatrixAndRhs_time = " << timeEndFillMatrixAndRhs - timeStartFillMatrixAndRhs << std::endl;

		if (parallel.myidWork == 0)
			SolveLinearSystem(getTimestat().timeSolveLinearSystem);

		//if (parallel.myidWork == 0)
		//{
		//for (int i = 0; i < sol.size(); ++i)
		//{
		//info('t') << i << " " << sol(i) << std::endl;
		//}
		//}


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

			
//			std::ostringstream ss;
//			ss << "solution_";
//			ss << currentStep;
//			std::ofstream solFile(ss.str());
//			VMlib::SaveToStream(sol, solFile);
//			solFile.close();
//			//exit(1);
//			

		}
		getTimestat().timeOther.second += omp_get_wtime();

	}

	getTimestat().timeOther.first += omp_get_wtime();
	wake->WakeSynchronize();

	for (int bou = 0; bou < boundary.size(); ++bou)
	{
		boundary[bou]->VirtualWakeSynchronize();
	}
	getTimestat().timeOther.second += omp_get_wtime();

	//Инициализация всех массивов для вычисления скорости и давления в точках
	if ((passport.timeDiscretizationProperties.saveVP > 0) && (!(currentStep % passport.timeDiscretizationProperties.saveVP)))
	{
		info('i') << "Saving VP to vtk-file " << std::endl;
		measureVP->Initialization();
	}


	//Вычисление скоростей вихрей: для тех, которые в следе, и виртуальных, а также в точках wakeVP
	CalcVortexVelo(dt, getTimestat().timeCalcVortexConvVelo, getTimestat().timeCalcVortexDiffVelo, getTimestat().timeVP);


	//Расчет и сохранение поля давления
	if ((passport.timeDiscretizationProperties.saveVP > 0) && (!(currentStep % passport.timeDiscretizationProperties.saveVP)))
	{
		measureVP->CalcPressure();
		measureVP->SaveVP(passport.dir, currentStep, getTimestat().timeVP);
	}





	/*	//Сохранятель слоев
		if (parallel.myidWork == 0)
		if (!(currentStep % passport.timeDiscretizationProperties.saveTXT))
		{
		std::string fname = fileNameStep("Bou", W.getPassport().timeDiscretizationProperties.nameLength, currentStep, "txt");

		std::ofstream outfile;
		outfile.open(passport.dir + "snapshots/" + fname);

		outfile << boundary.size() << '\n';
		for(size_t bou = 0; bou < boundary.size(); ++bou)
		outfile << airfoil[bou]->r.size() - 1 << ' ';
		outfile << '\n';

		for (size_t bou = 0; bou < boundary.size(); ++bou)
		for (size_t i = 0; i < airfoil[bou]->r.size() - 1; i++)
		{
		outfile << airfoil[bou]->r[i][0] << ' ' << airfoil[bou]->r[i][1] << ' ' << airfoil[bou]->r[i + 1][0] << ' ' << airfoil[bou]->r[i + 1][1] << ' ' << boundary[bou]->sheets.freeVortexSheet[i][0] << '\n';
		}//for i
		outfile.close();
		}
		*/
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




	//Вычисление сил, действующих на профиль и сохранение в файл	
	if (parallel.myidWork == 0)
	for (size_t mech = 0; mech < mechanics.size(); ++mech)
	{
		mechanics[mech]->GetHydroDynamForce(getTimestat().timeGetHydroDynamForce);

		getTimestat().timeOther.first += omp_get_wtime();
		mechanics[mech]->GenerateForcesString();
		mechanics[mech]->GeneratePositionString();
		getTimestat().timeOther.second += omp_get_wtime();
	}

	std::vector<Point2D> newPos;


//	//////////////
//	//Сохранятель пелены для проверки алгоритма протыкания
//	//до перемещения вихрей
//	if ((parallel.myidWork == 0) && (currentStep >= 0))
//	{
////		if (!(currentStep % 1))
//		{
//			int numVtx = wake->vtx.size();
//			numVtx += boundary[0]->virtualWake.vtx.size();
//
//			std::string fname = VMlib::fileNameStep("InsideOld", getPassport().timeDiscretizationProperties.nameLength, currentStep, "txt");
//			std::ofstream outfile;
//			VMlib::CreateDirectory(passport.dir, "inside");
//			outfile.open(passport.dir + "inside/" + fname);
//
//			outfile << numVtx << std::endl;
//			for (size_t i = 0; i < wake->vtx.size(); i++)
//			{
//				const Point2D& r = wake->vtx[i].r();
//				outfile << static_cast<int>(i) << " " << r[0] << " " << r[1] << " " <<\
//					velocity->wakeVortexesParams.convVelo[i][0] << " " << velocity->wakeVortexesParams.convVelo[i][1] << " " <<\
//					velocity->wakeVortexesParams.diffVelo[i][0] << " " << velocity->wakeVortexesParams.diffVelo[i][1] << std::endl;
//			}//for i	
//			for (size_t i = 0; i < boundary[0]->virtualWake.vtx.size(); i++)
//			{
//				const Point2D& r = boundary[0]->virtualWake.vtx[i].r();
//				outfile << static_cast<int>(i) << " " << r[0] << " " << r[1] << " " << \
//					velocity->virtualVortexesParams[0].convVelo[i][0] << " " << velocity->virtualVortexesParams[0].convVelo[i][1] << " " << \
//					velocity->virtualVortexesParams[0].diffVelo[i][0] << " " << velocity->virtualVortexesParams[0].diffVelo[i][1] << std::endl;
//			}//for i	
//			outfile.close();
//		}
//	}
//	//Сохранятель пелены для проверки алгоритма протыкания
//	//////////////
	
	//Движение вихрей (сброс вихрей + передвижение пелены)
	MoveVortexes(dt, newPos, getTimestat().timeMoveVortexes);


//	//////////////
//	//Сохранятель пелены для проверки алгоритма протыкания
//	//после перемещения вихрей
//	if ((parallel.myidWork == 0) && (currentStep >= 0))
//	{
////		if (!(currentStep % 1))
//		{
//			int numVtx = newPos.size();
//
//			std::string fname = VMlib::fileNameStep("InsideNew", getPassport().timeDiscretizationProperties.nameLength, currentStep, "txt");
//			std::ofstream outfile;
//			VMlib::CreateDirectory(passport.dir, "inside");
//			outfile.open(passport.dir + "inside/" + fname);
//
//			outfile << numVtx << std::endl;
//			for (size_t i = 0; i <newPos.size(); i++)
//			{
//				const Point2D& r = newPos[i];
//				outfile << static_cast<int>(i) << " " << r[0] << " " << r[1] << std::endl;
//			}//for i	
//			outfile.close();
//		}
//	}
//	//Сохранятель пелены для проверки алгоритма протыкания
//	//////////////


	getTimestat().timeOther.first += omp_get_wtime();
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
	{
		switch (passport.numericalSchemes.panelsType.second)
		{
		case 0:
			oldAirfoil.emplace_back(new AirfoilRect(*airfoil[afl]));
			break;
		case 1:

			//airfoil.emplace_back(new AirfoilCurv(GetPassport(), i, parallel));

			info('e') << "AirfoilCurv is not implemented now! " << std::endl;
			exit(1);
			break;
		}

		mechanics[afl]->Move();
	}//for

	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
		boundary[bou]->virtualWake.vtx.clear();
		boundary[bou]->virtualWake.vtx.resize(0);
	}
	getTimestat().timeOther.second += omp_get_wtime();

	CheckInside(newPos, getTimestat().timeCheckInside);

	getTimestat().timeOther.first += omp_get_wtime();


	//передача новых положений вихрей в пелену
	if (parallel.myidWork == 0)
	for (size_t i = 0; i < wake->vtx.size(); ++i)
		wake->vtx[i].r() = newPos[i];
	getTimestat().timeOther.second += omp_get_wtime();

	//Сохранятель пелены ("старый")
	//	if (parallel.myidWork == 0)
	//	if (!(currentStep % passport.timeDiscretizationProperties.saveTXT))
	//		wake.SaveKadr(wake.vtx, passport.dir, currentStep, timestat.timeSaveKadr);




	// Определение параметров, отвечающих за увеличение радиуса коллапса

	std::vector<double> rightBorder, horizSpan;
	rightBorder.reserve(airfoil.size());
	horizSpan.reserve(airfoil.size());

	for (size_t q = 0; q < airfoil.size(); ++q)
	{
		rightBorder.push_back(airfoil[q]->upRight[0]);
		horizSpan.push_back(airfoil[q]->upRight[0] - airfoil[q]->lowLeft[0]);
	}

	if (airfoil.size() > 0)
	{
		wake->collapseRightBorderParameter = *std::max_element(rightBorder.begin(), rightBorder.end());
		wake->collapseScaleParameter = *std::max_element(horizSpan.begin(), horizSpan.end());
	}
	else
	{
		wake->collapseRightBorderParameter = 0.0;
		wake->collapseScaleParameter = 1.0;
	}
#if defined(__CUDACC__) || defined(USE_CUDA)
	cuda.setCollapseCoeff(wake->collapseRightBorderParameter, wake->collapseScaleParameter);
#endif 

	wake->Restruct(getTimestat().timeRestruct);


	if ((parallel.myidWork == 0) && (currentStep >= 0))
	{
		if ((passport.timeDiscretizationProperties.saveTXT > 0) && (!(currentStep % passport.timeDiscretizationProperties.saveTXT)))
		{
			info('i') << "Saving kadr to txt-file " << std::endl;
			wake->SaveKadr(wake->vtx, passport.dir, currentStep, getTimestat().timeSaveKadr);
		}
		if ((passport.timeDiscretizationProperties.saveVTK > 0) && (!(currentStep % passport.timeDiscretizationProperties.saveVTK)))
		{
			info('i') << "Saving kadr to vtk-file " << std::endl;
			wake->SaveKadrVtk(wake->vtx, passport.dir, currentStep, getTimestat().timeSaveKadr);
		}
	}

	getTimestat().timeWakeSort.first = omp_get_wtime();
	
//Сортировка вихрей в пелене по абсциссе; было мнение, что позволяет оптимизировать счет на CUDA
/*
	std::sort(wake.vtx.begin(), wake.vtx.end(), 
		[](const Vortex2D &a, const Vortex2D &b) 
			{
				return a.r()[0] < b.r()[0];
			});
	*/

	getTimestat().timeWakeSort.second = omp_get_wtime();

getTimestat().timeOther.first += omp_get_wtime();
	oldAirfoil.clear();
	//oldAirfoil.resize(0);
getTimestat().timeOther.second += omp_get_wtime();


	//Засечка времени в конце шага
	time.second = omp_get_wtime();


	if (parallel.myidWork == 0)
	{
		info('i') << "Step = " << currentStep \
			<< " PhysTime = " << passport.physicalProperties.getCurrTime() \
			<< " StepTime = " << Times::dT(time) << std::endl;

		getTimestat().GenerateStatString();
	}

	passport.physicalProperties.addCurrTime(dt);
	currentStep++;

}//Step()


//Проверка проникновения вихрей внутрь профиля
void World2D::CheckInside(std::vector<Point2D>& newPos, timePeriod& time)
{
	time.first = omp_get_wtime();
		
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
		wake->Inside(newPos, *airfoil[afl], mechanics[afl]->isMoves, *oldAirfoil[afl]);

	time.second = omp_get_wtime();
}



//Решение системы линейных алгебраических уравнений
void World2D::SolveLinearSystem(timePeriod& time)
{
	time.first = omp_get_wtime();
	//sol = matr.partialPivLu().solve(rhs);


	if (currentStep == 0)
		invMatr = matr.inverse();


/*
	if (currentStep == 0)
	{
		invMatr = matr;
		std::ifstream fileMatrix("invMatrix4x3200.txt");
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

	time.second = omp_get_wtime();
}

//Заполнение матрицы системы для всех профилей
void World2D::FillMatrixAndRhs(timePeriod& time)
{
	time.first = omp_get_wtime();

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
			for (int i = 0; i < nVars; ++i)
			{
				if (currentStep == 0 || mechanics[bou]->isDeform)
				{
					for (int j = 0; j < nVars; ++j)
						matr(i + currentRow, j + currentRow) = locMatr(i, j);
					matr(currentRow + nVars, i + currentRow) = locLastLine(i);
					matr(i + currentRow, currentRow + nVars) = locLastCol(i);
				}
			}

			//\todo пока матрица считается только на первом шаге
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

	/*
	std::ostringstream ss;
	ss << "matrix_";
	ss << currentStep;
	std::ofstream matrFile(ss.str());
	VMlib::SaveToStream(matr, matrFile);
	matrFile.close();
	*/

	
	time.second = omp_get_wtime();
}//FillMatrixAndRhs(...)

//вычисляем размер матрицы и резервируем память под нее и под правую часть
void World2D::ReserveMemoryForMatrixAndRhs(timePeriod& time)
{
	time.first = omp_get_wtime();

	
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

	time.second = omp_get_wtime();
}//ReserveMemoryForMatrixAndRhs(...)


// Вычисляем скорости (и конвективные, и диффузионные) вихрей (в пелене и виртуальных), а также в точках вычисления VP 
void World2D::CalcVortexVelo(double dt, timePeriod& convTime, timePeriod& diffTime, timePeriod& vpTime)
{

	if (parallel.myidWork == 0)
	{
		//Обнуляем все скорости
		convTime.first += omp_get_wtime();
		velocity->wakeVortexesParams.convVelo.clear();
		velocity->wakeVortexesParams.convVelo.resize(wake->vtx.size(), { 0.0, 0.0 });
		convTime.second += omp_get_wtime();

		diffTime.first += omp_get_wtime();
		velocity->wakeVortexesParams.I0.clear();
		velocity->wakeVortexesParams.I0.resize(wake->vtx.size(), 0.0);

		velocity->wakeVortexesParams.I1.clear();
		velocity->wakeVortexesParams.I1.resize(wake->vtx.size(), 0.0);

		velocity->wakeVortexesParams.I2.clear();
		velocity->wakeVortexesParams.I2.resize(wake->vtx.size(), { 0.0, 0.0 });

		velocity->wakeVortexesParams.I3.clear();
		velocity->wakeVortexesParams.I3.resize(wake->vtx.size(), { 0.0, 0.0 });

		velocity->wakeVortexesParams.diffVelo.clear();
		velocity->wakeVortexesParams.diffVelo.resize(wake->vtx.size(), { 0.0, 0.0 });

		velocity->wakeVortexesParams.epsastWake.clear();
		velocity->wakeVortexesParams.epsastWake.resize(wake->vtx.size(), 0.0);
		diffTime.second += omp_get_wtime();

		//Создаем массивы под виртуальные вихри
		convTime.first += omp_get_wtime();
		velocity->virtualVortexesParams.clear();
		velocity->virtualVortexesParams.resize(boundary.size());
		convTime.second += omp_get_wtime();

		for (size_t bou = 0; bou < boundary.size(); ++bou)
		{
			convTime.first += omp_get_wtime();
			velocity->virtualVortexesParams[bou].convVelo.clear();
			velocity->virtualVortexesParams[bou].convVelo.resize(boundary[bou]->virtualWake.vtx.size(), { 0.0, 0.0 });
			convTime.second += omp_get_wtime();

			diffTime.first += omp_get_wtime();
			velocity->virtualVortexesParams[bou].I0.clear();
			velocity->virtualVortexesParams[bou].I0.resize(boundary[bou]->virtualWake.vtx.size(), 0.0);

			velocity->virtualVortexesParams[bou].I1.clear();
			velocity->virtualVortexesParams[bou].I1.resize(boundary[bou]->virtualWake.vtx.size(), 0.0);

			velocity->virtualVortexesParams[bou].I2.clear();
			velocity->virtualVortexesParams[bou].I2.resize(boundary[bou]->virtualWake.vtx.size(), { 0.0, 0.0 });

			velocity->virtualVortexesParams[bou].I3.clear();
			velocity->virtualVortexesParams[bou].I3.resize(boundary[bou]->virtualWake.vtx.size(), { 0.0, 0.0 });

			velocity->virtualVortexesParams[bou].diffVelo.clear();
			velocity->virtualVortexesParams[bou].diffVelo.resize(boundary[bou]->virtualWake.vtx.size(), { 0.0, 0.0 });

			velocity->virtualVortexesParams[bou].epsastWake.clear();
			velocity->virtualVortexesParams[bou].epsastWake.resize(boundary[bou]->virtualWake.vtx.size(), 0.0);
			diffTime.second += omp_get_wtime();
		}
	}


	/*
	info('t') << "st=1: "
		      << boundary[0]->virtualWake.vtx[0].r() << " "
		      << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
		      << velocity->virtualVortexesParams[0].convVelo[0] << " "  
		      << getPassport().physicalProperties.V0() << std::endl;
			  */
	

	convTime.first += omp_get_wtime();

	//Конвективные скорости всех вихрей (и в пелене, и виртуальных), индуцируемые вихрями в пелене
	//Подготовка CUDA
#if defined(__CUDACC__) || defined(USE_CUDA)
	cuda.RefreshWake();
	cuda.RefreshAfls();	
	cuda.RefreshVirtualWakes();
#endif

	convTime.second += omp_get_wtime();
	
	timePeriod convWakeTime, convVPTime;	

	//Вычисление конвективных скоростей вихрей (и в пелене, и виртуальных), а также в точках wakeVP
	velocity->CalcConvVelo(convWakeTime, convVPTime);

	convTime.first += convWakeTime.first;
	convTime.second += convWakeTime.second;

	vpTime.first += convVPTime.first;
	vpTime.second += convVPTime.second;

	convTime.first += omp_get_wtime();

	//std::cout << "Conv Velo time = " << tConvVeloFinish - tConvVeloStart << std::endl;

	/*
	info('t') << "st=2: "
		      << boundary[0]->virtualWake.vtx[0].r() << " "
		      << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
		      << velocity->virtualVortexesParams[0].convVelo[0] << " "  
		      << getPassport().physicalProperties.V0() << std::endl;
			  */
	

	//std::ostringstream sss;
	//sss << "epsast";
	//std::ofstream epsastFile(sss.str());
	//for (size_t i = 0; i < velocity->wakeVortexesParams.epsastWake.size(); ++i)
	//	epsastFile << velocity->wakeVortexesParams.epsastWake[i] << std::endl;
	//epsastFile.close();
	

	//POLARA (убрали вычисление конвективных скоростей по закону Био-Савара от слоев)
	//for (size_t bou = 0; bou < boundary.size(); ++bou)
	//
	//	boundary[bou]->GetConvVelocityToSetOfPoints(wake.vtx, velocity->wakeVortexesParams.convVelo);
	//		
	//	for (size_t targetBou = 0; targetBou < boundary.size(); ++targetBou)
	//		boundary[bou]->GetConvVelocityToSetOfPoints(boundary[targetBou]->virtualWake, velocity->virtualVortexesParams[targetBou].convVelo);
	//}

	convTime.second += omp_get_wtime();

	/*
	info('t') << "st=3: "
		      << boundary[0]->virtualWake.vtx[0].r() << " "
		      << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
		      << velocity->virtualVortexesParams[0].convVelo[0] << " "  
		      << getPassport().physicalProperties.V0() << std::endl;
			  */
	

	//std::ostringstream sss;
	//sss << "epsastVirtual";
	//std::ofstream epsastFile(sss.str());
	//for (size_t i = 0; i < velocity->virtualVortexesParams[0].epsastWake.size(); ++i)
	//	epsastFile << velocity->virtualVortexesParams[0].epsastWake[i] << std::endl;
	//epsastFile.close();

	//exit(1);

	//std::ostringstream sss;
	//sss << "bouVirtVelo_";
	//std::ofstream bouVirtVeloFile(sss.str());
	//for (size_t i = 0; i < velocity->virtualVortexesParams[0].convVelo.size(); ++i)
	//	bouVirtVeloFile << velocity->virtualVortexesParams[0].convVelo[i][0] << " " << velocity->virtualVortexesParams[0].convVelo[i][1] << std::endl;
	//bouVirtVeloFile.close();

	//exit(1);

	//std::ostringstream ss;
	//ss << "bouVelo_";
	//std::ofstream bouVeloFile(ss.str());
	//bouVeloFile << velocity->wakeVortexesParams.convVelo << std::endl;
	//bouVeloFile.close();


	diffTime.first += omp_get_wtime();

	//Вычисление диффузионных скоростей вихрей в следе и виртуальных вихрей, обусловленных влиянием вихрей (тех, что в следе, и виртуальных вихрей)
	
	if (getPassport().physicalProperties.nu > 0.0)
	{

		velocity->CalcDiffVelo();

		/*
		info('t') << "st=4: "
		<< boundary[0]->virtualWake.vtx[0].r() << " "
		<< velocity->virtualVortexesParams[0].diffVelo[0] << " "
		<< velocity->virtualVortexesParams[0].convVelo[0] << " "
		<< getPassport().physicalProperties.V0() << std::endl;
		*/

		//Вычисление диффузионных скоростей вихрей в следе и виртуальных вихрей, обусловленных влиянием поверхностей


		for (size_t afl = 0; afl < airfoil.size(); ++afl)
		{
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I0I3))		
			airfoil[afl]->GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(*wake, velocity->wakeVortexesParams.epsastWake, velocity->wakeVortexesParams.I0, velocity->wakeVortexesParams.I3);
#else
			airfoil[afl]->GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(*wake, velocity->wakeVortexesParams.epsastWake, velocity->wakeVortexesParams.I0, velocity->wakeVortexesParams.I3);
#endif

			for (size_t bou = 0; bou < boundary.size(); ++bou)
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I0I3))			
				airfoil[afl]->GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(boundary[bou]->virtualWake, velocity->virtualVortexesParams[bou].epsastWake, velocity->virtualVortexesParams[bou].I0, velocity->virtualVortexesParams[bou].I3);
#else
				airfoil[afl]->GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(boundary[bou]->virtualWake, velocity->virtualVortexesParams[bou].epsastWake, velocity->virtualVortexesParams[bou].I0, velocity->virtualVortexesParams[bou].I3);
#endif			
		}

		//влияние поверхности
		Point2D I3;
		double I0;


		for (size_t vt = 0; vt < velocity->wakeVortexesParams.diffVelo.size(); ++vt)
		{
			velocity->wakeVortexesParams.I0[vt] *= velocity->wakeVortexesParams.epsastWake[vt];
			velocity->wakeVortexesParams.I0[vt] += DPI * sqr(velocity->wakeVortexesParams.epsastWake[vt]);

			I3 = velocity->wakeVortexesParams.I3[vt];
			I0 = velocity->wakeVortexesParams.I0[vt];

			if (fabs(I0) > 1.e-8)
				velocity->wakeVortexesParams.diffVelo[vt] += I3 * (1.0 / I0);
		}


		for (size_t targetBou = 0; targetBou < boundary.size(); ++targetBou)
		for (size_t vt = 0; vt < velocity->virtualVortexesParams[targetBou].diffVelo.size(); ++vt)
		{
			velocity->virtualVortexesParams[targetBou].I0[vt] *= velocity->virtualVortexesParams[targetBou].epsastWake[vt];
			velocity->virtualVortexesParams[targetBou].I0[vt] += DPI * sqr(velocity->virtualVortexesParams[targetBou].epsastWake[vt]);

			I3 = velocity->virtualVortexesParams[targetBou].I3[vt];
			I0 = velocity->virtualVortexesParams[targetBou].I0[vt];

			if (fabs(I0) > 1.e-8)
				velocity->virtualVortexesParams[targetBou].diffVelo[vt] += I3 * (1.0 / I0);
		}

		diffTime.second += omp_get_wtime();

		/*
		info('t') << "st=5: "
		<< boundary[0]->virtualWake.vtx[0].r() << " "
		<< velocity->virtualVortexesParams[0].diffVelo[0] << " "
		<< velocity->virtualVortexesParams[0].convVelo[0] << " "
		<< GetPassport().physicalProperties.V0() << std::endl;
		*/



		diffTime.first += omp_get_wtime();

		//контроль застрелов диффузионной скорости
		if (parallel.myidWork == 0)
		for (size_t i = 0; i < velocity->wakeVortexesParams.diffVelo.size(); ++i)
		{
			Point2D& diffV = velocity->wakeVortexesParams.diffVelo[i];

			diffV *= passport.physicalProperties.nu;

			if (diffV.length() > 1.5*passport.physicalProperties.vRef)
				diffV.normalize(1.5*passport.physicalProperties.vRef);
		}

		/*
		info('t') << "st=6: "
		<< boundary[0]->virtualWake.vtx[0].r() << " "
		<< velocity->virtualVortexesParams[0].diffVelo[0] << " "
		<< velocity->virtualVortexesParams[0].convVelo[0] << " "
		<< GetPassport().physicalProperties.V0() << std::endl;
		*/

		//контроль "застрелов" диффузионной скорости
		if (parallel.myidWork == 0)
		for (size_t bou = 0; bou < velocity->virtualVortexesParams.size(); ++bou)
		for (size_t i = 0; i < velocity->virtualVortexesParams[bou].diffVelo.size(); ++i)
		{
			Point2D& diffV = velocity->virtualVortexesParams[bou].diffVelo[i];

			diffV *= passport.physicalProperties.nu;
			if (diffV.length() > 1.5*passport.physicalProperties.vRef)
				diffV.normalize(1.5*passport.physicalProperties.vRef);
		}


		if (airfoil.size() > 0)
		{
			for (size_t afl = 0; afl < airfoil.size(); ++afl)
			for (size_t i = 0; i < airfoil[afl]->viscousStress.size(); ++i)
				airfoil[afl]->viscousStress[i] *= passport.physicalProperties.nu;
		}

//		if((passport.timeDiscretizationProperties.saveVTK > 0) &&  !(currentStep % passport.timeDiscretizationProperties.saveVTK))
		if((passport.timeDiscretizationProperties.saveVTK > 0) &&  !(currentStep % 10))
		{
			std::string fname = VMlib::fileNameStep("VisStress", getPassport().timeDiscretizationProperties.nameLength, currentStep, "txt");
			std::ofstream outfile;
			VMlib::CreateDirectory(getPassport().dir, "visStress");
			outfile.open(getPassport().dir + "visStress/" + fname);

			outfile << airfoil[0]->viscousStress.size() << std::endl; //Сохранение числа вихрей в пелене

			for (size_t i = 0; i < airfoil[0]->viscousStress.size(); i++)
			{
				const Point2D& r = 0.5 * (airfoil[0]->getR(i + 1) + airfoil[0]->getR(i));
				double gi = airfoil[0]->viscousStress[i];
				outfile << static_cast<int>(i) << " " << r[0] << " " << r[1] << " " << gi << std::endl;
			}//for i	
			outfile.close();
		}
	}
	diffTime.second += omp_get_wtime();


	//std::ostringstream sss;
	//sss << "velo_";
	//std::ofstream veloFile(sss.str());
	//for (int i = 0; i < velocity->diffVirtualVelo[0].size(); ++i)
	//	veloFile << velocity->diffVirtualVelo[0][i] << std::endl;
	//veloFile.close();

 /*   info('t') << "st=7: "
	          << boundary[0]->virtualWake.vtx[0].r() << " "
	          << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
	          << velocity->virtualVortexesParams[0].convVelo[0] << " "  
	          << getPassport().physicalProperties.V0() << std::endl;
	*/

	//std::ostringstream sss;
	//sss << "velo_";
	//std::ofstream veloFile(sss.str());
	//for (int i = 0; i < velocity->diffVirtualVelo[0].size(); ++i)
	//	veloFile << velocity->diffVirtualVelo[0][i] << std::endl;
	//veloFile.close();


#pragma warning (push)
#pragma warning (disable: 4010)
	//info('t') << "Times: " << \
		tTemp[1] - tTemp[0] << " " << \
		tTemp[2] - tTemp[1] << " " << \
		tTemp[3] - tTemp[2] << " " << \
		tTemp[4] - tTemp[3] << " " << \
		tTemp[5] - tTemp[4] << " " << \
		tTemp[6] - tTemp[5] << " " << \
		tTemp[7] - tTemp[6] << " " << \
		tTemp[8] - tTemp[7] << " " << \
		tTemp[9] - tTemp[8] << " " << \
		std::endl;
#pragma warning (pop)

/*
	std::ofstream outfileEps;
	outfileEps.open(passport.dir + "/veloWake.txt");
	for (int i = 0; i < velocity->wakeVortexesParams.convVelo.size(); i++)
		outfileEps << velocity->wakeVortexesParams.convVelo[i][0] << " " << velocity->wakeVortexesParams.convVelo[i][1] << std::endl;
	outfileEps.close();

	std::ofstream outfileEpsVirt;
	outfileEpsVirt.open(passport.dir + "/veloVirtWake.txt");
	for (int j = 0; j < velocity->virtualVortexesParams.size(); j++)
		for (int i = 0; i < velocity->virtualVortexesParams[j].convVelo.size(); i++)
			outfileEpsVirt << velocity->virtualVortexesParams[j].convVelo[i][0] << " " << velocity->virtualVortexesParams[j].convVelo[i][1] << std::endl;
	outfileEpsVirt.close();
*/

}//CalcVortexVelo(...)

//Вычисляем новые положения вихрей (в пелене и виртуальных)
void World2D::MoveVortexes(double dt, std::vector<Point2D>& newPos, timePeriod& time)
{
	

	size_t nvt = wake->vtx.size();
	for (size_t i = 0; i < getNumberOfBoundary(); ++i)
		nvt += boundary[i]->virtualWake.vtx.size();

	newPos.clear();
	newPos.reserve(nvt);

	time.first = omp_get_wtime();

	if (parallel.myidWork == 0)
		for (size_t i = 0; i < wake->vtx.size(); ++i)
		{
			newPos.push_back(wake->vtx[i].r() \
				+ (velocity->wakeVortexesParams.convVelo[i] + 
				   velocity->wakeVortexesParams.diffVelo[i] + 
				   passport.physicalProperties.V0())*passport.timeDiscretizationProperties.dt);
			//	info('t') << *(newPos.end() - 1) << std::endl;
		}

	if (parallel.myidWork == 0)
		for (size_t bou = 0; bou < boundary.size(); ++bou)
			for (size_t i = 0; i < boundary[bou]->virtualWake.vtx.size(); ++i)
			{
				wake->vtx.push_back(boundary[bou]->virtualWake.vtx[i]);
				newPos.push_back(boundary[bou]->virtualWake.vtx[i].r() \
					+ (velocity->virtualVortexesParams[bou].diffVelo[i] +
					   velocity->virtualVortexesParams[bou].convVelo[i] +
					   passport.physicalProperties.V0())*passport.timeDiscretizationProperties.dt);
			}

	time.second = omp_get_wtime();
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
			//\todo пока матрица считается только на первом шаге
			if (currentStep == 0)
			{
				for (size_t oth = 0; oth < boundary.size(); ++oth)
				{
					size_t nVarsOther = boundary[oth]->GetUnknownsSize();

					if (bou != oth)
					{
						boundary[bou]->FillIQFromOther(*boundary[oth], IQ[bou][oth]);

					}// if (bou != oth)
				}// for oth
			}// if (currentStep == 0 || mechanics[oth]->isMoves)
		}// if (parallel.myidWork == 0)
	}// for bou
}//FillIQ()