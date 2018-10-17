/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#include "World2D.h"

//Конструктор
World2D::World2D(const Passport& passport_, const Parallel& parallel_) :
	passport(passport_),
	parallel(parallel_),
	wake(*this),
	timestat(*this),
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

	// загрузка пелены из файла
	if (passport.wakeDiscretizationProperties.fileWake != "")
		wake.ReadFromFile(passport.wakesDir); //Считываем из каталога с пеленой
	
	//airfoil.resize(GetPassport().airfoilParams.size());
	//boundary.resize(GetPassport().airfoilParams.size());
	//mechanics.resize(GetPassport().airfoilParams.size());

	switch (passport.numericalSchemes.velocityComputation)
	{
	case 0:
		velocity.reset(new VelocityBiotSavart(*this));
		break;
	case 1:
		/*
		velocity.reset(new VelocityFourier(*this));
		*/	
		info('e') << "VelocityFourier is not implemented now! " << std::endl;
		exit(1);
		break;
	}

	//считываем массив точек для подсчета и вывода поля скоростей и давлений
	if (passport.timeDiscretizationProperties.saveVP == 0)
		measureVP.reset(nullptr);
	else
	{
		measureVP.reset(new MeasureVP(*this));
		measureVP->ReadPointsFromFile(passport.dir);
	}

	velocity->virtualVortexesParams.resize(passport.airfoilParams.size());

	for (size_t i = 0; i < passport.airfoilParams.size(); ++i)
	{
		switch (passport.airfoilParams[i].panelsType)
		{
		case 0:			
			airfoil.emplace_back(new AirfoilRect(*this, i));
			break;
		case 1:
			/*
			airfoil.emplace_back(new AirfoilCurv(*this, i));
			*/
			info('e') << "AirfoilCurv is not implemented now! " << std::endl;
			exit(1);
			break;
		}

		airfoil[i]->ReadFromFile(passport.airfoilsDir);	//Считываем из каталога с коллекцией профилей
		
		switch (passport.airfoilParams[i].boundaryCondition)
		{
		case 0:
			/*
			boundary.emplace_back(new BoundaryMDV(*this, i));
			*/
			info('e') << "BoundaryMDV is not implemented now! " << std::endl;
			exit(1);
			break;
			
		case 1:
			/*
			boundary.emplace_back(new BoundaryVortColl(*this, i));
			*/
			info('e') << "BoundaryVortColl is not implemented now! " << std::endl;
			exit(1);
			break;
			
		case 2:
			boundary.emplace_back(new BoundaryConstLayerAver(*this, i));
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
			mechanics.emplace_back(new MechanicsRigidOscillMon(*this, i));
			break;
		case 4:
			mechanics.emplace_back(new MechanicsRigidRotateMon(*this, i));
			break;
		}

	}

	info.endl();
	info('i') << "Start solving problem " << passport.problemName << std::endl;
	info.endl();
}//World2D(...)


//Функция, возвращающая признак завершения счета в решаемой задаче
bool World2D::isFinished() const  //Проверка условия, что счет завершен
{
	return (passport.physicalProperties.getCurrTime() >= passport.timeDiscretizationProperties.timeStop);
}//isFinished()



//Основная функция выполнения одного шага по времени
void World2D::Step()
{
	timePeriod& time = timestat.timeWholeStep;

	//Очистка статистики
	timestat.ToZero();
	
	//Засечка времени в начале шага
	time.first = omp_get_wtime();

	const double& dt = passport.timeDiscretizationProperties.dt;


timestat.timeOther.first += omp_get_wtime();
	//вычисляем скорости панелей
	for (size_t i = 0; i < airfoil.size(); ++i)
		mechanics[i]->VeloOfAirfoilPanels(currentStep * dt);
		
	//вычисляем интенсивности присоединенных слоев
	for (size_t i = 0; i < airfoil.size(); ++i)
		boundary[i]->ComputeAttachedSheetsIntensity();
		
timestat.timeOther.second += omp_get_wtime();


	if (airfoil.size() > 0)
	{
		if (parallel.myidWork == 0)
			ReserveMemoryForMatrixAndRhs(timestat.timeReserveMemoryForMatrixAndRhs);
			
#if defined(__CUDACC__) || defined(USE_CUDA)
		cuda.RefreshWake();
		cuda.RefreshAfls();	
#endif

		FillMatrixAndRhs(timestat.timeFillMatrixAndRhs);

		if (parallel.myidWork == 0)
			SolveLinearSystem(timestat.timeSolveLinearSystem);

		/*if (parallel.myidWork == 0)
		{
			for (int i = 0; i < sol.size(); ++i)
			{
				info('t') << i << " " << sol(i) << std::endl;
			}
		}*/

timestat.timeOther.first += omp_get_wtime();
		if (parallel.myidWork == 0)
		{
			size_t currentRow = 0;
			for (size_t bou = 0; bou < boundary.size(); ++bou)
			{
				size_t nVars = boundary[bou]->GetUnknownsSize();
				Eigen::VectorXd locSol, locMechSol;
				locSol.resize(nVars);
				locMechSol.resize(mechanics[bou]->degOfFreedom);
				for (size_t i = 0; i < nVars; ++i)
					locSol(i) = sol(currentRow + i);

				for (size_t i = 0; i < mechanics[bou]->degOfFreedom; ++i)
					locMechSol(i) = sol(currentRow + nVars + 1 + i);

				boundary[bou]->SolutionToFreeVortexSheetAndVirtualVortex(locSol);
				mechanics[bou]->SolutionToMechanicalSystem(locMechSol);
				currentRow += nVars + 1 + mechanics[bou]->degOfFreedom;
			}
			
			/*
			std::ostringstream ss;
			ss << "solution_";
			ss << currentStep;
			std::ofstream solFile(ss.str());
			SaveToStream(sol, solFile);
			solFile.close();
			*/
		}
timestat.timeOther.second += omp_get_wtime();

	}

timestat.timeOther.first += omp_get_wtime();
	wake.WakeSynchronize();

	for (int bou = 0; bou < boundary.size(); ++bou)
	{
		boundary[bou]->VirtualWakeSynchronize();
	}
timestat.timeOther.second += omp_get_wtime();
		
	//Вычисление скоростей вихрей: для тех, которые в следе, и виртуальных
	CalcVortexVelo(dt, timestat.timeCalcVortexConvVelo, timestat.timeCalcVortexDiffVelo);
		
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	//Сохранение пелены до сдвига - для отладки
	size_t wakevirtualsize = wake.vtx.size();
	for (size_t bou = 0; bou < boundary.size(); bou++)
		wakevirtualsize += boundary[bou]->virtualWake.vtx.size();

timestat.timeOther.first += omp_get_wtime();

	std::vector<Vortex2D> vtx2;
	vtx2.reserve(wakevirtualsize);
	if (parallel.myidWork == 0)
	for (size_t i = 0; i < wake.vtx.size(); ++i)
	{
		vtx2.push_back(wake.vtx[i]);
	}
	if (parallel.myidWork == 0)
	for (size_t bou = 0; bou < boundary.size(); bou++)
	for (size_t i = 0; i < boundary[bou]->virtualWake.vtx.size(); ++i)
	{
		vtx2.push_back(boundary[bou]->virtualWake.vtx[i]);
	}
timestat.timeOther.second += omp_get_wtime();

	if (parallel.myidWork == 0)
	{
		if ( (passport.timeDiscretizationProperties.saveTXT>0) && (!(currentStep % passport.timeDiscretizationProperties.saveTXT)) )
			wake.SaveKadr(vtx2, passport.dir, currentStep, timestat.timeSaveKadr);
		if ( (passport.timeDiscretizationProperties.saveVTK>0) && (!(currentStep % passport.timeDiscretizationProperties.saveVTK)) )
			wake.SaveKadrVtk(vtx2, passport.dir, currentStep, timestat.timeSaveKadr);
	}

/*	//Сохранятель слоев
	if (parallel.myidWork == 0)
	if (!(currentStep % passport.timeDiscretizationProperties.saveTXT))
	{
		std::string fname = "Bou";
		if (currentStep < 10) fname += "0";
		if (currentStep < 100) fname += "0";
		if (currentStep < 1000) fname += "0";
		if (currentStep < 10000) fname += "0";

		std::ostringstream ss;
		ss << currentStep;
		fname += ss.str();
		fname += ".txt";

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
		
	//Расчет и сохранение поля давления
	if ( (passport.timeDiscretizationProperties.saveVP>0) && (!(currentStep % passport.timeDiscretizationProperties.saveVP)) )
		measureVP->CalcSaveVP(passport.dir, currentStep, timestat.timeVP);
		
	//Вычисление сил, действующих на профиль и сохранение в файл	
	if (parallel.myidWork == 0)
	for (size_t mech = 0; mech < mechanics.size(); ++mech)
	{
		mechanics[mech]->GetHydroDynamForce(timestat.timeGetHydroDynamForce);
		
timestat.timeOther.first += omp_get_wtime();
		mechanics[mech]->GenerateForcesString();
		mechanics[mech]->GeneratePositionString();
timestat.timeOther.second += omp_get_wtime();
	}

	std::vector<Point2D> newPos;
	
	//Движение вихрей (сброс вихрей + передвижение пелены)
	MoveVortexes(dt, newPos, timestat.timeMoveVortexes);
	

timestat.timeOther.first += omp_get_wtime();
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
	{
		switch (passport.airfoilParams[afl].panelsType)
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
timestat.timeOther.second += omp_get_wtime();

	CheckInside(newPos, timestat.timeCheckInside);
		
timestat.timeOther.first += omp_get_wtime();
	//передача новых положений вихрей в пелену
	if (parallel.myidWork == 0)
	for (size_t i = 0; i < wake.vtx.size(); ++i)
		wake.vtx[i].r() = newPos[i];
timestat.timeOther.second += omp_get_wtime();

	//Сохранятель пелены ("старый")
//	if (parallel.myidWork == 0)
//	if (!(currentStep % passport.timeDiscretizationProperties.saveTXT))
//		wake.SaveKadr(wake.vtx, passport.dir, currentStep, timestat.timeSaveKadr);


	wake.Restruct(timestat.timeRestruct);

	timestat.timeWakeSort.first = omp_get_wtime();
	
//Сортировка вихрей в пелене по абсциссе; было мнение, что позволяет оптимизтровать счет на CUDA
/*
	std::sort(wake.vtx.begin(), wake.vtx.end(), 
		[](const Vortex2D &a, const Vortex2D &b) 
			{
				return a.r()[0] < b.r()[0];
			});
	*/

	timestat.timeWakeSort.second = omp_get_wtime();

	timestat.timeOther.first += omp_get_wtime();
	oldAirfoil.clear();
	//oldAirfoil.resize(0);
	timestat.timeOther.second += omp_get_wtime();


	//Засечка времени в конце шага
	time.second = omp_get_wtime();


	if (parallel.myidWork == 0)
	{
		info('i') << "Step = " << currentStep \
			<< " PhysTime = " << passport.physicalProperties.getCurrTime() \
			<< " StepTime = " << Times::dT(time) << std::endl;

		timestat.GenerateStatString();
	}

	passport.physicalProperties.addCurrTime(dt);
	currentStep++;

}//Step()


//Проверка проникновения вихрей внутрь профиля
void World2D::CheckInside(std::vector<Point2D>& newPos, timePeriod& time)
{
	time.first = omp_get_wtime();
		
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
		wake.Inside(newPos, *airfoil[afl], mechanics[afl]->isMoves, *oldAirfoil[afl]);

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
	Eigen::VectorXd locRhs;
	std::vector <Eigen::MatrixXd> mechRows;	//строки с уравнениями для аэроупругой системы (слагаемые, отвечающие за интенсивность вихревого слоя, в механическом уравнении)
	std::vector <Eigen::MatrixXd> mechCols;	//столбцы с уравнениями для аэроупругой системы (слагаемые, отвечающие за скорость, в граничном условии)
	std::vector <Eigen::MatrixXd> mechCross; //пересечение строк и столбцов, относящхся к механической системе (слагаемое, отвечаюшее за скорость, в механическом уравнении)
	std::vector <std::vector<double>> mechRhs;	//компоненты правой части уравнений для аэроупругой системы


	double locLastRhs = 0.0;

	//обнуляем матрицу на первом шаге расчета
	if (currentStep == 0)
	{
		for (int i = 0; i < matr.rows(); ++i)
		for (int j = 0; j < matr.cols(); ++j)
			matr(i, j) = 0.0;
	}

	mechRows.resize(mechanics.size());
	mechCols.resize(mechanics.size());
	mechCross.resize(mechanics.size());
	mechRhs.resize(mechanics.size());

	for (size_t mech = 0; mech < boundary.size(); ++mech)
	{
		mechRows[mech].resize(mechanics[mech]->degOfFreedom, boundary[mech]->GetUnknownsSize());
		for (int i = 0; i < mechRows[mech].rows(); ++i)
		for (int j = 0; j < mechRows[mech].cols(); ++j)
			mechRows[mech](i, j) = 0.0;

		mechCols[mech].resize(boundary[mech]->GetUnknownsSize(), mechanics[mech]->degOfFreedom);
		for (int i = 0; i < mechCols[mech].rows(); ++i)
		for (int j = 0; j < mechCols[mech].cols(); ++j)
			mechCols[mech](i, j) = 0.0;

		mechCross[mech].resize(mechanics[mech]->degOfFreedom, mechanics[mech]->degOfFreedom);
		for (int i = 0; i < mechCross[mech].rows(); ++i)
		for (int j = 0; j < mechCross[mech].cols(); ++j)
			mechCross[mech](i, j) = 0.0;
			
		mechRhs[mech].resize(mechanics[mech]->degOfFreedom);
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
			locRhs.resize(nVars);
		}
		
		if (currentStep == 0 || mechanics[bou]->isDeform)
		{
			if (parallel.myidWork == 0)
			{
				boundary[bou]->FillMatrixSelf(locMatr, locLastLine, locLastCol);
				mechanics[bou]->FillMechanicsRowsAndCross(mechRows[bou], mechCross[bou]);
			}
		}

		boundary[bou]->FillRhs(passport.physicalProperties.V0(), locRhs, &locLastRhs, mechanics[bou]->isMoves, mechanics[bou]->isDeform);
		mechanics[bou]->FillMechanicsRhs(mechRhs[bou]);

		Eigen::MatrixXd attCols;
		Eigen::MatrixXd attRhs;
		attCols.resize(mechCols[bou].rows(), mechCols[bou].cols());
		attRhs.resize(mechCols[bou].rows(), 1);


		for (int i = 0; i < attCols.rows(); ++i)
		{
			for (int j = 0; j < attCols.cols(); ++j)
				attCols(i, j) = 0.0;

			attRhs(i) = 0.0;
		}

		mechanics[bou]->FillAtt(attCols, attRhs);

	    /// \todo Влияние присоединенных слоев от других профилей пока выключено (для поступательного движения вроде не нужно)
/*		for (size_t i = 0; i < airfoil.size(); ++i)
		{
			if (i != bou)
			if (mechanics[i]->isDeform || mechanics[i]->isMoves)
				boundary[bou]->FillRhsFromOther(*airfoil[i], locRhs);
		}
*/

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

				rhs(i + currentRow) = locRhs(i);
			}
			rhs(currentRow + nVars) = locLastRhs;

			//размазываем механическую систему
			for (int i = 0; i < mechanics[bou]->degOfFreedom; ++i)
			{
				//\todo пока матрица считается только на первом шаге
				if (currentStep == 0)
				{
					for (int j = 0; j < nVars; ++j)//строки и столбцы
					{
						matr(currentRow + nVars + 1 + i, currentRow + j) = mechRows[bou](i, j);
						matr(currentRow + j, currentRow + nVars + 1 + i) = mechCols[bou](j, i);
					}

					for (int j = 0; j < mechanics[bou]->degOfFreedom; ++j)//пересечения строк и столбцов
					{
						matr(currentRow + nVars + 1 + i, currentRow + nVars + 1 + j) = mechCross[bou](i, j);
					}
				}
				rhs(currentRow + nVars + 1 + i) = mechRhs[bou][i];
			}


			//добавляем влияние присоединенных слоев
			for (int i = 0; i < nVars; ++i)
			{
				//\todo пока матрица считается только на первом шаге
				if (currentStep == 0)
				{
					for (int j = 0; j < mechanics[bou]->degOfFreedom; ++j)	//строки и столбцы
						matr(currentRow + i, currentRow + nVars + 1 + j) += attCols(i, j);
				}

				rhs(currentRow + i) += attRhs(i);
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
					currentCol += nVarsOther + 1 + mechanics[oth]->degOfFreedom;
				}// for oth
			}// if (currentStep == 0 || mechanics[oth]->isMoves)

			currentRow += nVars + 1 + mechanics[bou]->degOfFreedom;

		}// if (parallel.myidWork == 0)
	}// for bou


	/*
	std::ostringstream ss;
	ss << "matrix_";
	ss << currentStep;
	std::ofstream matrFile(ss.str());
	SaveToStream(matr, matrFile);
	matrFile.close();
	
	
	
	std::ostringstream sss;
	sss << "rhs_";
	sss << currentStep;
	std::ofstream rhsFile(sss.str());
	SaveToStream(rhs, rhsFile);
	rhsFile.close();
	*/
	
	
	time.second = omp_get_wtime();
}//FillMatrixAndRhs(...)

//вычисляем размер матрицы и резервируем память под нее и под правую часть
void World2D::ReserveMemoryForMatrixAndRhs(timePeriod& time)
{
	time.first = omp_get_wtime();

	if (currentStep == 0)
	{
		size_t matrSize = boundary.size();

		for (auto it = boundary.begin(); it != boundary.end(); ++it)
			matrSize += (*it)->GetUnknownsSize();

		//добавляем уравнения (их количество совпадает с количеством степеней свободы для профилей) для решения механической системы
		for (auto it = mechanics.begin(); it != mechanics.end(); ++it)
			matrSize += (*it)->degOfFreedom;

		matr.resize(matrSize, matrSize);
		matr.setZero();

		rhs.resize(matrSize);
	}
	rhs.setZero();

	time.second = omp_get_wtime();
}//ReserveMemoryForMatrixAndRhs(...)


// Вычисляем скорости (и конвективные, и диффузионные) вихрей (в пелене и виртуальных) 
void World2D::CalcVortexVelo(double dt, timePeriod& convTime, timePeriod& diffTime)
{
	double tTemp[10];

	tTemp[0] = omp_get_wtime();

	if (parallel.myidWork == 0)
	{
		//Обнуляем все скорости
		convTime.first += omp_get_wtime();
		velocity->wakeVortexesParams.convVelo.clear();
		velocity->wakeVortexesParams.convVelo.resize(wake.vtx.size(), { 0.0, 0.0 });
		convTime.second += omp_get_wtime();

		diffTime.first += omp_get_wtime();
		velocity->wakeVortexesParams.I0.clear();
		velocity->wakeVortexesParams.I0.resize(wake.vtx.size(), 0.0);

		velocity->wakeVortexesParams.I1.clear();
		velocity->wakeVortexesParams.I1.resize(wake.vtx.size(), 0.0);

		velocity->wakeVortexesParams.I2.clear();
		velocity->wakeVortexesParams.I2.resize(wake.vtx.size(), { 0.0, 0.0 });

		velocity->wakeVortexesParams.I3.clear();
		velocity->wakeVortexesParams.I3.resize(wake.vtx.size(), { 0.0, 0.0 });

		velocity->wakeVortexesParams.diffVelo.clear();
		velocity->wakeVortexesParams.diffVelo.resize(wake.vtx.size(), { 0.0, 0.0 });

		velocity->wakeVortexesParams.epsastWake.clear();
		velocity->wakeVortexesParams.epsastWake.resize(wake.vtx.size(), 0.0);
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



	tTemp[1] = omp_get_wtime();

	/*
	info('t') << "st=1: "
		      << boundary[0]->virtualWake.vtx[0].r() << " "
		      << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
		      << velocity->virtualVortexesParams[0].convVelo[0] << " "  
		      << GetPassport().physicalProperties.V0() << std::endl;
	*/

	convTime.first += omp_get_wtime();

	//Конвективные скорости всех вихрей (и в пелене, и виртуальных), индуцируемые вихрями в пелене
	//Подготовка CUDA

	//double timeStartRefresh = omp_get_wtime();

#if defined(__CUDACC__) || defined(USE_CUDA)
	cuda.RefreshWake();
	cuda.RefreshAfls();	
#endif

	//double timeEndRefresh = omp_get_wtime();

	//info('t') << "Refresh_time = " << timeEndRefresh - timeStartRefresh << std::endl;

	//Скорости всех вихрей (и в пелене, и виртуальных), индуцируемые вихрями в пелене
	velocity->CalcConvVelo();

	tTemp[2] = omp_get_wtime();

	/*
	info('t') << "st=2: "
		      << boundary[0]->virtualWake.vtx[0].r() << " "
		      << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
		      << velocity->virtualVortexesParams[0].convVelo[0] << " "  
		      << GetPassport().physicalProperties.V0() << std::endl;
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

	//POLARA
	//вычисление конвективных скоростей по закону Био-Савара от виртуальных вихрей
	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONVVIRT))
		boundary[bou]->GPUGetConvVelocityToSetOfPointsFromVirtualVortexes(wake, velocity->wakeVortexesParams.convVelo);		
#else
		boundary[bou]->GetConvVelocityToSetOfPointsFromVirtualVortexes(wake, velocity->wakeVortexesParams.convVelo);		
#endif
		
		for (size_t targetBou = 0; targetBou < boundary.size(); ++targetBou)
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONVVIRT))
			boundary[bou]->GPUGetConvVelocityToSetOfPointsFromVirtualVortexes(boundary[targetBou]->virtualWake, velocity->virtualVortexesParams[targetBou].convVelo);
#else
			boundary[bou]->GetConvVelocityToSetOfPointsFromVirtualVortexes(boundary[targetBou]->virtualWake, velocity->virtualVortexesParams[targetBou].convVelo);
#endif	
	}

	tTemp[3] = omp_get_wtime();

	convTime.second += omp_get_wtime();

	/*
	info('t') << "st=3: "
		      << boundary[0]->virtualWake.vtx[0].r() << " "
		      << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
		      << velocity->virtualVortexesParams[0].convVelo[0] << " "  
		      << GetPassport().physicalProperties.V0() << std::endl;
	*/

	//std::ostringstream sss;
	//sss << "bouVirtVelo_";
	//std::ofstream bouVirtVeloFile(sss.str());
	//for (size_t i = 0; i < velocity->virtualVortexesParams[0].convVelo.size(); ++i)
	//	bouVirtVeloFile << velocity->virtualVortexesParams[0].convVelo[i][0] << " " << velocity->virtualVortexesParams[0].convVelo[i][1] << std::endl;
	//bouVirtVeloFile.close();

	//std::ostringstream ss;
	//ss << "bouVelo_";
	//std::ofstream bouVeloFile(ss.str());
	//bouVeloFile << velocity->wakeVortexesParams.convVelo << std::endl;
	//bouVeloFile.close();


	diffTime.first += omp_get_wtime();



	//Вычисление диффузионных скоростей вихрей в следе и виртуальных вихрей, обусловленных влиянием вихрей (тех, что в следе, и виртуальных вихрей)
	velocity->CalcDiffVelo();
		
	tTemp[4] = omp_get_wtime();
	
	/*
	info('t') << "st=4: "
	      << boundary[0]->virtualWake.vtx[0].r() << " "
	      << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
	      << velocity->virtualVortexesParams[0].convVelo[0] << " "  
	      << GetPassport().physicalProperties.V0() << std::endl;
	*/

	//Вычисление диффузионных скоростей вихрей в следе и виртуальных вихрей, обусловленных влиянием поверхностей


	for (size_t afl = 0; afl < airfoil.size(); ++afl)
	{
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I0I3))		
		airfoil[afl]->GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(wake, velocity->wakeVortexesParams.epsastWake, velocity->wakeVortexesParams.I0, velocity->wakeVortexesParams.I3);
#else
		airfoil[afl]->GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(wake, velocity->wakeVortexesParams.epsastWake, velocity->wakeVortexesParams.I0, velocity->wakeVortexesParams.I3);
#endif

		for (size_t bou = 0; bou < boundary.size(); ++bou)
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I0I3))			
			airfoil[afl]->GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(boundary[bou]->virtualWake, velocity->virtualVortexesParams[bou].epsastWake, velocity->virtualVortexesParams[bou].I0, velocity->virtualVortexesParams[bou].I3);
#else
			airfoil[afl]->GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(boundary[bou]->virtualWake, velocity->virtualVortexesParams[bou].epsastWake, velocity->virtualVortexesParams[bou].I0, velocity->virtualVortexesParams[bou].I3);
#endif			
	}


	tTemp[5] = omp_get_wtime();

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

	tTemp[6] = omp_get_wtime();

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

	tTemp[7] = omp_get_wtime();
	
	diffTime.second += omp_get_wtime();
	
	/*
	info('t') << "st=5: "
		  << boundary[0]->virtualWake.vtx[0].r() << " "
		  << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
		  << velocity->virtualVortexesParams[0].convVelo[0] << " "  
		  << GetPassport().physicalProperties.V0() << std::endl;
	*/

	
	//std::ofstream stressFile("Stresses", std::ios::app);
	//stressFile << currentStep << "	" << GetPassport().physicalProperties.getCurrTime() << "	";
	//for (size_t i = 0; i < airfoil[0]->viscousStress.size(); ++i)
	//	stressFile << airfoil[0]->viscousStress[i] << " ";
	//stressFile << "\n";
	//stressFile.close();

	
	diffTime.first += omp_get_wtime();
	
	/// \todo Сделать качественный контроль "застрелов" диффузионной скорости
		
	if (parallel.myidWork == 0)
	for (size_t i = 0; i < velocity->wakeVortexesParams.diffVelo.size(); ++i)
	{
		Point2D& diffV = velocity->wakeVortexesParams.diffVelo[i];

		diffV *= passport.physicalProperties.nu;


		if (passport.physicalProperties.V0().length2() > 0.0)
		{
			if (diffV.length2() > 0.25*passport.physicalProperties.V0().length2())
			{
				diffV.normalize(0.5*passport.physicalProperties.V0().length());
			}
		}
		/// \todo переделать относительную скорость
		else
		if (diffV.length2() > 0.25*1.0)
		{
			diffV.normalize(0.5*1.0);
		}
	}
	
	/*
    info('t') << "st=6: "
	          << boundary[0]->virtualWake.vtx[0].r() << " "
	          << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
	          << velocity->virtualVortexesParams[0].convVelo[0] << " "  
	          << GetPassport().physicalProperties.V0() << std::endl;
	*/

	/// \todo Сделать качественный контроль "застрелов" диффузионной скорости
	if (parallel.myidWork == 0)
	for (size_t bou = 0; bou < velocity->virtualVortexesParams.size();  ++bou)
	for (size_t i = 0; i < velocity->virtualVortexesParams[bou].diffVelo.size(); ++i)
	{
		Point2D& diffV = velocity->virtualVortexesParams[bou].diffVelo[i];

		diffV *= passport.physicalProperties.nu;
		if (passport.physicalProperties.V0().length2() > 0.0)
		{
			if (diffV.length2() > 0.25*passport.physicalProperties.V0().length2())
			{
				diffV.normalize(0.5*passport.physicalProperties.V0().length());
			}
		}
		/// \todo переделать относительную скорость
		else
		if (diffV.length2() > 0.25*1.0)
		{
			diffV.normalize(0.5*1.0);
		}
	}

	
	

	for (size_t i = 0; i < airfoil[0]->viscousStress.size(); ++i)
		airfoil[0]->viscousStress[i] *= passport.physicalProperties.nu;

	diffTime.second += omp_get_wtime();


	tTemp[8] = omp_get_wtime();

	/*
    info('t') << "st=7: "
	          << boundary[0]->virtualWake.vtx[0].r() << " "
	          << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
	          << velocity->virtualVortexesParams[0].convVelo[0] << " "  
	          << GetPassport().physicalProperties.V0() << std::endl;
	*/

	//std::ostringstream sss;
	//sss << "velo_";
	//std::ofstream veloFile(sss.str());
	//for (int i = 0; i < velocity->diffVirtualVelo[0].size(); ++i)
	//	veloFile << velocity->diffVirtualVelo[0][i] << std::endl;
	//veloFile.close();

	tTemp[9] = omp_get_wtime();

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

}//CalcVortexVelo(...)

//Вычисляем новые положения вихрей (в пелене и виртуальных)
void World2D::MoveVortexes(double dt, std::vector<Point2D>& newPos, timePeriod& time)
{
	/////////K_ZH 05.03.18
	newPos.clear();
	newPos.resize(0);


	time.first = omp_get_wtime();

	if (parallel.myidWork == 0)
		for (size_t i = 0; i < wake.vtx.size(); ++i)
		{
			newPos.push_back(wake.vtx[i].r() + (velocity->wakeVortexesParams.convVelo[i] + velocity->wakeVortexesParams.diffVelo[i] + passport.physicalProperties.V0())*passport.timeDiscretizationProperties.dt);
			//	info('t') << *(newPos.end() - 1) << std::endl;
		}

	if (parallel.myidWork == 0)
		for (size_t bou = 0; bou < boundary.size(); ++bou)
			for (size_t i = 0; i < boundary[bou]->virtualWake.vtx.size(); ++i)
			{
				wake.vtx.push_back(boundary[bou]->virtualWake.vtx[i]);
				newPos.push_back(boundary[bou]->virtualWake.vtx[i].r() \
					+ (velocity->virtualVortexesParams[bou].diffVelo[i] +
						velocity->virtualVortexesParams[bou].convVelo[i] +
						passport.physicalProperties.V0())*passport.timeDiscretizationProperties.dt);
			}

	time.second = omp_get_wtime();
}//MoveVortexes(...)