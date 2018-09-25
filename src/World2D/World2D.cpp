/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
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
\version 1.1
\date 2 апреля 2018 г.
*/

#include "World2D.h"
#include "Queue.h"


//Конструктор
World2D::World2D(const Queue& _queue, std::ostream& _telefile) :
	teleFile(_telefile),
	queue(_queue),
	parallel(_queue.parallel),
	wake(GetPassport().wakeDiscretizationProperties, parallel, cuda, airfoil, boundary),
	timestat(GetPassport()),
	cuda(gpu(boundary, wake, _queue.parallel))
{
	//myidWork;
	
	GetPassport().physicalProperties.setCurrTime(GetPassport().timeDiscretizationProperties.timeStart);
	currentStep = 0;

	// загрузка пелены из файла
	if (GetPassport().wakeDiscretizationProperties.fileWake != "")
		wake.ReadFromFile(GetPassport().wakesDir); //Считываем из каталога с пеленой

	//airfoil.resize(GetPassport().airfoilParams.size());
	//boundary.resize(GetPassport().airfoilParams.size());
	//mechanics.resize(GetPassport().airfoilParams.size());

	switch (GetPassport().numericalSchemes.velocityComputation)
	{
	case 0:
		velocity.reset(new VelocityBiotSavart(parallel, cuda, wake, boundary));
		break;
	case 1:
		/*
		velocity.reset(new VelocityFourier(parallel, wake, boundary));
		*/	
		teleFile << "VelocityFourier is not implemented now! " << std::endl;
		exit(1);
		break;
	}


	velocity->virtualVortexesParams.resize(GetPassport().airfoilParams.size());

	for (size_t i = 0; i < GetPassport().airfoilParams.size(); ++i)
	{
		switch (GetPassport().airfoilParams[i].panelsType)
		{
		case 0:			
				airfoil.emplace_back(new AirfoilRect(GetPassport(), i, parallel));
			break;
		case 1:
			/*
			airfoil.emplace_back(new AirfoilCurv(GetPassport(), i, parallel));
			*/
			teleFile << "AirfoilCurv is not implemented now! " << std::endl;
			exit(1);
			break;
		}

		airfoil[i]->ReadFromFile(GetPassport().airfoilsDir);	//Считываем из каталога с коллекцией профилей

		switch (GetPassport().airfoilParams[i].boundaryCondition)
		{
		case 0:
			/*
			boundary.emplace_back(new BoundaryMDV(GetPassport(), *(airfoil[i]), boundary, wake, parallel, cuda));
			*/
			teleFile << "BoundaryMDV is not implemented now! " << std::endl;
			exit(1);
			break;
			
		case 1:
			/*
			boundary.emplace_back(new BoundaryVortColl(GetPassport(), *(airfoil[i]), boundary, wake, parallel, cuda));
			*/
			teleFile << "BoundaryVortColl is not implemented now! " << std::endl;
			exit(1);
			break;
			
		case 2:
			boundary.emplace_back(new BoundaryConstLayerAver(GetPassport(), *(airfoil[i]), boundary, wake, parallel, cuda));
			break;
		}


		switch (GetPassport().airfoilParams[i].mechanicalSystem)
		{
		case 0:
			mechanics.emplace_back(new MechanicsRigidImmovable(GetPassport(), *(airfoil[i]), *(boundary[i]), velocity->virtualVortexesParams[i], parallel));
			break;
		case 1:
			mechanics.emplace_back(new MechanicsRigidGivenLaw (GetPassport(), *(airfoil[i]), *(boundary[i]), velocity->virtualVortexesParams[i], parallel));
			break;
		case 2:
			mechanics.emplace_back(new MechanicsRigidOscillPart(GetPassport(), *(airfoil[i]), *(boundary[i]), velocity->virtualVortexesParams[i], parallel));
			break;
		}
	}



	teleFile << "time = " << GetPassport().physicalProperties.getCurrTime() << std::endl;
}//World2D(...)


//Функция, возвращающая константную ссылку на паспорт конкретного расчета
const Passport& World2D::GetPassport() const
{
	return queue.task[problemNumber()].passport;
}//GetPassport()


//Функция, возвращающая номер текущей решаемой задачи(для справки)
size_t World2D::problemNumber() const //Номер текущей решаемой задачи (для справки)
{
	return queue.myProcState;
}//problemNumber()


//Функция, возвращающая признак завершения счета в решаемой задаче
bool World2D::isFinished() const  //Проверка условия, что счет завершен
{
	return (GetPassport().physicalProperties.getCurrTime() >= GetPassport().timeDiscretizationProperties.timeStop);
}//isFinished()



//Основная функция выполнения одного шага по времени
void World2D::Step(timePeriod& time)
{ 	
	//Очистка статистики
	timestat.ToZero();
	
	//Засечка времени в начале шага
	time.first = omp_get_wtime();
	
	const double& dt = GetPassport().timeDiscretizationProperties.dt;	

	//вычисляем скорости панелей
	for (size_t i = 0; i < airfoil.size(); ++i)
		mechanics[i]->VeloOfAirfoilPanels(currentStep * dt);

	//вычисляем интенсивности присоединенных слоев
	for (size_t i = 0; i < airfoil.size(); ++i)
		boundary[i]->ComputeAttachedSheetsIntensity();

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
	}

	wake.WakeSynchronize();

	for (int bou = 0; bou < boundary.size(); ++bou)
	{
		boundary[bou]->VirtualWakeSynchronize();
	}
	
	//Вычисление скоростей вихрей: для тех, которые в следе, и виртуальных
	CalcVortexVelo(dt, timestat.timeCalcVortexConvVelo, timestat.timeCalcVortexDiffVelo);

	//Вычисление сил, действующих на профиль и сохранение в файл	
	if (parallel.myidWork == 0)
	for (size_t mech = 0; mech < mechanics.size(); ++mech)
	{
		mechanics[mech]->GetHydroDynamForce(timestat.timeGetHydroDynamForce);
		mechanics[mech]->GenerateForcesString(currentStep, mech);
		mechanics[mech]->GeneratePositionString(currentStep, mech);
	}

	std::vector<Point2D> newPos;

	//Движение вихрей (сброс вихрей + передвижение пелены)
	MoveVortexes(dt, newPos, timestat.timeMoveVortexes);

	/*//Сохранение пелены после сдвига, но до коллапса и протыкания - для отладки
	std::vector<Vortex2D> vtx2(wake.vtx.size());
	//передача новых положений вихрей в пелену
	if (parallel.myidWork == 0)
	for (size_t i = 0; i < wake.vtx.size(); ++i)
	{
		vtx2[i].r() = newPos[i];
		vtx2[i].g() = wake.vtx[i].g();
	}

	//Переставила сохранятель пелены
	//Сохранятель пелены ("старый")
	if (parallel.myidWork == 0)
		if (!(currentStep % GetPassport().timeDiscretizationProperties.deltacntText))
			wake.SaveKadr(vtx2, GetPassport().dir, currentStep, timestat.timeSaveKadr);
	*/
	
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
	{
		switch (GetPassport().airfoilParams[afl].panelsType)
		{
		case 0:
			oldAirfoil.emplace_back(new AirfoilRect(*airfoil[afl]));
			break;
		case 1:
			
			//airfoil.emplace_back(new AirfoilCurv(GetPassport(), i, parallel));
			
			teleFile << "AirfoilCurv is not implemented now! " << std::endl;
			exit(1);
			break;
		}

		mechanics[afl]->Move();
	}//for

	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
		boundary[bou]->virtualWake.clear();
		boundary[bou]->virtualWake.resize(0);
	}

	CheckInside(newPos, timestat.timeCheckInside);

	//передача новых положений вихрей в пелену
	if (parallel.myidWork == 0)
	for (size_t i = 0; i < wake.vtx.size(); ++i)
		wake.vtx[i].r() = newPos[i];

	//Сохранятель пелены ("старый")
	if (parallel.myidWork == 0)
	if (!(currentStep % GetPassport().timeDiscretizationProperties.deltacntText))
		wake.SaveKadr(wake.vtx, GetPassport().dir, currentStep, timestat.timeSaveKadr);

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


	//Засечка времени в конце шага
	time.second = omp_get_wtime();

	if (parallel.myidWork == 0)
	{
		std::cout << "Step = " << currentStep \
			<< " PhysTime = " << GetPassport().physicalProperties.getCurrTime() \
			<< " StepTime = " << Times::dT(time) << std::endl;

		timestat.GenerateStatString(currentStep, wake.vtx.size());
	}
	
	GetPassport().physicalProperties.addCurrTime(dt);
	currentStep++;
	oldAirfoil.clear();
	//oldAirfoil.resize(0);

}//Step()


//Проверка проникновения вихрей внутрь профиля
void World2D::CheckInside(std::vector<Point2D>& newPos, timePeriod& time)
{
	time.first = omp_get_wtime();
		
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
	{	
		if (mechanics[afl]->isMoves)
			wake.InsideMovingBoundary(newPos, *oldAirfoil[afl], *airfoil[afl]);
		else
			wake.Inside(newPos, *airfoil[afl]);		
	}

	time.second = omp_get_wtime();
}

//Проверка проникновения вихрей внутрь профиля
void World2D::CheckInsideMovingBoundary(std::vector<Point2D>& newPos, timePeriod& time)
{
	time.first = omp_get_wtime();

	for (size_t afl = 0; afl < airfoil.size(); ++afl)
	{
		wake.InsideMovingBoundary(newPos, *oldAirfoil[afl] , * airfoil[afl]);
	}
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
	std::vector <Eigen::MatrixXd> mechRows;	//строки с уравнениями для аэроупругой системы
	std::vector <Eigen::MatrixXd> mechCols;	//столбцы с уравнениями для аэроупругой системы
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
	mechRhs.resize(mechanics.size());

	for (size_t mech = 0; mech < boundary.size(); ++mech)
	{
		mechRows[mech].resize(mechanics[mech]->degOfFreedom, matr.cols());
		for (int i = 0; i < mechRows[mech].rows(); ++i)
		for (int j = 0; j < mechRows[mech].cols(); ++j)
			mechRows[mech](i, j) = 0.0;

		mechCols[mech].resize(matr.rows(), mechanics[mech]->degOfFreedom);
		for (int i = 0; i < mechCols[mech].rows(); ++i)
		for (int j = 0; j < mechCols[mech].cols(); ++j)
			mechCols[mech](i, j) = 0.0;
			
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
		
		
		if (parallel.myidWork == 0)
			boundary[bou]->FillMatrixSelf(locMatr, locLastLine, locLastCol);


		boundary[bou]->FillRhs(GetPassport().physicalProperties.V0(), locRhs, &locLastRhs, mechanics[bou]->isMoves, mechanics[bou]->isDeform);
		for (size_t i = 0; i < airfoil.size(); ++i)
		{
			if (i != bou)
			if (mechanics[i]->isDeform || mechanics[i]->isMoves)
				boundary[bou]->FillRhsFromOther(*airfoil[i], locRhs);
		}

		mechanics[bou]->FillMechanicsRowsAndCols(mechRows[bou], mechCols[bou]);
		mechanics[bou]->FillMechanicsRhs(mechRhs[bou]);

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

			for (int i = 0; i < mechanics[bou]->degOfFreedom; ++i)
			{
				for (int j = 0; j < matr.cols(); ++j)
				{
					matr(currentStep + nVars + 1 + i, j) = mechRows[bou](i, j);
					matr(j, currentStep + nVars + 1 + i) = mechCols[bou](j, i);
				}
				rhs(currentStep + nVars + 1 + i) = mechRhs[bou][i];
			}

			size_t currentCol = 0;
			for (size_t oth = 0; oth < boundary.size(); ++oth)
			{
				size_t nVarsOther = boundary[oth]->GetUnknownsSize();
				if (currentStep == 0 || mechanics[oth]->isMoves)
				{
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
				}// if (currentStep == 0 || mechanics[oth]->isMoves)

				currentCol += nVarsOther + 1 + mechanics[oth]->degOfFreedom;
			}// for oth

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
	//time.first = omp_get_wtime();

	double tTemp[10];

	tTemp[0] = omp_get_wtime();

	if (parallel.myidWork == 0)
	{
		//Обнуляем все скорости
		velocity->wakeVortexesParams.convVelo.clear();
		velocity->wakeVortexesParams.convVelo.resize(wake.vtx.size(), { 0.0, 0.0 });

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

		//Создаем массивы под виртуальные вихри
		//TODO профиль пока строго 1

		velocity->virtualVortexesParams.clear();
		velocity->virtualVortexesParams.resize(airfoil.size());

		for (size_t afl = 0; afl < airfoil.size(); ++afl)
		{
			velocity->virtualVortexesParams[afl].convVelo.clear();
			velocity->virtualVortexesParams[afl].convVelo.resize(airfoil[afl]->np, { 0.0, 0.0 });

			velocity->virtualVortexesParams[afl].I0.clear();
			velocity->virtualVortexesParams[afl].I0.resize(airfoil[afl]->np, 0.0);

			velocity->virtualVortexesParams[afl].I1.clear();
			velocity->virtualVortexesParams[afl].I1.resize(airfoil[afl]->np, 0.0);

			velocity->virtualVortexesParams[afl].I2.clear();
			velocity->virtualVortexesParams[afl].I2.resize(airfoil[afl]->np, { 0.0, 0.0 });

			velocity->virtualVortexesParams[afl].I3.clear();
			velocity->virtualVortexesParams[afl].I3.resize(airfoil[afl]->np, { 0.0, 0.0 });

			velocity->virtualVortexesParams[afl].diffVelo.clear();
			velocity->virtualVortexesParams[afl].diffVelo.resize(airfoil[afl]->np, { 0.0, 0.0 });

			velocity->virtualVortexesParams[afl].epsastWake.clear();
			velocity->virtualVortexesParams[afl].epsastWake.resize(airfoil[afl]->np, 0.0);
		}
	}



	tTemp[1] = omp_get_wtime();

	/*
	std::cout << "st=1: "
		      << boundary[0]->virtualWake[0].r() << " "
		      << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
		      << velocity->virtualVortexesParams[0].convVelo[0] << " "  
		      << GetPassport().physicalProperties.V0() << std::endl;
	*/

	convTime.first = omp_get_wtime();

	//Конвективные скорости всех вихрей (и в пелене, и виртуальных), индуцируемые вихрями в пелене
	//Подготовка CUDA

	//double timeStartRefresh = omp_get_wtime();

#if defined(__CUDACC__) || defined(USE_CUDA)
	cuda.RefreshWake();
	cuda.RefreshAfls();	
#endif

	//double timeEndRefresh = omp_get_wtime();

	//std::cout << "Refresh_time = " << timeEndRefresh - timeStartRefresh << std::endl;

	//Скорости всех вихрей (и в пелене, и виртуальных), индуцируемые вихрями в пелене
	velocity->CalcConvVelo();

	tTemp[2] = omp_get_wtime();

	/*
	std::cout << "st=2: "
		      << boundary[0]->virtualWake[0].r() << " "
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
		cuda.ExpGetConvVelocityToSetOfPointsFromVirtualVortexes(wake.vtx.size(), wake.devWakePtr, boundary[bou]->virtualWake.size(), cuda.host_ptr_ptr_pnl[bou], velocity->wakeVortexesParams.convVelo, cuda.vels, cuda.dev_ptr_vel, wake.param.eps2);
#else
		boundary[bou]->GetConvVelocityToSetOfPointsFromVirtualVortexes(wake.vtx, velocity->wakeVortexesParams.convVelo);		
#endif
		
		for (size_t targetBou = 0; targetBou < boundary.size(); ++targetBou)
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONVVIRT))
			cuda.ExpGetConvVelocityToSetOfPointsFromVirtualVortexes(boundary[targetBou]->virtualWake.size(), cuda.host_ptr_ptr_pnl[targetBou], boundary[bou]->virtualWake.size(), cuda.host_ptr_ptr_pnl[bou], velocity->virtualVortexesParams[targetBou].convVelo, cuda.virtvels[targetBou], cuda.host_ptr_ptr_vel[targetBou], wake.param.eps2);
#else
			boundary[bou]->GetConvVelocityToSetOfPointsFromVirtualVortexes(boundary[targetBou]->virtualWake, velocity->virtualVortexesParams[targetBou].convVelo);
#endif	
	}

	tTemp[3] = omp_get_wtime();

	convTime.second = omp_get_wtime();

	/*
	    std::cout << "st=3: "
		      << boundary[0]->virtualWake[0].r() << " "
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


	diffTime.first = omp_get_wtime();

	//Вычисление диффузионных скоростей вихрей в следе и виртуальных вихрей, обусловленных влиянием вихрей (тех, что в следе, и виртуальных вихрей)
	velocity->CalcDiffVelo();

	
	tTemp[4] = omp_get_wtime();
	
	/*
	std::cout << "st=4: "
	      << boundary[0]->virtualWake[0].r() << " "
	      << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
	      << velocity->virtualVortexesParams[0].convVelo[0] << " "  
	      << GetPassport().physicalProperties.V0() << std::endl;
	*/

	//Вычисление диффузионных скоростей вихрей в следе и виртуальных вихрей, обусловленных влиянием поверхностей
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
	{
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I0I3))
		cuda.ExpGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(wake.vtx.size(), wake.devWakePtr, cuda.dev_ptr_rad, boundary[afl]->afl.r.size(), cuda.host_ptr_ptr_r[afl], velocity->wakeVortexesParams.I0, velocity->wakeVortexesParams.I3, cuda.i0, cuda.i3, cuda.dev_ptr_i0, cuda.dev_ptr_i3);
#else
		airfoil[afl]->GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(wake.vtx, velocity->wakeVortexesParams.epsastWake, velocity->wakeVortexesParams.I0, velocity->wakeVortexesParams.I3);
#endif

		for (size_t bou = 0; bou < boundary.size(); ++bou)
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I0I3))
			cuda.ExpGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(boundary[bou]->virtualWake.size(), cuda.host_ptr_ptr_pnl[bou], cuda.host_ptr_ptr_rad[bou], boundary[afl]->afl.r.size(), cuda.host_ptr_ptr_r[afl], velocity->virtualVortexesParams[bou].I0, velocity->virtualVortexesParams[bou].I3, cuda.virti0[bou], cuda.virti3[bou], cuda.host_ptr_ptr_i0[bou], cuda.host_ptr_ptr_i3[bou]);
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
	
	diffTime.second = omp_get_wtime();
	
	/*
	std::cout << "st=5: "
		  << boundary[0]->virtualWake[0].r() << " "
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

	/// \todo Сделать качественный контроль "застрелов" диффузионной скорости
	if (parallel.myidWork == 0)
	for (size_t i = 0; i < velocity->wakeVortexesParams.diffVelo.size(); ++i)
	{
		Point2D& diffV = velocity->wakeVortexesParams.diffVelo[i];

		diffV *= GetPassport().physicalProperties.nu;


		if (GetPassport().physicalProperties.V0().length2() > 0.0)
		{
			if (diffV.length2() > 0.25*GetPassport().physicalProperties.V0().length2())
			{
				diffV.normalize(0.5*GetPassport().physicalProperties.V0().length());
			}
		}
		/// \todo переделать отночительную скорость
		else
		if (diffV.length2() > 0.25*1.0)
		{
			diffV.normalize(0.5*1.0);
		}
	}
	
	/*
        std::cout << "st=6: "
	          << boundary[0]->virtualWake[0].r() << " "
	          << velocity->virtualVortexesParams[0].diffVelo[0] << " " 
	          << velocity->virtualVortexesParams[0].convVelo[0] << " "  
	          << GetPassport().physicalProperties.V0() << std::endl;
	*/

	if (parallel.myidWork == 0)
	for (size_t bou = 0; bou < velocity->virtualVortexesParams.size();  ++bou)
	for (size_t i = 0; i < velocity->virtualVortexesParams[bou].diffVelo.size(); ++i)
	{
		Point2D& diffV = velocity->virtualVortexesParams[bou].diffVelo[i];

		diffV *= GetPassport().physicalProperties.nu;
		if (GetPassport().physicalProperties.V0().length2() > 0.0)
		{
			if (diffV.length2() > 0.25*GetPassport().physicalProperties.V0().length2())
			{
				diffV.normalize(0.5*GetPassport().physicalProperties.V0().length());
			}
		}
		else
		if (diffV.length2() > 0.25*1.0)
		{
			diffV.normalize(0.5*1.0);
		}
	}

	


	for (size_t i = 0; i < airfoil[0]->viscousStress.size(); ++i)
		airfoil[0]->viscousStress[i] *= GetPassport().physicalProperties.nu;


	tTemp[8] = omp_get_wtime();

	/*
        std::cout << "st=7: "
	          << boundary[0]->virtualWake[0].r() << " "
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

	//std::cout << "Times: " << \
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
		newPos.push_back(wake.vtx[i].r() + (velocity->wakeVortexesParams.convVelo[i] + velocity->wakeVortexesParams.diffVelo[i]  + GetPassport().physicalProperties.V0())*GetPassport().timeDiscretizationProperties.dt);
	//	std::cout << *(newPos.end() - 1) << std::endl;
	}

	if (parallel.myidWork == 0)	
	for (size_t bou = 0; bou < boundary.size(); ++bou)
	for (size_t i = 0; i < boundary[bou]->virtualWake.size(); ++i)
	{
		wake.vtx.push_back(boundary[bou]->virtualWake[i]);
		newPos.push_back(boundary[bou]->virtualWake[i].r() \
		 + ( velocity->virtualVortexesParams[bou].diffVelo[i] + 
		     velocity->virtualVortexesParams[bou].convVelo[i] + 
		     GetPassport().physicalProperties.V0())*GetPassport().timeDiscretizationProperties.dt);
	}
	

	time.second = omp_get_wtime();
}//MoveVortexes(...)

