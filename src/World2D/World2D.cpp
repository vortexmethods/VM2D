/*!
\file
\brief Файл кода с описанием класса World2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include "World2D.h"
#include "Queue.h"


//Конструктор
World2D::World2D(const Queue& _queue, std::ostream& _telefile) :
	teleFile(_telefile),
	queue(_queue),
	parallel(_queue.parallel),
	wake(Passport().wakeDiscretizationProperties, parallel, airfoil),
	timestat(Passport())
{
	//myidWork;
	
	Passport().physicalProperties.setCurrTime(Passport().timeDiscretizationProperties.timeStart);
	currentStep = 0;

	// загрузка пелены из файла
	if (Passport().wakeDiscretizationProperties.fileWake != "")
		wake.ReadFromFile(Passport().wakesDir); //Считываем из каталога с пеленой

	//airfoil.resize(Passport().airfoilParams.size());
	//boundary.resize(Passport().airfoilParams.size());
	//mechanics.resize(Passport().airfoilParams.size());

	switch (Passport().numericalSchemes.velocityComputation)
	{
	case 0:
		velocity.reset(new VelocityBiotSavart(parallel, wake, boundary));
		break;
	case 1:
		velocity.reset(new VelocityFourier(parallel, wake, boundary));
		break;
	}


	velocity->virtualVortexesParams.resize(Passport().airfoilParams.size());

	for (size_t i = 0; i < Passport().airfoilParams.size(); ++i)
	{
		switch (Passport().airfoilParams[i].panelsType)
		{
		case 0:			
				airfoil.emplace_back(new AirfoilRect(Passport(), i, parallel));
			break;
		case 1:
			/*
			airfoil.emplace_back(new AirfoilCurv(Passport(), i, parallel));
			*/
			break;
		}

		airfoil[i]->ReadFromFile(Passport().airfoilsDir);	//Считываем из каталога с коллекцией профилей

		switch (Passport().airfoilParams[i].boundaryCondition)
		{
		case 0:
			//boundary.emplace_back(new BoundaryMDV(Passport(), *(airfoil[i]), boundary, wake, parallel));
			break;
		case 1:
			//boundary.emplace_back(new BoundaryVortColl(Passport(), *(airfoil[i]), boundary, wake, parallel));
			break;
		case 2:
			boundary.emplace_back(new BoundaryConstLayerAver(Passport(), *(airfoil[i]), boundary, wake, parallel));
			break;
		}


		switch (Passport().airfoilParams[i].mechanicalSystem)
		{
		case 0:
			mechanics.emplace_back(new MechanicsRigidImmovable(Passport(), *(airfoil[i]), *(boundary[i]), velocity->virtualVortexesParams[i], parallel));
			break;
		}
	}



	teleFile << "time = " << Passport().physicalProperties.getCurrTime() << std::endl;
}//World2D(...)


//Функция, возвращающая константную ссылку на паспорт конкретного расчета
const Passport& World2D::Passport() const
{
	return queue.task[problemNumber()].passport;
}//Passport()


//Функция, возвращающая номер текущей решаемой задачи(для справки)
int World2D::problemNumber() const //Номер текущей решаемой задачи (для справки)
{
	return queue.myProcState;
}//problemNumber()


//Функция, возвращающая признак завершения счета в решаемой задаче
bool World2D::isFinished() const  //Проверка условия, что счет завершен
{
	return (Passport().physicalProperties.getCurrTime() >= Passport().timeDiscretizationProperties.timeStop);
}//isFinished()


//Основная функция выполнения одного шага по времени
void World2D::Step(timePeriod& time)
{ 			
	//Очистка статистики
	timestat.ToZero();
	
	//Засечка времени в начале шага
	time.first = omp_get_wtime();
	
	const double& dt = Passport().timeDiscretizationProperties.dt;	
	
	//std::cout << "N_PROCESSORS = " << parallel.nProcWork << std::endl;


	if (airfoil.size() > 0)
	{
		if (parallel.myidWork == 0)
			ReserveMemoryForMatrixAndRhs(timestat.timeReserveMemoryForMatrixAndRhs);

		FillMatrixAndRhs(timestat.timeFillMatrixAndRhs);

		if (parallel.myidWork == 0)
			SolveLinearSystem(timestat.timeSolveLinearSystem);


		//virtualWake.virtVtx.resize(0);

		if (parallel.myidWork == 0)
		{
			int currentRow = 0;
			for (int bou = 0; bou < boundary.size(); ++bou)
			{
				int nVars = boundary[bou]->GetUnknownsSize();
				Eigen::VectorXd locSol;
				locSol.resize(nVars);
				for (int i = 0; i < nVars; ++i)
					locSol(i) = sol(currentRow + i);

				boundary[bou]->SolutionToFreeVortexSheetAndVirtualVortex(locSol);
				currentRow += nVars;
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

	//POLARA
	//TODO: Пока нет слоев, сразу сбрасываются вихри; 
	/*if (parallel.myidWork == 0)
	for (int bou = 0; bou < boundary.size(); ++bou)
	{
		for (size_t i = 0; i < boundary[bou]->virtualWake.size(); ++i)
			wake.vtx.push_back(boundary[bou]->virtualWake[i]);
	}*/

	//std::ostringstream ss1;
	//ss1 << "wakeFile_";
	//ss1 << currentStep;
	//std::ofstream wakefile(ss1.str());
	//for (int i = 0; i < wake.vtx.size(); i++)
	//	wakefile << wake.vtx[i].r()[0] << " " << wake.vtx[i].r()[1] << std::endl;
	//wakefile.close();

	////POLARA
	////TODO: Пока нет слоев, сразу сбрасываются вихри; 
	//if (parallel.myidWork == 0)
	//for (int bou = 0; bou < boundary.size(); ++bou)
	//{
	//	for (size_t i = 0; i < boundary[bou]->virtualWake.size(); ++i)
	//		boundary[bou]->virtualWake[i].g() = 0.0;
	//}

	wake.WakeSynchronize();
	
	//Переставила сохранятель пелены
	//Сохранятель пелены ("старый")
	if (parallel.myidWork == 0)
	if (!(currentStep % Passport().timeDiscretizationProperties.deltacntText))
		wake.SaveKadr(Passport().dir, currentStep, timestat.timeSaveKadr);

	//Движение вихрей (сброс вихрей + передвижение пелены)
	CalcVortexVelo(dt, timestat.timeCalcVortexVelo);

	//std::ostringstream ss2;
	//ss2 << "convVelo_";
	//std::ofstream convVeloFile(ss2.str());
	//convVeloFile << velocity->wakeVortexesParams.convVelo << std::endl;
	//convVeloFile.close();

	//std::ostringstream ss3;
	//ss3 << "diffVelo_";
	//std::ofstream diffVeloFile(ss3.str());
	//for (int i = 0; i < velocity->wakeVortexesParams.diffVelo.size(); i++)
	//	diffVeloFile << velocity->wakeVortexesParams.diffVelo[i][0] << " " << velocity->wakeVortexesParams.diffVelo[i][1] << std::endl;
	//diffVeloFile.close();

	//Вычисление сил, действующих на профиль и сохранение в файл	
	if (parallel.myidWork == 0)
	for (size_t mech = 0; mech < mechanics.size(); ++mech)
	{
		mechanics[mech]->GetHydroDynamForce(timestat.timeGetHydroDynamForce);
		mechanics[mech]->GenerateForcesString(currentStep, mech);
	}

	//if (parallel.myidWork == 0)
	//{
	//	std::ostringstream ss;
	//	ss << "convVelo";
	//	ss << currentStep;
	//	std::ofstream convVeloFile(ss.str());
	//	SaveToStream(velocity->wakeVortexesParams.convVelo, convVeloFile);
	//	convVeloFile.close();
	//}

	//if (parallel.myidWork == 0)
	//{
	//	std::ostringstream ss;
	//	ss << "diffVelo";
	//	ss << currentStep;
	//	std::ofstream diffVeloFile(ss.str());
	//	SaveToStream(velocity->wakeVortexesParams.diffVelo, diffVeloFile);
	//	diffVeloFile.close();
	//}

	/*std::ostringstream ss1;
	ss1 << "wakeFile";
	ss1 << currentStep;
	std::ofstream wakefile(ss1.str());
	for (int i = 0; i < wake.vtx.size(); i++)
		wakefile << wake.vtx[i].r() << std::endl;
	wakefile.close();*/

	std::vector<Point2D> newPos;

	MoveVortexes(dt, newPos, timestat.timeMoveVortexes);

	/*
	if (parallel.myidWork == 0)
	{
		std::ostringstream ss;
		ss << "newPos_";
		ss << currentStep;
		std::ofstream newPosFile(ss.str());
		SaveToStream(newPos, newPosFile);
		newPosFile.close();
	}
	*/


	for (int bou = 0; bou < boundary.size(); ++bou)
	{
		boundary[bou]->virtualWake.clear();
		boundary[bou]->virtualWake.resize(0);
	}
	
	CheckInside(newPos, timestat.timeCheckInside);

	//передача новых положений вихрей в пелену
	if (parallel.myidWork == 0)
	for (size_t i = 0; i < wake.vtx.size(); ++i)
		wake.vtx[i].r() = newPos[i];

	wake.Restruct(timestat.timeRestruct);
	
	//Сохранятель пелены ("старый")
	//if (parallel.myidWork == 0)
	//if (!(currentStep % Passport().timeDiscretizationProperties.deltacntText))
	//	wake.SaveKadr(Passport().dir, currentStep, timestat.timeSaveKadr);

	//Засечка времени в конце шага
	time.second = omp_get_wtime();

	if (parallel.myidWork == 0)
	{
		std::cout << "Step = " << currentStep \
			<< " PhysTime = " << Passport().physicalProperties.getCurrTime() \
			<< " StepTime = " << Times::dT(time) << std::endl;

		timestat.GenerateStatString(currentStep, wake.vtx.size());
	}
	
	Passport().physicalProperties.addCurrTime(dt);
	currentStep++;
}//Step()


//Проверка проникновения вихрей внутрь профиля
void World2D::CheckInside(std::vector<Point2D>& newPos, timePeriod& time)
{
	time.first = omp_get_wtime();
		
	for (size_t afl = 0; afl < airfoil.size(); ++afl)
	{		
		Airfoil& aflRef = *(wake.airfoils[afl]);
		wake.Inside(newPos, aflRef);
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
	double locLastRhs = 0.0;

	if (currentStep == 0)
	{
		for (int i = 0; i < matr.rows(); ++i)
		for (int j = 0; j < matr.cols(); ++j)
			matr(i, j) = 0.0;
	}

	int totalSolVars = matr.rows() - boundary.size();
	int currentRow = 0;

	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
		int nVars;	
	
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

		if (currentStep == 0)
		if (parallel.myidWork == 0)
			boundary[bou]->FillMatrixSelf(locMatr, locLastLine, locLastCol);

		boundary[bou]->FillRhs(Passport().physicalProperties.V0(), locRhs, &locLastRhs);

		//размазываем матрицу
		if (parallel.myidWork == 0)
		{
			for (int i = 0; i < nVars; ++i)
			{
				if (currentStep == 0)
				{
					for (int j = 0; j < nVars; ++j)
						matr(i + currentRow, j + currentRow) = locMatr(i, j);
					matr(totalSolVars + bou, i + currentRow) = locLastLine(i);
					matr(i + currentRow, totalSolVars + bou) = locLastCol(i);
				}

				rhs(i + currentRow) = locRhs(i);
			}
			rhs(totalSolVars + bou) = locLastRhs;

			if (currentStep == 0)
			{
				int currentCol = 0;
				for (size_t oth = 0; oth < boundary.size(); ++oth)
				{
					int nVarsOther = boundary[oth]->GetUnknownsSize();
					if (bou != oth)
					{
						otherMatr.resize(nVars, nVarsOther);
						boundary[bou]->FillMatrixFromOther(*boundary[oth], otherMatr);

						//размазываем матрицу
						for (int i = 0; i < nVars; ++i)
						{
							for (int j = 0; j < nVarsOther; ++j)
								matr(i + currentRow, j + currentCol) = otherMatr(i, j);
						}
					}

					currentCol += nVarsOther;
				}

				currentRow += nVars;
			}
		}
	}


	/*std::ostringstream ss;
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
	rhsFile.close();*/

	time.second = omp_get_wtime();
}//FillMatrixAndRhs(...)

//вычисляем размер матрицы и резервируем память под нее и под правую часть
void World2D::ReserveMemoryForMatrixAndRhs(timePeriod& time)
{
	time.first = omp_get_wtime();

	if (currentStep == 0)
	{
		int matrSize = boundary.size();

		for (auto it = boundary.begin(); it != boundary.end(); ++it)
			matrSize += (*it)->GetUnknownsSize();

		matr.resize(matrSize, matrSize);
		matr.setZero();

		rhs.resize(matrSize);
	}
	rhs.setZero();

	time.second = omp_get_wtime();
}//ReserveMemoryForMatrixAndRhs(...)


// Вычисляем скорости вихрей (в пелене и виртуальных)
void World2D::CalcVortexVelo(double dt, timePeriod& time)
{
	time.first = omp_get_wtime();

	if (parallel.myidWork == 0)
	{
		//Обнуляем все скорости
		velocity->wakeVortexesParams.convVelo.clear();
		velocity->wakeVortexesParams.convVelo.resize(wake.vtx.size(), { 0.0, 0.0 });

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

			velocity->virtualVortexesParams[afl].diffVelo.clear();
			velocity->virtualVortexesParams[afl].diffVelo.resize(airfoil[afl]->np, { 0.0, 0.0 });

			velocity->virtualVortexesParams[afl].epsastWake.clear();
			velocity->virtualVortexesParams[afl].epsastWake.resize(airfoil[afl]->np, 0.0);
		}
	}


	//Скорости всех вихрей (и в пелене, и виртуальных), индуцируемые вихрями в пелене
	velocity->CalcConvVelo();

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
		boundary[bou]->GetConvVelocityToSetOfPointsFromVirtualVortexes(wake.vtx, velocity->wakeVortexesParams.convVelo);
			
		for (size_t targetBou = 0; targetBou < boundary.size(); ++targetBou)
			boundary[bou]->GetConvVelocityToSetOfPointsFromVirtualVortexes(boundary[targetBou]->virtualWake, velocity->virtualVortexesParams[targetBou].convVelo);
	}

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

	velocity->CalcDiffVelo();

	for (size_t afl = 0; afl < airfoil.size(); ++afl)
	{
		airfoil[afl]->GetDiffVelocityToSetOfPointsAndViscousStresses(wake.vtx, velocity->wakeVortexesParams.epsastWake, velocity->wakeVortexesParams.diffVelo);
		
		for (size_t bou = 0; bou < boundary.size(); ++bou)
			airfoil[afl]->GetDiffVelocityToSetOfPointsAndViscousStresses(boundary[bou]->virtualWake, velocity->virtualVortexesParams[bou].epsastWake, velocity->virtualVortexesParams[bou].diffVelo);
	}

	for (size_t i = 0; i < airfoil[0]->viscousStress.size(); ++i)
		airfoil[0]->viscousStress[i] *= Passport().physicalProperties.nu;

	//std::ofstream stressFile("Stresses", std::ios::app);
	//stressFile << currentStep << "	" << Passport().physicalProperties.getCurrTime() << "	";
	//for (size_t i = 0; i < airfoil[0]->viscousStress.size(); ++i)
	//	stressFile << airfoil[0]->viscousStress[i] << " ";
	//stressFile << "\n";
	//stressFile.close();

	/// \todo Сделать качественный контроль "застрелов" диффузионной скорости
	if (parallel.myidWork == 0)
	for (size_t i = 0; i < velocity->wakeVortexesParams.diffVelo.size(); ++i)
	{
		Point2D& diffV = velocity->wakeVortexesParams.diffVelo[i];

		diffV *= Passport().physicalProperties.nu;
		if (diffV.length2() > 0.25)
		{
			diffV.normalize(0.5);
		}

	}

	if (parallel.myidWork == 0)
	for (size_t bou = 0; bou < velocity->virtualVortexesParams.size();  ++bou)
	for (size_t i = 0; i < velocity->virtualVortexesParams[bou].diffVelo.size(); ++i)
	{
		Point2D& diffV = velocity->virtualVortexesParams[bou].diffVelo[i];

		diffV *= Passport().physicalProperties.nu;
		if (diffV.length2() > 0.25)
		{
			diffV.normalize(0.5);
		}
	}


	//std::ostringstream sss;
	//sss << "velo_";
	//std::ofstream veloFile(sss.str());
	//for (int i = 0; i < velocity->diffVirtualVelo[0].size(); ++i)
	//	veloFile << velocity->diffVirtualVelo[0][i] << std::endl;
	//veloFile.close();

	time.second = omp_get_wtime();
}//CalcVortexVelo(...)

//Вычисляем новые положения вихрей (в пелене и виртуальных)
void World2D::MoveVortexes(double dt, std::vector<Point2D>& newPos, timePeriod& time)
{	
	time.first = omp_get_wtime();

	if (parallel.myidWork == 0)
	for (size_t i = 0; i < wake.vtx.size(); ++i)
	{
		newPos.push_back(wake.vtx[i].r() + (velocity->wakeVortexesParams.convVelo[i] + velocity->wakeVortexesParams.diffVelo[i]  + Passport().physicalProperties.V0())*Passport().timeDiscretizationProperties.dt);
	//	std::cout << *(newPos.end() - 1) << std::endl;
	}

	
	for (size_t bou = 0; bou< boundary.size(); bou++)
	for (size_t i = 0; i < boundary[bou]->virtualWake.size(); ++i)
	{
		wake.vtx.push_back(boundary[bou]->virtualWake[i]);
		newPos.push_back(boundary[bou]->virtualWake[i].r()  \
			+ (velocity->virtualVortexesParams[bou].diffVelo[i] + velocity->virtualVortexesParams[bou].convVelo[i] + Passport().physicalProperties.V0())*Passport().timeDiscretizationProperties.dt);
	}
	

	time.second = omp_get_wtime();
}//MoveVortexes(...)

