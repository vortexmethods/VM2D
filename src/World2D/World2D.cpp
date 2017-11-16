/*!
\file
\brief Файл кода с описанием класса World2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include <ctime>

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
		wake.ReadFromFile("./Wakes/"); //Считываем из каталога с пеленой

	airfoil.resize(Passport().airfoilParams.size());
	boundary.resize(Passport().airfoilParams.size());
	mechanics.resize(Passport().airfoilParams.size());

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

	for (size_t i = 0; i < airfoil.size(); ++i)
	{
		switch (Passport().airfoilParams[i].panelsType)
		{
		case 0:			
				airfoil[i].reset(new AirfoilRect(Passport(), i, parallel));
			break;
		case 1:
			/*
			airfoil[i].reset(new AirfoilCurv(Passport(), i, parallel));
			*/
			break;
		}

		airfoil[i]->ReadFromFile("./Airfoils/");	//Считываем из каталога с коллекцией профилей

		switch (Passport().airfoilParams[i].boundaryCondition)
		{
		case 0:
			boundary[i].reset(new BoundaryMDV(Passport(), *(airfoil[i]), boundary, wake, parallel));
			break;
		case 1:
			boundary[i].reset(new BoundaryVortColl(Passport(), *(airfoil[i]), boundary, wake, parallel));
			break;
		case 2:
			boundary[i].reset(new BoundaryConstLayerAver(Passport(), *(airfoil[i]), boundary, wake, parallel));
			break;
		}


		switch (Passport().airfoilParams[i].mechanicalSystem)
		{
		case 0:
			mechanics[i].reset(new MechanicsRigidImmovable(Passport(), *(airfoil[i]), *(boundary[i]), velocity->virtualVortexesParams[i], parallel));
			break;
		}
	}



	teleFile << "time = " << Passport().physicalProperties.getCurrTime() << std::endl;
}//World2D(...)


//Функция, возвраящающая константную ссылку на паспорт конкретного расчета
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
	time.first = clock();	
	
	const double& dt = Passport().timeDiscretizationProperties.dt;	
	
	//std::cout << "N_PROCESSORS = " << parallel.nProcWork << std::endl;

	if (parallel.myidWork == 0)
	{
		//TODO профиль пока строго 1 или строго 0
		if (airfoil.size() > 0 )
		{
			ReserveMemoryForMatrixAndRhs(timestat.timeReserveMemoryForMatrixAndRhs);

			FillMatrixAndRhs(timestat.timeFillMatrixAndRhs);

			SolveLinearSystem(timestat.timeSolveLinearSystem);

			//virtualWake.virtVtx.resize(0);

			//TODO профиль пока строго 1
			int nVars = boundary[0]->GetUnknownsSize();
			Eigen::VectorXd locSol;
			locSol.resize(nVars);
			for (int i = 0; i < nVars; ++i)
				locSol(i) = sol(i);

			boundary[0]->SolutionToFreeVortexSheetAndVirtualVortex(locSol);	

			std::ostringstream ss;
			ss << "solution_";
			ss << currentStep;
			ss << ".txt";
			std::ofstream solFile(ss.str());
			SaveToStream(sol, solFile);
			solFile.close();
			
		}
	}

	//Движение вихрей (сброс вихрей + передвижение пелены)
	CalcVortexVelo(dt, timestat.timeCalcVortexVelo);

	//Вычисление сил, действующих на профиль и сохранение в файл
	// TODO пока профиль ровно 1
	mechanics[0]->GetHydroDynamForce(timestat.timeGetHydroDynamForce);
	mechanics[0]->GenerateForcesString(currentStep, 0);

	std::vector<Point2D> newPos;

	MoveVortexes(dt, newPos, timestat.timeMoveVortexes);

	//std::ostringstream ss;
	//ss << "newPos_";
	//ss << currentStep;
	//ss << ".txt";
	//std::ofstream newPosFile(ss.str());
	//SaveToStream(newPos, newPosFile);
	//newPosFile.close();

	//TODO профиль 1
	if (airfoil.size() > 0)
	boundary[0]->virtualWake.resize(0);

	CheckInside(newPos, timestat.timeCheckInside);

	//передача новых положений вихрей в пелену
	for (size_t i = 0; i < wake.vtx.size(); ++i)
		wake.vtx[i].r() = newPos[i];

	wake.Restruct(timestat.timeRestruct);
	
	//Сохранятель пелены ("старый")
	if (parallel.myidWork == 0)
	if (!(currentStep % Passport().timeDiscretizationProperties.deltacntText))
		wake.SaveKadr(Passport().dir, currentStep, timestat.timeSaveKadr);

	//Засечка времени в конце шага
	time.second = clock();

	std::cout << "Step = " << currentStep \
		<< " PhysTime = " << Passport().physicalProperties.getCurrTime() \
		<< " StepTime = " << Times::dT(time) << std::endl;  

	timestat.GenerateStatString(currentStep, wake.vtx.size());

	Passport().physicalProperties.addCurrTime(dt);
	currentStep++;
}//Step()


//Проверка проникновения вихрей внутрь профиля
void World2D::CheckInside(std::vector<Point2D>& newPos, timePeriod& time)
{
	time.first = clock();

	//TODO профиль 1
	if (airfoil.size() > 0)
	{
		Airfoil& afl = *(wake.airfoils[0]);
		wake.Inside(newPos, afl);
	}

	time.second = clock();
}


//Решение системы линейных алгебраических уравнений
void World2D::SolveLinearSystem(timePeriod& time)
{
	time.first = clock();
	sol = matr.partialPivLu().solve(rhs);
	time.second = clock();
}

//Заполнение матрицы системы для всех профилей
void World2D::FillMatrixAndRhs(timePeriod& time)
{
	time.first = clock();

	Eigen::MatrixXd locMatr;
	Eigen::VectorXd locLastLine, locLastCol;
	Eigen::VectorXd locRhs;
	double locLastRhs = 0.0;

	int nVars = boundary[0]->GetUnknownsSize();
	locMatr.resize(nVars, nVars);
	locLastLine.resize(nVars);
	locLastCol.resize(nVars);
	locRhs.resize(nVars);

	boundary[0]->FillMatrixSelf(locMatr, locLastLine, locLastCol);
	boundary[0]->FillRhs(Passport().physicalProperties.V0(), locRhs, &locLastRhs);

	//размазываем матрицу
	//TODO: переделать для общего случая
	for (int i = 0; i < nVars; ++i)
	{
		for (int j = 0; j < nVars; ++j)
			matr(i, j) = locMatr(i, j);
		matr(nVars, i) = locLastLine(i);
		matr(i, nVars) = locLastCol(i);

		rhs(i) = locRhs(i);
	}

	matr(nVars, nVars) = 0.0;
	rhs(nVars) = locLastRhs;

	//std::ostringstream ss;
	//ss << "matrix_";
	//ss << currentStep;
	//ss << ".txt";
	//std::ofstream matrFile(ss.str());
	//SaveToStream(matr, matrFile);
	//matrFile.close();

	//std::ostringstream sss;
	//sss << "rhs_";
	//sss << currentStep;
	//sss << ".txt";
	//std::ofstream rhsFile(sss.str());
	//SaveToStream(rhs, rhsFile);
	//rhsFile.close();

	time.second = clock();
}//FillMatrixAndRhs(...)

//вычисляем размер матрицы и резервируем память под нее и под правую часть
void World2D::ReserveMemoryForMatrixAndRhs(timePeriod& time)
{
	time.first = clock();

	int matrSize = boundary.size();

	for (auto it = boundary.begin(); it != boundary.end(); ++it)
		matrSize += (*it)->GetUnknownsSize();

	matr.resize(matrSize, matrSize);
	matr.setZero();

	rhs.resize(matrSize);
	rhs.setZero();

	time.second = clock();
}//ReserveMemoryForMatrixAndRhs(...)


// Вычисляем скорости вихрей (в пелене и виртуальных)
void World2D::CalcVortexVelo(double dt, timePeriod& time)
{
	time.first = clock();

	//Обнуляем все скорости
	velocity->wakeVortexesParams.convVelo.clear();
	velocity->wakeVortexesParams.convVelo.resize(wake.vtx.size(), { 0.0, 0.0 });

	velocity->wakeVortexesParams.diffVelo.clear();
	velocity->wakeVortexesParams.diffVelo.resize(wake.vtx.size(), { 0.0, 0.0 });

	velocity->wakeVortexesParams.epsastWake.clear();
	velocity->wakeVortexesParams.epsastWake.resize(wake.vtx.size(), 0.0);

	//Создаем массивы под виртуальные вихри
	//TODO профиль пока строго 1
	if (airfoil.size() > 0)
	{
		velocity->virtualVortexesParams.clear();
		velocity->virtualVortexesParams.resize(1);

		velocity->virtualVortexesParams[0].convVelo.clear();
		velocity->virtualVortexesParams[0].convVelo.resize(airfoil[0]->np, { 0.0, 0.0 });

		velocity->virtualVortexesParams[0].diffVelo.clear();
		velocity->virtualVortexesParams[0].diffVelo.resize(airfoil[0]->np, { 0.0, 0.0 });

		velocity->virtualVortexesParams[0].epsastWake.clear();
		velocity->virtualVortexesParams[0].epsastWake.resize(airfoil[0]->np, 0.0 );
	}

	//Скорости всех вихрей (и в пелене, и виртуальных), индуцируемые вихрями в пелене
	velocity->CalcConvVelo();

	

	//std::ostringstream sss;
	//sss << "VirtVelo_";
	//sss << ".txt";
	//std::ofstream VirtVeloFile(sss.str());
	//for (size_t i = 0; i < velocity->virtualVortexesParams[0].convVelo.size(); ++i)
	//	VirtVeloFile << velocity->virtualVortexesParams[0].convVelo[i][0] << " " << velocity->virtualVortexesParams[0].convVelo[i][1] << std::endl;
	//VirtVeloFile.close();

	//Влияние профиля
	//TODO профиль пока строго 1 или строго 0
	if (airfoil.size() > 0)
	{
		
		boundary[0]->GetConvVelocityToSetOfPoints(wake.vtx, velocity->wakeVortexesParams.convVelo);
		boundary[0]->GetConvVelocityToSetOfPoints(boundary[0]->virtualWake, velocity->virtualVortexesParams[0].convVelo);

		//std::ostringstream sss;
		//sss << "bouVirtVelo_";
		//sss << ".txt";
		//std::ofstream bouVirtVeloFile(sss.str());
		//for (size_t i = 0; i < velocity->virtualVortexesParams[0].convVelo.size(); ++i)
		//	bouVirtVeloFile << velocity->virtualVortexesParams[0].convVelo[i][0] << " " << velocity->virtualVortexesParams[0].convVelo[i][1] << std::endl;
		//bouVirtVeloFile.close();

		//std::ostringstream ss;
		//ss << "bouVelo_";
		//ss << ".txt";
		//std::ofstream bouVeloFile(ss.str());
		//for (size_t i = 0; i < wakeVelo.size(); ++i)
		//	bouVeloFile << wakeVelo[i] << std::endl;
		//bouVeloFile.close();

		//TODO Профиль ровно 1
		}
	
	velocity->CalcDiffVelo(Passport().physicalProperties.nu);

	/*boundary[0]->GetConvVelocityToSetOfPoints(boundary[0]->virtualWake, velocity->virtualVortexesParams[0].convVelo);
	std::ostringstream sss;
	sss << "bouVirtVelo_";
	sss << ".txt";
	std::ofstream bouVirtVeloFile(sss.str());
	for (size_t i = 0; i < velocity->virtualVortexesParams[0].convVelo.size(); ++i)
		bouVirtVeloFile << velocity->virtualVortexesParams[0].convVelo[i] << std::endl;
	bouVirtVeloFile.close();*/

	if (airfoil.size() > 0)
	{
		airfoil[0]->GetDiffVelocityToSetOfPointsAndViscousStresses(wake.vtx, velocity->wakeVortexesParams.epsastWake, velocity->wakeVortexesParams.diffVelo, Passport().wakeDiscretizationProperties.epscol);
		airfoil[0]->GetDiffVelocityToSetOfPointsAndViscousStresses(boundary[0]->virtualWake, velocity->virtualVortexesParams[0].epsastWake, velocity->virtualVortexesParams[0].diffVelo, Passport().wakeDiscretizationProperties.epscol);
	}

	for (size_t i = 0; i < airfoil[0]->viscousStress.size(); ++i)
		airfoil[0]->viscousStress[i] *= Passport().physicalProperties.nu;

	std::ofstream stressFile("Stresses.txt", std::ios::app);
	stressFile << currentStep << "	" << Passport().physicalProperties.getCurrTime() << "	";
	for (size_t i = 0; i < airfoil[0]->viscousStress.size(); ++i)
		stressFile << airfoil[0]->viscousStress[i] << " ";
	stressFile << "\n";
	stressFile.close();

	for (size_t i = 0; i < velocity->wakeVortexesParams.diffVelo.size(); ++i)
	{
		Point2D& diffV = velocity->wakeVortexesParams.diffVelo[i];

		diffV *= Passport().physicalProperties.nu;
		if (diffV.length2() > 0.25)
		{
			diffV.normalize(0.5);
		}

	}

	for (size_t i = 0; i < velocity->virtualVortexesParams[0].diffVelo.size(); ++i)
	{
		Point2D& diffV = velocity->virtualVortexesParams[0].diffVelo[i];

		diffV *= Passport().physicalProperties.nu;
		if (diffV.length2() > 0.25)
		{
			diffV.normalize(0.5);
		}

	}
	

	//std::ostringstream sss;
	//sss << "velo_";
	//sss << ".txt";
	//std::ofstream veloFile(sss.str());
	//for (int i = 0; i < velocity->diffVirtualVelo[0].size(); ++i)
	//	veloFile << velocity->diffVirtualVelo[0][i] << std::endl;
	//veloFile.close();

	time.second = clock();
}//CalcVortexVelo(...)

//Вычисляем новые положения вихрей (в пелене и виртуальных)
void World2D::MoveVortexes(double dt, std::vector<Point2D>& newPos, timePeriod& time)
{	
	time.first = clock();

	if (parallel.myidWork == 0)
	for (size_t i = 0; i < wake.vtx.size(); ++i)
	{
		newPos.push_back(wake.vtx[i].r() + (velocity->wakeVortexesParams.convVelo[i] + velocity->wakeVortexesParams.diffVelo[i]  + Passport().physicalProperties.V0())*Passport().timeDiscretizationProperties.dt);
	//	std::cout << *(newPos.end() - 1) << std::endl;
	}



	//TODO профиль 1
	if (airfoil.size() > 0)
	for (size_t i = 0; i < boundary[0]->virtualWake.size(); ++i)
	{
		wake.vtx.push_back(boundary[0]->virtualWake[i]);
		newPos.push_back(boundary[0]->virtualWake[i].r()  \
			+ (velocity->virtualVortexesParams[0].diffVelo[i] + velocity->virtualVortexesParams[0].convVelo[i] + Passport().physicalProperties.V0())*Passport().timeDiscretizationProperties.dt);
	}

	time.second = clock();
}//MoveVortexes(...)

