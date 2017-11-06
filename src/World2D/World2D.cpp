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
	wake(Passport().wakeDiscretizationProperties, parallel)
{
	//myidWork;
	
	currentTime = Passport().timeDiscretizationProperties.timeStart;
	currentStep = 0;

	// загрузка пелены из файла
	if (Passport().wakeDiscretizationProperties.fileWake != "")
		wake.ReadFromFile("./Wakes/"); //Считываем из каталога с пеленой

	airfoil.resize(Passport().airfoilParams.size());
	boundary.resize(Passport().airfoilParams.size());

	for (size_t i = 0; i < airfoil.size(); ++i)
	{
		switch (Passport().airfoilParams[i].panelsType)
		{
		case 0:			
			for (size_t i = 0; i < airfoil.size(); ++i)
				airfoil[i].reset(new AirfoilRect());
			break;
		case 1:
			/*
			for (size_t i = 0; i < airfoil.size(); ++i)
			airfoil[i].reset(new AirfoilCurv());
			*/
			break;
		}

		airfoil[i]->ReadFromFile("./Airfoils/", Passport().airfoilParams[i]);	//Считываем из каталога с коллекцией профилей

		switch (Passport().airfoilParams[i].boundaryCondition)
		{
		case 0:
			boundary[i].reset(new BoundaryVortColl(airfoil[i], wake, parallel));
			break;
		case 1:
			boundary[i].reset(new BoundaryConstLayerAver(airfoil[i], wake, parallel));
			break;
		}
	}

	switch (Passport().numericalSchemes.velocityComputation)
	{
	case 0:
		velocity.reset(new VelocityDirect(parallel, wake));
		break;
	case 1:
		velocity.reset(new VelocityFourier(parallel, wake));
		break;
	}

	teleFile << "time = " << currentTime << std::endl;
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
	return (currentTime >= Passport().timeDiscretizationProperties.timeStop);
}//isFinished()


//Основная функция выполнения одного шага по времени
void World2D::Step() 
{ 		
	const double& dt = Passport().timeDiscretizationProperties.dt;	
	
	std::cout << "N_PROCESSORS = " << parallel.nProcWork << std::endl;

	if (parallel.myidWork == 0)
	{
		//вычисляем размер матрицы и резервируем память под нее и под правую часть
		int matrSize = boundary.size();

		for (auto it = boundary.begin(); it != boundary.end(); ++it)
			matrSize += (*it)->GetUnknownsSize();

		matr.resize(matrSize, matrSize);
		matr.setZero();

		rhs.resize(matrSize);
		rhs.setZero();

		//TODO профиль пока строго 1 или строго 0
		if (airfoil.size() > 0 )
		{

			int nVars = boundary[0]->GetUnknownsSize();

			Eigen::MatrixXd locMatr;
			Eigen::VectorXd locLastLine, locLastCol;
			Eigen::VectorXd locRhs;

			locMatr.resize(nVars, nVars);
			locLastLine.resize(nVars);
			locLastCol.resize(nVars);
			locRhs.resize(nVars);

			boundary[0]->FillMatrixSelf(locMatr, locLastLine, locLastCol);
			boundary[0]->FillRhs(Passport().physicalProperties.V0, locRhs);

			//размазываем матрицу
			//TODO: переделать для общего случая
			for (int i = 0; i < nVars; ++i)
			{
				for (int j = 0; j < nVars; ++j)
					matr(i, j) = locMatr(i, j);
				matr(nVars, i) = locLastLine(i);
				matr(i, nVars) = locLastCol(i);

				rhs(i) = locRhs(i);
				//TODO totalGamma сейчас жестко 0
			}

			/*ofstream matrixFile("matrix.txt");
			SaveToStream(matr, matrixFile);
			matrixFile.close();

			ofstream rhsFile("rhs.txt");
			SaveToStream(rhs, rhsFile);
			rhsFile.close();*/

			sol = matr.partialPivLu().solve(rhs);

			//TODO профиль пока строго 1
			Eigen::VectorXd locSol;
			locSol.resize(nVars);
			for (int i = 0; i < nVars; ++i)
				locSol(i) = sol(i);
			boundary[0]->SolutionToFreeVortexSheet(locSol);

			//*
			std::ofstream solFile("solution.txt");
			SaveToStream(sol, solFile);
			solFile.close();
			//*/
		}
	}

	//TODO: Сброс вихрей

	//Сохранятель пелены ("старый")
	if (parallel.myidWork == 0)
	if (!(currentStep % Passport().timeDiscretizationProperties.deltacntText))
		wake.SaveKadr(Passport().dir, currentStep);

	//Движение вихрей (конвективное, след на себя, без вязкости)
	//Скорости
	velocity->CalcConvVelo(dt);
	
	//Влияние профиля
	std::vector<Point2D> SheetsVelo;
	//TODO профиль пока строго 1 или строго 0
	if (airfoil.size() > 0)
	{
		boundary[0]->GetWakeVelocity(SheetsVelo, dt);
	}

	//Смещения
	//TODO: "Погрузить" ее в какой-нибудь класс или в сам World2D
	std::vector<Point2D> newPos;
	if (parallel.myidWork == 0)
	for (int i = 0; i < wake.nv; ++i)
	{
		newPos.push_back(wake.vtx[i].r() + (velocity->convVelo[i] + SheetsVelo[i] + Passport().physicalProperties.V0*Passport().timeDiscretizationProperties.dt));
	}

	//TODO профиль пока строго 1 или строго 0
	if (airfoil.size() > 0)
	{
		std::vector<double> gamma;
		gamma = wake.Inside(newPos, *airfoil[0]);
		for (int i = 0; i < wake.nv; ++i)
			wake.vtx[i].r() = newPos[i];
	}

	std::cout << "tm = " << currentTime << std::endl;

	currentTime += dt;
	currentStep++;

	wake.WakeSinchronize();
}//Step()
