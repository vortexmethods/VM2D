/*!
\file
\brief Файл кода с описанием класса Queue
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include "Queue.h"


//Конструктор
Queue::Queue(int& argc, char**& argv, std::ostream& _timeFile, void (*_CreateMpiTypes)())
	: timeFile(_timeFile)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcAll);
	MPI_Comm_rank(MPI_COMM_WORLD, &myidAll);
	MPI_Comm_group(MPI_COMM_WORLD, &groupAll);
	_CreateMpiTypes();	

	//готовим систему к запуску
	//(выполняется на глобальном нулевом процессоре)
	if (myidAll == 0)
	{
		//TODO: Переделать вывод в файл (вывод силы)
		std::ofstream forcesFile("Forces.txt");
		forcesFile << "currentStep	currentTime	Fx	Fy\n";
		forcesFile.close();

		//TODO: Переделать вывод в файл (вывод силы)
		std::ofstream stressFile("Stresses.txt");
		stressFile.close();

		//Устанавливаем флаг занятости в состояние "свободен" всем процессорам, 	
		procState.resize(nProcAll, MPI_UNDEFINED);
		procStateVar.resize(nProcAll, MPI_UNDEFINED);

		numberOfTask.solving  = 0; //число решаемых в данный момент задач
		numberOfTask.prepared = 0; //число подготовленных к запуску задач (в том числе)
		numberOfTask.finished = 0; //число уже отсчитанных задач	

		//текущий номер кванта времени
		currentKvant = -1;
	}//if myid==0 
}//Queue()


//Деструктор
Queue::~Queue()
{
	MPI_Finalize();
}//~Queue()


//Процедура постановка новых задач на отсчет и занятие процессоров
void Queue::TaskSplit()
{	
	//Пересылаем информацию о состоянии процессоров на все компьютеры
	MPI_Scatter(procState.data(), 1, MPI_INT, &myProcState, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (myProcState == MPI_UNDEFINED)
		world2D.reset(nullptr);

	if (myidAll == 0)
	{
		currentKvant++;

		int nfree = 0; //число свободных процессоров

		//Обнуляем число готовых к старту задач
		numberOfTask.prepared = 0;

		//считаем число свободных процессоров
		for (int i = 0; i < nProcAll; i++)
		if (procState[i] == MPI_UNDEFINED)
			nfree++;
	
		//формируем номер следующей задачи, которая возможно будет поставлена на отсчет:
		size_t taskFol = numberOfTask.finished + numberOfTask.solving;

		//проверяем, достаточно ли свободных процессоров для запуска еще одной задачи
		//при условии, что не все задачи уже решены
		while ((taskFol < task.size()) && (nfree >= task[taskFol].nProc))
		{
			int p = 0; //число найденных свободных процессоров под задачу
			int j = 0; //текущий счетчик процессоров

			//подбираем свободные процессоры
			do
			{
				if (procState[j] == MPI_UNDEFINED)
				{
					//состояние процессора устанавливаем в "занят текущей задачей"
					procState[j] = taskFol;

					//номер процессора посылаем в перечень процессоров, решающих текущую задачу
					task[taskFol].proc[p] = j;

					//увеличиваем число найденных процессоров на единицу
					p++;
				}//if (ProcState[j]==MPI_UNDEFINED)

				//переходим к следующему процессору
				j++;
			} while (p < task[taskFol].nProc);

			//состояние задачи устанавливаем в режим "стартует"
			task[taskFol].state = TaskState::starting;

			//отмечаем номер кванта, когда задача начала считаться
			task[taskFol].startEndKvant.first = currentKvant;

			//изменяем счетчики числа задач
			numberOfTask.prepared++;
			numberOfTask.solving++;

			//изменяем счетчик свободных процессоров
			nfree -= task[taskFol].nProc;

			//переходим к следующей задаче
			taskFol++;

		}//while ((nfree >= Task[task_fol].nproc)&&(task_fol<Task.size()))

	}//if (myidAll == 0)

	//Пересылаем информацию о состоянии процессоров на все компьютеры
	MPI_Scatter(procState.data(), 1, MPI_INT, &myProcState, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//список головных процессоров на тех задачах, 
	//которые только что поставлены на счет:
	std::vector<int> prepList;
	
	//список паспортов для стартующих задач
	//vector<Passport> pspList;

	//формируем списки prepList и pspList
	//(выполняется на глобальном нулевом процессоре)
	if (myidAll == 0)
	{
		for (size_t i = 0; i<task.size(); ++i)
		//находим задачи в состоянии "стартует"
		if (task[i].state == TaskState::starting)
		{
			//запоминаем номер того процессора, что является там головным
			prepList.push_back(task[i].proc[0]);
			
			//состояние задачи вереводим в "считает"
			task[i].state = TaskState::running;
		}
	}//if myid==0
		

	//Рассылаем количество стартующих задач в данном кванте на все машины
	MPI_Bcast(&numberOfTask.prepared, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (myidAll > 0)
	{
		prepList.resize(numberOfTask.prepared);
	}//if myid==0

	//следующий фрагмент выполняется только если в данном кванте стартуют новые задачи
	if (numberOfTask.prepared > 0)
	{
		//пересылаем список головных машин стартующих задач на все процессоры
		MPI_Bcast(prepList.data(), numberOfTask.prepared, MPI_INT, 0, MPI_COMM_WORLD);
		
		//формируем группу и коммуникатор стартующих головных процессоров
		commStarting = MPI_COMM_NULL;
		MPI_Group_incl(groupAll, numberOfTask.prepared, prepList.data(), &groupStarting);
		MPI_Comm_create(MPI_COMM_WORLD, groupStarting, &commStarting);
					
		//следующий фрагмент кода только для стартующих процессов, 
		//которые входят в коммуникатор comm_starting
		if (commStarting != MPI_COMM_NULL)
		{
			//получаем номер процесса в списке стартующих
			int myidStarting;
			MPI_Comm_rank(commStarting, &myidStarting);
			
			//TODO
			//<Тут будет подготовка стартующих задач>
			world2D.reset(new World2D(*this, std::cout));

			//Коммуникатор головных процессоров стартующих задач
			//выполнил свое черное дело и будет удален
			MPI_Comm_free(&commStarting);
		} //if(comm_starting != MPI_COMM_NULL)

		//их группа тоже больше без надобности
		MPI_Group_free(&groupStarting);
	}//if (task_prepared>0)

	//формируем группу и коммуникатор считающих головных процессоров
	//в том числе тех, которые стартуют, и тех, которые ранее входили в commStarting
	std::vector<int> solvList; //список номеров головных процессов

	if (myidAll == 0)
	{
		//Если нулевой просессор свободен 
		//(не явлется головным в решении задачи в данном кванте) - 
		//- формально присоединяем его в группу головных
		//необходимо для корректного обмена данными
		if (procState[0] == MPI_UNDEFINED)
			solvList.push_back(0);


		//далее ищем считающие задачи (в том числе только что стартующие)
		//и собираем номера их головных процессоров
		//(собираем в порядке просмотра номеров процессоров))
		for (int s = 0; s<nProcAll; s++)
		if ((procState[s] != MPI_UNDEFINED) && (task[procState[s]].state == TaskState::running) && (task[procState[s]].proc[0] == s))
			solvList.push_back(s);

		sizeCommSolving = solvList.size();

		flagFinish.clear();
		flagFinish.resize(solvList.size());

		if (myidAll == 0)
		for (size_t q = 0; q < solvList.size(); q++)
			std::cout << "solvList[" << q << "] = " << solvList[q] << std::endl;
	}//if myidAll==0


	//число sizeCommSolving на единицу больше, чем taskSolving, если нулевой процессор свободен
	//и равно taskSolving, если он занят решением задачи
	//(если он занят решением задачи, то непременно является головным)

	//пересылаем число головных процессоров и их список на все компьютеры
	//(в том числе головной процессор - независимо от его состояния)
	MPI_Bcast(&sizeCommSolving, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (myidAll > 0) 
		solvList.resize(sizeCommSolving);
	
	MPI_Bcast(solvList.data(), sizeCommSolving, MPI_INT, 0, MPI_COMM_WORLD);

	commSolving = MPI_COMM_NULL;

	//формируем группу и коммуникатор головных процессоров
	MPI_Group_incl(groupAll, sizeCommSolving, solvList.data(), &groupSolving);
	MPI_Comm_create(MPI_COMM_WORLD, groupSolving, &commSolving);

	if (myidAll == 0)
	{
		//составление неубывающего массива номеров для верного
		//распределения задач по процессоров, с тем чтобы задачи
		//с прошлого кванта оказались в тех же группах, что и раньше
		ConstructProcStateVar();	

		if (myidAll == 0)
		for (int q = 0; q < nProcAll; q++)
			std::cout << "procStateVar[" << q << "] = " << procStateVar[q] << std::endl;
	}
		
	//Пересылаем информацию о состоянии процессоров на все компьютеры
	MPI_Scatter(procStateVar.data(), 1, MPI_INT, &myProcStateVar, 1, MPI_INT, 0, MPI_COMM_WORLD);

	std::cout << "id = " << myidAll << ", myProcStateVar = " << myProcStateVar << std::endl;

	//Расщепляем весь комммуникатор MPI_COMM_WORLD (т.е. вообще все процессоры) 
	//на коммуникаторы решаемых задач
	//В результате все процессоры, решающие конкретную задачу,
	//объединяются в коммуникатор comm_work
	//Несмотря на то, что он называется везде одинаково, у каждой задачи он свой
	//и объединяет ровно те процессоры, которые надо
	//Процессоры, состояние которых установлено в "свободен", т.е.
	//ProcState=MPI_UNDEFINED (=-32766) ни в один 
	//новый коммуникатор commWork не попадут
	MPI_Comm_split(MPI_COMM_WORLD, myProcStateVar, 0, &parallel.commWork);

	//Рассылаем величину кванта времени
	//MPI_Bcast(&kvantTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//следующий код выполняется процессорами, участвующими в решении задач, 
	//следовательно, входящими в коммуникаторы comm_work
	
	std::cout << "myid = " << myidAll << std::endl;
	
	if (parallel.commWork != MPI_COMM_NULL)
	{
		std::cout << "inside_myid = " << myidAll << std::endl;
		//оперделяем "локальный" номер данного процессора в коммуникаторе commWork и число процессов в нем
		MPI_Comm_size(parallel.commWork, &parallel.nProcWork);
		MPI_Comm_rank(parallel.commWork, &parallel.myidWork);

		//World2D равен nullptr только в случае, когда выполняется первый шаг
		//в этом случае он создан только на главном процессоре группы;
		//на остальных происходит его создание
		if (world2D == nullptr)
			world2D.reset(new World2D(*this, *(defaults::defaultPtele)));

		//Синхронизация параметров
		MPI_Bcast(&(world2D->Passport().physicalProperties.currTime), 1, MPI_DOUBLE, 0, parallel.commWork);

		//TODO
		//Синхронизация паспортов
		//MPI_Bcast(...)
	}
}//TaskSplit()


//Процедура обновления состояния задач и процессоров
void Queue::TaskUpdate()//Обновление состояния задач и процессоров
{
	if (parallel.commWork != MPI_COMM_NULL)
	{
		//К этому месту все процессоры, решающие задачи, приходят уже выполнив все расчеты
		//в рамках одного кванта времени
		//Поэтому коммуникаторы решаемых задач нам уже ни к чему - уничтожаем их
		MPI_Comm_free(&parallel.commWork);
	}

	std::cout << "id = " << myidAll << ", commWork cleared!" << std::endl;

	//Код только для головных процессов, которые входят в коммуникатор commSolving
	//т.е. выполняется только на головных процессорах, решающих задачи.
	//Кажется, что можно было бы с тем же эфектом написать здесь if (myidWork==0),
	//но это не так, поскольку нужно включить в себя еще нулевой процессор,
	//который в любом случае присоединен к коммуникатору comm_solving,
	//даже если он не решает задачу (а если решает - то он всегда головной)
	if (commSolving != MPI_COMM_NULL)
	{
		//определяем номер процессора в данном коммуникаторе - 
		// - вдруг потребуется на будущее!
		int myidSolving;
		MPI_Comm_rank(commSolving, &myidSolving);

		//Алгоритм возвращения результатов расчета конкретной задачи		
		//Если срабатывает признак того, что решение задачи можно прекращать
		//или не решается ни одна задача (а в комуникатор входит
		//только 0-й процессор, присоединеннй туда насильно)
		int stopSignal = ((numberOfTask.solving == 0) || (world2D->isFinished())) ? 1 : 0;

		//пересылаем признаки прекращения решения задач на нулевой процессор 
		//коммуникатора comm_solving - т.е. на глобальный процессор с номером 0 ---
		// --- вот для чего его насильно присоединяли к этому коммуникатору
		//если сам по себе он туда не входил		
		MPI_Gather(&stopSignal, 1, MPI_INT, flagFinish.data(), 1, MPI_INT, 0, commSolving);
		
		//после пересылки информации коммуникатор comm_solving уничтожается
		MPI_Comm_free(&commSolving);

	} //if(comm_solving != MPI_COMM_NULL)

	//и соотвествующая ему группа тоже удаляется
	MPI_Group_free(&groupSolving);

	//Обновление состояния задач и высвобождение процессоров из отсчитавших задач
	//(выполняется на глобальном нулевом процессоре)
	if (myidAll == 0)
	{
		//Список номеров задач, которые сейчас считаются
		std::vector<int> taskList;

		//если состояние "считается" запоминаем эту задачу
		//задачи запоминаются в порядке просмотра процессоров
		for (int s = 0; s<nProcAll; s++)
		if ( (procState[s] != MPI_UNDEFINED) && (task[procState[s]].state == TaskState::running) && (task[procState[s]].proc[0] == s) )
			taskList.push_back(procState[s]);
			
		//Тем задачам, которые во флаге прекращения счета вернули "1", 
		//ставим состояние завершения счета.
		//Конструкция +sizeCommSolving-taskSolving введена для смещения в массиве на единицу
		//если нулевой процессор свободен - тогда он присоединяется формально в нулевом
		//элементе массива flagFinish, а содержательная часть массива смещается на 1
		//(в этом случае sizeCommSolving как раз будет на 1 больше, чем taskSolving, см. выше)
		//если же нулевой процесс занят решением задачи, то смещения нет,
		//поскольку в этом случае sizeCommSolving и taskSolving равны
		for (int i = 0; i < numberOfTask.solving; i++)
		if (flagFinish[i + sizeCommSolving - numberOfTask.solving] == 1)
			task[taskList[i]].state = TaskState::finishing;
	} //if myid==0

	if (myidAll == 0)
	{
		for (size_t i = 0; i < task.size(); i++)
		{
			//освобождаем процессоры от сосчитавшихся задач
			if (task[i].state == TaskState::finishing)
			{
				//устанавливаем состояние задачи в "отсчитано"
				task[i].state = TaskState::done;

				//отмечаем последний квант, когда задача считалась
				task[i].startEndKvant.second = currentKvant;

				//отметка в общий файл статистики о времени счета данной задачи (с точностью до кванта)
				//timeFile << "Task " << task[i].passport.GetFilePassport() << " time (rounded upto kvantTime) = " << (task[i].endKvant - task[i].startKvant + 1)*kvantTime << endl;

				//изменяем счетчики
				numberOfTask.solving--;
				numberOfTask.finished++;

				//состояние процессоров, решавших данную задачу, устанавливаем в "свободен"
				for (int p = 0; p < task[i].nProc; ++p)
				{
					procState[task[i].proc[p]] = MPI_UNDEFINED;
				}//for p
			}//if (Task[i].state == 3)
		}//for i

		//определяем, надо ли еще выделять квант времени или можно заканчивать работу
		//путем сравнения числа отсчитанных задач с общим числом задач
		nextKvant = (numberOfTask.finished < static_cast<int>(task.size())) ? 1 : 0;
	}//if (myid == 0)

	MPI_Bcast(&nextKvant, 1, MPI_INT, 0, MPI_COMM_WORLD);
}//TaskUpdate()


//Процедура, нумерующая задачи в возрастающем порядке
void Queue::ConstructProcStateVar()
{
	std::vector<bool> prFlag;
	int number = -1;

	prFlag.resize(nProcAll, true); //их надо просматривать

	for (int i = 0; i<nProcAll; i++)
	if (procState[i] == MPI_UNDEFINED)
	{
		prFlag[i] = false;  //уже просмотрели
		procStateVar[i] = MPI_UNDEFINED;
	} //if (ProcState[i]==MPI_UNDEFINED)

	for (int i = 0; i<nProcAll; i++)
	if (prFlag[i])
	{
		prFlag[i] = false;
		number++;
		
		for (int s = i; s<nProcAll; s++)
		if (procState[s] == procState[i])
		{
			procStateVar[s] = number;
			prFlag[s] = false;
		}//if (ProcState[s]==ProcState[i])
	} //if (pr_flag[i])									
}//ConstructProcStateVar()


// Запуск вычислительного конвейера (в рамках кванта времени)
void Queue::RunConveyer()
{
	if (parallel.commWork != MPI_COMM_NULL)
	{
		double kvantStartWallTime; //Время на главном процессоре, когда начался очередной квант
		double deltaWallTime;

		if (parallel.myidWork == 0)
			kvantStartWallTime = MPI_Wtime();

		std::cout << "id = " << myidAll << " - conv. started!" << std::endl;

		do
		{
			if (!world2D->isFinished())
				world2D->Step();
			
			if (parallel.myidWork == 0)
				deltaWallTime = MPI_Wtime() - kvantStartWallTime;
			MPI_Bcast(&deltaWallTime, 1, MPI_DOUBLE, 0, parallel.commWork);
		}
		//Проверка окончания кванта по времени (или завершения задачи)
		while (deltaWallTime < kvantTime);

		std::cout << "id = " << myidAll << " - conv. finished!" << std::endl;
	}//if (commWork != MPI_COMM_NULL)
}//RunConveyer()


//Добавление задачи в список
void Queue::AddTask(int _nProc, const Passport& _passport)
{
	Task newTask(_passport);
	newTask.nProc = _nProc;
	newTask.state = TaskState::waiting;
	newTask.proc.resize(_nProc);

	task.push_back(newTask);
}//AddTask(...)


//Загрузка списка задач
void Queue::LoadTasksList(const std::string& _tasksFile, const std::string& _defaultsFile, const std::string& _switchersFile)
{
	std::ifstream tasksFile(_tasksFile);
	std::ifstream defaultsFile(_defaultsFile);
	std::ifstream switchersFile(_switchersFile);
	std::vector<std::string> taskFolder;
		
	//создаем парсер и связываем его с нужными потоками, на выход получаем список папок с задачами
	std::unique_ptr<StreamParser> parserTask;
	parserTask.reset(new StreamParser(tasksFile, defaultsFile, taskFolder));
	
	for (size_t q = 0; q < taskFolder.size(); ++q)
	{
		*defaults::defaultPinfo << "----- Problem #" << q << " is being loaded -----" << std::endl;
		
		//считываем параметры задачи, если они есть
		std::vector<std::string> taskLine;
		std::string dir = taskFolder[q];
		parserTask->get(dir, taskLine);

		std::string pspFile;
		int np;
		
		std::unique_ptr<StreamParser> parserOneTask;
		std::stringstream paramStream(StreamParser::VectorStringToString(taskLine));
		parserOneTask.reset(new StreamParser(paramStream, defaultsFile));
		
		parserOneTask->get("pspfile", pspFile, &defaults::defaultPspFile);
		parserOneTask->get("np", np, &defaults::defaultNp);
		
		Passport psp("./" + dir + "/", pspFile, _defaultsFile, _switchersFile, taskLine);
		AddTask(np, psp);

		*defaults::defaultPinfo << "-------- Problem #" << q << " is loaded --------" << std::endl;
		*defaults::defaultPinfo << std::endl;
	}
	tasksFile.close();
	tasksFile.clear();
	
	defaultsFile.close();
	defaultsFile.clear();

	switchersFile.close();
	switchersFile.clear();
}//Queue::LoadTasks()