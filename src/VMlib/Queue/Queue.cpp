/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.9    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2020/07/22     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Queue.cpp                                                        |
| Info: Source code of VMlib                                                  |
|                                                                             |
| This file is part of VMlib.                                                 |
| VMLib is free software: you can redistribute it and/or modify it            |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VMlib is distributed in the hope that it will be useful, but WITHOUT        |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VMlib.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Файл кода с описанием класса Queue
\author Марчевский Илья Константинович
\version 1.9   
\date 22 июля 2020 г.
*/

#if defined(_WIN32)
 #include <direct.h>
#endif

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "Queue.h"

#include "Preprocessor.h"
#include "StreamParser.h"
#include "WorldGen.h"

#ifdef CODE2D
	#include "Airfoil2D.h"
	#include "Boundary2D.h"
	#include "MeasureVP2D.h"
	#include "Mechanics2D.h"
	#include "Passport2D.h"
	#include "Velocity2D.h"
	#include "Wake2D.h"
	#include "WakeDataBase2D.h"
	#include "World2D.h"
	#include "Tree2D.h"
#endif

#ifdef CODE3D
	#include "Body3D.h"
	#include "Passport3D.h"
	#include "Velocity3D.h"
	#include "Wake3D.h"
	#include "World3D.h"
#endif

using namespace VMlib;

//Конструктор
Queue::Queue(int& argc, char**& argv, void (*_CreateMpiTypes)())
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
		PrintUniversalLogoToStream(*defaults::defaultWorld2DLogStream);

		info.assignStream(defaults::defaultQueueLogStream, "queue");

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
	info('i') << "Goodbye!" << std::endl;
	MPI_Finalize();
}//~Queue()


//Процедура постановки новых задач на отсчет и занятие процессоров
void Queue::TaskSplit()
{		
	//Пересылаем информацию о состоянии процессоров на все компьютеры
	MPI_Scatter(procState.data(), 1, MPI_INT, &myProcState, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
	if (myProcState == MPI_UNDEFINED)
		world.reset(nullptr);
	
	if (myidAll == 0)
	{
		currentKvant++;

		int nfree = 0; //число свободных процессоров

		//Обнуляем число готовых к старту задач
		numberOfTask.prepared = 0;

		//считаем число свободных процессоров
		for (int i = 0; i < nProcAll; ++i)
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
					procState[j] = static_cast<int>(taskFol);

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
	
	//Вывод информации о занятости процессоров задачами
	if (myidAll == 0)
	{
		info.endl();
		info('i') << "ProcStates: " << std::endl;
		for (int i = 0; i < nProcAll; ++i)
			info('-') << "proc[" << i << "] <=> problem[" << procState[i] << "]" << std::endl;
		info.endl();
	}

	//Пересылаем информацию о состоянии процессоров на все компьютеры
	MPI_Scatter(procState.data(), 1, MPI_INT, &myProcState, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Синхронизируем информацию
	MPI_Bcast(&numberOfTask.solving, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&numberOfTask.prepared, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&numberOfTask.finished, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	//список головных процессоров на тех задачах, 
	//которые только что поставлены на счет:
	std::vector<int> prepList;
	
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
	}//if myidAll==0
		

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
			
			//подготовка стартующих задач
			parallel.myidWork = 0;
						
#ifdef CODE2D			
				world.reset(new VM2D::World2D(task[myProcState].getPassport(), parallel));
#endif

#ifdef CODE3D						
				world.reset(new VM3D::World3D(task[myProcState].getPassport(), parallel));
#endif			

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
		for (int s = 0; s<nProcAll; ++s)
		if ((procState[s] != MPI_UNDEFINED) && (task[procState[s]].state == TaskState::running) && (task[procState[s]].proc[0] == s))
			solvList.push_back(s);

		sizeCommSolving = static_cast<int>(solvList.size());

		flagFinish.clear();
		flagFinish.resize(solvList.size());

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
	}
		
	//Пересылаем информацию о состоянии процессоров на все компьютеры
	MPI_Scatter(procStateVar.data(), 1, MPI_INT, &myProcStateVar, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
	if (parallel.commWork != MPI_COMM_NULL)
	{
		//оперделяем "локальный" номер данного процессора в коммуникаторе commWork и число процессов в нем
		MPI_Comm_size(parallel.commWork, &parallel.nProcWork);
		MPI_Comm_rank(parallel.commWork, &parallel.myidWork);

		//World равен nullptr только в случае, когда выполняется первый шаг
		//в этом случае он создан только на главном процессоре группы;
		//на остальных происходит его создание
#ifdef CODE2D
			if (world == nullptr)
				world.reset(new VM2D::World2D(task[myProcState].getPassport(), parallel));
#endif

#ifdef CODE3D
			if (world == nullptr)
				world.reset(new VM3D::World3D(task[myProcState].getPassport(), parallel));
#endif
		

		//Синхронизация параметров
		MPI_Bcast(&(world->getPassportGen().timeDiscretizationProperties.currTime), 1, MPI_DOUBLE, 0, parallel.commWork);

		
		if ((parallel.myidWork == 0) && (world->getCurrentStep() == 0))
		{
#ifdef CODE2D					
			VM2D::World2D& world2D = dynamic_cast<VM2D::World2D&>(*world);
			//Создание файлов для записи сил
			for (size_t q = 0; q < world2D.getPassport().airfoilParams.size(); ++q)
				world2D.GenerateMechanicsHeader(q);
			//Создание файла для записи временной статистики
			world2D.getTimestat().GenerateStatHeader();
#endif

#ifdef CODE3D			
			VM3D::World3D& world3D = dynamic_cast<VM3D::World3D&>(*world);
			//Создание файлов для записи сил
			//for (size_t q = 0; q < world3D.getPassport().airfoilParams.size(); ++q)
			// world3D.GenerateMechanicsHeader(q);
			
			//Создание файла для записи временной статистики
			world3D.getTimestat().GenerateStatHeader();
#endif
		}
		//TODO
		//Синхронизация паспортов
		//MPI_Bcast(...)
	}

}//TaskSplit()


//Процедура обновления состояния задач и процессоров
void Queue::TaskUpdate()//Обновление состояния задач и процессоров
{
	info('i') << "------ Kvant finished ------" << std::endl;
	
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
		int stopSignal = ( (parallel.commWork == MPI_COMM_NULL) || (world->isFinished())) ? 1 : 0;
		
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


	if (parallel.commWork != MPI_COMM_NULL)
	{
		//К этому месту все процессоры, решающие задачи, приходят уже выполнив все расчеты
		//в рамках одного кванта времени
		//Поэтому коммуникаторы решаемых задач нам уже ни к чему - уничтожаем их
		MPI_Comm_free(&parallel.commWork);
	}

	//Обновление состояния задач и высвобождение процессоров из отсчитавших задач
	//(выполняется на глобальном нулевом процессоре)
	if (myidAll == 0)
	{
		//Список номеров задач, которые сейчас считаются
		std::vector<int> taskList;

		//если состояние "считается" запоминаем эту задачу
		//задачи запоминаются в порядке просмотра процессоров
		for (int s = 0; s<nProcAll; ++s)
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
		for (int i = 0; i < numberOfTask.solving; ++i)
		if (flagFinish[i + sizeCommSolving - numberOfTask.solving] == 1)
			task[taskList[i]].state = TaskState::finishing;
	} //if myid==0

	if (myidAll == 0)
	{
		for (size_t i = 0; i < task.size(); ++i)
		{
			//освобождаем процессоры от сосчитавшихся задач
			if (task[i].state == TaskState::finishing)
			{
				//устанавливаем состояние задачи в "отсчитано"
				task[i].state = TaskState::done;

				//отмечаем последний квант, когда задача считалась
				task[i].startEndKvant.second = currentKvant;
				
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

	for (int i = 0; i<nProcAll; ++i)
	if (procState[i] == MPI_UNDEFINED)
	{
		prFlag[i] = false;  //уже просмотрели
		procStateVar[i] = MPI_UNDEFINED;
	} //if (ProcState[i]==MPI_UNDEFINED)

	for (int i = 0; i<nProcAll; ++i)
	if (prFlag[i])
	{
		prFlag[i] = false;
		number++;
		
		for (int s = i; s<nProcAll; ++s)
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
		
		do
		{
			if (!world->isFinished())
				world->Step();
					
			if (parallel.myidWork == 0)
				deltaWallTime = MPI_Wtime() - kvantStartWallTime;
			MPI_Bcast(&deltaWallTime, 1, MPI_DOUBLE, 0, parallel.commWork);		
		}
		//Проверка окончания кванта по времени (или завершения задачи)
		while (deltaWallTime < kvantTime);

	}//if (commWork != MPI_COMM_NULL)
}//RunConveyer()


//Добавление задачи в список
void Queue::AddTask(int _nProc, std::unique_ptr<PassportGen> _passport)
{
	//Task newTask(std::move(_passport));
	//newTask.nProc = _nProc;
	//newTask.state = TaskState::waiting;
	//newTask.proc.resize(_nProc);
	//task.emplace_back(newTask);

	task.resize(task.size() + 1);
	(task.end() - 1)->nProc = _nProc;
	(task.end() - 1)->passport = std::move(_passport);
	(task.end() - 1)->proc.resize(_nProc);
	(task.end() - 1)->state = TaskState::waiting;



}//AddTask(...)


//Загрузка списка задач
void Queue::LoadTasksList(const std::string& _tasksFile, const std::string& _mechanicsFile, const std::string& _defaultsFile, const std::string& _switchersFile)
{
	if (
	fileExistTest(_tasksFile, info) && 
	fileExistTest(_defaultsFile, info) &&
	fileExistTest(_switchersFile, info)
	)
	{
		std::stringstream tasksFile(Preprocessor(_tasksFile).resultString);
		std::stringstream defaultsFile(Preprocessor(_defaultsFile).resultString);
		std::stringstream switchersFile(Preprocessor(_switchersFile).resultString);
		std::vector<std::string> taskFolder;

		std::unique_ptr<StreamParser> parserTaskList;
		parserTaskList.reset(new StreamParser(info, "parser", tasksFile, '(', ')'));

		std::vector<std::string> alltasks;
		parserTaskList->get("problems", alltasks);

		std::ptrdiff_t nTasks = std::count_if(alltasks.begin(), alltasks.end(), [](const std::string& a) {return (a.length() > 0);});
		info.endl();
		info('i') << "Number of problems to be solved: " << nTasks << std::endl;
		
		for (size_t i = 0; i < alltasks.size(); ++i)
		{
			if (alltasks[i].length() > 0)
			{
				//делим имя задачи + выражение в скобках на 2 подстроки
				std::pair<std::string, std::string> taskLine = StreamParser::SplitString(info, alltasks[i], false);

				std::string dir = taskLine.first;

				info.endl();
				info('i') << "-------- Loading problem #" << i << " (" << dir << ") --------" << std::endl;

				//вторую подстроку разделяем на вектор из строк по запятым, стоящим вне фигурных скобок		
				std::vector<std::string> vecTaskLineSecond = StreamParser::StringToVector(taskLine.second, '{', '}');

				//создаем парсер и связываем его с параметрами профиля
				std::stringstream tskStream(StreamParser::VectorStringToString(vecTaskLineSecond));
				std::unique_ptr<StreamParser> parserTask;

				parserTask.reset(new StreamParser(info, "problem parameters", tskStream, defaultsFile, {"pspfile", "np", "copyPspFile"}));

				std::string pspFile;
				int np;				

				//считываем нужные параметры с учетом default-значений
				parserTask->get("pspfile", pspFile, &defaults::defaultPspFile);
				parserTask->get("np", np, &defaults::defaultNp);				

				if (np > nProcAll)
				{
					info.endl();
					info('e') << "problem #" << i << " (" << dir << \
						") requires " << np << " processors (" << nProcAll << " available)" << std::endl;
					exit(-1);
				}
				else
				{					
					info('i') << "number of processors np = " << np << std::endl;
				}

				std::string copyPspFile;
				parserTask->get("copyPspFile", copyPspFile, &defaults::defaultCopyPspFile);
				if (copyPspFile.length() > 0)
				{
					//Копировать паспорт поручаем только одному процессу
					if (myidAll == 0)
					{
						size_t lastSlash = copyPspFile.find_last_of("\\/");
						std::string initDir(copyPspFile.begin(), copyPspFile.begin()+lastSlash);
						std::string command;

#if defined(_WIN32)
						std::replace(dir.begin(), dir.end(), '/', '\\');
						std::replace(initDir.begin(), initDir.end(), '/', '\\');
#endif

						VMlib::CreateDirectory(dir.c_str(), "");

#if defined(_WIN32)		
						command = "copy \"" + initDir + "\\*.*\" "+ dir + "\\";
						//system("dir");
#else
						command = "cp ./" + initDir + "/* " + dir + "/";
#endif

						std::cout << "Copying files from folder \"" << initDir << "\" to \"" << dir << "\"" << std::endl;
						
						int systemRet = system(command.c_str());
						if(systemRet == -1)
						{
							// The system method failed
							info('e') << "problem #" << i << " (" << dir << \
										 ") copying passport system method failed" << std::endl;
							exit(-1);
						}					
						std::cout << "Copying OK " << std::endl << std::endl;

						//copyFile(pspFile, "./" + dir + "/" + pspFile);
					}
				}

				MPI_Barrier(MPI_COMM_WORLD);

				std::unique_ptr<PassportGen> ptrPsp;
				
#ifdef CODE2D				
				ptrPsp.reset(new VM2D::Passport(info, dir, i, pspFile, _mechanicsFile, _defaultsFile, _switchersFile, vecTaskLineSecond));
#endif

#ifdef CODE3D			
				ptrPsp.reset(new VM3D::Passport(info, dir, i, pspFile, _mechanicsFile, _defaultsFile, _switchersFile, vecTaskLineSecond));
#endif				
								
				AddTask(np, std::move(ptrPsp));

				info('i') << "-------- Problem #" << i << " (" << dir << ") is loaded --------" << std::endl << std::endl;
			} 
		} //for i

		//tasksFile.close();
		tasksFile.clear();
	
		//defaultsFile.close();
		defaultsFile.clear();

		//switchersFile.close();
		switchersFile.clear();

	}
	else
	{
		exit(1);
	}

}//Queue::LoadTasks()