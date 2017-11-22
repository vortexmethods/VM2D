/*!
\file
\brief Основной файл программы VM2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

/*! 
\mainpage Вихревые методы для решения двумерных задач
Данный программный комплекс реализует вихревые методы решения двумерных задач гидродинамики и гидроупругости
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\date 01 декабря 2017
\version 1.0
*/

/// \defgroup Parallel Параллельные функции

#include "Queue.h"

void CreateMpiTypes()
{
	Vortex2D::CreateMpiType();
	Point2D::CreateMpiType();
}

void Initializers()
{
	
}


int main(int argc, char** argv)
{	
	Initializers();
	
	Queue queue(argc, argv, std::cout, CreateMpiTypes); //третий аргумент - timeFile

	queue.LoadTasksList("problems", "defaults", "switchers");

	do
	{
		queue.TaskSplit();
		
		queue.RunConveyer();

		queue.TaskUpdate();		
	} 
	while (queue.nextKvant);

	//cin.get();
}

