/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.8    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/03/09     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: VM2D.cpp                                                         |
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
\brief Основной файл программы VM2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.8
\date 09 марта 2020 г.
*/

/*! 
\mainpage Вихревые методы для решения двумерных задач
Данный программный комплекс реализует вихревые методы решения двумерных задач гидродинамики и гидроупругости
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\date 09 марта 2020 г.
\version 1.8
*/

/// \defgroup Parallel Параллельные функции

#include "Queue.h"

using namespace VMlib;

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

	Queue queue(argc, argv, CreateMpiTypes);

	queue.LoadTasksList("problems", "mechanics", "defaults", "switchers");	
	
	do
	{
		queue.TaskSplit();
		
		queue.RunConveyer();

		queue.TaskUpdate();		
	} 
	while (queue.nextKvant);

	//cin.get();
}

