/*--------------------------------*- VMlib -*-----------------*---------------*\
| ##  ## ##   ## ##     #####   |                            | Version 1.5    |
| ##  ## ### ### ##     ##  ##  |  VMlib: VM2D/VM3D Library   | 2019/02/20     |
| ##  ## ## # ## ##     #####   |  Open Source Code          *----------------*
|  ####  ##   ## ##     ##  ##  |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: VM2D.cpp                                                         |
| Info: Source code of VMlib                                                   |
|                                                                             |
| This file is part of VMlib.                                                  |
| VMlib is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VMlib is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VMlib.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Основной файл программы VM2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.5
\date 20 февраля 2019 г.
*/

/*! 
\mainpage Вихревые методы для решения двумерных задач
Данный программный комплекс реализует вихревые методы решения двумерных задач гидродинамики и гидроупругости
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\date 20 февраля 2019 г.
\version 1.5
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

