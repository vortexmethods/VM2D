/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.7    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2019/11/22     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Task.h                                                           |
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
\brief Заголовочный файл с описанием класса Task
\author Марчевский Илья Константинович
\version 1.7   
\date 22 ноября 2019 г.
*/

#ifndef TASK_H
#define TASK_H

#include "PassportGen.h"

namespace VMlib
{

	/*!
	\brief Тип-перечисление, опеделяющий состояние задачи в процессе работы программы

	Может принимать 5 различных значений:
	 - waiting (0)   --- ожидает запуска;
	 - starting (1)  --- стартует;
	 - running (2)   --- считает;
	 - finishing (3) --- останавливается;
	 - done (4)      --- отсчитано;

	\author Марчевский Илья Константинович
	\version 1.7
	\date 22 ноября 2019 г.
	*/
	enum class TaskState
	{
		/// задача ожидает запуска
		waiting = 0,

		/// задача стартует
		starting = 1,

		/// задача решается
		running = 2,

		/// задача финиширует
		finishing = 3,

		/// задача решена
		done = 4
	};


	/*!
	\brief Класс, опеделяющий описание каждой задачи в смысле прохождения процесса ее решения
	\author Марчевский Илья Константинович
	\version 1.7
	\date 22 ноября 2019 г.
	*/
	class Task
	{
	public:

		/// \brief Конструктор
		///
		/// \param _passport константная ссылка на паспорт задачи
		//Task(std::unique_ptr<PassportGen> _passport)
		//{
		//	passport = std::move(_passport);
		//};

		/// Деструктор
		//~Task(){};

		//далее - данные, требуемые для обработки задачи очередью

		/// Число процессоров
		int nProc;

		/// Номера процессоров, занятых решением данной задачи
		std::vector<int> proc;

		/// Флаг состояния задачи	
		TaskState state;

		/// Номера квантов времени, в которых начался и завершился счет задачи
		std::pair<int, int> startEndKvant;

		/// \brief Время в секундах, потраченное на решение данной задачи
		///
		/// Заполняется после окончания решения задачи
		double wTime;

		/// Паспорт задачи	
		std::shared_ptr<PassportGen> passport;


		const PassportGen& getPassport() const
		{
			return *passport;
		}
	};

}//namespace VMlib

#endif