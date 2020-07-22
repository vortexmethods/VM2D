/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.9    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2020/07/22     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: TimesGen.h                                                       |
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
\brief Заголовочный файл с описанием класса TimesGen
\author Марчевский Илья Константинович
\version 1.9   
\date 22 июля 2020 г.
*/ 

#ifndef TIMESGEN_H
#define TIMESGEN_H

#include "defs.h"

namespace VMlib
{

	/*!
	\brief Класс для сбора статистики времени исполнения основных шагов алгоритма и вывода ее в файл
	\author Марчевский Илья Константинович
	\version 1.9
	\date 22 июля 2020 г.
	*/
	class TimesGen
	{
	protected:
		/// Обнуление одного временного периода 
		/// \param[out] period промежуток времени, начало и конец которого будут обнулены
		static void ToZeroPeriod(timePeriod& period)
		{
			period.first = 0;
			period.second = 0;
		}//ToZero(...)

	public:
		/// Конструктор
		TimesGen() {};

		/// Деструктор
		virtual ~TimesGen() {};


		/// Генерация заголовка файла временной статистики
		virtual void GenerateStatHeader() const = 0;

		/// Сохранение строки со статистикой в файл временной статистики
		virtual void GenerateStatString() const = 0;
		
		/// Обнуление состояния временной статистики
		virtual void ToZero() = 0;

		/// Вычисление разницы во времени для пары засечек в секундах
		/// \param[in] t константная ссылка на пару засечек времени
		/// \return разницу в секундах
		static double dT(const timePeriod& t)
		{
			return (t.second - t.first);
		}//dT(...)
	};

}//namespace VMlib

#endif