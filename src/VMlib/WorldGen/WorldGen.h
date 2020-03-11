/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.8    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2020/03/09     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: WorldGen.h                                                       |
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
\brief Заголовочный файл с описанием класса WorldGen
\author Марчевский Илья Константинович
\version 1.8   
\date 09 марта 2020 г.
*/

#ifndef WORLDGEN_H
#define WORLDGEN_H

#include <memory>

#include "TimesGen.h"

namespace VMlib
{
	class PassportGen;
	class Parallel;
	class TimesGen;
	
	/*!
	\brief Класс, опеделяющий текущую решаемую задачу
	\author Марчевский Илья Константинович
	\version 1.8
	\date 09 марта 2020 г.
	*/
	class WorldGen
	{
	protected:
		/// Поток для вывода логов и сообщений об ошибках
		mutable LogStream info;

		/// Константная ссылка на паспорт конкретного расчета
		const PassportGen& passportGen;

		/// Константная ссылка на параметры исполнения задачи в параллельном MPI-режиме
		const VMlib::Parallel& parallel;

		/// Сведения о временах выполнения основных операций
		std::unique_ptr<TimesGen> timestat;

	public:
		/// Текущий номер шага в решаемой задаче
		size_t currentStep;
		/// \brief Возврат ссылки на объект LogStream
		/// Используется в техничеcких целях для организации вывода
		///
		/// \return ссылку на объект LogStream
		VMlib::LogStream& getInfo() const { return info; };

		/// \brief Возврат ссылки на поток вывода информации
		/// Необходимо для вывода телеметрической информации, информации об ошибках и т.п.
		///
		/// \param[in] x символ, определяющий стиль вывода сообщения
		/// \return ссылку на поток вывода информации
		std::ostream& getInfo(char x) const { return info(x); };

		/// \brief Возврат константной ссылки на параметры распараллеливания по MPI
		///
		/// \return константную ссылку на параметры распараллеливания по MPI
		const Parallel& getParallel() const { return parallel; };

		/// \brief Возврат номера текущего временного шага
		///
		/// \return номера текущего временного шага
		size_t getCurrentStep() const { return currentStep; };

		const PassportGen& getPassportGen() { return passportGen; };

		/// \brief Функция, возвращающая признак завершения счета в решаемой задаче
		///
		/// true, если задача решена и выполнен признак останова; false если требуется еще выполнять шаги по времени
		bool isFinished() const;

		/// \brief Конструктор
		///
		/// \param[in] passport_ константная ссылка на паспорт расчета
		/// \param[in] parallel_ коенстантная ссылка на параметры исполнения задачи в параллельном MPI-режиме
		WorldGen(const VMlib::PassportGen& passport_, const VMlib::Parallel& parallel_);

		/// Основная функция выполнения одного шага по времени	
		virtual void Step() = 0;
	};

}//namespace VMlib

#endif

