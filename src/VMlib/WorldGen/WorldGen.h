/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: WorldGen.h                                                       |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Заголовочный файл с описанием класса WorldGen
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef WORLDGEN_H
#define WORLDGEN_H

#include "TimesGen.h"

namespace VMlib
{
	class PassportGen;
	class TimesGen;
	
	/*!
	\brief Класс, опеделяющий текущую решаемую задачу
	\author Марчевский Илья Константинович
	\Version 1.14
	\date 6 марта 2026 г.
	*/
	class WorldGen
	{
	protected:
		/// Поток для вывода логов и сообщений об ошибках
		mutable LogStream info;

		/// Константная ссылка на паспорт конкретного расчета
		const PassportGen& passportGen;

		/// Сведения о временах выполнения основных операций
		std::unique_ptr<TimersGen> timers;

		/// Текущий номер шага в решаемой задаче
		size_t currentStep;

		/// Текущее время в решаемой задаче
		double currentTime;

	public:

		size_t nVtxBeforeMerging;

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
		/// const Parallel& getParallel() const { return parallel; };

		/// \brief Возврат номера текущего временного шага
		///
		/// \return номера текущего временного шага
		size_t getCurrentStep() const { return currentStep; };
		double getCurrentTime() const { return currentTime; };

		const PassportGen& getPassportGen() const { return passportGen; };

		/// \brief Функция, возвращающая признак завершения счета в решаемой задаче
		///
		/// true, если задача решена и выполнен признак останова; false если требуется еще выполнять шаги по времени
		bool isFinished() const;

		/// \brief Конструктор
		///
		/// \param[in] passport_ константная ссылка на паспорт расчета
		WorldGen(const VMlib::PassportGen& passport_);

		/// Деструктор
		virtual ~WorldGen() {};

		/// Функция выполнения предварительного шага
		//virtual void ZeroStep() = 0;

		/// Основная функция выполнения одного шага по времени	
		virtual void Step() = 0;
	};

}//namespace VMlib

#endif

