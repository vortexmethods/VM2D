/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.8    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/03/09     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Times2D.h                                                        |
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
\brief Заголовочный файл с описанием класса Times
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.8   
\date 09 марта 2020 г.
*/ 

#ifndef TIMES_H
#define TIMES_H

#include "TimesGen.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс для сбора статистики времени исполнения основных шагов алгоритма и вывода ее в файл
	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.8
	\date 09 марта 2020 г.
	*/
	class Times : public VMlib::TimesGen
	{
	private:
		/// Константная ссылка на решаемую задачу
		const World2D& W;

	public:
		/// Конструктор
		Times(const World2D& W_)
			: W(W_) {};

		/// Деструктор
		~Times() {};
		
		/// Начало и конец процесса выполнения шага целиком
		timePeriod timeWholeStep;

		/// Начало и конец процесса выделения памяти под матрицу и правую часть
		timePeriod timeReserveMemoryForMatrixAndRhs;

		/// Начало и конец процесса заполнения матрицы и формирования правой части
		timePeriod timeFillMatrixAndRhs;

		/// Начало и конец процесса решения системы линейных алгебраических уравнений
		timePeriod timeSolveLinearSystem;

		/// Начало и конец процесса вычисления конвективных скоростей вихрей
		timePeriod timeCalcVortexConvVelo;

		/// Начало и конец процесса вычисления диффузионных скоростей вихрей
		timePeriod timeCalcVortexDiffVelo;

		/// Начало и конец процесса вычисления нагрузок
		timePeriod timeGetHydroDynamForce;

		/// Начало и конец процесса перемещения вихрей
		timePeriod timeMoveVortexes;

		/// Начало и конец процесса контроля протыкания
		timePeriod timeCheckInside;

		/// Начало и конец процесса реструктуризации пелены
		timePeriod timeRestruct;

		/// Начало и конец процесса сортировки вихревого следа
		timePeriod timeWakeSort;

		/// Начало и конец процесса сохранения кадра в файл
		timePeriod timeSaveKadr;

		/// Начало и конец процесса подсчета полей скоростей и давления и сохранения их в файл
		timePeriod timeVP;

		/// Все прочее
		timePeriod timeOther;
		
		/// Генерация заголовка файла временной статистики
		virtual void GenerateStatHeader() const override;

		/// Сохранение строки со статистикой в файл временной статистики
		virtual void GenerateStatString() const override;
		
		//Обнуление состояния временной статистики
		virtual void ToZero() override;

	};

}//namespace VM2D

#endif