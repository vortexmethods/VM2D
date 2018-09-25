/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.3    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/09/26     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Times.h                                                          |
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
\brief Заголовочный файл с описанием класса Times
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/ 

#ifndef TIMES_H
#define TIMES_H

#include "defs.h"

class World2D;

/*!
\brief Класс для сбора статистики времени исполнения основных шагов алгоритма и вывода ее в файл
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/
class Times
{
private:
	/// Константная ссылка на решаемую задачу
	const World2D& W;

	/// Обнуление одного временного периода 
	/// \param[out] period промежуток времени, начало и конец которого будут обнулены
	void ToZero(timePeriod& period)
	{
		period.first = 0;
		period.second = 0;
	}//ToZero(...)

public:
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
	timePeriod timeVelocityPressure;
	
	/// Конструктор
	Times(const World2D& W_)
		: W(W_) {};

	/// Деструктор
	~Times() {};
	
	/// Генерация заголовка файла временной статистики
	void GenerateStatHeader() const;

	/// Сохранение строки со статистикой в файл временной статистики
	void GenerateStatString() const;
	
	/// Обнуление состояния временной статистики
	void ToZero();

	/// Вычисление разницы во времени для пары засечек в секундах
	/// \param[in] t константная ссылка на пару засечек времени
	/// \return разницу в секундах
	static double dT(const timePeriod& t)
	{
		return (t.second - t.first);
	}//dT(...)
};

#endif