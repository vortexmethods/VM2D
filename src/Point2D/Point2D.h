/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Point2D.h                                                        |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.	                                      |
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
\brief Заголовочный файл с описанием класса Point2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/

#ifndef POINT2D_H_
#define POINT2D_H_

#include <cmath>

#include "mpi.h"

#include "numvector.h"

/*!
\brief Класс, опеделяющий двумерный вектор

Наследуется от numvector<double, 2>, имеет дополнительные возможности:
- поворота на заданный угол против часовой стрелки;
- генерируется MPI-описатель для возможности его пересылки как единичного объекта.

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/
class Point2D 
	: public numvector<double, 2>
{
public:	
	/// MPI-описатель типа
	static MPI_Datatype mpiPoint2D;
	
	/// Пустой конструктор
	Point2D() { };

	/// \brief Конструктор и приведение типа из numvector<double, 2>
	///
	/// \param[in] _r константная ссылка на копируемый объект типа numvector<double, 2>
	Point2D(const numvector<double, 2>& _r);
	
	/// \brief Конструктор копирования
	///
	/// \param[in] _r константная ссылка на копируемый вектор
	Point2D(const Point2D& _r);
	
#if !defined(__CUDACC__)
	/// \brief Конструктор инициализации списком
	///
	/// \param[in] z константная ссылка на список инициализации из чисел типа double
	/// \warning Длина списка инициализации не проверяется, от него берутся только 2 первых элемента
	Point2D(const std::initializer_list<double>& z);
#endif

	/// Деструктор
	~Point2D() { };

	/// \brief Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
	///
	/// \param[in] angle угол поворота в радианах (по умолчанию \f$ \frac{\pi}{2} \f$)
	/// \return новый вектор, полученный поворотом старого 
	Point2D rotated(const double angle = 1.5707963267948966192313216916398) const;
	
	/// Cоздание MPI-описателя типа
	static void CreateMpiType();	
};


#endif
 