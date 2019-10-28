/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.6    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2019/10/28     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Point2D.h                                                        |
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
\brief Заголовочный файл с описанием класса Point2D
\author Марчевский Илья Константинович
\version 1.6   
\date 28 октября 2019 г.
*/

#ifndef POINT2D_H_
#define POINT2D_H_

#ifndef __CUDACC__
	#include "mpi.h"
#endif

#include "numvector.h"

namespace VMlib
{

	/*!
	\brief Класс, опеделяющий двумерный вектор

	Наследуется от numvector<double, 2>, имеет дополнительные возможности:
	- поворота на заданный угол против часовой стрелки;
	- генерируется MPI-описатель для возможности его пересылки как единичного объекта.

	\author Марчевский Илья Константинович
	\version 1.6
	\date 28 октября 2019 г.
	*/
	class Point2D
		: public numvector<double, 2>
	{
	public:

#ifndef __CUDACC__
		/// MPI-описатель типа
		static MPI_Datatype mpiPoint2D;
#endif
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

#ifndef __CUDACC__
		/// Cоздание MPI-описателя типа
		static void CreateMpiType();
#endif
	};

}//namespace VMlib

using VMlib::Point2D;

#endif
 