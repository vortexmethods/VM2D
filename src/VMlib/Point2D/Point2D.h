/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.12   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2024/01/14     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2024 Ilia Marchevsky                                     |
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
\Version 1.12
\date 14 января 2024 г.
*/

#ifndef POINT2D_H_
#define POINT2D_H_

#include "numvector.h"

#if defined(__CUDACC__)
#define HD __host__ __device__
#else
#define HD 
#endif


namespace VMlib
{
	struct TParticleCode
	{
		/// Мортоновский код частицы
		unsigned int key;

		/// Индекс частицы в глобальном массиве частиц
		int originNumber;
	};


	/*!
	\brief Класс, опеделяющий двумерный вектор

	Наследуется от numvector<double, 2>, имеет дополнительные возможности:
	- поворота на заданный угол против часовой стрелки;
	- генерируется MPI-описатель для возможности его пересылки как единичного объекта.

	\author Марчевский Илья Константинович
	\Version 1.12
	\date 14 января 2024 г.
	*/
	
	class Point2Df
		: public numvector<float, 2>
	{
	public:

		//#ifndef __CUDACC__
		//		/// MPI-описатель типа
		//		static MPI_Datatype mpiPoint2D;
		//#endif
				/// Пустой конструктор
		HD Point2Df() { };

		/// \brief Конструктор и приведение типа из numvector<double, 2>
		///
		/// \param[in] _r константная ссылка на копируемый объект типа numvector<double, 2>
		HD Point2Df(const numvector<float, 2>& _r) {
			data[0] = _r[0];
			data[1] = _r[1];
		};

		/// \brief Конструктор копирования
		///
		/// \param[in] _r константная ссылка на копируемый вектор
		HD Point2Df(const Point2Df& _r) {
			data[0] = _r[0];
			data[1] = _r[1];
		}

#if !defined(__CUDACC__)
		/// \brief Конструктор инициализации списком
		///
		/// \param[in] z константная ссылка на список инициализации из чисел типа double
		/// \warning Длина списка инициализации не проверяется, от него берутся только 2 первых элемента
		Point2Df(const std::initializer_list<float>& z);
#endif

		/// Деструктор
		HD ~Point2Df() { };

		/// \brief Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
		///
		/// \param[in] angle угол поворота в радианах (по умолчанию \f$ \frac{\pi}{2} \f$)
		/// \return новый вектор, полученный поворотом старого 
		HD Point2Df rotated(const float angle = 1.5707963267948966192313216916398f) const {
			Point2Df res;
			float cosa = cosf(angle);
			float sina = sinf(angle);

			res[0] = data[0] * cosa - data[1] * sina;
			res[1] = data[0] * sina + data[1] * cosa;
			return res;
		}

		//#ifndef __CUDACC__
		//		/// Cоздание MPI-описателя типа
		//		static void CreateMpiType();
		//#endif

				//operator numvector<double, 2>&() 
				//{
				//	return *this;
				//}
	};

	
	class Point2D
		: public numvector<double, 2>
	{
	public:

//#ifndef __CUDACC__
//		/// MPI-описатель типа
//		static MPI_Datatype mpiPoint2D;
//#endif
		/// Пустой конструктор
		HD Point2D() { };

		/// \brief Конструктор и приведение типа из numvector<double, 2>
		///
		/// \param[in] _r константная ссылка на копируемый объект типа numvector<double, 2>
		HD Point2D(const numvector<double, 2>& _r) {
			data[0] = _r[0];
			data[1] = _r[1];
		};

		/// \brief Конструктор копирования
		///
		/// \param[in] _r константная ссылка на копируемый вектор
		HD Point2D(const Point2D& _r) {
			data[0] = _r[0];
			data[1] = _r[1];
		}

#if !defined(__CUDACC__)
		/// \brief Конструктор инициализации списком
		///
		/// \param[in] z константная ссылка на список инициализации из чисел типа double
		/// \warning Длина списка инициализации не проверяется, от него берутся только 2 первых элемента
		Point2D(const std::initializer_list<double>& z);
#endif

		/// Деструктор
		HD ~Point2D() { };

		/// \brief Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
		///
		/// \param[in] angle угол поворота в радианах (по умолчанию \f$ \frac{\pi}{2} \f$)
		/// \return новый вектор, полученный поворотом старого 
		HD Point2D rotated(const double angle = 1.5707963267948966192313216916398) const {
			Point2D res;
			double cosa = cos(angle);
			double sina = sin(angle);

			res[0] = data[0] * cosa - data[1] * sina;
			res[1] = data[0] * sina + data[1] * cosa;
			return res;
		}

//#ifndef __CUDACC__
//		/// Cоздание MPI-описателя типа
//		static void CreateMpiType();
//#endif

		//operator numvector<double, 2>&() 
		//{
		//	return *this;
		//}
	};

}//namespace VMlib

using VMlib::Point2D;
using VMlib::Point2Df;
using VMlib::TParticleCode;


#endif
 