/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: v3D.h                                                            |
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
\brief Заголовочный файл с описанием класса v3D
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef V3D_H_
#define V3D_H_

#include "numvector.h"

namespace VMlib
{

	/*!
	\brief Класс, опеделяющий двумерный вектор

	Наследуется от numvector<double, 3>, имеет дополнительные возможности:
	- поворота на заданный угол вокруг заданного направления;
	- генерируется MPI-описатель для возможности его пересылки как единичного объекта.

	\author Марчевский Илья Константинович
	\Version 1.14
	\date 6 марта 2026 г.
	*/
	class v3D : public numvector<double, 3>
	{
	public:
		
//#ifndef __CUDACC__		
//		/// MPI-описатель типа
//		static MPI_Datatype mpiv3D;
//#endif

		/// Пустой конструктор
		v3D() { };

		/// \brief Конструктор и приведение типа из numvector<double, 3>
		///
		/// \param[in] _r константная ссылка на копируемый объект типа numvector<double, 3>
		v3D(const numvector<double, 3>& _r) 
		{
			data[0] = _r[0];
			data[1] = _r[1];
			data[2] = _r[2];
		};

		/// \brief Конструктор копирования
		///
		/// \param[in] _r константная ссылка на копируемый вектор
		v3D(const v3D& _r) 
		{
			data[0] = _r[0];
			data[1] = _r[1];
			data[2] = _r[2];
		};

#ifndef __CUDACC__
		/// \brief Конструктор инициализации списком
		///
		/// \param[in] z константная ссылка на список инициализации из чисел типа double
		/// \warning Длина списка инициализации не проверяется, от него берутся только 3 первых элемента
		v3D(const std::initializer_list<double>& z) 
		{
			for (size_t i = 0; i < 3; ++i)
				data[i] = *(z.begin() + i);
		};
#endif

		/// Деструктор
		~v3D() { };

		/// \brief Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
		///
		/// \param[in] angle угол поворота в радианах
		/// \param[in] axis константная ссылка на v3D, задающий ось поворота
		/// \return новый вектор, полученный поворотом старого 
		v3D rotated(const double angle, const v3D& axis) const 
		{
			/// \todo Реализовать!
			/// \warning Пока возвращает себя!

			v3D res = *this;
			//double cosa = cos(angle);
			//double sina = sin(angle);
			//
			//res[0] = r[0] * cosa - r[1] * sina;
			//res[1] = r[0] * sina + r[1] * cosa;
			return res;
		}//rotated(...);
	};


    /// \todo Исследовать целесообразность наличия явной перегрузки оператора умножения
	inline v3D operator*(double c, const v3D& x)
	{
		v3D res(x);
		for (size_t i = 0; i < 3; ++i)
			res[i] *= c;
		return res;
	}//operator*(...)


	inline std::pair<v3D, v3D>& operator+=(std::pair<v3D, v3D>& a, const std::pair<v3D, v3D>& b)
	{
		a.first += b.first;
		a.second += b.second;
		return a;
	}//operator+=(...)

	inline std::pair<v3D, v3D>& operator*=(std::pair<v3D, v3D>& a, double c)
	{
		a.first *= c;
		a.second *= c;
		return a;
	}//operator+=(...)

	inline std::pair<v3D, v3D> operator*(std::pair<v3D, v3D>& a, double c)
	{		
		return {a.first * c, a.second * c};
	}//operator+=(...)

}//namespace VMlib

using VMlib::v3D;

#endif
 