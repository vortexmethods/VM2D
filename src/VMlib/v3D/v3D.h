/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.5    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2019/02/20     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: v3D.h                                                            |
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
\brief Заголовочный файл с описанием класса v3D
\author Марчевский Илья Константинович
\version 1.5
\date 20 февраля 2019 г.
*/

#ifndef V3D_H_
#define V3D_H_

#include "mpi.h"

#include "numvector.h"

namespace VMlib
{

	/*!
	\brief Класс, опеделяющий двумерный вектор

	Наследуется от numvector<double, 3>, имеет дополнительные возможности:
	- поворота на заданный угол вокруг заданного направления;
	- генерируется MPI-описатель для возможности его пересылки как единичного объекта.

	\author Марчевский Илья Константинович
	\version 1.5
	\date 20 февраля 2019 г.
	*/
	class v3D : public numvector<double, 3>
	{
	public:
		/// MPI-описатель типа
		static MPI_Datatype mpiv3D;

		/// Пустой конструктор
		v3D() { };

		/// \brief Конструктор и приведение типа из numvector<double, 3>
		///
		/// \param[in] _r константная ссылка на копируемый объект типа numvector<double, 3>
		v3D(const numvector<double, 3>& _r);

		/// \brief Конструктор копирования
		///
		/// \param[in] _r константная ссылка на копируемый вектор
		v3D(const v3D& _r);

#if !defined(__CUDACC__)
		/// \brief Конструктор инициализации списком
		///
		/// \param[in] z константная ссылка на список инициализации из чисел типа double
		/// \warning Длина списка инициализации не проверяется, от него берутся только 3 первых элемента
		v3D(const std::initializer_list<double>& z);
#endif

		/// Деструктор
		~v3D() { };

		/// \brief Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
		///
		/// \param[in] angle угол поворота в радианах
		/// \param[in] axis константная ссылка на v3D, задающий ось поворота
		/// \return новый вектор, полученный поворотом старого 
		v3D rotated(const double angle, const v3D& axis) const;


		/// Cоздание MPI-описателя типа
		static void CreateMpiType();
	};


	inline std::pair<v3D, v3D>& operator+=(std::pair<v3D, v3D>& a, const std::pair<v3D, v3D>& b)
	{
		a.first += b.first;
		a.second += b.second;
		return a;
	}//operator+=(...)

}//namespace VMlib

using VMlib::v3D;

#endif
 