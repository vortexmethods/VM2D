/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.8    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2020/03/09     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: PairInt.h                                                        |
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
\brief Файл кода с описанием класса PairInt
\author Марчевский Илья Константинович
\version 1.8   
\date 09 марта 2020 г.
*/


/*!
\file
\brief Заголовочный файл с описанием класса pairInt
\author Марчевский Илья Константинович
\version 1.8
\date 09 марта 2020 г.
*/

#ifndef PAIRINT_H_
#define PAIRINT_H_

#include "mpi.h"

#include "numvector.h"

namespace VMlib
{

	/*!
	\brief Класс, опеделяющий пару целых чисел

	Наследуется от numvector<int, 2>, имеет дополнительные возможности:
	- поворота на заданный угол вокруг заданного направления;
	- генерируется MPI-описатель для возможности его пересылки как единичного объекта.

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.8
	\date 09 марта 2020 г.
	*/
	class PairInt : public numvector<int, 2>
	{
	public:
		/// MPI-описатель типа
		static MPI_Datatype mpiPairInt;

		/// Пустой конструктор
		PairInt() { };

		/// \brief Конструктор и приведение типа из numvector<int, 2>
		///
		/// \param[in] _r константная ссылка на копируемый объект типа numvector<int, 2>
		PairInt(const numvector<int, 2>& _r);

		/// \brief Конструктор копирования
		///
		/// \param[in] _r константная ссылка на копируемый вектор
		PairInt(const PairInt& _r);

#if !defined(__CUDACC__)
		/// \brief Конструктор инициализации списком
		///
		/// \param[in] z константная ссылка на список инициализации из чисел типа int
		/// \warning Длина списка инициализации не проверяется, от него берутся только 2 первых элемента
		PairInt(const std::initializer_list<int>& z);
#endif

		/// Деструктор
		~PairInt() { };

		/// Cоздание MPI-описателя типа
		static void CreateMpiType();
	};


	inline std::pair<PairInt, PairInt>& operator+=(std::pair<PairInt, PairInt>& a, const std::pair<PairInt, PairInt>& b)
	{
		a.first += b.first;
		a.second += b.second;
		return a;
	}//operator+=(...)

}//namespace VMlib

using VMlib::PairInt;

#endif
 