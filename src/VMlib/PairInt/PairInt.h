/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: PairInt.h                                                        |
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
\brief Заголовочный файл с описанием класса pairInt
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef PAIRINT_H_
#define PAIRINT_H_

#include "numvector.h"

namespace VMlib
{
	/*!
	\brief Класс, опеделяющий пару целых чисел

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна
	\Version 1.14
	\date 6 марта 2026 г.
	*/
	class PairInt : public numvector<int, 2>
	{
	public:
		/// Пустой конструктор
		PairInt() { };

		/// \brief Конструктор и приведение типа из numvector<int, 2>
		///
		/// \param[in] _r константная ссылка на копируемый объект типа numvector<int, 2>
		PairInt(const numvector<int, 2>& _r) { data[0] = _r[0];data[1] = _r[1]; };

		/// \brief Конструктор копирования
		///
		/// \param[in] _r константная ссылка на копируемый вектор
		PairInt(const PairInt& _r) { data[0] = _r[0];data[1] = _r[1]; };

#if !defined(__CUDACC__)
		/// \brief Конструктор инициализации списком
		///
		/// \param[in] z константная ссылка на список инициализации из чисел типа int
		/// \warning Длина списка инициализации не проверяется, от него берутся только 2 первых элемента
		PairInt(const std::initializer_list<int>& z) {
			for (size_t i = 0; i < 2; ++i)
				data[i] = *(z.begin() + i);
		};
#endif

		/// Деструктор
		~PairInt() { };
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
 