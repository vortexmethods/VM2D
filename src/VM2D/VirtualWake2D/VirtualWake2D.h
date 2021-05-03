/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: VirtualWake2D.h                                                  |
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
\brief Заголовочный файл с описанием класса VirtualWake
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
*/

#ifndef VIRTUALWAKE_H
#define VIRTUALWAKE_H

#include "WakeDataBase2D.h"

namespace VM2D
{

	class World2D;
	class Boundary;

	/*!
	\brief Класс, опеделяющий вихревой след (пелену)
	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.10
	\date 17 мая 2021 г.
	*/
	class VirtualWake : public WakeDataBase
	{
	private:

	public:

		/// Константная ссылка на границу, которой принадлежит данный виртуальный след
		const Boundary& bnd;

		/// Скорость вихрей виртуального следа конкретного профиля (равна Gamma/2) используется для расчета давления
		std::vector<Point2D> vecHalfGamma;

		/// Пара чисел: номер профиля и номер панели, на которой рожден виртуальный вихрь
		std::vector<std::pair<size_t, size_t>> aflPan;

		/// \brief Конструктор инициализации
		///
		/// \param[in] W_ константная ссылка на решаемую задачу	
		/// \param[in] bnd_ константная ссылка на границу	
		VirtualWake(const World2D& W_, const Boundary& bnd_)
			: WakeDataBase(W_), bnd(bnd_) { };

		/// Деструктор
		~VirtualWake() { };

	};

} //namespace VM2D

#endif
