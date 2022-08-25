/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.11   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2022/08/07     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: PointsCopy2D.h                                                   |
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
\brief Заголовочный файл с описанием класса PointsCopy
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.11
\date 07 августа 2022 г.
*/


#ifndef POINTSCOPY_H_
#define POINTSCOPY_H_

#include "Point2D.h"
#include "Vortex2D.h"
#include "defs.h"

namespace VM2D
{
	/*!
	\brief Класс, опеделяющий копию структуры данных следа (пелены), источников или "виртуальных" вихрей или источников
	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.11
	\date 07 августа 2022 г.
	*/
	class PointsCopy : public Vortex2D
	{
	public:
		Point2D veloCopy;
		PointType type;
		double domainRadius;

		///Для панелей
		/*
		Point2D veloCopyLin;
		std::vector<Point2D> i00, i01, i10, i11;

		Point2D a, c; //для предобуславливателя в Т0
		Point2D a1, c1; //для предобуславливателя в Т1

		Point2D panBegin, panEnd;
		/// Мультиполные моменты панелей
		Point2D dipP, quaP, octP, hexP;

		Point2D tau;
		double len;
		double gamLin;
		///\brief Список из пар значений - номер рассматриваемого профиля и номер панели на нем
		std::vector<std::pair<int, int>> aflPnl;
		*/

		PointsCopy(const Vortex2D& vtx_) :
			Vortex2D(vtx_), veloCopy({ 0.0, 0.0 })/*, veloCopyLin({ 0.0, 0.0 })*/
		{};
		~PointsCopy() {};
	};
}


using VM2D::PointsCopy;
#endif