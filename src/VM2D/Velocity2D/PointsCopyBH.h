/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: PointsCopyBH.h                                                   |
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
\brief Заголовок класса-обертки для точек и панелей для метода Барнса - Хата для CPU
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/


#ifndef POINTSCOPY_H
#define POINTSCOPY_H

#include "Vortex2D.h"
#include "ParamsBH.h"
#include "omp.h"

namespace BH
{
/*!
\brief Класс-обертка для точек и панелей
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

	class PointsCopy : public Vortex2D
	{
	public:

		/// Вычисленная скорость
		Point2D veloCopy; 

		/// Вычисленное расстояние до ближайших вихрей
		double epsast;
		
		/// Инициализирующий конструктор
		PointsCopy() :
			Vortex2D({ 0.0, 0.0 }, 0.0)			
		{};

		/// Конструктор копирования из точки
		PointsCopy(const Point2D& r) :
			Vortex2D(r, 0.0)
		{};

		/// Конструктор копирования из вихря
		PointsCopy(const Vortex2D& vtx_) :
			Vortex2D(vtx_), veloCopy({ 0.0, 0.0 }), epsast(0)
		{};

		/// Деструктор
		~PointsCopy() {};
	};

}//namespace BH

#endif
