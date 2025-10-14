/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.5    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2024/06/19     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2024 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: PointsCopy.h                                                     |
| Info: Source code of BH                                                     |
|                                                                             |
| This file is part of BH.                                                    |
| BH is free software: you can redistribute it and/or modify it               |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| BHcu is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with BH.  If not, see <http://www.gnu.org/licenses/>.                 |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Заголовок класса-обертки для точек и панелей
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 19 июня 2024 г.
*/


#ifndef POINTSCOPY_H_
#define POINTSCOPY_H_

#include "Vortex2D.h"
#include "Params.h"
#include "omp.h"

namespace BH
{
/*!
\brief Класс-обертка для точек и панелей
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 19 июня 2024 г.
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
