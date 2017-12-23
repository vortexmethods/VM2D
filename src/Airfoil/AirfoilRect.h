/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2017/12/01     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina       |
*-----------------------------------------------------------------------------*
| File name: AirfoilRect.h                                                    |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.	                                      |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.	                                                          |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Заголовочный файл с описанием класса AirfoilRect
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef AIRFOILRECT_H
#define AIRFOILRECT_H

#include "Airfoil.h"

/*!
\brief Класс, определяющий тип обтекаемого профиля

Тип профиля:
- профиль с прямолинейными панелями.

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.0
\date 1 декабря 2017 г.
*/
class AirfoilRect 
	: public Airfoil
{
public:

	/// Конструктор
	AirfoilRect(const Passport& passport_, const int numberInPassport_, const Parallel& parallel_)
		:Airfoil(passport_, numberInPassport_, parallel_)
	{ };

	AirfoilRect(const Airfoil& afl) : Airfoil(afl){};

	/// Деструктор
	virtual ~AirfoilRect() { };

	//далее -- реализация виртуальной функции
	virtual void ReadFromFile(const std::string& dir);

	virtual void GetDiffVelocityToSetOfPointsAndViscousStresses(const std::vector<Vortex2D>& points, std::vector<double>& domainRadius, std::vector<Point2D>& velo);
};

#endif