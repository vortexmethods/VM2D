/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Airfoil2D.cpp                                                    |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VM is distributed in the hope that it will be useful, but WITHOUT           |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Файл кода с описанием класса Airfoil
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7   
\date 22 ноября 2019 г.
*/

#include "Airfoil2D.h"

#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "WakeDataBase2D.h"

using namespace VM2D;

// Конструктор
Airfoil::Airfoil(const World2D& W_, const size_t numberInPassport_)
	: W(W_), numberInPassport(numberInPassport_)
{ }



//Проверка, идет ли вершина i следом за вершиной j
bool Airfoil::isAfter(size_t i, size_t j) const
{
	return ((i == j + 1) || (i == 0 && j == r_.size() - 1));
}//isAfter(...)







//Определяет, находится ли точка с радиус-вектором r внутри габаритного прямоугольника профиля
bool Airfoil::isInsideGabarits(const Point2D& r) const
{
	return (r[0] <= upRight[0] && (r[0] >= lowLeft[0] && r[1] >= lowLeft[1] && r[1] <= upRight[1]));
}//isInsideGabarits(...)


//Определяет, находится ли точка с радиус-вектором r вне габаритного прямоугольника профиля
bool Airfoil::isOutsideGabarits(const Point2D& r) const
{
	return (r[0] > upRight[0] || (r[0] < lowLeft[0] || r[1] < lowLeft[1] || r[1] > upRight[1]));
}//isOutsideGabarits(...)

