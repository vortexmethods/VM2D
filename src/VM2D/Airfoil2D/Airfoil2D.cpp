/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.9    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/07/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
\brief Файл кода с описанием класса Airfoil
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.9   
\date 22 июля 2020 г.
*/

#include "Airfoil2D.h"

#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Tree2D.h"
#include "Velocity2D.h"
#include "WakeDataBase2D.h"
#include "World2D.h"
#include "Wake2D.h"

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

//Вычисление средних значений eps на панелях
void Airfoil::calcMeanEpsOverPanel()
{	
	meanEpsOverPanel.clear();
	meanEpsOverPanel.resize(getNumberOfPanels());

	if (W.getParallel().myidWork == 0)
	{
		double midEps;

		const Boundary& bnd = W.getBoundary(numberInPassport);
		VortexesParams virtVortParams = W.getVelocity().virtualVortexesParams[numberInPassport];
				
		for (size_t i = 0; i < getNumberOfPanels(); ++i)
		{
			midEps = 0.0;

			for (int j = bnd.vortexBeginEnd[i].first; j < bnd.vortexBeginEnd[i].second; ++j)
				midEps += virtVortParams.epsastWake[j];
			
			midEps /= (bnd.vortexBeginEnd[i].second - bnd.vortexBeginEnd[i].first);

			meanEpsOverPanel[i] = midEps;
		}
	}//if id==0

	MPI_Bcast(meanEpsOverPanel.data(), (int)meanEpsOverPanel.size(), MPI_DOUBLE, 0, W.getParallel().commWork);
}