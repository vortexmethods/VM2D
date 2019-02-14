/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.5    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/02/20     |
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
\version 1.5   
\date 20 февраля 2019 г.
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
	: np(0), W(W_), numberInPassport(numberInPassport_)
{ }


// Вычисление нормалей, касательных и длин панелей по текущему положению вершин
void Airfoil::CalcNrmTauLen()		
{
	if (nrm.size() != np)
	{
		nrm.resize(np);
		tau.resize(np);
		len.resize(np);
	}

	Point2D rpan;

	for (size_t i = 0; i < np; ++i)
	{
		rpan = (r[i + 1] - r[i]);
		len[i] = rpan.length();
		tau[i] = rpan.unit();
			
		nrm[i] = { tau[i][1], -tau[i][0] };
	}

	//закольцовываем
	nrm.push_back(nrm[0]);
	tau.push_back(tau[0]);
	len.push_back(len[0]);
}//CalcNrmTauLen()


//Проверка, идет ли вершина i следом за вершиной j
bool Airfoil::isAfter(size_t i, size_t j) const
{
	return ((i == j + 1) || (i == 0 && j == np - 1));
}//isAfter(...)


//Перемещение профиля 
void Airfoil::Move(const Point2D& dr)	//перемещение профиля как единого целого на вектор dr
{
	for (size_t i = 0; i < np + 1; i++)
		r[i] += dr;
	rcm += dr;
	CalcNrmTauLen();
	GetGabarits();
}//Move(...)


//Поворот профиля 
void Airfoil::Rotate(double alpha)	//поворот профиля на угол alpha вокруг центра масс
{
	phiAfl += alpha;
	numvector<numvector<double, 2>, 2> rotMatrix = { { cos(alpha), sin(alpha) }, { -sin(alpha), cos(alpha) } };

	for (size_t i = 0; i < np + 1; i++)
		r[i] = rcm + dot(rotMatrix, r[i] - rcm);

	CalcNrmTauLen();
	GetGabarits();
}//Rotate(...)


//Масштабирование профиля
void Airfoil::Scale(double factor)	//масштабирование профиля на коэффициент factor относительно центра масс
{
	for (size_t i = 0; i < np + 1; i++)
		r[i] = rcm + factor*(r[i] - rcm);

	CalcNrmTauLen(); //строго говоря, меняются при этом только длины, нормали и касательные - без изменения
	GetGabarits();
}//Scale(...)


//Вычисляет габаритный прямоугольник профиля
void Airfoil::GetGabarits(double gap)	//определение габаритного прямоугольника
{
	lowLeft = {  1E+10,  1E+10 };
	upRight = { -1E+10, -1E+10 };

	for (size_t i = 0; i < np + 1; i++)
	{
		lowLeft[0] = std::min(lowLeft[0], r[i][0]);
		lowLeft[1] = std::min(lowLeft[1], r[i][1]);

		upRight[0] = std::max(upRight[0], r[i][0]);
		upRight[1] = std::max(upRight[1], r[i][1]);
	}

	Point2D size = upRight - lowLeft;
	lowLeft -= size*gap;
	upRight += size*gap;
}//GetGabarits(...)


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

