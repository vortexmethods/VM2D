/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MechanicsRigidOscillPart.cpp                                     |
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
\brief Файл кода с описанием класса MechanicsRigidOscillPart
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/

#include "MechanicsRigidOscillPart.h"


//Вычисление гидродинамической силы, действующей на профиль
void MechanicsRigidOscillPart::GetHydroDynamForce(timePeriod& time)
{
	time.first = omp_get_wtime();

	hydroDynamForce = { 0.0, 0.0 };

	Point2D hDFGam = { 0.0, 0.0 };	//гидродинамические силы, обусловленные Gamma_k
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные delta_k

	for (size_t i = 0; i < afl.np; ++i)
	{
		double GamK = boundary.virtualWake[i].g();
		double deltaK = boundary.sheets.freeVortexSheet[i][0] * afl.len[i] - afl.gammaThrough[i];  		 //afl.gammaThrough[i];
		Point2D VelK = virtVortParams.convVelo[i] + passport.physicalProperties.V0();
		Point2D rK = 0.5 * (afl.r[i + 1] + afl.r[i]);	//boundary.virtualWake[i].r();		

		hDFGam += GamK * Point2D({ VelK[1], -VelK[0] });
		hDFdelta += deltaK * Point2D({ -rK[1], rK[0] });
	}

	hydroDynamForce = hDFGam + (1.0 / passport.timeDiscretizationProperties.dt) * hDFdelta;

	time.second = omp_get_wtime();
}// GetHydroDynamForce(...)

// Вычисление скорости центра масс
Point2D MechanicsRigidOscillPart::VeloOfAirfoilRcm(double currTime)
{
	return {0.0, u};
	//return{ 0.0, 0.0 };
}//VeloOfAirfoilRcm(...)

// Вычисление положения центра масс
Point2D MechanicsRigidOscillPart::PositionOfAirfoilRcm(double currTime)
{
	return{ 0.0, y };
	//return{ 0.0, 0.0 };
}//PositionOfAirfoilRcm(...)

// Вычисление скоростей начал панелей
void MechanicsRigidOscillPart::VeloOfAirfoilPanels(double currTime)
{
	Point2D veloRcm = VeloOfAirfoilRcm(currTime);
	for (size_t i = 0; i < afl.v.size(); ++i)
		afl.v[i] = veloRcm;
		//afl.v[i] = { 0.0, 0.0 };
}//VeloOfAirfoilPanels(...)


void MechanicsRigidOscillPart::Move()
{
	
	double step = passport.timeDiscretizationProperties.dt;
	double dy, du;
	numvector<double, 2> kk[4];
	kk[0] = { u, (hydroDynamForce[1]-b*u-k*y)/m};
	kk[1] = { u + 0.5*step*kk[0][1], (hydroDynamForce[1] - b*(u + 0.5*step*kk[0][1]) - k*(y + 0.5*step*kk[0][0])) / m };
	kk[2] = { u + 0.5*step*kk[1][1], (hydroDynamForce[1] - b*(u + 0.5*step*kk[1][1]) - k*(y + 0.5*step*kk[1][0])) / m };
	kk[3] = { u + step*kk[2][1], (hydroDynamForce[1] - b*(u + step*kk[2][1]) - k*(y + step*kk[2][0])) / m };
	
	dy = step * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
	du = step * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	
	//afl.Move({0.0, 0.0});
	afl.Move({ 0.0, dy });

	//y += step*0.0;
	//u = 0.0;

	y += dy;
	u += du;

	//std::cout << "du = " << du << ", dy = " << dy << std::endl;
}//Move()

