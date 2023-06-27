/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DRigidImmovable.cpp                                    |
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
\brief Файл кода с описанием класса MechanicsRigidImmovable
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#include "Mechanics2DRigidImmovable.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;


MechanicsRigidImmovable::MechanicsRigidImmovable(const World2D & W_, size_t numberInPassport_)
	: Mechanics(W_, numberInPassport_, false, false)
{
	ReadSpecificParametersFromDictionary();
	Initialize({ 0.0, 0.0 }, W_.getAirfoil(numberInPassport_).rcm, 0.0, W_.getAirfoil(numberInPassport_).phiAfl);
};

//Вычисление скорости центра масс профиля
Point2D MechanicsRigidImmovable::VeloOfAirfoilRcm(double currTime)
{
	return{ 0.0, 0.0 };
}//VeloOfAirfoilRcm(...)

//Вычисление положения центра масс профиля
Point2D MechanicsRigidImmovable::PositionOfAirfoilRcm(double currTime)
{
	return afl.rcm;
}//PositionOfAirfoilRcm(...)

double MechanicsRigidImmovable::AngularVelocityOfAirfoil(double currTime)
{
	return 0.0;
}//AngularVelocityOfAirfoil(...)

//Вычисление угла поворота профиля
double MechanicsRigidImmovable::AngleOfAirfoil(double currTime)
{
	return afl.phiAfl;
}//AngleOfAirfoil(...)

// Вычисление скоростей начал панелей
void MechanicsRigidImmovable::VeloOfAirfoilPanels(double currTime)
{
	Point2D veloRcm = VeloOfAirfoilRcm(currTime);
	afl.setV(veloRcm);	
}//VeloOfAirfoilPanels(...)

//Вычисление гидродинамической силы, действующей на профиль
void MechanicsRigidImmovable::GetHydroDynamForce()
{
	W.getTimestat().timeGetHydroDynamForce.first += omp_get_wtime();

	const double& dt = W.getPassport().timeDiscretizationProperties.dt;

	hydroDynamForce = { 0.0, 0.0 };
	hydroDynamMoment = 0.0;

	viscousForce = { 0.0, 0.0 };
	viscousMoment = 0.0;

	Point2D hDFGam = { 0.0, 0.0 };	//гидродинамические силы, обусловленные Gamma_k
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные delta_k
	double hDMdelta = 0.0;

	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		/// \todo Учитываем только нулевой момент решения. Надо ли учитывать остальные?
		double deltaK =  boundary.sheets.freeVortexSheet(i, 0) * afl.len[i] - afl.gammaThrough[i];
		Point2D rK = 0.5 * (afl.getR(i + 1) + afl.getR(i)) - afl.rcm;

		hDFdelta += deltaK * Point2D({ -rK[1], rK[0] });
		hDMdelta += 0.5 * deltaK * rK.length2();	
	}

	const double rho = W.getPassport().physicalProperties.rho;

	hydroDynamForce = (rho / dt) * hDFdelta;
	hydroDynamMoment = (rho / dt) * hDMdelta;
		
	if (W.getPassport().physicalProperties.nu > 0.0)
		for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
		{
			Point2D rK = 0.5 * (afl.getR(i + 1) + afl.getR(i)) - afl.rcm;
			viscousForce += rho * afl.viscousStress[i] * afl.tau[i];
			viscousMoment += rho * (afl.viscousStress[i] * afl.tau[i]) & rK;
		}

	W.getTimestat().timeGetHydroDynamForce.second += omp_get_wtime();
}//GetHydroDynamForce()
