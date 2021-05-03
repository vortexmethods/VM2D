/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DRigidOscillPart.cpp                                   |
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
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
*/

#include "mpi.h"

#include "Mechanics2DRigidOscillPart.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Tree2D.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

//Вычисление гидродинамической силы, действующей на профиль
void MechanicsRigidOscillPart::GetHydroDynamForce()
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


	Point2D deltaVstep = { 0.0, u - uOld };

	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		/// \todo Учитываем только нулевой момент решения. Надо ли учитывать остальные?
		double deltaK = boundary.sheets.freeVortexSheet(i, 0) * afl.len[i] - afl.gammaThrough[i] + (deltaVstep & afl.tau[i]) * afl.len[i];
		Point2D rK = 0.5 * (afl.getR(i + 1) + afl.getR(i)) - afl.rcm;

		hDFdelta += deltaK * Point2D({ -rK[1], rK[0] });
		hDMdelta += 0.5 * deltaK * rK.length2();
	}

	hydroDynamForce = hDFdelta * (1.0 / dt);
	hydroDynamMoment = hDMdelta / dt;

	if (W.getPassport().physicalProperties.nu > 0.0)
		for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
		{
			Point2D rK = 0.5 * (afl.getR(i + 1) + afl.getR(i)) - afl.rcm;
			viscousForce += afl.viscousStress[i] * afl.tau[i];
			viscousMoment += (afl.viscousStress[i] * afl.tau[i]) & rK;
		}
	
	W.getTimestat().timeGetHydroDynamForce.second += omp_get_wtime();
}// GetHydroDynamForce()

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
	afl.setV(veloRcm);	
}//VeloOfAirfoilPanels(...)


void MechanicsRigidOscillPart::Move()
{
	uOld = u;
	yOld = y;
	
	double dy, du;

	//W.getInfo('t') << "k = " << k << std::endl;

	if (W.getParallel().myidWork == 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;		
		numvector<double, 2> kk[4];
		kk[0] = { u, (hydroDynamForce[1] - b*u - k*y) / m };
		kk[1] = { u + 0.5*dt*kk[0][1], (hydroDynamForce[1] - b*(u + 0.5*dt*kk[0][1]) - k*(y + 0.5*dt*kk[0][0])) / m };
		kk[2] = { u + 0.5*dt*kk[1][1], (hydroDynamForce[1] - b*(u + 0.5*dt*kk[1][1]) - k*(y + 0.5*dt*kk[1][0])) / m };
		kk[3] = { u + dt*kk[2][1], (hydroDynamForce[1] - b*(u + dt*kk[2][1]) - k*(y + dt*kk[2][0])) / m };

		dy = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		du = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}

	MPI_Bcast(&dy, 1, MPI_DOUBLE, 0, W.getParallel().commWork);
	MPI_Bcast(&du, 1, MPI_DOUBLE, 0, W.getParallel().commWork);

	afl.Move({ 0.0, dy });
	Vcm[1] += du; 
		
	y += dy;
	u += du;	
}//Move()


#ifdef BRIDGE
void MechanicsRigidOscillPart::ReadSpecificParametersFromDictionary()
{
        mechParamsParser->get("m", m);

        W.getInfo('i') << "mass " << "m = " << m << std::endl;

        mechParamsParser->get("b", b);

        W.getInfo('i') << "damping " << "b = " << b << std::endl;


        mechParamsParser->get("k", k);

//      std::vector<double> sh;
//      mechParamsParser->get("sh", sh);

//      k = m * 4.0 * PI * PI * sh[1] * sh[1] * W.getPassport().physicalProperties.vInf.length2();

        W.getInfo('i') << "rigidity k = " << k << std::endl;
}//ReadSpecificParametersFromDictionary()
#endif

#ifdef INITIAL
void MechanicsRigidOscillPart::ReadSpecificParametersFromDictionary()
{
	mechParamsParser->get("m", m);
	
	W.getInfo('i') << "mass " << "m = " << m << std::endl;
		
	std::vector<double> sh;
	mechParamsParser->get("sh", sh);
	
	k = m * 4.0 * PI * PI * sh[1] * sh[1] * W.getPassport().physicalProperties.vInf.length2();

	W.getInfo('i') << "rigidity k = " << k << std::endl;
}//ReadSpecificParametersFromDictionary()
#endif