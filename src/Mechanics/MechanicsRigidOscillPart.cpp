/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.2    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/06/14     |
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
\version 1.2
\date 14 июня 2018 г.
*/

#include "mpi.h"

#include "MechanicsRigidOscillPart.h"
#include "World2D.h"

//Вычисление гидродинамической силы, действующей на профиль
void MechanicsRigidOscillPart::GetHydroDynamForce(timePeriod& time)
{
	time.first = omp_get_wtime();

	const double& dt = W.getPassport().timeDiscretizationProperties.dt;

	hydroDynamForce = { 0.0, 0.0 };
	hydroDynamMoment = 0.0;

	Point2D hDFGam = { 0.0, 0.0 };	//гидродинамические силы, обусловленные Gamma_k
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные delta_k
	double hDMdelta = 0.0;


	Point2D deltaVstep = { 0.0, u - uOld };

	for (size_t i = 0; i < afl.np; ++i)
	{

		double deltaK = boundary.sheets.freeVortexSheet[i][0] * afl.len[i] - afl.gammaThrough[i] + deltaVstep * afl.tau[i] * afl.len[i];
		Point2D rK = 0.5 * (afl.r[i + 1] + afl.r[i]) - afl.rcm;

		hDFdelta += deltaK * Point2D({ -rK[1], rK[0] });
		hDMdelta += 0.5 * deltaK * (rK * rK);
	}

	hydroDynamForce = hDFdelta * (1.0 / dt);
	hydroDynamMoment = hDMdelta / dt;

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
	uOld = u;
	yOld = y;
	
	double dy, du;

	std::cout << "k = " << k << std::endl;

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
		
	y += dy;
	u += du;	
}//Move()


void MechanicsRigidOscillPart::ReadSpecificParametersFromDictionary()
{
	mechParamsParser->get("m", m);
	
	std::cout << "m = " << m << std::endl;
		
	std::vector<double> sh;
	mechParamsParser->get("sh", sh);
	
	k = m * 4.0 * PI * PI * sh[1] * sh[1] * W.getPassport().physicalProperties.vInf.length2();

	std::cout << "constr k = " << k << std::endl;
}

void MechanicsRigidOscillPart::FillAtt(Eigen::MatrixXd& col, Eigen::MatrixXd& rhs)
{
	for (int i = 0; i < afl.np; ++i)
		rhs(i, 0) = /*12*/ afl.v[i] * afl.tau[i];
}
