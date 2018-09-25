/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.2    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/06/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MechanicsRigidOscillStronglyCoupled.cpp                          |
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
\brief Файл кода с описанием класса MechanicsRigidOscillStronglyCoupled
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.2
\date 14 июня 2018 г.
*/

#include "mpi.h"

#include "MechanicsRigidOscillStronglyCoupled.h"
#include "World2D.h"

//Вычисление гидродинамической силы, действующей на профиль
void MechanicsRigidOscillStronglyCoupled::GetHydroDynamForce(timePeriod& time)
{
	time.first = omp_get_wtime();

	const double& dt = W.getPassport().timeDiscretizationProperties.dt;

	hydroDynamForce = { 0.0, 0.0 };
	hydroDynamMoment = 0.0;

	Point2D hDFGam = { 0.0, 0.0 };	//гидродинамические силы, обусловленные Gamma_k
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные delta_k
	double hDMdelta = 0.0;


	Point2D deltaVstep = { 0.0, Vcm[1] - VcmOld[1] };

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
Point2D MechanicsRigidOscillStronglyCoupled::VeloOfAirfoilRcm(double currTime)
{
	return Vcm;
	//return{ 0.0, 0.0 };
}//VeloOfAirfoilRcm(...)

// Вычисление положения центра масс
Point2D MechanicsRigidOscillStronglyCoupled::PositionOfAirfoilRcm(double currTime)
{
	return Rcm;
	//return{ 0.0, 0.0 };
}//PositionOfAirfoilRcm(...)

// Вычисление скоростей начал панелей
void MechanicsRigidOscillStronglyCoupled::VeloOfAirfoilPanels(double currTime)
{
	Point2D veloRcm = VeloOfAirfoilRcm(currTime);
	for (size_t i = 0; i < afl.v.size(); ++i)
		afl.v[i] = veloRcm;
	//afl.v[i] = { 0.0, 0.0 };
}//VeloOfAirfoilPanels(...)


void MechanicsRigidOscillStronglyCoupled::Move()
{
	VcmOld = Vcm;
	
	Point2D dRcm = Vcm * W.getPassport().timeDiscretizationProperties.dt;
	Rcm += dRcm;
	afl.Move(dRcm);
}//Move()

void MechanicsRigidOscillStronglyCoupled::ReadSpecificParametersFromDictionary()
{
	mechParamsParser->get("m", m);
	
	std::cout << "m = " << m << std::endl;
	std::vector<double> sh;
	mechParamsParser->get("sh", sh);

	k =  m * 4.0 * PI * PI * sh[1] * sh[1] * W.getPassport().physicalProperties.vInf.length2();

//	std::cout << "constr k = " << k << std::endl;
}

/// \todo пока считается, что вихревой слой состоит из отдельных изолированных вихрей
void MechanicsRigidOscillStronglyCoupled::FillMechanicsRowsAndCross(Eigen::MatrixXd& row, Eigen::MatrixXd& cross)
{
	cross(0, 0) = -(/*3*/ m + /*4*/ W.getPassport().timeDiscretizationProperties.dt * b);

	Point2D riC;

	for (size_t loc_i = 0; loc_i < afl.np; ++loc_i)
	{
		riC = 0.5 * (afl.r[loc_i + 1] + afl.r[loc_i]) - afl.rcm;
		row(0, loc_i) = /*1*/ riC.kcross()[1] * afl.len[loc_i];
		cross(0, 0) += /*2*/ riC.kcross()[1] * afl.tau[loc_i][1] * afl.len[loc_i];
	}
}
//FillMechanicsRowsAndCols(...)


/// \todo пока считается, что вихревой слой состоит из отдельных изолированных вихрей
void MechanicsRigidOscillStronglyCoupled::FillMechanicsRhs(std::vector<double>& rhs)
{
	rhs[0] = /*7*/-m * Vcm[1] + /*8*/ W.getPassport().timeDiscretizationProperties.dt * k * Rcm[1];

	for (size_t loc_i = 0; loc_i < afl.np; ++loc_i)
	{
		Point2D riC = 0.5 * (afl.r[loc_i + 1] + afl.r[loc_i]) - afl.rcm;
		rhs[0] += riC.kcross()[1] * (/*5*/afl.gammaThrough[loc_i] + /*6*/Vcm[1] * afl.tau[loc_i][1] * afl.len[loc_i]);
	}
}

void MechanicsRigidOscillStronglyCoupled::FillAtt(Eigen::MatrixXd& col, Eigen::MatrixXd& rhs)
{
	for (size_t i = 0; i < afl.np; ++i)
		col(i, 0) = /*12*/ -afl.tau[i][1];
}

void MechanicsRigidOscillStronglyCoupled::SolutionToMechanicalSystem(Eigen::VectorXd& col)
{
	VcmOld = Vcm;

	Vcm[1] = col(0);

	for (size_t i = 0; i < afl.np + 1; ++i)
		afl.v[i] = Vcm;
}