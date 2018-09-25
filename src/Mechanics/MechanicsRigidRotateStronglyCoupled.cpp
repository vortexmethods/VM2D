/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.3    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/09/26     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MechanicsRigidRotateStronglyCoupled.cpp                          |
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
\brief Файл кода с описанием класса MechanicsRigidRotateStronglyCoupled
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/

#include "mpi.h"

#include "MechanicsRigidRotateStronglyCoupled.h"
#include "World2D.h"

//Вычисление гидродинамической силы, действующей на профиль
void MechanicsRigidRotateStronglyCoupled::GetHydroDynamForce(timePeriod& time)
{
	time.first = omp_get_wtime();

	const double& dt = W.getPassport().timeDiscretizationProperties.dt;

	hydroDynamForce = { 0.0, 0.0 };
	hydroDynamMoment = 0.0;

	Point2D hDFGam = { 0.0, 0.0 };		//гидродинамические силы, обусловленные Gamma_k
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные delta_k
	Point2D hDFQ = { 0.0, 0.0 };		//гидродинамические силы, обусловленные Q_k
	double hDMGam = 0.0;				//гидродинамический момент, обусловленный Gamma_k
	double hDMdelta = 0.0;				//гидродинамический момент, обусловленный Gamma_k
	double hDMQ = 0.0;					//гидродинамический момент, обусловленный Gamma_k

	for (size_t i = 0; i < afl.np; ++i)
	{
		Point2D rK = 0.5 * (afl.r[i + 1] + afl.r[i]) - afl.rcm;
		Point2D velK = { -Wcm * rK[1], Wcm * rK[0] };

		double gAtt = Wcm * (rK ^ afl.tau[i]);
		double gAttOld = WcmOld * (rK ^ afl.tau[i]);
		double deltaGAtt = gAtt - gAttOld;

		double qAtt = Wcm * (rK ^ afl.nrm[i]);
		
		/*1*/
		double deltaK = boundary.sheets.freeVortexSheet[i][0] * afl.len[i] - afl.gammaThrough[i] + deltaGAtt * afl.len[i];
		hDFdelta += deltaK * Point2D({ -rK[1], rK[0] });
		hDMdelta += 0.5 * deltaK * (rK * rK);
		
		/*2*/
		hDFGam += 0.5 * afl.v[i].kcross() * gAtt * afl.len[i];
		hDMGam += 0.5 * rK ^ afl.v[i].kcross() * gAtt * afl.len[i];

		/*3*/
		hDFQ -= 0.5 * afl.v[i] * qAtt * afl.len[i];
		hDMQ -= 0.5 * rK ^ afl.v[i] * qAtt * afl.len[i];
	}

	hydroDynamForce = hDFGam + hDFdelta * (1.0 / dt) + hDFQ;
	hydroDynamMoment = hDMGam + hDMdelta / dt + hDMQ;

	time.second = omp_get_wtime();
}// CalcHydroDynamForce(...)

// Вычисление скорости центра масс
Point2D MechanicsRigidRotateStronglyCoupled::VeloOfAirfoilRcm(double currTime)
{
	return{ 0.0, 0.0 };
	//return{ 0.0, 0.0 };
}//VeloOfAirfoilRcm(...)

// Вычисление положения центра масс
Point2D MechanicsRigidRotateStronglyCoupled::PositionOfAirfoilRcm(double currTime)
{
	return afl.rcm;
}//PositionOfAirfoilRcm(...)

// Вычисление скоростей начал панелей
void MechanicsRigidRotateStronglyCoupled::VeloOfAirfoilPanels(double currTime)
{
	for (size_t i = 0; i < afl.v.size(); ++i)
	{
		afl.v[i][0] = afl.r[i].kcross()[0] * Wcm;
		afl.v[i][1] = afl.r[i].kcross()[1] * Wcm;
	}
}//VeloOfAirfoilPanels(...)


void MechanicsRigidRotateStronglyCoupled::Move()
{
	PhiOld = Phi;

	double dphi = Wcm * W.getPassport().timeDiscretizationProperties.dt;
	MPI_Bcast(&dphi, 1, MPI_DOUBLE, 0, W.getParallel().commWork);

	Phi += dphi;
	afl.Rotate(dphi);
}//Move()

void MechanicsRigidRotateStronglyCoupled::ReadSpecificParametersFromDictionary()
{
	mechParamsParser->get("J", J);
	mechParamsParser->get("k", k);
	mechParamsParser->get("tRotateAccel", tRotateAccel);
	mechParamsParser->get("tMomentAccel", tMomentAccel);
	mechParamsParser->get("Mz", Mz);

	std::cout << "J = " << J << std::endl;
	std::cout << "k = " << k << std::endl;
	std::cout << "tRotateAccel = " << tRotateAccel << std::endl;
	std::cout << "tMomentAccel = " << tMomentAccel << std::endl;
	std::cout << "Mz = " << Mz << std::endl;

//	std::vector<double> sh;
//	mechParamsParser->get("sh", sh);
//	k = 0.0 * m * 4.0 * PI * PI * sh[1] * sh[1] * W.getPassport().physicalProperties.V0().length2();
}

/// \todo пока считается, что вихревой слой состоит из отдельных изолированных вихрей
void MechanicsRigidRotateStronglyCoupled::FillMechanicsRowsAndCross(Eigen::MatrixXd& row, Eigen::MatrixXd& cross)
{
	cross(0, 0) = -(/*3*/ J + /*4*/ W.getPassport().timeDiscretizationProperties.dt * b);

	Point2D riC;

	for (size_t loc_i = 0; loc_i < afl.np; ++loc_i)
	{
		riC = 0.5 * (afl.r[loc_i + 1] + afl.r[loc_i]) - afl.rcm;

		row(0, loc_i) = /*1*/ 0.5 * riC * riC * afl.len[loc_i];
		cross(0, 0) += /*2*/ 0.5 * riC * riC * (riC ^ afl.tau[loc_i]) * afl.len[loc_i];
	}
}
//FillMechanicsRowsAndCols(...)


/// \todo пока считается, что вихревой слой состоит из отдельных изолированных вихрей
void MechanicsRigidRotateStronglyCoupled::FillMechanicsRhs(std::vector<double>& rhs)
{
	rhs[0] = /*7*/-J * Wcm + /*8*/ W.getPassport().timeDiscretizationProperties.dt * k * Phi -\
				/*11*/GetMz(W.getCurrentStep() * W.getPassport().timeDiscretizationProperties.dt)* W.getPassport().timeDiscretizationProperties.dt;

	for (size_t loc_i = 0; loc_i < afl.np; ++loc_i)
	{
		Point2D riC = 0.5 * (afl.r[loc_i + 1] + afl.r[loc_i]) - afl.rcm;

		rhs[0] += /*5*/ 0.5 * riC * riC * afl.gammaThrough[loc_i];
		rhs[0] += /*6*/ 0.5 * riC * riC * Wcm * (riC ^ afl.tau[loc_i]) * afl.len[loc_i];
		rhs[0] += /*9*/ 0.5 * afl.len[loc_i] * (Wcm * riC ^ afl.tau[loc_i]) * (-Wcm * riC) ^ riC * W.getPassport().timeDiscretizationProperties.dt;
		rhs[0] -= /*10*/ -0.5 * Wcm * riC * riC * afl.len[loc_i] * (riC ^ afl.nrm[loc_i]) * Wcm * W.getPassport().timeDiscretizationProperties.dt;
	}
}

void MechanicsRigidRotateStronglyCoupled::FillAtt(Eigen::MatrixXd& col, Eigen::MatrixXd& rhs)
{
	Point2D riC;
	for (size_t i = 0; i < afl.np; ++i)
	{
		riC = 0.5 * (afl.r[i + 1] + afl.r[i]) - afl.rcm;

		col(i, 0) = /*12 + 91*/ -0.5 * (riC ^ afl.tau[i]);
		for (size_t j = 0; j < afl.np; ++j)
		{
			col(i, 0) += /*92*/ IDPI * (Alpha(riC - afl.r[j + 1], riC - afl.r[j]) * afl.tau[j]
									 - Lambda(riC - afl.r[j + 1], riC - afl.r[j]) * afl.nrm[j]
									   )* afl.tau[i] * (riC ^ afl.tau[j]);

			col(i, 0) += /*93*/ IDPI * (Alpha(riC - afl.r[j + 1], riC - afl.r[j]) * afl.nrm[j] +
									   Lambda(riC - afl.r[j + 1], riC - afl.r[j]) * afl.tau[j]
									   )* afl.tau[i] * (riC ^ afl.nrm[j]);
		}
	}
		
}

void MechanicsRigidRotateStronglyCoupled::SolutionToMechanicalSystem(Eigen::VectorXd& col)
{
	WcmOld = Wcm;

	Wcm = col(0);

	for (size_t i = 0; i < afl.np + 1; ++i)
		afl.v[i] = afl.r[i].kcross() * Wcm;
}