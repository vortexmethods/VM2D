/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DRigidRotatePart.cpp                                   |
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
\brief Файл кода с описанием класса Mechanics2DRigidRotatePart
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#include "Mechanics2DRigidRotatePart.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;


MechanicsRigidRotatePart::MechanicsRigidRotatePart(const World2D& W_, size_t numberInPassport_)
	: Mechanics(W_, numberInPassport_, true, false)
	//w0(0.0), 
	//phi0(W_.getAirfoil(numberInPassport_).phiAfl)
{
	Vcm0 = { 0.0, 0.0 };
	Rcm0 = { W_.getAirfoil(numberInPassport_).rcm[0], W_.getAirfoil(numberInPassport_).rcm[1] };
	Vcm = Vcm0;
	Rcm = Rcm0;
	VcmOld = Vcm0;
	RcmOld = Rcm;
		
	ReadSpecificParametersFromDictionary();
	Initialize({ 0.0, 0.0 }, W_.getAirfoil(numberInPassport_).rcm, 0.0, W_.getAirfoil(numberInPassport_).phiAfl);
};

//Вычисление гидродинамической силы, действующей на профиль
void MechanicsRigidRotatePart::GetHydroDynamForce()
{
	W.getTimestat().timeGetHydroDynamForce.first += omp_get_wtime();

	const double& dt = W.getPassport().timeDiscretizationProperties.dt;

	hydroDynamForce = { 0.0, 0.0 };
	hydroDynamMoment = 0.0;

	viscousForce = { 0.0, 0.0 };
	viscousMoment = 0.0;

	Point2D hDFGam = { 0.0, 0.0 };	    //гидродинамические силы, обусловленные присоед.завихренностью
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные приростом завихренности
	Point2D hDFQ = { 0.0, 0.0 };		//гидродинамические силы, обусловленные присоед.источниками

	double hDMGam = 0.0;				//гидродинамический момент, обусловленный присоед.завихренностью
	double hDMdelta = 0.0;				//гидродинамический момент, обусловленный приростом завихренности
	double hDMQ = 0.0;					//гидродинамический момент, обусловленный присоед.источниками


	

	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		Point2D rK = 0.5 * (afl.getR(i + 1) + afl.getR(i)) - afl.rcm;	

		Point2D velK = { -Wcm * rK[1], Wcm * rK[0] };

		double gAtt = Wcm * (rK ^ afl.tau[i]);
		double gAttOld = WcmOld * (rK ^ afl.tau[i]);
		double deltaGAtt = gAtt - gAttOld;

		double qAtt = Wcm * (rK ^ afl.nrm[i]);

		/// \todo Учитываем только нулевой момент решения. Надо ли учитывать остальные?
		double deltaK = boundary.sheets.freeVortexSheet(i, 0) * afl.len[i] - afl.gammaThrough[i] + deltaGAtt * afl.len[i];
		
		/*1*/
		hDFdelta += deltaK * Point2D({ -rK[1], rK[0] });
		hDMdelta += 0.5 * deltaK * rK.length2();

		/*2*/
		hDFGam += 0.25 * (afl.getV(i) + afl.getV(i + 1)).kcross() * gAtt * afl.len[i];
		hDMGam += 0.25 * rK ^ (afl.getV(i) + afl.getV(i + 1)).kcross() * gAtt * afl.len[i];

		/*3*/
		hDFQ -= 0.25 * (afl.getV(i) + afl.getV(i + 1)) * qAtt * afl.len[i];
		hDMQ -= 0.25 * rK ^ (afl.getV(i) + afl.getV(i + 1)) * qAtt * afl.len[i];
	}

	const double rho = W.getPassport().physicalProperties.rho;

	hydroDynamForce = rho * (hDFGam + hDFdelta * (1.0 / dt) + hDFQ);
	hydroDynamMoment = rho * (hDMGam + hDMdelta / dt + hDMQ);

	if ((W.getPassport().physicalProperties.nu > 0.0)/* && (W.currentStep > 0)*/)
		for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
		{
			Point2D rK = 0.5 * (afl.getR(i + 1) + afl.getR(i)) - afl.rcm;
			viscousForce += rho * afl.viscousStress[i] * afl.tau[i];
			viscousMoment += rho * (afl.viscousStress[i] * afl.tau[i]) & rK;
		}

	W.getTimestat().timeGetHydroDynamForce.second += omp_get_wtime();
}// GetHydroDynamForce()

// Вычисление скорости центра масс
Point2D MechanicsRigidRotatePart::VeloOfAirfoilRcm(double currTime)
{
	return Point2D{0.0, 0.0};
}//VeloOfAirfoilRcm(...)

// Вычисление положения центра масс
Point2D MechanicsRigidRotatePart::PositionOfAirfoilRcm(double currTime)
{
	return Rcm;
}//PositionOfAirfoilRcm(...)

double MechanicsRigidRotatePart::AngularVelocityOfAirfoil(double currTime)
{
	return Wcm;
}//AngularVelocityOfAirfoil(...)

double MechanicsRigidRotatePart::AngleOfAirfoil(double currTime)
{
	return afl.phiAfl;
}//AngleOfAirfoil(...)

// Вычисление скоростей начал панелей
void MechanicsRigidRotatePart::VeloOfAirfoilPanels(double currTime)
{
	Point2D veloRcm = VeloOfAirfoilRcm(currTime);

	std::vector<Point2D> veloW(afl.getNumberOfPanels());
	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
		veloW[i] = veloRcm + Wcm * (afl.getR(i) - Rcm).kcross();

	afl.setV(veloW);

	circulationOld = circulation;
	circulation = 2.0 * afl.area * Wcm;
}//VeloOfAirfoilPanels(...)


void MechanicsRigidRotatePart::Move()
{
	double Jeff = J;

	PhiOld = Phi;
	WcmOld = Wcm;
	

	Point2D dr, dV;
	double dphi, dw;

	//W.getInfo('t') << "k = " << k << std::endl;


	dr[1] = 0.0;
	dV[1] = 0.0;

	dr[0] = 0.0;
	dV[0] = 0.0;
	
	double bw = 0.0;
	double kw = 0.0;
		

	double dt = W.getPassport().timeDiscretizationProperties.dt;
	double t = W.getPassport().physicalProperties.getCurrTime();
	
	
	double addMom = (t > 3.0 * tAccel) ? externalTorque : 0.0;

	if (t > 1.0 * tAccel)
	{
		Point2D kk[4];
		kk[0] = { Wcm, (hydroDynamMoment - 2.0 * bw * Wcm - kw * Phi - addMom) / Jeff };
		kk[1] = { Wcm + 0.5 * dt * kk[0][1], (hydroDynamMoment - 2.0 * bw * (Wcm + 0.5 * dt * kk[0][1]) - kw * (Phi + 0.5 * dt * kk[0][0]) - addMom) / Jeff };
		kk[2] = { Wcm + 0.5 * dt * kk[1][1], (hydroDynamMoment - 2.0 * bw * (Wcm + 0.5 * dt * kk[1][1]) - kw * (Phi + 0.5 * dt * kk[1][0]) - addMom) / Jeff };
		kk[3] = { Wcm + dt * kk[2][1], (hydroDynamMoment - 2.0 * bw * (Wcm + dt * kk[2][1]) - kw * (Phi + dt * kk[2][0]) - addMom) / Jeff };

		dphi = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		dw = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}
	else 
	{
		double a = wAccel / tAccel;
		dw = dt * a;
		dphi = 0.5 * a * sqr(t) - 0.5 * a * sqr(t - dt);
	}

	
	afl.Move(dr);
	afl.Rotate(dphi);

	Rcm += dr;
	Vcm += dV;

	Phi += dphi;

	//std::cout << "Phi = " << Phi << ", afl.Phi = " << afl.phiAfl << std::endl;

	Wcm += dw;

}//Move()



void MechanicsRigidRotatePart::ReadSpecificParametersFromDictionary()
{
	mechParamsParser->get("J", J);

	W.getInfo('i') << "moment of inertia " << "J = " << J << std::endl;

	mechParamsParser->get("wAccel", wAccel);

	W.getInfo('i') << "wAccel = " << wAccel << std::endl;

	mechParamsParser->get("tAccel", tAccel);

	W.getInfo('i') << "tAccel = " << tAccel << std::endl;

	mechParamsParser->get("externalTorque", externalTorque);

	W.getInfo('i') << "externalTorque = " << externalTorque << std::endl;

}//ReadSpecificParametersFromDictionary()
