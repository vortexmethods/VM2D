/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.11   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2022/08/07     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
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
\version 1.11
\date 07 августа 2022 г.
*/

#include "mpi.h"

#include "Mechanics2DRigidRotatePart.h"

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


MechanicsRigidRotatePart::MechanicsRigidRotatePart(const World2D& W_, size_t numberInPassport_)
	: Mechanics(W_, numberInPassport_, 1, true, false, true), 
	w0(0.0), 
	phi0(W_.getAirfoil(numberInPassport_).phiAfl)
{
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
//		Point2D Vcm = VeloOfAirfoilRcm(W.getCurrentStep() * dt);
//		Point2D VcmOld = VeloOfAirfoilRcm((std::max(W.getCurrentStep(), (size_t)1) - 1) * dt);

//		double Wcm = AngularVelocityOfAirfoil(W.getCurrentStep() * dt);
//		double WcmOld = AngularVelocityOfAirfoil((std::max(W.getCurrentStep(), (size_t)1) - 1) * dt);

		Point2D rK = 0.5 * (afl.getR(i + 1) + afl.getR(i));

		Point2D dr = rK - afl.rcm;

		Point2D velK = { -Wcm * dr[1], Wcm * dr[0] };

		double gAtt = Wcm * (dr ^ afl.tau[i]);
		double gAttOld = WcmOld * (dr ^ afl.tau[i]);
		double deltaGAtt = gAtt - gAttOld;

		double qAtt = Wcm * (dr ^ afl.nrm[i]);

		/// \todo Учитываем только нулевой момент решения. Надо ли учитывать остальные?
		double deltaK = boundary.sheets.freeVortexSheet(i, 0) * afl.len[i] - afl.gammaThrough[i] + deltaGAtt * afl.len[i];

		/*1*/
		hDFdelta += deltaK * Point2D({ -dr[1], dr[0] });
		hDMdelta += 0.5 * deltaK * dr.length2();

		/*2*/
		hDFGam += 0.25 * (afl.getV(i) + afl.getV(i + 1)).kcross() * gAtt * afl.len[i];
		hDMGam += 0.25 * dr ^ (afl.getV(i) + afl.getV(i + 1)).kcross() * gAtt * afl.len[i];

		/*3*/
		hDFQ -= 0.25 * (afl.getV(i) + afl.getV(i + 1)) * qAtt * afl.len[i];
		hDMQ -= 0.25 * dr ^ (afl.getV(i) + afl.getV(i + 1)) * qAtt * afl.len[i];
	}

	const double rho = W.getPassport().physicalProperties.rho;

	hydroDynamForce = rho * (hDFGam + hDFdelta * (1.0 / dt) + hDFQ);
	hydroDynamMoment = rho * (hDMGam + hDMdelta / dt + hDMQ);

	if (W.getPassport().physicalProperties.nu > 0.0)
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
	return { 0.0, 0.0 };
}//VeloOfAirfoilRcm(...)

// Вычисление положения центра масс
Point2D MechanicsRigidRotatePart::PositionOfAirfoilRcm(double currTime)
{
	return afl.rcm;
}//PositionOfAirfoilRcm(...)

// Вычисление угловой скорости профиля вокруг центра масс
double MechanicsRigidRotatePart::AngularVelocityOfAirfoil(double currTime)
{
	return w0;
}//AngleOfAirfoil(...)

// Вычисление угла поворота профиля вокруг центра масс
double MechanicsRigidRotatePart::AngleOfAirfoil(double currTime)
{
	return afl.phiAfl;
}//AngleOfAirfoil(...)

// Вычисление скоростей начал панелей
void MechanicsRigidRotatePart::VeloOfAirfoilPanels(double currTime)
{
	Point2D veloRcm = VeloOfAirfoilRcm(currTime);

	std::vector<Point2D> vel(afl.getNumberOfPanels(), { 0.0, 0.0 });
	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		vel[i][0] = (afl.getR(i) - afl.rcm).kcross()[0] * Wcm;
		vel[i][1] = (afl.getR(i) - afl.rcm).kcross()[1] * Wcm;
	}

	afl.setV(vel);
}//VeloOfAirfoilPanels(...)


void MechanicsRigidRotatePart::Move()
{
	WcmOld = Wcm;
	PhiOld = Phi;

	//До момента времени tAccel ротор равномерно ускоряется до скорости wTilde
	
	double t = W.getPassport().physicalProperties.getCurrTime();
	double dt = W.getPassport().timeDiscretizationProperties.dt;

	if (t < tAccel)
	{
		Wcm = t * wAccel / tAccel;
		Phi = phi0 + Wcm * t / 2.;

		W.getInfo('i') << "Wcm * t / 2. " << Wcm * t / 2. << std::endl;
	}
	else if (t < 3.0 * tAccel)
	{
		Wcm = WcmOld + (hydroDynamMoment + 0.0 * viscousMoment) * dt / J;
		Phi = PhiOld + Wcm * dt;
		
	}
	else
	{
		Wcm = WcmOld + (hydroDynamMoment + 0.0 * viscousMoment - externalTorque) * dt / J;
		Phi = PhiOld + Wcm * dt;

//		W.getInfo('i') << "Wcm * dt " << Wcm * dt << std::endl;
	}

	afl.Rotate(Phi - PhiOld);

	

//	W.getInfo('i') << "hydroDynamMoment " << hydroDynamMoment << std::endl;

//	W.getInfo('i') << "externalTorque " << externalTorque << std::endl;

//	W.getInfo('i') << "J " << J << std::endl;

//	W.getInfo('i') << "Wcm "  << Wcm << std::endl;

//	W.getInfo('i') << "Phi " << Phi << std::endl;

//	W.getInfo('i') << "PhiOld " << PhiOld << std::endl;

//	W.getInfo('i') << "Phi - PhiOld "  << Phi - PhiOld << std::endl;

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
