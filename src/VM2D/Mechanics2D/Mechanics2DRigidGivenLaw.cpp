/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.9    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/07/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DRigidGivenLaw.cpp                                     |
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
\brief Файл кода с описанием класса MechanicsRigidGivenLaw
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.9   
\date 22 июля 2020 г.
*/

#include <algorithm>

#include "Mechanics2DRigidGivenLaw.h"

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
void MechanicsRigidGivenLaw::GetHydroDynamForce()
{
	W.getTimestat().timeGetHydroDynamForce.first += omp_get_wtime();

	const double& dt = W.getPassport().timeDiscretizationProperties.dt;

	hydroDynamForce = { 0.0, 0.0 };
	hydroDynamMoment = 0.0;

	viscousForce = { 0.0, 0.0 };
	viscousMoment = 0.0;

	Point2D hDFGam = { 0.0, 0.0 };	    //гидродинамические силы, обусловленные присоед.завихренностью
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные проростом завихренности
	Point2D hDFQ = { 0.0, 0.0 };		//гидродинамические силы, обусловленные присоед.источниками
	
	double hDMGam = 0.0;				//гидродинамический момент, обусловленный присоед.завихренностью
	double hDMdelta = 0.0;				//гидродинамический момент, обусловленный проростом завихренности
	double hDMQ = 0.0;					//гидродинамический момент, обусловленный присоед.источниками

		
	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		Point2D Vcm = VeloOfAirfoilRcm(W.getCurrentStep() * dt);
		Point2D VcmOld = VeloOfAirfoilRcm((std::max(W.getCurrentStep(), (size_t)1) - 1) * dt);
		
		double Wcm = AngularVelocityOfAirfoil(W.getCurrentStep() * dt);
		double WcmOld = AngularVelocityOfAirfoil((std::max(W.getCurrentStep(), (size_t)1) - 1) * dt);

		Point2D rK = 0.5 * (afl.getR(i + 1) + afl.getR(i));
		
		Point2D dr = rK - afl.rcm;

		Point2D velK = { -Wcm * dr[1], Wcm * dr[0] };
		velK += Vcm;

		double gAtt = Wcm * (dr ^ afl.tau[i]) + (Vcm & afl.tau[i]);
		double gAttOld = WcmOld * (dr ^ afl.tau[i]) + (VcmOld & afl.tau[i]);
		double deltaGAtt = gAtt - gAttOld;

		double qAtt = Wcm * (dr ^ afl.nrm[i]) + (Vcm & afl.nrm[i]);
	
		/// \todo Учитываем только нулевой момент решения. Надо ли учитывать остальные?
		double deltaK = boundary.sheets.freeVortexSheet(i, 0) * afl.len[i] - afl.gammaThrough[i] + deltaGAtt * afl.len[i];
		
		/*1*/
		hDFdelta += deltaK * Point2D({ -dr[1], dr[0] });
		hDMdelta += 0.5 * deltaK * dr.length2();

		/*2*/
		hDFGam += 0.25 * (afl.getV(i) + afl.getV(i+1)).kcross() * gAtt * afl.len[i];
		hDMGam += 0.25 * dr ^ (afl.getV(i) + afl.getV(i + 1)).kcross() * gAtt * afl.len[i];

		/*3*/
		hDFQ -= 0.25 * (afl.getV(i) + afl.getV(i + 1)) * qAtt * afl.len[i];
		hDMQ -= 0.25 * dr ^ (afl.getV(i) + afl.getV(i + 1)) * qAtt * afl.len[i];		
	}

	hydroDynamForce = hDFGam + hDFdelta * (1.0 / dt) + hDFQ;
	hydroDynamMoment = hDMGam + hDMdelta / dt + hDMQ;

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
Point2D MechanicsRigidGivenLaw::VeloOfAirfoilRcm(double currTime)
{
	return VelocityOfCenterOfMass(currTime);
}//VeloOfAirfoilRcm(...)


// Вычисление положения центра масс
Point2D MechanicsRigidGivenLaw::PositionOfAirfoilRcm(double currTime)
{	
	return PositionOfCenterOfMass(currTime);
}//PositionOfAirfoilRcm(...)


// Вычисление угла поворота профиля вокруг центра масс
double MechanicsRigidGivenLaw::AngleOfAirfoil(double currTime)
{
	return RotationAngle(currTime);	
}//AngleOfAirfoil(...)

// Вычисление угловой скорости профиля вокруг центра масс
double MechanicsRigidGivenLaw::AngularVelocityOfAirfoil(double currTime)
{
	return AngularVelocity(currTime);
}//AngularVelocityOfAirfoil(...)





// Вычисление скоростей начал панелей
void MechanicsRigidGivenLaw::VeloOfAirfoilPanels(double currTime)
{
	Point2D veloRcm = VeloOfAirfoilRcm(currTime);

	std::vector<Point2D> vel(afl.getNumberOfPanels(), { 0.0, 0.0 });
	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		vel[i][0] = (afl.getR(i) - afl.rcm).kcross()[0] * Wcm;
		vel[i][1] = (afl.getR(i) - afl.rcm).kcross()[1] * Wcm;
		vel[i] += veloRcm;
	}
	
	afl.setV(vel);
}//VeloOfAirfoilPanels(...)

void MechanicsRigidGivenLaw::Move()
{
	double t = W.getPassport().physicalProperties.getCurrTime();
	double dt = W.getPassport().timeDiscretizationProperties.dt;

	//Point2D airfoilVelo = VeloOfAirfoilRcm(t);
	Point2D aflRcmOld = afl.rcm;
	double aflPhiOld = afl.phiAfl;
	afl.Move(PositionOfAirfoilRcm(t + dt) - aflRcmOld);
	afl.Rotate(AngleOfAirfoil(t + dt) - aflPhiOld);

	Vcm = VeloOfAirfoilRcm(t + dt);
	Phi = AngleOfAirfoil(t + dt);
	Wcm = AngularVelocityOfAirfoil(t + dt);
}//Move()


void MechanicsRigidGivenLaw::ReadSpecificParametersFromDictionary()
{
   /*
	mechParamsParser->get("Comega", Comega);
    
    W.getInfo('i') << "frequency " << "Comega = " << Comega << std::endl;
	
    mechParamsParser->get("CA", CA);
		
    W.getInfo('i') << "Amplitude CA = " << CA << std::endl;
	*/
}//ReadSpecificParametersFromDictionary()