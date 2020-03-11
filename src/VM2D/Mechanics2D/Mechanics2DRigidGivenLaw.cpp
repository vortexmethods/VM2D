/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.8    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/03/09     |
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
\brief Файл кода с описанием класса MechanicsRigidGivenLaw
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.8   
\date 09 марта 2020 г.
*/

#include <algorithm>

#include "Mechanics2DRigidGivenLaw.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"


#include "MeasureVP2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
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

	Point2D hDFGam = { 0.0, 0.0 };	//гидродинамические силы, обусловленные Gamma_k
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные delta_k
	double hDMdelta = 0.0;

	Point2D deltaVstep = VeloOfAirfoilRcm(W.getCurrentStep() * dt) - VeloOfAirfoilRcm((std::max(W.getCurrentStep(), (size_t)1) - 1) * dt);
		
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


const double accelStop = 1000.0;


// Вычисление скорости центра масс
Point2D MechanicsRigidGivenLaw::VeloOfAirfoilRcm(double currTime)
{
//	if(Comega * currTime < accelStop * 2. * PI)
//	{
//	    return{ CA * Comega * cos(Comega * currTime), 0.0 };
//	}
//	else return{ CA * Comega * cos(accelStop * 2. * PI), 0.0 };
	
	return{ CA*Comega*cos(Comega*currTime), 0.0 };
	//return { 1.0, 0.0 };
	
	//Тестовый вариант 
	//if (currTime < 10.0)
	//	return { -currTime / 10.0, 0.0 };
	//else return { -1.0, 0.0};


}//VeloOfAirfoilRcm(...)

// Вычисление положения центра масс
Point2D MechanicsRigidGivenLaw::PositionOfAirfoilRcm(double currTime)
{

//	if(Comega * currTime < accelStop * 2. * PI)
//	{
//	    return{ CA * sin(Comega * currTime), 0.0 };
//	}
//	else return{ CA * sin(accelStop * 2. * PI) + CA * Comega * cos(accelStop * 2. * PI) * (currTime - accelStop * 2. * PI / Comega), 0.0 };
	return{ CA*sin(Comega*currTime), 0.0};
	
	//return { 1.0*currTime, 0.0 };

	//Тестовый вариант
	//if (currTime < 10.0)
	//	return { -currTime / 10.0 * currTime * 0.5, 0.0 };
	//else return { -5.0 - (currTime - 10.0), 0.0 };

}//PositionOfAirfoilRcm(...)

// Вычисление скоростей начал панелей
void MechanicsRigidGivenLaw::VeloOfAirfoilPanels(double currTime)
{
	Point2D veloRcm = VeloOfAirfoilRcm(currTime);
	afl.setV(veloRcm);	
}//VeloOfAirfoilPanels(...)

void MechanicsRigidGivenLaw::Move()
{
	//Point2D airfoilVelo = VeloOfAirfoilRcm(W.getPassport().physicalProperties.getCurrTime());
	Point2D aflRcmOld = afl.rcm;
	afl.Move(PositionOfAirfoilRcm(W.getPassport().physicalProperties.getCurrTime() + W.getPassport().timeDiscretizationProperties.dt) - aflRcmOld);
}//Move()


void MechanicsRigidGivenLaw::ReadSpecificParametersFromDictionary()
{
    mechParamsParser->get("Comega", Comega);
    
    W.getInfo('i') << "frequency " << "Comega = " << Comega << std::endl;
	
    mechParamsParser->get("CA", CA);
		
    W.getInfo('i') << "Amplitude CA = " << CA << std::endl;
}//ReadSpecificParametersFromDictionary()