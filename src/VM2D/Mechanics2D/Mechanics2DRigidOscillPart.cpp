/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
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
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#include "Mechanics2DRigidOscillPart.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

MechanicsRigidOscillPart::MechanicsRigidOscillPart(const World2D& W_, size_t numberInPassport_)
	:
	Mechanics(W_, numberInPassport_, true, false)
	//, V0({ 0.0, 0.0 })
	//, r0({ W_.getAirfoil(numberInPassport_).rcm[0], W_.getAirfoil(numberInPassport_).rcm[1] })
{
	Vcm0 = { 0.0, 0.0 };
	Rcm0 = { W_.getAirfoil(numberInPassport_).rcm[0], W_.getAirfoil(numberInPassport_).rcm[1] };
	Vcm = Vcm0;
	Rcm = Rcm0;
	VcmOld = Vcm0;
	RcmOld = Rcm;

	strongCoupling = false;
		
	ReadSpecificParametersFromDictionary();
	Initialize(initVelocity, W_.getAirfoil(numberInPassport_).rcm + initDisplacement, initAngularVelocity, W_.getAirfoil(numberInPassport_).phiAfl + initAngularDisplacement);
};

//Вычисление гидродинамической силы, действующей на профиль
void MechanicsRigidOscillPart::GetHydroDynamForce()
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

		Point2D velK = 0.5 * (afl.getV(i) + afl.getV(i + 1));
		double gAtt = (velK & afl.tau[i]);

		double gAttOld = 0.0;
		if (W.currentStep > 0)
		{
			auto oldAfl = W.getOldAirfoil(numberInPassport);
			gAttOld = ((0.5 * (oldAfl.getV(i) + oldAfl.getV(i + 1))) & oldAfl.tau[i]);
		}

		double deltaGAtt = gAtt - gAttOld;

		double qAtt = (velK & afl.nrm[i]);

		/// \todo Учитываем только нулевой момент решения. Надо ли учитывать остальные?
		double deltaK = boundary.sheets.freeVortexSheet(i, 0) * afl.len[i] - afl.gammaThrough[i] + deltaGAtt * afl.len[i];

		/*1*/
		hDFdelta += deltaK * Point2D({ -rK[1], rK[0] });
		hDMdelta += 0.5 * deltaK * rK.length2();

		/*2*/
		hDFGam += 0.5 * velK.kcross() * gAtt * afl.len[i];
		hDMGam += 0.5 * (rK ^ velK.kcross()) * gAtt * afl.len[i];

		/*3*/
		hDFQ -= 0.5 * velK * qAtt * afl.len[i];
		hDMQ -= 0.5 * (rK ^ velK) * qAtt * afl.len[i];
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
Point2D MechanicsRigidOscillPart::VeloOfAirfoilRcm(double currTime)
{
	return Vcm;
}//VeloOfAirfoilRcm(...)

// Вычисление положения центра масс
Point2D MechanicsRigidOscillPart::PositionOfAirfoilRcm(double currTime)
{
	return Rcm;
}//PositionOfAirfoilRcm(...)

double MechanicsRigidOscillPart::AngularVelocityOfAirfoil(double currTime)
{
	return Wcm;
}//AngularVelocityOfAirfoil(...)

double MechanicsRigidOscillPart::AngleOfAirfoil(double currTime)
{
	if (afl.phiAfl != Phi)
	{
		std::cout << "afl.phiAfl != Phi" << std::endl;
		exit(100600);
	}

	return afl.phiAfl;
}//AngleOfAirfoil(...)

// Вычисление скоростей начал панелей
void MechanicsRigidOscillPart::VeloOfAirfoilPanels(double currTime)
{
	Point2D veloRcm = VeloOfAirfoilRcm(currTime);

	std::vector<Point2D> veloW(afl.getNumberOfPanels());
	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
		veloW[i] = veloRcm + Wcm * (afl.getR(i) - Rcm).kcross();

	afl.setV(veloW);

	circulationOld = circulation;
	circulation = 2.0 * afl.area * Wcm;

}//VeloOfAirfoilPanels(...)


void MechanicsRigidOscillPart::Move()
{
	Point2D meff;
	//if (W.getPassport().airfoilParams[numberInPassport].addedMass.length2() > 0)
	//	meff = Point2D{ m + W.getPassport().airfoilParams[numberInPassport].addedMass[0], m + W.getPassport().airfoilParams[numberInPassport].addedMass[1] };
	//else
		meff = Point2D{ m, m };

		double Jeff = J;

	VcmOld = Vcm;
	RcmOld = Rcm;
	PhiOld = Phi;
	WcmOld = Wcm;	

	Point2D dr, dV;
	double dphi, dw;

	//W.getInfo('t') << "k = " << k << std::endl;

	if (k[0] > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		Point2D kk[4];

		kk[0] = { Vcm[0], (hydroDynamForce[0] - 2.0 * b[0] * Vcm[0] - k[0] * Rcm[0]) / meff[0] };
		kk[1] = { Vcm[0] + 0.5 * dt * kk[0][1], (hydroDynamForce[0] - 2.0 * b[0] * (Vcm[0] + 0.5 * dt * kk[0][1]) - k[0] * (Rcm[0] + 0.5 * dt * kk[0][0])) / meff[0] };
		kk[2] = { Vcm[0] + 0.5 * dt * kk[1][1], (hydroDynamForce[0] - 2.0 * b[0] * (Vcm[0] + 0.5 * dt * kk[1][1]) - k[0] * (Rcm[0] + 0.5 * dt * kk[1][0])) / meff[0] };
		kk[3] = { Vcm[0] + dt * kk[2][1], (hydroDynamForce[0] - 2.0 * b[0] * (Vcm[0] + dt * kk[2][1]) - k[0] * (Rcm[0] + dt * kk[2][0])) / meff[0] };

		dr[0] = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		dV[0] = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dr[0] = 0.0;
		dV[0] = 0.0;
	}



	if (k[1] > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		Point2D kk[4];
		kk[0] = { Vcm[1], (hydroDynamForce[1] - 2.0 * b[1] * Vcm[1] - k[1] * Rcm[1]) / meff[1]};
		kk[1] = { Vcm[1] + 0.5 * dt * kk[0][1], (hydroDynamForce[1] - 2.0 * b[1] * (Vcm[1] + 0.5 * dt * kk[0][1]) - k[1] * (Rcm[1] + 0.5 * dt * kk[0][0])) / meff[1]};
		kk[2] = { Vcm[1] + 0.5 * dt * kk[1][1], (hydroDynamForce[1] - 2.0 * b[1] * (Vcm[1] + 0.5 * dt * kk[1][1]) - k[1] * (Rcm[1] + 0.5 * dt * kk[1][0])) / meff[1]};
		kk[3] = { Vcm[1] + dt * kk[2][1], (hydroDynamForce[1] - 2.0 * b[1] * (Vcm[1] + dt * kk[2][1]) - k[1] * (Rcm[1] + dt * kk[2][0])) / meff[1]};
				
		dr[1] = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		dV[1] = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dr[1] = 0.0;
		dV[1] = 0.0;
	}



	if (kw > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		Point2D kk[4];

		kk[0] = { Wcm, (hydroDynamMoment - 2.0 * bw * Wcm - kw * Phi) / Jeff };
		kk[1] = { Wcm + 0.5 * dt * kk[0][1], (hydroDynamMoment - 2.0 * bw * (Wcm + 0.5 * dt * kk[0][1]) - kw * (Phi + 0.5 * dt * kk[0][0])) / Jeff };
		kk[2] = { Wcm + 0.5 * dt * kk[1][1], (hydroDynamMoment - 2.0 * bw * (Wcm + 0.5 * dt * kk[1][1]) - kw * (Phi + 0.5 * dt * kk[1][0])) / Jeff };
		kk[3] = { Wcm + dt * kk[2][1], (hydroDynamMoment - 2.0 * bw * (Wcm + dt * kk[2][1]) - kw * (Phi + dt * kk[2][0])) / Jeff };

		dphi = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		dw = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dphi = 0.0;
		dw = 0.0;
	}

	afl.Move(dr);
	afl.Rotate(dphi);

	Rcm += dr;
	Vcm += dV;

	Phi += dphi;
	Wcm += dw;
}//Move()



void MechanicsRigidOscillPart::MoveKinematic()
{
	VcmOld = Vcm;
	RcmOld = Rcm;

	Point2D dr, dV;
	double dt = W.getPassport().timeDiscretizationProperties.dt;

	if (k[1] > 0)
	{			
		dr[1] = Vcm[1] * dt;
		dV[1] = 0.0;
	}
	else
	{
		dr[1] = 0.0;
		dV[1] = 0.0;
	}


	if (k[0] > 0)
	{		
		dr[0] = Vcm[0] * dt;
		dV[0] = 0.0;
	}
	else
	{
		dr[0] = 0.0;
		dV[0] = 0.0;
	}

	afl.Move(dr);	
			
	Rcm += dr;
	Vcm += dV;
}//MoveKinematic()


void MechanicsRigidOscillPart::MoveOnlyVelo()
{
	Point2D meff;
	if (W.getPassport().airfoilParams[numberInPassport].addedMass.length2() > 0)
		meff = Point2D{ m + W.getPassport().airfoilParams[numberInPassport].addedMass[0], m + W.getPassport().airfoilParams[numberInPassport].addedMass[1] };
	else
		meff = Point2D{ m, m };

	Point2D dr, dV;	

	if (k[1] > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		Point2D kk[4];
	
		kk[0] = { Vcm[1], (hydroDynamForce[1] - 2.0 * b[1] * Vcm[1] - k[1] * Rcm[1]) / meff[1]};
		kk[1] = { Vcm[1] + 0.5 * dt * kk[0][1], (hydroDynamForce[1] - 2.0 * b[1] * (Vcm[1] + 0.5 * dt * kk[0][1]) - k[1] * (Rcm[1] + 0.5 * dt * kk[0][0])) / meff[1]};
		kk[2] = { Vcm[1] + 0.5 * dt * kk[1][1], (hydroDynamForce[1] - 2.0 * b[1] * (Vcm[1] + 0.5 * dt * kk[1][1]) - k[1] * (Rcm[1] + 0.5 * dt * kk[1][0])) / meff[1]};
		kk[3] = { Vcm[1] + dt * kk[2][1], (hydroDynamForce[1] - 2.0 * b[1] * (Vcm[1] + dt * kk[2][1]) - k[1] * (Rcm[1] + dt * kk[2][0])) / meff[1]};

		dr[1] = 0.0;// dt* (kk[0][0] + 2.0 * kk[1][0] + 2.0 * kk[2][0] + kk[3][0]) / 6.0;
		dV[1] = dt * (kk[0][1] + 2.0 * kk[1][1] + 2.0 * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dr[1] = 0.0;
		dV[1] = 0.0;
	}


	if (k[0] > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		Point2D kk[4];
		
		kk[0] = { Vcm[0], (hydroDynamForce[0] - 2.0 * b[0] * Vcm[0] - k[0] * Rcm[0]) / meff[0]};
		kk[1] = { Vcm[0] + 0.5 * dt * kk[0][1], (hydroDynamForce[0] - 2.0 * b[0] * (Vcm[0] + 0.5 * dt * kk[0][1]) - k[0] * (Rcm[0] + 0.5 * dt * kk[0][0])) / meff[0]};
		kk[2] = { Vcm[0] + 0.5 * dt * kk[1][1], (hydroDynamForce[0] - 2.0 * b[0] * (Vcm[0] + 0.5 * dt * kk[1][1]) - k[0] * (Rcm[0] + 0.5 * dt * kk[1][0])) / meff[0]};
		kk[3] = { Vcm[0] + dt * kk[2][1], (hydroDynamForce[0] - 2.0 * b[0] * (Vcm[0] + dt * kk[2][1]) - k[0] * (Rcm[0] + dt * kk[2][0])) / meff[0]};

		dr[0] = 0.0; // dt* (kk[0][0] + 2.0 * kk[1][0] + 2.0 * kk[2][0] + kk[3][0]) / 6.0;
		dV[0] = dt * (kk[0][1] + 2.0 * kk[1][1] + 2.0 * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dr[0] = 0.0;
		dV[0] = 0.0;
	}

	//afl.Move({ dx, dy });
	Vcm += dV;	
	
}//MoveOnlyVelo()



#if defined(INITIAL) || defined(BRIDGE) 
void MechanicsRigidOscillPart::ReadSpecificParametersFromDictionary()
{
	/*
	mechParamsParser->get("m", m);	
	W.getInfo('i') << "mass " << "m = " << m << std::endl;

	mechParamsParser->get("J", J);
	W.getInfo('i') << "moment of inertia " << "J = " << J << std::endl;

	
	mechParamsParser->get("k", k);
	mechParamsParser->get("kw", kw);

	W.getInfo('i') << "linear rigidity (kx, ky) = " << k << std::endl;	
	W.getInfo('i') << "rotational rigidity kw = " << kw << std::endl;	

	Point2D c;
	mechParamsParser->get("c", c);//Логарифмический декремент
	double cw;
	mechParamsParser->get("cw", cw);//Логарифмический декремент

	b[0] = c[0] / 2.0;
	b[1] = c[1] / 2.0;
	bw = cw / 2.0;

	W.getInfo('i') << "linear damping (bx, by) = " << b << std::endl;
	W.getInfo('i') << "rotational damping bw = " << bw << std::endl;

	mechParamsParser->get("initDisplacement", initDisplacement, &defaults::defaultInitDisplacement);
	mechParamsParser->get("initAngularDisplacement", initAngularDisplacement, &defaults::defaultInitAngularDisplacement);
	W.getInfo('i') << "initial displacement: " << "translational = " << initDisplacement << ", rotational = " << initAngularDisplacement << std::endl;

	mechParamsParser->get("initVelocity", initVelocity, &defaults::defaultInitVelocity);
	mechParamsParser->get("initAngularVelocity", initAngularVelocity, &defaults::defaultInitAngularVelocity);
	W.getInfo('i') << "initial velocity: " << "translational = " << initVelocity << ", rotational = " << initAngularVelocity << std::endl;

	//*/

	//*
	mechParamsParser->get("m", m);	
	W.getInfo('i') << "mass " << "m = " << m << std::endl;

	mechParamsParser->get("J", J);
	W.getInfo('i') << "moment of inertia " << "J = " << J << std::endl;
		
	Point2D sh;
	mechParamsParser->get("sh", sh);	
	k[0] = m * sqr(2.0 * PI * sh[0] / W.getPassport().airfoilParams[numberInPassport].chord) * W.getPassport().physicalProperties.vInf.length2();
	k[1] = m * sqr(2.0 * PI * sh[1] / W.getPassport().airfoilParams[numberInPassport].chord) * W.getPassport().physicalProperties.vInf.length2();
	
	double shw;
	mechParamsParser->get("shw", shw);
	kw = J * sqr(2.0 * PI * shw / W.getPassport().airfoilParams[numberInPassport].chord) * W.getPassport().physicalProperties.vInf.length2();
	
	W.getInfo('i') << "linear rigidity (kx, ky) = " << k << std::endl;	
	W.getInfo('i') << "rotational rigidity kw = " << kw << std::endl;	

	Point2D zeta;
	mechParamsParser->get("zeta", zeta);//Логарифмический декремент
	double zetaw;
	mechParamsParser->get("zetaw", zetaw);//Логарифмический декремент

	b[0] = zeta[0] / (2.0 * PI) * sqrt(k[0] * m);
	b[1] = zeta[1] / (2.0 * PI) * sqrt(k[1] * m);
	bw = zetaw / (2.0 * PI) * sqrt(kw * J);

	W.getInfo('i') << "linear damping (bx, by) = " << b << std::endl;
	W.getInfo('i') << "rotational damping bw = " << bw << std::endl;


	mechParamsParser->get("initDisplacement", initDisplacement, &defaults::defaultInitDisplacement);
	mechParamsParser->get("initAngularDisplacement", initAngularDisplacement, &defaults::defaultInitAngularDisplacement);
	W.getInfo('i') << "initial displacement: " << "translational = " << initDisplacement << ", rotational = " << initAngularDisplacement << std::endl;

	mechParamsParser->get("initVelocity", initVelocity, &defaults::defaultInitVelocity);
	mechParamsParser->get("initAngularVelocity", initAngularVelocity, &defaults::defaultInitAngularVelocity);
	W.getInfo('i') << "initial velocity: " << "translational = " << initVelocity << ", rotational = " << initAngularVelocity << std::endl;
	//*/

}//ReadSpecificParametersFromDictionary()
#endif