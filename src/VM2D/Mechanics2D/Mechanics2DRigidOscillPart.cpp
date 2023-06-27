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
	Mechanics(W_, numberInPassport_, 0, true, false, false), 
	Vx0(0.0), 
	Vy0(0.0), 
	x0(W_.getAirfoil(numberInPassport_).rcm[0]), 
	y0(W_.getAirfoil(numberInPassport_).rcm[1]) 
	//bx(0.0 * 0.731),
	//by(0.0 * 0.731)
{
	Vx = Vx0;
	Vy = Vy0;

	x = x0;
	y = y0;

	VxOld = Vx0;
	VyOld = Vy0;
	
	xOld = x0;
	yOld = y0;
		
	ReadSpecificParametersFromDictionary();
	Initialize({ 0.0, 0.0 }, W_.getAirfoil(numberInPassport_).rcm, 0.0, W_.getAirfoil(numberInPassport_).phiAfl);
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
	
	Point2D hDFGam = { 0.0, 0.0 };	//гидродинамические силы, обусловленные Gamma_k
	Point2D hDFdelta = { 0.0, 0.0 };	//гидродинамические силы, обусловленные delta_k
	double hDMdelta = 0.0;


	Point2D deltaVstep;
	if (W.getPassport().airfoilParams[numberInPassport].addedMass.length2() > 0)
		deltaVstep = { 0.0, 0.0 };  //Для итерационной процедуры
	else
		deltaVstep = { Vx - VxOld, Vy - VyOld }; //Для безытерационной процедуры	

	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		/// \todo Учитываем только нулевой момент решения. Надо ли учитывать остальные?
		double deltaK = boundary.sheets.freeVortexSheet(i, 0) * afl.len[i] - afl.gammaThrough[i] + (deltaVstep & afl.tau[i]) * afl.len[i];
		Point2D rK = 0.5 * (afl.getR(i + 1) + afl.getR(i)) - afl.rcm;

		hDFdelta += deltaK * Point2D({ -rK[1], rK[0] });
		hDMdelta += 0.5 * deltaK * rK.length2();
	}

	const double rho = W.getPassport().physicalProperties.rho;

	hydroDynamForce = (rho / dt) * hDFdelta;
	hydroDynamMoment = (rho / dt) * hDMdelta;

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
	return { Vx, Vy};
	//return{ 0.0, 0.0 };
}//VeloOfAirfoilRcm(...)

// Вычисление положения центра масс
Point2D MechanicsRigidOscillPart::PositionOfAirfoilRcm(double currTime)
{
	return{ x, y };
	//return{ 0.0, 0.0 };
}//PositionOfAirfoilRcm(...)

double MechanicsRigidOscillPart::AngularVelocityOfAirfoil(double currTime)
{
	return 0.0;
}//AngularVelocityOfAirfoil(...)

double MechanicsRigidOscillPart::AngleOfAirfoil(double currTime)
{
	return afl.phiAfl;
}//AngleOfAirfoil(...)

// Вычисление скоростей начал панелей
void MechanicsRigidOscillPart::VeloOfAirfoilPanels(double currTime)
{
	Point2D veloRcm = VeloOfAirfoilRcm(currTime);
	afl.setV(veloRcm);	
}//VeloOfAirfoilPanels(...)


void MechanicsRigidOscillPart::Move()
{
	Point2D meff;
	if (W.getPassport().airfoilParams[numberInPassport].addedMass.length2() > 0)
		meff = Point2D{ m + W.getPassport().airfoilParams[numberInPassport].addedMass[0], m + W.getPassport().airfoilParams[numberInPassport].addedMass[1] };
	else
		meff = Point2D{ m, m };

	VxOld = Vx;
	VyOld = Vy;

	xOld = x;
	yOld = y;

	double dx, dVx, dy, dVy;

	//W.getInfo('t') << "k = " << k << std::endl;


	if (ky > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		numvector<double, 2> kk[4];
		kk[0] = { Vy, (hydroDynamForce[1] - 2.0 * by * Vy - ky * y) / meff[1]};
		kk[1] = { Vy + 0.5 * dt * kk[0][1], (hydroDynamForce[1] - 2.0 * by * (Vy + 0.5 * dt * kk[0][1]) - ky * (y + 0.5 * dt * kk[0][0])) / meff[1]};
		kk[2] = { Vy + 0.5 * dt * kk[1][1], (hydroDynamForce[1] - 2.0 * by * (Vy + 0.5 * dt * kk[1][1]) - ky * (y + 0.5 * dt * kk[1][0])) / meff[1
		]};
		kk[3] = { Vy + dt * kk[2][1], (hydroDynamForce[1] - 2.0 * by * (Vy + dt * kk[2][1]) - ky * (y + dt * kk[2][0])) / meff[1]};

		dy = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		dVy = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dy = 0.0;
		dVy = 0.0;
	}


	if (kx > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		numvector<double, 2> kk[4];
		kk[0] = { Vx, (hydroDynamForce[0] - 2.0 * bx * Vx - kx * x) / meff[0]};
		kk[1] = { Vx + 0.5 * dt * kk[0][1], (hydroDynamForce[0] - 2.0 * bx * (Vx + 0.5 * dt * kk[0][1]) - kx * (x + 0.5 * dt * kk[0][0])) / meff[0]};
		kk[2] = { Vx + 0.5 * dt * kk[1][1], (hydroDynamForce[0] - 2.0 * bx * (Vx + 0.5 * dt * kk[1][1]) - kx * (x + 0.5 * dt * kk[1][0])) / meff[0]};
		kk[3] = { Vx + dt * kk[2][1], (hydroDynamForce[0] - 2.0 * bx * (Vx + dt * kk[2][1]) - kx * (x + dt * kk[2][0])) / meff[0]};

		dx = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		dVx = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dx = 0.0;
		dVx = 0.0;
	}

	afl.Move({ dx, dy });
	Vcm[0] += dVx;
	Vcm[1] += dVy;

	x += dx;
	y += dy;

	Vx += dVx;
	Vy += dVy;

}//Move()

void MechanicsRigidOscillPart::UpdateU()
{
	Point2D meff;
	if (W.getPassport().airfoilParams[numberInPassport].addedMass.length2() > 0)
		//#ifdef addm 
		meff = Point2D{ m + W.getPassport().airfoilParams[numberInPassport].addedMass[0], m + W.getPassport().airfoilParams[numberInPassport].addedMass[1] };
	//#else
	else
		meff = Point2D{ m, m };
	//double meff = m;
	//#endif



	VxOld = Vx;
	VyOld = Vy;

	xOld = x;
	yOld = y;

	double dx, dVx, dy, dVy;

	//W.getInfo('t') << "k = " << k << std::endl;


	if (ky > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		numvector<double, 2> kk[4];
		kk[0] = { Vy, (hydroDynamForce[1] - 2.0 * by * Vy - ky * y) / meff[1]};
		kk[1] = { Vy + 0.5 * dt * kk[0][1], (hydroDynamForce[1] - 2.0 * by * (Vy + 0.5 * dt * kk[0][1]) - ky * (y + 0.5 * dt * kk[0][0])) / meff[1]};
		kk[2] = { Vy + 0.5 * dt * kk[1][1], (hydroDynamForce[1] - 2.0 * by * (Vy + 0.5 * dt * kk[1][1]) - ky * (y + 0.5 * dt * kk[1][0])) / meff[1]};
		kk[3] = { Vy + dt * kk[2][1], (hydroDynamForce[1] - 2.0 * by * (Vy + dt * kk[2][1]) - ky * (y + dt * kk[2][0])) / meff[1]};

		dy = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		dVy = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dy = 0.0;
		dVy = 0.0;
	}


	if (kx > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		numvector<double, 2> kk[4];
		kk[0] = { Vx, (hydroDynamForce[0] - 2.0 * bx * Vx - kx * x) / meff[0]};
		kk[1] = { Vx + 0.5 * dt * kk[0][1], (hydroDynamForce[0] - 2.0 * bx * (Vx + 0.5 * dt * kk[0][1]) - kx * (x + 0.5 * dt * kk[0][0])) / meff[0]};
		kk[2] = { Vx + 0.5 * dt * kk[1][1], (hydroDynamForce[0] - 2.0 * bx * (Vx + 0.5 * dt * kk[1][1]) - kx * (x + 0.5 * dt * kk[1][0])) / meff[0]};
		kk[3] = { Vx + dt * kk[2][1], (hydroDynamForce[0] - 2.0 * bx * (Vx + dt * kk[2][1]) - kx * (x + dt * kk[2][0])) / meff[0]};

		dx = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		dVx = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dx = 0.0;
		dVx = 0.0;
	}

	//afl.Move({ dx, dy });
	Vcm[0] += dVx;
	Vcm[1] += dVy;

	//x += dx;
	//y += dy;

	Vx += dVx;
	Vy += dVy;

}//Move()


/*
void MechanicsRigidOscillPart::RecalcU(Point2D forcePrev) //ИК
{
	double dx, dVx, dy, dVy;

#ifdef addm 
	double meff = m;// +W.getPassport().physicalProperties.rho * PI * 0.5 * 0.5;
#else
	double meff = m;
#endif

		//W.getInfo('t') << "k = " << k << std::endl;

	Point2D force = 1.0 * forcePrev;// +0.5 * Qiter;

	if (ky > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		numvector<double, 2> kk[4];
		kk[0] = { Vy, (force[1] - 2.0*by * Vy - ky * y) / meff };
		kk[1] = { Vy + 0.5 * dt * kk[0][1], (force[1] - 2.0 * by * (Vy + 0.5 * dt * kk[0][1]) - ky * (y + 0.5 * dt * kk[0][0])) / meff };
		kk[2] = { Vy + 0.5 * dt * kk[1][1], (force[1] - 2.0 * by * (Vy + 0.5 * dt * kk[1][1]) - ky * (y + 0.5 * dt * kk[1][0])) / meff };
		kk[3] = { Vy + dt * kk[2][1], (force[1] - 2.0 * by * (Vy + dt * kk[2][1]) - ky * (y + dt * kk[2][0])) / meff };

		dy = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		dVy = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dy = 0.0;
		dVy = 0.0;
	}


	if (kx > 0)
	{
		double dt = W.getPassport().timeDiscretizationProperties.dt;
		numvector<double, 2> kk[4];
		kk[0] = { Vx, (force[0] - 2.0 * bx * Vx - kx * x) / meff };
		kk[1] = { Vx + 0.5 * dt * kk[0][1], (force[0] - 2.0 * bx * (Vx + 0.5 * dt * kk[0][1]) - kx * (x + 0.5 * dt * kk[0][0])) / meff };
		kk[2] = { Vx + 0.5 * dt * kk[1][1], (force[0] - 2.0 * bx * (Vx + 0.5 * dt * kk[1][1]) - kx * (x + 0.5 * dt * kk[1][0])) / meff };
		kk[3] = { Vx + dt * kk[2][1], (force[0] - 2.0 * bx * (Vx + dt * kk[2][1]) - kx * (x + dt * kk[2][0])) / meff };

		dx = dt * (kk[0][0] + 2. * kk[1][0] + 2. * kk[2][0] + kk[3][0]) / 6.0;
		dVx = dt * (kk[0][1] + 2. * kk[1][1] + 2. * kk[2][1] + kk[3][1]) / 6.0;
	}
	else
	{
		dx = 0.0;
		dVx = 0.0;
	}



	//afl.Move({ 0.0, dy });
	//Vcm[1] += du;

	//y += dy;
	Vx = VxIter + 0.5 * dVx;
	Vy = VyIter + 0.5 * dVy;
}//RecalcU()
*/


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

//template <typename T>
//inline T sqr(T x) 
//{ return x * x; }

#ifdef INITIAL
void MechanicsRigidOscillPart::ReadSpecificParametersFromDictionary()
{
	mechParamsParser->get("m", m);
	
	W.getInfo('i') << "mass " << "m = " << m << std::endl;
		
	std::vector<double> sh;
	mechParamsParser->get("sh", sh);	
	mechParamsParser->get("sh", sh);	
	


	kx = m * sqr(2.0 * PI * sh[0] / W.getPassport().airfoilParams[numberInPassport].chord) * W.getPassport().physicalProperties.vInf.length2();
	ky = m * sqr(2.0 * PI * sh[1] / W.getPassport().airfoilParams[numberInPassport].chord) * W.getPassport().physicalProperties.vInf.length2();

	W.getInfo('i') << "rigidity kx = " << kx << std::endl;
	W.getInfo('i') << "rigidity ky = " << ky << std::endl;

	std::vector<double> zeta;
	mechParamsParser->get("zeta", zeta);//Логарифмический декремент

	bx = zeta[0] / (2.0 * PI) * sqrt(kx * m);
	by = zeta[1] / (2.0 * PI) * sqrt(ky * m);

	W.getInfo('i') << "damping bx = " << bx << std::endl;
	W.getInfo('i') << "damping by = " << by << std::endl;


}//ReadSpecificParametersFromDictionary()
#endif