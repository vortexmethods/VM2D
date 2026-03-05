/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DDeformable.cpp                                        |
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
\brief Файл кода с описанием класса MechanicsDeformable
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\author Серебровская Екатерина Александровна
\Version 1.14
\date 6 марта 2026 г.
*/

#include "Mechanics2DDeformable.h"

#include "Airfoil2D.h"
#include "Airfoil2DDeformable.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;


Beam::Beam(const World2D& W_, bool fsi_, double x0_, double L_, int R_) :
	W(W_),
	fsi(fsi_),
	rho(10000.0),
	F(0.02),
	EJ(1.1703239289446188),
	R(R_),
	x0(x0_),
	L(L_)
{
	qCoeff.resize(R);
	currentPhi.resize(R);
	currentDPhi.resize(R);
	presLastSteps.reserve(nLastSteps);
};

double Beam::phi(int n, double t) const
{
	return currentPhi[n];
}

void Beam::solveDU(int n, double dt)
{
	double cm = rho * F;
	double ck = EJ * sqr(sqr(unitLambda[n] / L));
	double cq = qCoeff[n];
	//std::cout << "ck = " << ck << std::endl;
	//std::cout << "lam = " << unitLambda[n] << std::endl;
	//std::cout << "lam = " << unitLambda[n] / L << std::endl;

	//double omega = sqrt(ck / cm);
	//double cDamp = 0.0; //0.005;

	double phiAst = currentPhi[n] + 0.5 * dt * currentDPhi[n];
	double psiAst = currentDPhi[n] + 0.5 * dt * (-ck * currentPhi[n] /*- cDamp * omega * currentDPhi[n]*/ - cq) / cm;

	currentPhi[n] += dt * psiAst;
	currentDPhi[n] += dt * (-ck * phiAst /*- cDamp * omega * psiAst*/ - cq) / cm;
}

void Beam::solveDU_RK(int n, double dt) {
	double cm = rho * F;
	double ck = EJ * sqr(sqr(unitLambda[n] / L));
	double cq = qCoeff[n];

	double omega = sqrt(ck / cm);
	double cDamp = 0.0;

	auto acceleration = [&](double phi, double dphi) {
		return (-ck * phi - cDamp * omega * dphi - cq) / cm;
		};

	double phi = currentPhi[n];
	double dphi = currentDPhi[n];

	double k1_phi = dphi;
	double k1_dphi = acceleration(phi, dphi);

	double k2_phi = dphi + 0.5 * dt * k1_dphi;
	double k2_dphi = acceleration(phi + 0.5 * dt * k1_phi, dphi + 0.5 * dt * k1_dphi);

	double k3_phi = dphi + 0.5 * dt * k2_dphi;
	double k3_dphi = acceleration(phi + 0.5 * dt * k2_phi, dphi + 0.5 * dt * k2_dphi);

	double k4_phi = dphi + dt * k3_dphi;
	double k4_dphi = acceleration(phi + dt * k3_phi, dphi + dt * k3_dphi);

	currentPhi[n] = phi + (dt / 6.0) * (k1_phi + 2 * k2_phi + 2 * k3_phi + k4_phi);
	currentDPhi[n] = dphi + (dt / 6.0) * (k1_dphi + 2 * k2_dphi + 2 * k3_dphi + k4_dphi);
	
	std::cout << "phi = " << currentPhi[n] << std::endl;
}


double Beam::getTotalDisp(double x, double t) const
{
	double result = 0.0;
	for (int i = 0; i < R; ++i)
		result += phi(i, t) * shape(i, x);
	return result;
}

double Beam::getGivenLaw(double x, double t, double deformParam) const //имитатор деформации упругой линии
{
	const double c1 = -0.825;
	const double c2 = 1.625;
	double alpha;
    if (fsi) //Turek
		alpha = 0.0; //Turek
	else
		alpha = 0.1; //Fish
	

	const double lambda = 1.0;
	const double length = 1.0;
	const double f = deformParam;

	auto A = [c1, c2](double xi) {return 1.0 + (xi - 1.0) * c1 + (xi * xi - 1.0) * c2;};

	//double result = alpha * A(x + 0.5) * sin(DPI * ((x + 0.5) / (length * lambda) - f * t));
	double result = W.getPassport().physicalProperties.accelCft(t) * alpha * A(x + 0.5) * sin(DPI * ((x + 0.5) / (length * lambda) - f * t));

	return result;
}




MechanicsDeformable::MechanicsDeformable(const World2D& W_, size_t numberInPassport_)
	:
	Mechanics(W_, numberInPassport_, true, true)
{
	const auto& afl = W_.getAirfoil(numberInPassport_);

	Vcm0 = { 0.0, 0.0 };
	Rcm0 = { afl.rcm[0], afl.rcm[1] };
	Vcm = Vcm0;
	Rcm = Rcm0;
	VcmOld = Vcm0;
	RcmOld = Rcm;

		
	ReadSpecificParametersFromDictionary();

	Point2D zero = { 0.0, 0.0 };

	Initialize(zero, afl.rcm + zero, 0.0, afl.phiAfl + 0.0);

	if (afl.phiAfl != 0)
	{
		W.getInfo('e') << "Airfoil rotation for Turek problem is not allowed" << std::endl;
		exit(2345);
	}

	if (fsi)
	{
		//Turek
		//Выделение упругой хорды
		indexOfUpperRightAngle = 0;
		double x0 = afl.getR(indexOfUpperRightAngle)[0];
		double x1 = x0;
		while (fabs(x1 - x0) < 1e-12)
		{
			++indexOfUpperRightAngle;
			x0 = x1;
			x1 = afl.getR(indexOfUpperRightAngle)[0];
		}
		--indexOfUpperRightAngle;

		upperShifts.resize(indexOfUpperRightAngle - 1);
		for (size_t i = 0; i < upperShifts.size(); ++i)
			upperShifts[i] = afl.getR(i + 1)[1] - afl.getR(0)[1];


		indexOfUpperLeftAngle = indexOfUpperRightAngle;
		double y0 = afl.getR(indexOfUpperLeftAngle)[1];
		double y1 = y0;
		while (fabs(y1 - y0) < 1e-12)
		{
			++indexOfUpperLeftAngle;
			y0 = y1;
			y1 = afl.getR(indexOfUpperLeftAngle)[1];
		}
		--indexOfUpperLeftAngle;


		indexOfLowerRightAngle = afl.getNumberOfPanels() + 1;
		x0 = afl.getR(indexOfLowerRightAngle)[0];
		x1 = x0;
		while (fabs(x1 - x0) < 1e-12)
		{
			--indexOfLowerRightAngle;
			x0 = x1;
			x1 = afl.getR(indexOfLowerRightAngle)[0];
		}
		++indexOfLowerRightAngle;


		lowerShifts.resize(afl.getNumberOfPanels() - indexOfLowerRightAngle - 1);
		for (size_t i = 0; i < lowerShifts.size(); ++i)
			lowerShifts[i] = afl.getR(afl.getNumberOfPanels() - 1 - i)[1] - afl.getR(0)[1];

		indexOfLowerLeftAngle = indexOfLowerRightAngle;
		y0 = afl.getR(indexOfLowerLeftAngle)[1];
		y1 = y0;
		while (fabs(y1 - y0) < 1e-12)
		{
			--indexOfLowerLeftAngle;
			y0 = y1;
			y1 = afl.getR(indexOfLowerLeftAngle)[1];
		}
		++indexOfLowerLeftAngle;

		//W.getInfo('i') << "UR: " << afl.getR(indexOfUpperRightAngle) << std::endl;
		//W.getInfo('i') << "UL: " << afl.getR(indexOfUpperLeftAngle) << std::endl;
		//W.getInfo('i') << "LR: " << afl.getR(indexOfLowerRightAngle) << std::endl;
		//W.getInfo('i') << "LL: " << afl.getR(indexOfLowerLeftAngle) << std::endl;

		if (indexOfUpperLeftAngle - indexOfUpperRightAngle != indexOfLowerRightAngle - indexOfLowerLeftAngle)
		{
			W.getInfo('e') << "indexOfUpperLeftAngle - indexOfUpperRightAngle != indexOfLowerRightAngle - indexOfLowerLeftAngle" << std::endl;
			exit(2346);
		}

		chord.resize(indexOfUpperLeftAngle - indexOfUpperRightAngle);

		for (size_t i = 0; i < indexOfUpperLeftAngle - indexOfUpperRightAngle; ++i)
		{
			size_t idxUp = indexOfUpperLeftAngle - 1 - i;
			size_t idxDn = indexOfLowerLeftAngle + i;
			Point2D rUpLeft = afl.getR(idxUp + 1);
			Point2D rUpRight = afl.getR(idxUp);
			Point2D rDnLeft = afl.getR(idxDn);
			Point2D rDnRight = afl.getR(idxDn + 1);

			if ((fabs(rUpLeft[0] - rDnLeft[0]) > 0.01 * (rUpRight[0] - rUpLeft[0])) || (fabs(rUpLeft[0] - rDnLeft[0]) > 0.01 * (rDnRight[0] - rDnLeft[0])) ||
				(fabs(rUpRight[0] - rDnRight[0]) > 0.01 * (rUpRight[0] - rUpLeft[0])) || (fabs(rUpLeft[0] - rDnLeft[0]) > 0.01 * (rDnRight[0] - rDnLeft[0])))
			{
				W.getInfo('e') << "x_up != x_dn" << std::endl;
				exit(2347);
			}

			chord[i].beg = 0.5 * (rUpLeft + rDnLeft);
			chord[i].end = 0.5 * (rUpRight + rDnRight);
			chord[i].infPanels = { idxUp, idxDn };
			chord[i].rightSemiWidth = 0.5 * (rUpRight - rDnRight)[1];
		}
		initialChord = chord;
		beam = std::make_unique<Beam>(W, fsi, chord[0].beg[0], initialChord.back().end[0] - initialChord[0].beg[0], 3); //Beam
	}
	else //Fish
	{
		int np = (int)afl.getNumberOfPanels();
		chord.resize(np / 2);

		chord[0].beg = afl.getR(np / 2);
		chord[np / 2 - 1].end = afl.getR(0);
		chord[np / 2 - 1].rightSemiWidth = 0.0;

		for (size_t i = 0; i < np / 2; ++i)
		{
			if (i != 0)
				chord[i].beg = 0.5 * (afl.getR(np / 2 - i) + afl.getR(np / 2 + i));

			if (i != np / 2 - 1)
				chord[i].end = 0.5 * (afl.getR(np / 2 - i - 1) + afl.getR(np / 2 + i + 1));

			chord[i].infPanels = { np / 2 - i, np / 2 + i + 1 };

			if (i != np / 2 - 1)
				chord[i].rightSemiWidth = (afl.getR(np / 2 - i - 1) - afl.getR(np / 2 + i + 1)).length() * 0.5;
		}
		initialChord = chord;
	}
	
	//std::ofstream of("chord.txt");		
	//for (size_t i = 0; i < chord.size(); ++i)
	//	of << chord[i].beg[0] << " " << chord[i].beg[1] << " " << chord[i].end[0] << " " << chord[i].end[1] << " " << chord[i].infPanels.first << " " << chord[i].infPanels.second << std::endl;
	//of.close();	

	

	initialPossibleWays = afl.possibleWays;
};

//Вычисление гидродинамической силы, действующей на профиль
void MechanicsDeformable::GetHydroDynamForce()
{
	W.getTimers().start("Force");	

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
		if (W.getCurrentStep() > 0)
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

	W.getTimers().stop("Force");
}// GetHydroDynamForce()

// Вычисление скорости центра масс
Point2D MechanicsDeformable::VeloOfAirfoilRcm(double currTime)
{
	return Vcm;
}//VeloOfAirfoilRcm(...)

// Вычисление положения центра масс
Point2D MechanicsDeformable::PositionOfAirfoilRcm(double currTime)
{
	return Rcm;
}//PositionOfAirfoilRcm(...)

double MechanicsDeformable::AngularVelocityOfAirfoil(double currTime)
{
	return Wcm;
}//AngularVelocityOfAirfoil(...)

double MechanicsDeformable::AngleOfAirfoil(double currTime)
{
	if (afl.phiAfl != Phi)
	{
		std::cout << "afl.phiAfl != Phi" << std::endl;
		exit(100600);
	}

	return afl.phiAfl;
}//AngleOfAirfoil(...)

// Вычисление скоростей начал панелей
void MechanicsDeformable::VeloOfAirfoilPanels(double currTime)
{
	std::vector<Point2D> veloW(afl.getNumberOfPanels(), {0.0, 0.0});

	//if (W.getCurrentStep() == 0)
	//	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	//		veloW[i] = { 0.0, 0.0 };//afl.getR(i).kcross();

	if (W.getCurrentStep() > 0)
		for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
			veloW[i] = (1.0 / W.getPassport().timeDiscretizationProperties.dt) * (afl.getR(i) - W.getOldAirfoil(0).getR(i));

	afl.setV(veloW);


	//Циркуляция
	circulationOld = circulation;
	circulation = 0.0;
	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
		circulation += 0.5 * afl.len[i] * ((afl.getV(i) + afl.getV(i + 1)) & afl.tau[i]);

}//VeloOfAirfoilPanels(...)


void MechanicsDeformable::Move()
{
	double x0 = chord[0].beg[0];
	double t = W.getCurrentTime();


	if (fsi) //turek
	{

		for (int i = 0; i < beam->R; ++i)
			beam->solveDU(i, W.getPassport().timeDiscretizationProperties.dt);

		for (size_t i = 0; i < chord.size(); ++i)
		{
			chord[i].beg[1] = beam->getTotalDisp(chord[i].beg[0], t);
			chord[i].end[1] = beam->getTotalDisp(chord[i].end[0], t);
		}

		std::vector<Point2D> upperPoints(chord.size() + upperShifts.size());
		std::vector<Point2D> lowerPoints(chord.size() + lowerShifts.size());

		for (size_t i = 0; i < chord.size() - 1; ++i)
		{
			const Point2D& begIp1 = chord[i + 1].beg;
			const Point2D& endI = chord[i].end;
			Point2D normI = (endI - chord[i].beg).unit().kcross();
			Point2D normIp1 = (chord[i + 1].end - begIp1).unit().kcross();
			upperPoints[i] = endI + 0.5 * chord[i].rightSemiWidth * (normI + normIp1);
			lowerPoints[i] = endI - 0.5 * chord[i].rightSemiWidth * (normI + normIp1);
		}
		Point2D normBack = (chord.back().end - chord.back().beg).unit().kcross();
		upperPoints[chord.size() - 1] = chord.back().end + chord.back().rightSemiWidth * normBack;
		lowerPoints[chord.size() - 1] = chord.back().end - chord.back().rightSemiWidth * normBack;

		for (size_t i = 0; i < upperShifts.size(); ++i)
			upperPoints[chord.size() + i] = chord.back().end + upperShifts[upperShifts.size() - 1 - i] * normBack;
		for (size_t i = 0; i < lowerShifts.size(); ++i)
			lowerPoints[chord.size() + i] = chord.back().end + lowerShifts[lowerShifts.size() - 1 - i] * normBack;

		afl.setR(0) = { chord.back().end[0], beam->getTotalDisp(chord.back().end[0], t) };
		for (size_t i = 0; i < upperPoints.size(); ++i)
			afl.setR(indexOfUpperLeftAngle - 1 - i) = upperPoints[i];
		for (size_t i = 0; i < lowerPoints.size(); ++i)
			afl.setR(indexOfLowerLeftAngle + 1 + i) = lowerPoints[i];
	}
	else //Fish
	{
		for (size_t i = 0; i < chord.size(); ++i)
		{
			chord[i].beg[1] = beam->getGivenLaw(chord[i].beg[0], t, deformParam);
			chord[i].end[1] = beam->getGivenLaw(chord[i].end[0], t, deformParam);
		}

		std::vector<Point2D> upperPoints(chord.size());
		std::vector<Point2D> lowerPoints(chord.size());

		for (size_t i = 0; i < chord.size() - 1; ++i)
		{
			const Point2D& begIp1 = chord[i + 1].beg;
			const Point2D& endI = chord[i].end;
			Point2D normI = (endI - chord[i].beg).unit().kcross();
			Point2D normIp1 = (chord[i + 1].end - begIp1).unit().kcross();
			upperPoints[i] = endI + 0.5 * chord[i].rightSemiWidth * (normI + normIp1);
			lowerPoints[i] = endI - 0.5 * chord[i].rightSemiWidth * (normI + normIp1);
		}

		Point2D normFront = (chord[0].end - chord[0].beg).unit().kcross();
		Point2D normBack = (chord.back().end - chord.back().beg).unit().kcross();

		int nph = (int)chord.size();

		afl.setR(0) = chord.back().end;
		afl.setR(nph) = chord[0].beg;
		for (size_t i = 0; i < nph - 1; ++i)
		{
			afl.setR(nph - 1 - i) = upperPoints[i];
			afl.setR(nph + 1 + i) = lowerPoints[i];
		}
	}
	
	afl.CalcNrmTauLen();
	afl.GetGabarits();

	afl.possibleWays.clear();
	afl.possibleWays.resize(initialPossibleWays.size());

	for (size_t w = 0; w < initialPossibleWays.size(); ++w)
	{
		std::vector<Point2D> initialWay = initialPossibleWays[w];
		afl.possibleWays[w].resize(initialWay.size());
		for (size_t p = 0; p < initialWay.size(); ++p)
		{
			double x = initialWay[p][0];
			double y;

			if (fsi) //turek
				y = beam->getTotalDisp(x, t); 
			else //fish
				y = beam->getGivenLaw(x, t, deformParam); 
			
			afl.possibleWays[w][p] = Point2D{ x,y };
		}	
	}

	afl.lightningTest();
}//Move()







#if defined(INITIAL) || defined(BRIDGE) 
void MechanicsDeformable::ReadSpecificParametersFromDictionary()
{
	mechParamsParser->get("deformParam", deformParam);
	W.getInfo('i') << "mechDeformable: " << "deformParam = " << deformParam << std::endl;

	mechParamsParser->get("fsi", fsi);
	W.getInfo('i') << "mechDeformable: " << "fsi = " << fsi << std::endl;
}//ReadSpecificParametersFromDictionary()
#endif