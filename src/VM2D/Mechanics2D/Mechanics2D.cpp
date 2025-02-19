/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Mechanics2D.cpp                                                  |
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
\brief Файл кода с описанием класса Mechanics.cpp
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#include "Mechanics2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"



using namespace VM2D;

//Конструктор
Mechanics::Mechanics(const World2D& W_, size_t numberInPassport_, bool isMoves_, bool isDeform_)
	:
	W(W_),
	numberInPassport(numberInPassport_),
	afl(W_.getNonConstAirfoil(numberInPassport_)),
	boundary(W_.getBoundary(numberInPassport_)),
	virtVortParams(W_.getVelocity().virtualVortexesParams[numberInPassport_]),

	Vcm0({ 0.0, 0.0 }),
	Wcm0(0.0),
	Rcm0({ 0.0, 0.0 }),
	Phi0(0.0),

	/*Vcm0(Vcm0_),
	Wcm0(Wcm0_),
	Rcm0(Rcm0_),
	Phi0(Phi0_),*/

	isMoves(isMoves_), 
	isDeform(isDeform_), 

	hydroDynamForce({0.0, 0.0}),
	hydroDynamMoment(0.0),
	viscousForce({ 0.0, 0.0 }),
	viscousMoment(0.0)

{
	ReadParametersFromDictionary();
}; 

//Задание начального положения и начальной скорости
void Mechanics::Initialize(Point2D Vcm0_, Point2D Rcm0_, double Wcm0_, double Phi0_)
{
	Vcm0 = Vcm0_;
	Wcm0 = Wcm0_;
	Rcm0 = Rcm0_;
	Phi0 = Phi0_;

	Vcm = Vcm0;
	Wcm = Wcm0;
	Rcm = Rcm0;
	Phi = Phi0;
	VcmOld = Vcm0;
	WcmOld = Wcm0;
	RcmOld = Rcm0;
	PhiOld = Phi0;

	circulation = 2.0 * afl.area * Wcm;
	circulationOld = circulation;

	afl.Move(Rcm - W.getAirfoil(numberInPassport).rcm);
	afl.Rotate(Phi - W.getAirfoil(numberInPassport).phiAfl);
}//Initialize(...)

//Парсинг списка параметров механической системы
void Mechanics::ReadParametersFromDictionary()
{
	std::stringstream ss(W.getPassport().airfoilParams[afl.numberInPassport].mechanicalSystemParameters);
	
	mechParamsParser.reset(new VMlib::StreamParser(W.getInfo(), "mechanical parser", ss));
}//ReadParametersFromDictionary()


// Генерация заголовка файла нагрузок
void Mechanics::GenerateForcesHeader()
{
	std::stringstream forceFileName, forceFileNameCsv;
	forceFileName << W.getPassport().dir << "forces-airfoil-" << numberInPassport;
	forceFileNameCsv << W.getPassport().dir << "forces-airfoil-" << numberInPassport << ".csv";

	std::ofstream newForcesFile(forceFileName.str());
	std::ofstream newForcesFileCsv(forceFileNameCsv.str());

	VMlib::PrintLogoToTextFile(newForcesFile, forceFileName.str(), "Hydrodynamic loads for the airfoil " + W.getPassport().airfoilParams[numberInPassport].fileAirfoil);

	if (!W.getPassport().calcCoefficients)
		VMlib::PrintHeaderToTextFile(newForcesFile, "currentStep     currentTime     Fx     Fy     Mz     Ftaux     Ftauy     Mtau");
	else
		VMlib::PrintHeaderToTextFile(newForcesFile, "currentStep     currentTime     CFx    CFy    CMz    CFtaux    CFtauy    CMtau");

	if (!W.getPassport().calcCoefficients)
		newForcesFileCsv << "step,time,Fx,Fy,Mz,Ftaux,Ftauy,Mtau";
	else
		newForcesFileCsv << "step,time,CFx,CFy,CMz,CFtaux,CFtauy,CMtau";

	newForcesFile.close();
	newForcesFile.clear();

	newForcesFileCsv.close();
	newForcesFileCsv.clear();

}//GenerateForcesHeader()


//Генерация заголовка файла положения профиля
void Mechanics::GeneratePositionHeader()
{
	if (isMoves)
	{
		std::stringstream positionFileName, positionFileNameCsv;
		positionFileName << W.getPassport().dir << "position-airfoil-" << numberInPassport;
		positionFileNameCsv << W.getPassport().dir << "position-airfoil-" << numberInPassport << ".csv";

		std::ofstream newPositionFile(positionFileName.str());
		std::ofstream newPositionFileCsv(positionFileNameCsv.str());

		VMlib::PrintLogoToTextFile(newPositionFile, positionFileName.str(), "Position of the airfoil " + W.getPassport().airfoilParams[numberInPassport].fileAirfoil);

		VMlib::PrintHeaderToTextFile(newPositionFile, "currentStep     currentTime     x     y     phi     Vx     Vy     w");

		newPositionFileCsv << "step,time,x,y,phi,Vx,Vy,w";

		newPositionFile.close();
		newPositionFile.clear();

		newPositionFileCsv.close();
		newPositionFileCsv.clear();
	}
}//GeneratePositionHeader()


//Сохранение строки со статистикой в файл нагрузок
void Mechanics::GenerateForcesString()
{
	W.getTimestat().timeOther.first += omp_get_wtime();

	std::stringstream forceFileName, forceFileNameCsv;
	forceFileName << W.getPassport().dir << "forces-airfoil-" << numberInPassport;
	forceFileNameCsv << W.getPassport().dir << "forces-airfoil-" << numberInPassport << ".csv";

	//double cShock = (W.getPassport().physicalProperties.getCurrTime() > W.getPassport().physicalProperties.timeAccel + 2.0 * W.getPassport().timeDiscretizationProperties.dt) ? 1.0 : 0.0;
	double cShock = 1.0;

	Point2D cartesianHydroForce = cShock * hydroDynamForce;
	Point2D cartesianViscoForce = cShock * viscousForce;

	double cartesianHydroMoment = cShock * hydroDynamMoment;
	double cartesianViscoMoment = cShock * viscousMoment;

	double rotationAngle = W.getPassport().airfoilParams[numberInPassport].angle;
	if (W.getPassport().geographicalAngles)
	{
		cartesianHydroForce = cartesianHydroForce.rotated(-0.5 * PI - W.getPassport().airfoilParams[numberInPassport].angle);
		cartesianViscoForce = cartesianViscoForce.rotated(-0.5 * PI - W.getPassport().airfoilParams[numberInPassport].angle);
		//rotationAngle += 0.5 * PI;
	}

	//if (W.getPassport().rotateForces)
	//{
	//	cartesianHydroForce = cartesianHydroForce.rotated(rotationAngle);
	//	cartesianViscoForce = cartesianViscoForce.rotated(rotationAngle);
	//}

	double flowPress = 0.5 * W.getPassport().physicalProperties.rho * sqr(W.getPassport().physicalProperties.vRef);
	if (W.getPassport().calcCoefficients)
	{
		const double& chord = W.getPassport().airfoilParams[numberInPassport].chord;
		cartesianHydroForce /= (flowPress * chord);
		cartesianViscoForce /= (flowPress * chord);

		cartesianHydroMoment /= (flowPress * sqr(chord));
		cartesianViscoMoment /= (flowPress * sqr(chord));
	}


	std::ofstream forcesFile(forceFileName.str(), std::ios::app);
	forcesFile << std::endl << W.getCurrentStep() << "	" << W.getPassport().physicalProperties.getCurrTime() << "	" << cartesianHydroForce[0] << "	" << cartesianHydroForce[1] << "	" << cartesianHydroMoment << "	" << cartesianViscoForce[0] << "	" << cartesianViscoForce[1] << "	" << cartesianViscoMoment;
	forcesFile.close();

	std::ofstream forcesFileCsv(forceFileNameCsv.str(), std::ios::app);
	forcesFileCsv << std::endl << W.getCurrentStep() << "," << W.getPassport().physicalProperties.getCurrTime() << "," << cartesianHydroForce[0] << "," << cartesianHydroForce[1] << "," << cartesianHydroMoment << "," << cartesianViscoForce[0] << "," << cartesianViscoForce[1] << "," << cartesianViscoMoment;
	forcesFileCsv.close();

	W.getTimestat().timeOther.second += omp_get_wtime();
}//GenerateForcesString()


//Сохранение строки со статистикой в файл положения
void Mechanics::GeneratePositionString()
{
	W.getTimestat().timeOther.first += omp_get_wtime();

	if (isMoves)
	{
		std::stringstream positionFileName, positionFileNameCsv;
		positionFileName << W.getPassport().dir << "position-airfoil-" << numberInPassport;
		positionFileNameCsv << W.getPassport().dir << "position-airfoil-" << numberInPassport << ".csv";

		std::ofstream positionFile(positionFileName.str(), std::ios::app);
		positionFile << std::endl << W.getCurrentStep() << "	" << W.getPassport().physicalProperties.getCurrTime() << "	" << afl.rcm[0] << "	" << afl.rcm[1] << "	" << Phi << "	" << Vcm[0] << "	" << Vcm[1] << "	" << Wcm;
		positionFile.close();

		std::ofstream positionFileCsv(positionFileNameCsv.str(), std::ios::app);
		positionFileCsv << std::endl << W.getCurrentStep() << "," << W.getPassport().physicalProperties.getCurrTime() << "," << afl.rcm[0] << "," << afl.rcm[1] << "," << Phi << "," << Vcm[0] << "," << Vcm[1] << "," << Wcm;
		positionFileCsv.close();
	}

	W.getTimestat().timeOther.second += omp_get_wtime();
}//GeneratePositionString()