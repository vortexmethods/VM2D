/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Mechanics.cpp                                                    |
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
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#include "Mechanics.h"
#include "World2D.h"

//Конструктор
Mechanics::Mechanics(const World2D& W_, size_t numberInPassport_, int degOfFreedom_, bool isMoves_, bool isDeform_, bool isRotate_, Point2D Vcm0_, Point2D Rcm0_, double Wcm0_, double Phi0_)
: W(W_), numberInPassport(numberInPassport_), afl(W_.getNonConstAirfoil(numberInPassport_)), boundary(W_.getBoundary(numberInPassport_)), virtVortParams(W_.getVelocity().virtualVortexesParams[numberInPassport_]), degOfFreedom(degOfFreedom_), isMoves(isMoves_), isDeform(isDeform_), isRotate(isRotate_),\
Vcm0(Vcm0_), Rcm0(Rcm0_), Wcm0(Wcm0_), Phi0(Phi0_)
{
	Vcm = Vcm0;
	Wcm = Wcm0;
	Rcm = Rcm0;
	Phi = Phi0;
	VcmOld = Vcm0;
	WcmOld = Wcm0;
	RcmOld = Rcm0;
	PhiOld = Phi0;
	ReadParametersFromDictionary();
}; 


//Парсинг списка параметров механической системы
void Mechanics::ReadParametersFromDictionary()
{
	std::stringstream ss(W.getPassport().airfoilParams[afl.numberInPassport].mechanicalSystemParameters);
	
	mechParamsParser.reset(new StreamParser(W.getInfo(), "mechanical parser", ss));
}//ReadParametersFromDictionary()


// Генерация заголовка файла нагрузок
void Mechanics::GenerateForcesHeader()
{
	std::stringstream forceFileName, forceFileNameCsv;
	forceFileName << W.getPassport().dir << "forces-airfoil-" << numberInPassport;
	forceFileNameCsv << W.getPassport().dir << "forces-airfoil-" << numberInPassport << ".csv";

	std::ofstream newForcesFile(forceFileName.str());
	std::ofstream newForcesFileCsv(forceFileNameCsv.str());

	PrintLogoToTextFile(newForcesFile, forceFileName.str(), "Hydrodynamic loads for the airfoil " + W.getPassport().airfoilParams[numberInPassport].fileAirfoil);

	PrintHeaderToTextFile(newForcesFile, "currentStep     currentTime     Fx     Fy     Mz");

	newForcesFileCsv << "t,Fx,Fy,Mz";

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

		PrintLogoToTextFile(newPositionFile, positionFileName.str(), "Position of the airfoil " + W.getPassport().airfoilParams[numberInPassport].fileAirfoil);

		PrintHeaderToTextFile(newPositionFile, "currentStep     currentTime     x     y     phi     Vx     Vy     w");

		newPositionFileCsv << "t,x,y,phi,Vx,Vy,w";

		newPositionFile.close();
		newPositionFile.clear();

		newPositionFileCsv.close();
		newPositionFileCsv.clear();
	}
}//GeneratePositionHeader()


//Сохранение строки со статистикой в файл нагрузок
void Mechanics::GenerateForcesString()
{
	std::stringstream forceFileName, forceFileNameCsv;
	forceFileName << W.getPassport().dir << "forces-airfoil-" << numberInPassport;
	forceFileNameCsv << W.getPassport().dir << "forces-airfoil-" << numberInPassport << ".csv";

	std::ofstream forcesFile(forceFileName.str(), std::ios::app);
	forcesFile << std::endl << W.getCurrentStep() << "	" << W.getPassport().physicalProperties.getCurrTime() << "	" << hydroDynamForce[0] << "	" << hydroDynamForce[1] << "	" << hydroDynamMoment;
	forcesFile.close();

	std::ofstream forcesFileCsv(forceFileNameCsv.str(), std::ios::app);
	forcesFileCsv << std::endl << W.getPassport().physicalProperties.getCurrTime() << "," << hydroDynamForce[0] << "," << hydroDynamForce[1] << "," << hydroDynamMoment;
	forcesFileCsv.close();
}//GenerateForcesString()


//Сохранение строки со статистикой в файл положения
void Mechanics::GeneratePositionString()
{
	if (isMoves)
	{
		std::stringstream positionFileName, positionFileNameCsv;
		positionFileName << W.getPassport().dir << "position-airfoil-" << numberInPassport;
		positionFileNameCsv << W.getPassport().dir << "position-airfoil-" << numberInPassport << ".csv";

		std::ofstream forcesFile(positionFileName.str(), std::ios::app);
		forcesFile << std::endl << W.getCurrentStep() << "	" << W.getPassport().physicalProperties.getCurrTime() << "	" << afl.rcm[0] << "	" << afl.rcm[1] << "	" << Phi << "	" << Vcm[0] << "	" << Vcm[1] << "	" << Wcm;
		forcesFile.close();

		std::ofstream forcesFileCsv(positionFileNameCsv.str(), std::ios::app);
		forcesFileCsv << std::endl << W.getPassport().physicalProperties.getCurrTime() << "," << afl.rcm[0] << "," << afl.rcm[1] << "," << Phi << "," << Vcm[0] << "," << Vcm[1] << "," << Wcm;
		forcesFileCsv.close();
	}
}//GeneratePositionString()