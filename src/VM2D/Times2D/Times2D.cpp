/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Times2D.cpp                                                      |
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
\brief Файл кода с описанием класса Times
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7   
\date 22 ноября 2019 г.
*/

#include "Times2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

//Генерация заголовка файла временной статистики
void Times::GenerateStatHeader() const
{
	std::stringstream timeStatFileName;
	timeStatFileName << W.getPassport().dir << "timestat";

	std::ofstream timeStatFile(timeStatFileName.str());
	VMlib::PrintLogoToTextFile(timeStatFile, timeStatFileName.str(), "Time statistics");

	VMlib::PrintHeaderToTextFile(timeStatFile, "step Time N tStep tMem tMatRhs tSolve tConvVelo tDiffVelo tForce tVelPres tMove tInside tRestr tWakeSort tSave tOther");

	timeStatFile.close();
	timeStatFile.clear();
}//GenerateStatHeader()



//Сохранение строки со статистикой в файл временной статистики
void Times::GenerateStatString() const
{
	std::ofstream timestatFile(W.getPassport().dir + "timestat", std::ios::app);

	timestatFile << std::endl
		<< W.getCurrentStep() << "\t"
		<< W.getPassport().physicalProperties.getCurrTime() << "\t"
		<< W.getWake().vtx.size() << "\t"
		<< dT(timeWholeStep) << "\t"
		<< dT(timeReserveMemoryForMatrixAndRhs) << "\t"
		<< dT(timeFillMatrixAndRhs) << "\t"
		<< dT(timeSolveLinearSystem) << "\t"
		<< dT(timeCalcVortexConvVelo) << "\t"
		<< dT(timeCalcVortexDiffVelo) << "\t"
		<< dT(timeGetHydroDynamForce) << "\t"
		<< dT(timeVP) << "\t"
		<< dT(timeMoveVortexes) << "\t"
		<< dT(timeCheckInside) << "\t"
		<< dT(timeRestruct) << "\t"
		<< dT(timeWakeSort) << "\t"
		<< dT(timeSaveKadr) << "\t"
		<< dT(timeOther);

	timestatFile.close();
}//GenerateStatString()


//Обнуление состояния временной статистики
void Times::ToZero()
{
	ToZeroPeriod(timeWholeStep);
	ToZeroPeriod(timeReserveMemoryForMatrixAndRhs);
	ToZeroPeriod(timeFillMatrixAndRhs);
	ToZeroPeriod(timeSolveLinearSystem);
	ToZeroPeriod(timeCalcVortexConvVelo);
	ToZeroPeriod(timeCalcVortexDiffVelo);
	ToZeroPeriod(timeGetHydroDynamForce);
	ToZeroPeriod(timeVP);
	ToZeroPeriod(timeMoveVortexes);
	ToZeroPeriod(timeCheckInside);
	ToZeroPeriod(timeRestruct);
	ToZeroPeriod(timeWakeSort);
	ToZeroPeriod(timeSaveKadr);
	ToZeroPeriod(timeOther);	
}// ToZero()