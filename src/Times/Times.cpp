/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.2    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/06/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Times.cpp                                                        |
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
\brief Файл кода с описанием класса Times
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.2
\date 14 июня 2018 г.
*/


#include "Times.h"
#include "World2D.h"


//Генерация заголовка файла временной статистики
void Times::GenerateStatHeader() const
{
	std::stringstream timeStatFileName;
	timeStatFileName << W.getPassport().dir << "timestat";

	std::ofstream timeStatFile(timeStatFileName.str());
	PrintLogoToTextFile(timeStatFile, timeStatFileName.str(), "Time statistics");

	PrintHeaderToTextFile(timeStatFile, "step Time N tStep tMem tMatRhs tSolve tConvVelo tDiffVelo tForce tMove tInside tRestr tWakeSort tSave");

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
		<< dT(timeMoveVortexes) << "\t"
		<< dT(timeCheckInside) << "\t"
		<< dT(timeRestruct) << "\t"
		<< dT(timeWakeSort) << "\t"
		<< dT(timeSaveKadr);

	timestatFile.close();
}//GenerateStatString()


//Обнуление состояния временной статистики
void Times::ToZero()
{
	ToZero(timeWholeStep);
	ToZero(timeReserveMemoryForMatrixAndRhs);
	ToZero(timeFillMatrixAndRhs);
	ToZero(timeSolveLinearSystem);
	ToZero(timeCalcVortexConvVelo);
	ToZero(timeCalcVortexDiffVelo);
	ToZero(timeGetHydroDynamForce);
	ToZero(timeMoveVortexes);
	ToZero(timeCheckInside);
	ToZero(timeRestruct);
	ToZero(timeWakeSort);
	ToZero(timeSaveKadr);
}// ToZero()