/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Velocity2DBarnesHut.cpp                                          |
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
\brief Файл кода с описанием класса VelocityBarnesHut
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7
\date 22 ноября 2019 г.
*/

#include "Velocity2DBarnesHut.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
//#include "Tree2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

// Конструктор
VelocityBarnesHut::VelocityBarnesHut(const World2D& W_) :
	Velocity(W_)
{
	
};

// Деструктор
VelocityBarnesHut::~VelocityBarnesHut()
{
};

// Создание массива указателей на массив точек (для сортировки), используется для массивов pointsCopyVP и pointsCopyWake
void VelocityBarnesHut::CreatePointsCopy(std::vector<PointsCopy>& pointsCopy, const WakeDataBase& points, const PointType type)
{
	pointsCopy.reserve(pointsCopy.size() + points.vtx.size());

	for (int i = 0; i < (int)points.vtx.size(); ++i) 
		pointsCopy.emplace_back(PointsCopy(points.vtx[i], -1, i, type));
}//CreatePointsCopy(...)


// Создание массива указателей на массив точек (для сортировки), используется для массива sheetsCopy
void VelocityBarnesHut::CreateSheetsCopy(std::vector<PointsCopy>& pointsCopy, const PointType type)
{
	size_t nSh = 0;
	for (size_t i = 0; i < W.getNumberOfBoundary(); ++i)
		nSh += W.getBoundary(i).sheets.getSheetSize();
	
	pointsCopy.reserve(pointsCopy.size() + nSh);

	switch (type)
	{
	case PointType::sheetGam:
		for (int i = 0; i < (int)(W.getNumberOfBoundary()); ++i)
			for (int j = 0; j < (int)(W.getBoundary(i).sheets.getSheetSize()); ++j)
				pointsCopy.emplace_back(PointsCopy(Vortex2D{0.5 * (W.getAirfoil(i).getR(j) + W.getAirfoil(i).getR(j + 1)), 0.0}, i, j, type ));
		break;
	case PointType::source:
		for (int i = 0; i < (int)(W.getNumberOfBoundary()); ++i)
			for (int j = 0; j < (int)(W.getBoundary(i).sheets.getSheetSize()); ++j)
				pointsCopy.emplace_back(PointsCopy( Vortex2D{0.5 * (W.getAirfoil(i).getR(j) + W.getAirfoil(i).getR(j + 1)), W.getBoundary(i).sheets.attachedSourceSheet(j, 0) * W.getAirfoil(i).len[j]}, i, j, type ));
		break;
	}
}//CreateSheetsCopy(...)

//Вычисление конвективных скоростей вихрей и виртуальных вихрей в вихревом следе, а также в точках wakeVP
void VelocityBarnesHut::CalcConvVelo(timePeriod& convWakeTime, timePeriod& convVPTime)
{

	/// \todo Реализовать засечки времени

	//CreatePointsCopy(pointsCopyWake, W.getWake(), PointType::wake);
	//CreatePointsCopy(pointsCopyWake, W.getSource(), PointType::sourceWake);

	
	double tCPUSTART, tCPUEND;

	tCPUSTART = omp_get_wtime();
	treeWake->BuildTree(pointsCopyWake);

	tCPUEND = omp_get_wtime();
	std::cout << "BuildTree: " << tCPUEND - tCPUSTART << std::endl;

	tCPUSTART = omp_get_wtime();
	treeWake->rootCell->CalculateCellsParams();
	tCPUEND = omp_get_wtime();
	std::cout << "CalculateCellsParams: " << tCPUEND - tCPUSTART << std::endl;


	double timeCoeff, timeBiot, timeSum, tCPUA, tCPUB;

	timeCoeff = timeBiot = timeSum = 0.0;

	// расчет влияния на pointsCopyAll
	for (size_t i = 0; i < treeWake->lowCells.size(); ++i)
	{
		tCPUA = omp_get_wtime();
		treeWake->lowCells[i]->CalcCoeffToLowLevel(treeWake->rootCell);
		tCPUB = omp_get_wtime();

		timeCoeff += tCPUB - tCPUA;

		treeWake->lowCells[i]->CalcConvVeloByBiotSavart();
		tCPUA = omp_get_wtime();

		timeBiot += tCPUA - tCPUB;
		
		//treeWake->lowCells[i]->CalcConvVeloByBiotSavart(timeBiot);
		//tCPUA = omp_get_wtime();
		
		// \todo ????? Переделать
		//treeWake->lowCells[i]->CalcConvVeloToLowLevelCell();

		tCPUB = omp_get_wtime();

		timeSum += tCPUB - tCPUA;
	}
	
	std::cout << "CalcCoeffToLowLevel: " << timeCoeff << std::endl;
	std::cout << "CalcConvVeloByBiotSavart: " << timeBiot << std::endl;
	std::cout << "CalcConvVeloToLowLevelCell: " << timeSum << std::endl;


	treeSheetsGam->BuildTree(sheetsGamCopy);


	// расчет скоростей самих виртуальных вихрей и их влияния на wake
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		std::vector<Point2D>& Vel = W.getNonConstVelocity().virtualVortexesParams[bou].convVelo;
		for (int i = 0; i < Vel.size(); ++i)
			Vel[i] = W.getBoundary(bou).afl.getV(W.getBoundary(bou).virtualWake.aflPan[i].second) \
			+ W.getBoundary(bou).virtualWake.vecHalfGamma[i] \
			- W.getPassport().physicalProperties.V0(); // V0 потом прибавляется ко всем скоростям в функции MoveVortexes

		W.getBoundary(bou).GetConvVelocityToSetOfPointsFromVirtualVortexes(W.getWake(), wakeVortexesParams.convVelo);
	}


	// если нужно, расчет влияния на pointsCopyVP
	if ((W.getPassport().timeDiscretizationProperties.saveVP > 0) && (!(W.getCurrentStep() % W.getPassport().timeDiscretizationProperties.saveVP)))
	{

		for (size_t i = 0; i < treeVP->lowCells.size(); ++i)
		{
			treeVP->lowCells[i]->CalcCoeffToLowLevel(treeWake->rootCell);
			treeVP->lowCells[i]->CalcConvVeloByBiotSavart();
			//treeVP->lowCells[i]->CalcConvVeloToLowLevelCell();
		}
		// расчет влияния виртуальных вихрей 
		for (size_t i = 0; i < W.getNumberOfBoundary(); i++)
			W.getNonConstBoundary(i).GetConvVelocityToSetOfPointsFromVirtualVortexes(W.getMeasureVP().getWakeVP(), W.getNonConstMeasureVP().getNonConstVelocity());
	}
}//CalcConvVelo()

void VelocityBarnesHut::BuildTrees(PointType type)
{
	switch (type)
	{
		case PointType::wake: 
		{
			pointsCopyWake.clear();
			CreatePointsCopy(pointsCopyWake, W.getWake(), PointType::wake);
			treeWake.reset(new Tree(W));
			treeWake->BuildTree(pointsCopyWake);
			break;
		}
		case PointType::sheetGam:
		{
			sheetsGamCopy.clear();
			CreateSheetsCopy(sheetsGamCopy, PointType::sheetGam);
			treeSheetsGam.reset(new Tree(W));
			treeSheetsGam->BuildTree(sheetsGamCopy);
			break;
		}
		case PointType::wakeVP:
		{
			pointsCopyVP.clear();
			CreatePointsCopy(pointsCopyVP, W.getMeasureVP().getWakeVP(), PointType::wakeVP);
			treeVP.reset(new Tree(W));
			treeVP->BuildTree(pointsCopyVP);
			break;
		}
		case PointType::sourceWake :
		{
			sourcesCopy.clear();
			CreateSheetsCopy(sourcesCopy, PointType::source);
			CreatePointsCopy(sourcesCopy, W.getSource(), PointType::sourceWake);
			treeSources->BuildTree(sourcesCopy);
			break;
		}
	default:
		break;
	}
}//BuildTrees(...)


void VelocityBarnesHut::GetWakeInfluenceToRhs(std::vector<double>& wakeRhs) const
{
	//Здесь предполагается, что тип схемы у всех профилей одинаковы
	size_t shDim = W.getBoundary(0).sheetDim;
	std::vector<double> velI(shDim, 0.0);
	std::vector<double> sumVelI(shDim, 0.0);

	Point2D approxVel;

	for (size_t i = 0; i < treeSheetsGam->lowCells.size(); ++i)
	{
		treeSheetsGam->lowCells[i]->closeCells.clear();
		treeSheetsGam->lowCells[i]->CalcCoeffToLowLevel(treeWake->rootCell);

		for (auto it = treeSheetsGam->lowCells[i]->itBegin; it != treeSheetsGam->lowCells[i]->itEnd; ++it) // цикл по панелям
		{
			int aflN = it->aflN;
			int panNum = it->num;
			for (size_t k = 0; k < shDim; ++k)
				sumVelI[k] = 0;

			for (size_t j = 0; j < treeSheetsGam->lowCells[i]->closeCells.size(); ++j)
			{
				W.getAirfoil(aflN).GetInfluenceFromVorticesToPanel(panNum, &(*treeSheetsGam->lowCells[i]->closeCells[j]->itBegin), std::distance(treeSheetsGam->lowCells[i]->closeCells[j]->itBegin, treeSheetsGam->lowCells[i]->closeCells[j]->itEnd), velI);

				for (size_t k = 0; k < shDim; ++k)
					sumVelI[k] += velI[k];
			}

			// результат приближенного влияния завписывается в veloCopy
			treeSheetsGam->lowCells[i]->CalcInfluenceFromVortexFarCells(it);

//			sumVelI[0] += it->veloCopy & W.getAirfoil(aflN).tau[panNum];

			for (size_t k = 0; k < shDim; ++k)
				wakeRhs[W.getDispBoundaryInSystem(aflN) + panNum + k * W.getAirfoil(aflN).getNumberOfPanels()] = -sumVelI[k];
		}

		treeSheetsGam->lowCells[i]->closeCells.clear();
		treeSheetsGam->lowCells[i]->CalcCoeffToLowLevel(treeSources->rootCell);

		for (auto it = treeSheetsGam->lowCells[i]->itBegin; it != treeSheetsGam->lowCells[i]->itEnd; ++it) // цикл по панелям
		{
			int aflN = it->aflN;
			int panNum = it->num;
			for (size_t k = 0; k < shDim; ++k)
				sumVelI[k] = 0;

			for (size_t j = 0; j < treeSheetsGam->lowCells[i]->closeCells.size(); ++j)
			{
				W.getAirfoil(aflN).GetInfluenceFromSourcesToPanel(panNum, &(*treeSheetsGam->lowCells[i]->closeCells[j]->itBegin), std::distance(treeSheetsGam->lowCells[i]->closeCells[j]->itBegin, treeSheetsGam->lowCells[i]->closeCells[j]->itEnd), velI);

				for (size_t k = 0; k < shDim; ++k)
					sumVelI[k] += velI[k];
			}

			treeSheetsGam->lowCells[i]->CalcInfluenceFromSourceFarCells(it);

			sumVelI[0] += it->veloCopy & W.getAirfoil(aflN).tau[panNum];

			for (size_t k = 0; k < shDim; ++k)
				wakeRhs[W.getDispBoundaryInSystem(aflN) + panNum + k * W.getAirfoil(aflN).getNumberOfPanels()] += -sumVelI[k];
		}

	}

}//GetWakeInfluenceToRhs(...)



void VelocityBarnesHut::FillRhs(Eigen::VectorXd& rhs) const
{
	treeWake->rootCell->CalculateCellsParams();


	
}//FillRhs(...)