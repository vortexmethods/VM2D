/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.8    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/03/09     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
\version 1.8
\date 09 марта 2020 г.
*/

#include "Velocity2DBarnesHut.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
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
void VelocityBarnesHut::CreatePointsCopy(PointsCopy& pointsCopy, const WakeDataBase& points, const PointType type)
{
	pointsCopy.reserve(pointsCopy.size() + points.vtx.size());

	for (int i = 0; i < (int)points.vtx.size(); ++i) 
		pointsCopy.emplace_back(points.vtx[i], {-1, -1}, i, type);
}//CreatePointsCopy(...)


// Создание массива указателей на массив точек (для сортировки), используется для массива sheetsCopy
void VelocityBarnesHut::CreateSheetsCopy(PointsCopy& pointsCopy, const PointType type)
{
	size_t nSh = 0;
	for (size_t i = 0; i < W.getNumberOfBoundary(); ++i)
		nSh += W.getBoundary(i).sheets.getSheetSize();
	
	pointsCopy.reserve(pointsCopy.size() + nSh);
	int num;
	switch (type)
	{
	case PointType::sheetGam:
		num = 0;
		for (int i = 0; i < (int)(W.getNumberOfBoundary()); ++i)
			for (int j = 0; j < (int)(W.getBoundary(i).sheets.getSheetSize()); ++j)
			{
				pointsCopy.emplace_back(Vortex2D{ 0.5 * (W.getAirfoil(i).getR(j) + W.getAirfoil(i).getR(j + 1)), W.getBoundary(i).sheets.attachedVortexSheet(j, 0) * W.getAirfoil(i).len[j] }, { i, j }, num, type);
				num++;
			}
		break;
	case PointType::source:
		num = 0;
		for (int i = 0; i < (int)(W.getNumberOfBoundary()); ++i)
			for (int j = 0; j < (int)(W.getBoundary(i).sheets.getSheetSize()); ++j)
			{
				pointsCopy.emplace_back(Vortex2D{ 0.5 * (W.getAirfoil(i).getR(j) + W.getAirfoil(i).getR(j + 1)), W.getBoundary(i).sheets.attachedSourceSheet(j, 0) * W.getAirfoil(i).len[j] }, { i, j }, num, type);
				num++;
			}
		break;
	}
}//CreateSheetsCopy(...)

//Вычисление конвективных скоростей вихрей и виртуальных вихрей в вихревом следе, а также в точках wakeVP
void VelocityBarnesHut::CalcConvVelo()
{

	/// \todo Реализовать засечки времени и разобраться с вычисление скоростей в самих виртуальных вихрях (надо вызывать функцию CalcConvVelocityAtVirtualVortexes)
	
	W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
	
	int bou, pan;
	for (size_t i = 0; i < sheetsGamCopy.size(); ++i)
	{
		bou = sheetsGamCopy.aflPnl[i].first;
		pan = sheetsGamCopy.aflPnl[i].second;
		sheetsGamCopy.vtx[i].g() += W.getBoundary(bou).sheets.freeVortexSheet(pan, 0) * W.getAirfoil(bou).len[pan];
	}

	if (treeSheetsGam)
		treeSheetsGam->rootCell.CalculateCellsParams();

	double timeCoeff, timeBiot, timeSum;

	timeCoeff = timeBiot = timeSum = 0.0;

	// WAKE:
	if (treeWake)
	{
#pragma omp parallel for reduction(+:timeCoeff, timeBiot, timeSum) schedule(dynamic, 16)
		for (int i = 0; i < (int)(treeWake->lowCells.size()); ++i)
		{
			std::vector<numvector<double, 3>> savedEe2;
			auto& lowCellI = *(treeWake->lowCells[i]);
			
			lowCellI.ClearABCD();

			// расчет влияния на wake от wake 
			lowCellI.CalcABCDandCloseCellsToLowLevel(treeWake->rootCell, true);
			lowCellI.CalcConvVeloByBiotSavartFromVortices(true, savedEe2);

			lowCellI.closeCells.clear();

			// расчет влияния на wake от sheets 
			if (treeSheetsGam)
			{
				lowCellI.CalcABCDandCloseCellsToLowLevel(treeSheetsGam->rootCell, true);
				lowCellI.CalcConvVeloByBiotSavartFromSheets(true, savedEe2);
			}
			
			lowCellI.SetDomainRadius(savedEe2);

			lowCellI.closeCells.clear();
			
			// расчет влияния на wake от sources 
			if (treeSources)
			{
				lowCellI.CalcABCDandCloseCellsToLowLevel(treeSources->rootCell, true);
				lowCellI.CalcConvVeloByBiotSavartFromSources();
			}
			lowCellI.closeCells.clear();
			
			for (auto it = lowCellI.itBegin; it != lowCellI.itEnd; ++it)
				lowCellI.CalcInfluenceFromVortexFarCells(it);

		}
		//treeWake->rootCell->PrintTree();

		//std::cout << "CalcCoeffToLowLevel: " << timeCoeff << std::endl;
		//std::cout << "CalcConvVeloByBiotSavart: " << timeBiot << std::endl;
		//std::cout << "CalcConvVeloToLowLevelCell: " << timeSum << std::endl;
	}

	// расчет скоростей самих виртуальных вихрей и их влияния на wake
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
		W.getBoundary(bou).CalcConvVelocityAtVirtualVortexes(virtualVortexesParams[bou].convVelo);
	
	W.getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();

	// вычисление eps* для виртуальных вихрей 
	W.getTimestat().timeCalcVortexDiffVelo.first += omp_get_wtime();
	if(treeSheetsGam)
	for (size_t i = 0; i < treeSheetsGam->lowCells.size(); ++i)
	{		
		//if(treeWake)
			//treeSheetsGam->lowCells[i]->CalcABCDandCloseCellsToLowLevel(treeWake->rootCell, false); // закомментировать!!! 23.01
		treeSheetsGam->lowCells[i]->CalcDomainRadiusForVirtualWake();
	}
	W.getTimestat().timeCalcVortexDiffVelo.second += omp_get_wtime();

	W.getTimestat().timeVP.first += omp_get_wtime();


	// WAKEVP:
	// если нужно, расчет влияния на pointsCopyVP
	//if ((W.getPassport().timeDiscretizationProperties.saveVP > 0) && (!(W.getCurrentStep() % W.getPassport().timeDiscretizationProperties.saveVP)))
	//{
	//	if (treeVP)
	//	{
	//		std::vector<numvector<double, 3>> savedEe2;

	//		for (size_t i = 0; i < treeVP->lowCells.size(); ++i)
	//		{
	//			treeVP->lowCells[i]->ClearABCD();

	//			// расчет влияния на wakeVP от wake
	//			if (treeWake)
	//			{
	//				treeVP->lowCells[i]->CalcABCDandCloseCellsToLowLevel(treeWake->rootCell, true);
	//				treeVP->lowCells[i]->CalcConvVeloByBiotSavartFromVortices(false, savedEe2);
	//				treeVP->lowCells[i]->closeCells.clear();
	//			}

	//			// расчет влияния на wakeVP от sources 
	//			if (treeSources)
	//			{
	//				treeVP->lowCells[i]->CalcABCDandCloseCellsToLowLevel(treeSources->rootCell, true);
	//				treeVP->lowCells[i]->CalcConvVeloByBiotSavartFromSources();
	//				treeVP->lowCells[i]->closeCells.clear();
	//			}
	//			// расчет влияния на wakeVP от sheets 
	//			if (treeSheetsGam)
	//			{
	//				treeVP->lowCells[i]->CalcABCDandCloseCellsToLowLevel(treeSheetsGam->rootCell, true);
	//				treeVP->lowCells[i]->CalcConvVeloByBiotSavartFromSheets(false, savedEe2);
	//				treeVP->lowCells[i]->closeCells.clear();
	//			}

	//			for (auto it = treeVP->lowCells[i]->itBegin; it != treeVP->lowCells[i]->itEnd; ++it)
	//				treeVP->lowCells[i]->CalcInfluenceFromVortexFarCells(it);
	//		}
	//	}//if (treeVP)
	//}
	
	/// \todo Нужно раздать скорости в точках VP в оригинальные значения 
	W.getTimestat().timeVP.second += omp_get_wtime();

	W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
	
	//раздача копий в оригинальные значения 
	wakeVortexesParams.convVelo.swap(pointsCopyWake.velo);
	wakeVortexesParams.epsastWake.swap(pointsCopyWake.domainRadius);

	W.getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();
	
}//CalcConvVelo()

void VelocityBarnesHut::BuildTrees(PointType type)
{
	switch (type)
	{
		case PointType::wake: 
		{
			pointsCopyWake.clear();
			CreatePointsCopy(pointsCopyWake, W.getWake(), PointType::wake);
			if (pointsCopyWake.size() > 0)
				treeWake.reset(new Tree(W, pointsCopyWake));
			else 
				treeWake = nullptr;
			break;
		}
		case PointType::sheetGam:
		{
			sheetsGamCopy.clear();
			CreateSheetsCopy(sheetsGamCopy, PointType::sheetGam);
			if (sheetsGamCopy.size() > 0)
				treeSheetsGam.reset(new Tree(W, sheetsGamCopy));
			else
				treeSheetsGam = nullptr;
			break;
		}
		case PointType::wakeVP:
		{
			pointsCopyVP.clear();
			CreatePointsCopy(pointsCopyVP, W.getMeasureVP().getWakeVP(), PointType::wakeVP);
			if (pointsCopyVP.size() > 0)
				treeVP.reset(new Tree(W, pointsCopyVP));
			else
				treeVP = nullptr;
			break;
		}
		case PointType::sourceWake :
		{
			sourcesCopy.clear();
			CreateSheetsCopy(sourcesCopy, PointType::source);
			CreatePointsCopy(sourcesCopy, W.getSource(), PointType::sourceWake);
			if (sourcesCopy.size() > 0)
				treeSources.reset(new Tree(W, sourcesCopy));
			else 
				treeSources = nullptr;
			break;
		}
	default:
		break;
	}
}//BuildTrees(...)


void VelocityBarnesHut::Initialization()
{
	W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();

	BuildTrees(PointType::wake);
	if (treeWake)
		treeWake->rootCell.CalculateCellsParams();

	BuildTrees(PointType::sourceWake);
	if (treeSources)
		treeSources->rootCell.CalculateCellsParams();

	if (W.getNumberOfAirfoil() > 0)
		BuildTrees(PointType::sheetGam);

	W.getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();

	W.getTimestat().timeVP.first += omp_get_wtime();
	if ((W.getPassport().timeDiscretizationProperties.saveVP > 0) && (!(W.getCurrentStep() % W.getPassport().timeDiscretizationProperties.saveVP)))
		BuildTrees(PointType::wakeVP);
	W.getTimestat().timeVP.second += omp_get_wtime();
}//Initialization()

void VelocityBarnesHut::GetWakeInfluenceToRhs(std::vector<double>& wakeRhs) const
{
	//Здесь предполагается, что тип схемы у всех профилей одинаковый
	size_t shDim = W.getBoundary(0).sheetDim;
	std::vector<double> velI(shDim, 0.0);

	Point2D approxVel;
	for (size_t i = 0; i < treeSheetsGam->lowCells.size(); ++i)
	{

		if (treeSources)
		{
			treeSheetsGam->lowCells[i]->ClearABCD();
			treeSheetsGam->lowCells[i]->closeCells.clear();
			treeSheetsGam->lowCells[i]->CalcABCDandCloseCellsToLowLevel(treeSources->rootCell, true);

			for (auto it = treeSheetsGam->lowCells[i]->itBegin; it != treeSheetsGam->lowCells[i]->itEnd; ++it) // цикл по панелям
			{
				int aflN = it.getAflPnl().first;
				int panNum = it.getAflPnl().second;

				velI.assign(shDim, 0.0);

				for (size_t j = 0; j < treeSheetsGam->lowCells[i]->closeCells.size(); ++j)
					W.getAirfoil(aflN).GetInfluenceFromSourcesToPanel(panNum, &(treeSheetsGam->lowCells[i]->closeCells[j]->itBegin.getVtx()), \
						std::distance(treeSheetsGam->lowCells[i]->closeCells[j]->itBegin, treeSheetsGam->lowCells[i]->closeCells[j]->itEnd), velI);
				
				for (size_t j = 0; j < shDim; ++j)
					velI[j] *= IDPI / W.getAirfoil(aflN).len[panNum];

				// результат приближенного влияния записывается в veloCopy конкретной панели (там уже есть накопленные влияния от дальних источников)
				treeSheetsGam->lowCells[i]->CalcInfluenceFromSourceFarCells(it);

				velI[0] += treeSheetsGam->allPnt.velo[it.getNum()] & W.getAirfoil(aflN).tau[panNum];

				for (size_t k = 0; k < shDim; ++k)
					wakeRhs[W.getDispBoundaryInSystem(aflN) + panNum + k * W.getAirfoil(aflN).getNumberOfPanels()] = velI[k];
			}
		}//if(treeSources)

		// аналогично для wake		
		if (treeWake)
		{
			treeSheetsGam->lowCells[i]->ClearABCD();
			treeSheetsGam->lowCells[i]->closeCells.clear();
			treeSheetsGam->lowCells[i]->CalcABCDandCloseCellsToLowLevel(treeWake->rootCell, true);
			int count = 0; 
			for (auto it = treeSheetsGam->lowCells[i]->itBegin; it != treeSheetsGam->lowCells[i]->itEnd; ++it) // цикл по панелям
			{
				//if (count == 0)
					//std::cout << treeSheetsGam->lowCells[i]->itEnd - treeSheetsGam->lowCells[i]->itBegin << std::endl;
				int aflN = it.getAflPnl().first;
				int panNum = it.getAflPnl().second;

				velI.assign(shDim, 0.0);

				for (size_t j = 0; j < treeSheetsGam->lowCells[i]->closeCells.size(); ++j)
					W.getAirfoil(aflN).GetInfluenceFromVorticesToPanel(panNum, &(treeSheetsGam->lowCells[i]->closeCells[j]->itBegin.getVtx()), \
						std::distance(treeSheetsGam->lowCells[i]->closeCells[j]->itBegin, treeSheetsGam->lowCells[i]->closeCells[j]->itEnd), velI);
				
				for (size_t j = 0; j < shDim; ++j)
					velI[j] *= IDPI / W.getAirfoil(aflN).len[panNum];

				// результат приближенного влияния записывается в velo конкретной панели (там уже есть накопленные влияния от дальних источников)
				treeSheetsGam->lowCells[i]->CalcInfluenceFromVortexFarCells(it);

				velI[0] += treeSheetsGam->allPnt.velo[it.getNum()] & W.getAirfoil(aflN).tau[panNum];
				
				for (size_t k = 0; k < shDim; ++k)
					wakeRhs[W.getDispBoundaryInSystem(aflN) + panNum + k * W.getAirfoil(aflN).getNumberOfPanels()] += velI[k];
				
				//std::cout << count << std::endl;
				//count++;
			}
		}//if (treeWake)
	}

}//GetWakeInfluenceToRhs(...)

#if defined(USE_CUDA)
//Генерация вектора влияния вихревого следа на профиль
void VelocityBarnesHut::GPUGetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeVelo) const
{
	const int& id = W.getParallel().myidWork;

	const size_t& nvt = W.getWake().vtx.size();
	const size_t& nsr = W.getSource().vtx.size();

	if (afl.numberInPassport == 0)
	{
		size_t nTotPan = 0;
		for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s)
			nTotPan += W.getAirfoil(s).getNumberOfPanels();

		double*& dev_ptr_pt = afl.devRPtr;
		double*& dev_ptr_vt = W.getWake().devVtxPtr;
		double*& dev_ptr_sr = W.getSource().devVtxPtr;
		double*& dev_ptr_rhs = afl.devRhsPtr;
		std::vector<double> locrhs(nTotPan);

		VMlib::parProp par = W.getParallel().SplitMPI(nTotPan, true);

		if ((nvt > 0) || (nsr > 0))
		{
			cuCalculateRhs(par.myDisp, par.myLen, nTotPan, dev_ptr_pt, nvt, dev_ptr_vt, nsr, dev_ptr_sr, dev_ptr_rhs);

			W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_rhs, (double*)&locrhs[0]);

			std::vector<double> newRhs;
			if (id == 0)
				newRhs.resize(nTotPan);

			MPI_Gatherv(locrhs.data(), par.myLen, MPI_DOUBLE, newRhs.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

			if (id == 0)
			{
				size_t curGlobPnl = 0;
				for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s)
				{
					std::vector<double>& tmpRhs = W.getNonConstAirfoil(s).tmpRhs;
					const size_t& np = W.getAirfoil(s).getNumberOfPanels();
					tmpRhs.resize(0);
					tmpRhs.insert(tmpRhs.end(), newRhs.begin() + curGlobPnl, newRhs.begin() + curGlobPnl + np);
					curGlobPnl += np;
				}
			}
		}
	}

	if (id == 0)
	{
		if ((nvt > 0) || (nsr > 0))
			wakeVelo = std::move(afl.tmpRhs);
		else
			wakeVelo.resize(afl.getNumberOfPanels(), 0.0);
	}

}//GPUGetWakeInfluenceToRhs(...)
#endif

///Заполнение правой части СЛАУ 
void VelocityBarnesHut::FillRhs(Eigen::VectorXd& rhs) const
{
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_RHS))
	Eigen::VectorXd locRhs;
	std::vector<double> lastRhs(W.getNumberOfBoundary());

	size_t currentRow = 0;

	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		//double tt0 = omp_get_wtime();

		const Airfoil& afl = W.getAirfoil(bou);
		size_t np = afl.getNumberOfPanels();

		size_t nVars;

		if (W.getParallel().myidWork == 0)
		{
			nVars = W.getBoundary(bou).GetUnknownsSize();

			locRhs.resize(nVars);
		}


		std::vector<double> wakeRhs;

		GPUGetWakeInfluenceToRhs(afl, wakeRhs);

		std::vector<double> vInfRhs;
		afl.GetInfluenceFromVInfToPanel(vInfRhs);

		if (W.getParallel().myidWork == 0)
		{
#pragma omp parallel for \
	default(none) \
	shared(locRhs, afl, bou, wakeRhs, vInfRhs, np) \
	schedule(dynamic, DYN_SCHEDULE)
			for (int i = 0; i < (int)afl.getNumberOfPanels(); ++i)
			{
				/// \todo ВАЖНО!!!
				//Здесь учет скорости профиля и присоединенных слоев источников сделан только для прямолинейных панелей
				locRhs(i) = -vInfRhs[i] - wakeRhs[i] + 0.25 * ((afl.getV(i) + afl.getV(i + 1)) & afl.tau[i]); //0.25 * (afl.getV(i) + afl.getV(i + 1))*afl.tau[i] - прямолинейные
				if (W.getBoundary(bou).sheetDim > 1)
					locRhs(np + i) = -vInfRhs[np + i] - wakeRhs[np + i];

				//влияние присоединенных слоев от самого себя и от других профилей				
				for (size_t q = 0; q < W.getNumberOfBoundary(); ++q)
				{
					const auto& sht = W.getBoundary(q).sheets;
					const auto& iq = W.getIQ(bou, q);

					const Airfoil& aflOther = W.getAirfoil(q);
					if (W.getMechanics(q).isDeform || W.getMechanics(q).isMoves)
					{
						for (size_t j = 0; j < aflOther.getNumberOfPanels(); ++j)
						{
							if ((i != j) || (bou != q))
							{
								locRhs(i) += -iq.first(i, j) * sht.attachedVortexSheet(j, 0);
								locRhs(i) += -iq.second(i, j) * sht.attachedSourceSheet(j, 0); // getIQ(bou, q).second(i, j) пока забито нулем для криволинейных
							}//if (i != j)
						}//for j
					}
				}//for q
			}//for i

			lastRhs[bou] = 0.0;

			for (size_t q = 0; q < afl.gammaThrough.size(); ++q)
				lastRhs[bou] += afl.gammaThrough[q];

			//размазываем правую часть		
			for (size_t i = 0; i < nVars; ++i)
				rhs(i + currentRow) = locRhs(i);

			rhs(currentRow + nVars) = lastRhs[bou];

			currentRow += nVars + 1;

		}// if (parallel.myidWork == 0)

	}// for bou
#else

	std::vector<double> wakeRhs;
	wakeRhs.resize(rhs.size());


//#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_RHS))
	//GPUGetWakeInfluenceToRhs(wakeRhs); //\todo это задел 
//#else
	GetWakeInfluenceToRhs(wakeRhs);
//#endif
	
	size_t currentRow = 0;

	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		size_t nPanLast = 0;

		const Airfoil& afl = W.getAirfoil(bou);

		size_t nVars = W.getBoundary(bou).GetUnknownsSize();

		if (W.getParallel().myidWork == 0)
		{
			std::vector<double> vInfRhs;
			afl.GetInfluenceFromVInfToPanel(vInfRhs);

			//запись компонент, отвечающих за константные проекционные функции
#pragma omp parallel for default(none) shared(afl, bou, wakeRhs, vInfRhs, currentRow, rhs) schedule(static, 1)  //schedule(dynamic, DYN_SCHEDULE)
			for (int i = 0; i < afl.getNumberOfPanels(); ++i)
			{
				
				/// \todo ВАЖНО!!!
				//Здесь учет скорости профиля и присоединенных слоев источников сделан только для прямолинейных панелей
				rhs(currentRow + i) = -vInfRhs[i] - wakeRhs[currentRow + i] + 0.25 * ((afl.getV(i) + afl.getV(i + 1)) & afl.tau[i]);

				//влияние присоединенных слоев от самого себя и от других профилей				
				for (int q = 0; q < W.getNumberOfBoundary(); ++q)
				{
					const auto& sht = W.getBoundary(q).sheets;
					const auto& iq = W.getIQ(bou, q);

					const Airfoil& aflOther = W.getAirfoil(q);
					if (W.getMechanics(q).isDeform || W.getMechanics(q).isMoves)
					{
						for (size_t j = 0; j < aflOther.getNumberOfPanels(); ++j)
						{
							if ((i != j) || (bou != q))
							{

								rhs(currentRow + i) += -iq.first(i, j) * sht.attachedVortexSheet(j, 0);
								rhs(currentRow + i) += -iq.second(i, j) * sht.attachedSourceSheet(j, 0);
							}//if (i != j)
						}//for j
					}
				}//for q
			}//for i

			currentRow += afl.getNumberOfPanels();

			//запись компонент, отвечающих за линейные проекционные функции
			if (W.getBoundary(bou).sheetDim > 1)
			{
#pragma omp parallel for \
	default(none) \
	shared(afl, bou, wakeRhs, currentRow) \
	schedule(dynamic, DYN_SCHEDULE)
				for (int i = 0; i < afl.getNumberOfPanels(); ++i)
				{
					size_t iLin = afl.getNumberOfPanels() + i;

					rhs(currentRow + i) = -wakeRhs[currentRow + i];

					//влияние присоединенных слоев от самого себя и от других профилей				
					for (int q = 0; q < W.getNumberOfBoundary(); ++q)
					{
						const Airfoil& aflOther = W.getAirfoil(q);
						if (W.getMechanics(q).isDeform || W.getMechanics(q).isMoves)
						{
							for (int j = 0; j < aflOther.getNumberOfPanels(); ++j)
							{
								size_t jLin = aflOther.getNumberOfPanels() + j;
								if ((i != j) || (bou != q))
								{
									rhs(currentRow + i) += -W.getIQ(bou, q).first(iLin, jLin) * W.getBoundary(q).sheets.attachedVortexSheet(j, 0);
									rhs(currentRow + i) += -W.getIQ(bou, q).second(iLin, jLin) * W.getBoundary(q).sheets.attachedSourceSheet(j, 0);
								}//if (i != j)
							}//for j
						}
					}//for q
				}//for i

				currentRow += afl.getNumberOfPanels();
			}//if (W.getBoundary(bou).sheetDim > 1)
			

			//запись компонент, отвечающих за квадратичные проекционные функции
			if (W.getBoundary(bou).sheetDim > 2)
			{
#pragma omp parallel for \
	default(none) \
	shared(afl, bou, wakeRhs) \
	schedule(dynamic, DYN_SCHEDULE)
				for (int i = 0; i < afl.getNumberOfPanels(); ++i)
				{
					size_t iQuad = 2 * afl.getNumberOfPanels() + i;

					rhs(currentRow + i) = -wakeRhs[currentRow + i];

					//влияние присоединенных слоев от самого себя и от других профилей				
					for (int q = 0; q < W.getNumberOfBoundary(); ++q)
					{
						const Airfoil& aflOther = W.getAirfoil(q);
						if (W.getMechanics(q).isDeform || W.getMechanics(q).isMoves)
						{
							for (size_t j = 0; j < aflOther.getNumberOfPanels(); ++j)
							{
								size_t jQuad = 2 * aflOther.getNumberOfPanels() + j;
								if ((i != j) || (bou != q))
								{
									rhs(currentRow + i) += -W.getIQ(bou, q).first(iQuad, jQuad) * W.getBoundary(q).sheets.attachedVortexSheet(j, 0);
									rhs(currentRow + i) += -W.getIQ(bou, q).second(iQuad, jQuad) * W.getBoundary(q).sheets.attachedSourceSheet(j, 0);
								}//if (i != j)
							}//for j
						}
					}//for q
				}//for i
				currentRow += afl.getNumberOfPanels();
			}//if (W.getBoundary(bou).sheetDim > 2)

			rhs[currentRow] = 0.0;

			for (size_t q = 0; q < afl.gammaThrough.size(); ++q)
				rhs[currentRow] += afl.gammaThrough[q];

			currentRow++;
		}// if (parallel.myidWork == 0)

	}// for bou
#endif
}//FillRhs(...)