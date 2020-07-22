/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.9    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/07/22     |
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
\brief Файл кода с описанием класса VelocityBarnesHut
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.9
\date 22 июля 2020 г.
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


//Вычисление конвективных скоростей вихрей и виртуальных вихрей в вихревом следе, а также в точках wakeVP
void VelocityBarnesHut::CalcConvVelo()
{
	/// \todo Реализовать засечки времени и разобраться с вычисление скоростей в самих виртуальных вихрях (надо вызывать функцию CalcConvVelocityAtVirtualVortexes)
	
	W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
	
	Tree& treeWakeRef = W.getNonConstTree(PointType::wake);
	Tree& treeSheetsGamRef = W.getNonConstTree(PointType::sheetGam);
	//Tree& treeVPRef = W.getNonConstTree(PointType::wakeVP);
	Tree& treeSheetsSourceRef = W.getNonConstTree(PointType::source);
	Tree& treeSourceWakeRef = W.getNonConstTree(PointType::sourceWake);


	int bou, pan;

	if (W.treeSheetsGam)
	for (size_t i = 0; i < treeSheetsGamRef.allPnt.size(); ++i)
	{
		bou = treeSheetsGamRef.allPnt.aflPnl[i].first;
		pan = treeSheetsGamRef.allPnt.aflPnl[i].second;
		treeSheetsGamRef.allPnt.vtx[i].g() += W.getBoundary(bou).sheets.freeVortexSheet(pan, 0) * W.getAirfoil(bou).len[pan];
	}

	if (W.treeSheetsGam)
		treeSheetsGamRef.rootCell.CalculateCellsParams();

	double timeCoeff, timeBiot, timeSum;

	timeCoeff = timeBiot = timeSum = 0.0;

	// WAKE:
	if (W.treeWake)
	{
#pragma omp parallel for  reduction(+:timeCoeff, timeBiot, timeSum) schedule(dynamic, 16)
		for (int i = 0; i < (int)(treeWakeRef.lowCells.size()); ++i)
		{
			std::vector<numvector<double, 3>> savedEe2;
			auto& lowCellI = *(treeWakeRef.lowCells[i]);
			
			lowCellI.ClearABCD();

			// расчет влияния на wake от wake 
			lowCellI.CalcABCDandCloseCellsToLowLevel(treeWakeRef.rootCell, true);
			lowCellI.CalcConvVeloByBiotSavartFromVortices(true, savedEe2);

			lowCellI.closeCells.clear();

			// расчет влияния на wake от sheets 
			if (W.treeSheetsGam)
			{
				lowCellI.CalcABCDandCloseCellsToLowLevel(treeSheetsGamRef.rootCell, true);
				lowCellI.CalcConvVeloByBiotSavartFromSheets(true, savedEe2);
			}
			
			lowCellI.SetDomainRadius(savedEe2);

			lowCellI.closeCells.clear();
			
			for (auto it = lowCellI.itBegin; it != lowCellI.itEnd; ++it)
				lowCellI.CalcInfluenceFromVortexFarCells(it);

			lowCellI.ClearABCD();

			// расчет влияния на wake от sources 
			if (W.treeSheetsSource)
			{
				lowCellI.CalcABCDandCloseCellsToLowLevel(treeSheetsSourceRef.rootCell, true);
				lowCellI.CalcConvVeloByBiotSavartFromSources();
			}
			lowCellI.closeCells.clear();
			
			// расчет влияния на wake от sourcesWake 
			if (W.treeSourcesWake)
			{
				lowCellI.CalcABCDandCloseCellsToLowLevel(treeSourceWakeRef.rootCell, true);
				lowCellI.CalcConvVeloByBiotSavartFromSources();
			}
			lowCellI.closeCells.clear();

			for (auto it = lowCellI.itBegin; it != lowCellI.itEnd; ++it)
				lowCellI.CalcInfluenceFromSourceFarCells(it);
		}
	}//if(W.treeWake)

	
	// расчет скоростей самих виртуальных вихрей
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
		W.getBoundary(bou).CalcConvVelocityAtVirtualVortexes(virtualVortexesParams[bou].convVelo);
	
	W.getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();

	// вычисление eps* для виртуальных вихрей 
	W.getTimestat().timeCalcVortexDiffVelo.first += omp_get_wtime();

	if(W.treeSheetsGam)
		for (size_t i = 0; i < treeSheetsGamRef.lowCells.size(); ++i)
		{
			treeSheetsGamRef.lowCells[i]->closeCells.clear();
			if (W.treeWake)
				treeSheetsGamRef.lowCells[i]->CalcABCDandCloseCellsToLowLevel(treeWakeRef.rootCell, false);
			treeSheetsGamRef.lowCells[i]->CalcDomainRadiusForVirtualWake();
		}

#if (defined(USE_CUDA))
	{
		for (size_t afl = 0; afl < W.getNumberOfAirfoil(); ++afl)
		{
			auto& ref = W.getNonConstVelocity().virtualVortexesParams[afl].epsastWake;
			cuCopyFixedArray(W.getBoundary(afl).virtualWake.devRadPtr, ref.data(), sizeof(double) * ref.size());
		}
	}
#endif

	W.getTimestat().timeCalcVortexDiffVelo.second += omp_get_wtime();

	W.getTimestat().timeVP.first += omp_get_wtime();

	/// \todo сделать по аналогии с Wake 
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
	if (W.treeWake)
	{
		wakeVortexesParams.convVelo.swap(W.treeWake->allPnt.velo);
		wakeVortexesParams.epsastWake.swap(W.treeWake->allPnt.domainRadius);

#if (defined(USE_CUDA))
		{
			auto& ref = wakeVortexesParams.epsastWake;
			cuCopyFixedArray(W.getWake().devRadPtr, ref.data(), sizeof(double) * ref.size());
		}
#endif
	}
	W.getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();
	
}//CalcConvVelo()

void VelocityBarnesHut::GetWakeInfluenceToRhs(std::vector<double>& wakeRhs) const
{
	Tree& treeWakeRef = W.getNonConstTree(PointType::wake);
	Tree& treeSheetsGamRef = W.getNonConstTree(PointType::sheetGam);
	Tree& treeSheetsSourceRef = W.getNonConstTree(PointType::source);
	//Tree& treeSourceWakeRef = W.getNonConstTree(PointType::sourceWake);

	//Здесь предполагается, что тип схемы у всех профилей одинаковый
	size_t shDim = W.getBoundary(0).sheetDim;
	std::vector<double> velI(shDim, 0.0);

	Point2D approxVel;
	for (size_t i = 0; i < treeSheetsGamRef.lowCells.size(); ++i)
	{

		if (W.treeSheetsSource)
		{
			treeSheetsGamRef.lowCells[i]->ClearABCD();
			treeSheetsGamRef.lowCells[i]->closeCells.clear();
			treeSheetsGamRef.lowCells[i]->CalcABCDandCloseCellsToLowLevel(treeSheetsSourceRef.rootCell, true);

			for (auto it = treeSheetsGamRef.lowCells[i]->itBegin; it != treeSheetsGamRef.lowCells[i]->itEnd; ++it) // цикл по панелям
			{
				int aflN = it.getAflPnl().first;
				int panNum = it.getAflPnl().second;

				velI.assign(shDim, 0.0);

				for (size_t j = 0; j < treeSheetsGamRef.lowCells[i]->closeCells.size(); ++j)
					W.getAirfoil(aflN).GetInfluenceFromSourcesToPanel(panNum, &(treeSheetsGamRef.lowCells[i]->closeCells[j]->itBegin.getVtx()), \
						std::distance(treeSheetsGamRef.lowCells[i]->closeCells[j]->itBegin, treeSheetsGamRef.lowCells[i]->closeCells[j]->itEnd), velI);
				
				for (size_t j = 0; j < shDim; ++j)
					velI[j] *= IDPI / W.getAirfoil(aflN).len[panNum];

				// результат приближенного влияния записывается в veloCopy конкретной панели (там уже есть накопленные влияния от дальних источников)
				treeSheetsGamRef.lowCells[i]->CalcInfluenceFromSourceFarCells(it);

				velI[0] += treeSheetsGamRef.allPnt.velo[it.getNum()] & W.getAirfoil(aflN).tau[panNum];

				for (size_t k = 0; k < shDim; ++k)
					wakeRhs[W.getDispBoundaryInSystem(aflN) + panNum + k * W.getAirfoil(aflN).getNumberOfPanels()] = velI[k];
			}
		}//if(treeSources)

		treeSheetsGamRef.lowCells[i]->ClearABCD();
		treeSheetsGamRef.lowCells[i]->closeCells.clear();

		// аналогично для wake		
		if (W.treeWake)
		{			
			treeSheetsGamRef.lowCells[i]->CalcABCDandCloseCellsToLowLevel(treeWakeRef.rootCell, true);
			//int count = 0; 
			for (auto it = treeSheetsGamRef.lowCells[i]->itBegin; it != treeSheetsGamRef.lowCells[i]->itEnd; ++it) // цикл по панелям
			{
				//if (count == 0)
					//std::cout << treeSheetsGam->lowCells[i]->itEnd - treeSheetsGam->lowCells[i]->itBegin << std::endl;
				int aflN = it.getAflPnl().first;
				int panNum = it.getAflPnl().second;

				velI.assign(shDim, 0.0);

				for (size_t j = 0; j < treeSheetsGamRef.lowCells[i]->closeCells.size(); ++j)
					W.getAirfoil(aflN).GetInfluenceFromVorticesToPanel(panNum, &(treeSheetsGamRef.lowCells[i]->closeCells[j]->itBegin.getVtx()), \
						std::distance(treeSheetsGamRef.lowCells[i]->closeCells[j]->itBegin, treeSheetsGamRef.lowCells[i]->closeCells[j]->itEnd), velI);
				
				for (size_t j = 0; j < shDim; ++j)
					velI[j] *= IDPI / W.getAirfoil(aflN).len[panNum];

				// результат приближенного влияния записывается в velo конкретной панели (там уже есть накопленные влияния от дальних источников)
				treeSheetsGamRef.lowCells[i]->CalcInfluenceFromVortexFarCells(it);

				velI[0] += treeSheetsGamRef.allPnt.velo[it.getNum()] & W.getAirfoil(aflN).tau[panNum];
				
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

void VelocityBarnesHut::GetWakeInfluenceToRhsBS(const Airfoil& afl, std::vector<double>& wakeRhs) const
{
	size_t np = afl.getNumberOfPanels();
	int id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(np);

	std::vector<double> locVeloWake, locVeloWakeLin;
	locVeloWake.resize(par.myLen);

	size_t shDim = W.getBoundary(afl.numberInPassport).sheetDim;

	if (shDim != 1)
		locVeloWakeLin.resize(par.myLen);

	//локальные переменные для цикла	
	std::vector<double> velI(shDim, 0.0);

#pragma omp parallel for default(none) shared(locVeloWake, par, shDim, afl, locVeloWakeLin, IDPI) private(velI)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI.assign(shDim, 0.0);

		if (W.getWake().vtx.size() > 0)
		{
			//Учет влияния следа
			afl.GetInfluenceFromVorticesToPanel(par.myDisp + i, &(*W.getWake().vtx.begin()), std::distance(W.getWake().vtx.begin(), W.getWake().vtx.end()), velI);
		}

		if (W.getSource().vtx.size() > 0)
		{
			//Учет влияния источников
			afl.GetInfluenceFromSourcesToPanel(par.myDisp + i, &(*W.getSource().vtx.begin()), std::distance(W.getSource().vtx.begin(), W.getSource().vtx.end()), velI);
		}

		for (size_t j = 0; j < shDim; ++j)
			velI[j] *= IDPI / afl.len[par.myDisp + i];

		locVeloWake[i] = velI[0];

		if (shDim != 1)
			locVeloWakeLin[i] = velI[1];
	}

	if (id == 0)
		wakeRhs.resize(W.getBoundary(afl.numberInPassport).GetUnknownsSize());

	MPI_Gatherv(locVeloWake.data(), par.myLen, MPI_DOUBLE, wakeRhs.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

	if (shDim != 1)
		MPI_Gatherv(locVeloWakeLin.data(), par.myLen, MPI_DOUBLE, wakeRhs.data() + np, par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);
}//GetWakeInfluenceToRhsBS(...)

void VelocityBarnesHut::FillRhs(Eigen::VectorXd& rhs) const
{
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


#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_RHS))
		GPUGetWakeInfluenceToRhs(afl, wakeRhs);
#else
		GetWakeInfluenceToRhsBS(afl, wakeRhs);
#endif				

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
}

///Заполнение правой части СЛАУ 
/*void VelocityBarnesHut::FillRhs(Eigen::VectorXd& rhs) const
{

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
				rhs(currentRow + i) = -vInfRhs[i] -wakeRhs[currentRow + i] + 0.25 * ((afl.getV(i) + afl.getV(i + 1)) & afl.tau[i]);

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
	shared(afl, bou, wakeRhs, currentRow, rhs) \
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
	shared(afl, bou, wakeRhs, rhs, currentRow) \
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
}//FillRhs(...)
*/
//Вычисление числителей и знаменателей диффузионных скоростей в вихревом следе
void VelocityBarnesHut::CalcDiffVeloI1I2ToWakeFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	if (W.treeWake)
	{
		const int& id = W.getParallel().myidWork;
		if (id == 0)
		{
			double tCPUSTART, tCPUEND;

			tCPUSTART = omp_get_wtime();

			std::vector<double> selfI1;
			std::vector<Point2D> selfI2;

			selfI1.resize(pointsDb.vtx.size(), 0.0);
			selfI2.resize(pointsDb.vtx.size(), { 0.0, 0.0 });

#pragma warning (push)
#pragma warning (disable: 4101)
			//Локальные переменные для цикла
			Point2D Rij;
			double rij;
			double diffRadius, domRad;
			double posJx;
			double prd;
#pragma warning (pop)

			auto& tree = *(W.treeWake);
			std::vector<double> domRadinI;
			double distance;

#pragma omp parallel for private(Rij, rij, prd, diffRadius, domRad, posJx, domRadinI, distance) schedule(dynamic, 16)
			for (int s = 0; s < (int)tree.lowCells.size(); ++s)
			{
				domRadinI.resize(tree.lowCells[s]->itEnd - tree.lowCells[s]->itBegin);

				int cnt = 0;
				for (auto it = tree.lowCells[s]->itBegin; it != tree.lowCells[s]->itEnd; ++it)
				{
					int i = it.getNum();
					domRadinI[cnt] = std::max(domainRadius[i], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());
					++cnt;
				}
				/// \todo Понять природу магической константы 8.0 и синхронизировать с GPU
				distance = 8.0 * *std::max_element(domRadinI.begin(), domRadinI.end());

				tree.lowCells[s]->CalcCloseZone(tree.rootCell, distance);

				for (size_t k = 0; k < tree.lowCells[s]->closeCells.size(); ++k)
				{
					int cntr = 0;
					for (auto it = tree.lowCells[s]->itBegin; it != tree.lowCells[s]->itEnd; ++it)
					{
						const Vortex2D& vtxI = it.getVtx();
						int i = it.getNum();

						domRad = domRadinI[cntr];
						++cntr;

						/// \todo Понять природу магической константы 8.0 и синхронизировать с GPU
						diffRadius = 8.0 * domRad;


						for (auto itCl = tree.lowCells[s]->closeCells[k]->itBegin; itCl != tree.lowCells[s]->closeCells[k]->itEnd; ++itCl)
						{
							const Vortex2D& vtxJ = itCl.getVtx();
							posJx = vtxJ.r()[0];

							Rij = vtxI.r() - vtxJ.r();
							rij = Rij.length();

							if (rij < diffRadius && rij > 1.e-10)
							{
								prd = vtxJ.g() * exp(-rij / domRad);
								selfI2[i] += (prd / rij) * Rij;
								selfI1[i] += prd;
							}
						}//for j
					}
				}
			}

			for (size_t i = 0; i < I1.size(); ++i)
			{
				I1[i] += selfI1[i];
				I2[i] += selfI2[i];
			}

			tCPUEND = omp_get_wtime();
		}
	}
}//CalcDiffVeloI1I2ToWakeFromWake(...)

void VelocityBarnesHut::CalcDiffVeloI1I2ToWakeFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	if (W.treeWake)
	{
		const int& id = W.getParallel().myidWork;
		if (id == 0)
		{
			double tCPUSTART, tCPUEND;

			tCPUSTART = omp_get_wtime();

			std::vector<double> selfI1;
			std::vector<Point2D> selfI2;

			selfI1.resize(pointsDb.vtx.size(), 0.0);
			selfI2.resize(pointsDb.vtx.size(), { 0.0, 0.0 });

			//синхронизация свободного вихревого слоя ???
			bnd.sheets.FreeSheetSynchronize();

#pragma warning (push)
#pragma warning (disable: 4101)
			//Локальные переменные для цикла
			Point2D Rij;
			double rij;
			double diffRadius;
			//double left;
			//double right;
			double posJx;
			double domRad;
			double prd;
#pragma warning (pop)

			auto& tree = *(W.treeWake);
			std::vector<double> domRadinI;
			double distance;

#pragma omp parallel for private(Rij, rij, prd, diffRadius, domRad, posJx, domRadinI, distance) schedule(dynamic, 16)
			for (int s = 0; s < (int)tree.lowCells.size(); ++s)
			{
				domRadinI.resize(tree.lowCells[s]->itEnd - tree.lowCells[s]->itBegin);

				int cnt = 0;
				for (auto it = tree.lowCells[s]->itBegin; it != tree.lowCells[s]->itEnd; ++it)
				{
					int i = it.getNum();
					domRadinI[cnt] = std::max(domainRadius[i], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());
					++cnt;
				}
				/// \todo Понять природу магической константы 8.0 и синхронизировать с GPU		
				distance = 8.0 * *std::max_element(domRadinI.begin(), domRadinI.end());

				tree.lowCells[s]->closeCells.clear();
				tree.lowCells[s]->CalcCloseZone(W.treeSheetsGam->rootCell, distance);

				for (size_t k = 0; k < tree.lowCells[s]->closeCells.size(); ++k)
				{
					int cntr = 0;
					for (auto it = tree.lowCells[s]->itBegin; it != tree.lowCells[s]->itEnd; ++it)
					{
						const Vortex2D& vtxI = it.getVtx();
						int i = it.getNum();

						domRad = domRadinI[cntr];
						++cntr;

						/// \todo Понять природу магической константы 8.0 и синхронизировать с GPU
						diffRadius = 8.0 * domRad;

						for (auto itCl = tree.lowCells[s]->closeCells[k]->itBegin; itCl != tree.lowCells[s]->closeCells[k]->itEnd; ++itCl)
						{
							/// \todo Сделать переменной и синхронизировать с GPU
							const int nQuadPt = 3;
							int j = itCl.getAflPnl().second;
							auto& bndCl = W.getBoundary(itCl.getAflPnl().first);
							if (&bndCl == &bnd)
								/// \todo Учитываем пока только нулевой момент решения
							{
								const double ptG = bndCl.sheets.freeVortexSheet(j, 0) * bndCl.afl.len[j] / nQuadPt;

								for (int q = 0; q < nQuadPt; ++q)
								{
									const Point2D& ptJ = bndCl.afl.getR(j) + bndCl.afl.tau[j] * (q + 0.5) * bndCl.afl.len[j] * (1.0 / nQuadPt);  // vorticesDb.vtx[j];
									posJx = ptJ[0];

									Rij = vtxI.r() - ptJ;
									rij = Rij.length();
									if (rij < diffRadius && rij > 1.e-10)
									{
										prd = ptG * exp(-rij / domRad);
										selfI2[i] += (prd / rij) * Rij;
										selfI1[i] += prd;
									}
								}//for q
							}
						}//for itCL
					}//for it
				}//for k
			}//for s


			for (size_t i = 0; i < I1.size(); ++i)
			{
				I1[i] += selfI1[i];
				I2[i] += selfI2[i];
			}

			tCPUEND = omp_get_wtime();


		}//if(id == 0)
	}
}//CalcDiffVeloI1I2ToWakeFromSheets(...)
