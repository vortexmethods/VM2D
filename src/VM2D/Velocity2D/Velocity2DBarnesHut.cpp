/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
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
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
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
	W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
	
	Tree& treeWakeRef = W.getNonConstTree(PointType::wake);
	Tree& treeSheetsGamRef = W.getNonConstTree(PointType::sheetGam);
	//Tree& treeVPRef = W.getNonConstTree(PointType::wakeVP);
	Tree& treeSheetsSourceRef = W.getNonConstTree(PointType::source);
	Tree& treeSourceWakeRef = W.getNonConstTree(PointType::sourceWake);


	// WAKE:
	if (W.treeWake)
	{
#pragma omp parallel for /*reduction(+:timeCoeff, timeBiot, timeSum)*/ schedule(dynamic, 16)
		for (int i = 0; i < (int)(treeWakeRef.lowCells.size()); ++i)
		{
			std::vector<numvector<double, 3>> savedEe2;
			//todo
			// расчет влияния на wake от wake 
			// расчет влияния на wake от sheets 
			// расчет влияния на wake от sources 
			// расчет влияния на wake от sourcesWake 
		}
	}//if(W.treeWake)

	
	// расчет скоростей самих виртуальных вихрей
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
		W.getBoundary(bou).CalcConvVelocityAtVirtualVortexes(virtualVortexesParams[bou].convVelo);
	
	W.getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();

	// вычисление eps* для виртуальных вихрей 
	W.getTimestat().timeCalcVortexDiffVelo.first += omp_get_wtime();
	//todo

	//синхронизация с CUDA после расчетов eps* для виртуальных вихрей
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
	//todo Аналогично для wakeVP
	W.getTimestat().timeVP.second += omp_get_wtime();

	W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
	
	//раздача копий в оригинальные значения 
	if (W.treeWake)
	{
		//todo

		//синхронизация с CUDA после расчетов eps* для следа
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

	//Здесь предполагается, что тип схемы у всех профилей одинаковый
	for (size_t i = 0; i < treeSheetsGamRef.lowCells.size(); ++i)
	{

		if (W.treeSheetsSource)
		{
			//todo
			// расчет влияния на панели от sources 
		}//if(treeSources)


		// аналогично для wake		
		if (W.treeWake)
		{			
			//todo
			// расчет влияния на панели от wake 
		}//if (treeWake)
	}

}//GetWakeInfluenceToRhs(...)

///Заполнение правой части СЛАУ 
void VelocityBarnesHut::FillRhs(Eigen::VectorXd& rhs) const
{

	std::vector<double> wakeRhs;
	wakeRhs.resize(rhs.size());
//\todo реализовать GPU
	GetWakeInfluenceToRhs(wakeRhs);

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

//Вычисление числителей и знаменателей диффузионных скоростей в вихревом следе
void VelocityBarnesHut::CalcDiffVeloI1I2ToWakeFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	//todo
}//CalcDiffVeloI1I2ToWakeFromWake(...)

void VelocityBarnesHut::CalcDiffVeloI1I2ToWakeFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	//todo
}//CalcDiffVeloI1I2ToWakeFromSheets(...)
