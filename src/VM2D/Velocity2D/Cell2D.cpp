/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Cell2D.cpp                                                       |
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
\brief Файл кода с описанием класса Cell
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
*/

#include "Cell2D.h"
#include <algorithm>

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Tree2D.h"
#include "Velocity2D.h"
#include "Velocity2DBarnesHut.h"
#include "Wake2D.h"
#include "World2D.h"
#include <fstream>
using std::ifstream;
using std::ofstream;


using namespace VM2D;

/// Конструктор
/*
Cell::Cell(Tree& Tree_)
	: tree(Tree_)
{
	leftDown = { 0.0, 0.0 };
	rightUp = { 0.0, 0.0 };

	for (int i = 0; i < 2; ++i)
	{
		sumGam[i] = 0.0;
		centerGam[i] = { 0.0, 0.0 };
		mChildren[i] = nullptr;
		level = 0;
	}//for i

	a = 0.0;
	b = 0.0;
	c = 0.0;
	d = 0.0;

};

/// Конструктор
Cell::Cell(Tree& Tree_, const Point2D leftDown_, const Point2D rightUp_)
	: leftDown(leftDown_), rightUp(rightUp_), tree(Tree_)
{

}

/// Деструктор
Cell::~Cell()
{

};

void Cell::CreateRootCell(PointsCopy& pointsCopy)
{
	auto spanY = std::minmax_element(pointsCopy.begin(), pointsCopy.end(), PointsCopy_iterator::LessYrr);
	auto spanX = std::minmax_element(pointsCopy.begin(), pointsCopy.end(), PointsCopy_iterator::LessXrr);

	leftDown[1] = spanY.first.getVtx().r()[1];
	rightUp[1] = spanY.second.getVtx().r()[1];

	leftDown[0] = spanX.first.getVtx().r()[0];
	rightUp[0] = spanX.second.getVtx().r()[0];

	itBegin = pointsCopy.begin();
	itEnd = pointsCopy.end();
}//CreateRootCell(...)

// Построение иерархической структуры ячеек (дерева)
void Cell::CreateAllCells(PointsCopy::iterator& itBegin, PointsCopy::iterator& itEnd)
{
	double w = rightUp[0] - leftDown[0];
	double h = rightUp[1] - leftDown[1];

	if (std::distance(itBegin, itEnd) > 1)
	{
		if (w > h)
		{
			double midX = (rightUp[0] + leftDown[0]) / 2.0;

			mChildren[0].reset(new Cell(tree, leftDown, { midX, rightUp[1] }));
			mChildren[1].reset(new Cell(tree, { midX, leftDown[0] }, rightUp));

			// сортировка по х
			// НЕСОВМЕСТИМОСТЬ
			//std::sort(itBegin, itEnd, PointsCopy_iterator::LessXcc);

			PointsCopy::iterator itDop(itBegin);

			mChildren[0]->itBegin = itBegin;
			mChildren[1]->itEnd = itEnd;

			bool flag = true;

			while (flag)
			{
				if (itDop.getVtx().r()[0] < midX)
					itDop++;
				else
				{
					flag = false;
					mChildren[0]->itEnd = itDop;
					mChildren[1]->itBegin = itDop;

					//обрезаем по вихрям
					mChildren[0]->rightUp[0] = (itDop - 1).getVtx().r()[0];
					mChildren[1]->leftDown[0] = itDop.getVtx().r()[0];
				}
			}//while (flag)

			// для границ по горизонтали
			for (int i = 0; i < 2; ++i)
			{
				auto spanY = std::minmax_element(mChildren[i]->itBegin, mChildren[i]->itEnd, PointsCopy_iterator::LessYrr);

				mChildren[i]->leftDown[1] = spanY.first.getVtx().r()[1];
				mChildren[i]->rightUp[1] = spanY.second.getVtx().r()[1];
			}// for i
		} //if (w > h)
		else
		{
			double midY = (leftDown[1] + rightUp[1]) / 2.0;

			mChildren[0].reset(new Cell(tree, leftDown, { rightUp[0], midY }));
			mChildren[1].reset(new Cell(tree, { leftDown[0], midY }, rightUp));

			//сортировка по у 
			// НЕСОВМЕСТИМОСТЬ
			//std::sort(itBegin, itEnd, PointsCopy_iterator::LessYcc);

			PointsCopy::iterator itDop = itBegin;

			mChildren[0]->itBegin = itBegin;
			mChildren[1]->itEnd = itEnd;

			bool flag = true;

			while (flag)
			{
				if (itDop.getVtx().r()[1] < midY)
					itDop++;
				else
				{
					flag = false;
					mChildren[0]->itEnd = itDop;
					mChildren[1]->itBegin = itDop;

					//обрезаем по вихрям
					mChildren[0]->rightUp[1] = (itDop - 1).getVtx().r()[1];
					mChildren[1]->leftDown[1] = itDop.getVtx().r()[1];
				}
			}//while (flag)


			//для границ по вертикали
			for (int i = 0; i < 2; ++i)
			{
				auto spanX = std::minmax_element(mChildren[i]->itBegin, mChildren[i]->itEnd, PointsCopy_iterator::LessXrr);

				mChildren[i]->leftDown[0] = spanX.first.getVtx().r()[0];
				mChildren[i]->rightUp[0] = spanX.second.getVtx().r()[0];
			}
		} //else

		// номер уровня потомков 
		for (int j = 0; j < 2; ++j)
		{
			mChildren[j]->level = level + 1;

			if ((mChildren[j]->level < tree.numOfLevels) && ((mChildren[j]->itEnd - mChildren[j]->itBegin) > tree.minNumOfVort) && \
				(std::max(fabs(mChildren[j]->rightUp[0] - mChildren[j]->leftDown[0]), fabs(mChildren[j]->rightUp[1] - mChildren[j]->leftDown[1])) > tree.minDist))
				mChildren[j]->CreateAllCells(mChildren[j]->itBegin, mChildren[j]->itEnd);
			else
			{
				mChildren[j]->lowLevel = true;
				tree.lowCells.push_back(mChildren[j].get());
			}
		}
	}
	else
	{
		lowLevel = true;
		tree.lowCells.push_back(this);
	}
}//CreateAllCells(...)

//Вычисление параметров всех ячеек дерева (циркуляций и центров завихренности)	
void Cell::CalculateCellsParams()
{
	double gg[2];
	Point2D rr[2];

	for (int i = 0; i < 2; ++i)
	{
		gg[i] = 0.0;
		rr[i].toZero();
	}

	if (!lowLevel)
	{
		mChildren[0]->CalculateCellsParams();
		mChildren[1]->CalculateCellsParams();

		for (int k = 0; k < 2; ++k)
		{
			gg[k] = mChildren[0]->sumGam[k] + mChildren[1]->sumGam[k];

			if (fabs(gg[k]) > 1e-10) //когда все вихри внутри области одной циркуляции
			{
				rr[k] += mChildren[0]->sumGam[k] * mChildren[0]->centerGam[k] + mChildren[1]->sumGam[k] * mChildren[1]->centerGam[k];
				rr[k] /= gg[k];
			}
			sumGam[k] = gg[k];
			centerGam[k] = rr[k];
			centerGam[k] = rr[k];
		}// for k
	}//if(!lowLevel)

	else
	{
		for(auto it = itBegin; it != itEnd; ++it)
		{
			double gDop = it.getVtx().g();
			int k = (int)(gDop < 0.0);
			gg[k] += gDop;
			rr[k] += gDop * it.getVtx().r();
		}

		for (int k = 0; k < 2; ++k)
		{
			if (fabs(gg[k]) > 1e-10)
				rr[k] /= gg[k];

			sumGam[k] = gg[k];
			centerGam[k] = rr[k];			
		}// for k
	}//else
}//CalculateCellsParams()

//Обнуление коэффициентов a,b,c,d для ячейки нижнего уровня 
void Cell::ClearABCD()
{
	a = b = c = d = 0.0;
}// ClearABCD()

void Cell::CalcCloseZone(Cell& infCell, double distance)
{
	if (&infCell != this)
	{
		//сумма габаритов ячейки, на которую считается влияние (h0), и влияющей (h)
		double h0, h;
		h0 = fabs(rightUp[0] - leftDown[0]) + fabs(rightUp[1] - leftDown[1]);
		h = fabs(infCell.rightUp[0] - infCell.leftDown[0]) + fabs(infCell.rightUp[1] - infCell.leftDown[1]);

		//r0 --- нижний уровень, r1 --- влияющий
		Point2D r0, r1;

		// центры прямоугольников
		r0 = 0.5 * (rightUp + leftDown);
		r1 = 0.5 * (infCell.rightUp + infCell.leftDown);

		double dc = dist(r0, r1);

		// если выполнен критерий дальности => считаем коэффициенты
		if (dc - 0.5 * (h + h0) > distance)
		{
			
		}//if crit
		else // если не выполнен критерий, то рекурсия 
		{
			if (!(infCell.lowLevel))
			{
				CalcCloseZone(*(infCell.mChildren[0]), distance);
				CalcCloseZone(*(infCell.mChildren[1]), distance);
			}
			else
				closeCells.push_back(&infCell);
		}
	}//if (lowTree != this)
	else closeCells.push_back(&infCell); 

}//CalcCloseZone(...)

//Расчет коэффициентов a,b,c,d для ячейки нижнего уровня
void Cell::CalcABCDandCloseCellsToLowLevel(Cell& infCell, bool calcCoef)
{
	if (&infCell != this)
	{
		//сумма габаритов ячейки, на которую считается влияние (h0), и влияющей (h)
		double h0, h;
		h0 = fabs(rightUp[0] - leftDown[0]) + fabs(rightUp[1] - leftDown[1]);
		h = fabs(infCell.rightUp[0] - infCell.leftDown[0]) + fabs(infCell.rightUp[1] - infCell.leftDown[1]);

		//r0 --- нижний уровень, r1 --- влияющий, p --- с плюсом, m --- с минусом  
		Point2D r0, r1, r1p, r1m;

		// центры прямоугольников
		r0 = 0.5 * (rightUp + leftDown);
		r1 = 0.5 * (infCell.rightUp + infCell.leftDown);

		double crit = dist(r0, r1);
		//double crit = sqrt(sqr(r0[0] - r1[0]) + sqr(r0[1] - r1[1]));
		
		// если выполнен критерий дальности => считаем коэффициенты
		if ((crit >= (h0 + h + tree.W.getPassport().wakeDiscretizationProperties.eps) / tree.theta))
		{ 
			if (calcCoef)
			{
				if ((itEnd - itBegin) <= 3)
					VMlib::ModifyE2(closestCellDist.data(), crit);

				//центры завихренностей
				r1p = infCell.centerGam[0];
				r1m = infCell.centerGam[1];

				//суммарные циркуляции
				double gam[2];
				gam[0] = infCell.sumGam[0];
				gam[1] = infCell.sumGam[1];

				Point2D  rrP, rrM;
				rrP = r0 - r1p;
				rrM = r0 - r1m;

				double nrmP2 = rrP.length2();
				double nrmM2 = rrM.length2();

				a += -IDPI * (gam[0] * rrP[1] / nrmP2 + gam[1] * rrM[1] / nrmM2);
				b += IDPI * (gam[0] * rrP[0] / nrmP2 + gam[1] * rrM[0] / nrmM2);

				double nrmP4 = nrmP2 * nrmP2;
				double nrmM4 = nrmM2 * nrmM2;

				c += 2.0 * IDPI * (gam[0] * rrP[0] * rrP[1] / nrmP4 + \
					gam[1] * rrM[0] * rrM[1] / nrmM4);

				d += IDPI * (gam[0] * (rrP[1] * rrP[1] - rrP[0] * rrP[0]) / nrmP4 + \
					gam[1] * (rrM[1] * rrM[1] - rrM[0] * rrM[0]) / nrmM4);
			}
		}//if crit
		else // если не выполнен критерий, то рекурсия 
		{
			if (!(infCell.lowLevel))
			{
				CalcABCDandCloseCellsToLowLevel(*(infCell.mChildren[0]), calcCoef);
				CalcABCDandCloseCellsToLowLevel(*(infCell.mChildren[1]), calcCoef);
			}
			else
				closeCells.push_back(&infCell);
		}
	}//if (lowTree != this)
	else closeCells.push_back(&infCell); //себя тоже добавляем в ближнюю зону 
	
}//CalcABCDandCloseCellsToLowLevel(...)

//Расчет конвективных скоростей вихрей внутри одной ячейки от вихрей внутри всех ближних ячеек
void Cell::CalcConvVeloByBiotSavartFromVortices(bool calcRadius, std::vector<numvector<double, 3>>& savedEe2)
{
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2eps, dst2;

	if(calcRadius)
		savedEe2.reserve((int)(itEnd - itBegin));

	for (auto it = itBegin; it != itEnd; ++it)
	{
		numvector<double, 3> ee2;
		ee2[0] = ee2[1] = ee2[2] = 10000.0;

		velI.toZero();

		const Point2D &posI = it.getVtx().r();

		for (size_t k = 0; k < closeCells.size(); ++k)		
		for (auto it2 = closeCells[k]->itBegin; it2 != closeCells[k]->itEnd; ++it2)
		{
			const Point2D& posJ = it2.getVtx().r();
						
			dst2 = dist2(posI, posJ);

			//Модифицируем массив квадратов расстояний до ближайших  вихрей из wake
			if (calcRadius)
				VMlib::ModifyE2(&ee2[0], dst2);

			const double& gamJ = it2.getVtx().g();

			tempVel.toZero();
			dst2eps = VMlib::boundDenom(dst2, tree.W.getPassport().wakeDiscretizationProperties.eps2); //Сглаживать надо!!!

			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2eps);
			velI += tempVel;			
		}//for it2		
	
		velI *= IDPI;

		tree.allPnt.velo[it.getNum()] += velI;
		
		if (calcRadius)
		{
			for (size_t i = 0; i < ee2.size(); ++i)
				VMlib::ModifyE2(ee2.data(), sqr(closestCellDist[i]));

			savedEe2.push_back(ee2);
		}
	}//for it
}//CalcConvVeloByBiotSavartFromVortices(...)

//Расчет eps* для виртуальных вихрей от всех вихрей (в пелене и виртуальных) внутри всех ближних ячеек
void Cell::CalcDomainRadiusForVirtualWake()
{
	double dst2;

	std::vector<numvector<double,3>> savedEe2; 

	int numVtx = 0;
	for (auto it = itBegin; it != itEnd; ++it)
	{
		numVtx += tree.W.getBoundary(it.getAflPnl().first).vortexBeginEnd[it.getAflPnl().second].second - \
			tree.W.getBoundary(it.getAflPnl().first).vortexBeginEnd[it.getAflPnl().second].first;
	}
	savedEe2.reserve(numVtx);


	for (auto it = itBegin; it != itEnd; ++it) //цикл по панелям
	{
		const Boundary& bou = tree.W.getBoundary(it.getAflPnl().first);
		for (int virt = bou.vortexBeginEnd[it.getAflPnl().second].first; virt < bou.vortexBeginEnd[it.getAflPnl().second].second; ++virt)
		{
			numvector<double, 3> tempEe2;
			tempEe2[0] = tempEe2[1] = tempEe2[2] = 10000.;
			const Point2D &posI = bou.virtualWake.vtx[virt].r();

			// В ближней зоне ячеек treeSheetsGam уже есть wake
			for(size_t k = 0; k < closeCells.size(); ++k)
				for (auto it2 = closeCells[k]->itBegin; it2 != closeCells[k]->itEnd; ++it2)
				{
					const Point2D& posJ = it2.getVtx().r();

					dst2 = dist2(posI, posJ);

					//Модифицируем массив квадратов расстояний до ближайших вихрей из wake
						VMlib::ModifyE2(&tempEe2[0], dst2);
				}//for it2

			for (size_t i = 0; i < tempEe2.size(); ++i)
				VMlib::ModifyE2(tempEe2.data(), sqr(closestCellDist[i]));

			savedEe2.push_back(tempEe2);			
		}// for virt
	}// for it

	closeCells.clear();

	CalcABCDandCloseCellsToLowLevel(tree.rootCell, false);

	int currI = 0;
	for (auto it = itBegin; it != itEnd; ++it) //цикл по панелям
	{
		const Boundary& bou = tree.W.getBoundary(it.getAflPnl().first);
		for (int virt = bou.vortexBeginEnd[it.getAflPnl().second].first; virt < bou.vortexBeginEnd[it.getAflPnl().second].second; ++virt)
		{
			const Point2D &posI = bou.virtualWake.vtx[virt].r();

			for (size_t k = 0; k < closeCells.size(); ++k)
				for (auto it2 = closeCells[k]->itBegin; it2 != closeCells[k]->itEnd; ++it2)
				{
					const Boundary& bou2 = tree.W.getBoundary(it2.getAflPnl().first);
					for (int virt2 = bou2.vortexBeginEnd[it2.getAflPnl().second].first; virt2 < bou2.vortexBeginEnd[it2.getAflPnl().second].second; ++virt2)
					{
						const Point2D& posJ = bou2.virtualWake.vtx[virt2].r();

						dst2 = dist2(posI, posJ);

						//Модифицируем массив квадратов расстояний до ближайших вихрей из wake
						VMlib::ModifyE2(&savedEe2[currI][0], dst2);

					}
				}//for it2			
			
			for (size_t i = 0; i < savedEe2[currI].size(); ++i)
				VMlib::ModifyE2(savedEe2[currI].data(), sqr(closestCellDist[i]));			
			
			tree.W.getNonConstVelocity().virtualVortexesParams[it.getAflPnl().first].epsastWake[virt] = \
				sqrt((savedEe2[currI][0] + savedEe2[currI][1] + savedEe2[currI][2]) / 3.0);
			currI++;
		}// for virt
	}// for it

}//CalcDomainRadiusForVirtualWake()


//Расчет конвективных скоростей вихрей внутри одной ячейки от источников внутри всех ближних ячеек
void Cell::CalcConvVeloByBiotSavartFromSources()
{
	//Локальные переменные для цикла
	Point2D velI, velSheet;
	Point2D tempVel;
	double dst2eps, dst2;

	for (auto it = itBegin; it != itEnd; ++it)
	{
		velI.toZero();

		const Point2D &posI = it.getVtx().r();

		for (size_t k = 0; k < closeCells.size(); ++k)
			for (auto it2 = closeCells[k]->itBegin; it2 != closeCells[k]->itEnd; ++it2)
			{
				const Point2D& posJ = it2.getVtx().r();

				dst2 = dist2(posI, posJ);
				
				if (closeCells[k]->tree.allPnt.type == PointType::sourceWake)
				{
					const double& gamJ = it2.getVtx().g();

					tempVel.toZero();
					dst2eps = VMlib::boundDenom(dst2, tree.W.getPassport().wakeDiscretizationProperties.eps2); //Сглаживать надо!!!

					tempVel = { posI[0] - posJ[0], posI[1] - posJ[1] };
					tempVel *= (gamJ / dst2eps);
					velI += tempVel;
				}
				else
				{
					tree.W.getAirfoil(it2.getAflPnl().first).GetInfluenceFromSourceSheetToVortex(it2.getAflPnl().second, it.getVtx(), velSheet);
					velI += velSheet;
				}
			}//for it2


		velI *= IDPI;
		tree.allPnt.velo[it.getNum()] += velI; 
		
	}//for it
}//CalcConvVeloByBiotSavartFromSources()

// Расчет конвективных скоростей вихрей внутри одной ячейки от вихревых слоев (присоединенный + свободный) внутри всех ближних ячеек
void Cell::CalcConvVeloByBiotSavartFromSheets(bool calcRadius, std::vector<numvector<double, 3>>& savedEe2)
{
	//Локальные переменные для цикла
	Point2D velI, velSheet;
	Point2D tempVel;
	double dst2;

	int currI = 0;

	for (auto it = itBegin; it != itEnd; ++it)
	{
		velI.toZero();

		const Point2D &posI = it.getVtx().r();	

		for (size_t k = 0; k < closeCells.size(); ++k)
			for (auto it2 = closeCells[k]->itBegin; it2 != closeCells[k]->itEnd; ++it2)//цикл по панелям 
			{
				const Boundary& bou = tree.W.getBoundary(it2.getAflPnl().first);


				tree.W.getAirfoil(it2.getAflPnl().first).GetInfluenceFromVortexSheetToVortex(it2.getAflPnl().second, it.getVtx(), velSheet);
				velI += velSheet;		
				
				if(calcRadius)
				for (int virt = bou.vortexBeginEnd[it2.getAflPnl().second].first; virt < bou.vortexBeginEnd[it2.getAflPnl().second].second; ++virt)
				{
					dst2 = dist2(posI, bou.virtualWake.vtx[virt].r());

					//Модифицируем массив квадратов расстояний до ближайших вихрей из virtual wake
					VMlib::ModifyE2(savedEe2[currI].data(), dst2);

				}

			}//for it2

		velI *= IDPI;

		tree.allPnt.velo[it.getNum()] += velI;

		currI++;
	}//for it
}//CalcConvVeloByBiotSavartFromSheets()


void Cell::SetDomainRadius(const std::vector<numvector<double, 3>>& savedEe2)
{
	int i = 0;
	for (auto it = itBegin; it != itEnd; ++it, ++i)
		tree.allPnt.domainRadius[it.getNum()] = sqrt((savedEe2[i][0] + savedEe2[i][1] + savedEe2[i][2]) / 3.0);
	
}

// Расчет приближенного влияния от вихрей на элемент it внутри ячейки нижнего уровня
void Cell::CalcInfluenceFromVortexFarCells(PointsCopy::iterator it)
{
	Point2D Rc;
	Rc[0] = 0.5*(rightUp[0] + leftDown[0]);
	Rc[1] = 0.5*(rightUp[1] + leftDown[1]);

	Point2D deltaPos = it.getVtx().r() - Rc;

	tree.allPnt.velo[it.getNum()][0] += a + c * deltaPos[0] + d * deltaPos[1];
	tree.allPnt.velo[it.getNum()][1] += b + d * deltaPos[0] - c * deltaPos[1];

}//CalcInfluenceFromFarCells(...)

// Расчет приближенного влияния от источников на элемент it внутри ячейки нижнего уровня
void Cell::CalcInfluenceFromSourceFarCells(PointsCopy::iterator it)
{
	Point2D Rc;
	Rc[0] = 0.5*(rightUp[0] + leftDown[0]);
	Rc[1] = 0.5*(rightUp[1] + leftDown[1]);

	Point2D deltaPos = it.getVtx().r() - Rc;

	//\todo проверить формулы 
	tree.allPnt.velo[it.getNum()][0] +=  b + d * deltaPos[0] - c * deltaPos[1];
	tree.allPnt.velo[it.getNum()][1] += -a - c * deltaPos[0] - d * deltaPos[1];

}//CalcInfluenceFromFarCells(...)

/// Печать всего дерева в текстовый файл
void Cell::PrintTree()
{
	ofstream outfile;
	if (level == 0)
		outfile.open(tree.W.getPassport().dir+"tree.txt");
	else outfile.open(tree.W.getPassport().dir + "tree.txt", std::ios_base::app);
	
	if (lowLevel)
		outfile << leftDown[0] << " " << leftDown[1] << " " << rightUp[0] << " " << rightUp[1] << std::endl;

	outfile.close();
	if (!lowLevel)
	{
		mChildren[0]->PrintTree();
		mChildren[1]->PrintTree();
		return;
	}//if (!lowLevel)
	else
		return;
}//PrintAllTree()
*/
