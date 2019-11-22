/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
\brief Файл кода с описанием класса Cell
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7
\date 22 ноября 2019 г.
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
Cell::Cell(const Point2D leftDown_, const Point2D rightUp_, Tree& Tree_)
	: leftDown(leftDown_), rightUp(rightUp_), tree(Tree_)
{

}

/// Деструктор
Cell::~Cell()
{

};

// Построение корня дерева на основе заданных вихрей
void Cell::MakeRootCell(std::vector<PointsCopy>& pointsCopy) 
{
	//сортировка по у
	std::sort(pointsCopy.begin(), pointsCopy.end(), [](const PointsCopy &a, const PointsCopy &b) { return a.r()[1]  < b.r()[1]; });

	size_t numElem = pointsCopy.size();

	leftDown[1] = pointsCopy[0].r()[1];
	rightUp[1] = pointsCopy[numElem - 1].r()[1];

	//сортировка по x
	std::sort(pointsCopy.begin(), pointsCopy.end(), [](const PointsCopy &a, const PointsCopy &b) { return a.r()[0] < b.r()[0]; });

	leftDown[0] = pointsCopy[0].r()[0];
	rightUp[0] = pointsCopy[numElem - 1].r()[0];

	itBegin = pointsCopy.begin();
	itEnd = pointsCopy.end();

	level = 0;
}//MakeRootCell(...)

// Построение иерархической структуры ячеек (дерева)
void Cell::CreateAllCells(std::vector<PointsCopy>::iterator &itBegin, std::vector<PointsCopy>::iterator &itEnd)
{
	double w = rightUp[0] - leftDown[0];
	double h = rightUp[1] - leftDown[1];

	if (w > h)
	{
		double midX = (rightUp[0] + leftDown[0]) / 2.0;

		mChildren[0].reset(new Cell(leftDown, { midX, rightUp[1] }, tree));
		mChildren[1].reset(new Cell({ midX, leftDown[0] }, rightUp, tree));

		// сортировка по х
		std::sort(itBegin, itEnd, [](const PointsCopy &a, const PointsCopy &b) { return a.r()[0] < b.r()[0]; });

		std::vector<PointsCopy>::iterator itDop(itBegin);

		mChildren[0]->itBegin = itBegin;
		mChildren[1]->itEnd = itEnd;

		bool flag = true;

		while (flag)
		{
			if ((*itDop).r()[0] < midX)
				itDop++;
			else
			{
				flag = false;
				mChildren[0]->itEnd = itDop;
				mChildren[1]->itBegin = itDop;

				//обрезаем по вихрям
				mChildren[0]->rightUp[0] = (*(itDop - 1)).r()[0];
				mChildren[1]->leftDown[0] = (*itDop).r()[0];
			}
		}//while (flag)

		// для границ по горизонтали
		for (int i = 0; i < 2; ++i)
		{
			//сортировка по у 
			std::sort(mChildren[i]->itBegin, mChildren[i]->itEnd, [](const PointsCopy &a, const PointsCopy &b) { return a.r()[1] < b.r()[1]; });

			mChildren[i]->leftDown[1] = mChildren[i]->itBegin->r()[1];
			mChildren[i]->rightUp[1] = ((mChildren[i]->itEnd) - 1)->r()[1];

		}// for i
	} //if (w > h)
	else
	{
		double midY = (leftDown[1] + rightUp[1]) / 2.0;

		mChildren[0].reset(new Cell(leftDown, { rightUp[0], midY }, tree));
		mChildren[1].reset(new Cell({ leftDown[0], midY }, rightUp, tree));

		//сортировка по у 
		std::sort(itBegin, itEnd, [](const PointsCopy &a, const PointsCopy &b) { return a.r()[1] < b.r()[1]; });

		std::vector<PointsCopy>::iterator itDop = itBegin;

		mChildren[0]->itBegin = itBegin;
		mChildren[1]->itEnd = itEnd;

		bool flag = true;
		
		while (flag)
		{
			if ((*itDop).r()[1] < midY)
				itDop++;		
			else
			{
				flag = false;
				mChildren[0]->itEnd = itDop;
				mChildren[1]->itBegin = itDop;

				//обрезаем по вихрям
				mChildren[0]->rightUp[1] = (*(itDop - 1)).r()[1];
				mChildren[1]->leftDown[1] = (*itDop).r()[1];
			}
		}//while (flag)


		//для границ по вертикали
		for (int i = 0; i < 2; ++i)
		{
			// сортировка по х
			std::sort(mChildren[i]->itBegin, mChildren[i]->itEnd, [](const PointsCopy &a, const PointsCopy &b) { return a.r()[0] < b.r()[0]; });

			mChildren[i]->leftDown[0] = mChildren[i]->itBegin->r()[0];
			mChildren[i]->rightUp[0] = ((mChildren[i]->itEnd) - 1)->r()[0];
		}
	} //else

	// номер уровня потомков 
	for (int j = 0; j < 2; ++j)
		mChildren[j]->level = level + 1;



	if ((mChildren[0]->level < tree.numOfLevels) && ((mChildren[0]->itEnd - mChildren[0]->itBegin) > tree.minNumOfVort) && \
		(fabs(mChildren[0]->rightUp[0] - mChildren[0]->leftDown[0]) + fabs(mChildren[0]->rightUp[1] - mChildren[0]->leftDown[1]) > tree.minDist))
			mChildren[0]->CreateAllCells(mChildren[0]->itBegin, mChildren[0]->itEnd);
	else
	{
		mChildren[0]->lowLevel = true;
		tree.lowCells.push_back(mChildren[0]);
	}

	if ((mChildren[1]->level < tree.numOfLevels) && ((mChildren[1]->itEnd - mChildren[1]->itBegin) > tree.minNumOfVort) && \
		(fabs(mChildren[1]->rightUp[0] - mChildren[1]->leftDown[0]) + fabs(mChildren[1]->rightUp[1] - mChildren[1]->leftDown[1]) > tree.minDist))
			mChildren[1]->CreateAllCells(mChildren[1]->itBegin, mChildren[1]->itEnd);
	else
	{
		mChildren[1]->lowLevel = true;
		tree.lowCells.push_back(mChildren[1]);
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
			double gDop = it->g();
			int k = (int)(gDop < 0.0);
			gg[k] += gDop;
			rr[k] += gDop * it->r();
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

//Расчет коэффициентов a,b,c,d для ячейки нижнего уровня
void Cell::CalcCoeffToLowLevel(std::shared_ptr<Cell> infCell)
{
	if (&(*infCell) != this)
	{
		//сумма габаритов ячейки, на которую считается влияние (h0), и влияющей (h)
		double h0, h; 
		h0 = fabs(rightUp[0] - leftDown[0]) + fabs(rightUp[1] - leftDown[1]);
		h = fabs(infCell->rightUp[0] - infCell->leftDown[0]) + fabs(infCell->rightUp[1] - infCell->leftDown[1]);

		//r0 --- нижний уровень, r1 --- влияющий, p --- с плюсом, m --- с минусом  
		Point2D r0, r1, r1p, r1m;

		// центры прямоугольников
		r0 = 0.5*(rightUp + leftDown);
		r1 = 0.5*(infCell->rightUp + infCell->leftDown);

		//центры завихренностей
		r1p = infCell->centerGam[0];
		r1m = infCell->centerGam[1];

		//суммарные циркуляции
		double gam[2];
		gam[0] = infCell->sumGam[0];
		gam[1] = infCell->sumGam[1];

		double crit = dist(r0, r1);

		// если выполнен критерий дальности => считаем коэффициенты
		if (crit >= (h0 + h + tree.W.getPassport().wakeDiscretizationProperties.eps) / tree.theta)
		{
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
		}//if crit
		else // если не выполнен критерий, то рекурсия 
		{
			if (!(infCell->lowLevel))
			{
				CalcCoeffToLowLevel(infCell->mChildren[0]);
				CalcCoeffToLowLevel(infCell->mChildren[1]);
			}
			else
				closeCells.push_back(infCell);
		}
	}//if (lowTree != this)
	else closeCells.push_back(infCell); //себя тоже добавляем в ближнюю зону 
}//CalcCoeffToLowLevel(...)

//Расчет конвективных скоростей вихрей внутри одной ячейки от вихрей внутри всех ближних ячеек
void Cell::CalcConvVeloByBiotSavart()
{
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2eps, dst2;

	for (auto it = itBegin; it != itEnd; ++it)
	{
		double ee2[3] = { 10000.0, 10000.0, 10000.0 };

		velI.toZero();

		const Point2D &posI = it->r();

		for (size_t k = 0; k < closeCells.size(); ++k)		
		for (auto it2 = closeCells[k]->itBegin; it2 != closeCells[k]->itEnd; ++it2)
		{
			const Point2D& posJ = it2->r();
						
			dst2 = dist2(posI, posJ);

			//Модифицируем массив квадратов расстояний до ближайших вихрей из wake
			VMlib::ModifyE2(ee2, dst2);

			//if ((it->type != PointType::virtwake) && (it2->type != PointType::virtwake))
			{
				const double& gamJ = it2->g();

				tempVel.toZero();
				dst2eps = VMlib::boundDenom(dst2, tree.W.getPassport().wakeDiscretizationProperties.eps2); //Сглаживать надо!!!

				tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
				tempVel *= (gamJ / dst2eps);
				velI += tempVel;
			}				
		}//for it2
				
		//if (it->type != PointType::virtwake)
			for (size_t j = 0; j < tree.W.getSource().vtx.size(); ++j)
			{
				const Point2D& posJ = tree.W.getSource().vtx[j].r();
				const double& gamJ = tree.W.getPassport().physicalProperties.accelCft() * tree.W.getSource().vtx[j].g();

				tempVel.toZero();

				dst2 = dist2(posI, posJ);
				dst2eps = VMlib::boundDenom(dst2, tree.W.getPassport().wakeDiscretizationProperties.eps2); //Сглаживать надо!!!

				tempVel = { posI[0] - posJ[0], posI[1] - posJ[1] };
				tempVel *= (gamJ / dst2eps);
				velI += tempVel;
			}

		//if (it->type != PointType::virtwake)
		{
			velI *= IDPI;
			it->veloCopy += velI;
		}
		
		it->domainRadiusCopy = std::max(sqrt((ee2[0] + ee2[1] + ee2[2]) / 3.0), 2.0 * tree.W.getPassport().wakeDiscretizationProperties.epscol);

	}//for it
}//CalcConvByBiotSavart()

void Cell::CalcInfluenceFromVortexFarCells(std::vector<PointsCopy>::iterator it)
{
	Point2D Rc;
	Rc[0] = 0.5*(rightUp[0] + leftDown[0]);
	Rc[1] = 0.5*(rightUp[1] + leftDown[1]);

	Point2D deltaPos = it->r() - Rc;

	it->veloCopy[0] += a + c * deltaPos[0] + d * deltaPos[1];
	it->veloCopy[1] += b + d * deltaPos[0] - c * deltaPos[1];
	
}//CalcInfluenceFromFarCells(...)

void Cell::CalcInfluenceFromSourceFarCells(std::vector<PointsCopy>::iterator it)
{
	Point2D Rc;
	Rc[0] = 0.5*(rightUp[0] + leftDown[0]);
	Rc[1] = 0.5*(rightUp[1] + leftDown[1]);

	Point2D deltaPos = it->r() - Rc;

	//\todo проверить формулы 
	it->veloCopy[0] +=  b + d * deltaPos[0] - c * deltaPos[1];  
	it->veloCopy[1] += -a - c * deltaPos[0] - d * deltaPos[1];

}//CalcInfluenceFromFarCells(...)

/// Печать всего дерева в текстовый файл
void Cell::PrintTree()
{
	ofstream outfile;
	if (level == 0)
		outfile.open("tree.txt");
	else outfile.open("tree.txt", std::ios_base::app);

	//if (lowLevel)
	outfile << level << " " << leftDown[0] << " " << leftDown[1] << " " << rightUp[0] << " " << rightUp[1] << std::endl;

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

void Cell::PrintMyVortexes()
{
	ofstream outfile;
	outfile.open("vortexes.txt", std::ios_base::app);

	outfile << level << " " << itEnd - itBegin << std::endl;

	for (auto it = itBegin; it != itEnd; ++it)
		outfile << it->r()[0] << " " << it->r()[1] << " " << it->g() << std::endl;
	outfile.close();
}