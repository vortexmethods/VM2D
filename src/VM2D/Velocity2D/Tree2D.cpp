/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.11   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2022/08/07     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Tree2D.cpp                                                       |
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
\brief Файл кода с описанием класса Tree
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.11
\date 07 августа 2022 г.
*/

#include "Tree2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Velocity2DBarnesHut.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

/*
/// Конструктор
Tree::Tree(const World2D& W_, const WakeDataBase& points, const PointType type, int numOfLevels_) : W(W_), numOfLevels(numOfLevels_), rootCell(*this)
{
	minDist = 2.0 * W.getPassport().wakeDiscretizationProperties.eps;
	allPnt.clear();

	switch (type)
	{
	case PointType::wake:
		CreatePointsCopy(W.getWake(), PointType::wake);
		break;

	case PointType::wakeVP:
		CreatePointsCopy(W.getMeasureVP().getWakeVP(), PointType::wakeVP);
		break;

	case PointType::sourceWake:
		CreatePointsCopy(W.getSource(), PointType::sourceWake);
		break;

	default:
		break;
	}
	
	rootCell.CreateRootCell(allPnt);
	rootCell.CreateAllCells(rootCell.itBegin, rootCell.itEnd);
};

/// Конструктор
Tree::Tree(const World2D& W_, const PointType type, int numOfLevels_) : W(W_), numOfLevels(numOfLevels_), rootCell(*this)
{
	minDist = 2.0 * W.getPassport().wakeDiscretizationProperties.eps;
	allPnt.clear();

	switch (type)
	{
	case PointType::sheetGam:
		CreateSheetsCopy(PointType::sheetGam);
		break;

	case PointType::source:
		CreateSheetsCopy(PointType::source);
		break;

	default:
		break;
	}

	rootCell.CreateRootCell(allPnt);
	rootCell.CreateAllCells(rootCell.itBegin, rootCell.itEnd);
};

/// Деструктор
Tree::~Tree()
{

};

// Создание массива указателей на массив точек (для сортировки), используется для массивов pointsCopyVP и pointsCopyWake
void Tree::CreatePointsCopy(const WakeDataBase& points, const PointType type)
{
	allPnt.reserve(points.vtx.size());

	for (int i = 0; i < (int)points.vtx.size(); ++i)
		allPnt.emplace_back(points.vtx[i], { -1, -1 }, i, type);
}//CreatePointsCopy(...)

// Создание массива указателей на массив точек (для сортировки), используется для массива sheetsCopy
void Tree::CreateSheetsCopy(const PointType type)
{
	size_t nSh = 0;
	for (size_t i = 0; i < W.getNumberOfBoundary(); ++i)
		nSh += W.getBoundary(i).sheets.getSheetSize();

	allPnt.reserve(nSh);
	int num;
	switch (type)
	{
	case PointType::sheetGam:
		num = 0;
		for (int i = 0; i < (int)(W.getNumberOfBoundary()); ++i)
			for (int j = 0; j < (int)(W.getBoundary(i).sheets.getSheetSize()); ++j)
			{
				allPnt.emplace_back(Vortex2D{ 0.5 * (W.getAirfoil(i).getR(j) + W.getAirfoil(i).getR(j + 1)), W.getBoundary(i).sheets.attachedVortexSheet(j, 0) * W.getAirfoil(i).len[j] }, { i, j }, num, type);
				num++;
			}
		break;
	case PointType::source:
		num = 0;
		for (int i = 0; i < (int)(W.getNumberOfBoundary()); ++i)
			for (int j = 0; j < (int)(W.getBoundary(i).sheets.getSheetSize()); ++j)
			{
				allPnt.emplace_back(Vortex2D{ 0.5 * (W.getAirfoil(i).getR(j) + W.getAirfoil(i).getR(j + 1)), W.getBoundary(i).sheets.attachedSourceSheet(j, 0) * W.getAirfoil(i).len[j] }, { i, j }, num, type);
				num++;
			}
		break;
	default:
		break;
	}
}//CreateSheetsCopy(...)
*/