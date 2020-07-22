/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.9    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/07/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: PointsCopy2D.cpp                                                 |
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
\brief Файл кода с описанием класса PointsCopy
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.9   
\date 22 июля 2020 г.
*/



#include "PointsCopy2D.h"



void PointsCopy::reserve(size_t size)
{
	vtx.reserve(size);
	aflPnl.reserve(size);
	num.reserve(size);
	velo.reserve(size);
	domainRadius.reserve(size);
}//reserve(...)

void PointsCopy::emplace_back(Vortex2D _vtx, std::pair<int, int> _aflPnl, int _num, PointType _type, Point2D _velo, double _domainRadius)
{
	vtx.emplace_back(_vtx);
	aflPnl.emplace_back(_aflPnl);
	num.emplace_back(_num);
	type = _type;
	velo.emplace_back(_velo);
	domainRadius.emplace_back(_domainRadius);
}//emplace_back(...)

void PointsCopy::push_back(Vortex2D _vtx, std::pair<int, int> _aflPnl, int _num)
{
	vtx.push_back(_vtx);
	aflPnl.push_back(_aflPnl);
	num.push_back(_num);
}//push_back(...)

void PointsCopy::clear()
{
	vtx.clear();
	aflPnl.clear();
	num.clear();
	velo.clear();
	domainRadius.clear();
}//clear()

size_t PointsCopy::size() const
{
	return vtx.size();
}//size()


PointsCopy::iterator PointsCopy::begin()
{
	return { vtx.begin(), aflPnl.begin(), num.begin() };
}//begin()

PointsCopy::iterator PointsCopy::end()
{
	return { vtx.end(), aflPnl.end(), num.end() };
}//end()


//////////////////////////////////////

PointsCopy_iterator::all_reference& PointsCopy_iterator::all_reference::operator= (all_copy const& p)
{
	vtx = p.vtx;
	aflPnl = p.aflPnl;
	num = p.num;

	return *this;
}

PointsCopy_iterator::all_reference& PointsCopy_iterator::all_reference::operator= (all_copy&& p)
{
	vtx = std::move(p.vtx);
	aflPnl = std::move(p.aflPnl);
	num = std::move(p.num);

	return *this;
}

bool PointsCopy_iterator::LessXrr(all_reference const& lhs, all_reference const& rhs)
{
	return lhs.vtx.r()[0] < rhs.vtx.r()[0];
}

bool PointsCopy_iterator::LessYrr(all_reference const& lhs, all_reference const& rhs)
{
	return lhs.vtx.r()[1] < rhs.vtx.r()[1];
}

bool PointsCopy_iterator::LessXcr(all_copy const& lhs, all_reference const& rhs)
{
	return lhs.vtx.r()[0] < rhs.vtx.r()[0];
}

bool PointsCopy_iterator::LessYcr(all_copy const& lhs, all_reference const& rhs)
{
	return lhs.vtx.r()[1] < rhs.vtx.r()[1];
}

bool PointsCopy_iterator::LessXrc(all_reference const& lhs, all_copy const& rhs)
{
	return lhs.vtx.r()[0] < rhs.vtx.r()[0];
}

bool PointsCopy_iterator::LessYrc(all_reference const& lhs, all_copy const& rhs)
{
	return lhs.vtx.r()[1] < rhs.vtx.r()[1];
}

bool PointsCopy_iterator::LessXcc(all_copy const& lhs, all_copy const& rhs)
{
	return lhs.vtx.r()[0] < rhs.vtx.r()[0];
}

bool PointsCopy_iterator::LessYcc(all_copy const& lhs, all_copy const& rhs)
{
	return lhs.vtx.r()[1] < rhs.vtx.r()[1];
}
