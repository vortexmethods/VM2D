/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Boundary.cpp                                                     |
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
\brief Файл кода с описанием класса Boundary
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#include "Boundary.h"
#include "World2D.h"

//Конструктор
Boundary::Boundary(const World2D& W_, size_t numberInPassport_, int sheetDim_) : 
	W(W_ ), 
	numberInPassport(numberInPassport_),
	sheetDim(sheetDim_), 
	afl(W_.getAirfoil(numberInPassport_)),
	CC(W_.getAirfoil(numberInPassport_).r),
	virtualWake(W_, *this)
{
	sheets.SetLayersDim(afl.np, sheetDim);
	oldSheets.SetLayersDim(afl.np, sheetDim);
	virtualWake.vtx.resize(0);
	virtualWake.vecHalfGamma.resize(0);
}//Boundary(...)


//MPI-синхронизация вихревого следа
void Boundary::VirtualWakeSynchronize()
{
	int nV;
	if (W.getParallel().myidWork == 0)
		nV = static_cast<int>(virtualWake.vtx.size());
	MPI_Bcast(&nV, 1, MPI_INT, 0, W.getParallel().commWork);

	if (W.getParallel().myidWork > 0)
		virtualWake.vtx.resize(nV);

	MPI_Bcast(virtualWake.vtx.data(), nV, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);
}//VirtualWakeSinchronize()