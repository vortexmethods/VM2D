/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
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
| (at your option) any later version.	                                      |
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
\version 1.1
\date 2 апреля 2018 г.
*/


#include "Boundary.h"



//Конструктор
Boundary::Boundary(const Passport& passport_, const Airfoil& afl_, const std::vector<std::unique_ptr<Boundary>>& allBoundary_, int sheetDim_, const Wake& wake_, const Parallel& parallel_, gpu& cuda_)
	: passport(passport_ ), afl(afl_), allBoundary(allBoundary_), sheetDim(sheetDim_), wake(wake_), parallel(parallel_), cuda(cuda_), CC(afl.r)
{
	sheets.SetLayersDim(afl.np, sheetDim);
	virtualWake.resize(0);
}//Boundary(...)


//MPI-синхронизация вихревого следа
void Boundary::VirtualWakeSynchronize()
{
	int nV;
	if (parallel.myidWork == 0)
		nV = (int)virtualWake.size();
	MPI_Bcast(&nV, 1, MPI_INT, 0, parallel.commWork);

	if (parallel.myidWork > 0)
		virtualWake.resize(nV);

	MPI_Bcast(virtualWake.data(), nV, Vortex2D::mpiVortex2D, 0, parallel.commWork);
}//VirtualWakeSinchronize()