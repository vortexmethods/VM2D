/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.6    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/10/28     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Sheet2D.cpp                                                      |
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
\brief Файл кода с описанием класса Sheet
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.6   
\date 28 октября 2019 г.
*/

#include "mpi.h"

#include "Sheet2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

Sheet::Sheet(const World2D& W_, int dim_)
	: dim(dim_), W(W_)
{ };


void Sheet::SetLayers(size_t np)
{
	freeVortexSheet_.resize(np * dim, 0.0);
	attachedVortexSheet_.resize(np * dim, 0.0);
	attachedSourceSheet_.resize(np * dim, 0.0);
};//SetLayers(...)


void Sheet::FreeSheetSynchronize() const
{
	//Синхронизировать размер не нужно, он и так синхронизирован
	MPI_Bcast(const_cast<double*>(freeVortexSheet_.data()), (int)freeVortexSheet_.size(), MPI_DOUBLE, 0, W.getParallel().commWork);
}
