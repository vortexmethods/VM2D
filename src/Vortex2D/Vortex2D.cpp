/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Vortex2D.cpp                                                     |
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
\brief Файл кода с описанием класса Vortex2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/

#include "Vortex2D.h"


MPI_Datatype Vortex2D::mpiVortex2D;
size_t Vortex2D::offsPos; 
size_t Vortex2D::offsGam;


void Vortex2D::CreateMpiType()
{
	int          len[4] = { 1, 2, 1, 1 };
	MPI_Aint     pos[4] = { 0, offsetof(Vortex2D, pos), offsetof(Vortex2D, gam), sizeof(Vortex2D) };
	MPI_Datatype typ[4] = { MPI_LB, MPI_DOUBLE, MPI_DOUBLE, MPI_UB };

	MPI_Type_create_struct(4, len, pos, typ, &mpiVortex2D);
	MPI_Type_commit(&mpiVortex2D);

	offsPos = offsetof(Vortex2D, pos);
	offsGam = offsetof(Vortex2D, gam);
}//CreateMpiType()