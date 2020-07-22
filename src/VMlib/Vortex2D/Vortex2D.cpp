/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.9    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2020/07/22     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Vortex2D.cpp                                                     |
| Info: Source code of VMlib                                                  |
|                                                                             |
| This file is part of VMlib.                                                 |
| VMLib is free software: you can redistribute it and/or modify it            |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VMlib is distributed in the hope that it will be useful, but WITHOUT        |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VMlib.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Файл кода с описанием класса Vortex2D
\author Марчевский Илья Константинович
\version 1.9   
\date 22 июля 2020 г.
*/

#include <stddef.h>

#include "Vortex2D.h"

using namespace VMlib;

MPI_Datatype Vortex2D::mpiVortex2D;
size_t Vortex2D::offsPos; 
size_t Vortex2D::offsGam;


void Vortex2D::CreateMpiType()
{
	int          len[2] = { 2, 1 };
	MPI_Aint     pos[2] = { offsetof(Vortex2D, pos), offsetof(Vortex2D, gam) };
	MPI_Datatype typ[2] = { MPI_DOUBLE, MPI_DOUBLE };

	MPI_Type_create_struct(2, len, pos, typ, &mpiVortex2D);
	MPI_Type_commit(&mpiVortex2D);

	offsPos = offsetof(Vortex2D, pos);
	offsGam = offsetof(Vortex2D, gam);

	/*
	int          len[4] = { 1, 2, 1, 1 };
	MPI_Aint     pos[4] = { 0, offsetof(Vortex2D, pos), offsetof(Vortex2D, gam), sizeof(Vortex2D) };
	MPI_Datatype typ[4] = { MPI_LB, MPI_DOUBLE, MPI_DOUBLE, MPI_UB };

	MPI_Type_create_struct(4, len, pos, typ, &mpiVortex2D);
	MPI_Type_commit(&mpiVortex2D);

	offsPos = offsetof(Vortex2D, pos);
	offsGam = offsetof(Vortex2D, gam);
	*/
}//CreateMpiType()
