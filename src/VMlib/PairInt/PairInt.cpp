/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.10   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2021/05/17     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: PairInt.cpp                                                      |
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
\brief Файл кода с описанием класса v3D
\author Марчевский Илья Константинович
\version 1.10
\date 17 мая 2021 г.
*/

#include "PairInt.h"

using namespace VMlib;

MPI_Datatype PairInt::mpiPairInt;


//Конструктор и приведение типа из numvector<int, 2>
PairInt::PairInt(const numvector<int, 2>& _r)
{
	data()[0] = _r[0];
	data()[1] = _r[1];
}//PairInt(...)


//Конструктор копирования
PairInt::PairInt(const PairInt& _r)
{
	data()[0] = _r[0];
	data()[1] = _r[1];
}//PairInt(...)


//Конструктор инициализации списком
PairInt::PairInt(const std::initializer_list<int>& z)
{
	for (size_t i = 0; i < 2; ++i)
		data()[i] = *(z.begin() + i);
}//PairInt(...)


// Cоздание MPI-описателя типа
void PairInt::CreateMpiType()
{
	int          len[1] = { 2 };
	MPI_Aint     pos[1] = { 0 };
	MPI_Datatype typ[1] = { MPI_INT };

	MPI_Type_create_struct(1, len, pos, typ, &mpiPairInt);
	MPI_Type_commit(&mpiPairInt);


	/*
	int          len[3] = { 1, 2, 1 };
	//MPI_Aint     pos[3] = { 0, offsetof(PairInt, data()), sizeof(PairInt) };
	MPI_Aint     pos[3] = { 0, 0, sizeof(PairInt) };
	MPI_Datatype typ[3] = { MPI_LB, MPI_INT, MPI_UB };

	MPI_Type_create_struct(3, len, pos, typ, &mpiPairInt);
	MPI_Type_commit(&mpiPairInt);
	*/
}//CreateMpiType()
