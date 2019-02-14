/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.5    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2019/02/20     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: v3D.cpp                                                          |
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
\version 1.5
\date 20 февраля 2019 г.
*/

#include "v3D.h"

using namespace VMlib;

MPI_Datatype v3D::mpiv3D;


//Конструктор и приведение типа из numvector<double, 3>
v3D::v3D(const numvector<double, 3>& _r)
{
	r[0] = _r[0];
	r[1] = _r[1];
	r[2] = _r[2];
}//v3D(...)


//Конструктор копирования
v3D::v3D(const v3D& _r)
{
	r[0] = _r[0];
	r[1] = _r[1];
	r[2] = _r[2];
}//v3D(...)


//Конструктор инициализации списком
v3D::v3D(const std::initializer_list<double>& z)
{
	for (size_t i = 0; i < 3; ++i)
		r[i] = *(z.begin() + i);
}//v3D(...)


//Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
v3D v3D::rotated(const double angle, const v3D& axis) const
{
	/// \todo Реализовать!
	/// \warning Пока возвращает себя!
	
	v3D res = *this;
	//double cosa = cos(angle);
	//double sina = sin(angle);
	//
	//res[0] = r[0] * cosa - r[1] * sina;
	//res[1] = r[0] * sina + r[1] * cosa;
	return res;
}//rotated(...)


// Cоздание MPI-описателя типа
void v3D::CreateMpiType()
{
	int          len[3] = { 1, 3, 1 };
	MPI_Aint     pos[3] = { 0, offsetof(v3D, r), sizeof(v3D) };
	MPI_Datatype typ[3] = { MPI_LB, MPI_DOUBLE, MPI_UB };

	MPI_Type_create_struct(3, len, pos, typ, &mpiv3D);
	MPI_Type_commit(&mpiv3D);
}//CreateMpiType()