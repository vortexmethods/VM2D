/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.7    |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2019/11/22     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Point2D.cpp                                                      |
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
\brief Файл кода с описанием класса Point2D
\author Марчевский Илья Константинович
\version 1.7   
\date 22 ноября 2019 г.
*/

#include <stddef.h>

#include "Point2D.h"

using namespace VMlib;

MPI_Datatype Point2D::mpiPoint2D;


//Конструктор и приведение типа из numvector<double, 2>
Point2D::Point2D(const numvector<double, 2>& _r)
{
	r[0] = _r[0];
	r[1] = _r[1];
}//Point2D(...)


//Конструктор копирования
Point2D::Point2D(const Point2D& _r)
{
	r[0] = _r[0];
	r[1] = _r[1];
}//Point2D(...)


//Конструктор инициализации списком
Point2D::Point2D(const std::initializer_list<double>& z)
{
	for (size_t i = 0; i < 2; ++i)
		r[i] = *(z.begin() + i);
}//Point2D(...)


//Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
Point2D Point2D::rotated(const double angle) const
{
	Point2D res;
	double cosa = cos(angle);
	double sina = sin(angle);

	res[0] = r[0] * cosa - r[1] * sina;
	res[1] = r[0] * sina + r[1] * cosa;
	return res;
}//rotated(...)


// Cоздание MPI-описателя типа
void Point2D::CreateMpiType()
{
	int          len[3] = { 1, 2, 1 };
	MPI_Aint     pos[3] = { 0, offsetof(Point2D, r), sizeof(Point2D) };
	MPI_Datatype typ[3] = { MPI_LB, MPI_DOUBLE, MPI_UB };

	MPI_Type_create_struct(3, len, pos, typ, &mpiPoint2D);
	MPI_Type_commit(&mpiPoint2D);
}//CreateMpiType()
