/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.2    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/06/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: VelocityFourier.h                                                |
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
\brief Заголовочный файл с описанием класса WakeDataBase
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.2
\date 14 июня 2018 г.
*/

#ifndef WAKEDATABASE_H
#define WAKEDATABASE_H

#include <vector>

#include "gpudefs.h"
#include "Point2D.h"
#include "Vortex2D.h"


class WakeDataBase
{
public:
	/// Список вихревых элементов
	std::vector<Vortex2D> vtx;

#if defined(USE_CUDA)
	mutable double* devWakePtr;
	mutable size_t devNWake;

	mutable std::vector<Point2D> tmpVels;
	mutable double* devVelsPtr;

	mutable std::vector<double> tmpRads;
	mutable double* devRadsPtr;

	mutable std::vector<double> tmpI0;
	mutable double* devI0Ptr;

	mutable std::vector<double> tmpI1;
	mutable double* devI1Ptr;

	mutable std::vector<Point2D> tmpI2;
	mutable double* devI2Ptr;

	mutable std::vector<Point2D> tmpI3;
	mutable double* devI3Ptr;

	//Данные для коллапса
	mutable std::vector<int> tmpNei;
	mutable int* devMeshPtr;
	mutable int* devNeiPtr;

	mutable std::vector<numvector<int, 2>> mesh;
#endif

	WakeDataBase() {};
	virtual ~WakeDataBase() {}; 
};

#endif

