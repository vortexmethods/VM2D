/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.5    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/02/20     |
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
\version 1.5   
\date 20 февраля 2019 г.
*/

#include "Sheet2D.h"

using namespace VM2D;

//Установка pазмерностей всех векторов и их обнуление
void Sheet::SetLayersDim(size_t np, size_t layerDim)
{
	freeVortexSheet.resize(np);
	attachedVortexSheet.resize(np);
	attachedSourceSheet.resize(np);

	for (size_t j = 0; j < np; ++j)
	{
		freeVortexSheet[j].resize(layerDim, 0.0);
		attachedVortexSheet[j].resize(layerDim, 0.0);
		attachedSourceSheet[j].resize(layerDim, 0.0);
	}
};//SetLayersDim(...)