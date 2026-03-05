/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: cuKnn.cuh                                                        |
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
\brief Заголовки интерфейса поиска ближайших соседей на CUDA
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\author Серебровская Екатерина Александровна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef CUKNN_CUH
#define CUKNN_CUH

#include <vector>
#include "Vortex2D.h"

static const int knbForRestruct = 5;

struct prDoubleInt
{
	double first;
	size_t second;
};

class VectorsForKnn
{
public:
	VectorsForKnn();
	~VectorsForKnn();

	//указатель на вихри
	Vortex2D* vtxPtr;

	//мортоновские коды вихрей
	int* mcdataPtr;
	int* mcdata_unsortedPtr;

	//порядок сортировки
	int* indexPtr;
	int* index_unsortedPtr;

	prDoubleInt* initdistPtr;       // указатели на структуру данных соседей = { расстояние; до какого элемента (в отсортированном массиве) }
	prDoubleInt* initdistPtrSdvig;  // то же, но после сдвига мортоновских кодов, потом мерджится initdistPtr

	void* sortBuffer;
	int sortBufferSizeInBytes;

	void ResizeAll(int k, int nvtx);

private:
	int reservedVtx;

};


template<int k>
double kNNcuda(Point2D minr, Point2D maxr, int blocks, const std::vector<Vortex2D>& vtx,
	std::vector<std::pair<double, size_t>>& initdist,
	VectorsForKnn*& vecForKnn,
	double cSP, double cRBP, double maxG, double epsCol, int typ);




#endif //KNN_CUH