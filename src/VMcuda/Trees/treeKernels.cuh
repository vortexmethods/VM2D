/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: treeKernels.cuh                                                  |
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
\brief Заголовки функций-оберток для работы с деревом на CUDA
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 14 января 2024 г.
*/

#ifndef TREEKERNELS_CUH
#define TREEKERNELS_CUH

#include "Point2D.h"
#include "cudaTreeInfo.h"

#ifdef USE_CUDA

#define WARPSIZE 32

//было 28 - и иногда падало
#define MAXDEPTH 32

namespace BHcu
{
	void CudaTestError(const char* msg);
	float treeBoundingBoxWrapper(CudaTreeInfo& treeInfo);
	float treeMortonCodesWrapper(CudaTreeInfo& treeInfo, bool sort);
	float treeMortonInternalNodesWrapper(CudaTreeInfo& treeInfo, bool sort);

	float treeClearKernelWrapper(CudaTreeInfo& treeInfo, const int order);
	float treeSummarizationWrapper(CudaTreeInfo& treeInfo, const int order, bool calcAABB);
	float treeCalcAABBWrapper(CudaTreeInfo& treeInfo);

	float treePanelsToPointsCalculationWrapper(const CudaTreeInfo& treePanelsInfo, const CudaTreeInfo& controlTreeInfo,
											   int order, double itolsq2, Point2D* __restrict velD);

	float treeClosestPanelToPointsCalculationWrapper(const CudaTreeInfo& treePanelsInfo, const CudaTreeInfo& controlTreeInfo, \
													 int* __restrict closePnlD, \
													 bool findOnlyInside = false, double* pseudoNormals = nullptr);

	float treePanelsSegmentsIntersectionCalculationWrapper(const CudaTreeInfo& treePanelsInfo, const CudaTreeInfo& controlTreeInfo, \
														   int* __restrict closePnlD);

	float treeMatrToVecNoCalculationWrapper(CudaTreeInfo& treePanelsInfo, double itolsq);

	
	__global__
		void treeI1I2CalculationKernel(
			const int nnodesd, const int nbodiesd,
			const double minRd,
			const int2* __restrict Mchildd,
			const double2* __restrict momsd,
			const double4* __restrict vtxd,
			const int* __restrict MmortonCodesIdxd,
			const double2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
			const double4* __restrict Mlowerupperd,
			double* __restrict I1d,
			double2* __restrict I2d,
			const double* __restrict epsast);

	__global__
		void ClearI0I3(int nbodiesd, float* __restrict I0d, float2* __restrict I3d);

	__global__
		void treeI0I3CalculationKernel(
			const int nnodesd, const int nbodiesd,
			const double minRd,
			const int2* __restrict Mchildd,
			const double2* __restrict momsd,
			const double4* __restrict vtxd,
			const int* __restrict MmortonCodesIdxd,
			const double2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
			const double4* __restrict Mlowerupperd, int* range,
			float* __restrict I0d,
			float2* __restrict I3d,
			const double* __restrict epsast, const double* __restrict meanEpsd, const int npan, const double* __restrict pansd, double* __restrict visstrd);
}

#endif

#endif