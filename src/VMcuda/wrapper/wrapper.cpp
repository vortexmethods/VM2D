/*--------------------------------*- BHgpu -*----------------*---------------*\
| #####   ##  ##                |                            | Version 1.5    |
| ##  ##  ##  ##   ####  ##  ## |  BHgpu: Barnes-Hut method  | 2023/08/29     |
| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*
| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |
| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2023 I. Marchevsky, E. Ryatina, A. Kolganova             |
| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |
*-----------------------------------------------------------------------------*
| File name: main.cpp                                                         |
| Info: Source code of BHgpu                                                  |
|                                                                             |
| This file is part of BHgpu.                                                 |
| BHcu is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| BHcu is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with BHgpu.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/

/*
 * Portions of this program were originally released under the following license
 *
 * CUDA BarnesHut v3.1: Simulation of the gravitational forces
 * in a galactic cluster using the Barnes-Hut n-body algorithm
 *
 * Copyright (c) 2013, Texas State University-San Marcos. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted for academic, research, experimental, or personal use provided that
 * the following conditions are met:
 *
 *    * Redistributions of source code must retain the above copyright notice,
 *      this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright notice,
 *      this list of conditions and the following disclaimer in the documentation
 *      and/or other materials provided with the distribution.
 *    * Neither the name of Texas State University-San Marcos nor the names of its
 *      contributors may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 *
 * For all other uses, please contact the Office for Commercialization and Industry
 * Relations at Texas State University-San Marcos <http://www.txstate.edu/ocir/>.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Martin Burtscher <burtscher@txstate.edu>
 *
 */

 /*!
 \file
 \brief Barnes-Hut method (CUDA) for 2D vortex particles + Morton tree
 \author Марчевский Илья Константинович
 \author Рятина Евгения Павловна
 \author Колганова Александра Олеговна
 \version 1.5
 \date 29 августа 2023 г.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include "omp.h"

#include "types.cuh"
#include "cuKernels.cuh"
#include "cuLib2D.cuh"
#include "wrapper.h"


namespace BHcu
{

	const real IDPI = (real)0.15915494309189534;

	/******************************************************************************/
	/******************************************************************************/

	struct CudaCalcGab
	{
		realPoint* maxpt;
		realPoint* minpt;
		int blocks;

		CudaCalcGab() //(realPoint*& maxpt, realPoint*& minpt, int blocks)
			//: maxpt_(maxpt), minpt_(minpt), blocks_(blocks) 
		{
			CudaSelect(0);
			setBlocks(blocks);
			
			maxpt = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));
			minpt = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));
		};

		float calc(int npoints, const realVortex* pointsl)
		{
			float time;
			time = McuBoundingBoxKernelFree(nullptr, maxpt, minpt, npoints, (void*)pointsl, sizeof(realVortex), (int)realVortex::offsPos);
			return time;
		}

		~CudaCalcGab()
		{
			cudaDelete(maxpt, 1);
			cudaDelete(minpt, 2);
		}
	};

	struct CudaSorter
	{
		int* pointsMortonCodesKeyl;
		int* pointsMortonCodesKeyUnsortl;

		int* pointsMortonCodesIdxl;
		int* pointsMortonCodesIdxUnsortl;

		int npoints_;
		CudaCalcGab gab;
		const realVortex* pointsl_;

		CudaSorter(int npoints, const realVortex* pointsl)
			: npoints_(npoints), pointsl_(pointsl)
		{
			pointsMortonCodesKeyl = (int*)cudaNew(npoints, sizeof(int));
			pointsMortonCodesKeyUnsortl = (int*)cudaNew(npoints, sizeof(int));

			pointsMortonCodesIdxl = (int*)cudaNew(npoints, sizeof(int));
			pointsMortonCodesIdxUnsortl = (int*)cudaNew(npoints, sizeof(int));
		};

		float calc()
		{
			float timeGab, timeCodes;
			timeGab = gab.calc(npoints_, pointsl_);
			
			timeCodes = McuMortonCodesKernelFree(gab.maxpt, gab.minpt,
				pointsMortonCodesKeyUnsortl, pointsMortonCodesIdxUnsortl,
				pointsMortonCodesKeyl, pointsMortonCodesIdxl, nullptr,
				npoints_, pointsl_, (int)sizeof(realVortex), (int)Vortex2D::offsPos, true);
			return timeGab + timeCodes;
		}

		~CudaSorter()
		{
			cudaDelete(pointsMortonCodesKeyl, 3);
			cudaDelete(pointsMortonCodesKeyUnsortl, 4);
			cudaDelete(pointsMortonCodesIdxl, 5);
			cudaDelete(pointsMortonCodesIdxUnsortl, 6);
		}
	};


	void rebuildBaseTree(CUDApointers& ptrs, const int nbodies, /*const realVortex* vtxl*/ const void* ptrl, int sizeOfElement, int offsetOfPointInElement, int objType, int nnodes, int order, double* timing, 
		bool bodyZeroSize, const double* XYXY)
	{
		timing[0] += cuInitializationKernel();
		
		double tBBK	= McuBoundingBoxKernel(ptrs, nbodies, ptrl, sizeOfElement, offsetOfPointInElement);
		timing[1] += tBBK;
		
		double tTBK = 0;
		tTBK += McuMortonCodesKernel(ptrs, nbodies, ptrl, sizeOfElement, offsetOfPointInElement /*(int)sizeof(realVortex), (int)Vortex2D::offsPos*/);
		
		//std::cout << nbodies << std::endl;
		
		tTBK += McuMortonInternalNodesKernel(ptrs, nbodies);
		tTBK += McuMortonInternalCellsGeometryKernel(ptrs, nbodies, nnodes);
		timing[2] += tTBK;
		
		double tCLK = cuClearKernel2(ptrs, order, nnodes, nbodies);
		timing[3] += tCLK;
		
		double tSKK = 0;
		tSKK += cuAABBKernel2(ptrs, nnodes, nbodies, ptrl, sizeOfElement, offsetOfPointInElement, bodyZeroSize, XYXY);
		tSKK += cuClearKernel2(ptrs, order, nnodes, nbodies);
		tSKK += cuSummarizationKernel2(ptrs, order, nnodes, nbodies, (double*)ptrl, objType);
		
		timing[4] += tSKK;


		//std::cout << " BBK = " << tBBK << " TBK = " << tTBK << " CLK = " << tCLK << " SKK = " << tSKK << std::endl;
	}


	double memoryAllocate(CUDApointers& ptrs, int nnodes, int nbodies, int nbodiesOld, int blocks, int order)
	{
		double starttime, endtime;
		starttime = omp_get_wtime();

		if (nbodiesOld > 0)
		{
			//std::cout << "BHgpu: free CUDA-memory" << std::endl;
			cudaDelete(ptrs.massl, 7);

			cudaDelete(ptrs.momsl, 8);
			cudaDelete(ptrs.El, 9);

			cudaDelete(ptrs.maxrl, 10);
			cudaDelete(ptrs.minrl, 11);

			///For Morton tree
			cudaDelete(ptrs.MmortonCodesKeyUnsortl, 12);
			cudaDelete(ptrs.MmortonCodesIdxUnsortl, 13);
			cudaDelete(ptrs.MmortonCodesKeyl, 14);
			cudaDelete(ptrs.MmortonCodesIdxl, 15);

			cudaDelete(ptrs.Mposl, 16);
			cudaDelete(ptrs.Mlowerl, 17);
			cudaDelete(ptrs.Mupperl, 18);
			cudaDelete(ptrs.Mparentl, 19);
			cudaDelete(ptrs.Mchildl, 20);
			cudaDelete(ptrs.Mrangel, 21);

			cudaDelete(ptrs.MlevelUnsortl, 22);
			cudaDelete(ptrs.MlevelSortl, 23);
			cudaDelete(ptrs.MindexUnsortl, 24);
			cudaDelete(ptrs.MindexSortl, 25);
			cudaDelete(ptrs.MindexSortTl, 26);
		}

		//std::cout << "BHgpu: allocation GPU-memory: nbodies = " << nbodies << ", nnodes = " << nnodes << ", order = " << order << std::endl;

		//unsigned long long int mem = 0;
		ptrs.massl = (int*)cudaNew(nbodies - 1, sizeof(int));
		//mem += (nbodies - 1) * sizeof(int);

		ptrs.momsl = (realPoint*)cudaNew((nbodies - 1) * order, sizeof(realPoint));

		//printf("ALLOCATED for MOMS = %d bytes for %d bodies, order = %d, sizeof = %d\n", int((nbodies - 1) * order * sizeof(realPoint)), nbodies - 1, order, sizeof(realPoint));

		ptrs.El = nullptr;
		//mem += (nbodies - 1) * order * sizeof(realPoint);

		ptrs.maxrl = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));
		ptrs.minrl = (realPoint*)cudaNew(blocks * FACTOR1, sizeof(realPoint));
		//mem += 2 * blocks * FACTOR1 * sizeof(realPoint);

		///For MortonTree
		ptrs.MmortonCodesKeyUnsortl = (int*)cudaNew(nbodies, sizeof(int));
		ptrs.MmortonCodesKeyl = (int*)cudaNew(nbodies, sizeof(int));
		ptrs.MmortonCodesIdxUnsortl = (int*)cudaNew(nbodies, sizeof(int));
		ptrs.MmortonCodesIdxl = (int*)cudaNew(nbodies, sizeof(int));
		//mem += 4 * nbodies * sizeof(int);

		ptrs.Mposl = (realPoint*)cudaNew(nbodies - 1, sizeof(realPoint));
		ptrs.Mlowerl = (realPoint*)cudaNew(nbodies - 1, sizeof(realPoint));
		ptrs.Mupperl = (realPoint*)cudaNew(nbodies - 1, sizeof(realPoint));
		//mem += 3 * (nbodies - 1) * sizeof(realPoint);

		ptrs.Mparentl = (int*)cudaNew(nnodes, sizeof(int));
		//mem += nnodes * sizeof(int);

		ptrs.Mchildl = (intPair*)cudaNew(nbodies - 1, sizeof(intPair));
		//mem += (nbodies - 1) * sizeof(intPair);

		ptrs.Mrangel = (intPair*)cudaNew(nnodes, sizeof(intPair)); //Нужно ли для всех?
		//std::cout << "Mrangel.size = " << nnodes << ", " << nnodes * sizeof(intPair) << " bytes" << std::endl;
		//mem += nnodes * sizeof(intPair); 

		ptrs.MlevelUnsortl = (int*)cudaNew(nbodies - 1, sizeof(int));
		ptrs.MlevelSortl = (int*)cudaNew(nbodies - 1, sizeof(int));
		ptrs.MindexUnsortl = (int*)cudaNew(nbodies - 1, sizeof(int));
		ptrs.MindexSortl = (int*)cudaNew(nbodies - 1, sizeof(int));
		ptrs.MindexSortTl = (int*)cudaNew(nbodies - 1, sizeof(int));
		//mem += 5 * (nbodies - 1) * sizeof(int);

		endtime = omp_get_wtime();
		return endtime - starttime;
	}



	// wake to points
	double wrapperInfluenceToPoints(
		const realVortex* vtxl, const realVortex* pointsl, int sizeOfArrayElement, int offsetOfPointInElement,		
		realPoint* vell, real* epsastl,
		CUDApointers& ptrs, bool rebuild, int nbodies, int npoints, double* timing, real eps, real theta,
		size_t& nbodiesOld, int nbodiesUp, int order,
		size_t nAfls, size_t* nVtxs, double** ptrVtxs, bool calcVelo, bool calcRadius)
	{
		double starttime, endtime;
		starttime = omp_get_wtime();

		//Число мультипроцессоров, заполняется функцией  setBlocks(blocks)
		int blocks;

		//Число ячеек дерева и тел
		int nnodes, nnodesUp;

		//Радиус вихря и параметр близости и их квадраты
		real epssq = (real)(eps * eps);
		real itolsq = (real)(1 / (theta * theta));

		CudaSelect(0);
		setBlocks(blocks); //"достает" число блоков, равное числу мультипроцессоров (blocks - по ссылке)		

		nnodes = nbodies * 2;
		if (nnodes < 1024 * blocks)
			nnodes = 1024 * blocks;
		while ((nnodes & (32 - 1)) != 0)  // 32 - это размер варпа
			nnodes++;
		nnodes--;

		if (rebuild)
		{
			nnodesUp = nbodiesUp * 2;
			if (nnodesUp < 1024 * blocks)
				nnodesUp = 1024 * blocks;
			while ((nnodesUp & (32 - 1)) != 0)  // 32 - это размер варпа
				nnodesUp++;
			nnodesUp--;
		}

		KernelsOptimization();


		for (int i = 0; i < 6; i++)
			timing[i] = 0;

		if (rebuild)
		{
			if (nbodiesUp > nbodiesOld)
				timing[1] += memoryAllocate(ptrs, nnodesUp, nbodiesUp, (int)nbodiesOld, blocks, order);

			nbodiesOld = nbodiesUp;			
			rebuildBaseTree(ptrs, nbodies, vtxl, (int)sizeof(realVortex), (int)Vortex2D::offsPos, 0, nnodes, order, timing, true, nullptr);
		}


		bool calcEpsAst = calcRadius;

		if (pointsl != vtxl)
		{
			TMortonCodesCalculator mCodes(npoints, (void*)pointsl, sizeOfArrayElement, offsetOfPointInElement);
			timing[5] += mCodes.timeOfWork;			
			timing[5] += cuForceCalculationKernel2points(ptrs, order, nnodes, nbodies, itolsq, epssq, vtxl, mCodes.getSortedIndices(), npoints, pointsl, vell, calcEpsAst, epsastl, nAfls, nVtxs, ptrVtxs);
		}
		else
			timing[5] += cuForceCalculationKernel2points(ptrs, order, nnodes, nbodies, itolsq, epssq, vtxl, ptrs.MmortonCodesIdxl, npoints, pointsl, vell, calcEpsAst, epsastl, nAfls, nVtxs, ptrVtxs);

		timing[6] = timing[1] + timing[2] + timing[3] + timing[4] + timing[5];

		endtime = omp_get_wtime();
		return endtime - starttime;
	}


	double wrapperInfluenceFromPanelsToPoints(
		double* dev_ptr_r,  //начала и концы панелей
		double* dev_ptr_freeVortexSheet, double* dev_ptr_freeVortexSheetLin,
		double* dev_ptr_attachedVortexSheet, double* dev_ptr_attachedVortexSheetLin,
		double* dev_ptr_attachedSourceSheet, double* dev_ptr_attachedSourceSheetLin,
		Vortex2D* dev_ptr_pt, //вихри в следе
		double*& dev_ptr_vel,  //куда сохранить результат 
		CUDApointers& ptrsi,  //указатель на дерево на i-м профиле
		bool rebuild,         //признак перестроения дерева вихрей
		int npt,			  //число вихрей в следе
		int npnli,			  //общее число панелей на профиле
		double* timings,      //засечки времени
		double eps,           //eps
		double theta,         //theta
		int order,            //order
		int schemeType)
	{		
		double starttime, endtime;
		starttime = omp_get_wtime();
		
		///*
		//Число мультипроцессоров, заполняется функцией  setBlocks(blocks)
		int blocks;

		//Число ячеек дерева и тел
		int nnodes, nnodesUp;

		//Радиус вихря и параметр близости и их квадраты
		real epssq = (real)(eps * eps);
		real itolsq = (real)(1 / (theta * theta));

		CudaSelect(0);
		setBlocks(blocks); //"достает" число блоков, равное числу мультипроцессоров (blocks - по ссылке)		

		//nnodes = nbodies * 2;
		nnodes = npnli * 2;
		if (nnodes < 1024 * blocks)
			nnodes = 1024 * blocks;
		while ((nnodes & (32 - 1)) != 0)  // 32 - это размер варпа
			nnodes++;
		nnodes--;

		if (rebuild)
		{
			//nnodesUp = nbodiesUp * 2;
			nnodesUp = npnli * 2;

			if (nnodesUp < 1024 * blocks)
				nnodesUp = 1024 * blocks;
			while ((nnodesUp & (32 - 1)) != 0)  // 32 - это размер варпа
				nnodesUp++;
			nnodesUp--;
		}

		KernelsOptimization();

		for (int i = 0; i < 6; i++)
			timings[i] = 0;

		double* panelPoints = (double*)cudaNew(12 * npnli, sizeof(double));//Для кусочно-постоянной схемы

		McuVerticesAndSheetsToPanelPoints(npnli, (double*)dev_ptr_r, 
			(double*)dev_ptr_freeVortexSheet, (double*)dev_ptr_freeVortexSheetLin,
			(double*)dev_ptr_attachedVortexSheet, (double*)dev_ptr_attachedVortexSheetLin,
			(double*)dev_ptr_attachedSourceSheet, (double*)dev_ptr_attachedSourceSheetLin,
			(double*)panelPoints, schemeType);

		if (rebuild)
		{
			//if (nbodiesUp > nbodiesOld)
			//timing[1] += memoryAllocate(ptrs, nnodesUp, nbodiesUp, (int)nbodiesOld, blocks, order);
			timings[1] += memoryAllocate(ptrsi, nnodesUp, npnli, npnli, blocks, order);

			//nbodiesOld = nbodiesUp;
			//rebuildBaseTree(ptrs, nbodies, vtxl, (int)sizeof(realVortex), (int)Vortex2D::offsPos, nnodes, order, timing);

			rebuildBaseTree(ptrsi, npnli, panelPoints, 12*(int)sizeof(double), 0, schemeType, nnodes, order, timings, false, dev_ptr_r);
		}

		if (npt > 0)
		{
			TMortonCodesCalculator mCodes(npt, (void*)dev_ptr_pt, sizeof(Vortex2D), 0);
			timings[5] += mCodes.timeOfWork;

			//timing[5] += cuForceCalculationKernel2points(ptrs, order, nnodes, nbodies, itolsq, epssq, vtxl, ptrs.MmortonCodesIdxl, npoints, pointsl, vell, true, epsastl, nAfls, nVtxs, ptrVtxs);
			timings[5] += cuForceCalculationKernelFromPanels2points(ptrsi, order, nnodes, npnli, itolsq, epssq, (double*)dev_ptr_r, panelPoints, 12 * (int)sizeof(double), 0, mCodes.getSortedIndices(), npt, dev_ptr_pt, (realPoint*)dev_ptr_vel, schemeType);
		}

			cudaDelete(panelPoints, 28);
						

		timings[6] = timings[1] + timings[2] + timings[3] + timings[4] + timings[5];
		//*/

		endtime = omp_get_wtime();
		return endtime - starttime;
	}



	double wrapperInfluenceToRHS(
		const realVortex* dev_ptr_vt,  //вихри в следе
		const double* dev_ptr_pt,      //начала и концы панелей
		double* dev_ptr_rhs,           //куда сохранить результат (для T0 и верхней половины T1)
		double* dev_ptr_rhslin,        //куда сохранить результат (для нижней половины T1)

		CUDApointers& ptrs,     //указатели на дерево вихрей
		bool rebuild,           //признак перестроения дерева вихрей

		int nvt,               //число вихрей в следе
		int nTotPan,           //общее число панелей на всех профилях
		double* timingsToRHS,  //засечки времени
		double theta,          //theta
		size_t& nbodiesOld, int nbodiesUp, int order, int scheme)
	{
		double starttime, endtime;
		starttime = omp_get_wtime();
		
		//Число мультипроцессоров, заполняется функцией  setBlocks(blocks)
		int blocks;

		//Число ячеек дерева и тел
		int nnodes, nnodesUp;

		KernelsOptimization();

		for (int i = 0; i < 6; i++)
			timingsToRHS[i] = 0;

		//Радиус вихря и параметр близости и их квадраты
		real itolsq = (real)(1 / (theta * theta));

		CudaSelect(0);
		setBlocks(blocks); //"достает" число блоков, равное числу мультипроцессоров (blocks - по ссылке)		

		nnodes = nvt * 2;
		if (nnodes < 1024 * blocks)
			nnodes = 1024 * blocks;
		while ((nnodes & (32 - 1)) != 0)  // 32 - это размер варпа
			nnodes++;
		nnodes--;

		if (rebuild)
		{
			nnodesUp = nbodiesUp * 2;
			if (nnodesUp < 1024 * blocks)
				nnodesUp = 1024 * blocks;
			while ((nnodesUp & (32 - 1)) != 0)  // 32 - это размер варпа
				nnodesUp++;
			nnodesUp--;
		}

		if (rebuild)
		{
			if (nbodiesUp > nbodiesOld)
				timingsToRHS[1] += memoryAllocate(ptrs, nnodesUp, nbodiesUp, (int)nbodiesOld, blocks, order);

			nbodiesOld = nbodiesUp;
			rebuildBaseTree(ptrs, nvt, dev_ptr_vt, (int)sizeof(realVortex), (int)Vortex2D::offsPos, 0, nnodes, order, timingsToRHS, true, nullptr);

			//std::cout << " BBK = " << timingsToRHS[1] << " TBK = " << timingsToRHS[2] << " CLK = " << timingsToRHS[3] << " SKK = " << timingsToRHS[4] << std::endl;
		}

		Point2D* controlPoints = (Point2D*)cudaNew(nTotPan, sizeof(Point2D));
		realPoint* El = (realPoint*)cudaNew(nTotPan * order, sizeof(realPoint));

		McuVerticesToControlPoints(nTotPan, (double*)dev_ptr_pt, (double*)controlPoints);

		//CudaSorter srt(nTotPan, controlPoints);
		//timingsToRHS[5] += srt.calc();

		TMortonCodesCalculator mCodes(nTotPan, (void*)controlPoints, sizeof(Point2D), 0);
		timingsToRHS[5] += mCodes.timeOfWork;

		double* ptrToLin = nullptr;
		if (scheme == 2)
			ptrToLin = dev_ptr_rhslin;

		timingsToRHS[5] += cuRhsCalculationKernel(ptrs, order, nnodes, nvt, itolsq, dev_ptr_vt, 
			mCodes.getSortedIndices(), El,
			nTotPan, dev_ptr_pt, (const real*)controlPoints, dev_ptr_rhs, ptrToLin);

		cudaDelete(El, 27);	
		cudaDelete(controlPoints, 28);

		//std::cout << " RhsCK = " << timingsToRHS[5] << std::endl;


		timingsToRHS[6] = timingsToRHS[1] + timingsToRHS[2] + timingsToRHS[3] + timingsToRHS[4] + timingsToRHS[5];

		endtime = omp_get_wtime();

		return endtime - starttime;
		
	}

	wrapperMatrixToVector::wrapperMatrixToVector(
		const double* dev_ptr_pt_,      //начала и концы панелей
		double* dev_ptr_rhs_,           //куда сохранить результат (для T0 и верхней половины T1)
		double* dev_ptr_rhslin_,        //куда сохранить результат (для нижней половины T1)

		CUDApointers& ptrs_,     //указатели на дерево вихрей
		bool rebuild_,           //признак перестроения дерева вихрей

		int nTotPan_,           //общее число панелей на всех профилях
		double* timingsMatr_,  //засечки времени
		double theta_,          //theta
		int order_, int scheme_)
		:
		dev_ptr_pt(dev_ptr_pt_), dev_ptr_rhs(dev_ptr_rhs_),
		dev_ptr_rhslin(dev_ptr_rhslin_), ptrs(ptrs_), rebuild(rebuild_),
		nTotPan(nTotPan_), timingsMatr(timingsMatr_), theta(theta_), order(order_), scheme(scheme_)
	{
		panelPoints = (double*)cudaNew(12 * nTotPan, sizeof(double));//Для кусочно-постоянной схемы
		controlPoints = (Point2D*)cudaNew(nTotPan, sizeof(Point2D));
		El = (realPoint*)cudaNew(nTotPan * order, sizeof(realPoint));

		McuVerticesToControlPoints(nTotPan, (double*)dev_ptr_pt, (double*)controlPoints);
		mCodesPtr = new TMortonCodesCalculator(nTotPan, (void*)controlPoints, sizeof(Point2D), 0);
		timingsMatr[5] += mCodesPtr->timeOfWork;
	};

	wrapperMatrixToVector::~wrapperMatrixToVector()
	{
		delete mCodesPtr;
		cudaDelete(panelPoints, 25);
		cudaDelete(El, 27);
		cudaDelete(controlPoints, 28);
	}


		double wrapperMatrixToVector::calculate(
			double* dev_ptr_freeVortexSheet, 
			double* dev_ptr_freeVortexSheetLin)
		{
			double starttime, endtime;
			starttime = omp_get_wtime();

			//Число мультипроцессоров, заполняется функцией  setBlocks(blocks)
			int blocks;

			//Число ячеек дерева и тел
			int nnodes, nnodesUp;

			KernelsOptimization();

			for (int i = 0; i < 6; i++)
				timingsMatr[i] = 0;

			//Радиус вихря и параметр близости и их квадраты
			real itolsq = (real)(1 / (theta * theta));

			CudaSelect(0);
			setBlocks(blocks); //"достает" число блоков, равное числу мультипроцессоров (blocks - по ссылке)		

			nnodes = nTotPan * 2;
			if (nnodes < 1024 * blocks)
				nnodes = 1024 * blocks;
			while ((nnodes & (32 - 1)) != 0)  // 32 - это размер варпа
				nnodes++;
			nnodes--;

			if (rebuild)
			{
				nnodesUp = nTotPan * 2;
				if (nnodesUp < 1024 * blocks)
					nnodesUp = 1024 * blocks;
				while ((nnodesUp & (32 - 1)) != 0)  // 32 - это размер варпа
					nnodesUp++;
				nnodesUp--;
			}

			//realVortex* controlPointsVortexes = (realVortex*)cudaNew(nTotPan, sizeof(realVortex));
			//McuVerticesToControlPointsVortexes(nTotPan, (double*)dev_ptr_pt, (double*)controlPointsVortexes, dev_gam);


			
			McuVerticesAndSheetsToPanelPoints(nTotPan, (double*)dev_ptr_pt,
				(double*)dev_ptr_freeVortexSheet, (double*)dev_ptr_freeVortexSheetLin,
				nullptr, nullptr,
				nullptr, nullptr,
				(double*)panelPoints, scheme);



			if (rebuild)
			{
				//if (nbodiesUp > nbodiesOld)
				timingsMatr[1] += memoryAllocate(ptrs, nnodesUp, nTotPan, nTotPan, blocks, order);

				//nbodiesOld = nbodiesUp;
				//rebuildBaseTree(ptrs, nTotPan, controlPointsVortexes, (int)sizeof(realVortex), (int)Vortex2D::offsPos, 0, nnodes, order, timingsToRHS, true, nullptr);
				rebuildBaseTree(ptrs, nTotPan, panelPoints, 12*(int)sizeof(double), 0, scheme, nnodes, order, timingsMatr, false, dev_ptr_pt);
				//std::cout << " BBK = " << timingsToRHS[1] << " TBK = " << timingsToRHS[2] << " CLK = " << timingsToRHS[3] << " SKK = " << timingsToRHS[4] << std::endl;
			}

			


			double* ptrToLin = nullptr;
			if (scheme == 2)
				ptrToLin = dev_ptr_rhslin;


			timingsMatr[5] += cuMatrixMulVectorCalculationKernel(ptrs, order, nnodes, nTotPan, itolsq, dev_ptr_pt, 
				dev_ptr_freeVortexSheet,
				dev_ptr_freeVortexSheetLin,
				mCodesPtr->getSortedIndices(), El,
				nTotPan, dev_ptr_pt, (const real*)controlPoints, dev_ptr_rhs, ptrToLin);


			//cudaDelete(controlPointsVortexes, 29);

			//std::cout << " RhsCK = " << timingsToRHS[5] << std::endl;


			timingsMatr[6] = timingsMatr[1] + timingsMatr[2] + timingsMatr[3] + timingsMatr[4] + timingsMatr[5];

			endtime = omp_get_wtime();

			return endtime - starttime;

		}





	
	double wrapperDiffusiveVeloI1I2(const realVortex* vtxl,
		real* i1l,
		realPoint* i2l,
		real* epsastl,
		CUDApointers& ptrs,
		bool rebuild,
		int nbodies,
		double* timing,
		real minRd,
		size_t& nbodiesOld,
		int nbodiesUp,
		int order)
	{
		double starttime, endtime;
		starttime = omp_get_wtime();

		//Число мультипроцессоров, заполняется функцией  setBlocks(blocks)
		int blocks;

		//Число ячеек дерева и тел
		int nnodes, nnodesUp;

		//Радиус вихря и параметр близости и их квадраты
		//real epssq = (real)(eps * eps);

		CudaSelect(0);
		setBlocks(blocks); //"достает" число блоков, равное числу мультипроцессоров (blocks - по ссылке)		

		nnodes = nbodies * 2;
		if (nnodes < 1024 * blocks)
			nnodes = 1024 * blocks;
		while ((nnodes & (32 - 1)) != 0)  // 32 - это размер варпа
			nnodes++;
		nnodes--;

		if (rebuild)
		{
			nnodesUp = nbodiesUp * 2;
			if (nnodesUp < 1024 * blocks)
				nnodesUp = 1024 * blocks;
			while ((nnodesUp & (32 - 1)) != 0)  // 32 - это размер варпа
				nnodesUp++;
			nnodesUp--;
		}

		KernelsOptimization();

		for (int i = 0; i < 6; i++)
			timing[i] = 0;

		if (rebuild)
		{
			if (nbodiesUp > nbodiesOld)
				timing[1] += memoryAllocate(ptrs, nnodesUp, nbodiesUp, (int)nbodiesOld, blocks, order);

			nbodiesOld = nbodiesUp;
			rebuildBaseTree(ptrs, nbodies, vtxl, (int)sizeof(realVortex), (int)Vortex2D::offsPos, 0, nnodes, order, timing, true, nullptr);
		}
		timing[5] += cuI1I2CalculationKernel2(ptrs, order, nnodes, nbodies, minRd, vtxl, i1l, i2l, epsastl);
		timing[6] = timing[1] + timing[2] + timing[3] + timing[4] + timing[5];



		endtime = omp_get_wtime();
		return endtime - starttime;
	}


	double wrapperDiffusiveVeloI0I3(const realVortex* vtxl, float* i0l, floatPoint* i3l, real* epsastl, real* ptr_r, CUDApointers& ptrs, bool rebuild, int nbodies, size_t nPan, real* visstr, double* timing, real* meanEps, real minRd, size_t& nbodiesOld, int nbodiesUp, int order)
	{
		double starttime, endtime;
		starttime = omp_get_wtime();

		//Число мультипроцессоров, заполняется функцией  setBlocks(blocks)
		int blocks;

		//Число ячеек дерева и тел
		int nnodes, nnodesUp;

		CudaSelect(0);
		setBlocks(blocks); //"достает" число блоков, равное числу мультипроцессоров (blocks - по ссылке)		

		nnodes = nbodies * 2;
		if (nnodes < 1024 * blocks)
			nnodes = 1024 * blocks;
		while ((nnodes & (32 - 1)) != 0)  // 32 - это размер варпа
			nnodes++;
		nnodes--;

		if (rebuild)
		{
			nnodesUp = nbodiesUp * 2;
			if (nnodesUp < 1024 * blocks)
				nnodesUp = 1024 * blocks;
			while ((nnodesUp & (32 - 1)) != 0)  // 32 - это размер варпа
				nnodesUp++;
			nnodesUp--;
		}

		KernelsOptimization();

		for (int i = 0; i < 6; i++)
			timing[i] = 0;

		if (rebuild)
		{
			if (nbodiesUp > nbodiesOld)
				timing[1] += memoryAllocate(ptrs, nnodesUp, nbodiesUp, (int)nbodiesOld, blocks, order);

			nbodiesOld = nbodiesUp;
			rebuildBaseTree(ptrs, nbodies, vtxl, (int)sizeof(realVortex), (int)Vortex2D::offsPos, 0, nnodes, order, timing, true, nullptr);
		}
		
		timing[5] += cuI0I3CalculationKernel2(ptrs, order, nnodes, nbodies, minRd, vtxl, i0l, i3l, epsastl, meanEps, (int)nPan, ptr_r, visstr);
		
		timing[6] = timing[1] + timing[2] + timing[3] + timing[4] + timing[5];

		endtime = omp_get_wtime();
		return endtime - starttime;
	}

}