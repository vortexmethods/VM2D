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
			time = McuBoundingBoxKernelFree(nullptr, maxpt, minpt, npoints, pointsl);
			return time;
		}

		~CudaCalcGab()
		{
			cudaDelete(maxpt);
			cudaDelete(maxpt);
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
				npoints_, pointsl_);
			return timeGab + timeCodes;
		}

		~CudaSorter()
		{
			cudaDelete(pointsMortonCodesKeyl);
			cudaDelete(pointsMortonCodesKeyUnsortl);
			cudaDelete(pointsMortonCodesIdxl);
			cudaDelete(pointsMortonCodesIdxUnsortl);
		}
	};


	void rebuildBaseTree(CUDApointers& ptrs, const int nbodies, const realVortex* vtxl, int nnodes, int order, double* timing)
	{
		timing[0] += cuInitializationKernel();
		timing[1] += McuBoundingBoxKernel(ptrs, nbodies, vtxl);
		
		timing[2] += McuMortonCodesKernel(ptrs, nbodies, vtxl);
		timing[2] += McuMortonInternalNodesKernel(ptrs, nbodies);
		timing[2] += McuMortonInternalCellsGeometryKernel(ptrs, nbodies, nnodes);

		timing[3] += cuClearKernel2(ptrs, order, nnodes, nbodies);

		timing[4] += cuAABBKernel2(ptrs, nnodes, nbodies, vtxl);

		timing[4] += cuClearKernel2(ptrs, order, nnodes, nbodies);

		timing[4] += cuSummarizationKernel2(ptrs, order, nnodes, nbodies, vtxl);
	}


	double memoryAllocate(CUDApointers& ptrs, int nnodes, int nbodies, int nbodiesOld, int blocks, int order)
	{
		double starttime, endtime;
		starttime = omp_get_wtime();

		if (nbodiesOld > 0)
		{
			//std::cout << "BHgpu: free CUDA-memory" << std::endl;
			cudaDelete(ptrs.massl);

			cudaDelete(ptrs.momsl);
			cudaDelete(ptrs.El);

			cudaDelete(ptrs.maxrl);
			cudaDelete(ptrs.minrl);

			///For Morton tree
			cudaDelete(ptrs.MmortonCodesKeyUnsortl);
			cudaDelete(ptrs.MmortonCodesIdxUnsortl);
			cudaDelete(ptrs.MmortonCodesKeyl);
			cudaDelete(ptrs.MmortonCodesIdxl);

			cudaDelete(ptrs.Mposl);
			cudaDelete(ptrs.Mlowerl);
			cudaDelete(ptrs.Mupperl);
			cudaDelete(ptrs.Mparentl);
			cudaDelete(ptrs.Mchildl);
			cudaDelete(ptrs.Mrangel);

			cudaDelete(ptrs.MlevelUnsortl);
			cudaDelete(ptrs.MlevelSortl);
			cudaDelete(ptrs.MindexUnsortl);
			cudaDelete(ptrs.MindexSortl);
			cudaDelete(ptrs.MindexSortTl);
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


	// N body <=> wake to wake
	double wrapperInfluence(const realVortex* vtxl, realPoint* vell, 
		real* epsastl, CUDApointers& ptrs, 
		int nbodies, double* timing, real eps, real theta, 
		size_t& nbodiesOld, int nbodiesUp, 
		int order,
		size_t nAfls, size_t* nVtxs, double** ptrVtxs)
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


		nnodesUp = nbodiesUp * 2;
		if (nnodesUp < 1024 * blocks)
			nnodesUp = 1024 * blocks;
		while ((nnodesUp & (32 - 1)) != 0)  // 32 - это размер варпа
			nnodesUp++;
		nnodesUp--;

		KernelsOptimization();

		for (int i = 0; i < 6; i++)
			timing[i] = 0;

		if (nbodiesUp > nbodiesOld)
			timing[1] += memoryAllocate(ptrs, nnodesUp, nbodiesUp, (int)nbodiesOld, blocks, order);

		nbodiesOld = nbodiesUp;
		rebuildBaseTree(ptrs, nbodies, vtxl, nnodes, order, timing);

		timing[5] += cuForceCalculationKernel2points(ptrs, order, nnodes, nbodies, itolsq, epssq, vtxl, 
			ptrs.MmortonCodesIdxl, nbodies, vtxl,  vell, true, epsastl, nAfls, nVtxs, ptrVtxs);
		timing[6] = timing[1] + timing[2] + timing[3] + timing[4] + timing[5];

		endtime = omp_get_wtime();
		return endtime - starttime;

	}


	// wake to points
	double wrapperInfluenceToPoints(
		const realVortex* vtxl, const realVortex* pointsl, realPoint* vell, real* epsastl,
		CUDApointers& ptrs, bool rebuild, int nbodies, int npoints, double* timing, real eps, real theta,
		size_t& nbodiesOld, int nbodiesUp, int order,
		size_t nAfls, size_t* nVtxs, double** ptrVtxs)
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
			rebuildBaseTree(ptrs, nbodies, vtxl, nnodes, order, timing);
		}

		CudaSorter srt(npoints, pointsl);
		timing[5] += srt.calc();

		timing[5] += cuForceCalculationKernel2points(ptrs, order, nnodes, nbodies, itolsq, epssq, vtxl, srt.pointsMortonCodesIdxl, npoints, pointsl, vell, true, epsastl,
			nAfls, nVtxs, ptrVtxs);


		timing[6] = timing[1] + timing[2] + timing[3] + timing[4] + timing[5];

		endtime = omp_get_wtime();
		return endtime - starttime;

	}



	double wrapperInfluenceToRHS(
		const realVortex* dev_ptr_vt,  //вихри в следе
		const double* dev_ptr_pt,      //начала и концы панелей
		double* dev_ptr_rhs,           //куда сохранить результат (для T0 и верхней половины T1)
		double* dev_ptr_rhslin,        //куда сохранить результат (для нижней половины T1)

		CUDApointers& ptrs,     //указатели на делево вихрей
		bool rebuild,           //признак перестроения делева вихрей

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

		KernelsOptimization();

		for (int i = 0; i < 6; i++)
			timingsToRHS[i] = 0;

		if (rebuild)
		{
			if (nbodiesUp > nbodiesOld)
				timingsToRHS[1] += memoryAllocate(ptrs, nnodesUp, nbodiesUp, (int)nbodiesOld, blocks, order);

			nbodiesOld = nbodiesUp;
			rebuildBaseTree(ptrs, nvt, dev_ptr_vt, nnodes, order, timingsToRHS);
		}

		Vortex2D* pointsl = (Vortex2D*)cudaNew(nTotPan, sizeof(Vortex2D));
		realPoint* El = (realPoint*)cudaNew(nTotPan * order, sizeof(realPoint));

		McuVerticesToControlPoints(nTotPan, (double*)dev_ptr_pt, (double*)pointsl);

		CudaSorter srt(nTotPan, pointsl);
		timingsToRHS[5] += srt.calc();

		double* ptrToLin = nullptr;
		if (scheme == 1)
			ptrToLin = dev_ptr_rhslin;

		timingsToRHS[5] += cuRhsCalculationKernel(ptrs, order, nnodes, nvt, itolsq, dev_ptr_vt, 
			srt.pointsMortonCodesIdxl, El,
			nTotPan, dev_ptr_pt, (const real*)pointsl, dev_ptr_rhs, ptrToLin);

		cudaDelete(El);
		cudaDelete(pointsl);

		timingsToRHS[6] = timingsToRHS[1] + timingsToRHS[2] + timingsToRHS[3] + timingsToRHS[4] + timingsToRHS[5];

		endtime = omp_get_wtime();
		return endtime - starttime;
		
	}

	double wrapperDiffusiveVelo(const realVortex* vtxl, real* i1l, realPoint* i2l, real* epsastl, CUDApointers& ptrs, bool rebuild, int nbodies, double* timing, real eps, real theta, size_t& nbodiesOld, int nbodiesUp, int order,
		size_t nAfls, size_t* nVtxs, double** ptrVtxs)
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
			rebuildBaseTree(ptrs, nbodies, vtxl, nnodes, order, timing);
		}

		timing[5] += cuDiffVelCalculationKernel2(ptrs, order, nnodes, nbodies, itolsq, epssq, vtxl, i1l, i2l, true, epsastl, nAfls, nVtxs, ptrVtxs);
		timing[6] = timing[1] + timing[2] + timing[3] + timing[4] + timing[5];

		endtime = omp_get_wtime();
		return endtime - starttime;
	}




}