/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: cuKnn.cu                                                         |
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
\brief Реализация поиска ближайших соседей на CUDA
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\author Серебровская Екатерина Александровна
\Version 1.14
\date 6 марта 2026 г.
*/


#include "cuKnn.cuh"

#include <utility>
#include "Vortex2D.h"
#include <omp.h>
#include <algorithm>
#include <iostream>

#include "cudaTreeInfo.h"


VectorsForKnn::VectorsForKnn() : reservedVtx(0), sortBufferSizeInBytes(0) {};

VectorsForKnn::~VectorsForKnn() {
	cudaFree(initdistPtr);
	cudaFree(initdistPtrSdvig);

	cudaFree(vtxPtr);

	cudaFree(mcdataPtr);
	cudaFree(mcdata_unsortedPtr);

	cudaFree(indexPtr);
	cudaFree(index_unsortedPtr);

	cudaFree(sortBuffer);
};

void VectorsForKnn::ResizeAll(int k, int nvtx)
{
	if ((nvtx > 0) && (nvtx > reservedVtx))
	{
		if (reservedVtx > 0)
		{
			cudaFree(vtxPtr);
			cudaFree(mcdataPtr);
			cudaFree(mcdata_unsortedPtr);
			cudaFree(indexPtr);
			cudaFree(index_unsortedPtr);
			cudaFree(initdistPtr);
			cudaFree(initdistPtrSdvig);
			cudaFree(sortBuffer);
		}

		reservedVtx = nvtx * 2;
		//std::cout << "knn: reservedVtx = " << reservedVtx << ", nvtx = " << nvtx << std::endl;

		cudaMalloc((void**)&vtxPtr, reservedVtx * sizeof(Vortex2D));		
		
		cudaMalloc((void**)&mcdataPtr, reservedVtx * sizeof(int));
		cudaMalloc((void**)&mcdata_unsortedPtr, reservedVtx * sizeof(int));		
		
		cudaMalloc((void**)&indexPtr, reservedVtx * sizeof(int));
		cudaMalloc((void**)&index_unsortedPtr, reservedVtx * sizeof(int));

		cudaMalloc((void**)&initdistPtr, reservedVtx * 2 * k * sizeof(prDoubleInt));
		cudaMalloc((void**)&initdistPtrSdvig, reservedVtx * 2 * k * sizeof(prDoubleInt));

		cudaMalloc((void**)&sortBuffer, reservedVtx * 64 * sizeof(char));
		sortBufferSizeInBytes = reservedVtx * 64;
	}
}


const int twoPowCodeLengthVar = (1 << codeLength);

//"Разрежение" двоичного представления беззнакового целого, вставляя по одному нулику между всеми битами
__host__ __device__ unsigned int ExpandBitscuda(unsigned int v)
{
	// вставит 1 нуль
	v = (v | (v << 8)) & 0x00FF00FF;      //  00000000`00000000`abcdefgh`ijklmnop 
	//                                      | 00000000`abcdefgh`ijklmnop`00000000
	//                                      = 00000000`abcdefgh`XXXXXXXX`ijklmnop
	//                                      & 00000000`11111111`00000000`11111111
	//                                      = 00000000`abcdefgh`00000000`ijklmnop

	v = (v | (v << 4)) & 0x0F0F0F0F;      //  00000000`abcdefgh`00000000`ijklmnop 
	//                                      | 0000abcd`efgh0000`0000ijkl`mnop0000
	//                                      = 0000abcd`XXXXefgh`0000ijkl`XXXXmnop
	//                                      & 00001111`00001111`00001111`00001111
	//                                      = 0000abcd`0000efgh`0000ijkl`0000mnop

	v = (v | (v << 2)) & 0x33333333;      //  0000abcd`0000efgh`0000ijkl`0000mnop 
	//                                      | 00abcd00`00efgh00`00ijkl00`00mnop00
	//                                      = 00abXXcd`00efXXgh`00ijXXkl`00mnXXop
	//                                      & 00110011`00110011`00110011`00110011
	//                                      = 00ab00cd`00ef00gh`00ij00kl`00mn00op

	v = (v | (v << 1)) & 0x55555555;      //  00ab00cd`00ef00gh`00ij00kl`00mn00op 
	//                                      | 0ab00cd0`0ef00gh0`0ij00kl0`0mn00op0
	//                                      = 0aXb0cXd`0eXf0gXh`0iXj0kXl`0mXn0oXp
	//                                      & 01010101`01010101`01010101`01010101
	//                                      = 0a0b0c0d`0e0f0g0h`0i0j0k0l`0m0n0o0p
	return v;
}


//Мортоновский код для пары из чисел типа double
//Исходное число - строго в диапазоне [0, 1 / (1.0 + 1/(2^codeLength - 1) ))
__host__ __device__ int Morton2Dcuda(const Point2D& r)
{
	const Point2D& rscale = twoPowCodeLengthVar * r;
	const unsigned int& xx = ExpandBitscuda((unsigned int)(rscale[0]));
	const unsigned int& yy = ExpandBitscuda((unsigned int)(rscale[1]));
	return (int)(yy | (xx << 1));
}


template <int k>
__device__ int BinSearchcuda(int iii, double2* beg, double x)
{
	int mid = -1;

	int low = 0;
	int high = k - 1;

	double2 buf[k];

#pragma unroll
	for (int i = 0; i < k; ++i)
		buf[i] = beg[i];

	if (x > buf[high].x)
	{
		return high + 1;
	}

	while (low <= high) {
		mid = (low + high) / 2;
				
		if (buf[mid].x == x)
		{
			return mid + 1; //выход из цикла
		}
		if (buf[mid].x < x)
			low = mid + 1;
		else
			high = mid - 1;
	}

	return mid;
}


template <int k>
__device__ void newSortcudaBuf(double2* beg) 
{	
	size_t cnt = 0;
	double2 buf[k];
	double2 result[k];

#pragma unroll
	for (size_t i = 0; i < k; i++)
		buf[i] = beg[i];

	for (size_t i = 0; i < k; i++) {
		double elem = buf[i].x;

		cnt = 0;
		for (size_t j = 0; j < k; j++) {
			cnt += (buf[j].x < elem);
		}
		//if (!AssumeDistinct) 
		for (size_t j = 0; j < i; j++)
			cnt += (buf[j].x == elem);
		result[cnt] = buf[i];
	}

#pragma unroll
	for (int i = 0; i < k; ++i)	
		beg[i] = result[i];
}


template <int k>
void /*__host__*/ __device__ newMergecuda(
	int iii,
	double2* beg,
	double2* candidateNN
)
{
	int locBuf[2*k];
	int counterBuf[2*k];
	int counterScanBuf[2*k];
	int offsetBuf[2*k];
	double2 updateNNBuf[2*k];
	double2 candidateNNBuf[2*k];
	double2 buf[k];

#pragma unroll
	for (size_t j = 0; j < k; ++j)
		buf[j] = beg[j];

#pragma unroll
	for (size_t j = 0; j < 2 * k; ++j)
		candidateNNBuf[j] = candidateNN[j];

#pragma unroll
	for (size_t j = 0; j < 2 * k; ++j)
		locBuf[j] = BinSearchcuda<k>(iii, buf, candidateNNBuf[j].x);
	
#pragma unroll
	for (int j = 0; j < 2 * k; ++j)
		counterBuf[j] = j & 1;

	for (int j = 0; j < 2 * k; ++j)
	{
		if (( locBuf[j] > 0) && ((locBuf[k] == k) 
			|| (__double_as_longlong(candidateNNBuf[j].y) == __double_as_longlong(buf[locBuf[j] - 1].y))))
			offsetBuf[j] = k + 1;
		else
			offsetBuf[j] = (counterBuf[2 * locBuf[j]])++;
	}

	counterScanBuf[0] = 0;

	for (int j = 1; j < 2 * k; ++j)
		counterScanBuf[j] = counterScanBuf[j - 1] + counterBuf[j - 1];

	size_t index;
	for (int j = 0; j < k; ++j)
	{
		index = counterScanBuf[2 * j + 1];
		if (index < k)
			updateNNBuf[index] = buf[j];
	}

	for (int j = 0; j < 2 * k; ++j)
	{
		if (2 * locBuf[j] < 2 * k)
		{
			index = counterScanBuf[2 * locBuf[j]] + offsetBuf[j];
			if (index < k)
				updateNNBuf[index] = candidateNNBuf[j];
		}
	}

#pragma unroll
	for (int s = 0; s < k; ++s)
		beg[s] = updateNNBuf[s];
	
}

#include "cuSort.cuh"


__device__ void calcCheckMergeVortices(const Vortex2D& vtxi, const Vortex2D& vtxk, double maxG, double cSP, double cRBP, double epsCol, int type, double& d2, int& result1, int& result2)
{
	int flagExit = false;
	int check = false;

	d2 = (vtxi.r() - vtxk.r()).length2();
	if (type == -1)
	{	
		result1 = false;
		result2 = true;
		return;
	}

	//линейное увеличение радиуса коллапса				    
	double mnog = std::max(1.0, /* 2.0 * */ (vtxi.r()[0] - cRBP) / cSP);
	double r2test = (epsCol * mnog) * (epsCol * mnog);

	if (type == 1)
		r2test *= 4.0; //Увеличение радиуса коллапса в 2 раза для коллапса вихрей разных знаков		
	
	if (d2 > 1.0 * r2test)
	{
		flagExit = true;
		result1 = flagExit;
		result2 = false;
		return;
	}

	const double gi = vtxi.g();
	const double gk = vtxk.g();

	if (d2 < r2test)
	{
		switch (type)
		{
		case 0:		
			check = (fabs(gi * gk) != 0.0) && (fabs(gi + gk) < (mnog * mnog) * maxG);
			break;
		case 1:
			check = (gi * gk < 0.0);
			break;
		case 2:
			check = (gi * gk > 0.0) && (fabs(gi + gk) < (mnog * mnog) * maxG);
			break;
		}
	}//if r2 < r2_test

	result1 = flagExit;
	result2 = check;
	return;
}



template <int k>
__global__ void findNeib(
	double2* initdistPtr,
	double2* initdistPtrSdvig,
	int* index_dev_raw, 
	Vortex2D* devVtxPtrV2D, 
	int n, 
	int sdvig, 
	double maxG, 
	double cSP,
	double cRBP, 
	double epsCol, 
	int type
	)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	double2* curPosition = nullptr;

	int flagExit, check;
	double d2;

	int cntr = 0;
	int search = i - 1;

	if (i < n)
	{
		const int id = index_dev_raw[i];
		const Vortex2D vtxi = devVtxPtrV2D[id];
		

		while ((cntr < k) && (search >= 0))
		{
			const Vortex2D vtxk = devVtxPtrV2D[index_dev_raw[search]];
			calcCheckMergeVortices(vtxi, vtxk, maxG, cSP, cRBP, epsCol, type, d2, flagExit, check);
			if (flagExit)
				break;

			if ((type==-1) || ((index_dev_raw[search] > id) && check))
			{
				//__syncthreads();
				curPosition = (sdvig ? initdistPtrSdvig : initdistPtr) + (id * 2 * k) + cntr;
				curPosition->x = d2;
				curPosition->y = __longlong_as_double(index_dev_raw[search]);
				++cntr;
			}
			--search;
		}

		__threadfence();
		__syncthreads();


#pragma unroll
		for (int w = cntr; w < k; ++w)			
		{
			curPosition = (sdvig ? initdistPtrSdvig : initdistPtr) + (id * 2 * k) + cntr;
			curPosition->x = 100000000.0;
			curPosition->y = __longlong_as_double(0);
			++cntr;			
		}
		__threadfence();
		__syncthreads();

		search = i + 1;

		while ((cntr < 2 * k) && (search < n))
		{
			const Vortex2D& vtxk = devVtxPtrV2D[index_dev_raw[search]];
			calcCheckMergeVortices(vtxi, vtxk, maxG, cSP, cRBP, epsCol, type, d2, flagExit, check);
			if (flagExit)
				break;

			if ((type==-1) || ((index_dev_raw[search] > id) && check))
			{
				curPosition = (sdvig ? initdistPtrSdvig : initdistPtr) + (id * 2 * k) + cntr;
				curPosition->x = d2;
				curPosition->y = __longlong_as_double(index_dev_raw[search]);
				++cntr;
			}
			++search;
		}
		__threadfence();
		__syncthreads();

#pragma unroll
		for (int w = cntr; w < 2 * k; ++w)
		{
			curPosition = (sdvig ? initdistPtrSdvig : initdistPtr) + (id * 2 * k) + cntr;
			curPosition->x = 100000000.0;
			curPosition->y = __longlong_as_double(0);
			++cntr;			
		}
		__threadfence();
		__syncthreads();

		if (sdvig == 0)					
			newSortcudaBuf<2*k>(initdistPtr + (id * 2 * k));
		else
		{			
			newMergecuda<k>(id, initdistPtr + (id * 2 * k), initdistPtrSdvig + (id * 2 * k));
			newSortcudaBuf<2*k>(initdistPtr + (id * 2 * k));
		}
		__threadfence();
		__syncthreads();
	}
};


__global__ void Fill123(int* ptr, int n)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < n)	
		ptr[i] = i;
};



__global__ void CalcMortonCodes(Vortex2D* vtxPtr, int n, Point2D LL, double scale, int sdvig, Point2D shift005, int* mcdata_unsortedPtr)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n)
	{
		Point2D sh = ((vtxPtr[i].r() - LL) * (0.75 / scale) + sdvig * shift005);
		mcdata_unsortedPtr[i] = Morton2Dcuda(sh);
	};
}

template <int k>
__global__ void reordNeib(int n,
	prDoubleInt* initdistPtr,
	prDoubleInt* initdistPtrUpdate)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < n)
	{
#pragma unroll
		for (int q = 0; q < k; ++q)
			initdistPtrUpdate[i * k + q] = initdistPtr[i * (2 * k) + q];

	}//if n
}

template <int k>
double kNNcuda(Point2D minr, Point2D maxr, int blocks, const std::vector<Vortex2D>& vtx,
	std::vector<std::pair<double, size_t>>& initdist,
	VectorsForKnn*& vecForKnn,
	double cSP, double cRBP, double maxG, double epsCol, int type
)
{
	const int nSdvig = 5;

	double t1 = -omp_get_wtime();

	const size_t n = vtx.size();

	if (vecForKnn == nullptr)
		vecForKnn = new VectorsForKnn;

	vecForKnn->ResizeAll((int)k, (int)n);
	cudaMemcpy(vecForKnn->vtxPtr, vtx.data(), n * sizeof(Vortex2D), cudaMemcpyHostToDevice);
	t1 += omp_get_wtime();

	double tAll = -omp_get_wtime();
	
	double t2 = -omp_get_wtime();

	double scale = std::max(maxr[0] - minr[0], maxr[1] - minr[1]);
	
	
	Point2D shift005;
	shift005[0] = shift005[1] = 0.05;

	dim3 blockA, threadA = 32;
	blockA.x = ((int)n + threadA.x - 1) / threadA.x;

	t2 += omp_get_wtime();

	double t3 = 0.0;
	double t4 = 0.0;
	
	for (size_t sdvig = 0; sdvig < nSdvig; ++sdvig)
	{
		//масштабирование координат вихрей в [0;0.75)^2, поиск их мортоновских кодов 	
		//Сдвигаем все точки и считаем новые коды		
		t3 -= omp_get_wtime();

		CalcMortonCodes<<<blockA, threadA>>> (vecForKnn->vtxPtr, (int)n, minr, scale, (int)sdvig, shift005, vecForKnn->mcdata_unsortedPtr);

		Fill123<<<blockA, threadA>>>(vecForKnn->index_unsortedPtr, (int)n);

		BHcu::RadixSortFromCUBReservedMem(
			vecForKnn->mcdata_unsortedPtr,
			vecForKnn->mcdataPtr,
			vecForKnn->index_unsortedPtr,
			vecForKnn->indexPtr,
			(int)(n), 0, 2 * codeLength,
			vecForKnn->sortBuffer,
			vecForKnn->sortBufferSizeInBytes
			);
		t3 += omp_get_wtime();

		t4 -= omp_get_wtime();

		findNeib<k><<<blockA, threadA>>>(
			(double2*)vecForKnn->initdistPtr,
			(double2*)vecForKnn->initdistPtrSdvig,
			vecForKnn->indexPtr,
			vecForKnn->vtxPtr,			
			(int)n,
			(int)sdvig,
			maxG,
			cSP,
			cRBP,
			epsCol,
			type
			);
		cudaDeviceSynchronize();
		t4 += omp_get_wtime();
	}//for sdvig

	double t5 = -omp_get_wtime();
	//std::cout << "tB-tA = " << tB-tA << std::endl;

	reordNeib<k> << <blockA, threadA >> > ((int)n, vecForKnn->initdistPtr, vecForKnn->initdistPtrSdvig);
	cudaDeviceSynchronize();
	t5 += omp_get_wtime();

	double t6 = -omp_get_wtime();
	double tMem = -omp_get_wtime();
	cudaMemcpy(initdist.data(), (char*)vecForKnn->initdistPtrSdvig, n * k * sizeof(prDoubleInt), cudaMemcpyDeviceToHost);

	t6 += omp_get_wtime();
	tAll += omp_get_wtime();
	
	//std::cout << "t1 = " << t1 * 1000 << ", t2 = " << t2 * 1000 << ", t3 = " << t3 * 1000 << ", t4 = " << t4 * 1000 << ", t5 = " << t5 * 1000 << ", t6 = " << t6 * 1000 << std::endl;
	return tAll;
}


template
double kNNcuda<3>(Point2D minr, Point2D maxr, int blocks, const std::vector<Vortex2D>& vtx,
	std::vector<std::pair<double, size_t>>& initdist,
	VectorsForKnn*& vecForKnn,
	double cSP, double cRBP, double maxG, double epsCol, int typ);

template
double kNNcuda<4>(Point2D minr, Point2D maxr, int blocks, const std::vector<Vortex2D>& vtx,
	std::vector<std::pair<double, size_t>>& initdist,
	VectorsForKnn*& vecForKnn,
	double cSP, double cRBP, double maxG, double epsCol, int typ);

template
double kNNcuda<5>(Point2D minr, Point2D maxr, int blocks, const std::vector<Vortex2D>& vtx,
	std::vector<std::pair<double, size_t>>& initdist,
	VectorsForKnn*& vecForKnn,
	double cSP, double cRBP, double maxG, double epsCol, int typ);

template
double kNNcuda<6>(Point2D minr, Point2D maxr, int blocks, const std::vector<Vortex2D>& vtx,
	std::vector<std::pair<double, size_t>>& initdist,
	VectorsForKnn*& vecForKnn,
	double cSP, double cRBP, double maxG, double epsCol, int typ);