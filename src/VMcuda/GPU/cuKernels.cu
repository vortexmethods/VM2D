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
| File name: cuKernels.cu                                                     |
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

/*!
\file
\brief Реализация CUDA-ядер
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 29 августа 2023 г.
*/

#include "cuKernels.cuh"

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <cuda.h>
#include "operations.cuh"
#include "cuSort.cuh"

#include "Point2D.h"

#define WARPSIZE 32

//было 28 - и иногда падало
#define MAXDEPTH 32 

#define BLOCKD 32

#define codeLength 14
#define twoPowCodeLength (1 << codeLength)
#define rbound (1 - (real)1.0 / twoPowCodeLength)

#define idpi ((real)0.15915494309189534)

namespace BHcu
{
    int blocks;
__device__ volatile int bottomd;
__device__ unsigned int blkcntd;
__device__ real iQuadSideFactor;
__device__ real2 shift;

void setBlocks(int& blocks_)
{
     blocks_ = blocks;
}

void CudaSelect(int dev)
{
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0) {
        fprintf(stderr, "There is no device supporting CUDA\n");
        exit(-1);
    }

    if ((dev < 0) || (deviceCount <= dev)) {
        fprintf(stderr, "There is no device %d\n", dev);
        exit(-1);
    }
    cudaSetDevice(dev);


    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);

    blocks = deviceProp.multiProcessorCount;

/*
    printf("\n");
    printf("                          GPU Device Properties                         \n");
    printf("------------------------------------------------------------------------\n");
    printf("Name:                                  %s\n", properties.name); 
    printf("CUDA driver/runtime version:           %d.%d/%d.%d\n", driverVersion / 1000, (driverVersion % 100) / 10, runtimeVersion / 1000, (runtimeVersion % 100) / 10);
    printf("CUDA compute capabilitiy:              %d.%d\n", properties.major, properties.minor);
    printf("Number of multiprocessors:             %d\n", properties.multiProcessorCount);
    
    if (false)
    {
        printf("GPU clock rate:                        %d (MHz)\n", properties.clockRate / fact);
        printf("Memory clock rate:                     %d (MHz)\n", properties.memoryClockRate / fact);
        printf("Memory bus width:                      %d-bit\n", properties.memoryBusWidth);
        printf("Theoretical memory bandwidth:          %d (GB/s)\n", (properties.memoryClockRate / fact * (properties.memoryBusWidth / 8) * 2) / fact);
        printf("Device global memory:                  %d (MB)\n", (int)(properties.totalGlobalMem / (fact * fact)));
        printf("Shared memory per block:               %d (KB)\n", (int)(properties.sharedMemPerBlock / fact));
        printf("Constant memory:                       %d (KB)\n", (int)(properties.totalConstMem / fact));
        printf("Maximum number of threads per block:   %d\n", properties.maxThreadsPerBlock);
        printf("Maximum thread dimension:              [%d, %d, %d]\n", properties.maxThreadsDim[0], properties.maxThreadsDim[1], properties.maxThreadsDim[2]);
        printf("Maximum grid size:                     [%d, %d, %d]\n", properties.maxGridSize[0], properties.maxGridSize[1], properties.maxGridSize[2]);
    }
    printf("------------------------------------------------------------------------\n");  
*/
    
}


void cudaDelete(void* cudaPtr)
{
	cudaFree(cudaPtr);
}


void* cudaNew(int n, size_t sizeType)
{
	void* cudaPtr;
	cudaMalloc(&cudaPtr, sizeType * n);
	CudaTest("couldn't allocate device memory");

	return cudaPtr;
}

void cudaCopyVecToDevice(void* hostPtr, void* cudaPtr, size_t n, size_t typeSize)
{
	cudaMemcpy(cudaPtr, hostPtr, typeSize * n, cudaMemcpyHostToDevice);
	CudaTest("couldn't copy data from host to device");
}

void cudaCopyVecFromDevice(void* cudaPtr, void* hostPtr, size_t n, size_t typeSize)
{
	cudaMemcpy(hostPtr, cudaPtr, typeSize * n, cudaMemcpyDeviceToHost);
	CudaTest("couldn't copy data from device to host");
}


//////////////////
/// Error TEST
//////////////////


void CudaTest(const char* msg)
{
    cudaError_t e;

    //cudaThreadSynchronize();
    cudaDeviceSynchronize();
    if (cudaSuccess != (e = cudaGetLastError())) {
        fprintf(stderr, "%s: %d\n", msg, e);
        fprintf(stderr, "%s\n", cudaGetErrorString(e));
        exit(-1);
    }
}


//////////////////
/// CUDA Kernels
//////////////////


/******************************************************************************/
/*** initialize memory ********************************************************/
/******************************************************************************/

__global__ void InitializationKernel()
{
    blkcntd = 0;
}


/******************************************************************************/
/*** compute center and radius ************************************************/
/******************************************************************************/
__global__
__launch_bounds__(THREADS1, FACTOR1)
void MBoundingBoxKernel(
    const int nbodiesd, 
    const real3* __restrict vtxd, 
    real2* __restrict Mposd, 
    volatile real2* __restrict maxrd, 
    volatile real2* __restrict minrd)
{
    register int i, j, k, inc;
    register real2 val;
    register real2 minr, maxr;
    __shared__ volatile real2 sminr[THREADS1], smaxr[THREADS1];

    // initialize with valid data (in case #bodies < #threads)
    minr.x = maxr.x = vtxd[0].x;
    minr.y = maxr.y = vtxd[0].y;

    // scan all bodies
    i = threadIdx.x;
    inc = THREADS1 * gridDim.x;
    
    for (j = i + blockIdx.x * THREADS1; j < nbodiesd; j += inc) 
    {
        val.x = vtxd[j].x;
        val.y = vtxd[j].y;

        minr.x = realmin(minr.x, val.x);
        maxr.x = realmax(maxr.x, val.x);

        minr.y = realmin(minr.y, val.y);
        maxr.y = realmax(maxr.y, val.y);
    }

    // reduction in shared memory
    sminr[i].x = minr.x;
    smaxr[i].x = maxr.x;
    sminr[i].y = minr.y;
    smaxr[i].y = maxr.y;

    for (j = THREADS1 / 2; j > 0; j /= 2) {
        __syncthreads();
        if (i < j) {
            k = i + j;
            sminr[i].x = minr.x = realmin(minr.x, sminr[k].x);
            smaxr[i].x = maxr.x = realmax(maxr.x, smaxr[k].x);
            sminr[i].y = minr.y = realmin(minr.y, sminr[k].y);
            smaxr[i].y = maxr.y = realmax(maxr.y, smaxr[k].y);
        }
    }

    // write block result to global memory
    if (i == 0) {
        k = blockIdx.x;
        minrd[k].x = minr.x;
        maxrd[k].x = maxr.x;
        minrd[k].y = minr.y;
        maxrd[k].y = maxr.y;

        __threadfence();

        inc = gridDim.x - 1;
        if (inc == atomicInc(&blkcntd, inc)) {

            // I'm the last block, so combine all block results
            for (j = 0; j <= inc; j++) {
                minr.x = realmin(minr.x, minrd[j].x);
                maxr.x = realmax(maxr.x, maxrd[j].x);
                minr.y = realmin(minr.y, minrd[j].y);
                maxr.y = realmax(maxr.y, maxrd[j].y);
            }

            // create root node 
            if (Mposd != nullptr)
            {
                Mposd[0].x = (minr.x + maxr.x) / 2;
                Mposd[0].y = (minr.y + maxr.y) / 2;
            }

            minrd[0].x = minr.x;
            minrd[0].y = minr.y;
            maxrd[0].x = maxr.x;
            maxrd[0].y = maxr.y;
        }
    }
    
}

/******************************************************************************/
/*** Morton codes *************************************************************/
/******************************************************************************/

__global__
void MMortonCodesKernel (
    const int nbodies, 
    const real3* __restrict vtxd, 
    int* __restrict MmortonCodesKeyUnsortd, 
    int* __restrict MmortonCodesIdxUnsortd,
    const real2* maxrd, 
    const real2* minrd)
{
	int bdy = blockDim.x * blockIdx.x + threadIdx.x;

    //if (bdy == 0)
    //    printf("x: %f -- %f, y: %f -- %f\n", minrd->x, maxrd->x, minrd->y, maxrd->y);

    register real lmax, lx0, ly0, quadSideFactor;
   
	if (bdy < nbodies)
	{
        lmax = realmax(maxrd->x - minrd->x, maxrd->y - minrd->y);
        lx0 = (minrd->x + maxrd->x) / 2; //координаты центра
        ly0 = (minrd->y + maxrd->y) / 2;

        quadSideFactor = rbound / lmax; //1;

        real x0 = (vtxd[bdy].x - lx0) * quadSideFactor + rbound/2;
        real y0 = (vtxd[bdy].y - ly0) * quadSideFactor + rbound/2;

       
        real x = twoPowCodeLength * x0;
        real y = twoPowCodeLength * y0;

		unsigned int xx = MExpandBits((unsigned int)x);
		unsigned int yy = MExpandBits((unsigned int)y);
		MmortonCodesKeyUnsortd[bdy] = yy | (xx << 1);
		MmortonCodesIdxUnsortd[bdy] = bdy;
	}

    if (bdy == 0)
    {
        iQuadSideFactor = 1 / quadSideFactor;
        shift.x = (rbound / 2) * iQuadSideFactor - lx0;
        shift.y = (rbound / 2) * iQuadSideFactor - ly0;
    }
}


/******************************************************************************/
/*** Morton Internal nodes tree build *****************************************/
/******************************************************************************/
__global__
void MMortonInternalNodesKernel(
    const int nbodies, 
    const int* __restrict MmortonCodesKeyd, 
    int* __restrict Mparentd, 
    int2* __restrict Mchildd, 
    int2* __restrict Mranged)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    
    if (i < nbodies - 1)
    {
        int d = sign(Delta(i, i + 1, nbodies, MmortonCodesKeyd) - Delta(i, i - 1, nbodies, MmortonCodesKeyd));
        int delta_min = Delta(i, i - d, nbodies, MmortonCodesKeyd);

        int Lmax = 2;
        while (Delta(i, i + Lmax * d, nbodies, MmortonCodesKeyd) > delta_min)
            Lmax *= 2;

        int L = 0;
        for (int t = (Lmax >> 1); t >= 1; t >>= 1)
            if (Delta(i, i + (L + t) * d, nbodies, MmortonCodesKeyd) > delta_min)
                L += t;

        int j = i + L * d;

        int delta_node = Delta(i, j, nbodies, MmortonCodesKeyd);

        int s = 0;
        for (int p = 1, t = ceilhalf(L); L > (1 << (p - 1)); ++p, t = ceilpow2(L, p))
        {
            int dl = Delta(i, i + (s + t) * d, nbodies, MmortonCodesKeyd);
            if (dl > delta_node)
                s += t;
        }//for p


        int gamma = i + s * d +   d * (d < 0);   //последнее слагаемое = std::min(d, 0);

        int Mmin = min(i, j);
        int Mmax = max(i, j);
        
        const int& left = gamma;
        const int& right = gamma + 1;

        // Левый потомок - лист или внутренний узел
        int childLeft = Mchildd[i].x = (Mmin == gamma) * nbodies + left;
        
        Mranged[childLeft].x = Mmin;
        Mranged[childLeft].y = gamma;
        Mparentd[childLeft] = i;

        // Правый потомок - лист или внутренний узел
        int childRight = Mchildd[i].y = (Mmax == gamma + 1) * nbodies + right;

        Mranged[childRight].x = gamma+1;
        Mranged[childRight].y = Mmax;
        Mparentd[childRight] = i;
    }
}

/******************************************************************************/
/*** Morton Internal nodes geometry calculation *******************************/
/******************************************************************************/
__global__
void MMortonInternalCellsGeometryKernel(
    const int nbodies,
    const int* __restrict MmortonCodesKeyd,   
    const int2* __restrict Mranged,
    int* __restrict MlevelUnsortd,
    int* __restrict MindexUnsortd,
    const real2* maxrd,
    const real2* minrd
)
{
    int cell = blockDim.x * blockIdx.x + threadIdx.x;
 
    if (cell < nbodies - 1)
    {
        int prLength = min(Delta(Mranged[cell].x, Mranged[cell].y, nbodies, MmortonCodesKeyd), 2 * codeLength);               
        prLength -= min(Delta(Mranged[0].x, Mranged[0].y, nbodies, MmortonCodesKeyd), 2 * codeLength);
        
        MlevelUnsortd[cell] = prLength;
        MindexUnsortd[cell] = cell;
    }

}//MMortonInternalCellsGeometryKernel(...)


/******************************************************************************/
/*** permutation list transposition *******************************************/
/******************************************************************************/
__global__
void MTransposeIndexKernel(
    const int nbodiesd, const int nnodesd,
    const int* __restrict MindexSortd, 
    int* __restrict MindexSortTd)
{
    register const int cell = blockDim.x * blockIdx.x + threadIdx.x;
    register const int newcell = MindexSortd[cell];

    if (cell < nbodiesd - 1)
        MindexSortTd[newcell] = cell;       

}//MTransposeIndexKernel



/******************************************************************************/
/*** build tree ***************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(1024, 1)
void ClearKernel2(
    const int nnodesd, const int nbodiesd, 
    volatile int* __restrict massd)
{
    register int k, inc, bottom;

    bottom = nnodesd - (nbodiesd - 1); //bottomd;
    inc = blockDim.x * gridDim.x;
    k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size
    if (k < bottom) k += inc;

// 0 1 ... (nb-1)  (nb+0) ... (nb+nb-2)
// --------------  --------------------
//     bodies              cells

    // iterate over all cells assigned to thread
    while (k < nnodesd) 
	{
        massd[nnodesd - 1 - k] = -1;
        k += inc;
    }
}



__global__
__launch_bounds__(THREADS3, FACTOR3)
void AABBKernel2(
    const int nnodesd, const int nbodiesd,
    const int2* __restrict Mchildd,
    volatile int* __restrict massd,    
    const real3* __restrict vtxd,     
    const int* __restrict MmortonCodesIdxd,
    real2* __restrict Mposd, 
    real2* __restrict Mlowerd,
    real2* __restrict Mupperd,    
    const int* __restrict MindexSortd, const int* __restrict MindexSortTd
)
{
    register int i, j, k, ch, inc, flag;

    register real2 lower[2];
    register real2 upper[2];



    register int m, cm;

    inc = blockDim.x * gridDim.x;
    k = ((nnodesd - (nbodiesd - 1)) & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size

    if (k < (nnodesd - (nbodiesd - 1)))
        k += inc;

    //MortonTree:
    // 0 1 2 ... (nb-2) x (nb+0) (nb+1) (nb+2) ... (nb+(nb-1))
    // ----------------   -----------------------------------
    //      cells                         bodies

    //Martin's tree:
    // 0 1 2 ... (nb-1) x x x x (nn-(nb-1)) ... (nn-2) (nn-1)
    // ----------------          ----------------------------
    //      bodies                 sorted and reversed cells

    flag = 0;
    j = 0;
    // iterate over all cells assigned to thread
    while (k < nnodesd)
    {
        if (massd[nnodesd - 1 - k] >= 0)
        {
            k += inc;
        }
        else
        {
            j = 2;
            for (i = 0; i < 2; i++) {

                //computation of child[k*2+i]
                register const int srt = MindexSortd[(nnodesd - 1) - k];
                int chd = i * Mchildd[srt].y + (1 - i) * Mchildd[srt].x;   // i==0 => .x;  i==1 => .y
                ch = (chd >= nbodiesd) ? chd - nbodiesd : (nnodesd - 1) - MindexSortTd[chd];

                if ((chd >= nbodiesd) || (massd[nnodesd - 1 - ch] >= 0))
                    j--;
            }

            if (j == 0) {
                // all children are ready
                const int kch = ((nnodesd - 1) - k);
                cm = 0;

                const register int sortedCell = MindexSortd[(nnodesd - 1) - k];


                const int2 chdPair = Mchildd[sortedCell];

                for (i = 0; i < 2; i++)
                {
                    //computation of ch = child[k*2+i]
                    const int chd = i * chdPair.y + (1 - i) * chdPair.x;
                    if (chd >= nbodiesd)
                    {
                        ch = chd - nbodiesd;
                        const register int sortedBody = MmortonCodesIdxd[ch];
                        lower[i] = upper[i] = real2{ vtxd[sortedBody].x, vtxd[sortedBody].y };

                        m = 1;
                    }
                    else
                    {
                        register const int srtT = MindexSortTd[chd];
                        ch = (nnodesd - 1) - srtT;
                        lower[i] = Mlowerd[chd];
                        upper[i] = Mupperd[chd];

                        m = massd[nnodesd - 1 - ch];
                    }
                    // add child's contribution

                    const int kchSort = MindexSortd[kch];
                    const real2 up = real2{ ::fmax(upper[0].x, upper[1].x), ::fmax(upper[0].y, upper[1].y) };
                    const real2 lo = real2{ ::fmin(lower[0].x, lower[1].x), ::fmin(lower[0].y, lower[1].y) };

                    Mupperd[kchSort] = up;                                        
                    Mlowerd[kchSort] = lo;                    
                    Mposd[kchSort] = (real)0.5 * (up + lo);                    
                    cm += m;
                }
                flag = 1;
            }
        }
        __threadfence();

        if (flag != 0) {
            massd[nnodesd - 1 - k] = cm;
            k += inc;
            flag = 0;
        }
    }
}

/******************************************************************************/
/*** compute center of mass ***************************************************/
/******************************************************************************/
#include "ShiftKernels/IncludeKer.cu"

/******************************************************************************/
/*** compute force ************************************************************/
/******************************************************************************/

__global__
//__launch_bounds__(THREADS5rhs, FACTOR5)
void ForceCalculationKernel2points(
    const int nnodesd, const int nbodiesd,
    const real itolsqd, const real epssqd,
    const int2* __restrict Mchildd,
    const int order, const real2* __restrict momsd,
    const real3* __restrict vtxd,
    const int* __restrict MmortonCodesIdxd,
    const real2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
    const int npointsd, const real3* __restrict pointsd,
    const int* __restrict MmortonCodesIdxPointsd,
    real2* __restrict veld,
    const real2* __restrict Mlowerd,
    const real2* __restrict Mupperd,
    bool calcEpsAst,
    real* __restrict epsast,
    size_t nAfls, size_t* nVtxs, double** ptrVtxs)
{
    register int j, k, n, depth, base, sbase, pd, nd;
    register real2 p, v, dr, ps;
    register real r2;
    register const real2* mom;

    register real d_1 = 1e+5;
    register real d_2 = 1e+5;
    register real d_3 = 1e+5;
    register real dst23, dst12;

    register real2 th;


    __shared__ volatile int pos[MAXDEPTH * THREADS5rhs / WARPSIZE], node[MAXDEPTH * THREADS5rhs / WARPSIZE];



    // figure out first thread in each warp (lane 0)
    base = threadIdx.x / WARPSIZE;
    sbase = base * WARPSIZE;
    j = base * MAXDEPTH;
    //diff = threadIdx.x - sbase;

    __syncthreads();
    __threadfence_block();

    // iterate over all bodies assigned to thread
    for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npointsd; k += blockDim.x * gridDim.x)
    {
        d_1 = d_2 = d_3 = 1e+5;

        const int indexOfPoint = MmortonCodesIdxPointsd[k]; //k;//MmortonCodesIdxd[k];
        //p = real2{ vtxd[indexInParticles].x, vtxd[indexInParticles].y };
        p = real2{ pointsd[indexOfPoint].x, pointsd[indexOfPoint].y };

        v.x = 0;
        v.y = 0;

        // initialize iteration stack, i.e., push root node onto stack
        depth = j;
        if (sbase == threadIdx.x)
        {
            pos[j] = 0;
            node[j] = nnodesd - 1;
        }

        do
        {
            // stack is not empty
            pd = pos[depth];
            nd = node[depth];

            register int2 chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];

            register real gm;
            register real2 sumSide2;
            bool isVortex;

            while (pd < 2)
            {
                // node on top of stack has more children to process

                // load child pointer
                //computation of n = childd[nd + pd] (pd = 0 или pd = 1)
                int chd = pd * chBoth.y + (1 - pd) * chBoth.x;
                ++pd;

                isVortex = (chd >= nbodiesd);

                if (isVortex)
                {
                    n = chd - nbodiesd;
                    ps = real2{ vtxd[MmortonCodesIdxd[n]].x, vtxd[MmortonCodesIdxd[n]].y };
                    gm = vtxd[MmortonCodesIdxd[n]].z;
                    sumSide2 = real2{ (real)0, (real)0 };
                }
                else
                {
                    register const int srtT = MindexSortTd[chd];
                    n = (nnodesd - 1) - srtT;
                    ps = Mposd[chd];
                    mom = momsd + (srtT * order);
                    gm = mom[0].x;
                    sumSide2 = Mupperd[chd] - Mlowerd[chd]; /*Msized[chd]*/;
                }

                dr = p - ps;
                r2 = (dr.x * dr.x + dr.y * dr.y);   // compute distance squared

                // check if all threads agree that cell is far enough away (or is a body)
                if (isVortex || __all_sync(0xffffffff, ((sumSide2.x + sumSide2.y) * (sumSide2.x + sumSide2.y) + epssqd) * itolsqd < r2))
                {
                    if (calcEpsAst)
                    {
                        if ((r2 < d_3) && (r2 > 0))
                        {
                            dst23 = realmin(r2, d_2);
                            d_3 = realmax(r2, d_2);

                            dst12 = realmin(dst23, d_1);
                            d_2 = realmax(dst23, d_1);

                            d_1 = dst12;
                        }
                    }

#ifdef CALCinDOUBLE
                    real f = gm / realmax(r2, epssqd);
#else
                    real f = fdividef(gm, realmax(r2, epssqd));
#endif
                    v += f * dr;

                    if ((!isVortex) && (order > 1))
                    {
#ifdef CALCinDOUBLE
                        real2 cftr = (r2 ? (1 / r2) : (real)0) * dr;
#else
                        real2 cftr = (r2 ? fdividef(idpi, r2) : 0.0f) * dr;
#endif
                        th = cftr;

                        for (int s = 1; s < order; ++s)
                        {
                            th = multz(th, cftr);
#ifdef CALCinFLOAT                                    
                            if (isinf(th.x) || isinf(th.y))
                            {
                                //printf("s = %d\n", s);
                                break;
                            }
#endif
                            v += multzA(th, mom[s]);
                        }
                    }
                }
                else
                {
                    // push cell onto stack
                    if (sbase == threadIdx.x)
                    {  // maybe don't push and inc if last child
                        pos[depth] = pd;
                        node[depth] = nd;
                    }
                    depth++;
                    pd = 0;
                    nd = n;

                    chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];
                }

            }
            depth--;  // done with this level
        } while (depth >= j);

        if (calcEpsAst)
        {
            for (size_t afl = 0; afl < nAfls; ++afl)
                for (size_t j = 0; j < nVtxs[afl]; ++j)
                {
                    dr.x = p.x - ptrVtxs[afl][(j) * 3 + 0 + 0];
                    dr.y = p.y - ptrVtxs[afl][(j) * 3 + 0 + 1];
                    r2 = dr.x * dr.x + dr.y * dr.y;


                    if ((d_3 > r2) && (r2 > 0))
                    {
                        dst23 = realmin(r2, d_2);
                        d_3 = realmax(r2, d_2);

                        dst12 = realmin(dst23, d_1);
                        d_2 = realmax(dst23, d_1);

                        d_1 = dst12;
                    }
                }
        }

        // update velocity

        real2 result = real2{ -idpi * v.y, idpi * v.x };
        veld[indexOfPoint] = result;

        if (calcEpsAst)
            epsast[indexOfPoint] = sqrt((d_1 + d_2 + d_3) / 3);
    }
}


__global__
__launch_bounds__(THREADS5diff, FACTOR5diff)
void DiffusiveForceCalculationKernel2(
    const int nnodesd, const int nbodiesd,
    const real epssqd,
    const int2* __restrict Mchildd,
    const int order, const real2* __restrict momsd,
    const real3* __restrict vtxd,
    const int* __restrict MmortonCodesIdxd,
    const real2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
    real* __restrict I1d,
    real2* __restrict I2d,
    const real2* __restrict Mlowerd,
    const real2* __restrict Mupperd,
    const real* __restrict epsast)

{
    register int j, k, n, depth, base, sbase, pd, nd;
    register real2 p, i2, dr, ps;
    register real r2sq, i1;
    register const real2* mom;

    register real rdi, diffRadius, diffRadiusMonopole;
    register real expr;



    __shared__ volatile int pos[MAXDEPTH * THREADS5diff / WARPSIZE], node[MAXDEPTH * THREADS5diff / WARPSIZE];



    // figure out first thread in each warp (lane 0)
    base = threadIdx.x / WARPSIZE;
    sbase = base * WARPSIZE;
    j = base * MAXDEPTH;
    //diff = threadIdx.x - sbase;

    __syncthreads();
    __threadfence_block();

    // iterate over all bodies assigned to thread
    for (k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd; k += blockDim.x * gridDim.x)
    {
        
        const int indexInParticles = MmortonCodesIdxd[k];

        p = real2{ vtxd[indexInParticles].x, vtxd[indexInParticles].y };
       
        i1 = 0;
        i2.x = 0;
        i2.y = 0;

        // initialize iteration stack, i.e., push root node onto stack
        depth = j;
        if (sbase == threadIdx.x)
        {
            pos[j] = 0;
            node[j] = nnodesd - 1;
        }

        do
        {            
            pd = pos[depth];
            nd = node[depth];

            register int2 chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];

            register real gm;
            register real2 sumSide2;
            bool isVortex;

            while (pd < 2)
            {
                // node on top of stack has more children to process

                // load child pointer
                //computation of n = childd[nd + pd] (pd = 0 или pd = 1)
                int chd = pd * chBoth.y + (1 - pd) * chBoth.x;
                ++pd;

                isVortex = (chd >= nbodiesd);

                if (isVortex)
                {
                    n = chd - nbodiesd;

                    ps = real2{ vtxd[MmortonCodesIdxd[n]].x, vtxd[MmortonCodesIdxd[n]].y };
                    gm = vtxd[MmortonCodesIdxd[n]].z;
                    sumSide2 = real2{ (real)0, (real)0 };
                }
                else
                {
                    register const int srtT = MindexSortTd[chd];
                    n = (nnodesd - 1) - srtT;

                    ps = Mposd[chd];
                    mom = momsd + (srtT * order);
                    gm = mom[0].x;

                    sumSide2 = Mupperd[chd] - Mlowerd[chd]; /*Msized[chd]*/;
                }

                dr = p - ps;
                r2sq = sqrt(dr.x * dr.x + dr.y * dr.y);   // compute distance 
                
                rdi = max(epsast[indexInParticles], sqrt(epssqd));
                diffRadius = 8.0 * rdi;
                diffRadiusMonopole = 4.0 * rdi;

                

                // check if all threads agree that cell is far enough away (or is a body)
                //if (isVortex || __all_sync(0xffffffff, ((sumSide2.x + sumSide2.y) * (sumSide2.x + sumSide2.y) + epssqd) * itolsqd < r2))
                if (__all_sync(0xffffffff, r2sq - 0.5 * (sumSide2.x + sumSide2.y) > diffRadius))
                {

                }
                else if (isVortex || __all_sync(0xffffffff, r2sq - 0.5*(sumSide2.x + sumSide2.y) > diffRadiusMonopole))
                {
                    if ((r2sq < diffRadius) && (r2sq > 1e-10))
                    {
                        expr = gm * exp(-r2sq / rdi);
                        i1 += expr;
                        i2 += expr / r2sq * dr;
                    }
                }
                else
                {
                    // push cell onto stack
                    if (sbase == threadIdx.x)
                    {  // maybe don't push and inc if last child
                        pos[depth] = pd;
                        node[depth] = nd;
                    }
                    depth++;
                    pd = 0;
                    nd = n;

                    chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];
                }

            }
            depth--;  // done with this level
        } while (depth >= j);
               

        // update velocity
        I1d[indexInParticles] = i1;
        I2d[indexInParticles] = i2; 
    }
}//DiffusiveVelo






__global__
void VerticesToControlPointsKernel(int nTotPan, real* dev_ptr_pt, real* pointsl)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nTotPan)
    {
        pointsl[3 * i + 0] = 0.5 * (dev_ptr_pt[4 * i + 0] + dev_ptr_pt[4 * i + 2]);
        pointsl[3 * i + 1] = 0.5 * (dev_ptr_pt[4 * i + 1] + dev_ptr_pt[4 * i + 3]);
        pointsl[3 * i + 2] = 0.0;
    }
}



__global__
__launch_bounds__(THREADS5rhs, FACTOR5rhs)
void RhsCalculationKernel
(
    const int nnodesd, const int nbodiesd,
    const real itolsqd, 
    const int2* __restrict Mchildd,
    const int order, const real2* __restrict momsd,
    const real3* __restrict vtxd,
    const int* __restrict MmortonCodesIdxd,
    const real2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
    const int npointsd, const real* __restrict dev_ptr_pt, const real3* __restrict pointsd, real2* __restrict Ed,
    const int* __restrict MmortonCodesIdxPointsd,
    real* __restrict veld,
    real* __restrict vellind,
    const real2* __restrict Mlowerd,
    const real2* __restrict Mupperd)
{
    register int j, k, n, depth, base, sbase, pd, nd;
    register real2 p, dr, ps;
    register real r2;
    register const real2* mom;
    register real val, vallin;

    __shared__ volatile int pos[MAXDEPTH * THREADS5rhs / WARPSIZE], node[MAXDEPTH * THREADS5rhs / WARPSIZE];

    // figure out first thread in each warp (lane 0)
    base = threadIdx.x / WARPSIZE;
    sbase = base * WARPSIZE;
    j = base * MAXDEPTH;

    __syncthreads();
    __threadfence_block();

    // iterate over all bodies assigned to thread
    for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npointsd; k += blockDim.x * gridDim.x)
    {
        const int indexOfPoint = MmortonCodesIdxPointsd[k];

        p = real2{ pointsd[indexOfPoint].x, pointsd[indexOfPoint].y };

        real2 beg, end;
        beg = real2{ dev_ptr_pt[4 * indexOfPoint + 0], dev_ptr_pt[4 * indexOfPoint + 1] };
        end = real2{ dev_ptr_pt[4 * indexOfPoint + 2], dev_ptr_pt[4 * indexOfPoint + 3] };

        real dlen2 = (end.x - beg.x) * (end.x - beg.x) + (end.y - beg.y) * (end.y - beg.y);
        real idlen = 1 / sqrt(dlen2);
        real2 tau = (end - beg) * idlen;

        val = 0;
        vallin = 0;

        // initialize iteration stack, i.e., push root node onto stack
        depth = j;
        if (sbase == threadIdx.x)
        {
            pos[j] = 0;
            node[j] = nnodesd - 1;
        }

        do
        {
            // stack is not empty
            pd = pos[depth];
            nd = node[depth];

            register int2 chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];

            register real gm;
            register real2 sumSide2;
            bool isVortex;

            while (pd < 2)
            {
                // node on top of stack has more children to process

                // load child pointer
                //computation of n = childd[nd + pd] (pd = 0 или pd = 1)
                int chd = pd * chBoth.y + (1 - pd) * chBoth.x;
                ++pd;

                isVortex = (chd >= nbodiesd);

                if (isVortex)
                {
                    n = chd - nbodiesd;
                   
                    ps = real2{ vtxd[MmortonCodesIdxd[n]].x, vtxd[MmortonCodesIdxd[n]].y };
                    gm = vtxd[MmortonCodesIdxd[n]].z;
                    sumSide2 = real2{ (real)0, (real)0 };
                }
                else
                {
                    register const int srtT = MindexSortTd[chd];
                    n = (nnodesd - 1) - srtT;
                    ps = Mposd[chd];
                    mom = momsd + (srtT * order);
                    gm = mom[0].x;
                    sumSide2 = Mupperd[chd] - Mlowerd[chd]; /*Msized[chd]*/;
                }

                //ps - положение вихря/кластера
                //p - центр панели
                dr = p - ps;
                r2 = (dr.x * dr.x + dr.y * dr.y);   // compute distance squared               

                // check if all threads agree that cell is far enough away (or is a body)
                if (isVortex || __all_sync(0xffffffff, ((sumSide2.x + sumSide2.y) * (sumSide2.x + sumSide2.y) + dlen2) * itolsqd < r2))
                {
                    if (isVortex)
                    {
                        real2 s = ps - beg;
                        real2 p = ps - end;

                        real alpha = atan2(p.x * s.y - p.y * s.x, p.x * s.x + p.y * s.y);

                        real tempVel = gm * alpha;
                        val -= tempVel;

                        if (vellind != nullptr)
                        {
                            
                            real2 u1;
                            u1.x = 0.5 * idlen * ((p.x + s.x) * tau.x * tau.x \
                                + 2.0 * (p.y + s.y) * tau.x * tau.y - (p.x + s.x) * tau.y * tau.y);
                            u1.y = 0.5 * idlen * (-(p.y + s.y) * tau.x * tau.x \
                                + 2.0 * (p.x + s.x) * tau.x * tau.y + (p.y + s.y) * tau.y * tau.y);

                            real lambda = 0.5 * log((s.x * s.x + s.y * s.y) / (p.x * p.x + p.y * p.y));

                            real tempVelLin = gm * (alpha * (u1.x * tau.x + u1.y * tau.y) + lambda * (-u1.y * tau.x + u1.x * tau.y));

                            vallin -= tempVelLin;
                           
                        }
                    }
                    else
                    {
                        real dist2 = (p.x - ps.x) * (p.x - ps.x) + (p.y - ps.y) * (p.y - ps.y);

                        // Это просто смещение указателя
                        real2* Eloc = Ed + order * indexOfPoint;

                        for (int q = 0; q < order; ++q)
                            Eloc[q].x = Eloc[q].y = 0;

                        real2 theta = (p - ps) / dist2;

                        int cftStart = 1;
                        int cftDiag;
                        real factorial;

                        for (int q = 0; q < order; ++q)
                        {
                            factorial = 1;
                            cftDiag = 1;
                            for (int s = q; s >= 0; --s)
                            {
                                Eloc[s] += (cftStart * cftDiag / factorial)  * multzA(theta, mom[q - s]);
                                
                                cftDiag = -cftDiag;
                                factorial *= (q - s + 1);
                            }
                            theta = ((q + 1) / dist2) * multz(theta, p - ps);
                            cftStart = -cftStart;
                        }



                        real2 v = Eloc[0];
                        real2 vL{0.0, 0.0};

                        real2 rPan = end - beg;
                        real panLen2 = (rPan.x * rPan.x + rPan.y * rPan.y);
                        real2 dPos2 = real2{ 0.0, 0.0 }; // т.к. разложение строится в центре панели

                        real2 kp, km;

                        real2 mulP = kp = dPos2 + rPan;
                        real2 mulM = km = dPos2 - rPan;

                        real2 taudL = (0.5 / panLen2) * rPan;
                        real2 taudLc = taudL;
                        
                        real iFact = 1;
                        for (int k = 1; k < order; ++k)
                        {
                            iFact /= (k+1);
                            mulP = multz(mulP, kp);
                            mulM = multz(mulM, km);
                            taudL = taudL * 0.5;
                            //if (fabs(Eloc[k].x) + fabs(Eloc[k].y) < 1e+10)
                            v += iFact * multz(Eloc[k], multzA(taudL, mulP - mulM));

                            
                            if (vellind != nullptr)
                                vL += (iFact / (k + 2)) * multz(Eloc[k], multzA(multz(taudL, taudLc), multz(mulP, (k + 1) * rPan - dPos2) + multz(mulM, (k + 1) * rPan + dPos2)));
                        }
                        
                        val += (-v.y * rPan.x + v.x * rPan.y);

                        
                        if (vellind != nullptr)
                            vallin += (-vL.y * rPan.x + vL.x * rPan.y);    
                        
                    }                    
                }
                else //идем глубже по дереву
                {
                    // push cell onto stack
                    if (sbase == threadIdx.x)
                    {  // maybe don't push and inc if last child
                        if (depth >= MAXDEPTH * THREADS5rhs / WARPSIZE)
                            printf("?????????????????????????????????????????????????????????????????\n");

                        pos[depth] = pd;
                        node[depth] = nd;
                    }
                    depth++;
                    pd = 0;
                    nd = n;

                    chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];
                }

            }
            depth--;  // done with this level
        } while (depth >= j);


        // update velocity
        veld[indexOfPoint] = val * idpi * idlen;
        
        if (vellind != nullptr)
            vellind[indexOfPoint] = vallin* idpi* idlen;
        
    }
}

/******************************************************************************/
/*** compute force (direct) ***************************************************/
/******************************************************************************/


__global__
//__launch_bounds__(THREADSD, FACTORD)
void ForceDirectCalculationKernel(
    const int nbodiesd,
    const real epssqd,
    const real3* __restrict vtxd,    
    real2* __restrict veld)
{
    __shared__ real3 shvtx[BLOCKD];    

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    real2 pt;
    pt.x = vtxd[i].x;
    pt.y = vtxd[i].y;

    real2 vel;
    vel.x = vel.y = 0;

    real2 dr;
    real dr2, izn;

    //vortices
    for (size_t j = 0; j < nbodiesd; j += BLOCKD)
    {
        shvtx[threadIdx.x] = vtxd[j + threadIdx.x];               

        __syncthreads();

        for (size_t q = 0; q < BLOCKD; ++q)
        {
            if (j + q < nbodiesd)
            {
                dr.x = pt.x - shvtx[q].x;
                dr.y = pt.y - shvtx[q].y;
                dr2 = dr.x * dr.x + dr.y * dr.y;

                izn = shvtx[q].z / realmax(dr2, epssqd);// / CUboundDenom(dr2, eps2); //РЎРіР»Р°Р¶РёРІР°С‚СЊ РЅР°РґРѕ!!!

                vel.x -= dr.y * izn;
                vel.y += dr.x * izn;

            }
        }
        __syncthreads();
    }

    if (i < nbodiesd)
    {
        veld[i].x = vel.x;// * iDPI;
        veld[i].y = vel.y;// * iDPI;
    }
    //*/
}


/******************************************************************************/


void KernelsOptimization()
{
    // set L1/shared memory configuration

    cudaFuncSetCacheConfig(ClearKernel2, cudaFuncCachePreferL1);
    cudaGetLastError();  // reset error value
        
    cudaFuncSetCacheConfig(ForceDirectCalculationKernel, cudaFuncCachePreferEqual); //d
    cudaGetLastError();  // reset error value
        
    cudaFuncSetCacheConfig(SummarizationKernel2_14, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_16, cudaFuncCachePreferL1);
    cudaGetLastError();  // reset error value

    cudaFuncSetCacheConfig(ForceCalculationKernel2points, cudaFuncCachePreferL1); //d
    cudaGetLastError();  // reset error value

}


/******************************************************************************/





//////////////////
/// Wrappers
//////////////////



    /******************************************************************************/
    /*** initialize memory ********************************************************/
    /******************************************************************************/

    float cuInitializationKernel()
    {
        //fprintf(stderr, "IKKernel\n");
        
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);

        cudaEventRecord(start, 0);
        InitializationKernel<<<1, 1>>> ();
        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);
        CudaTest("kernel 0 launch failed");
        
        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }


    /******************************************************************************/
    /*** compute center and radius ************************************************/
    /******************************************************************************/

    float McuBoundingBoxKernelFree(
        realPoint* __restrict Mposl,
        realPoint* __restrict maxrl,
        realPoint* __restrict minrl,
        int nbodiesd,
        const realVortex* __restrict vtxd)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        MBoundingBoxKernel << <blocks * FACTOR1, THREADS1 >> > (nbodiesd, (real3*)vtxd, (real2*)Mposl, (real2*)maxrl, (real2*)minrl);
        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("Mkernel 1 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


	float McuBoundingBoxKernel(
        CUDApointers ptr,
		int nbodiesd,
		const realVortex* __restrict vtxd)
	{
        float time;
        time = McuBoundingBoxKernelFree(ptr.Mposl, ptr.maxrl, ptr.minrl, nbodiesd, vtxd);
		return time;
	}

	/******************************************************************************/
	/*** Morton codes *************************************************************/
	/******************************************************************************/
    float McuMortonCodesKernelFree(
        const realPoint* __restrict maxrl,
        const realPoint* __restrict minrl,
        int* __restrict MmortonCodesKeyUnsortl, 
        int* __restrict MmortonCodesIdxUnsortl, 
        int* __restrict MmortonCodesKeyl, 
        int* __restrict MmortonCodesIdxl,
        const intPair* __restrict Mrangel,

        int nbodiesd,
        const realVortex* __restrict vtxd)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = (nbodiesd + 31) / 32;
        dim3 Mthreads = 32;
            

        MMortonCodesKernel << <Mblocks, Mthreads >> > (nbodiesd, (const real3*)vtxd, MmortonCodesKeyUnsortl, MmortonCodesIdxUnsortl, (const real2*)maxrl, (const real2*)minrl);

        ///RadixSort

        RadixSortFromCUB(
            MmortonCodesKeyUnsortl, MmortonCodesKeyl, \
            MmortonCodesIdxUnsortl, MmortonCodesIdxl, \
            nbodiesd, 0, 2 * codeLength);

        //Заполнение нулевой ячейки (диапазон для корня дерева)
        int totalRange[2] = { 0, nbodiesd - 1 };
        if (Mrangel != nullptr)
            cudaCopyVecToDevice(totalRange, (void*)Mrangel, 2, sizeof(int));

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("Mkernel 1 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    };


    float McuMortonCodesKernel(
        CUDApointers ptr,
        int nbodiesd,
        const realVortex* __restrict vtxd)
    {
        float time;
        time = McuMortonCodesKernelFree(
            ptr.maxrl, ptr.minrl, 
            ptr.MmortonCodesKeyUnsortl, ptr.MmortonCodesIdxUnsortl, ptr.MmortonCodesKeyl, ptr.MmortonCodesIdxl, 
            ptr.Mrangel, 
            nbodiesd, vtxd);
        
		return time;
	}

    /******************************************************************************/
    /*** Morton Internal nodes build **********************************************/
    /******************************************************************************/

    float McuMortonInternalNodesKernel(
        CUDApointers ptr,
        int nbodiesd)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = ((nbodiesd - 1) + 31) / 32;
        dim3 Mthreads = 32;

        MMortonInternalNodesKernel<<<Mblocks, Mthreads>>> (nbodiesd, ptr.MmortonCodesKeyl, ptr.Mparentl, (int2*)ptr.Mchildl, (int2*)ptr.Mrangel);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("Mkernel 2 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;

    }



    /******************************************************************************/
    /*** Morton Internal nodes geometry calculation *******************************/
    /******************************************************************************/
    float McuMortonInternalCellsGeometryKernel(
        CUDApointers ptr,
        int nbodiesd,
        int nnodesd)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = ((nbodiesd) + 31) / 32;
        dim3 Mthreads = 32;

		MMortonInternalCellsGeometryKernel<<<Mblocks, Mthreads>>>(nbodiesd, ptr.MmortonCodesKeyl, /*(real2*)ptr.Mposl, (real2*)ptr.Msizel, */ (int2*)ptr.Mrangel, ptr.MlevelUnsortl, ptr.MindexUnsortl, (real2*)ptr.maxrl, (real2*)ptr.minrl);

        RadixSortFromCUB( \
            ptr.MlevelUnsortl, ptr.MlevelSortl, \
            ptr.MindexUnsortl, ptr.MindexSortl, \
            nbodiesd-1, 0, 2 * codeLength);

        MTransposeIndexKernel << <Mblocks, Mthreads >> > (nbodiesd, nnodesd, ptr.MindexSortl, ptr.MindexSortTl);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("Mkernel 3 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }



	 	 

    /******************************************************************************/
    /*** build tree ***************************************************************/
    /******************************************************************************/   

    float cuClearKernel2(
        CUDApointers ptr,
        const int order,
        int nnodesd, int nbodiesd)
    {
        //fprintf(stderr, "CxKernel\n");
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        //printf("MEMSET: Trying to fill %d bytes for %d bodies, order = %d, sizeof = %d\n", (int)((nbodiesd - 1) * order * sizeof(realPoint)), nbodiesd - 1, order, sizeof(realPoint));

        cudaMemset((void*)ptr.momsl, 0, (nbodiesd-1) * order * sizeof(realPoint));

        //CudaTest("Memset in kernel clear2 launch failed");

        //printf("nnodesd = %d, nbodiesd = %d\n", nnodesd, nbodiesd);

        ClearKernel2 << <blocks * 1, 1024 >> > (nnodesd, nbodiesd, ptr.massl);
        
	       	

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel clear2 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }





    float cuAABBKernel2(
        CUDApointers ptr,
        const int nnodesd, const int nbodiesd,
        const realVortex* __restrict vtxd)
    {
        //fprintf(stderr, "AABBKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        const intPair* __restrict Mchildd = ptr.Mchildl;
        volatile int* massd = ptr.massl;
        const int* __restrict MmortonCodesIdxd = ptr.MmortonCodesIdxl;
        realPoint* __restrict Mposd = ptr.Mposl;
        realPoint* __restrict Mlowerd = ptr.Mlowerl;
        realPoint* __restrict Mupperd = ptr.Mupperl;
        
        const int* __restrict MindexSortd = ptr.MindexSortl;
        const int* __restrict MindexSortTd = ptr.MindexSortTl;

        AABBKernel2 << <blocks * FACTOR3, THREADS3 >> > (nnodesd, nbodiesd, (int2*)Mchildd, massd, (real3*)vtxd, MmortonCodesIdxd, (real2*)Mposd, (real2*)Mlowerd, (real2*)Mupperd, MindexSortd, MindexSortTd);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 3a launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }



    /******************************************************************************/
    /*** compute multipole moments for all the cells ******************************/
    /******************************************************************************/
    float cuSummarizationKernel2(
        CUDApointers ptr,
        const int order,
        const int nnodesd, const int nbodiesd,
        const realVortex* __restrict vtxd)
    {
        //fprintf(stderr, "SKKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);
        
        const intPair* __restrict Mchildd = ptr.Mchildl;
        volatile int* massd = ptr.massl;
        const realPoint* __restrict momsd = ptr.momsl;
        const int* __restrict MmortonCodesIdxd = ptr.MmortonCodesIdxl;
        const realPoint* __restrict Mposd = ptr.Mposl;
        const int* __restrict MindexSortd = ptr.MindexSortl; 
        const int* __restrict MindexSortTd = ptr.MindexSortTl;

#include "ShiftKernels/SwitchKer.cu"

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 3 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }

   

    /******************************************************************************/
    /*** compute force ************************************************************/
    /******************************************************************************/
    /*
    float cuForceCalculationKernel2(
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real itolsqd, real epssqd,
        const realVortex* __restrict vtxd,
        realPoint* __restrict veld,
        bool calcEpsAst, real* __restrict epsastd,
        size_t nAfls, size_t* nVtxs, double** ptrVtxs)
    {
        //fprintf(stderr, "FCKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        ForceCalculationKernel2 << <blocks * FACTOR5, THREADS5 >> > (
            nnodesd, nbodiesd, itolsqd, epssqd, (int2*)ptr.Mchildl, order, (real2*)ptr.momsl,
            (real3*)vtxd, ptr.MmortonCodesIdxl,
            (real2*)ptr.Mposl, ptr.MindexSortl, ptr.MindexSortTl,
            (real2*)veld, 
            (real2*)ptr.Mlowerl, (real2*)ptr.Mupperl,
            calcEpsAst, epsastd,
            nAfls, nVtxs, ptrVtxs);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 5 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }
    */

    /******************************************************************************/
    /*** compute diffusive velo****************************************************/
    /******************************************************************************/
    float cuDiffVelCalculationKernel2(
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real itolsqd, real epssqd,
        const realVortex* __restrict vtxd,
        //realPoint* __restrict veld,
        real* __restrict I1d,
        realPoint* __restrict I2d,
        bool calcEpsAst, real* __restrict epsastd,
        size_t nAfls, size_t* nVtxs, double** ptrVtxs)
    {
        //fprintf(stderr, "FCKernel\n");

        //cudaEvent_t start, stop;
        float time = 0;

        //cudaEventCreate(&start);  cudaEventCreate(&stop);
        //cudaEventRecord(start, 0);

        DiffusiveForceCalculationKernel2 << <blocks * FACTOR5diff, THREADS5diff >> > (
            nnodesd, nbodiesd, epssqd, (int2*)ptr.Mchildl, order, (real2*)ptr.momsl,
            (real3*)vtxd, ptr.MmortonCodesIdxl,
            (real2*)ptr.Mposl, ptr.MindexSortl, ptr.MindexSortTl,
             (real*)I1d,
            (real2*)I2d,
            (real2*)ptr.Mlowerl, (real2*)ptr.Mupperl,
            epsastd);


        //cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 5 diffvelo launch failed");

        //cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


     /******************************************************************************/
    /*** compute force ************************************************************/
    /******************************************************************************/
    float cuForceCalculationKernel2points(
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real itolsqd, real epssqd,
        const realVortex* __restrict vtxd,
        const int* __restrict MmortonCodesIdxl,
        int npointsd, const realVortex* pointsd,
        realPoint* __restrict veld,
        bool calcEpsAst, real* __restrict epsastd,
        size_t nAfls, size_t* nVtxs, double** ptrVtxs)
    {
        //fprintf(stderr, "FCKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        ForceCalculationKernel2points << <blocks * FACTOR5, THREADS5rhs >> > (
            nnodesd, nbodiesd, itolsqd, epssqd, (int2*)ptr.Mchildl, order, (real2*)ptr.momsl,
            (real3*)vtxd, ptr.MmortonCodesIdxl,
            (real2*)ptr.Mposl, ptr.MindexSortl, ptr.MindexSortTl,
            npointsd, (real3*)pointsd, MmortonCodesIdxl,
            (real2*)veld, 
            (real2*)ptr.Mlowerl, (real2*)ptr.Mupperl,
            calcEpsAst, epsastd,
            nAfls, nVtxs, ptrVtxs);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 5 convvelo launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }



    float McuVerticesToControlPoints(int nTotPan, double* dev_ptr_pt, double* pointsl)
    {
        cudaEvent_t startD, stopD;
        float timeD;

        cudaEventCreate(&startD);  cudaEventCreate(&stopD);
        cudaEventRecord(startD, 0);

        VerticesToControlPointsKernel <<<nTotPan, BLOCKD >> > (nTotPan, dev_ptr_pt, pointsl);
        cudaEventRecord(stopD, 0);  cudaEventSynchronize(stopD);  cudaEventElapsedTime(&timeD, startD, stopD);

        CudaTest("kernel direct launch failed");

        cudaEventDestroy(startD);  cudaEventDestroy(stopD);

        return timeD;
    }



    float cuRhsCalculationKernel
    (
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real itolsqd, 
        const realVortex* __restrict vtxd,
        //CUDApointers ptrPoints,
        const int* __restrict MmortonCodesIdxl,
        realPoint* __restrict El,
        int nTotPan, const real* dev_ptr_pt, const real* pointsd,
        real* __restrict veld,
        real* __restrict vellind)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        RhsCalculationKernel << <blocks * FACTOR5rhs, THREADS5rhs >> > (
            nnodesd, nbodiesd, itolsqd, 
            (int2*)ptr.Mchildl, order, (real2*)ptr.momsl,
            (real3*)vtxd, ptr.MmortonCodesIdxl,
            (real2*)ptr.Mposl, ptr.MindexSortl, ptr.MindexSortTl,
            nTotPan, (const real*)dev_ptr_pt, (const real3*)pointsd, (real2*)El, MmortonCodesIdxl,
            (real*)veld, (real*)vellind,
            (real2*)ptr.Mlowerl, (real2*)ptr.Mupperl);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 5 rhs launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


    /******************************************************************************/
    /*** compute force (direct) ***************************************************/
    /******************************************************************************/

    float cuForceDirectCalculationKernel(
        int nbodiesd,
        real epssqd,
        const realVortex* __restrict vtxd,
        volatile realPoint* __restrict veld)
    {
        //fprintf(stderr, "DFKernel\n");
        
        cudaEvent_t startD, stopD;
        float timeD;

        cudaEventCreate(&startD);  cudaEventCreate(&stopD);
        cudaEventRecord(startD, 0);
        
        ForceDirectCalculationKernel<<<(nbodiesd + BLOCKD - 1) / BLOCKD, BLOCKD>>> (nbodiesd, epssqd, (real3*)vtxd, (real2*)veld);
        cudaEventRecord(stopD, 0);  cudaEventSynchronize(stopD);  cudaEventElapsedTime(&timeD, startD, stopD);
        
        CudaTest("kernel direct launch failed");

        cudaEventDestroy(startD);  cudaEventDestroy(stopD);
        
        return timeD;
    }

}//namespace BHcu