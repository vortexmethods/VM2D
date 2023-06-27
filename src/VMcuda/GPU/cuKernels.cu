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
#define pi (3.1415926535897932384626433832795)
#define pif (3.1415926535897932384626433832795f)


namespace BHcu
{
    __device__ double myAtomicAdd(double* address, double val)






    {
        unsigned long long int* address_as_ull = (unsigned long long int*)address;
        unsigned long long int old = *address_as_ull, assumed;

        do {
            assumed = old;
            old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));

            // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
        } while (assumed != old);
        return __longlong_as_double(old);
    }

    __device__ double myAtomicExch(double* address, double val)
    {
        unsigned long long int* address_as_ull = (unsigned long long int*)address;
        unsigned long long int old = *address_as_ull, assumed;

        do {
            assumed = old;
            old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val));

            // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
        } while (assumed != old);
        return __longlong_as_double(old);
    }

    __device__ inline double myMax(double x, double y)
    {
        return (x > y) ? x : y;
    }

    __device__ inline double myMin(double x, double y)
    {
        return (x > y) ? y : x;
    }

    __device__ inline float myMax(float x, float y)
    {
        return (x > y) ? x : y;
    }

    __device__ inline float myMin(float x, float y)
    {
        return (x > y) ? y : x;
    }


    __device__ inline int sqr(int x)
    {
        return x * x;
    }

    __device__ inline double sqr(double x)
    {
        return x * x;
    }

    __device__ inline double sqrf(float x)
    {
        return x * x;
    }

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


void cudaDelete(void* cudaPtr, int code)
{
	cudaError_t err1 = cudaFree(cudaPtr);

    if (err1 != cudaSuccess)
        std::cout << cudaGetErrorString(err1) << " (cudaDelete), code = " << code << std::endl;
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
    const void* __restrict ptrd, 
    int sizeOfElement, 
    int offsetOfPointInElement,
    real2* __restrict Mposd, 
    volatile real2* __restrict maxrd, 
    volatile real2* __restrict minrd)
{
    register int i, j, k, inc;
    register real2 val;
    register real2 minr, maxr;
    __shared__ volatile real2 sminr[THREADS1], smaxr[THREADS1];



    // initialize with valid data (in case #bodies < #threads)
    real* rootPtrX = (real*)((char*)ptrd + 0 * sizeOfElement + offsetOfPointInElement);
    
    minr.x = maxr.x = *(rootPtrX+0); // vtxd[0].x;
    minr.y = maxr.y = *(rootPtrX+1); // vtxd[0].y;    
    
    // scan all bodies
    i = threadIdx.x;
    inc = THREADS1 * gridDim.x;
    
    for (j = i + blockIdx.x * THREADS1; j < nbodiesd; j += inc) 
    {
        real* ptrValX = (real*)((char*)ptrd + j * sizeOfElement + offsetOfPointInElement);
        
        val.x = *(ptrValX+0);
        val.y = *(ptrValX+1);
        //val.x = vtxd[j].x;
        //val.y = vtxd[j].y;

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
    const void* __restrict ptrl, 
    int sizeOfElement, 
    int offsetOfPointInElement,
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

        real* ptrX = (real*)((char*)ptrl + bdy * sizeOfElement + offsetOfPointInElement);

        real x0 = (*(ptrX + 0) - lx0) * quadSideFactor + rbound / 2;
        real y0 = (*(ptrX + 1) - ly0) * quadSideFactor + rbound / 2;
        //real x0 = (vtxd[bdy].x - lx0) * quadSideFactor + rbound/2;
        //real y0 = (vtxd[bdy].y - ly0) * quadSideFactor + rbound/2;

       
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
    

    //bool cond = ( (blockIdx.x == 113) && (threadIdx.x == 0) );

    if (i < nbodies - 1)
    {
        //if (cond)
        //{
        //    printf("MmortonCodesKeyd = %p\n", MmortonCodesKeyd);
        //    printf("Mparentd = %p\n", Mparentd);
        //    printf("Mchildd = %p\n", Mchildd);
        //    printf("Mranged = %p\n", Mranged);
        //}
        
        //if (cond)
        //    printf("i = %d\n", i);
        
        int d = sign(Delta(i, i + 1, nbodies, MmortonCodesKeyd) - Delta(i, i - 1, nbodies, MmortonCodesKeyd));

        //if (cond)
        //    printf("d = %d\n", d);

        int delta_min = Delta(i, i - d, nbodies, MmortonCodesKeyd);

        //if (cond)
        //    printf("delta_min = %d\n", delta_min);

        int Lmax = 2;
        while (Delta(i, i + Lmax * d, nbodies, MmortonCodesKeyd) > delta_min)
            Lmax *= 2;

        //if (cond)
        //    printf("Lmax = %d\n", Lmax);

        int L = 0;
        for (int t = (Lmax >> 1); t >= 1; t >>= 1)
            if (Delta(i, i + (L + t) * d, nbodies, MmortonCodesKeyd) > delta_min)
                L += t;

        //if (cond)
        //    printf("L = %d\n", L);

        int j = i + L * d;

        int delta_node = Delta(i, j, nbodies, MmortonCodesKeyd);

        //if (cond)
        //    printf("delta_node = %d\n", delta_node);

        int s = 0;
        for (int p = 1, t = ceilhalf(L); L > (1 << (p - 1)); ++p, t = ceilpow2(L, p))
        {
            int dl = Delta(i, i + (s + t) * d, nbodies, MmortonCodesKeyd);
            if (dl > delta_node)
                s += t;
        }//for p

        //if (cond)
        //    printf("s = %d\n", s);


        int gamma = i + s * d +   d * (d < 0);   //последнее слагаемое = std::min(d, 0);

        int Mmin = min(i, j);
        int Mmax = max(i, j);
        
        const int& left = gamma;
        const int& right = gamma + 1;

        // Левый потомок - лист или внутренний узел
        int childLeft = Mchildd[i].x = (Mmin == gamma) * nbodies + left;
        
        //if (cond)
        //    printf("childLeft = %d\n", childLeft);

        Mranged[childLeft].x = Mmin;
        Mranged[childLeft].y = gamma;
        Mparentd[childLeft] = i;

        // Правый потомок - лист или внутренний узел
        int childRight = Mchildd[i].y = (Mmax == gamma + 1) * nbodies + right;

        //if (cond)
        //    printf("childRight = %d\n", childRight);

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
    
    if (cell < nbodiesd - 1)
    {
        register const int newcell = MindexSortd[cell];
        MindexSortTd[newcell] = cell;
    }

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
    /*const real3* __restrict vtxd,*/
    const void* __restrict ptrd,
    int sizeOfElement,
    int offsetOfPointInElement,
    const int* __restrict MmortonCodesIdxd,
    real2* __restrict Mposd,
    real2* __restrict Mlowerd,
    real2* __restrict Mupperd,
    const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
    bool bodiesZeroSize,
    const double* XYXY
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
                        real* PtrX = (real*)((char*)ptrd + sortedBody * sizeOfElement + offsetOfPointInElement);

                        /*lower[i] = upper[i] = real2{vtxd[sortedBody].x, vtxd[sortedBody].y};*/
                        
                        if (bodiesZeroSize)
                            lower[i] = upper[i] = real2{ *(PtrX + 0), *(PtrX + 1) };
                        else
                        {
                            lower[i] = real2{ ::fmin(XYXY[4 * sortedBody + 0], XYXY[4 * sortedBody + 2]), ::fmin(XYXY[4 * sortedBody + 1], XYXY[4 * sortedBody + 3]) };
                            upper[i] = real2{ ::fmax(XYXY[4 * sortedBody + 0], XYXY[4 * sortedBody + 2]), ::fmax(XYXY[4 * sortedBody + 1], XYXY[4 * sortedBody + 3]) };
                        }

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
__launch_bounds__(THREADS5W2W, FACTOR5)
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


    __shared__ volatile int pos[MAXDEPTH * THREADS5W2W / WARPSIZE], node[MAXDEPTH * THREADS5W2W / WARPSIZE];



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
					if (depth < MAXDEPTH * THREADS5W2W / WARPSIZE)
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
					}//if depth <
                }

            }
            depth--;  // done with this level
        } while (depth >= j);

        
        // 14-05-2024
        /*
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
        */

        // update velocity

        real2 result = real2{ -idpi * v.y, idpi * v.x };
        veld[indexOfPoint] = result;

        if (calcEpsAst)
            epsast[indexOfPoint] = sqrt((d_1 + d_2 + d_3) / 3);
    }
}


__global__
__launch_bounds__(THREADS5S2W, FACTOR5)
void ForceCalculationKernelFromPanels2points(
    const int nnodesd, const int nbodiesd,
    const real itolsqd, const real epssqd,
    const int2* __restrict Mchildd,
    const int order, const real2* __restrict momsd,
    const double* __restrict dev_ptr_r, //начала и концы панелей
    const void* __restrict ptrd,
    int sizeOfElement,
    int offsetOfPointInElement,
    const int* __restrict MmortonCodesIdxd,
    const real2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
    const int npointsd, const real3* __restrict pointsd,
    const int* __restrict MmortonCodesIdxPointsd,
    real2* __restrict veld,
    const real2* __restrict Mlowerd,
    const real2* __restrict Mupperd, 
    int schemeType
)
{
    register int j, k, n, depth, base, sbase, pd, nd;
    register real2 p, v, dr, ps;
    register real r2;
    register const real2* mom;

    register real2 th;


    __shared__ volatile int pos[MAXDEPTH * THREADS5S2W / WARPSIZE], node[MAXDEPTH * THREADS5S2W / WARPSIZE];

   

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
        const int indexOfPoint = MmortonCodesIdxPointsd[k]; //k;//MmortonCodesIdxd[k];
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

            register real gm, gmFree, gmAttVortex, gmAttSource;
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

                real2 panBegin, panEnd;

                if (isVortex)
                {
                    n = chd - nbodiesd;
                    real* ptrValX = (real*)((char*)ptrd + MmortonCodesIdxd[n] * sizeOfElement + offsetOfPointInElement);
                    ps = real2{ *(ptrValX + 0), *(ptrValX + 1) };
                    panBegin = real2{ dev_ptr_r[MmortonCodesIdxd[n] * 4 + 0], dev_ptr_r[MmortonCodesIdxd[n] * 4 + 1] };
                    panEnd = real2{ dev_ptr_r[MmortonCodesIdxd[n] * 4 + 2], dev_ptr_r[MmortonCodesIdxd[n] * 4 + 3] };

                    gmFree = *(ptrValX + 6);
                    gmAttVortex = *(ptrValX + 7);
                    gmAttSource = *(ptrValX + 8);

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
                    
                    if (isVortex)
                    {   
                        //TODELETE//
                        //atomicAdd(&directCalc, 1);

                        real panLen = sqrt((panEnd.x - panBegin.x) * (panEnd.x - panBegin.x) + (panEnd.y - panBegin.y) * (panEnd.y - panBegin.y));
                        real2 vecS = p - panBegin; 
                        real2 vecP = p - panEnd;
                        real alpha = atan2(-vecP.y * vecS.x + vecP.x * vecS.y, vecP.x * vecS.x + vecP.y * vecS.y);
                        real2 u0 = real2{ (panEnd.x - panBegin.x) / panLen, (panEnd.y - panBegin.y) / panLen };
                        real lambda = 0.5 * log((vecS.x * vecS.x + vecS.y * vecS.y) / (vecP.x * vecP.x + vecP.y * vecP.y));
                        real2 skos0 = real2{ alpha * u0.y + lambda * u0.x, -alpha * u0.x + lambda * u0.y };

                        switch (schemeType)
                        {
                        case 1:
                        case 2:
                            v += ((gmFree + gmAttVortex) / panLen) * skos0;
                            break;

                        case -1:
                        case -2:
                            v += (gmAttSource / panLen) * skos0;
                            break;

                        };
                        

                        if ((schemeType == 2) || (schemeType == -2))
                        {
                            n = chd - nbodiesd;
                            real* ptrValX = (real*)((char*)ptrd + MmortonCodesIdxd[n] * sizeOfElement + offsetOfPointInElement);
                            real gmLinFree = *(ptrValX + 9);
                            real gmLinAttVortex = *(ptrValX + 10);
                            real gmLinAttSource = *(ptrValX + 11);
                            real2 vecA = vecP + vecS;
                            real2 vecB = u0;
                            real2 omega = (vecA.x * vecB.x + vecA.y * vecB.y) * vecB + real2{ vecA.y * vecB.x * vecB.y - vecA.x * vecB.y * vecB.y, -vecA.y * vecB.x * vecB.x + vecA.x * vecB.x * vecB.y }; // (vecA^ vecB)* (vecB.kcross());
                            real2 u1 = (0.5 / panLen) * omega;
                            real2 skos1 = real2{ alpha * u1.y + lambda * u1.x - u0.x, -alpha * u1.x + lambda * u1.y - u0.y };


                            switch (schemeType)
                            {             
                            case 2:
                                v += ((gmLinFree + gmLinAttVortex) / panLen) * skos1;
                                break;

                            case -2:
                                v += (gmLinAttSource / panLen) * skos1;
                                break;
                            };                            
                            
                        }
                    }
                    else                    
                    {
#ifdef CALCinDOUBLE
                        real f = gm / realmax(r2, epssqd);
#else
                        real f = fdividef(gm, realmax(r2, epssqd));
#endif
                        v += f * dr;
                    }

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
                    if (depth < MAXDEPTH * THREADS5S2W / WARPSIZE)
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
                    }//if depth <
                }

            }
            depth--;  // done with this level
        } while (depth >= j);


        /*if (calcEpsAst)
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
        }*/


        // update velocity
        real2 result = (schemeType < 0) ? real2{ idpi * v.x, idpi * v.y } : real2{ -idpi * v.y, idpi * v.x };
                
        //int newcntr = atomicAdd(&counterGPU, 1) + 1;
        
        //printf("k = %d, v = (%f, %f)\n", k, result.x, result.y);

        veld[indexOfPoint] = result;
        
        
        __syncthreads();

        //TODELETE//
        //if (newcntr == npointsd)
        //    printf("directCalc = %d\n", directCalc);

       /* if (calcEpsAst)
            epsast[indexOfPoint] = sqrt((d_1 + d_2 + d_3) / 3);*/
    }
}


__global__
__launch_bounds__(THREADS5I1I2, FACTOR5I1I2)
void I1I2CalculationKernel2(
    const int nnodesd, const int nbodiesd,
    const real minRd,
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



    __shared__ volatile int pos[MAXDEPTH * THREADS5I1I2 / WARPSIZE], node[MAXDEPTH * THREADS5I1I2 / WARPSIZE];



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
                
                rdi = max(epsast[indexInParticles], minRd);
                diffRadius = 8.0 * rdi;
                diffRadiusMonopole = 4.0 * rdi;

                

                // check if all threads agree that cell is far enough away (or is a body)
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
                    if (depth < MAXDEPTH * THREADS5I1I2 / WARPSIZE)
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
					}//if depth<
                }

            }
            depth--;  // done with this level
        } while (depth >= j);
               

        // update velocity
        I1d[indexInParticles] = i1;
        I2d[indexInParticles] = i2; 
    }
}//DiffusiveVeloI1I2

__global__
__launch_bounds__(THREADS5, FACTOR5)
void ClearI0I3(int nbodiesd, float3* __restrict i0i3Analyt, float* __restrict I0d, float2* __restrict I3d)
{
    for (int k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd; k += blockDim.x * gridDim.x)
    {
        i0i3Analyt[k].x = 0;
        i0i3Analyt[k].y = 0; 
        i0i3Analyt[k].z = 0;
        I0d[k] = 0;
        I3d[k].x = 0;
        I3d[k].y = 0;
    }
}

__global__
__launch_bounds__(THREADS5, FACTOR5)
void I0I3Recalc(int nbodiesd, const float3* __restrict i0i3Analyt, float* __restrict I0d, float2* __restrict I3d)
{
    for (int k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd; k += blockDim.x * gridDim.x)
    {
        if (i0i3Analyt[k].x != 0.0)
        {
            I0d[k] = i0i3Analyt[k].x;
            I3d[k].x = i0i3Analyt[k].y;
            I3d[k].y = i0i3Analyt[k].z;
        }
    }
}

__global__
__launch_bounds__(THREADS5I0I3, FACTOR5I0I3)
void I0I3CalculationKernel2(
    const int nnodesd, const int nbodiesd,
    const real minRd,
    const int2* __restrict Mchildd,
    const int order, const real2* __restrict momsd,
    const real3* __restrict vtxd,
    const int* __restrict MmortonCodesIdxd,
    const real2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
    float* __restrict I0d,
    float2* __restrict I3d,
    const real2* __restrict Mlowerd,
    const real2* __restrict Mupperd,
    const real* __restrict epsast, const real* __restrict meanEpsd, const int npan, const real* __restrict pansd, real* __restrict visstrd, float3* __restrict i0i3Analyt, int* range)

{
    register int j, k, n, depth, base, sbase, pd, nd;
    register float2 p, q, ps;
    register float i0;
    register float2 i3;
    register const double2* mom;

    __shared__ volatile int pos[MAXDEPTH * THREADS5I0I3 / WARPSIZE], node[MAXDEPTH * THREADS5I0I3 / WARPSIZE];

    // figure out first thread in each warp (lane 0)
    base = threadIdx.x / WARPSIZE;
    sbase = base * WARPSIZE;
    j = base * MAXDEPTH;

    float rdi; 
 
    float iDDomRad;

    float d;
    float2 beg, end;

    float lenj, lenj_m;
    float2 tau;
    float s;
    float2 norm;
    float vs, vsm;
    float meanepsj2;

    float v0x, v0y;
    float hx, hy;
    float xix, xiy, lxi;
    float expon;
    float mnx, mny;
    int new_n;
    float den;
    float xi_mx, xi_my, lxi_m;
    float mnog1;

    float d1, d2, d3, d4;
    float s1, s2, s3, s4;
 
    __syncthreads();
    __threadfence_block();

    
    // iterate over all bodies assigned to thread
    for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npan; k += blockDim.x * gridDim.x)
    {
        beg.x = (float)pansd[k * 4 + 0];
        beg.y = (float)pansd[k * 4 + 1];
        end.x = (float)pansd[k * 4 + 2];
        end.y = (float)pansd[k * 4 + 3];

        vs = 0.0f;
        //meanepsj2 = sqr(meanEpsd[k]); // commented 30-05
        meanepsj2 = sqrf(myMax((float)meanEpsd[k], (float)minRd));


        lenj = length(end - beg);

        tau = normalize(end - beg);        

        norm.x = tau.y;
        norm.y = -tau.x;

        p = 0.5f * (beg + end);


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

            register float gm;
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
                    i0 = 0.0f;
                    i3.x = 0.0f;
                    i3.y = 0.0f;

                    n = chd - nbodiesd;

                    ps = float2{ (float)vtxd[MmortonCodesIdxd[n]].x, (float)vtxd[MmortonCodesIdxd[n]].y };
                    gm = (float)(vtxd[MmortonCodesIdxd[n]].z);

                    rdi = myMax((float)epsast[MmortonCodesIdxd[n]], (float)minRd);
                    //rdi = epsast[MmortonCodesIdxd[n]]; // commented 30-05
                    iDDomRad = 1.0f / rdi;
                }
                else
                {
                    register const int srtT = MindexSortTd[chd];
                    n = (nnodesd - 1) - srtT;

                    ps.x = (float)Mposd[chd].x;
                    ps.y = (float)Mposd[chd].y;
                    mom = momsd + (srtT * order);
                    gm = (float)(mom[0].x);
                }

                q = ps - p;


                s = (q & tau);
                d = fabsf(q & norm);

                const float farDist = 50.0f;
                const float closeDist = 5.0f;
                float smax; 
                bool cond = false;
                if (!(isVortex))
                {
                    float2 low{ (float)Mlowerd[chd].x, (float)Mlowerd[chd].y };
                    float2 up{ (float)Mupperd[chd].x, (float)Mupperd[chd].y };

                    d1 = fabsf((low - p) & norm);
                    s1 = (low - p) & tau;
                    float2 crd;
                    crd.x = low.x;
                    crd.y = up.y;
                    d2 = fabsf((crd - p) & norm);
                    s2 = (crd - p) & tau;
                    crd.x = up.x;
                    crd.y = low.y;
                    d3 = fabsf((crd - p) & norm);
                    s3 = (crd - p) & tau;
                    d4 = fabsf((up - p) & norm);
                    s4 = (up - p) & tau;

                    d = myMin(myMin(myMin(d1, d2), d3), d4);
                    s = myMin(myMin(myMin(s1, s2), s3), s4);
                    smax  = myMax(myMax(myMax(s1, s2), s3), s4);

                    cond = (low.x < p.x) && (p.x < up.x) && (low.y < p.y) && (p.y < up.y);                  
                }

                                                                
                if (__all_sync(0xffffffff, ((d > farDist * lenj) || ( (fabsf(s) > farDist * lenj) && (s*smax > 0))) && !(cond)))
                {

                }
                else if (isVortex)
                {   

                    v0x = tau.x * lenj;
                    v0y = tau.y * lenj;

                    if ((d > closeDist * lenj) || (fabsf(s) > closeDist * lenj))
                    {
                        xix = q.x * iDDomRad;
                        xiy = q.y * iDDomRad;
                        lxi = sqrtf(xix * xix + xiy * xiy);

                        expon = expf(-lxi) * lenj;
                        mnx = norm.x * expon;
                        mny = norm.y * expon;

                        i0 = (xix * mnx + xiy * mny) * (lxi + 1.0) / (lxi * lxi);
                        i3.x = mnx;
                        i3.y = mny;                        
                        
                        vs += gm * expon / (pif * meanepsj2);

                    }
                    else if ((d > 0.01f * lenj) || (fabsf(s) > 0.45f * lenj))
                    {
                        den = (fabsf(s) < 0.5f * lenj) ? d : (fabsf(s) + d - 0.5f * lenj);
                        new_n = (int)min((int)max(ceilf(10.0f * lenj / den), 1.0f), 20);

                        hx = v0x / new_n;
                        hy = v0y / new_n;
                        vsm = 0;
                        for (int m = 0; m < new_n; ++m)
                        {
                            xi_mx = (ps.x - (beg.x + hx * (m + 0.5f))) * iDDomRad;
                            xi_my = (ps.y - (beg.y + hy * (m + 0.5f))) * iDDomRad;

                            lxi_m = sqrtf(xi_mx * xi_mx + xi_my * xi_my);

                            lenj_m = lenj / new_n;
                            expon = expf(-lxi_m) * lenj_m;

                            mnx = norm.x * expon;
                            mny = norm.y * expon;

                            
                            i0 += (xi_mx * mnx + xi_my * mny) * (lxi_m + 1.0) / (lxi_m * lxi_m);
                            i3.x += mnx;
                            i3.y += mny;                             

                            vsm += expon;
                        }//for m
                        
                        vsm *= gm / (pif * meanepsj2);
                        vs += vsm;
                    }
                    else
                    {
                        i0 = -pif * rdi;

                        mnog1 = 2.0f * rdi * (1.0f - expf(-lenj * 0.5f * iDDomRad) * coshf(fabsf(s) * iDDomRad));
                        i3.x = mnog1 * norm.x;
                        i3.y = mnog1 * norm.y;

                        vs += mnog1 * gm / (pif * meanepsj2);

                        i0i3Analyt[MmortonCodesIdxd[n]].x = i0;
                        i0i3Analyt[MmortonCodesIdxd[n]].y = i3.x;
                        i0i3Analyt[MmortonCodesIdxd[n]].z = i3.y;
                    }


//#ifdef  __CUDA_ARCH__
//#if __CUDA_ARCH__ < 600
                   //myAtomicAdd(I0d + MmortonCodesIdxd[n], i0);
                   //myAtomicAdd((float*)I3d + 2 * MmortonCodesIdxd[n], i3.x);
                   //myAtomicAdd((float*)I3d + 2 * MmortonCodesIdxd[n] + 1, i3.y);                    
//#else                   
                    //if ((ps.y > 0) && (ps.x > 0.18) && (ps.x < 0.19))
                    //if (MmortonCodesIdxd[n] == 140)
                    //    printf("vortex: %d, from_pan = %d, i0 = %f, i3 = {%f, %f}\n", MmortonCodesIdxd[n], k, i0, i3.x, i3.y);

                   atomicAdd(I0d + MmortonCodesIdxd[n], i0);
                   atomicAdd((float*)I3d + 2 * MmortonCodesIdxd[n], i3.x);
                   atomicAdd((float*)I3d + 2 * MmortonCodesIdxd[n] + 1, i3.y);
//#endif
//#endif //  __CUDA_ARCH__	
                }// if (isVortex)
                else //идем глубже по дереву
                {
                    
                    
                    
                    if (depth < MAXDEPTH * THREADS5I0I3 / WARPSIZE)
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
                    }//if depth<
                } 
            }
            depth--;  // done with this level
        } while (depth >= j);

        // update viscous stress 
        if (k < npan)
            visstrd[k] = (double)vs;
    }//for k
}//DiffusiveVeloI0I3





__global__
void VerticesToControlPointsKernel(int nTotPan, const real* dev_ptr_pt, real* pointsl)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nTotPan)
    {
        pointsl[2 * i + 0] = 0.5 * (dev_ptr_pt[4 * i + 0] + dev_ptr_pt[4 * i + 2]);
        pointsl[2 * i + 1] = 0.5 * (dev_ptr_pt[4 * i + 1] + dev_ptr_pt[4 * i + 3]);
    }
}


__global__
void VerticesToControlPointsVortexesKernel(int nTotPan, const real* dev_ptr_pt, real* pointsl, const real* gaml)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < nTotPan)
    {
        pointsl[3 * i + 0] = 0.5 * (dev_ptr_pt[4 * i + 0] + dev_ptr_pt[4 * i + 2]);
        pointsl[3 * i + 1] = 0.5 * (dev_ptr_pt[4 * i + 1] + dev_ptr_pt[4 * i + 3]);
        pointsl[3 * i + 2] = gaml[i];
    }
}


__global__
void VerticesAndSheetsToPanelPointsKernel(int npnli, const double* dev_ptr_pt, 
    const double* dev_ptr_freeVortexSheet, const double* dev_ptr_freeVortexSheetLin,
    const double* dev_ptr_attachedVortexSheet, const double* dev_ptr_attachedVortexSheetLin,
    const double* dev_ptr_attachedSourceSheet, const double* dev_ptr_attachedSourceSheetLin,
    double* panelPoints, int schemeType)
{
    size_t i = blockIdx.x * blockDim.x + threadIdx.x;
    double dx, dy, dl;
    if (i < npnli)
    {
        panelPoints[12 * i + 0] = 0.5 * (dev_ptr_pt[4 * i + 0] + dev_ptr_pt[4 * i + 2]);
        panelPoints[12 * i + 1] = 0.5 * (dev_ptr_pt[4 * i + 1] + dev_ptr_pt[4 * i + 3]);

        panelPoints[12 * i + 2] = dev_ptr_pt[4 * i + 0];
        panelPoints[12 * i + 3] = dev_ptr_pt[4 * i + 1];
        panelPoints[12 * i + 4] = dev_ptr_pt[4 * i + 2];
        panelPoints[12 * i + 5] = dev_ptr_pt[4 * i + 3];

        dx = dev_ptr_pt[4 * i + 2] - dev_ptr_pt[4 * i + 0];
        dy = dev_ptr_pt[4 * i + 3] - dev_ptr_pt[4 * i + 1];
        dl = sqrt(dx * dx + dy * dy);
           

        panelPoints[12 * i + 6] = dl * dev_ptr_freeVortexSheet[i];

        if (dev_ptr_attachedVortexSheet != nullptr)
            panelPoints[12 * i + 7] = dl * dev_ptr_attachedVortexSheet[i];
        else
            panelPoints[12 * i + 7] = 0.0;

        if (dev_ptr_attachedSourceSheet != nullptr)
            panelPoints[12 * i + 8] = dl * dev_ptr_attachedSourceSheet[i];
        else
            panelPoints[12 * i + 8] = 0.0;

        panelPoints[12 * i + 9] = 0.0;
        panelPoints[12 * i + 10] = 0.0;
        panelPoints[12 * i + 11] = 0.0;

        if (abs(schemeType) == 2)
        {
            panelPoints[12 * i + 9] = dl * dev_ptr_freeVortexSheetLin[i];
            
            if (dev_ptr_attachedVortexSheetLin != nullptr)
                panelPoints[12 * i + 10] = dl * dev_ptr_attachedVortexSheetLin[i];
            else
                panelPoints[12 * i + 10] = 0.0;
            
            if (dev_ptr_attachedSourceSheetLin != nullptr)
                panelPoints[12 * i + 11] = dl * dev_ptr_attachedSourceSheetLin[i];
            else
                panelPoints[12 * i + 11] = 0.0;
        }

        //if (i == 71 || i == 72)
        //    printf("%d: pos:(%f, %f) const:(%f, %f, %f), lin:(%f, %f, %f)\n", (int)i, panelPoints[12 * i + 0], panelPoints[12 * i + 1], panelPoints[12 * i + 6], panelPoints[12 * i + 7], panelPoints[12 * i + 8], panelPoints[12 * i + 9], panelPoints[12 * i + 10], panelPoints[12 * i + 11]);
    }
}



__global__
__launch_bounds__(THREADS5rhs, FACTOR5rhs)
void RhsCalculationKernel
(
    const int nnodesd,  // Число узлов дерева (по вихрям)
    const int nbodiesd, // Число вихрей в пелене

    const real itolsqd, // 1/theta^2
    const int2* __restrict Mchildd, //Массив потомков узлов дерева

    const int order, //порядок мультипольного разложения
    const real2* __restrict momsd, //моменты узлов дерева
    const real3* __restrict vtxd,  //вихри в пелене
    const int* __restrict MmortonCodesIdxd, //порядок сортировки вихрей в пелене

    const real2* __restrict Mposd, //массив координат центров ячеек дерева
    const int* __restrict MindexSortd, //порядок сортировки внутренних узлов дерева
    const int* __restrict MindexSortTd,//обратный порядок сортировки внутренних узлов дерева
     
    const int npointsd,                  //количество панелей на теле
    const real* __restrict dev_ptr_pt,   //массив четверок (beg.x, beg.y, end.x, end.y)
    const real2* __restrict pointsd,     //массив центров панелей
    
    real2* __restrict Ed,                //коэффициенты локального разложения
    const int* __restrict MmortonCodesIdxPointsd,//порядок сортировки центров панелей
    
    real* __restrict veld,               //куда сохранять ответ
    real* __restrict vellind,            //куда сохранять ответ

    const real2* __restrict Mlowerd,     //координаты нижнего левого угла ячейки дерева
    const real2* __restrict Mupperd      //координаты верхнего правого угла ячейки дерева
)
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
        p = pointsd[indexOfPoint];

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
                        real2 ss = ps - beg;
                        real2 pp = ps - end;

                        real alpha = atan2(pp.x * ss.y - pp.y * ss.x, pp.x * ss.x + pp.y * ss.y);

                        real tempVel = gm * alpha;
                        
                        if (r2 > 1e-20)
                            val -= tempVel;
                        //else
                        //    val -= 0.5 /idpi * gm;                       

                        if (vellind != nullptr)
                        {

                            real2 u1;
                            u1.x = 0.5 * idlen * ((pp.x + ss.x) * tau.x * tau.x \
                                + 2.0 * (pp.y + ss.y) * tau.x * tau.y - (pp.x + ss.x) * tau.y * tau.y);
                            u1.y = 0.5 * idlen * (-(pp.y + ss.y) * tau.x * tau.x \
                                + 2.0 * (pp.x + ss.x) * tau.x * tau.y + (pp.y + ss.y) * tau.y * tau.y);

                            real lambda = 0.5 * log((ss.x * ss.x + ss.y * ss.y) / (pp.x * pp.x + pp.y * pp.y));

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
                        
                        //printf("far!\n");
                        val += (-v.y * rPan.x + v.x * rPan.y);
                        
                        if (vellind != nullptr)
                            vallin += (-vL.y * rPan.x + vL.x * rPan.y);    
                        
                    }                    
                }
                else //идем глубже по дереву
                {
                    if (depth < MAXDEPTH * THREADS5rhs / WARPSIZE)
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
					}//if depth<
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



/// \brief Вспомогательная функция, которая определяет, находится ли панель itI на контуре после панели itJ
    ///
    /// \param[in] itI константная ссылка на обертку для второй ("правой") панели
    /// \param[in] itJ константная ссылка на обертку для первой ("левой") панели
__device__ bool isEqual(double2 paniBeg, const double2 panjEnd)
{
    return ((paniBeg.x - panjEnd.x) * (paniBeg.x - panjEnd.x) + (paniBeg.y - panjEnd.y) * (paniBeg.y - panjEnd.y) < 1e-20);
}

/// \brief Вспомогательная функция вычисления угла между векторами (в диапазоне (-pi...pi]) (со знаком, поворот от первого ко второму)
///
/// \param[in] p константная ссылка на первый вектор
/// \param[in] s константная ссылка на второй вектор 
__device__ double Alpha(double2 p, double2 s)
{
    return atan2(p.x * s.y - p.y * s.x, p.x * s.x + p.y * s.y);        
}

/// \brief  Вспомогательная функция вычисления логарифма отношения норм векторов
/// 
/// \param[in] p константная ссылка на первый вектор
/// \param[in] s константная ссылка на второй вектор 
__device__ double Lambda(double2 p, double2 s)
{
    return 0.5 * log((s.x * s.x + s.y * s.y) / (p.x * p.x + p.y * p.y));    
}


/// Вспомогательная функция вычисления величины \f$ (\vec a \cdot \vec b) \cdot \vec c + (\vec a \times \vec b) \times \vec c \f$
__device__ double2 Omega(double2 a, double2 b, double2 c)
{
    double adotb = a.x * b.x + a.y * b.y;
    double acrossb = a.x * b.y - a.y * b.x;
    return real2{ adotb * c.x - acrossb * c.y, adotb * c.y + acrossb * c.x};
};
    
__global__
__launch_bounds__(THREADS5slae, FACTOR5slae)
void MatrToVecCalculationKernel
(
    const int nnodesd,  // Число узлов дерева (по вихрям)
    const int nbodiesd, // Число вихрей в пелене

    const real itolsqd, // 1/theta^2
    const int2* __restrict Mchildd, //Массив потомков узлов дерева

    const int order, //порядок мультипольного разложения
    const real2* __restrict momsd, //моменты узлов дерева
    const real* __restrict rpnl,  //вихри в пелене
    const real* dev_ptr_freeVortexSheet,
    const real* dev_ptr_freeVortexSheetLin,

    const int* __restrict MmortonCodesIdxd, //порядок сортировки вихрей в пелене

    const real2* __restrict Mposd, //массив координат центров ячеек дерева
    const int* __restrict MindexSortd, //порядок сортировки внутренних узлов дерева
    const int* __restrict MindexSortTd,//обратный порядок сортировки внутренних узлов дерева

    const int npointsd,                  //количество панелей на теле
    const real* __restrict dev_ptr_pt,   //массив четверок (beg.x, beg.y, end.x, end.y)
    const real2* __restrict pointsd,     //массив центров панелей

    real2* __restrict Ed,                //коэффициенты локального разложения
    const int* __restrict MmortonCodesIdxPointsd,//порядок сортировки центров панелей

    real* __restrict veld,               //куда сохранять ответ
    real* __restrict vellind,            //куда сохранять ответ

    const real2* __restrict Mlowerd,     //координаты нижнего левого угла ячейки дерева
    const real2* __restrict Mupperd      //координаты верхнего правого угла ячейки дерева
)
{
    register int j, k, n, depth, base, sbase, pd, nd;
    register real2 p, dr, ps;
    register real r2;
    register const real2* mom;
    register real val, vallin;
    
    __shared__ volatile int pos[MAXDEPTH * THREADS5slae / WARPSIZE], node[MAXDEPTH * THREADS5slae / WARPSIZE];

    // figure out first thread in each warp (lane 0)
    base = threadIdx.x / WARPSIZE;
    sbase = base * WARPSIZE;
    j = base * MAXDEPTH;

    __syncthreads();
    __threadfence_block();

    // iterate over all bodies assigned to thread
    for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npointsd; k += blockDim.x * gridDim.x)
    {
        //if (k == 0)
        //    printf("mom0[187]: %f\n", momsd[187]);

        const int indexOfPoint = MmortonCodesIdxPointsd[k];
        p = pointsd[indexOfPoint];
        //if (indexOfPoint == 0)
        //    printf("x,y = %f, %f\n", p.x, p.y);

        real2 beg, end;
        beg = real2{ dev_ptr_pt[4 * indexOfPoint + 0], dev_ptr_pt[4 * indexOfPoint + 1] };
        end = real2{ dev_ptr_pt[4 * indexOfPoint + 2], dev_ptr_pt[4 * indexOfPoint + 3] };

        real2 di = real2{ end.x - beg.x, end.y - beg.y };

        real dlen2 = (end.x - beg.x) * (end.x - beg.x) + (end.y - beg.y) * (end.y - beg.y);
        real dilen = sqrt(dlen2);
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

                    ps = real2{ 0.5*(rpnl[MmortonCodesIdxd[n]*4 + 0] + rpnl[MmortonCodesIdxd[n] * 4 + 2]), 0.5 * (rpnl[MmortonCodesIdxd[n] * 4 + 1] + rpnl[MmortonCodesIdxd[n] * 4 + 3]) };
                    
                    double lj = sqrt((rpnl[MmortonCodesIdxd[n] * 4 + 2]- rpnl[MmortonCodesIdxd[n] * 4 + 0]) * (rpnl[MmortonCodesIdxd[n] * 4 + 2] - rpnl[MmortonCodesIdxd[n] * 4 + 0]) +
                        (rpnl[MmortonCodesIdxd[n] * 4 + 3] - rpnl[MmortonCodesIdxd[n] * 4 + 1]) * (rpnl[MmortonCodesIdxd[n] * 4 + 3] - rpnl[MmortonCodesIdxd[n] * 4 + 1]) );

                    gm = lj * dev_ptr_freeVortexSheet[MmortonCodesIdxd[n]];// vtxd[MmortonCodesIdxd[n]].z;

                    //if (indexOfPoint == 0)
                    //{
                    //    if (MmortonCodesIdxd[n] == 71 || MmortonCodesIdxd[n] == 72)
                    //    {
                    //        printf("gm[%d] = %f\n", (int)MmortonCodesIdxd[n], gm);
                    //        printf("freesh[%d] = %f\n", (int)MmortonCodesIdxd[n], dev_ptr_freeVortexSheet[MmortonCodesIdxd[n]]);
                    //    }
                    //}


                    sumSide2 = real2{ fabs(rpnl[MmortonCodesIdxd[n] * 4 + 2] - rpnl[MmortonCodesIdxd[n] * 4 + 0]),
                        fabs(rpnl[MmortonCodesIdxd[n] * 4 + 3] - rpnl[MmortonCodesIdxd[n] * 4 + 1]) };
                }
                else
                {
                    register const int srtT = MindexSortTd[chd];
                    //printf("srtT = %d\n", (int)srtT);
                    n = (nnodesd - 1) - srtT;
                    ps = Mposd[chd];
                    mom = momsd + (srtT * order);
                    gm = mom[0].x;
                    sumSide2 = Mupperd[chd] - Mlowerd[chd]; /*Msized[chd]*/;

                    //if (indexOfPoint == 0 && chd == 78)
                    //{
                    //    printf("srtT = %d, n = %d, ps = (%f, %f), mom0.x = %f, lower=(%f, %f), upper=(%f, %f)\n", srtT, (int)n, ps.x, ps.y, mom[0].x, Mlowerd[chd].x, Mlowerd[chd].y, Mupperd[chd].x, Mupperd[chd].y);
                    //}


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
                        real tempVelNew;

                        double2 ptpanBeg = real2{ rpnl[MmortonCodesIdxd[n] * 4 + 0], rpnl[MmortonCodesIdxd[n] * 4 + 1] };
                        double2 ptpanEnd = real2{ rpnl[MmortonCodesIdxd[n] * 4 + 2], rpnl[MmortonCodesIdxd[n] * 4 + 3] };

                        double lj = sqrt((rpnl[MmortonCodesIdxd[n] * 4 + 2] - rpnl[MmortonCodesIdxd[n] * 4 + 0]) * (rpnl[MmortonCodesIdxd[n] * 4 + 2] - rpnl[MmortonCodesIdxd[n] * 4 + 0]) +
                            (rpnl[MmortonCodesIdxd[n] * 4 + 3] - rpnl[MmortonCodesIdxd[n] * 4 + 1]) * (rpnl[MmortonCodesIdxd[n] * 4 + 3] - rpnl[MmortonCodesIdxd[n] * 4 + 1]));

                        double gmlin = 0.0;
                        if (dev_ptr_freeVortexSheetLin != nullptr)                      
                            gmlin = lj * dev_ptr_freeVortexSheetLin[MmortonCodesIdxd[n]];                        

                        if (r2 > 1e-20)
                        {

                            //bool cond = ((sumSide2.x + sumSide2.y) * (sumSide2.x + sumSide2.y) + dlen2)* itolsqd < r2;
                            //if ((int)indexOfPoint == 0)
                            //    printf("i = %d, j = %d, isVortex = %d, farZone = %d, S^2 = %f, dlen2 = %f, r2 = %f\n",
                            //        (int)indexOfPoint,
                            //        (int)MmortonCodesIdxd[n],
                            //        (int)isVortex,
                            //        (int)cond,
                            //        (sumSide2.x + sumSide2.y) * (sumSide2.x + sumSide2.y),
                            //        dlen2,
                            //        r2
                            //    );
                            double2 dj = real2{ ptpanEnd.x - ptpanBeg.x, ptpanEnd.y - ptpanBeg.y };
                            double djlen = sqrt(dj.x * dj.x + dj.y * dj.y);
                            double ilenj = 1.0 / djlen;
                            const double2 tauj = real2{ ilenj * dj.x, ilenj * dj.y };

                            double2 p1 = real2{ end.x - ptpanEnd.x, end.y - ptpanEnd.y };
                            double2 s1 = real2{ end.x - ptpanBeg.x, end.y - ptpanBeg.y };
                            double2 p2 = real2{ beg.x - ptpanEnd.x, beg.y - ptpanEnd.y };
                            double2 s2 = real2{ beg.x - ptpanBeg.x, beg.y - ptpanBeg.y };

                            double3 alpha = double3{ \
                                isEqual(ptpanBeg, end) ? 0.0 : Alpha(s2, s1), \
                                Alpha(s2, p1), \
                                isEqual(beg, ptpanEnd) ? 0.0 : Alpha(p1, p2) \
                            };

                            double3 lambda = double3{ \
                                isEqual(ptpanBeg, end) ? 0.0 : Lambda(s2, s1), \
                                Lambda(s2, p1), \
                                isEqual(beg, ptpanEnd) ? 0.0 : Lambda(p1, p2) \
                            };

                            double2 v00_0 = Omega(s1, tau, tauj);
                            double2 v00_1 = Omega(di, tau, tauj);
                            v00_1.x *= -1.0;
                            v00_1.y *= -1.0;

                            double2 v00_2 = Omega(p2, tau, tauj);
                            

                            double2 i00 = real2{ 
                                ilenj * (alpha.x * v00_0.x + alpha.y * v00_1.x + alpha.z * v00_2.x \
                                    - (lambda.x * v00_0.y + lambda.y * v00_1.y + lambda.z * v00_2.y)),
                                ilenj * (alpha.x * v00_0.y + alpha.y * v00_1.y + alpha.z * v00_2.y \
                                    + (lambda.x * v00_0.x + lambda.y * v00_1.x + lambda.z * v00_2.x))
                                };

                            tempVelNew = -gm * (i00.x * tau.x + i00.y * tau.y);
                            
                            //printf("%f, %f\n", tempVel, tempVelNew);
                            
                            val -= tempVelNew;

                            //if (indexOfPoint == 0)
                            //{
                            //    if (MmortonCodesIdxd[n] == 71 || MmortonCodesIdxd[n] == 72)
                            //        printf("contrib %d: %f\n", (int)MmortonCodesIdxd[n], tempVelNew);
                            //}


                            
                            if (vellind != nullptr)
                            {
                            double s1len2 = s1.x * s1.x + s1.y * s1.y;                            

                            double2 om1 = Omega(s1, tau, tauj);
                            double2 om2 = Omega(real2{ s1.x + p2.x, s1.y + p2.y }, tauj, tauj);
                            double sc = (p1.x + s1.x) * tauj.x + (p1.y + s1.y) * tauj.y;

                            double2 v01_0 = real2{ 0.5 / (djlen) * (sc * om1.x - s1len2 * tau.x), 0.5 / (djlen) * (sc * om1.y - s1len2 * tau.y) };
                            double2 v01_1 = real2{ -0.5 * dilen / djlen * om2.x, -0.5 * dilen / djlen * om2.y };
                            

                            double2 i01 = real2{
                                ilenj * ((alpha.x + alpha.z) * v01_0.x + (alpha.y + alpha.z) * v01_1.x\
                                    - (((lambda.x + lambda.z) * v01_0.y + (lambda.y + lambda.z) * v01_1.y) - 0.5 * dilen * tauj.y)),
                                ilenj * ((alpha.x + alpha.z) * v01_0.y + (alpha.y + alpha.z) * v01_1.y\
                                    + (((lambda.x + lambda.z) * v01_0.x + (lambda.y + lambda.z) * v01_1.x) - 0.5 * dilen * tauj.x))
                            };

                            double2 om3 = Omega(real2{ s1.x + p2.x, s1.y + p2.y }, tau, tau);

                            double2 v10_0 = real2{ -0.5 / dilen * (((s1 + s2) & tau) * om1.x - s1len2 * tauj.x), -0.5 / dilen * (((s1 + s2) & tau) * om1.y - s1len2 * tauj.y) };
                            double2 v10_1 = real2{ 0.5 * djlen / dilen * om3.x, 0.5 * djlen / dilen * om3.y };
                            

                            double2 i10 = real2{
                                ilenj * ((alpha.x + alpha.z) * v10_0.x + alpha.z * v10_1.x \
                                - (((lambda.x + lambda.z) * v10_0.y + lambda.z * v10_1.y) + 0.5 * djlen * tau.y)),
                                ilenj * ((alpha.x + alpha.z) * v10_0.y + alpha.z * v10_1.y \
                                + (((lambda.x + lambda.z) * v10_0.x + lambda.z * v10_1.x) + 0.5 * djlen * tau.x)),
                            };
       
                            double2 om4 = Omega(real2{ s1.x - 3.0 * p2.x, s1.y - 3.0 * p2.y }, tau, tauj);
                            double2 om5 = Omega(di, tauj, tauj);
                            double2 om6 = Omega(dj, tau, tau);
                            double sc4 = s1.x * om4.x + s1.y * om4.y;

                            

                            double2 v11_0 = real2{
                                1.0 / (12.0 * dilen * djlen) * (2.0 * sc4 * om1.x - s1len2 * (s1.x - 3.0 * p2.x)) - 0.25 * om1.x,
                                1.0 / (12.0 * dilen * djlen) * (2.0 * sc4 * om1.y - s1len2 * (s1.y - 3.0 * p2.y)) - 0.25 * om1.y
                            };
                            double2 v11_1 = real2{ -dilen / (12.0 * djlen) * om5.x, -dilen / (12.0 * djlen) * om5.y };
                            double2 v11_2 = real2{ -djlen / (12.0 * dilen) * om6.x, -djlen / (12.0 * dilen) * om6.y };
                                                        

                            double2 i11 = {
                                ilenj * ((alpha.x + alpha.z) * v11_0.x + (alpha.y + alpha.z) * v11_1.x + alpha.z * v11_2.x\
                                - ((lambda.x + lambda.z) * v11_0.y + (lambda.y + lambda.z) * v11_1.y + lambda.z * v11_2.y \
                                    + 1.0 / 12.0 * (djlen * tau.y + dilen * tauj.y - 2.0 * om1.y))),
                                ilenj * ((alpha.x + alpha.z) * v11_0.y + (alpha.y + alpha.z) * v11_1.y + alpha.z * v11_2.y\
                                + ((lambda.x + lambda.z) * v11_0.x + (lambda.y + lambda.z) * v11_1.x + lambda.z * v11_2.x \
                                    + 1.0 / 12.0 * (djlen * tau.x + dilen * tauj.x - 2.0 * om1.x)))
                            };
                            double tempVelNewB = -gmlin * (i01.x * tau.x + i01.y * tau.y);
                            //printf("gm = %f, gmlin = %f\n", gm, gmlin);

                            val -= tempVelNewB;


                            double tempVelNewC = -gm * (i10.x * tau.x + i10.y * tau.y);
                            double tempVelNewD = -gmlin * (i11.x * tau.x + i11.y * tau.y);
                            vallin -= (tempVelNewC + tempVelNewD);
                            }

                           
                        }//if (posI != posJ)



                    }
                    else
                    {
                        //if (indexOfPoint == 0)
                        //{
                        //    //printf("[%f, %f] x [%f, %f]\n", Mlowerd[chd].x, Mupperd[chd].x, Mlowerd[chd].y, Mupperd[chd].y);
                        //    printf("{{%f, %f}, {%f, %f}, mom.x=%f, gm=%f, SrtT=%d, chd=%d},\n", Mlowerd[chd].x, Mlowerd[chd].y, Mupperd[chd].x, Mupperd[chd].y, mom[0].x, gm, (int)(MindexSortTd[chd]), (int)chd);
                        //}
                        
                        
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
                                //if (sqrt(theta.x * theta.x + theta.y * theta.y) > 1e+16)
                                //{
                                //    printf("theta[%d] = (%f, %f)\n", q, theta.x, theta.y);
                                //    printf("mom[%d] = (%f, %f)\n", q - s, mom[q - s].x, mom[q - s].y);
                                //}
                                
                                Eloc[s] += (cftStart * cftDiag / factorial) * multzA(theta, mom[q - s]);

                                cftDiag = -cftDiag;
                                factorial *= (q - s + 1);
                            }
                            theta = ((q + 1) / dist2) * multz(theta, p - ps);
                            cftStart = -cftStart;
                        }



                        real2 v = Eloc[0];
                        real2 vL{ 0.0, 0.0 };


                        //if (indexOfPoint == 0 && chd == 78)
                        //{
                        //    //printf("order = %d\n", (int)order);
                        //    printf("val += {%f, %f}\n", v.x, v.y);
                        //}

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
                            iFact /= (k + 1);
                            mulP = multz(mulP, kp);
                            mulM = multz(mulM, km);
                            taudL = taudL * 0.5;
                            //if (fabs(Eloc[k].x) + fabs(Eloc[k].y) < 1e+10)
                            v += iFact * multz(Eloc[k], multzA(taudL, mulP - mulM));


                            if (vellind != nullptr)
                                vL += (iFact / (k + 2)) * multz(Eloc[k], multzA(multz(taudL, taudLc), multz(mulP, (k + 1) * rPan - dPos2) + multz(mulM, (k + 1) * rPan + dPos2)));
                        }

                        //printf("far!\n");



                        
                        //if (indexOfPoint == 0 && chd==78)
                        //    printf("val_old = %f\n", val);

                        //if (indexOfPoint == 0 && chd == 78)
                        //    printf("add, v={%f, %f}, rPan={%f, %f}\n", v.x, v.y, rPan.x, rPan.y);
                            
                        val += (-v.y * rPan.x + v.x * rPan.y);
                            
                        //if (indexOfPoint == 0 && chd == 78)
                        //    printf("val_added = %f\n", (-v.y * rPan.x + v.x * rPan.y));

                        if (vellind != nullptr)
                            vallin += (-vL.y * rPan.x + vL.x * rPan.y);
                    }
                }
                else //идем глубже по дереву
                {
                    if (depth < MAXDEPTH * THREADS5slae / WARPSIZE)
                    {

                        // push cell onto stack
                        if (sbase == threadIdx.x)
                        {  // maybe don't push and inc if last child
                            if (depth >= MAXDEPTH * THREADS5slae / WARPSIZE)
                                printf("?????????????????????????????????????????????????????????????????\n");

                            pos[depth] = pd;
                            node[depth] = nd;
                        }
                        depth++;
                        pd = 0;
                        nd = n;
                        chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];
                    }//if depth<
                }

            }
            depth--;  // done with this level
        } while (depth >= j);


        // update velocity
        veld[indexOfPoint] = val * idpi * idlen;

        if (vellind != nullptr)
            vellind[indexOfPoint] = vallin * idpi * idlen;

        //if (indexOfPoint == 0)
        //{
        //    //printf("order = %d\n", (int)order);
        //    printf("veld[0] = %f\n", val* idpi * idlen);
        //}

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
        
    cudaFuncSetCacheConfig(SummarizationKernel2_1, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_2, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_3, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_4, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_5, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_6, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_7, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_8, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_9, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_10, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_11, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_12, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_13, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_14, cudaFuncCachePreferL1);
    cudaFuncSetCacheConfig(SummarizationKernel2_15, cudaFuncCachePreferL1);
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
        int nbodies,
        const void* __restrict ptrl,
        int sizeOfElement,
        int offsetOfPointInElement)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        MBoundingBoxKernel << <blocks * FACTOR1, THREADS1 >> > (nbodies, ptrl, sizeOfElement, offsetOfPointInElement, (real2*)Mposl, (real2*)maxrl, (real2*)minrl);
        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("Mkernel 1a launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


	float McuBoundingBoxKernel(
        CUDApointers ptr,
		int nbodiesd,
		/*const realVortex* __restrict vtxd*/
        const void* __restrict vtxd,
        int sizeOfElement,
        int offsetOfPointInElement)
	{
        float time;
  time = McuBoundingBoxKernelFree(ptr.Mposl, ptr.maxrl, ptr.minrl, nbodiesd, vtxd, sizeOfElement, offsetOfPointInElement /*sizeof(realVortex), (int)realVortex::offsPos*/);
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
        const void* __restrict ptrl,
        int sizeOfElement,
        int offsetOfPointInElement,
        bool sort)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = (nbodiesd + 31) / 32;
        dim3 Mthreads = 32;
            

        MMortonCodesKernel << <Mblocks, Mthreads >> > (nbodiesd, (const void*)ptrl, sizeOfElement, offsetOfPointInElement, MmortonCodesKeyUnsortl, MmortonCodesIdxUnsortl, (const real2*)maxrl, (const real2*)minrl);

        CudaTest("MMortonCodesKernel launch failed");


        ///RadixSort
        if (sort)
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
        const void* __restrict vtxd,
        int sizeOfElement,
        int offsetOfPointInElement)
    {
        float time;
        time = McuMortonCodesKernelFree(
            ptr.maxrl, ptr.minrl, 
            ptr.MmortonCodesKeyUnsortl, ptr.MmortonCodesIdxUnsortl, ptr.MmortonCodesKeyl, ptr.MmortonCodesIdxl, 
            ptr.Mrangel, 
            nbodiesd, vtxd, sizeOfElement, offsetOfPointInElement, true);
        
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
        CudaTest("Mkernel 3a launch failed");

        RadixSortFromCUB( \
            ptr.MlevelUnsortl, ptr.MlevelSortl, \
            ptr.MindexUnsortl, ptr.MindexSortl, \
            nbodiesd-1, 0, 2 * codeLength);
        CudaTest("Mkernel 3b launch failed");

        MTransposeIndexKernel << <Mblocks, Mthreads >> > (nbodiesd, nnodesd, ptr.MindexSortl, ptr.MindexSortTl);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("Mkernel 3c launch failed");

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
        /*const realVortex* __restrict vtxd*/
        const void* __restrict ptrd,
        int sizeOfElement,
        int offsetOfPointInElement,
        bool bodyZeroSize,
        const double* XYXY
        )
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

       AABBKernel2 << <blocks * FACTOR3, THREADS3 >> > (nnodesd, nbodiesd, (int2*)Mchildd, massd, /*(real3*)vtxd,*/ ptrd, sizeOfElement, offsetOfPointInElement, MmortonCodesIdxd, (real2*)Mposd, (real2*)Mlowerd, (real2*)Mupperd, MindexSortd, MindexSortTd, bodyZeroSize, XYXY);

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
        /*const realVortex* __restrict vtxd*/
        const double* __restrict vtxd,
        int objType)
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
    float cuI1I2CalculationKernel2(
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real epssqd,
        const realVortex* __restrict vtxd,
        //realPoint* __restrict veld,
        real* __restrict I1d,
        realPoint* __restrict I2d,
        real* __restrict epsastd)
    {
        //fprintf(stderr, "FCKernel\n");

        //cudaEvent_t start, stop;
        float time = 0;

        //cudaEventCreate(&start);  cudaEventCreate(&stop);
        //cudaEventRecord(start, 0);

        I1I2CalculationKernel2 << <blocks * FACTOR5I1I2, THREADS5I1I2 >> > (
            nnodesd, nbodiesd, epssqd, (int2*)ptr.Mchildl, order, (real2*)ptr.momsl,
            (real3*)vtxd, ptr.MmortonCodesIdxl,
            (real2*)ptr.Mposl, ptr.MindexSortl, ptr.MindexSortTl,
             (real*)I1d,
            (real2*)I2d,
            (real2*)ptr.Mlowerl, (real2*)ptr.Mupperl,
            epsastd);


        //cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 5 I1I2 launch failed");

        //cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


    float cuI0I3CalculationKernel2(
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real epssqd,
        const realVortex* __restrict vtxd,
        float* __restrict I0d,
        floatPoint* __restrict I3d,
        real* __restrict epsastd, 
        const real* __restrict meanEpsd,
        int npans,
        real* __restrict pans,
        real* __restrict visstr)
    {
        float time = 0;
        
        float3* i0i3Analyt;

        cudaMalloc(&i0i3Analyt, nbodiesd * sizeof(float3));

        ClearI0I3 << <blocks * FACTOR5, THREADS5 >> > (nbodiesd, i0i3Analyt, (float*)I0d, (float2*)I3d);
        CudaTest("kernel 5a I0I3 launch failed");

        I0I3CalculationKernel2 << <blocks * FACTOR5I0I3, THREADS5I0I3 >> > (
            nnodesd, nbodiesd, epssqd, (int2*)ptr.Mchildl, order, (real2*)ptr.momsl,
            (real3*)vtxd, ptr.MmortonCodesIdxl,
            (real2*)ptr.Mposl, ptr.MindexSortl, ptr.MindexSortTl,
            (float*)I0d,
            (float2*)I3d,
            (real2*)ptr.Mlowerl, (real2*)ptr.Mupperl,
            epsastd, meanEpsd, npans, (real*)pans, (real*)visstr, (float3*)i0i3Analyt, (int*)ptr.Mrangel);
        CudaTest("kernel 5b I0I3 launch failed");

        I0I3Recalc << <blocks * FACTOR5, THREADS5 >> > (nbodiesd, (float3*)i0i3Analyt, (float*)I0d, (float2*)I3d);
        CudaTest("kernel 5c I0I3 launch failed");
       
        cudaFree(i0i3Analyt);

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

        ForceCalculationKernel2points << <blocks * FACTOR5, THREADS5W2W >> > (
            nnodesd, nbodiesd, itolsqd, epssqd, (int2*)ptr.Mchildl, order, (real2*)ptr.momsl,
            (real3*)vtxd, ptr.MmortonCodesIdxl,
            (real2*)ptr.Mposl, ptr.MindexSortl, ptr.MindexSortTl,
            npointsd, (real3*)pointsd, MmortonCodesIdxl,
            (real2*)veld, 
            (real2*)ptr.Mlowerl, (real2*)ptr.Mupperl,
            calcEpsAst, epsastd,
            nAfls, nVtxs, ptrVtxs);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 5a convvelo launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }


/******************************************************************************/
/*** compute force ************************************************************/
/******************************************************************************/
    float cuForceCalculationKernelFromPanels2points(
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real itolsqd, real epssqd,
        const double* __restrict dev_ptr_r, //начала и концы панелей
        const void* __restrict vtxd,
        int sizeOfElement,
        int offsetOfPointInElement,
        const int* __restrict MmortonCodesIdxl,
        int npointsd, const realVortex* pointsd,
        realPoint* __restrict veld, int schemeType
    )
    {
        //fprintf(stderr, "FCKernel\n");

        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        ForceCalculationKernelFromPanels2points << <blocks * FACTOR5, THREADS5S2W >> > (
            nnodesd, nbodiesd, itolsqd, epssqd, (int2*)ptr.Mchildl, order, (real2*)ptr.momsl, dev_ptr_r,
            vtxd, sizeOfElement, offsetOfPointInElement, ptr.MmortonCodesIdxl,
            (real2*)ptr.Mposl, ptr.MindexSortl, ptr.MindexSortTl,
            npointsd, (real3*)pointsd, MmortonCodesIdxl,
            (real2*)veld,
            (real2*)ptr.Mlowerl, (real2*)ptr.Mupperl, schemeType);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 5b convvelo launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }





    float McuVerticesAndSheetsToPanelPoints(int npnli, const double* dev_ptr_pt, 
        const double* dev_ptr_freeVortexSheet, const double* dev_ptr_freeVortexSheetLin,
        const double* dev_ptr_attachedVortexSheet, const double* dev_ptr_attachedVortexSheetLin,
        const double* dev_ptr_attachedSourceSheet, const double* dev_ptr_attachedSourceSheetLin,
        double* panelPoints, int schemeType)
    {
        cudaEvent_t startD, stopD;
        float timeD;

        cudaEventCreate(&startD);  cudaEventCreate(&stopD);
        cudaEventRecord(startD, 0);

        VerticesAndSheetsToPanelPointsKernel << <npnli, BLOCKD >> > (npnli, dev_ptr_pt, dev_ptr_freeVortexSheet, dev_ptr_freeVortexSheetLin, dev_ptr_attachedVortexSheet, dev_ptr_attachedVortexSheetLin, dev_ptr_attachedSourceSheet, dev_ptr_attachedSourceSheetLin, panelPoints, schemeType);
        cudaEventRecord(stopD, 0);  cudaEventSynchronize(stopD);  cudaEventElapsedTime(&timeD, startD, stopD);

        CudaTest("kernel McuVerticesAndSheetsToPanelPoints launch failed");

        cudaEventDestroy(startD);  cudaEventDestroy(stopD);

        return timeD;

    }

    float McuVerticesToControlPoints(int nTotPan, const double* dev_ptr_pt, double* pointsl)
    {
        cudaEvent_t startD, stopD;
        float timeD;

        cudaEventCreate(&startD);  cudaEventCreate(&stopD);
        cudaEventRecord(startD, 0);

        VerticesToControlPointsKernel <<<nTotPan, BLOCKD >> > (nTotPan, dev_ptr_pt, pointsl);
        cudaEventRecord(stopD, 0);  cudaEventSynchronize(stopD);  cudaEventElapsedTime(&timeD, startD, stopD);

        CudaTest("kernel McuVerticesToControlPoints launch failed");

        cudaEventDestroy(startD);  cudaEventDestroy(stopD);

        return timeD;
    }


    float McuVerticesToControlPointsVortexes(int nTotPan, const double* dev_ptr_pt, double* pointsl, const double* gaml)
    {
        cudaEvent_t startD, stopD;
        float timeD;

        cudaEventCreate(&startD);  cudaEventCreate(&stopD);
        cudaEventRecord(startD, 0);

        VerticesToControlPointsVortexesKernel << <nTotPan, BLOCKD >> > (nTotPan, dev_ptr_pt, pointsl, gaml);
        cudaEventRecord(stopD, 0);  cudaEventSynchronize(stopD);  cudaEventElapsedTime(&timeD, startD, stopD);

        CudaTest("kernel McuVerticesToControlPointsVortexes launch failed");

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
        int nTotPan, const real* dev_ptr_pt, 
        const real* pointsd,
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
            (real3*)vtxd, 
            ptr.MmortonCodesIdxl,
            (real2*)ptr.Mposl, 
            ptr.MindexSortl, ptr.MindexSortTl,
            nTotPan, (const real*)dev_ptr_pt, (const real2*)pointsd, (real2*)El, MmortonCodesIdxl,
            (real*)veld, (real*)vellind,
            (real2*)ptr.Mlowerl, (real2*)ptr.Mupperl);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel 5 rhs launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }




    float cuMatrixMulVectorCalculationKernel
    (
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real itolsqd,
        const real* __restrict rpnl,
        const real* dev_ptr_freeVortexSheet,
        const real* dev_ptr_freeVortexSheetLin,
        //CUDApointers ptrPoints,
        const int* __restrict MmortonCodesIdxl,
        realPoint* __restrict El,
        int nTotPan, const real* dev_ptr_pt,
        const real* pointsd,
        real* __restrict veld,
        real* __restrict vellind)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        MatrToVecCalculationKernel << <blocks * FACTOR5slae, THREADS5slae >> > (
            nnodesd, nbodiesd, itolsqd,
            (int2*)ptr.Mchildl, order, (real2*)ptr.momsl,
            (real*)rpnl,
            dev_ptr_freeVortexSheet,
            dev_ptr_freeVortexSheetLin,
            ptr.MmortonCodesIdxl,
            (real2*)ptr.Mposl,
            ptr.MindexSortl, ptr.MindexSortTl,
            nTotPan, (const real*)dev_ptr_pt, (const real2*)pointsd, (real2*)El, MmortonCodesIdxl,
            (real*)veld, (real*)vellind,
            (real2*)ptr.Mlowerl, (real2*)ptr.Mupperl);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTest("kernel MatrToVecCalculationKernel rhs launch failed");

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