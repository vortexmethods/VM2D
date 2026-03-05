/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: treeKernels.cu                                                   |
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
\brief Реализация функций работы с деревом на CUDA
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cuda.h>

#include "Gpudefs.h"
#include "operations.cuh"
#include "treeKernels.cuh"


namespace BHcu
{

    //////////////////////////////////////Kernels//////////////////////////////////////////////////// 

    __global__
        __launch_bounds__(THREADSgab, FACTORgab)
        void treeBoundingBoxKernel(
            unsigned int* dev_blkcnt,
            const int nbodiesd,
            const void* __restrict ptrd,
            int sizeOfElement,
            int offsetOfPointInElement,
            double2* __restrict Mposd,
            volatile double2* __restrict maxrd,
            volatile double2* __restrict minrd)
    {
        int i, j, k, inc;
        double2 val;
        double2 minr, maxr;
        __shared__ volatile double2 sminr[THREADSgab], smaxr[THREADSgab];

        // initialize with valid data (in case #bodies < #threads)
        if (sizeOfElement == 24)
        {
            double3 root = *(double3*)((char*)ptrd + 0 * sizeOfElement + offsetOfPointInElement);
            minr.x = maxr.x = root.x;
            minr.y = maxr.y = root.y;
        }
        else
        {
            minr = *(double2*)((char*)ptrd + 0 * sizeOfElement + offsetOfPointInElement);
            maxr = minr;
        }

        // scan all bodies
        i = threadIdx.x;
        inc = THREADSgab * gridDim.x;

        if (blockIdx.x == 0 && threadIdx.x == 0)
        {
            *dev_blkcnt = 0;
            __threadfence();
        }
        __syncthreads();

        for (j = i + blockIdx.x * THREADSgab; j < nbodiesd; j += inc)
        {
            if (sizeOfElement == 24)
            {
                double3 ValX = *(double3*)((char*)ptrd + j * sizeOfElement + offsetOfPointInElement);
                val.x = ValX.x;
                val.y = ValX.y;
            }
            else
                val = *(double2*)((char*)ptrd + j * sizeOfElement + offsetOfPointInElement);

            minr.x = min(minr.x, val.x);
            maxr.x = max(maxr.x, val.x);

            minr.y = min(minr.y, val.y);
            maxr.y = max(maxr.y, val.y);
        }

        // reduction in shared memory
        sminr[i].x = minr.x;
        smaxr[i].x = maxr.x;
        sminr[i].y = minr.y;
        smaxr[i].y = maxr.y;

        for (j = THREADSgab / 2; j > 0; j /= 2) {
            __syncthreads();
            if (i < j) {
                k = i + j;
                sminr[i].x = minr.x = min(minr.x, sminr[k].x);
                smaxr[i].x = maxr.x = max(maxr.x, smaxr[k].x);
                sminr[i].y = minr.y = min(minr.y, sminr[k].y);
                smaxr[i].y = maxr.y = max(maxr.y, smaxr[k].y);
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
            if (inc == atomicInc(dev_blkcnt, inc)) {

                // I'm the last block, so combine all block results
                for (j = 0; j <= inc; j++) {
                    minr.x = min(minr.x, minrd[j].x);
                    maxr.x = max(maxr.x, maxrd[j].x);
                    minr.y = min(minr.y, minrd[j].y);
                    maxr.y = max(maxr.y, maxrd[j].y);
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
    }//treeBoundingBoxKernel(...)


    __global__
        void treeMortonCodesKernel(
            const int nbodies,
            const void* __restrict ptrl,
            int sizeOfElement,
            int offsetOfPointInElement,
            int* __restrict mortonCodesKeyUnsortD,
            int* __restrict mortonCodesIdxUnsortD,
            const double2* maxrd,
            const double2* minrd)
    {
        int bdy = blockDim.x * blockIdx.x + threadIdx.x;

        register double lmax, lx0, ly0, quadSideFactor;

        if (bdy < nbodies)
        {
            double2 pmax = *maxrd;
            double2 pmin = *minrd;

            lmax = ::max(pmax.x - pmin.x, pmax.y - pmin.y);
            lx0 = (pmin.x + pmax.x) / 2; //координаты центра
            ly0 = (pmin.y + pmax.y) / 2;

            quadSideFactor = rbound / lmax; //1;

            double2 xy;
            if (sizeOfElement == 24)
            {
                double3 ptrX = *(double3*)((char*)ptrl + bdy * sizeOfElement + offsetOfPointInElement);
                xy.x = ptrX.x;
                xy.y = ptrX.y;
            }
            else
            {
                double2 ptrX = *(double2*)((char*)ptrl + bdy * sizeOfElement + offsetOfPointInElement);
                xy.x = ptrX.x;
                xy.y = ptrX.y;
            }

            double x0 = (xy.x - lx0) * quadSideFactor + rbound / 2;
            double y0 = (xy.y - ly0) * quadSideFactor + rbound / 2;

            double x = twoPowCodeLength * x0;
            double y = twoPowCodeLength * y0;

            unsigned int xx = MExpandBits((unsigned int)x);
            unsigned int yy = MExpandBits((unsigned int)y);
            mortonCodesKeyUnsortD[bdy] = yy | (xx << 1);
            mortonCodesIdxUnsortD[bdy] = bdy;           
        }
    }//treeMortonCodesKernel(...)


    __global__
        void treeMortonInternalNodesKernel(
            const int nbodies,
            const int* __restrict MmortonCodesKeyd,
            int* __restrict Mparentd,
            int2* __restrict Mchildd,
            int2* __restrict Mranged,
            int* __restrict MlevelUnsortd,
            int* __restrict MindexUnsortd)
    {
        int base = blockDim.x * blockIdx.x;
        int i = base + threadIdx.x;


        const int nLines = 9;
        const int nHalf = nLines / 2;
        __shared__ int cd[32 * nLines];

        for (int s = 0; s < nHalf; ++s)
            if (i >= (nHalf - s) * 32)
                cd[s * 32 + threadIdx.x] = MmortonCodesKeyd[i - (nHalf - s) * 32];

        if (i < nbodies)
            cd[nHalf * 32 + threadIdx.x] = MmortonCodesKeyd[i];

        for (int s = nHalf + 1; s < nLines; ++s)
            if (i + (s - nHalf) * 32 < nbodies)
                cd[s * 32 + threadIdx.x] = MmortonCodesKeyd[i + (s - nHalf) * 32];

        if (i < nbodies - 1)
        {
            int codei = cd[nHalf * 32 + threadIdx.x];

            int codeip1 = ((i + 1) >= nbodies) ? -1 : cd[nHalf * 32 + threadIdx.x + 1];
            int codeim1 = ((i - 1) < 0) ? -1 : cd[nHalf * 32 + threadIdx.x - 1];

            int Deltap1 = Delta(codei, codeip1, i, i + 1, nbodies);
            int Deltam1 = Delta(codei, codeim1, i, i - 1, nbodies);

            int d = sign(Deltap1 - Deltam1);

            int delta_min = (d > 0) ? Deltam1 : Deltap1;

            int Lmax = 2;
            int pos = i + Lmax * d;

            bool cond = ((pos >= base - nHalf * 32) && (pos < base + (nHalf + 1) * 32)) ? (Delta(codei, cd[nHalf * 32 + threadIdx.x + Lmax * d], i, pos, nbodies) > delta_min) :
                (Delta(codei, i, pos, nbodies, MmortonCodesKeyd) > delta_min);

            while (cond)
            {
                Lmax *= 2;
                pos = i + Lmax * d;

                cond = ((pos >= base - nHalf * 32) && (pos < base + (nHalf + 1) * 32)) ? (Delta(codei, cd[nHalf * 32 + threadIdx.x + Lmax * d], i, pos, nbodies) > delta_min) :
                    (Delta(codei, i, pos, nbodies, MmortonCodesKeyd) > delta_min);
            }

            int L = 0;
            for (int t = (Lmax >> 1); t >= 1; t >>= 1)
            {
                pos = i + (L + t) * d;

                cond = ((pos >= base - nHalf * 32) && (pos < base + (nHalf + 1) * 32)) ? (Delta(codei, cd[nHalf * 32 + threadIdx.x + (L + t) * d], i, pos, nbodies) > delta_min) :
                    (Delta(codei, i, pos, nbodies, MmortonCodesKeyd) > delta_min);

                if (cond)
                    L += t;
            }

            int j = i + L * d;
            pos = j;

            int delta_node = ((pos >= base - nHalf * 32) && (pos < base + (nHalf + 1) * 32)) ?
                Delta(codei, cd[nHalf * 32 + threadIdx.x + L * d], i, pos, nbodies) :
                Delta(codei, i, j, nbodies, MmortonCodesKeyd);

            MlevelUnsortd[i] = delta_node;
            MindexUnsortd[i] = i;

            int s = 0;
            for (int p = 1, t = ceilhalf(L); L > (1 << (p - 1)); ++p, t = ceilpow2(L, p))
            {
                pos = i + (s + t) * d;
               
                int dl = ((pos >= base - nHalf * 32) && (pos < base + (nHalf + 1) * 32)) ?
                    Delta(codei, cd[nHalf * 32 + threadIdx.x + (s + t) * d], i, pos, nbodies) :
                    Delta(codei, i, pos, nbodies, MmortonCodesKeyd);

                if (dl > delta_node)
                    s += t;
            }//for p

            int gamma = i + s * d + d * (d < 0);   //последнее слагаемое = std::min(d, 0);

            int Mmin = ::min(i, j);
            int Mmax = ::max(i, j);

            int left = gamma;
            int right = gamma + 1;

            // Левый потомок - лист или внутренний узел
            int childLeft = (Mmin == gamma) * nbodies + left;
            int2 pair = make_int2(Mmin, gamma);

            Mranged[childLeft] = pair;
            Mparentd[childLeft] = i;

            // Правый потомок - лист или внутренний узел
            int childRight = (Mmax == gamma + 1) * nbodies + right;

            pair = make_int2(gamma + 1, Mmax);

            Mranged[childRight] = pair;
            Mparentd[childRight] = i;

            Mchildd[i] = make_int2(childLeft, childRight);
        }
    }//treeMortonInternalNodesKernel(...)


    __global__
        void treeTransposeIndexKernel(
            const int nbodiesd,
            const int* __restrict MindexSortd,
            int* __restrict MindexSortTd)
    {
        register const int cell = blockDim.x * blockIdx.x + threadIdx.x;

        if (cell < nbodiesd - 1)
        {
            register const int newcell = MindexSortd[cell];
            MindexSortTd[newcell] = cell;
        }
    }//treeTransposeIndexKernel(...)


    __global__
        __launch_bounds__(1024, 1)
        void treeClearKernel(
            const int nnodesd, const int nbodiesd,
            volatile int* __restrict massd)
    {
        register int k, inc, bottom;

        bottom = nnodesd - (nbodiesd - 1);
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
    }//treeClearKernel(...)


    __global__
        __launch_bounds__(THREADSupward, FACTORupward)
        void treeSummarizationKernel12(
            const int nnodesd, const int nbodiesd,
            const int2* __restrict Mchildd,
            volatile int* __restrict massd,
            const int order, double2* __restrict momsd,  //momsd  - без volatile
            const double* __restrict vtxd, 
            object_T objectType, tree_T treeType,
            const int* __restrict MmortonCodesIdxd,
            double2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
            bool calcAABB,
            int sizeOfElement,
            int offsetOfPointInElement,
            double4* __restrict Mlowerupperd,
            const double4* __restrict XYXY,
            tree_T fromVortexOrSource,
            scheme_T constOrLin
        )
    {
        register int i, j, k, ch, inc, flag;

        register double4 lowerupper[2];

        register double2 mom0;
        register double2 mom1;
        register double2 mom2;
        register double2 mom3;
        register double2 mom4;
        register double2 mom5;
        register double2 mom6;
        register double2 mom7;
        register double2 mom8;
        register double2 mom9;
        register double2 mom10;
        register double2 mom11;

        register double2 cen, dr;

        register int cm;
        int m[2];


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
                register int srt = MindexSortd[(nnodesd - 1) - k];
                register int2 chdPair = Mchildd[srt];

                int chdSorted[2];

#pragma unroll
                for (i = 0; i < 2; i++) {
                    int chd = i * chdPair.y + (1 - i) * chdPair.x;   // i==0 => .x;  i==1 => .y

                    ch = (chd >= nbodiesd) ? chd - nbodiesd : (chdSorted[i] = (nnodesd - 1) - MindexSortTd[chd]);

                    if ((chd >= nbodiesd) || (massd[nnodesd - 1 - ch] >= 0))
                        j--;
                }

                if (j == 0)
                {
                    // all children are ready
                    const int kch = ((nnodesd - 1) - k) * orderAlignment;
                    cm = 0;
#pragma unroll
                    for (i = 0; i < 2; i++)
                    {
                        const int chd = i * chdPair.y + (1 - i) * chdPair.x;
                        if (chd >= nbodiesd)
                        {
                            ch = chd - nbodiesd;
                            const register int sortedBody = MmortonCodesIdxd[ch];

                            double4 xyAB = XYXY[sortedBody];
                            lowerupper[i] = make_double4(
                                ::fmin(xyAB.x, xyAB.z), ::fmin(xyAB.y, xyAB.w),
                                ::fmax(xyAB.x, xyAB.z), ::fmax(xyAB.y, xyAB.w)
                            );     
                        }
                        else
                        {
                            ch = (nnodesd - 1) - chdSorted[i];
                            lowerupper[i] = Mlowerupperd[chd];
                        }
                    }//for i

                    const double4 loup = make_double4(
                        ::fmin(lowerupper[0].x, lowerupper[1].x), ::fmin(lowerupper[0].y, lowerupper[1].y),
                        ::fmax(lowerupper[0].z, lowerupper[1].z), ::fmax(lowerupper[0].w, lowerupper[1].w)
                    );


                    Mlowerupperd[srt] = loup;
                    cen = Mposd[srt] = make_double2(0.5 * (loup.x + loup.z), 0.5 * (loup.y + loup.w));

                    const double2 zero = make_double2(0.0, 0.0);
                    register double2 momh0 = zero;
                    register double2 momh1 = zero;
                    register double2 momh2 = zero;
                    register double2 momh3 = zero;
                    register double2 momh4 = zero;
                    register double2 momh5 = zero;
                    register double2 momh6 = zero;
                    register double2 momh7 = zero;
                    register double2 momh8 = zero;
                    register double2 momh9 = zero;
                    register double2 momh10 = zero;
                    register double2 momh11 = zero;

                    for (i = 0; i < 2; i++)
                    {
                        //computation of ch = child[k*2+i]
                        const int chd = i * chdPair.y + (1 - i) * chdPair.x;
                        if (chd >= nbodiesd)
                        {
                            ch = chd - nbodiesd;
                            const register int sortedBody = MmortonCodesIdxd[ch];

                            if (objectType == object_T::point3)
                            {
                                mom0 = double2{ vtxd[sortedBody * 3 + 2], 0.0 };
                                mom1 = mom2 = mom3 = mom4 = mom5 = mom6 = mom7 = mom8 = mom9 = mom10 = mom11 = double2{ 0.0, 0.0 };
                                double3 pos = *(double3*)(&vtxd[sortedBody * 3]);
                                dr = make_double2(pos.x, pos.y) - cen;
                                m[i] = 1;
                            } //objectType==0
                            
                            if (objectType == object_T::panel)
                            {
                                double2 panBegin, panEnd;
                                panBegin = *(double2*)(&vtxd[sortedBody * 12 + 2]);
                                panEnd = *(double2*)(&vtxd[sortedBody * 12 + 4]);

                                double2 rcur, rd2Pow;
                                rcur = rd2Pow = multz(0.5 * (panEnd - panBegin), 0.5 * (panEnd - panBegin));
                                double gam;
                                double2 gammVort;
                                switch (fromVortexOrSource)
                                {
                                case tree_T::vortex:
                                    gammVort = *(double2*)(&vtxd[sortedBody * 12 + 6]);
                                    gam = gammVort.x + gammVort.y;
                                    break;

                                case tree_T::source:
                                    gam = vtxd[sortedBody * 12 + 8];
                                    break;
                                };

                                mom1 = mom3 = mom5 = mom7 = mom9 = mom11 = zero;
                                mom0 = double2{ gam, 0.0 };

                                mom2 = (gam / 3) * rcur;
                                rcur = multz(rcur, rd2Pow);
                                mom4 = (gam / 5) * rcur;
                                rcur = multz(rcur, rd2Pow);
                                mom6 = (gam / 7) * rcur;
                                rcur = multz(rcur, rd2Pow);
                                mom8 = (gam / 9) * rcur;
                                rcur = multz(rcur, rd2Pow);
                                mom10 = (gam / 11) * rcur;
                                
                                if (constOrLin == scheme_T::linScheme)
                                {
                                    double gamLin;

                                    switch (fromVortexOrSource)
                                    {
                                    case tree_T::vortex:
                                        gamLin = vtxd[sortedBody * 12 + 9] + vtxd[sortedBody * 12 + 10];                               
                                        break;
                                    case tree_T::source:
                                        gamLin = vtxd[sortedBody * 12 + 11];
                                        break;
                                    };

                                    rcur = 0.5 * (panEnd - panBegin);
                                    mom1 = gamLin * (0.5 / 3) * rcur;
                                    rcur = multz(rcur, rd2Pow);
                                    mom3 = gamLin * (0.5 / 5) * rcur;
                                    rcur = multz(rcur, rd2Pow);
                                    mom5 = gamLin * (0.5 / 7) * rcur;
                                    rcur = multz(rcur, rd2Pow);
                                    mom7 = gamLin * (0.5 / 9) * rcur;
                                    rcur = multz(rcur, rd2Pow);
                                    mom9 = gamLin * (0.5 / 11) * rcur;
                                    rcur = multz(rcur, rd2Pow);
                                    mom11 = gamLin * (0.5 / 13) * rcur;
                                }

                                double2 pos = *(double2*)(&vtxd[sortedBody * 12 + 0]);
                                dr = pos - cen;
                                m[i] = 1;
                            } //objectType == panel
                        }
                        else
                        {
                            register const int srtT = MindexSortTd[chd];
                            ch = (nnodesd - 1) - srtT;

                            const int nch = srtT * orderAlignment;

                            mom0 = double2{ momsd[nch + 0].x, (double)0 };
                            mom1 = momsd[nch + 1];
                            mom2 = momsd[nch + 2];
                            mom3 = momsd[nch + 3];
                            mom4 = momsd[nch + 4];
                            mom5 = momsd[nch + 5];
                            mom6 = momsd[nch + 6];
                            mom7 = momsd[nch + 7];
                            mom8 = momsd[nch + 8];
                            mom9 = momsd[nch + 9];
                            mom10 = momsd[nch + 10];
                            mom11 = momsd[nch + 11];

                            dr = Mposd[chd] - cen;
                            m[i] = massd[nnodesd - 1 - ch];
                        }
                        // add child's contribution

                        momh0 += mom0;
                        momh1 += mom1;
                        momh2 += mom2;
                        momh3 += mom3;
                        momh4 += mom4;
                        momh5 += mom5;
                        momh6 += mom6;
                        momh7 += mom7;
                        momh8 += mom8;
                        momh9 += mom9;
                        momh10 += mom10;
                        momh11 += mom11;

                        double2 z = dr;

                        momh1 += multz(mom0, z);
                        momh2 += 2 * multz(mom1, z);
                        momh3 += 3 * multz(mom2, z);
                        momh4 += 4 * multz(mom3, z);
                        momh5 += 5 * multz(mom4, z);
                        momh6 += 6 * multz(mom5, z);
                        momh7 += 7 * multz(mom6, z);
                        momh8 += 8 * multz(mom7, z);
                        momh9 += 9 * multz(mom8, z);
                        momh10 += 10 * multz(mom9, z);
                        momh11 += 11 * multz(mom10, z);

                        z = multz(z, dr);

                        momh2 += multz(mom0, z);
                        momh3 += 3 * multz(mom1, z);
                        momh4 += 6.0 * multz(mom2, z);
                        momh5 += 10.0 * multz(mom3, z);
                        momh6 += 15.0 * multz(mom4, z);
                        momh7 += 21.0 * multz(mom5, z);
                        momh8 += 28.0 * multz(mom6, z);
                        momh9 += 36.0 * multz(mom7, z);
                        momh10 += 45.0 * multz(mom8, z);
                        momh11 += 55.0 * multz(mom9, z);

                        z = multz(z, dr);

                        momh3 += multz(mom0, z);
                        momh4 += 4 * multz(mom1, z);
                        momh5 += 10.0 * multz(mom2, z);
                        momh6 += 20.0 * multz(mom3, z);
                        momh7 += 35.0 * multz(mom4, z);
                        momh8 += 56.0 * multz(mom5, z);
                        momh9 += 84.0 * multz(mom6, z);
                        momh10 += 120.0 * multz(mom7, z);
                        momh11 += 165.0 * multz(mom8, z);

                        z = multz(z, dr);

                        momh4 += multz(mom0, z);
                        momh5 += 5 * multz(mom1, z);
                        momh6 += 15.0 * multz(mom2, z);
                        momh7 += 35.0 * multz(mom3, z);
                        momh8 += 70.0 * multz(mom4, z);
                        momh9 += 126.0 * multz(mom5, z);
                        momh10 += 210.0 * multz(mom6, z);
                        momh11 += 330.0 * multz(mom7, z);

                        z = multz(z, dr);

                        momh5 += multz(mom0, z);
                        momh6 += 6 * multz(mom1, z);
                        momh7 += 21.0 * multz(mom2, z);
                        momh8 += 56.0 * multz(mom3, z);
                        momh9 += 126.0 * multz(mom4, z);
                        momh10 += 252.0 * multz(mom5, z);
                        momh11 += 462.0 * multz(mom6, z);

                        z = multz(z, dr);

                        momh6 += multz(mom0, z);
                        momh7 += 7 * multz(mom1, z);
                        momh8 += 28.0 * multz(mom2, z);
                        momh9 += 84.0 * multz(mom3, z);
                        momh10 += 210.0 * multz(mom4, z);
                        momh11 += 462.0 * multz(mom5, z);

                        z = multz(z, dr);

                        momh7 += multz(mom0, z);
                        momh8 += 8 * multz(mom1, z);
                        momh9 += 36.0 * multz(mom2, z);
                        momh10 += 120.0 * multz(mom3, z);
                        momh11 += 330.0 * multz(mom4, z);

                        z = multz(z, dr);

                        momh8 += multz(mom0, z);
                        momh9 += 9 * multz(mom1, z);
                        momh10 += 45.0 * multz(mom2, z);
                        momh11 += 165.0 * multz(mom3, z);

                        z = multz(z, dr);

                        momh9 += multz(mom0, z);
                        momh10 += 10 * multz(mom1, z);
                        momh11 += 55.0 * multz(mom2, z);

                        z = multz(z, dr);

                        momh10 += multz(mom0, z);
                        momh11 += 11 * multz(mom1, z);

                        z = multz(z, dr);

                        momh11 += multz(mom0, z);
                    }

                    momh0.y = 0;
                    momsd[kch + 0] = momh0;
                    momsd[kch + 1] = momh1;
                    momsd[kch + 2] = momh2;
                    momsd[kch + 3] = momh3;
                    momsd[kch + 4] = momh4;
                    momsd[kch + 5] = momh5;
                    momsd[kch + 6] = momh6;
                    momsd[kch + 7] = momh7;
                    momsd[kch + 8] = momh8;
                    momsd[kch + 9] = momh9;
                    momsd[kch + 10] = momh10;
                    momsd[kch + 11] = momh11;

                    cm += m[0] + m[1];

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
    }//treeSummarizationKernel12(...)


    __global__
        __launch_bounds__(THREADSupward, FACTORupward)
        void treeCalcAABBKernel(
            const int nnodesd, const int nbodiesd,
            const int2* __restrict Mchildd,
            volatile int* __restrict massd,
            const void* __restrict ptrd,
            int sizeOfElement,
            int offsetOfPointInElement,
            const int* __restrict MmortonCodesIdxd,
            double2* __restrict Mposd,
            double4* __restrict Mlowerupperd,
            const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
            const double4* __restrict XYXY
        )
    {
        register int i, j, k, ch, inc, flag;

        register double4 lowerupper[2];

        register int cm;
        int m[2];

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

                            double4 xyAB = XYXY[sortedBody];
                            lowerupper[i] = make_double4(
                                ::fmin(xyAB.x, xyAB.z), ::fmin(xyAB.y, xyAB.w),
                                ::fmax(xyAB.x, xyAB.z), ::fmax(xyAB.y, xyAB.w)
                            );

                            m[i] = 1;
                        }
                        else
                        {
                            register const int srtT = MindexSortTd[chd];
                            ch = (nnodesd - 1) - srtT;
                            lowerupper[i] = Mlowerupperd[chd];
                            m[i] = massd[nnodesd - 1 - ch];
                        }
                        // add child's contribution
                    }


                    const int kchSort = MindexSortd[kch];
                    const double4 loup = double4{
                        ::fmin(lowerupper[0].x, lowerupper[1].x), ::fmin(lowerupper[0].y, lowerupper[1].y),
                        ::fmax(lowerupper[0].z, lowerupper[1].z), ::fmax(lowerupper[0].w, lowerupper[1].w)
                    };

                    Mlowerupperd[kchSort] = loup;
                    Mposd[kchSort] = make_double2(0.5 * (loup.x + loup.z), 0.5 * (loup.y + loup.w));

                    cm += (m[0] + m[1]);

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
    }//treeCalcAABBKernel(...)


    /******************************************************************************/
    /*** compute force ************************************************************/
    /******************************************************************************/

    


    


    __global__
        __launch_bounds__(THREADSpanToPoint, FACTORpanToPoint)
        void treePanelsToPointsCalculationKernel(
            const int nnodesd, const int nbodiesd,
            const double itolsqd, const double epssqd,
            const int2* __restrict Mchildd,
            const int order, const double2* __restrict momsd,
            const double* __restrict dev_ptr_r, //начала и концы панелей
            const void* __restrict ptrd,
            int sizeOfElement,
            int offsetOfPointInElement,
            const int* __restrict MmortonCodesIdxd,
            const double2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
            const double4* __restrict Mlowerupperd,
            const int npointsd, const double3* __restrict pointsd,
            const int* __restrict MmortonCodesIdxPointsd,
            double2* __restrict veld,
            scheme_T schemeType, tree_T vortexOrSource
        )
    {
        register int j, k, n, depth, base, sbase, pd, nd;
        register double2 p, v, dr, ps;
        register double r2;
        register const double2* mom;

        register double2 th;


        __shared__ volatile int pos[MAXDEPTH * THREADSpanToPoint / WARPSIZE], node[MAXDEPTH * THREADSpanToPoint / WARPSIZE];

        // figure out first thread in each warp (lane 0)
        base = threadIdx.x / WARPSIZE;
        sbase = base * WARPSIZE;
        j = base * MAXDEPTH;

        __syncthreads();
        __threadfence_block();

        // iterate over all bodies assigned to thread
        for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npointsd; k += blockDim.x * gridDim.x)
        {
            const int indexOfPoint = MmortonCodesIdxPointsd[k]; //k;//MmortonCodesIdxd[k];
            p = double2{ pointsd[indexOfPoint].x, pointsd[indexOfPoint].y };

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

                register double gm, gmFree, gmAttVortex, gmAttSource;
                register double2 sumSide2;
                bool isVortex;

                while (pd < 2)
                {
                    // node on top of stack has more children to process

                    // load child pointer
                    //computation of n = childd[nd + pd] (pd = 0 или pd = 1)
                    int chd = pd * chBoth.y + (1 - pd) * chBoth.x;
                    ++pd;

                    isVortex = (chd >= nbodiesd);

                    double2 panBegin, panEnd;

                    if (isVortex)
                    {
                        n = chd - nbodiesd;
                        double* ptrValX = (double*)((char*)ptrd + MmortonCodesIdxd[n] * sizeOfElement + offsetOfPointInElement);
                        ps = double2{ *(ptrValX + 0), *(ptrValX + 1) };
                        panBegin = double2{ dev_ptr_r[MmortonCodesIdxd[n] * 4 + 0], dev_ptr_r[MmortonCodesIdxd[n] * 4 + 1] };
                        panEnd = double2{ dev_ptr_r[MmortonCodesIdxd[n] * 4 + 2], dev_ptr_r[MmortonCodesIdxd[n] * 4 + 3] };

                        gmFree = *(ptrValX + 6);
                        gmAttVortex = *(ptrValX + 7);
                        gmAttSource = *(ptrValX + 8);

                        sumSide2 = double2{ 0.0, 0.0 };
                    }
                    else
                    {
                        register const int srtT = MindexSortTd[chd];
                        n = (nnodesd - 1) - srtT;
                        ps = Mposd[chd];
                        mom = momsd + (srtT * orderAlignment);
                        gm = mom[0].x;
                        double4 gab = Mlowerupperd[chd];
                        sumSide2 = make_double2(gab.z - gab.x, gab.w - gab.y); /*Msized[chd]*/;
                    }

                    dr = p - ps;
                    r2 = (dr.x * dr.x + dr.y * dr.y);   // compute distance squared

                    // check if all threads agree that cell is far enough away (or is a body)
                    if (isVortex || __all_sync(0xffffffff, ((sumSide2.x + sumSide2.y) * (sumSide2.x + sumSide2.y) + epssqd) * itolsqd < r2))
                    {

                        if (isVortex)
                        {
                            double panLen = sqrt((panEnd.x - panBegin.x) * (panEnd.x - panBegin.x) + (panEnd.y - panBegin.y) * (panEnd.y - panBegin.y));
                            double2 vecS = p - panBegin;
                            double2 vecP = p - panEnd;
                            double alpha = atan2(-vecP.y * vecS.x + vecP.x * vecS.y, vecP.x * vecS.x + vecP.y * vecS.y);
                            double2 u0 = double2{ (panEnd.x - panBegin.x) / panLen, (panEnd.y - panBegin.y) / panLen };
                            double lambda = 0.5 * log((vecS.x * vecS.x + vecS.y * vecS.y) / (vecP.x * vecP.x + vecP.y * vecP.y));
                            double2 skos0 = double2{ alpha * u0.y + lambda * u0.x, -alpha * u0.x + lambda * u0.y };

                            switch (vortexOrSource)
                            {
                            case tree_T::vortex:
                                v += ((gmFree + gmAttVortex) / panLen) * skos0;
                                break;
                            case tree_T::source:
                                v += (gmAttSource / panLen) * skos0;
                                break;
                            };

                            if (schemeType == scheme_T::linScheme)
                            {
                                n = chd - nbodiesd;
                                double* ptrValX = (double*)((char*)ptrd + MmortonCodesIdxd[n] * sizeOfElement + offsetOfPointInElement);
                                double gmLinFree = *(ptrValX + 9);
                                double gmLinAttVortex = *(ptrValX + 10);
                                double gmLinAttSource = *(ptrValX + 11);
                                double2 vecA = vecP + vecS;
                                double2 vecB = u0;
                                double2 omega = (vecA.x * vecB.x + vecA.y * vecB.y) * vecB + double2{ vecA.y * vecB.x * vecB.y - vecA.x * vecB.y * vecB.y, -vecA.y * vecB.x * vecB.x + vecA.x * vecB.x * vecB.y }; // (vecA^ vecB)* (vecB.kcross());
                                double2 u1 = (0.5 / panLen) * omega;
                                double2 skos1 = double2{ alpha * u1.y + lambda * u1.x - u0.x, -alpha * u1.x + lambda * u1.y - u0.y };


                                switch (vortexOrSource)
                                {
                                case tree_T::vortex:
                                    v += ((gmLinFree + gmLinAttVortex) / panLen) * skos1;
                                    break;

                                case tree_T::source:
                                    v += (gmLinAttSource / panLen) * skos1;
                                    break;
                                };
                            }
                        }
                        else
                        {
                            double f = gm / max(r2, epssqd);
                            v += f * dr;
                        }

                        if ((!isVortex) && (order > 1))
                        {
                            double2 cftr = (r2 ? (1.0 / r2) : 0.0) * dr;
                            th = cftr;

                            for (int s = 1; s < order; ++s)
                            {
                                th = multz(th, cftr);
                                v += multzA(th, mom[s]);
                            }
                        }
                    }
                    else
                    {
                        if (depth < MAXDEPTH * THREADSpanToPoint / WARPSIZE)
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

            // update velocity
            double2 result = (vortexOrSource == tree_T::source) ? double2{ idpid * v.x, idpid * v.y } : double2{ -idpid * v.y, idpid * v.x };

            veld[indexOfPoint] = result;


            __syncthreads();

        }
    }//treePanelsToPointsCalculationKernel(...)


    //__device__ const float machineEps = 2.0e-7f;
    __device__ const float onePlusMachineEps = 1.0f + 2.0e-7f;

    template <typename T>
    __device__
        __forceinline void swap_device(T& a, T& b) { T c(a); a = b; b = c; }

    __device__
        __forceinline double mindist2(const double4& lowerUpper, const double2& rhs) noexcept
    {
        const double dx = ::fmin(lowerUpper.z, ::fmax(lowerUpper.x, rhs.x)) - rhs.x;
        const double dy = ::fmin(lowerUpper.w, ::fmax(lowerUpper.y, rhs.y)) - rhs.y;
        return dx * dx + dy * dy;
    }

    __device__
        __forceinline double minmaxdist2(const double4& lowerUpper, const double2& rhs) noexcept
    {
        double2 rm_sq = make_double2((lowerUpper.x - rhs.x) * (lowerUpper.x - rhs.x),
            (lowerUpper.y - rhs.y) * (lowerUpper.y - rhs.y));
        double2 rM_sq = make_double2((lowerUpper.z - rhs.x) * (lowerUpper.z - rhs.x),
            (lowerUpper.w - rhs.y) * (lowerUpper.w - rhs.y));

        if ((lowerUpper.z + lowerUpper.x) * 0.5 < rhs.x)
            swap_device(rm_sq.x, rM_sq.x);

        if ((lowerUpper.w + lowerUpper.y) * 0.5 < rhs.y)
            swap_device(rm_sq.y, rM_sq.y);

        const double dx = rm_sq.x + rM_sq.y;
        const double dy = rM_sq.x + rm_sq.y;

        return ::fmin(dx, dy);
    }

    __device__
        __forceinline double2 distance_calculator_point2segment(const double2& point, const double4& object)
    {
        int location;

        double a = object.w - object.y;
        double b = object.z - object.x;

        double distanceSegment;

        double2 dr = make_double2(point.x - object.x, point.y - object.y);

        double r_numerator = dr.x * b + dr.y * a;
        double r_denomenator = b * b + a * a;
        double r = r_numerator / r_denomenator;

        double s = (dr.x * a - dr.y * b);

        if ((r >= 0) && (r <= 1))
        {
            distanceSegment = s * s / r_denomenator;
            location = 0;
        }
        else
        {
            double dist1 = dr.x * dr.x + dr.y * dr.y;
            double dist2 = (point.x - object.z) * (point.x - object.z) + (point.y - object.w) * (point.y - object.w);
            if (dist1 < dist2)
            {
                distanceSegment = dist1;
                location = -1;
            }
            else
            {
                distanceSegment = dist2;
                location = 1;
            }
        }

        return make_double2(distanceSegment, location);
    }



    __global__
        __launch_bounds__(THREADSnear, FACTORnear)
        void treeClosestPanelToPointsCalculationKernel(
            const int nnodesd, const int nbodiesd,
            const int2* __restrict Mchildd,
            const double* __restrict dev_ptr_r, //начала и концы панелей
            const int* __restrict MmortonCodesIdxd,
            const double4* __restrict Mlowerupperd,
            const int npointsd, const double2* __restrict pointsd,
            const int* __restrict MmortonCodesIdxPointsd,
            int* __restrict near,            
            bool findOnlyInside,
            double* pseudonormals)
    {        
        int k = threadIdx.x + blockIdx.x * blockDim.x;
        {
            if (k >= npointsd)
                return;

            int idx = MmortonCodesIdxPointsd[k];
            double2 query = pointsd[idx];
            double2 stack[32];
            double2* stack_ptr = stack;

            double4 gab0 = Mlowerupperd[0];

            double md2 = mindist2(gab0, query);
            int2 nearestLocation{ -1, 0 };
            if ((findOnlyInside) && (md2 != 0))
            {
                near[idx] = -1;
                return;
            }

            *stack_ptr++ = make_double2(__longlong_as_double(0), md2);

            double dist_to_nearest_object = 1e+100;
            double4 pnlLeftGab, pnlRightGab;            

            do
            {                
                const double2 node = *--stack_ptr;
                if (node.y > dist_to_nearest_object)
                {
                    // if aabb mindist > already_found_mindist, it cannot have a nearest
                    continue;
                }

                const int2 LR_idx = Mchildd[__double_as_longlong(node.x)];

                bool isLeftLeaf = (LR_idx.x >= nbodiesd);
                bool isRightLeaf = (LR_idx.y >= nbodiesd);               

                double4 L_box, R_box;
                if (!isLeftLeaf)
                    L_box = Mlowerupperd[LR_idx.x];
                else
                {
                    int n = LR_idx.x - nbodiesd;
                    pnlLeftGab = *(double4*)(dev_ptr_r + MmortonCodesIdxd[n] * 4);
                    L_box = make_double4(
                        ::fmin(pnlLeftGab.x, pnlLeftGab.z), ::fmin(pnlLeftGab.y, pnlLeftGab.w),
                        ::fmax(pnlLeftGab.x, pnlLeftGab.z), ::fmax(pnlLeftGab.y, pnlLeftGab.w)
                    );
                }

                if (!isRightLeaf)
                    R_box = Mlowerupperd[LR_idx.y];
                else
                {
                    int n = LR_idx.y - nbodiesd;
                    pnlRightGab = *(double4*)(dev_ptr_r + MmortonCodesIdxd[n] * 4);
                    R_box = make_double4(
                        ::fmin(pnlRightGab.x, pnlRightGab.z), ::fmin(pnlRightGab.y, pnlRightGab.w),
                        ::fmax(pnlRightGab.x, pnlRightGab.z), ::fmax(pnlRightGab.y, pnlRightGab.w)
                    );
                }

                const double L_mindist2 = mindist2(L_box, query);
                const double R_mindist2 = mindist2(R_box, query);

                const double L_minmaxdist2 = minmaxdist2(L_box, query);
                const double R_minmaxdist2 = minmaxdist2(R_box, query);

                if (L_mindist2 <= R_minmaxdist2 * onePlusMachineEps) // L is worth considering
                {
                    if (isLeftLeaf) // leaf node
                    {
                        int n = LR_idx.x - nbodiesd;                       

                        double2 dist_code = distance_calculator_point2segment(query, pnlLeftGab);
                        if (dist_code.x <= dist_to_nearest_object)
                        {
                            dist_to_nearest_object = dist_code.x;
                            nearestLocation = make_int2(MmortonCodesIdxd[n], (int)dist_code.y);
                        }
                    }
                    else
                    {
                        *stack_ptr++ = make_double2(__longlong_as_double(LR_idx.x), L_mindist2);
                    }
                }

                if (R_mindist2 <= L_minmaxdist2 * onePlusMachineEps) // R is worth considering
                {                    
                    if (isRightLeaf) // leaf node
                    {
                        int n = LR_idx.y - nbodiesd;

                        double2 dist_code = distance_calculator_point2segment(query, pnlRightGab);
                        if (dist_code.x <= dist_to_nearest_object)
                        {
                            dist_to_nearest_object = dist_code.x;
                            nearestLocation = make_int2(MmortonCodesIdxd[n], (int)dist_code.y);
                        }
                    }
                    else
                    {                        
                        *stack_ptr++ = make_double2(__longlong_as_double(LR_idx.y), R_mindist2);
                    }
                }                
            } while (stack < stack_ptr); 
            
            //Проверка внутри-снаружи
            if (findOnlyInside)
            {
                int pnl = nearestLocation.x;
                int code = nearestLocation.y;
                double2 nrm = *(double2*)&(pseudonormals[6 * pnl + (code + 1) * 2]);
                double2 pt;
                switch (code)
                {
                case -1:
                    pt = *(double2*)&(dev_ptr_r[4 * pnl]);
                    break;                
                case 1:
                    pt = *(double2*)&(dev_ptr_r[4 * pnl + 2]);
                    break;
                default:
                    pt = 0.5 * (*(double2*)&(dev_ptr_r[4 * pnl]) + *(double2*)&(dev_ptr_r[4 * pnl + 2]));
                }
                pt -= query;
                if ((pt.x * nrm.x + pt.y * nrm.y) < 0) //значит, что снаружи
                    nearestLocation.x = -1;
            }
            
            near[idx] = nearestLocation.x;
        }
    }


    __host__ __device__
    bool ifAabbContainsSegment(const double2 start, const double2 finish, const double4& object) noexcept
    {
        // Completely outside.
        if ((start.x <= object.x && finish.x <= object.x) || (start.y <= object.y && finish.y <= object.y) ||
            (start.x >= object.z && finish.x >= object.z) || (start.y >= object.w && finish.y >= object.w))
            return false;

        double m = (finish.y - start.y) / (finish.x - start.x);

        double y = m * (object.x - start.x) + start.y;
        if (y > object.y && y < object.w) return true;

        y = m * (object.z - start.x) + start.y;
        if (y > object.y && y < object.z) return true;

        double x = (object.y - start.y) / m + start.x;
        if (x > object.x && x < object.z) return true;

        x = (object.w - start.y) / m + start.x;
        if (x > object.x && x < object.z) return true;

        // Start or end inside.
        if ((start.x > object.x && start.x < object.z && start.y > object.y && start.y < object.w)
            || (finish.x > object.x && finish.x < object.z && finish.y > object.y && finish.y < object.w)) return true;

        return false;
    }


    __device__ int getLineIntersection(double2 p0, double2 p1, double4 p23, double2& rint)

    {
        double2 s1, s2;
        s1 = p1 - p0;
        s2 = make_double2(p23.z - p23.x, p23.w - p23.y);

        double s, t;
        double den = (-s2.x * s1.y + s1.x * s2.y);

        double dx = (p0.x - p23.x) / den;
        double dy = (p0.y - p23.y) / den;

        s = (-s1.y * dx + s1.x * dy);
        t = (s2.x * dy - s2.y * dx);

        if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
        {
            rint = p0 + t * s1;
            return 1;
        }

        return 0; // No collision
    }


    //определяет расстояние по лучу до точки пересечения этим лучом отрезка object (или infinity)
    __device__ double  getRayDistanceToIntersectionLine(const double2 start, const double2 finish, const double4& object)
    {
        double2 rint;
        double infin = 1e+100;

        int intersect = getLineIntersection(start, finish, object, rint);
        if (intersect)
            return (rint.x - start.x) * (rint.x - start.x) + (rint.y - start.y) * (rint.y - start.y);
        else
            return infin;
    }

    //определяет расстояние по лучу до точки пересечения этим лучом прямоугольника (или infinity)
    __device__ double getRayDistanceToIntersectionRectangle(const double2& start, const double2& finish, const double4& object)
    {
        if ((object.x <= start.x) && (start.x <= object.z) && (object.y <= start.y) && (start.y <= object.w) &&
            (object.x <= finish.x) && (finish.x <= object.z) && (object.y <= finish.y) && (finish.y <= object.w))
            return 0;

        double d1, d2, d3, d4;
        double2 rint;
        double infin = 1e+100;
        double4 object1, object2, object3, object4;
        object1 = make_double4(object.z, object.w, object.z, object.y);
        object2 = make_double4(object.z, object.y, object.x, object.y);
        object3 = make_double4(object.x, object.y, object.x, object.w);
        object4 = make_double4(object.x, object.w, object.z, object.w);

        int intersect1 = getLineIntersection(start, finish, object1, rint);
        d1 = (intersect1) ? (rint.x - start.x) * (rint.y - start.x) + (rint.y - start.y) * (rint.y - start.y) : infin;

        int intersect2 = getLineIntersection(start, finish, object2, rint);
        d2 = (intersect2) ? (rint.x - start.x) * (rint.x - start.x) + (rint.y - start.y) * (rint.y - start.y) : infin;

        int intersect3 = getLineIntersection(start, finish, object3, rint);
        d3 = (intersect3) ? (rint.x - start.x) * (rint.x - start.x) + (rint.y - start.y) * (rint.y - start.y) : infin;

        int intersect4 = getLineIntersection(start, finish, object4, rint);
        d4 = (intersect4) ? (rint.x - start.x) * (rint.x - start.x) + (rint.y - start.y) * (rint.y - start.y) : infin;

        d1 = ::min(d1, d2);
        d1 = ::min(d1, d3);
        d1 = ::min(d1, d4);
        return d1;
    };
    
    __global__
    __launch_bounds__(THREADSsegintersect, FACTORsegintersect)
    void treePanelsSegmentsIntersectionCalculationKernel(
        const int nnodesd, const int nbodiesd,
        const int2* __restrict Mchildd,
        const double* __restrict dev_ptr_r, //начала и концы панелей
        const int* __restrict MmortonCodesIdxd,
        const double4* __restrict Mlowerupperd,
        const int nrays, const double4* __restrict q,
        const int* __restrict MmortonCodesIdxPointsd,
        int* __restrict near
        )
    {
        double infin = 1e+100;
        unsigned int hitted = 0xFFFFFFFF;
        double dist_to_hitted_object = infin;

        int k = threadIdx.x + blockIdx.x * blockDim.x;
        if (k >= nrays)
            return;

        double2 stack[32];
        double2* stack_ptr = stack;
        double4 gab0 = Mlowerupperd[0];

        int idx = MmortonCodesIdxPointsd[k];
        
        double4 ray = q[idx];

        bool crossRoot = ifAabbContainsSegment(make_double2(ray.x, ray.y), make_double2(ray.z, ray.w), gab0);
        if (crossRoot)
            *stack_ptr++ = make_double2(__longlong_as_double(0), 0);

        while (stack < stack_ptr)
        {
            const auto node = *--stack_ptr;
            if (node.y > dist_to_hitted_object)
                continue;

            const auto obj_idx = __double_as_longlong(node.x);

            if (obj_idx >= nbodiesd)
            {

                int n = obj_idx - nbodiesd;
                double4 pnlRightGab = *(double4*)(dev_ptr_r + MmortonCodesIdxd[n] * 4);
                double hit = getRayDistanceToIntersectionLine(make_double2(ray.x, ray.y), make_double2(ray.z, ray.w), pnlRightGab); // check_hit_obj - проверка пересечения лучом отрезка

                if (hit < dist_to_hitted_object)
                {
                    dist_to_hitted_object = hit;
                    hitted = obj_idx;
                }
            }
            else
            {
                const auto L_idx = Mchildd[__double_as_longlong(node.x)].x;
                const auto R_idx = Mchildd[__double_as_longlong(node.x)].y;

                double4 L_box, R_box;
                if (L_idx < nbodiesd)
                {
                    L_box = Mlowerupperd[L_idx];
                }
                else
                {
                    int n = L_idx - nbodiesd;
                    double4 pnlRightGab = *(double4*)(dev_ptr_r + MmortonCodesIdxd[n] * 4);
                    L_box = make_double4(
                        ::fmin(pnlRightGab.x, pnlRightGab.z), ::fmin(pnlRightGab.y, pnlRightGab.w),
                        ::fmax(pnlRightGab.x, pnlRightGab.z), ::fmax(pnlRightGab.y, pnlRightGab.w)
                    );
                }

                if (R_idx < nbodiesd)
                {
                    R_box = Mlowerupperd[R_idx];
                }
                else
                {
                    int n = R_idx - nbodiesd;
                    double4 pnlRightGab = *(double4*)(dev_ptr_r + MmortonCodesIdxd[n] * 4);
                    R_box = make_double4(
                        ::fmin(pnlRightGab.x, pnlRightGab.z), ::fmin(pnlRightGab.y, pnlRightGab.w),
                        ::fmax(pnlRightGab.x, pnlRightGab.z), ::fmax(pnlRightGab.y, pnlRightGab.w)
                    );
                }

                auto h_L = getRayDistanceToIntersectionRectangle(make_double2(ray.x, ray.y), make_double2(ray.z, ray.w), L_box);  ///////// check_hit_node - проверка пересечения лучом прямоугольника
                auto h_R = getRayDistanceToIntersectionRectangle(make_double2(ray.x, ray.y), make_double2(ray.z, ray.w), R_box);

                if ((h_L) == infin && (h_R) == infin)
                {
                }
                else if ((h_L) < infin && (h_R) == infin)
                {                    
                    *stack_ptr++ = make_double2(__longlong_as_double(L_idx), h_L);
                }
                else if ((h_R) < infin && (h_L) == infin)
                {                    
                    *stack_ptr++ = make_double2(__longlong_as_double(R_idx), h_R);
                }
                else
                {    
                    if (h_L > h_R)
                    {
                        *stack_ptr++ = make_double2(__longlong_as_double(L_idx), h_L);
                        *stack_ptr++ = make_double2(__longlong_as_double(R_idx), h_R);
                    }
                    else
                    {
                        *stack_ptr++ = make_double2(__longlong_as_double(R_idx), h_R);
                        *stack_ptr++ = make_double2(__longlong_as_double(L_idx), h_L);
                    }
                }
            }
        }

        near[idx] = (hitted == 0xFFFFFFFF) ? -1 : MmortonCodesIdxd[hitted - nbodiesd];
    }

    // Функция обхода дерева для подсчета ячеек в дальней зоне (без вычислений!)
    __global__
        __launch_bounds__(THREADSslae, FACTORslae)
        void treeMatrToVecNoCalculationKernel
        (
            const int nnodesd,  // Число узлов дерева (по вихрям)
            const int nbodiesd, // Число вихрей в пелене

            const double itolsqd, // 1/theta^2
            const int2* __restrict Mchildd, //Массив потомков узлов дерева

            const double* __restrict rpnl,  //начала и концы влияющих панелей

            const int* __restrict MmortonCodesIdxd, //порядок сортировки влияющих панелей

            const double2* __restrict Mposd, //массив координат центров ячеек дерева
            const int* __restrict MindexSortd, //порядок сортировки внутренних узлов дерева
            const int* __restrict MindexSortTd,//обратный порядок сортировки внутренних узлов дерева
            const double4* __restrict Mlowerupperd,//координаты нижнего левого и верхнего правого угла ячейки дерева

            const int npointsd,                  //количество панелей на теле
            const double* __restrict dev_ptr_pt,   //массив четверок (beg.x, beg.y, end.x, end.y)
            const int* __restrict MmortonCodesIdxPointsd,//порядок сортировки центров контрольных панелей

            int* __restrict nClosePanelsl,
            int* __restrict nFarCellsl
        )
    {
        register int j, k, n, depth, base, sbase, pd, nd;
        register double2 p, dr, ps;
        register double r2;

        __shared__ volatile int pos[MAXDEPTH * THREADSslae / WARPSIZE], node[MAXDEPTH * THREADSslae / WARPSIZE];

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
            
            nClosePanelsl[indexOfPoint] = 0;
            nFarCellsl[indexOfPoint] = 0;

            double2 beg, end;
            beg = double2{ dev_ptr_pt[4 * indexOfPoint + 0], dev_ptr_pt[4 * indexOfPoint + 1] };
            end = double2{ dev_ptr_pt[4 * indexOfPoint + 2], dev_ptr_pt[4 * indexOfPoint + 3] };
            p = 0.5 * (beg + end);

            double dlen2 = (end.x - beg.x) * (end.x - beg.x) + (end.y - beg.y) * (end.y - beg.y);


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

                register double2 sumSide2;
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

                        ps = double2{ 0.5 * (rpnl[MmortonCodesIdxd[n] * 4 + 0] + rpnl[MmortonCodesIdxd[n] * 4 + 2]), 0.5 * (rpnl[MmortonCodesIdxd[n] * 4 + 1] + rpnl[MmortonCodesIdxd[n] * 4 + 3]) };

                        sumSide2 = double2{ fabs(rpnl[MmortonCodesIdxd[n] * 4 + 2] - rpnl[MmortonCodesIdxd[n] * 4 + 0]),
                            fabs(rpnl[MmortonCodesIdxd[n] * 4 + 3] - rpnl[MmortonCodesIdxd[n] * 4 + 1]) };
                    }
                    else
                    {
                        register const int srtT = MindexSortTd[chd];
                        n = (nnodesd - 1) - srtT;
                        ps = Mposd[chd];
                        double4 gab = Mlowerupperd[chd];
                        sumSide2 = make_double2(gab.z - gab.x, gab.w - gab.y); /*Msized[chd]*/;
                    }

                    //ps - положение вихря/кластера
                    //p - центр панели
                    dr = p - ps;
                    r2 = (dr.x * dr.x + dr.y * dr.y);   // compute distance squared               

                    // check if all threads agree that cell is far enough away (or is a body)
                    if (isVortex || __all_sync(0xffffffff, ((sumSide2.x + sumSide2.y) * (sumSide2.x + sumSide2.y) + dlen2) * itolsqd < r2))
                    {
                        if (isVortex)
                            ++nClosePanelsl[indexOfPoint];
                        else
                            ++nFarCellsl[indexOfPoint];                        
                    }
                    else //идем глубже по дереву
                    {
                        if (depth < MAXDEPTH * THREADSslae / WARPSIZE)
                        {

                            // push cell onto stack
                            if (sbase == threadIdx.x)
                            {  // maybe don't push and inc if last child
                                if (depth >= MAXDEPTH * THREADSslae / WARPSIZE)
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
        }
    }//treeMatrToVecNoCalculationKernel(...)



    ////////// Auxilary functions for SLAE //////////

    /// \brief Вспомогательная функция, которая определяет, находится ли панель itI на контуре после панели itJ
    ///
    /// \param[in] itI константная ссылка на обертку для второй ("правой") панели
    /// \param[in] itJ константная ссылка на обертку для первой ("левой") панели
    __device__ __forceinline bool isEqual(const double2& paniBeg, const double2& panjEnd)
    {
        //return ((paniBeg.x - panjEnd.x) * (paniBeg.x - panjEnd.x) + (paniBeg.y - panjEnd.y) * (paniBeg.y - panjEnd.y) < 1e-20);
        return (fabs(paniBeg.x - panjEnd.x) + fabs(paniBeg.y - panjEnd.y) < 1e-10);
    }

    /// \brief Вспомогательная функция вычисления угла между векторами (в диапазоне (-pi...pi]) (со знаком, поворот от первого ко второму)
    ///
    /// \param[in] p константная ссылка на первый вектор
    /// \param[in] s константная ссылка на второй вектор 
    __device__ __forceinline double Alpha(const double2& p, const double2& s)
    {
        return atan2(p.x * s.y - p.y * s.x, p.x * s.x + p.y * s.y);
    }

    /// \brief  Вспомогательная функция вычисления логарифма отношения норм векторов
    /// 
    /// \param[in] p константная ссылка на первый вектор
    /// \param[in] s константная ссылка на второй вектор 
    __device__ __forceinline double Lambda(const double2& p, const double2& s)
    {
        return 0.5 * log((s.x * s.x + s.y * s.y) / (p.x * p.x + p.y * p.y));
    }


    /// Вспомогательная функция вычисления величины \f$ (\vec a \cdot \vec b) \cdot \vec c + (\vec a \times \vec b) \times \vec c \f$
    __device__ __forceinline double2 Omega(const double2& a, const double2& b, const double2& c)
    {
        double adotb = a.x * b.x + a.y * b.y;
        double acrossb = a.x * b.y - a.y * b.x;
        return double2{ adotb * c.x - acrossb * c.y, adotb * c.y + acrossb * c.x };
    };

    ////////// Auxilary functions for SLAE //////////




    template <int order>
    __global__
        __launch_bounds__(THREADSslae, FACTORslae)
        void treeMatrToVecZeroIterKernel12
        (
            const int nnodesd,  // Число узлов дерева (по вихрям)
            const int nbodiesd, // Число вихрей в пелене

            const double itolsqd, // 1/theta^2
            const int2* __restrict Mchildd, //Массив потомков узлов дерева

            const double2* __restrict momsd, //моменты узлов дерева
            const double* __restrict rpnl,  //начала и концы влияющих панелей
            const double* dev_pnl,
            /*const double* dev_ptr_freeVortexSheet,
            const double* dev_ptr_freeVortexSheetLin,*/

            const int* __restrict MmortonCodesIdxd, //порядок сортировки вихрей в пелене

            const double2* __restrict Mposd, //массив координат центров ячеек дерева
            const int* __restrict MindexSortd, //порядок сортировки внутренних узлов дерева
            const int* __restrict MindexSortTd,//обратный порядок сортировки внутренних узлов дерева
            const double4* __restrict Mlowerupperd,//координаты нижнего левого и верхнего правого угла ячейки дерева
            scheme_T schemeType,

            const int npointsd,                  //количество панелей на теле
            const double* __restrict dev_ptr_pt,   //массив четверок (beg.x, beg.y, end.x, end.y)

            double2* __restrict Ed,                //коэффициенты локального разложения
            const int* __restrict MmortonCodesIdxPointsd,//порядок сортировки центров панелей

            double* __restrict veld,               //куда сохранять ответ
            double* __restrict vellind,            //куда сохранять ответ

            int* __restrict closeCellsPfl,
            int* __restrict ClosePrefixSuml,
            int* __restrict farCellsPfl,
            int* __restrict FarPrefixSuml,
            double2* __restrict i00save,
            double2* __restrict i01save,
            double2* __restrict i10save,
            double2* __restrict i11save, int iter
        )
    {
        const int WIDTH5 = min(blockDim.x, WARPSIZE);

        register int base = threadIdx.x / WIDTH5;
        register int sbase = base * WIDTH5;
        register int j = base * MAXDEPTH;
        register int diff = threadIdx.x - sbase;
        register int depth;

        int npointsUp = ((npointsd + WIDTH5 - 1) / WIDTH5) * WIDTH5;

        register int k;
        register int nd;
        register int2 chBoth;
        register int pd;

        register const double2* mom;
        register int indexOfPoint;
        register int srtT;

        register double2 p;
        register double val, vallin;

        bool isVortex;
        register double2 ps;
        register double gm;

        register double sumSide2;
        register double2 dr;
        register double r2;

        register double2 beg, end;
        register double dlen2;
        register double2 tau;
        register double2 di;
        register double dilen, idlen;

        register double4 gab;

        register double2 Eloc[order];

        int closecntr, farcntr;

        __shared__ volatile int pos[MAXDEPTH * maxTHREADSslae / WARPSIZE], node[MAXDEPTH * maxTHREADSslae / WARPSIZE];
        __shared__ double2 momSh[std::max((maxTHREADSslae / WARPSIZE), 1) * orderAlignment];

        __syncthreads();
        __threadfence_block();

        // iterate over all bodies assigned to thread
        for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npointsUp; k += blockDim.x * gridDim.x)
        {
            if (k < npointsd)
            {
                indexOfPoint = MmortonCodesIdxPointsd[k];

                closecntr = 0;
                farcntr = 0;

                const double4 pnl = *(double4*)&(dev_ptr_pt[4 * indexOfPoint]);
                beg = make_double2(pnl.x, pnl.y);
                end = make_double2(pnl.z, pnl.w);
                p = 0.5 * (beg + end);

                di = make_double2(end.x - beg.x, end.y - beg.y);
                dlen2 = di.x * di.x + di.y * di.y;

                idlen = rsqrt(dlen2);
                dilen = 1.0 / idlen;
                tau = (end - beg) * idlen;

                val = 0;
                vallin = 0;
            }

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

                chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];

                register int mcn;

                while (pd < 2)
                {
                    // node on top of stack has more children to process

                    // load child pointer
                    //computation of n = childd[nd + pd] (pd = 0 или pd = 1)
                    register int chd = pd * chBoth.y + (1 - pd) * chBoth.x;

                    ++pd;

                    isVortex = (chd >= nbodiesd);

                    register int n;


                    if (isVortex)
                    {
                        n = chd - nbodiesd;

                        mcn = MmortonCodesIdxd[n];

                        gab = *(double4*)&(rpnl[mcn * 4]);
                        const double2 pnlVec = make_double2(gab.z - gab.x, gab.w - gab.y);
                        ps = double2{ 0.5 * (gab.x + gab.z), 0.5 * (gab.y + gab.w) };

                        double lj = sqrt(pnlVec.x * pnlVec.x + pnlVec.y * pnlVec.y);

                        sumSide2 = fabs(pnlVec.x) + fabs(pnlVec.y);
                        sumSide2 *= sumSide2;
                    }
                    else
                    {
                        srtT = MindexSortTd[chd];
                        n = (nnodesd - 1) - srtT;
                        ps = Mposd[chd];

                        gab = Mlowerupperd[chd];
                        sumSide2 = gab.z - gab.x + gab.w - gab.y;
                        sumSide2 *= sumSide2;
                    }
                
                    if (k >= npointsd)                
                        sumSide2 = dlen2 = 0.0;
                    
                    //ps - положение вихря/кластера
                    //p - центр панели
                    dr = p - ps;
                    r2 = (k < npointsd) ? (dr.x * dr.x + dr.y * dr.y) : 1000000.0;   // compute distance squared               
 
                    // check if all threads agree that cell is far enough away (or is a body)
                    if (isVortex || __all_sync(0xffffffff, (sumSide2 + dlen2) * itolsqd < r2))
                    {
                        if (!isVortex)
                        {
                            mom = momsd + (srtT * orderAlignment);
                            if (diff < order)
                                momSh[base * (orderAlignment)+diff] = mom[diff];
                        }
                        __threadfence_block();
                        __syncthreads();//added

                        if (k < npointsd)
                        {
                            if (isVortex)
                            {
                                closeCellsPfl[ClosePrefixSuml[indexOfPoint] + closecntr] = chd;

                                double tempVelNew;
                                const double2 ptpanBeg = make_double2(gab.x, gab.y);
                                const double2 ptpanEnd = make_double2(gab.z, gab.w);

                                const double2 ptpanVec = make_double2(ptpanEnd.x - ptpanBeg.x, ptpanEnd.y - ptpanBeg.y);
                                const double invlj = rsqrt(ptpanVec.x * ptpanVec.x + ptpanVec.y * ptpanVec.y);
                                const double lj = 1 / invlj;

                                gm = lj * dev_pnl[12 * mcn + 6];//dev_ptr_freeVortexSheet[mcn];
                                double gmlin = 0.0;
                                if (schemeType == scheme_T::linScheme)//if (dev_ptr_freeVortexSheetLin != nullptr)
                                    gmlin = lj * dev_pnl[12 * mcn + 9];//dev_ptr_freeVortexSheetLin[mcn];

                                if (r2 > 1e-20)
                                {
                                    double2 i00;
                                    double2 dj;
                                    double djlen;
                                    double ilenj;
                                    double2 tauj;
                                    double2 p1, s1, p2, s2;
                                    double3 alpha, lambda;
                                    bool condBefore, condAfter;
                                    if (iter == 0)
                                    {
                                        dj = make_double2(ptpanEnd.x - ptpanBeg.x, ptpanEnd.y - ptpanBeg.y);
                                        ilenj = rsqrt(dj.x * dj.x + dj.y * dj.y);
                                        djlen = 1 / ilenj;
                                        tauj = make_double2(ilenj * dj.x, ilenj * dj.y);

                                        p1 = make_double2(end.x - ptpanEnd.x, end.y - ptpanEnd.y);
                                        s1 = make_double2(end.x - ptpanBeg.x, end.y - ptpanBeg.y);
                                        p2 = make_double2(beg.x - ptpanEnd.x, beg.y - ptpanEnd.y);
                                        s2 = make_double2(beg.x - ptpanBeg.x, beg.y - ptpanBeg.y);

                                        condBefore = isEqual(beg, ptpanEnd);
                                        condAfter = isEqual(ptpanBeg, end);

                                        alpha = make_double3(\
                                            condAfter ? 0.0 : Alpha(s2, s1), \
                                            Alpha(s2, p1), \
                                            condBefore ? 0.0 : Alpha(p1, p2) \
                                        );

                                        lambda = make_double3(\
                                            condAfter ? 0.0 : Lambda(s2, s1), \
                                            Lambda(s2, p1), \
                                            condBefore ? 0.0 : Lambda(p1, p2) \
                                        );

                                        double2 v00_0 = Omega(s1, tau, tauj);
                                        double2 v00_1 = Omega(di, tau, tauj);

                                        double2 v00_2 = Omega(p2, tau, tauj);

                                        i00 = make_double2(
                                            ilenj * (alpha.x * v00_0.x - alpha.y * v00_1.x + alpha.z * v00_2.x \
                                                - (lambda.x * v00_0.y - lambda.y * v00_1.y + lambda.z * v00_2.y)),
                                            ilenj * (alpha.x * v00_0.y - alpha.y * v00_1.y + alpha.z * v00_2.y \
                                                + (lambda.x * v00_0.x - lambda.y * v00_1.x + lambda.z * v00_2.x))
                                        );

                                        i00save[ClosePrefixSuml[indexOfPoint] + closecntr] = i00;
                                    }
                                    else 
                                        i00 = i00save[ClosePrefixSuml[indexOfPoint] + closecntr];

                                    tempVelNew = -gm * (i00.x * tau.x + i00.y * tau.y);                                    
                                    val -= tempVelNew;

                                    if (schemeType == scheme_T::linScheme)
                                    {
                                        double2 i01, i10, i11;

                                        if (iter == 0)
                                        {
                                            double s1len2 = s1.x * s1.x + s1.y * s1.y;

                                            double2 om1 = Omega(s1, tau, tauj);
                                            double2 om2 = Omega(make_double2(s1.x + p2.x, s1.y + p2.y), tauj, tauj);
                                            double sc = (p1.x + s1.x) * tauj.x + (p1.y + s1.y) * tauj.y;

                                            const double hLj = 0.5 * ilenj;
                                            const double pA = hLj * sc;
                                            const double pB = hLj * s1len2;
                                            const double hLij = hLj * dilen;

                                            double2 v01_0 = make_double2(pA * om1.x - pB * tau.x, pA * om1.y - pB * tau.y);
                                            double2 v01_1 = make_double2(hLij * om2.x, hLij * om2.y);

                                            i01 = make_double2(\
                                                ilenj * ((alpha.x + alpha.z) * v01_0.x - (alpha.y + alpha.z) * v01_1.x\
                                                    - (((lambda.x + lambda.z) * v01_0.y - (lambda.y + lambda.z) * v01_1.y) - 0.5 * dilen * tauj.y)),
                                                ilenj * ((alpha.x + alpha.z) * v01_0.y - (alpha.y + alpha.z) * v01_1.y\
                                                    + (((lambda.x + lambda.z) * v01_0.x - (lambda.y + lambda.z) * v01_1.x) - 0.5 * dilen * tauj.x))
                                            );

                                            i01save[ClosePrefixSuml[indexOfPoint] + closecntr] = i01;

                                            double2 om3 = Omega(make_double2(s1.x + p2.x, s1.y + p2.y), tau, tau);

                                            const double hLi = 0.5 * idlen;
                                            const double pC = hLi * ((s1 + s2) & tau);
                                            const double pD = hLi * s1len2;
                                            const double hLji = hLi * djlen;

                                            double2 v10_0 = make_double2(pC * om1.x - pD * tauj.x, pC * om1.y - pD * tauj.y);
                                            double2 v10_1 = make_double2(hLji * om3.x, hLji * om3.y);

                                            i10 = make_double2(
                                                ilenj * (-(alpha.x + alpha.z) * v10_0.x + alpha.z * v10_1.x \
                                                    - ((-(lambda.x + lambda.z) * v10_0.y + lambda.z * v10_1.y) + 0.5 * djlen * tau.y)),
                                                ilenj * (-(alpha.x + alpha.z) * v10_0.y + alpha.z * v10_1.y \
                                                    + ((-(lambda.x + lambda.z) * v10_0.x + lambda.z * v10_1.x) + 0.5 * djlen * tau.x))
                                            );

                                            //11-01-2026
                                            i10save[ClosePrefixSuml[indexOfPoint] + closecntr] = i10;

                                            double2 om4 = Omega(make_double2(s1.x - 3.0 * p2.x, s1.y - 3.0 * p2.y), tau, tauj);
                                            double2 om5 = Omega(di, tauj, tauj);
                                            double2 om6 = Omega(dj, tau, tau);
                                            double sc4 = s1.x * om4.x + s1.y * om4.y;

                                            const double mn = idlen * ilenj / 12.0;
                                            const double qA = mn * 2.0 * sc4 - 0.25;
                                            const double qB = mn * s1len2;

                                            double2 v11_0 = make_double2(
                                                qA * om1.x - qB * (s1.x - 3 * p2.x),
                                                qA * om1.y - qB * (s1.y - 3 * p2.y)
                                            );

                                            double qC = hLij / 6.0;
                                            double qD = hLji / 6.0;

                                            double2 v11_1 = make_double2(qC * om5.x, qC * om5.y);
                                            double2 v11_2 = make_double2(qD * om6.x, qD * om6.y);

                                            i11 = make_double2(
                                                ilenj * ((alpha.x + alpha.z) * v11_0.x - (alpha.y + alpha.z) * v11_1.x - alpha.z * v11_2.x\
                                                    - ((lambda.x + lambda.z) * v11_0.y - (lambda.y + lambda.z) * v11_1.y - lambda.z * v11_2.y \
                                                        + 1.0 / 12.0 * (djlen * tau.y + dilen * tauj.y - 2.0 * om1.y))),
                                                ilenj * ((alpha.x + alpha.z) * v11_0.y - (alpha.y + alpha.z) * v11_1.y - alpha.z * v11_2.y\
                                                    + ((lambda.x + lambda.z) * v11_0.x - (lambda.y + lambda.z) * v11_1.x - lambda.z * v11_2.x \
                                                        + 1.0 / 12.0 * (djlen * tau.x + dilen * tauj.x - 2.0 * om1.x)))
                                            );

                                            i11save[ClosePrefixSuml[indexOfPoint] + closecntr] = i11;
                                        }
                                        else
                                        {
                                            i01 = i01save[ClosePrefixSuml[indexOfPoint] + closecntr];
                                            i10 = i10save[ClosePrefixSuml[indexOfPoint] + closecntr];
                                            i11 = i11save[ClosePrefixSuml[indexOfPoint] + closecntr];
                                        }
                                        double tempVelNewB = -gmlin * (i01.x * tau.x + i01.y * tau.y);
                                        val -= tempVelNewB;

                                        double tempVelNewC = -gm * (i10.x * tau.x + i10.y * tau.y);
                                        double tempVelNewD = -gmlin * (i11.x * tau.x + i11.y * tau.y);

                                        vallin -= (tempVelNewC + tempVelNewD);
                                    }
                                }//if (posI != posJ)

                                ++closecntr;

                            }
                            else
                            {
                                farCellsPfl[FarPrefixSuml[indexOfPoint] + farcntr] = chd;

                                double2 theta = (p - ps) / r2;
#pragma unroll
                                for (int q = 0; q < order; ++q)
                                    Eloc[q] = make_double2(0.0, 0.0);

#pragma unroll
                                for (int q = 0; q < order; ++q)
                                {
                                    for (int s = q; s >= 0; --s)
                                        Eloc[s] += (((2 * q - s) & 1) ? -ifac[q - s + 1] : ifac[q - s + 1]) * multzA(theta, momSh[base * (orderAlignment)+q - s]);
                                    theta = ((q + 1) / r2) * multz(theta, p - ps);
                                }

                                double2 v = 0.5 * Eloc[0];
                                double2 vL{ 0.0, 0.0 };

                                double2 rPan = end - beg;

                                double2 kp;

                                double2 mulP = kp = rPan;

                                double2 taudL = (0.5 / dlen2) * rPan;
                                double2 taudLc = taudL;
#pragma unroll
                                for (int k = 1; k < order; ++k)
                                {
                                    mulP = multz(mulP, kp);
                                    taudL *= 0.5;

                                    if (!(k & 1))
                                        v += ifac[k + 2] * multz(Eloc[k], multzA(taudL, mulP));
                                    else
                                        if (schemeType == scheme_T::linScheme)
                                            vL += ifac[k + 3] * multz(Eloc[k], multzA(multz(taudL, taudLc), multz(mulP, (k + 1) * rPan)));
                                }

                                v *= 2.0;
                                if (schemeType == scheme_T::linScheme)
                                    vL *= 2.0;

                                val += (-v.y * rPan.x + v.x * rPan.y);

                                if (schemeType == scheme_T::linScheme)
                                    vallin += (-vL.y * rPan.x + vL.x * rPan.y);

                                ++farcntr;
                            }
                        }
                        __threadfence_block();
                        __syncthreads();//added
                    }
                    else //идем глубже по дереву
                    {
                        if (depth < MAXDEPTH * maxTHREADSslae / WARPSIZE)
                        {
                            if (pd == 1)
                            {
                                if (sbase == threadIdx.x)
                                {
                                    pos[depth] = 1;
                                    node[depth] = nd;
                                }
                                ++depth;
                            }
                            pd = 0;
                            nd = n;

                            chBoth = Mchildd[MindexSortd[(nnodesd - 1) - nd]];
                        }//if depth<
                    }
                }
                --depth;

            } while (depth >= j);

            // update velocity
            if (k < npointsd)
            {
                register double cf = idpid * idlen;
                veld[indexOfPoint] = cf * val;

                if (schemeType == scheme_T::linScheme)
                    vellind[indexOfPoint] = cf * vallin;
            };
        }
    }//treeMatrToVecZeroIterKernel12(...)

    template <int order>
    __global__
        __launch_bounds__(THREADSslae, FACTORslae)
        void treeMatrToVecNonZeroIterKernel12 //new 
        (
            const int nnodesd,  // Число узлов дерева (по вихрям)
            const int nbodiesd, // Число вихрей в пелене

            const double2* __restrict momsd, //моменты узлов дерева
            const double* __restrict rpnl,  //вихри в пелене
            const double* dev_pnl,
            /*const double* dev_ptr_freeVortexSheet,
            const double* dev_ptr_freeVortexSheetLin,*/

            const int* __restrict MmortonCodesIdxd, //порядок сортировки вихрей в пелене

            const double2* __restrict Mposd, //массив координат центров ячеек дерева
            const int* __restrict MindexSortTd,//обратный порядок сортировки внутренних узлов дерева
            scheme_T schemeType,

            const int npointsd,                  //количество панелей на теле
            const double* __restrict dev_ptr_pt,   //массив четверок (beg.x, beg.y, end.x, end.y)

            double2* __restrict Ed,                //коэффициенты локального разложения
            const int* __restrict MmortonCodesIdxPointsd,//порядок сортировки центров панелей

            double* __restrict veld,               //куда сохранять ответ
            double* __restrict vellind,            //куда сохранять ответ

            int* __restrict nClosePanelsl,
            int* __restrict nFarCellsl,

            const int* __restrict closeCellsPfl,
            const int* __restrict ClosePrefixSuml,
            const int* __restrict farCellsPfl,
            const int* __restrict FarPrefixSuml,

            const double2* __restrict i00save,
            const double2* __restrict i01save,
            const double2* __restrict i10save,
            const double2* __restrict i11save
        )
    {
        const int WIDTH5 = min(blockDim.x, WARPSIZE);
        register int k, n;
        register int base = threadIdx.x / WIDTH5;
        register int sbase = base * WIDTH5;
        register int diff = threadIdx.x - sbase;
        register double2 p, ps;
        register const double2* mom;
        register double val, vallin;

        register int MmortonCodesIdxd_n = 0;

        register double2 Eloc[order];


        bool calcLin = (schemeType == scheme_T::linScheme);


        __syncthreads();
        __threadfence_block();

        // iterate over all bodies assigned to thread

        int npointsdUp = ((npointsd + 31) / 32) * 32;

        int indexOfPoint;
        double2 beg, end;
        double2 di, tau;
        double dlen2, idlen, lj, deltax, deltay;
        __shared__ int nClosePnls;
        __shared__ int nFarCells;

        for (k = threadIdx.x + blockIdx.x * blockDim.x; k < npointsdUp; k += blockDim.x * gridDim.x)
        {
            if (k < npointsd)
            {
                indexOfPoint = MmortonCodesIdxPointsd[k];

                double4 begend = *(double4*)(dev_ptr_pt + 4 * indexOfPoint);

                beg = make_double2(begend.x, begend.y);
                end = make_double2(begend.z, begend.w);
                p = 0.5 * (beg + end);

                di = make_double2(end.x - beg.x, end.y - beg.y);
                dlen2 = di.x * di.x + di.y * di.y;

                idlen = rhypot(di.x, di.y);

                tau = (end - beg) * idlen;

                lj = 0.0; // used only: if(isVortex)!

                val = 0;
                vallin = 0;

                nClosePnls = nClosePanelsl[indexOfPoint]; 
                nFarCells = nFarCellsl[indexOfPoint];
            }

            __threadfence();
            __syncthreads();

            register double gm, gmlin;


            int lock = k / 32 * 32;
            int indexOfFirstPointInWarp = MmortonCodesIdxPointsd[lock];
            int shift = ClosePrefixSuml[indexOfFirstPointInWarp];


            __shared__ double freeSheetIntensity[32];
            __shared__ double freeSheetIntensityLin[32];

            int nBlocks = (nClosePnls + 31) / 32;

            for (int blk = 0; blk < nBlocks; ++blk)
            {
                if (blk * 32 + threadIdx.x < nClosePnls)
                {
                    int chd = closeCellsPfl[shift + (blk * 32 + threadIdx.x)];
                    n = chd - nbodiesd;
                    MmortonCodesIdxd_n = MmortonCodesIdxd[n];

                    freeSheetIntensity[threadIdx.x] = dev_pnl[12 * MmortonCodesIdxd_n + 6];//dev_ptr_freeVortexSheet[MmortonCodesIdxd_n];
                    if (calcLin)
                        freeSheetIntensityLin[threadIdx.x] = dev_pnl[12 * MmortonCodesIdxd_n + 9];//dev_ptr_freeVortexSheetLin[MmortonCodesIdxd_n];
                }

                __syncthreads();

                if (k < npointsd)
                {
                    for (int q = 0; q < 32; ++q)
                    {
                        int pcntr = blk * 32 + q;
                        if (pcntr < nClosePnls)
                        {
                            int chd = closeCellsPfl[shift + pcntr];
                            n = chd - nbodiesd;

                            MmortonCodesIdxd_n = MmortonCodesIdxd[n];

                            double4 begendbuf = *(double4*)(rpnl + MmortonCodesIdxd_n * 4);
                            double2 begbuf = make_double2(begendbuf.x, begendbuf.y);
                            double2 endbuf = make_double2(begendbuf.z, begendbuf.w);

                            deltax = endbuf.x - begbuf.x;
                            deltay = endbuf.y - begbuf.y;
                            lj = hypot(deltax, deltay);

                            gm = lj * freeSheetIntensity[q]; //from sharedMem

                            double tempVelNew;

                            if (indexOfPoint != MmortonCodesIdxd_n)
                            {
                                int position = ClosePrefixSuml[indexOfPoint] + pcntr;

                                double2 i00;
                                i00 = i00save[position];

                                tempVelNew = -gm * (i00.x * tau.x + i00.y * tau.y);
                                val -= tempVelNew;

                                if (vellind)
                                {
                                    gmlin = lj * freeSheetIntensityLin[q]; //from sharedMem
                                    double2 i01, i10, i11;
                                    i01 = i01save[position];
                                    i10 = i10save[position];
                                    i11 = i11save[position];

                                    double tempVelNewB = -gmlin * (i01.x * tau.x + i01.y * tau.y);
                                    val -= tempVelNewB;


                                    double tempVelNewC = -gm * (i10.x * tau.x + i10.y * tau.y);
                                    double tempVelNewD = -gmlin * (i11.x * tau.x + i11.y * tau.y);
                                    vallin -= (tempVelNewC + tempVelNewD);
                                }
                            }//if (posI != posJ)
                        }
                    }
                }//if k

                __syncthreads();
            }
            //*/



            __shared__ double2 momSh[32 * orderAlignment];

            shift = FarPrefixSuml[indexOfFirstPointInWarp];

            nBlocks = (nFarCells + 31) / 32;

            for (int blk = 0; blk < nBlocks; ++blk)
            {
#pragma unroll
                for (int cluster = 0; cluster < 32; ++cluster)
                {
                    if (blk * 32 + cluster < nFarCells)
                    {
                        int chd = farCellsPfl[shift + (blk * 32 + cluster)];
                        register const int srtT = MindexSortTd[chd];
                        mom = momsd + (srtT * orderAlignment);
                        if (diff < order)
                            momSh[cluster * (orderAlignment)+diff] = mom[diff];
                    }
                }

                //*/
                __syncthreads();

                if (k < npointsd)
                {
#pragma unroll
                    for (int qq = 0; qq < 32; ++qq)
                    {
                        int pcntr = blk * 32 + qq;
                        if (pcntr < nFarCells)
                        {
                            int chd = farCellsPfl[shift + pcntr];
                            ps = Mposd[chd];
                            double dist2 = (p.x - ps.x) * (p.x - ps.x) + (p.y - ps.y) * (p.y - ps.y);
                            double2 theta = (p - ps) / dist2;

#pragma unroll
                            for (int q = 0; q < order; ++q)
                                Eloc[q] = make_double2(0.0, 0.0);

#pragma unroll
                            for (int q = 0; q < order; ++q)
                            {
                                for (int s = q; s >= 0; --s)
                                    Eloc[s] += (((2 * q - s) & 1) ? -ifac[q - s + 1] : ifac[q - s + 1]) * multzA(theta, momSh[qq * (orderAlignment)+q - s]);
                                theta = ((q + 1) / dist2) * multz(theta, p - ps);
                            }

                            double2 v = 0.5 * Eloc[0];
                            double2 vL{ 0.0, 0.0 };

                            double2 rPan = di;

                            double2 kp;

                            double2 mulP = kp = rPan;

                            double2 taudL = (0.5 / dlen2) * rPan;
                            double2 taudLc = taudL;
#pragma unroll
                            for (int k = 1; k < order; ++k)
                            {
                                mulP = multz(mulP, kp);
                                taudL *= 0.5;

                                if (!(k & 1))
                                    v += ifac[k + 2] * multz(Eloc[k], multzA(taudL, mulP));
                                else
                                    if (calcLin)
                                        vL += ifac[k + 3] * multz(Eloc[k], multzA(multz(taudL, taudLc), multz(mulP, (k + 1) * rPan)));
                            }
                            v *= 2.0;
                            if (calcLin)
                                vL *= 2.0;

                            val += (-v.y * rPan.x + v.x * rPan.y);

                            if (calcLin)
                                vallin += (-vL.y * rPan.x + vL.x * rPan.y);
                        }
                    }//for blk
                }//if k

                __syncthreads();
            }

            // update velocity
            if (k < npointsd)
            {
                double cft = idpid * idlen;
                veld[indexOfPoint] = cft * val;

                if (vellind)
                    vellind[indexOfPoint] = cft * vallin;
            }
        }
    }//treeMatrToVecNonZeroIterKernel12(...)

    /******************************************************************************/
    /*** compute diffusive velo ***************************************************/
    /******************************************************************************/


    __global__
        __launch_bounds__(THREADSI1I2, FACTORI1I2)
        void treeI1I2CalculationKernel(
            const int nnodesd, const int nbodiesd,
            const double minRd,
            const int2* __restrict Mchildd,
            const double2* __restrict momsd,
            const double3* __restrict vtxd,
            const int* __restrict MmortonCodesIdxd,
            const double2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
            const double4* __restrict Mlowerupperd,
            double* __restrict I1d,
            double2* __restrict I2d,
            const double* __restrict epsast)

    {
        register int j, k, n, depth, base, sbase, pd, nd;
        register double2 p, i2, dr, ps;
        register double r2sq, i1;
        register const double2* mom;

        register double rdi, diffRadius, diffRadiusMonopole;
        register double expr;

        __shared__ volatile int pos[MAXDEPTH * THREADSI1I2 / WARPSIZE], node[MAXDEPTH * THREADSI1I2 / WARPSIZE];

        // figure out first thread in each warp (lane 0)
        base = threadIdx.x / WARPSIZE;
        sbase = base * WARPSIZE;
        j = base * MAXDEPTH;

        __syncthreads();
        __threadfence_block();

        // iterate over all bodies assigned to thread
        for (k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd; k += blockDim.x * gridDim.x)
        {
            const int indexInParticles = MmortonCodesIdxd[k];

            p = double2{ vtxd[indexInParticles].x, vtxd[indexInParticles].y };

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

                register double gm;
                register double2 sumSide2;
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

                        ps = double2{ vtxd[MmortonCodesIdxd[n]].x, vtxd[MmortonCodesIdxd[n]].y };
                        gm = vtxd[MmortonCodesIdxd[n]].z;
                        sumSide2 = double2{ 0.0, 0.0 };
                    }
                    else
                    {
                        register const int srtT = MindexSortTd[chd];
                        n = (nnodesd - 1) - srtT;

                        ps = Mposd[chd];
                        mom = momsd + (srtT * orderAlignment);
                        gm = mom[0].x;

                        double4 gab = Mlowerupperd[chd];
                        sumSide2 = make_double2(gab.z - gab.x, gab.w - gab.y); /*Msized[chd]*/;
                    }

                    dr = p - ps;
                    r2sq = sqrt(dr.x * dr.x + dr.y * dr.y);   // compute distance 

                    rdi = max(epsast[indexInParticles], minRd);
                    diffRadius = 8.0 * rdi;
                    diffRadiusMonopole = 4.0 * rdi;

                    // check if all threads agree that cell is far enough away (or is a body)
                    if (__all_sync(0xffffffff, r2sq - 0.5 * (sumSide2.x + sumSide2.y) > diffRadius)) {}
                    else if (isVortex || __all_sync(0xffffffff, r2sq - 0.5 * (sumSide2.x + sumSide2.y) > diffRadiusMonopole))
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
                        if (depth < MAXDEPTH * THREADSI1I2 / WARPSIZE)
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
    }//treeI1I2CalculationKernel(...)


    __global__
        __launch_bounds__(THREADSforces, FACTORforces)
        void ClearI0I3(int nbodiesd, float* __restrict I0d, float2* __restrict I3d)
    {
        for (int k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd; k += blockDim.x * gridDim.x)
        {
            I0d[k] = 0;
            I3d[k].x = 0;
            I3d[k].y = 0;
        }
    }//ClearI0I3(...)


    __global__
        __launch_bounds__(THREADSI0I3, FACTORI0I3)
        void treeI0I3CalculationKernel(
            const int nnodesd, const int nbodiesd,
            const double minRd,
            const int2* __restrict Mchildd,
            const double2* __restrict momsd,
            const double3* __restrict vtxd,
            const int* __restrict MmortonCodesIdxd,
            const double2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd,
            const double4* __restrict Mlowerupperd, int* range,
            float* __restrict I0d,
            float2* __restrict I3d,
            const double* __restrict epsast, const double* __restrict meanEpsd, const int npan, const double* __restrict pansd, double* __restrict visstrd)
    {
        register int j, k, n, depth, base, sbase, pd, nd;
        register float2 p, q, ps;
        register float i0;
        register float2 i3;
        register const double2* mom;

        __shared__ volatile int pos[MAXDEPTH * THREADSI0I3 / WARPSIZE], node[MAXDEPTH * THREADSI0I3 / WARPSIZE];

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
                        mom = momsd + (srtT * orderAlignment);
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
                        double4 gab = Mlowerupperd[chd];
                        float2 low{ (float)gab.x, (float)gab.y };
                        float2 up{ (float)gab.z, (float)gab.w };

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
                        smax = myMax(myMax(myMax(s1, s2), s3), s4);

                        cond = (low.x < p.x) && (p.x < up.x) && (low.y < p.y) && (p.y < up.y);
                    }


                    if (__all_sync(0xffffffff, ((d > farDist * lenj) || ((fabsf(s) > farDist * lenj) && (s * smax > 0))) && !(cond))) {}
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

                            vs += gm * expon / (valPif * meanepsj2);
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

                            vsm *= gm / (valPif * meanepsj2);
                            vs += vsm;
                        }
                        else
                        {
                            i0 += -valPif * rdi;

                            float argcosh = fabsf(s) * iDDomRad;

                            float prod = (argcosh < 50.0f) ? expf(-lenj * 0.5f * iDDomRad) * coshf(fabsf(s) * iDDomRad) :
                                0.5f * expf((s - 0.5f * lenj) * iDDomRad);

                            mnog1 = 2.0f * rdi * (1.0f - prod);
                            i3.x += mnog1 * norm.x;
                            i3.y += mnog1 * norm.y;

                            vs += mnog1 * gm / (valPif * meanepsj2);
                        }

                        atomicAdd(I0d + MmortonCodesIdxd[n], i0);
                        atomicAdd((float*)I3d + 2 * MmortonCodesIdxd[n], i3.x);
                        atomicAdd((float*)I3d + 2 * MmortonCodesIdxd[n] + 1, i3.y);

                    }// if (isVortex)
                    else //идем глубже по дереву
                    {
                        if (depth < MAXDEPTH * THREADSI0I3 / WARPSIZE)
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
    }//treeI0I3CalculationKernel(...)



    //////////////////////////////////////Wrappers///////////////////////////////////////////////// 


    void CudaTestError(const char* msg)
    {
        //#ifdef DEBUG
        cudaError_t e;

        //cudaThreadSynchronize();
        cudaDeviceSynchronize();
        if (cudaSuccess != (e = cudaGetLastError())) {
            fprintf(stderr, "%s: %d\n", msg, e);
            fprintf(stderr, "%s\n", cudaGetErrorString(e));
            exit(-1);
        }
        //#endif // DEBUG
    }//CudaTestError(...)
    

    float treeBoundingBoxWrapper(CudaTreeInfo& treeInfo)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);
        unsigned int* dev_blkcntr;

        cudaMalloc(&dev_blkcntr, sizeof(unsigned int));
        treeBoundingBoxKernel << <treeInfo.nBlock * FACTORgab, THREADSgab >> > (dev_blkcntr, treeInfo.nObject, treeInfo.objectD, treeInfo.sizeOfElement, treeInfo.offsetOfPointInElement, (double2*)treeInfo.centerD, (double2*)treeInfo.maxrD, (double2*)treeInfo.minrD);
        cudaFree(dev_blkcntr);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTestError("treeBoundingBoxKernel launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }//treeBoundingBoxWrapper(...)


    float treeMortonCodesWrapper(CudaTreeInfo& treeInfo, bool sort)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = (treeInfo.nObject + 31) / 32;
        dim3 Mthreads = 32;

        treeMortonCodesKernel <<<Mblocks, Mthreads>>> (treeInfo.nObject, (const void*)treeInfo.objectD, treeInfo.sizeOfElement, \
                                                       treeInfo.offsetOfPointInElement, treeInfo.mortonCodesKeyUnsortD, treeInfo.mortonCodesIdxUnsortD, \
                                                       (const double2*)treeInfo.maxrD, (const double2*)treeInfo.minrD);

        CudaTestError("treeMortonCodesKernel launch failed");

        ///RadixSort
        if (sort)
            treeInfo.RadixSortMortonCodes();

        CudaTestError("RadixSortMortonCodes launch failed");

        //Заполнение нулевой ячейки (диапазон для корня дерева)
        int totalRange[2] = { 0, treeInfo.nObject - 1 };
        if (treeInfo.rangeD != nullptr)
            cudaMemcpy((void*)treeInfo.rangeD, &totalRange, 2 * sizeof(int), cudaMemcpyHostToDevice);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);
        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }//treeMortonCodesWrapper(...)


    float treeMortonInternalNodesWrapper(CudaTreeInfo& treeInfo, bool sort)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        dim3 Mblocks = ((treeInfo.nObject - 1) + 31) / 32;
        dim3 Mthreads = 32;

        treeMortonInternalNodesKernel << <Mblocks, Mthreads >> > (treeInfo.nObject,
            treeInfo.mortonCodesKeyD, treeInfo.parentD, (int2*)treeInfo.childD, (int2*)treeInfo.rangeD, treeInfo.levelUnsortD, treeInfo.indexUnsortD);
        CudaTestError("treeMortonInternalNodesKernel launch failed");
        if (sort)
        {
            treeInfo.RadixSortInternalCells();
            CudaTestError("RadixSortInternalCells launch failed");

            treeTransposeIndexKernel << <Mblocks, Mthreads >> > (treeInfo.nObject, treeInfo.indexSortD, treeInfo.indexSortTD);
            CudaTestError("treeTransposeIndexKernel launch failed");
        }

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);      
        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }//treeMortonInternalNodesWrapper(...)


    float treeClearKernelWrapper(CudaTreeInfo& treeInfo, const int order)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        if (order > 0)
            cudaMemset((void*)treeInfo.momsD, 0, (treeInfo.nObject - 1) * orderAlignment * sizeof(Point2D));

        treeClearKernel << < treeInfo.nBlock * 1, 1024 >> > (treeInfo.nNode, treeInfo.nObject, treeInfo.massD);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTestError("kernel clear2 launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }//treeClearKernelWrapper(...)
    

    float treeSummarizationWrapper(CudaTreeInfo& treeInfo, const int order, bool calcAABB)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        //case 12:
        treeSummarizationKernel12 << <treeInfo.nBlock * FACTORupward, THREADSupward >> > (\
            treeInfo.nNode, treeInfo.nObject, (int2*)treeInfo.childD, treeInfo.massD, order, (double2*)treeInfo.momsD, (double*)treeInfo.objectD, \
            treeInfo.objectType, treeInfo.treeType, treeInfo.mortonCodesIdxD, (double2*)treeInfo.centerD, treeInfo.indexSortD, treeInfo.indexSortTD, \
            true, treeInfo.sizeOfElement, treeInfo.offsetOfPointInElement, (double4*)treeInfo.lowerupperD, (double4*)treeInfo.gabForLeavesD, treeInfo.treeType, treeInfo.schemeType);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTestError("treeSummarizationKernel12 launch failed");
        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }//treeSummarizationWrapper(...)


    float treeCalcAABBWrapper(CudaTreeInfo& treeInfo)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        treeCalcAABBKernel << <treeInfo.nBlock * FACTORupward, THREADSupward >> > (treeInfo.nNode, treeInfo.nObject, (int2*)treeInfo.childD, treeInfo.massD, \
            treeInfo.objectD, treeInfo.sizeOfElement, treeInfo.offsetOfPointInElement, treeInfo.mortonCodesIdxD, (double2*)treeInfo.centerD, (double4*)treeInfo.lowerupperD, \
            treeInfo.indexSortD, treeInfo.indexSortTD, (double4*)treeInfo.gabForLeavesD);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTestError("treeCalcAABBKernel launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }//treeCalcAABBWrapper(...)

    ////DELETE
    //float treeVorticesToPointsCalculationWrapper(
    //    const CudaTreeInfo& treeVorticesInfo, const CudaTreeInfo& controlTreeInfo,
    //    const double itolsq, const double epssq,
    //    Point2D* __restrict velD, bool calcEpsAst, double* __restrict epsastD)
    //{
    //    cudaEvent_t start, stop;
    //    float time;

    //    cudaEventCreate(&start);  cudaEventCreate(&stop);
    //    cudaEventRecord(start, 0);

    //    treeVorticesToPointsCalculationKernel<12> << <treeVorticesInfo.nBlock * FACTORforces, THREADSforces >> > (
    //        treeVorticesInfo.nNode, treeVorticesInfo.nObject, itolsq, epssq, (int2*)treeVorticesInfo.childD, (double2*)treeVorticesInfo.momsD, \
    //        (double3*)treeVorticesInfo.objectD, treeVorticesInfo.mortonCodesIdxD, (double2*)treeVorticesInfo.centerD, treeVorticesInfo.indexSortD, \
    //        treeVorticesInfo.indexSortTD, (double4*)treeVorticesInfo.lowerupperD, \
    //        controlTreeInfo.nObject, (double3*)controlTreeInfo.objectD, controlTreeInfo.mortonCodesIdxD, (double2*)velD, calcEpsAst, epsastD);

    //    cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

    //    CudaTestError("treeVorticesInfluenceCalculationKernel launch failed");

    //    cudaEventDestroy(start);  cudaEventDestroy(stop);
    //    return time;
    //}//treeVorticesInfluenceCalculationWrapper(...)


    float treePanelsToPointsCalculationWrapper(const CudaTreeInfo& treePanelsInfo, const CudaTreeInfo& controlTreeInfo,
        int order, double itolsq2, const double epssq, Point2D* __restrict velD)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        treePanelsToPointsCalculationKernel << <treePanelsInfo.nBlock * FACTORpanToPoint, THREADSpanToPoint >> > (
            treePanelsInfo.nNode, treePanelsInfo.nObject, itolsq2, epssq, (int2*)treePanelsInfo.childD, order, (double2*)treePanelsInfo.momsD, \
            (double*)treePanelsInfo.gabForLeavesD, treePanelsInfo.objectD, treePanelsInfo.sizeOfElement, treePanelsInfo.offsetOfPointInElement, \
            treePanelsInfo.mortonCodesIdxD, (double2*)treePanelsInfo.centerD, treePanelsInfo.indexSortD, treePanelsInfo.indexSortTD, \
            (double4*)treePanelsInfo.lowerupperD, controlTreeInfo.nObject, (double3*)controlTreeInfo.objectD, controlTreeInfo.mortonCodesIdxD, (double2*)velD, \
            treePanelsInfo.schemeType, treePanelsInfo.treeType);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTestError("treePanelsToPointsCalculationKernel launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }//treePanelsToPointsCalculationWrapper(...)



    float treeClosestPanelToPointsCalculationWrapper(const CudaTreeInfo& treePanelsInfo, const CudaTreeInfo& controlTreeInfo, 
        int* __restrict closePnlD, bool findOnlyInside, double* pseudoNormals)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        treeClosestPanelToPointsCalculationKernel<<<(controlTreeInfo.nObject + THREADSnear - 1) / THREADSnear, THREADSnear>>>(
                treePanelsInfo.nNode,
                treePanelsInfo.nObject,
                (int2*)treePanelsInfo.childD,
                (double*)treePanelsInfo.gabForLeavesD,
                treePanelsInfo.mortonCodesIdxD,
                (double4*)treePanelsInfo.lowerupperD,

                controlTreeInfo.nObject,
                (double2*)controlTreeInfo.objectD,
                controlTreeInfo.mortonCodesIdxD,
                closePnlD,
                findOnlyInside,
                pseudoNormals
            );

        CudaTestError("treeClosestPanelToPointsCalculationKernel launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }//treePanelsToPointsCalculationWrapper(...)



    float treePanelsSegmentsIntersectionCalculationWrapper(
        const CudaTreeInfo& treePanelsInfo, const CudaTreeInfo& controlTreeInfo,
        int* __restrict closePnlD)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);
        
        treePanelsSegmentsIntersectionCalculationKernel << <(controlTreeInfo.nObject + THREADSsegintersect - 1) / THREADSsegintersect, THREADSsegintersect >>> (
            treePanelsInfo.nNode, 
            treePanelsInfo.nObject, 
            (int2*)treePanelsInfo.childD, 
            (double*)treePanelsInfo.gabForLeavesD, 
            treePanelsInfo.mortonCodesIdxD,
            (double4*)treePanelsInfo.lowerupperD,
            
            controlTreeInfo.nObject,
            (double4*)controlTreeInfo.gabForLeavesD,
            controlTreeInfo.mortonCodesIdxD,
            closePnlD
            );

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTestError("treePanelsSegmentsIntersectionCalculationKernel launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }//treePanelsSegmentsIntersectionCalculationWrapper(...)


    float treeMatrToVecNoCalculationWrapper(CudaTreeInfo& treePanelsInfo, double itolsq)
    {
        cudaEvent_t start, stop;
        float time;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        treeMatrToVecNoCalculationKernel << <treePanelsInfo.nBlock * FACTORslae, THREADSslae >> > (
            treePanelsInfo.nNode, treePanelsInfo.nObject, itolsq,
            (int2*)treePanelsInfo.childD, 
            (double*)treePanelsInfo.gabForLeavesD,
            treePanelsInfo.mortonCodesIdxD,
            (double2*)treePanelsInfo.centerD,
            treePanelsInfo.indexSortD, treePanelsInfo.indexSortTD, (double4*)treePanelsInfo.lowerupperD,
            treePanelsInfo.nObject,
            (double*)treePanelsInfo.gabForLeavesD,
            treePanelsInfo.mortonCodesIdxD,
            treePanelsInfo.matVecMulInfo.nClosePanelsD, treePanelsInfo.matVecMulInfo.nFarCellsD);

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        CudaTestError("kernel treeMatrToVecNoCalculationKernel launch failed");

        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }//treeMatrToVecNoCalculationWrapper(...)

    float treeMatrToVecWrapper(CudaTreeInfo& treePanelsInfo, int order, double itolsq, double* __restrict resD, double* __restrict resLinD, int iter)
    {
        cudaEvent_t start, stop;
        float time;


        int blocks = treePanelsInfo.nBlock;

        //int nTotPan = treePanelsInfo.nObject;
        //int th5slae = nTotPan / blocks;
        //if (th5slae < order)
        //    th5slae = order;
        //else
        //    th5slae = (int)std::ceil((double)nTotPan / blocks);
        //
        //int resultBlock = 1;
        //if (th5slae >= 64)  //Округляем вниз до предыдущей степени двойки, т.е.
        //{                                           // (1024...2047) -> 512
        //    while (resultBlock <= th5slae)           // ( 512...1023) -> 256
        //        resultBlock <<= 1;                  // ( 256...511 ) -> 128
        //    th5slae = min(1024, (resultBlock >> 2)); // ( 128...255 ) ->  64
        //}                                           // (  64...127 ) ->  32        
        //
        //if ((th5slae > 32) && (th5slae < 64))
        //    th5slae = 32;
        //
        //th5slae = 32;

        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        
        if (iter == 0)
        {
            treeMatrToVecZeroIterKernel12<12> << <blocks * FACTORslae, THREADSslae >> > (
                treePanelsInfo.nNode, treePanelsInfo.nObject, itolsq,
                (int2*)treePanelsInfo.childD, (double2*)treePanelsInfo.momsD,
                (double*)treePanelsInfo.gabForLeavesD,
                (double*)treePanelsInfo.objectD,
                treePanelsInfo.mortonCodesIdxD,
                (double2*)treePanelsInfo.centerD,
                treePanelsInfo.indexSortD, treePanelsInfo.indexSortTD, (double4*)treePanelsInfo.lowerupperD,
                treePanelsInfo.schemeType,
                treePanelsInfo.nObject, (double*)treePanelsInfo.gabForLeavesD, (double2*)treePanelsInfo.ED,
                treePanelsInfo.mortonCodesIdxD,
                resD, resLinD,
                treePanelsInfo.matVecMulInfo.closeCellsIdxD, treePanelsInfo.matVecMulInfo.closePrefixSumD,
                treePanelsInfo.matVecMulInfo.farCellsIdxD, treePanelsInfo.matVecMulInfo.farPrefixSumD,
                (double2*)treePanelsInfo.matVecMulInfo.i00D, (double2*)treePanelsInfo.matVecMulInfo.i01D, (double2*)treePanelsInfo.matVecMulInfo.i10D, (double2*)treePanelsInfo.matVecMulInfo.i11D, iter);

            CudaTestError("treeMatrToVecZeroIterKernel12 launch failed");
        }
        else
        {
            treeMatrToVecNonZeroIterKernel12<12> << <blocks * FACTORslae, THREADSslae >> > (
                treePanelsInfo.nNode, treePanelsInfo.nObject, (double2*)treePanelsInfo.momsD,
                (double*)treePanelsInfo.gabForLeavesD,
                (double*)treePanelsInfo.objectD,
                treePanelsInfo.mortonCodesIdxD,
                (double2*)treePanelsInfo.centerD, treePanelsInfo.indexSortTD,
                treePanelsInfo.schemeType,
                treePanelsInfo.nObject, (double*)treePanelsInfo.gabForLeavesD, (double2*)treePanelsInfo.ED,
                treePanelsInfo.mortonCodesIdxD,
                resD, resLinD,
                treePanelsInfo.matVecMulInfo.nClosePanelsD,
                treePanelsInfo.matVecMulInfo.nFarCellsD,
                treePanelsInfo.matVecMulInfo.closeCellsIdxD, treePanelsInfo.matVecMulInfo.closePrefixSumD,
                treePanelsInfo.matVecMulInfo.farCellsIdxD, treePanelsInfo.matVecMulInfo.farPrefixSumD,
                (double2*)treePanelsInfo.matVecMulInfo.i00D, (double2*)treePanelsInfo.matVecMulInfo.i01D, (double2*)treePanelsInfo.matVecMulInfo.i10D, (double2*)treePanelsInfo.matVecMulInfo.i11D);
         
            CudaTestError("treeMatrToVecNonZeroIterKernel12 launch failed");
        }             

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

        cudaEventDestroy(start);  cudaEventDestroy(stop);  

        return time;
    }//treeMatrToVecWrapper(...)
}


