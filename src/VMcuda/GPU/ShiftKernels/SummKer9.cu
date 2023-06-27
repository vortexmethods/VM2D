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
| File name: SummKer_n.cu                                                     |
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
\brief Сдвиг мультипольных моментов для схемы с order = 9
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 29 августа 2023 г.
*/


__global__ 
__launch_bounds__(THREADS3, FACTOR3)
void SummarizationKernel2_9(
    const int nnodesd, const int nbodiesd,
    const int2* __restrict Mchildd,
    volatile int* __restrict massd,
    const int order, real2* __restrict momsd,  //momsd  - без volatile
    const double* __restrict vtxd, int objectType, const int* __restrict MmortonCodesIdxd,
    const real2* __restrict Mposd, const int* __restrict MindexSortd, const int* __restrict MindexSortTd
)
{
    register int i, j, k, ch, inc, flag;

    register real2 mom0;
    register real2 mom1;
    register real2 mom2;
    register real2 mom3;
    register real2 mom4;
    register real2 mom5;
    register real2 mom6;
    register real2 mom7;
    register real2 mom8;

    register real2 cen, dr;

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
                int chd = i * Mchildd[srt].y + (1-i) * Mchildd[srt].x;   // i==0 => .x;  i==1 => .y
                ch = (chd >= nbodiesd) ? chd - nbodiesd : (nnodesd - 1) - MindexSortTd[chd];

                if ((chd >= nbodiesd) || (massd[nnodesd - 1 - ch] >= 0))
                    j--;
            }

            if (j == 0) {
                // all children are ready
                const int kch = ((nnodesd - 1) - k) * order;
                cm = 0;

                const register int sortedCell = MindexSortd[(nnodesd - 1) - k];

                cen = Mposd[sortedCell];
                const int2 chdPair = Mchildd[sortedCell];

                for (i = 0; i < 2; i++)
                {
                    //computation of ch = child[k*2+i]
                    const int chd = i * chdPair.y + (1-i) * chdPair.x;
                    if (chd >= nbodiesd)
                    {
                         ch = chd - nbodiesd;
                         const register int sortedBody = MmortonCodesIdxd[ch];
                         if (objectType == 0)
                         {
                              mom0 = real2{ vtxd[sortedBody*3+2], (real)0 };
                              mom1 = mom2 = mom3 = mom4 = mom5 = mom6 = mom7 = mom8 = real2{ 0, 0 };
                              dr = real2{vtxd[sortedBody*3+0], vtxd[sortedBody*3+1]} - cen;
                              m = 1;
                         } //objectType==0
                         if ((objectType == 1) || (objectType == 2) || (objectType == -1) || (objectType == -2))
                         {
                              real2 panBegin, panEnd;
                              panBegin = real2{vtxd[sortedBody * 12 + 2], vtxd[sortedBody * 12 + 3]};
                              panEnd = real2{vtxd[sortedBody * 12 + 4], vtxd[sortedBody * 12 + 5]};

                              real2 rcur, rd2Pow;
                              rcur = rd2Pow = multz(0.5 * (panEnd - panBegin), 0.5 * (panEnd - panBegin));
                              real gam;
                              switch (objectType)
                              {
                              case 1:
                              case 2:
                                  gam = vtxd[sortedBody * 12 + 6] + vtxd[sortedBody * 12 + 7];
                                  break;

                              case -1:
                              case -2:
                                  gam = vtxd[sortedBody * 12 + 8];
                                  break;
                              };

                              mom1 = mom3 = mom5 = mom7 = real2{ 0, 0 };
                              mom0 = real2{ gam, (real)0 };

                              mom2 = (gam / 3) * rcur; 
                              rcur = multz(rcur, rd2Pow);
                              mom4 = (gam / 5) * rcur; 
                              rcur = multz(rcur, rd2Pow);
                              mom6 = (gam / 7) * rcur; 
                              rcur = multz(rcur, rd2Pow);
                              mom8 = (gam / 9) * rcur; 
                         if ((objectType == 2) || (objectType == -2))
                         {
                              real gamLin;

                              switch (objectType)
                              {
                              case 2:
                                  gamLin = vtxd[sortedBody * 12 + 9] + vtxd[sortedBody * 12 + 10];
                                  break;
                              case -2:
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
                         }

                              dr = real2{vtxd[sortedBody*12+0], vtxd[sortedBody*12+1]} - cen;
                              m = 1;
                         } //objectType==1
                    }
                    else
                    {
                         register const int srtT = MindexSortTd[chd];
                         ch = (nnodesd - 1) - srtT;
                         const int nch = srtT * order;
                         mom0 = real2{ momsd[nch + 0].x, (real)0 };
                         mom1 = momsd[nch + 1];
                         mom2 = momsd[nch + 2];
                         mom3 = momsd[nch + 3];
                         mom4 = momsd[nch + 4];
                         mom5 = momsd[nch + 5];
                         mom6 = momsd[nch + 6];
                         mom7 = momsd[nch + 7];
                         mom8 = momsd[nch + 8];
                         //for (int s = 1; s < order; ++s)
                         //    mom[s] = momsd[ch * order + s];

                         dr = Mposd[chd] - cen;
                         m = massd[nnodesd - 1 - ch];
                     }
                     // add child's contribution
                     momsd[kch + 0].x += mom0.x;

                     register real2 momh1 = mom1;
                     register real2 momh2 = mom2;
                     register real2 momh3 = mom3;
                     register real2 momh4 = mom4;
                     register real2 momh5 = mom5;
                     register real2 momh6 = mom6;
                     register real2 momh7 = mom7;
                     register real2 momh8 = mom8;

                     //for (int p = 1; p < order; ++p)
                     //    momh[p] = mom[p];

                     real2 z = dr;

                     momh1 += multz(mom0, z);
                     momh2 += 2 * multz(mom1, z);
                     momh3 += 3 * multz(mom2, z);
                     momh4 += 4 * multz(mom3, z);
                     momh5 += 5 * multz(mom4, z);
                     momh6 += 6 * multz(mom5, z);
                     momh7 += 7 * multz(mom6, z);
                     momh8 += 8 * multz(mom7, z);

                     z = multz(z, dr);

                     momh2 += multz(mom0, z);
                     momh3 += 3 * multz(mom1, z);
                     momh4 += 6.0 * multz(mom2, z);
                     momh5 += 10.0 * multz(mom3, z);
                     momh6 += 15.0 * multz(mom4, z);
                     momh7 += 21.0 * multz(mom5, z);
                     momh8 += 28.0 * multz(mom6, z);

                     z = multz(z, dr);

                     momh3 += multz(mom0, z);
                     momh4 += 4 * multz(mom1, z);
                     momh5 += 10.0 * multz(mom2, z);
                     momh6 += 20.0 * multz(mom3, z);
                     momh7 += 35.0 * multz(mom4, z);
                     momh8 += 56.0 * multz(mom5, z);

                     z = multz(z, dr);

                     momh4 += multz(mom0, z);
                     momh5 += 5 * multz(mom1, z);
                     momh6 += 15.0 * multz(mom2, z);
                     momh7 += 35.0 * multz(mom3, z);
                     momh8 += 70.0 * multz(mom4, z);

                     z = multz(z, dr);

                     momh5 += multz(mom0, z);
                     momh6 += 6 * multz(mom1, z);
                     momh7 += 21.0 * multz(mom2, z);
                     momh8 += 56.0 * multz(mom3, z);

                     z = multz(z, dr);

                     momh6 += multz(mom0, z);
                     momh7 += 7 * multz(mom1, z);
                     momh8 += 28.0 * multz(mom2, z);

                     z = multz(z, dr);

                     momh7 += multz(mom0, z);
                     momh8 += 8 * multz(mom1, z);

                     z = multz(z, dr);

                     momh8 += multz(mom0, z);

                     //for (int s = 1; s < order; ++s)
                     //{
                     //    for (int p = s; p < order; ++p)
                     //        momh[p] += binomCft[p * order + s] * multz(mom[p - s], z);
                     //    z = multz(z, dr);
                     //}

                     momsd[kch + 1] += momh1;
                     momsd[kch + 2] += momh2;
                     momsd[kch + 3] += momh3;
                     momsd[kch + 4] += momh4;
                     momsd[kch + 5] += momh5;
                     momsd[kch + 6] += momh6;
                     momsd[kch + 7] += momh7;
                     momsd[kch + 8] += momh8;

                     //for (int p = 1; p < order; ++p)
                     //    momsd[k * (order)+p] += momh[p];

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