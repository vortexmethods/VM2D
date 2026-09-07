/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: cpuTreeInfo.cpp                                                  |
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
\brief Реализация класса дерева для реализации быстрых алгоритмов на CPU
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#include "cpuTreeInfo.h"
#include "Gpudefs.h"
#include <algorithm>
#include <omp.h>
#include <stdio.h>
#include "defs.h"
#include "TimesGen.h"

namespace VM2D
{
    //"Разрежение" двоичного представления беззнакового целого, вставляя по одному нулику между всеми битами
    inline unsigned int ExpandBits(unsigned int v)
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

    //Округление "в потолок" результата деления на степень двойки
    inline int ceilpow2(int x, int p) //  =ceil(x / 2^p)
    {
        return (x >> p) + !!(x & ((1 << p) - 1));
    }

    //Округление "в потолок" результата деления пополам
    inline int ceilhalf(int x) //  =ceil(x / 2), т.е. предыдущая функция при p=1
    {
        return (x >> 1) + (x & 1);
    }

    int CpuTreeInfo::Delta(int i, int j) const
    {
        if ((j < 0) || (j > (int)object.size() - 1))
            return -1;

        if (i > j)
            std::swap(i, j);

        //if ((i < 0) || (j > n-1))
        //    exit(111);

        const unsigned int ki = mortonCodesKey[i];
        const unsigned int kj = mortonCodesKey[j];

        //Поиск номера самого старшего ненулевого бита в числе c 
        int count = 0;
        for (unsigned int c = (ki ^ kj); c; c >>= 1, ++count);

        if ((!count) && (i != j))
        {
            int addCount = 0;
            //единички к номерам i и j добавлены для совместимости с Wolfram Mathematica, 
            //для кода они не важны, но дерево без них почти наверняка построится по-другому        
            for (unsigned int add = ((i + 1) ^ (j + 1)); add; add >>= 1, ++addCount);
            return 2 * codeLength + (2 * codeLength - addCount);
        }//if ((!count) && (i != j))

        return (2 * codeLength - count);
    }//Delta(...)




    CpuTreeInfo::CpuTreeInfo(tree_T treeType_, object_T objectType_, scheme_T schemeType_)
        :
        treeType(treeType_),
        objectType(objectType_),
        schemeType(schemeType_)
    {
        /*
        switch (objectType)
        {
        case object_T::point2:
            sizeOfElement = sizeof(Point2D);            
            break;

        case object_T::point4:
            sizeOfElement = sizeof(Vortex2D);
            break;

        case object_T::panel:
            switch (treeType)
            {
            case tree_T::aux:
                sizeOfElement = sizeof(double) * 6;
                break;
            case tree_T::contr:
                sizeOfElement = sizeof(Point2D);
                break;
            case tree_T::vortex:
            case tree_T::source:
                sizeOfElement = sizeof(double) * 12;
                break;
            }
            break;
        }
        */
    };

    CpuTreeInfo::~CpuTreeInfo() {};

    
    float CpuTreeInfo::Update(const std::vector<Vortex2D>& vtx, int cntrLev)
    {
        VMlib::vmTimer timer;
        timer.reset();

        timer.start();
        int nObject = (int)vtx.size();

        if (cntrLev == 0)
            controlLevel = std::max(4, (int)(log2(nObject) - 3));
        else 
            controlLevel = cntrLev;

        indexControlCells.resize(0);
        indexControlCells.reserve(nObject);

        bool inflTree = (treeType == tree_T::vortex || treeType == tree_T::source);

        object.resize(nObject);
        //if (objectType == object_T::point4 && inflTree) 
        {
            gamma.resize(nObject);
            sigma.resize(nObject);
        }

        gabForLeaves.resize(nObject);
        mortonCodesKeyUnsort.resize(nObject);
        mortonCodesIdxUnsort.resize(nObject);
        mortonCodesKey.resize(nObject);
        mortonCodesIdx.resize(nObject);

        levelUnsort.resize(nObject - 1);
        levelSort.resize(nObject - 1);
        indexUnsort.resize(nObject - 1);
        indexSort.resize(nObject - 1);
        indexSortT.resize(nObject - 1);
        range.resize(2 * nObject);
        parent.resize(2 * nObject);
        child.resize(nObject - 1);
        lowerupper.resize(nObject - 1);
        center.resize(nObject - 1);

        if (inflTree)
        {
            moms.resize((nObject - 1) * orderAlignment);
            mass.resize(nObject - 1);
        }


        for (int v = 0; v < nObject; ++v)
        {
            Point2D r = vtx[v].r();
            object[v] = r;
            if (objectType == object_T::point4 /* && inflTree*/) //влияющее дерево вихрей
            {
                //if (inflTree)
                {
                    gamma[v] = vtx[v].g();
                    sigma[v] = vtx[v].sigma();
                }
                gabForLeaves[v] = Point4D({ r[0] - sigma[v], r[1] - sigma[v], r[0] + sigma[v], r[1] + sigma[v] });
            }  
            else
                gabForLeaves[v] = Point4D({ r[0] - 0.0, r[1] - 0.0, r[0] + 0.0, r[1] + 0.0 });            
        }

        timer.stop();
        return (float)timer.duration();
    }



    float CpuTreeInfo::UpdatePanelGeometry(const std::vector<std::pair<Point2D, Point2D>>& panels, int cntrLev)
    {
        VMlib::vmTimer timer;
        timer.reset();

        timer.start();

        int nObject = (int)panels.size();

        if (cntrLev == 0)
            controlLevel = std::max(4, (int)(log2(nObject) - 3));
        else
            controlLevel = cntrLev;

        indexControlCells.resize(0);
        indexControlCells.reserve(nObject);

        object.resize(nObject);
        
        if (treeType != tree_T::aux && treeType != tree_T::contr)
        {
            gamma.resize(nObject);
            sigma.resize(nObject);
        }

        gabForLeaves.resize(nObject);
        mortonCodesKeyUnsort.resize(nObject);
        mortonCodesIdxUnsort.resize(nObject);
        mortonCodesKey.resize(nObject);
        mortonCodesIdx.resize(nObject);

        levelUnsort.resize(nObject - 1);
        levelSort.resize(nObject - 1);
        indexUnsort.resize(nObject - 1);
        indexSort.resize(nObject - 1);
        indexSortT.resize(nObject - 1);
        range.resize(2 * nObject);
        parent.resize(2 * nObject);
        child.resize(nObject - 1);
        lowerupper.resize(nObject - 1);
        center.resize(nObject - 1);

        mass.resize(nObject - 1);

        if (treeType != tree_T::contr)
        {
            moms.resize((nObject - 1) * orderAlignment);        
        }

        for (int v = 0; v < nObject; ++v)
        {
            Point2D r = 0.5 * (panels[v].first + panels[v].second);
            object[v] = r;

            if (treeType != tree_T::aux && treeType != tree_T::contr)
            {
            //    gamma[v] = vtx[v].g();
            //    sigma[v] = vtx[v].sigma();
            }



            gabForLeaves[v] = Point4D({ panels[v].first[0], panels[v].first[1], panels[v].second[0], panels[v].second[1] });
        }

        timer.stop();
        return (float)timer.duration();
    }




    //Сортировка листьев
    void CpuTreeInfo::RadixSortMortonCodes()
    {
        int nObject = (int)object.size();
        if (nObject <= 0)
            return;

        std::memcpy(mortonCodesKey.data(), mortonCodesKeyUnsort.data(), mortonCodesKeyUnsort.size() * sizeof(unsigned));
        std::memcpy(mortonCodesIdx.data(), mortonCodesIdxUnsort.data(), mortonCodesIdxUnsort.size() * sizeof(int));

        codesSorter.sort((int*)mortonCodesKey.data(), mortonCodesIdx.data(), 0, mortonCodesKeyUnsort.size());

        /*
        std::vector<std::pair<unsigned, int>> pr(mortonCodesKeyUnsort.size());
        for (size_t q = 0; q < mortonCodesKeyUnsort.size(); ++q)
            pr[q] = { mortonCodesKeyUnsort[q], mortonCodesIdxUnsort[q] };
        std::sort(pr.begin(), pr.end(), [](auto ka, auto kb) {return ka.first < kb.first;});

        for (size_t q = 0; q < mortonCodesKeyUnsort.size(); ++q)
        {
            mortonCodesKey[q] = pr[q].first;
            mortonCodesIdx[q] = pr[q].second;
        }
        */

    }//RadixSortMortonCodes()

    void CpuTreeInfo::RadixSortInternalCells()
    {
        int nObject = (int)object.size();
        int n = nObject - 1;

        if (n == 0)
            return;


        std::memcpy(levelSort.data(), levelUnsort.data(), levelUnsort.size() * sizeof(unsigned));
        std::memcpy(indexSort.data(), indexUnsort.data(), levelUnsort.size() * sizeof(int));

        levelSorter.sort((int*)levelSort.data(), indexSort.data(), 0, levelUnsort.size());

        /*
        std::vector<std::pair<unsigned, int>> pr(levelUnsort.size());
        for (size_t q = 0; q < levelUnsort.size(); ++q)
            pr[q] = { levelUnsort[q], indexUnsort[q] };
        std::sort(pr.begin(), pr.end(), [](auto ka, auto kb) {return ka.first < kb.first;});

        for (size_t q = 0; q < levelUnsort.size(); ++q)
        {
            levelSort[q] = pr[q].first;
            indexSort[q] = pr[q].second;
        }
        */

        for (int k = 0; k < n; ++k)
            indexSortT[indexSort[k]] = k;

    }//RadixSortInternalCells()


   

    void CpuTreeInfo::CalcAABB()
    {
        using double4 = Point4D;
        using double2 = Point2D;
        using int2 = std::pair<int, int>;

#pragma omp parallel
        {
            int i, j, ch, flag;

            double4 lu[2];

            //int cm;
            int m[2];


            const int nnodes = 2 * (int)object.size() - 1;
            const int nbodies = (int)object.size();
#pragma omp for
            for (int k = nbodies; k < nnodes; ++k)
            {
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
                while (flag == 0)
                {
                    j = 2;
                    const int kch = ((nnodes - 1) - k);
                    const int srt = indexSort[kch];
                    int2 chdPair = child[srt];

                    //cm = 0;

                    //int chdSorted[2];

                    for (i = 0; i < 2; i++)
                    {
                        int chd = i * chdPair.second + (1 - i) * chdPair.first;   // i==0 => .x;  i==1 => .y
                        ch = (chd >= nbodies) ? (chd - nbodies) : ((nnodes - 1) - indexSortT[chd]);
                        if ((chd >= nbodies) || (mass[nnodes - 1 - ch] >= 0))
                            j--;
                    }

                    if (j == 0)
                    {
                        for (i = 0; i < 2; i++)
                        {
                            const int chd = i * chdPair.second + (1 - i) * chdPair.first;
                            if (chd >= nbodies)
                            {
                                ch = chd - nbodies;
                                const int sortedBody = mortonCodesIdx[ch];

                                double4 xyAB = gabForLeaves[sortedBody];
                                lu[i] = double4{
                                    ::fmin(xyAB[0], xyAB[2]), ::fmin(xyAB[1], xyAB[3]),
                                    ::fmax(xyAB[0], xyAB[2]), ::fmax(xyAB[1], xyAB[3])
                                };
                                m[i] = 1;
                            }
                            else
                            {
                                const int srtT = indexSortT[chd];
                                lu[i] = lowerupper[chd];
                                m[i] = mass[srtT];
                            }
                        }

                        const double4 loup = double4{
                                                ::fmin(lu[0][0], lu[1][0]),
                                                ::fmin(lu[0][1], lu[1][1]),
                                                ::fmax(lu[0][2], lu[1][2]),
                                                ::fmax(lu[0][3], lu[1][3]) };

                        lowerupper[srt] = loup;

                        // Центр текущего узла = центр его AABB
                        center[srt] = double2{
                            0.5 * (loup[0] + loup[2]),
                            0.5 * (loup[1] + loup[3])
                        };
                        flag = 1;
                    }//if j==0

#pragma omp flush

                    if (flag != 0)
                        mass[nnodes - 1 - k] = m[0] + m[1];

                }//while flag
            }//for k
        }
    }//CalcAABB(...)



    float CpuTreeInfo::Build()
    {
        VMlib::vmTimer timer;
        timer.start();

        int nObject = (int)object.size();
        if (nObject > 0)
        {
            //treeBoundingBox;
            //to_parallel
            auto minmaxX = std::minmax_element(object.begin(), object.end(), Point2D::cmp<'x'>);
            auto minmaxY = std::minmax_element(object.begin(), object.end(), Point2D::cmp<'y'>);
            minr = Point2D({ (*minmaxX.first)[0], (*minmaxY.first)[1] });
            maxr = Point2D({ (*minmaxX.second)[0], (*minmaxY.second)[1] });

            //treeMortonCodes;
            double lmax, quadSideFactor;
            lmax = std::max(maxr[0] - minr[0], maxr[1] - minr[1]);
            Point2D rcen = 0.5 * (maxr + minr); //координаты центра
            quadSideFactor = rbound / lmax; //1;

#pragma omp parallel for
            for (int bdy = 0; bdy < nObject; ++bdy)
            {
                Point2D rScaled = twoPowCodeLength * ((object[bdy] - rcen) * quadSideFactor + 0.5 * Point2D{ rbound, rbound });

                unsigned int xx = ExpandBits((unsigned int)rScaled[0]);
                unsigned int yy = ExpandBits((unsigned int)rScaled[1]);
                mortonCodesKeyUnsort[bdy] = yy | (xx << 1);
                mortonCodesIdxUnsort[bdy] = bdy;
                //printf("{%d, %d},\n", bdy, mortonCodesKeyUnsort[bdy]);
            }



            RadixSortMortonCodes();
            //for (int i = 0; i < nObject; ++i)
                //printf("{%d, %d},\n",  mortonCodesKey[i], mortonCodesIdx[i]);

            //treeMortonInternalNodes
            //if (treeType != tree_T::contr) для CPU нужны внутренние ячейки для обоих деревьев
            {
#pragma omp parallel for
                for (int i = 0; i < nObject - 1; ++i)
                {
                    int codei = mortonCodesKey[i];

                    int Deltap1 = Delta(i, i + 1);
                    int Deltam1 = Delta(i, i - 1);

                    int d = (Deltap1 - Deltam1 > 0) ? 1 : -1;

                    int delta_min = (d > 0) ? Deltam1 : Deltap1;

                    int Lmax = 2;
                    int pos = i + Lmax * d;

                    while (Delta(i, pos) > delta_min)
                    {
                        Lmax *= 2;
                        pos = i + Lmax * d;
                    }

                    int L = 0;
                    for (int t = (Lmax >> 1); t >= 1; t >>= 1)
                    {
                        pos = i + (L + t) * d;

                        if (Delta(i, pos) > delta_min)
                            L += t;
                    }

                    int j = i + L * d;
                    pos = j;

                    int delta_node = Delta(i, j);

                    levelUnsort[i] = delta_node;
                    indexUnsort[i] = i;
                    //if (i == 2)
                    //    printf("!!! %d, %d\n", mortonCodesKey[i], mortonCodesKey[j]);

                    int s = 0;
                    for (int p = 1, t = ceilhalf(L); L > (1 << (p - 1)); ++p, t = ceilpow2(L, p))
                    {
                        pos = i + (s + t) * d;

                        int dl = Delta(i, pos);

                        if (dl > delta_node)
                            s += t;
                    }//for p

                    int gammaPos = i + s * d + d * (d < 0);   //                    = std::min(d, 0);

                    int Mmin = std::min(i, j);
                    int Mmax = std::max(i, j);

                    int left = gammaPos;
                    int right = gammaPos + 1;

                    //               -                         
                    int childLeft = (Mmin == gammaPos) * nObject + left;
                    range[childLeft] = { Mmin, gammaPos };
                    parent[childLeft] = i;

                    //                -                         
                    int childRight = (Mmax == gammaPos + 1) * nObject + right;
                    range[childRight] = { gammaPos + 1, Mmax };
                    parent[childRight] = i;

                    child[i] = { childLeft, childRight };
                    //printf("i:%d, %d, {%d, %d}\n", indexUnsort[i], levelUnsort[i], childLeft, childRight);
                }

                //if (treeType == tree_T::contr)
                {
                    for (int i = 0; i < nObject - 1; ++i)
                    {
                        if (levelUnsort[i] >= controlLevel && levelUnsort[parent[i]] < controlLevel)
                            indexControlCells.push_back(i);
                    }
                    for (int i = 0; i < nObject; ++i)
                    {
                        if (levelUnsort[parent[nObject + i]] < controlLevel)
                            indexControlCells.push_back(nObject + i);
                    }
                }
                RadixSortInternalCells();
            }// if contr
        }
        timer.stop();
        
        return (float)timer.duration();
    }//Build()


    float CpuTreeInfo::UpwardTraversal(int order)
    {
        VMlib::vmTimer timer;
        timer.start();

        if (object.size() > 0)
        {
            mass.assign(object.size() - 1, -1);

            if (treeType != tree_T::contr && treeType != tree_T::aux)
            {
#include <SummChoice.h>
            }
            else
                CalcAABB();
        }
        timer.stop();
        
        return (float)timer.duration();
    }//UpwardTraversal(...)


   
    inline double mindist2(const Point4D& lowerUpper, const Point2D& rhs) noexcept
    {        
        const double dx = ::fmin(lowerUpper[2], ::fmax(lowerUpper[0], rhs[0])) - rhs[0];
        const double dy = ::fmin(lowerUpper[3], ::fmax(lowerUpper[1], rhs[1])) - rhs[1];
        return dx * dx + dy * dy;
    }

    inline double minmaxdist2(const Point4D& lowerUpper, const Point2D& rhs) noexcept
    {
       std::pair<double, double> rm_sq = std::make_pair<double, double>((lowerUpper[0] - rhs[0]) * (lowerUpper[0] - rhs[0]),
            (lowerUpper[1] - rhs[1]) * (lowerUpper[1] - rhs[1]));
       std::pair<double, double> rM_sq = std::make_pair<double, double>((lowerUpper[2] - rhs[0]) * (lowerUpper[2] - rhs[0]),
            (lowerUpper[3] - rhs[1]) * (lowerUpper[3] - rhs[1]));

        if ((lowerUpper[2] + lowerUpper[0]) * 0.5 < rhs[0])
            std::swap(rm_sq.first, rM_sq.first);

        if ((lowerUpper[3] + lowerUpper[1]) * 0.5 < rhs[1])
            std::swap(rm_sq.second, rM_sq.second);

        const double dx = rm_sq.first + rM_sq.second;
        const double dy = rM_sq.first + rm_sq.second;

        return fmin(dx, dy);
    }

 
    inline std::pair<double, int> distance_calculator_point2segment(const Point2D& point, const Point4D& object)
    {
        int location;

        double a = object[3] - object[1];
        double b = object[2] - object[0];

        double distanceSegment;

        std::pair<double, double> dr = std::make_pair<double, double>(point[0] - object[0], point[1] - object[1]);

        double r_numerator = dr.first * b + dr.second * a;
        double r_denomenator = b * b + a * a;
        double r = r_numerator / r_denomenator;

        double s = (dr.first * a - dr.second * b);

        if ((r >= 0) && (r <= 1))
        {
            distanceSegment = s * s / r_denomenator;
            location = 0;
        }
        else
        {
            double dist1 = dr.first * dr.first + dr.second * dr.second;
            double dist2 = sqr(point[0] - object[2]) + sqr(point[1] - object[3]);
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

        return std::make_pair(distanceSegment, location);
    }

  inline double MinDist2AABB(const Point4D& gabI, const Point4D& gabJ)
  {
      double dx = 0.0, dy = 0.0;
      
      if (gabI[0] > gabJ[2])
          dx = gabI[0] - gabJ[2];
      else if (gabJ[0] > gabI[2])
          dx = gabJ[0] - gabI[2];

      if (gabI[1] > gabJ[3])
          dy = gabI[1] - gabJ[3];
      else if (gabJ[1] > gabI[3])
          dy = gabJ[1] - gabI[3];

      return dx * dx + dy * dy;
  }


  const double onePlusMachineEps = 1.0f + 2.0e-7f;

  float CpuTreeInfo::DownwardTraversalClosestPanelToPoints(CpuTreeInfo& cntrTree, std::vector<std::pair<int, double>>& indexPnlDist, bool findOnlyInside, double* pseudonormals)
  {      
      using double4 = Point4D;
      const int nbodies = (int)object.size(); //количество вихрей
      const int nnodes = 2 * nbodies - 1; //количество узлов дерева вихрей (листья-вихри + внутренние узлы)
      const int npoints = (int)cntrTree.object.size(); //количество точек наблюдения

      int indexOfPoint;
      const int maxDepth = 32;
      
      VMlib::vmTimer timer;
      timer.start();
      
#pragma omp parallel for schedule(dynamic, 1)
      for (int indexK = 0; indexK < (int)cntrTree.indexControlCells.size(); ++indexK)
      {
          Point2D query;

          const int k = cntrTree.indexControlCells[indexK];
          
          Point4D gabCntrl;
          bool isCntrLeaf = (k >= cntrTree.object.size());

          if (isCntrLeaf)
          {
              const int leafIndex = cntrTree.mortonCodesIdx[k - npoints];
              indexOfPoint = leafIndex;
              query = cntrTree.object[leafIndex];
              Point4D pnlCntrGab = cntrTree.gabForLeaves[cntrTree.mortonCodesIdx[leafIndex]];
              gabCntrl = {
                  fmin(pnlCntrGab[0], pnlCntrGab[2]), fmin(pnlCntrGab[1], pnlCntrGab[3]),
                  fmax(pnlCntrGab[0], pnlCntrGab[2]), fmax(pnlCntrGab[1], pnlCntrGab[3])
              };
          }
          else
          {
              query = cntrTree.center[k];
              gabCntrl = cntrTree.lowerupper[k];
          }

          std::pair<int, double> stack[maxDepth];
          int depth;

          double md2 = mindist2(lowerupper[0], query);
          std::pair<int, int> nearestLocation{ -1, 0 };

          if ((findOnlyInside) && (md2 != 0))
          {
              indexPnlDist[indexK] = { -1, -1.0 };
              continue;
          }

          depth = 0;
          stack[0] = { 0, md2 };

          double dist_to_nearest_object = 1e+100;

          while (depth >= 0)
          {
              std::pair<int, double> nd = stack[depth];
              --depth;

              if (nd.second > dist_to_nearest_object)
                  continue;

              std::pair<int, int> chBoth = child[nd.first];

              bool isLeftLeaf = (chBoth.first >= object.size());
              bool isRightLeaf = (chBoth.second >= object.size());

              double4 L_box, R_box;
              double4 pnlLeftGab, pnlRightGab;

              if (!isLeftLeaf)
                  L_box = lowerupper[chBoth.first];
              else
              {
                  int n = chBoth.first - (int)object.size();
                  pnlLeftGab = gabForLeaves[mortonCodesIdx[n]];
                  L_box = {
                      fmin(pnlLeftGab[0], pnlLeftGab[2]), fmin(pnlLeftGab[1], pnlLeftGab[3]),
                      fmax(pnlLeftGab[0], pnlLeftGab[2]), fmax(pnlLeftGab[1], pnlLeftGab[3])
                  };
              }

              if (!isRightLeaf)
                  R_box = lowerupper[chBoth.second];
              else
              {
                  int n = chBoth.second - (int)object.size();
                  pnlRightGab = gabForLeaves[mortonCodesIdx[n]];
                  R_box = {
                      fmin(pnlRightGab[0], pnlRightGab[2]), fmin(pnlRightGab[1], pnlRightGab[3]),
                      fmax(pnlRightGab[0], pnlRightGab[2]), fmax(pnlRightGab[1], pnlRightGab[3])
                  };
              }

              double L_mindist2 = mindist2(L_box, query);
              double R_mindist2 = mindist2(R_box, query);

              double L_minmaxdist2 = minmaxdist2(L_box, query);
              double R_minmaxdist2 = minmaxdist2(R_box, query);

              if (L_mindist2 <= R_minmaxdist2 * onePlusMachineEps) // L is worth considering
              {
                  if (isLeftLeaf) // leaf node
                  {
                      int n = chBoth.first - (int)object.size();

                      std::pair<double, int> dist_code = distance_calculator_point2segment(query, pnlLeftGab);
                      if (dist_code.first <= dist_to_nearest_object)
                      {
                          dist_to_nearest_object = dist_code.first;
                          nearestLocation = { mortonCodesIdx[n], dist_code.second };
                      }
                  }
                  else
                  {
                      ++depth;
                      stack[depth] = { chBoth.first, L_mindist2 };
                  }
              }

              if (R_mindist2 <= L_minmaxdist2 * onePlusMachineEps) // R is worth considering
              {
                  if (isRightLeaf) // leaf node
                  {
                      int n = chBoth.second - (int)object.size();

                      std::pair<double, int> dist_code = distance_calculator_point2segment(query, pnlRightGab);
                      if (dist_code.first <= dist_to_nearest_object)
                      {
                          dist_to_nearest_object = dist_code.first;
                          nearestLocation = { mortonCodesIdx[n], (int)dist_code.second };
                      }
                  }
                  else
                  {
                      ++depth;
                      stack[depth] = { chBoth.second, R_mindist2 };
                  }
              }
          }//while depth


          if (isCntrLeaf)
          {
              indexPnlDist[cntrTree.mortonCodesIdx[k - npoints]] = { nearestLocation.first, sqrt(dist_to_nearest_object) };
              continue;
          }

          double radius = sqrt(dist_to_nearest_object) + 0.5 * sqrt(sqr(gabCntrl[2] - gabCntrl[0]) + sqr(gabCntrl[3] - gabCntrl[1]));
          

          int stackCandidates[maxDepth];
          depth = 0;          
          stackCandidates[0] = 0;

          std::vector<int> candidates;
          candidates.reserve(100);

          while (depth >= 0)
          {
              int nd = stackCandidates[depth];
              --depth;

              std::pair<int, int> chBoth = child[nd];

              bool isLeftLeaf = (chBoth.first >= object.size());
              bool isRightLeaf = (chBoth.second >= object.size());

              Point4D gabL, gabR;

              if (isLeftLeaf)
              {
                  int n = chBoth.first - (int)object.size();
                  Point4D pnlLeftGab = gabForLeaves[mortonCodesIdx[n]];
                  gabL = {
                      fmin(pnlLeftGab[0], pnlLeftGab[2]), fmin(pnlLeftGab[1], pnlLeftGab[3]),
                      fmax(pnlLeftGab[0], pnlLeftGab[2]), fmax(pnlLeftGab[1], pnlLeftGab[3])
                  };
              }
              else
                  gabL = lowerupper[chBoth.first];
                        

              if (isRightLeaf)          
              {
                  int n = chBoth.second - (int)object.size();
                  Point4D pnlRightGab = gabForLeaves[mortonCodesIdx[n]];
                  gabR = {
                      fmin(pnlRightGab[0], pnlRightGab[2]), fmin(pnlRightGab[1], pnlRightGab[3]),
                      fmax(pnlRightGab[0], pnlRightGab[2]), fmax(pnlRightGab[1], pnlRightGab[3])
                  };
              }
              else
                  gabR = lowerupper[chBoth.second];
              

              double dist2L = MinDist2AABB(gabCntrl, gabL);
              double dist2R = MinDist2AABB(gabCntrl, gabR);
                        
              if (dist2L < sqr(radius))
              {              
                  if (isLeftLeaf)
                      candidates.push_back(mortonCodesIdx[chBoth.first - object.size()]);
                  else                  
                      stackCandidates[++depth] = chBoth.first;                  
              }

              if (dist2R < sqr(radius))
              {                  
                  if (isRightLeaf)
                      candidates.push_back(mortonCodesIdx[chBoth.second - object.size()]);
                  else                      
                      stackCandidates[++depth] = chBoth.second;                  
              }

          }//while depth

          for (int particle = cntrTree.range[k].first; particle <= cntrTree.range[k].second; ++particle)
          {          
              Point2D observPoint = cntrTree.object[cntrTree.mortonCodesIdx[particle]];

              double currentMinDist2 = 1e+100;
              int indexMinDist = -1;
              for (int i = 0; i < candidates.size(); ++i)
              {
                  Point4D pnl = gabForLeaves[candidates[i]];                 
                  double dst2 = distance_calculator_point2segment(observPoint, gabForLeaves[candidates[i]]).first;
                  if (dst2 < currentMinDist2)
                  {
                      currentMinDist2 = dst2;
                      indexMinDist = candidates[i];
                  }
              }
              indexPnlDist[cntrTree.mortonCodesIdx[particle]] = { indexMinDist, sqrt(currentMinDist2) };
          }
      }//for k

      //*/
            
      timer.stop();      
      return (float)timer.duration();
  }

      
    float CpuTreeInfo::DownwardTraversalVorticesToPoints(CpuTreeInfo& cntrTree, std::vector<Point2D>& vel, std::vector<double>& epsast, double theta, int order, bool calcRadius)
    {
        float t1 = (float)omp_get_wtime();

        using double4 = Point4D;
        using double2 = Point2D;
        using int2 = std::pair<int, int>;

        if (object.size() > 0)
        {
#pragma omp parallel
            {
                double itolsq = 1.0 / (theta * theta);

                const int nbodies = (int)object.size(); //количество вихрей
                const int nnodes = 2 * nbodies - 1; //количество узлов дерева вихрей (листья-вихри + внутренние узлы)
                const int npoints = (int)cntrTree.object.size(); //количество точек наблюдения

                int nd;      // индекс родительской ячейки дерева, которую обходим, и которая находится на верхушке стека; 
                int2 chBoth; //индексы обоих потомков узла nd
                int pd;      // 0 или 1 --- какого потомка ячейки nd обходим (левого или правого)   

                // для вычисления квадратов расстояний до трех ближайших вихрей
                double d_1, d_2, d_3, dst23, dst12;
                int indexOfPoint; //истинный индекс точки наблюдения	
                int srtT;         //истинный индекс внутреннего узла дерева

                double2 p;        //координаты точки наблюдения
                double2 v{ 0.0,0.0 };        //результат расчета скорости в точке наблюдения

                bool isVortex;             //признак того, что обрабатываемая вершина - лист
                double2 ps;       //координаты влияющего вихря, если обрабатываемая вершина --- лист, или центра влияющей ячейки, если обрабатывается внутренняя ячейка дерева
                double gm;        //циркуляция вихря, если обрабатываемая вершина --- лист, или 0-й мультипольный момент (суммарная циркуляция вихрей) если обрабатывается внутренняя ячейка дерева
                double sgm2;
                double sgm2Contr;
                double sumSide2;  //габариты обрабатываемой ячейки
                double2 dr;       //радиус-вектор из центра влияющей ячейки в точку наблюдения    
                double r2;        //квадрат модуля предыдущего

                const int maxDepth = 32;

                int posStack[maxDepth];
                int nodeStack[maxDepth];
                int depth;
                const Point2D* mom = nullptr;

#pragma omp for schedule(dynamic,1) 
                for (int indexK = 0; indexK < (int)cntrTree.indexControlCells.size(); ++indexK)
                {
                    const int k = cntrTree.indexControlCells[indexK];
                    std::vector<Point2D> vParticles;
                    std::vector<numvector<double, 3>> epsastParticles;

                    v.toZero();
                    d_1 = d_2 = d_3 = 1e+5;
                    std::vector<Point2D> Ek(order, { 0.0, 0.0 });

                    bool isCntrLeaf = (k >= npoints);

                    if (isCntrLeaf)
                    {
                        const int leafIndex = cntrTree.mortonCodesIdx[k - npoints];
                        indexOfPoint = leafIndex;
                        p = cntrTree.object[leafIndex];
                        sgm2Contr = sqr(sigma[leafIndex]);
                    }
                    else
                    {
                        p = cntrTree.center[k];
                        const Point4D& gab = cntrTree.lowerupper[k];
                        sgm2Contr = sqr((gab[2] - gab[0] + gab[3] - gab[1]));
                        vParticles.resize(cntrTree.range[k].second - cntrTree.range[k].first + 1, { 0.0, 0.0 });
                        epsastParticles.resize(cntrTree.range[k].second - cntrTree.range[k].first + 1, { 1e+5, 1e+5, 1e+5 });
                    }

                    depth = 0;
                    posStack[0] = 0;
                    nodeStack[0] = nnodes - 1;

                    while (depth >= 0)
                    {
                        pd = posStack[depth];
                        nd = nodeStack[depth];

                        chBoth = child[indexSort[(nnodes - 1) - nd]];

                        while (pd < 2)
                        {
                            const int chd = (pd == 0) ? chBoth.first : chBoth.second;

                            ++pd;
                            posStack[depth] = pd;

                            isVortex = (chd >= nbodies);

                            int n;
                            if (isVortex)
                            {
                                n = chd - nbodies;
                                const int vortexIndex = mortonCodesIdx[n]; // истинный индекс вихря

                                ps = object[vortexIndex];
                                gm = gamma[vortexIndex];
                                sgm2 = sqr(sigma[vortexIndex]);
                                sumSide2 = 0.0;                 // ???
                            }
                            else
                            {
                                srtT = indexSortT[chd];
                                n = (nnodes - 1) - srtT;

                                ps = center[chd];

                                const Point4D& gab = lowerupper[chd];
                                sumSide2 = sqr((gab[2] - gab[0] + gab[3] - gab[1]));
                            }

                            dr = p - ps;
                            r2 = dr[0] * dr[0] + dr[1] * dr[1];
                            if ((isVortex && isCntrLeaf) || ((sumSide2 + sgm2Contr) * itolsq < r2)) // если выполнен критерий дальности (может быть как влияние кластера, так и одного дальнего вихря)
                            {
                                // Если это ячейка, берём её mm
                                if (!isVortex)
                                {
                                    mom = moms.data() + srtT * orderAlignment;
                                    gm = mom[0][0]; // нулевой момент = суммарная циркуляция
                                    sgm2 = 0.0;
                                }

                                if ((calcRadius) && (isCntrLeaf))
                                {
                                    if ((r2 < d_3) && (r2 > 0.0))
                                    {
                                        dst23 = std::fmin(r2, d_2);
                                        d_3 = std::fmax(r2, d_2);

                                        dst12 = std::fmin(dst23, d_1);
                                        d_2 = std::fmax(dst23, d_1);

                                        d_1 = dst12;
                                    }
                                }
                                const double f = gm / std::fmax(r2, sgm2);

                                if (isCntrLeaf)
                                    v += f * dr;

                                if (!isCntrLeaf)
                                    Ek[0] += f * dr;  // = theta[m] * m[m]

                                // Вклад высших мм если влияет ячейка
                                if ((order > 1) && r2 > 0.0)
                                {
                                    Point2D cftr = (1.0 / r2) * dr;
                                    Point2D th = cftr;

                                    if (isVortex)
                                    {
                                        if (!isCntrLeaf)
                                            for (int s = 1; s < order; ++s)
                                            {
                                                th = s * multz(th, cftr);
                                                Ek[s] += (s % 2 ? -1.0 : 1.0) * th * gm;
                                            }
                                    }
                                    else
                                    {
                                        for (int s = 1; s < order; ++s)
                                        {
                                            th = s * multz(th, cftr);

                                            if (isCntrLeaf)
                                            {
                                                Point2D add = ifac[s + 1] * multzA(th, mom[s]);
                                                v += add;
                                            }
                                            else
                                            {
                                                for (int q = 0; q <= s; ++q)
                                                {
                                                    //if (s == q)
                                                    //    std::cout << "E[" << q << "] += " << (q % 2 ? -1.0 : 1.0) << " * " << ifac[s - q + 1] << " * " << th << " * " << mom[s - q] << "\n";
                                                    Ek[q] += ((q % 2 ? -1.0 : 1.0) * ifac[s - q + 1]) * multzA(th, mom[s - q]);
                                                }

                                                if (s == order - 1)
                                                    Ek[0] += ifac[s + 1] * multzA(multz(th, cftr), mom[s]);
                                            }
                                        }
                                    }
                                }
                            }
                            else if (!isVortex)
                            {
                                // Ячейка слишком близка -> спускаемся ниже по дереву
                                if (depth + 1 < maxDepth)
                                {
                                    if (pd == 1)//если это была обработка левого потомка, и по нему пришлось идти вниз по дереву
                                    {
                                        posStack[depth] = 1;
                                        nodeStack[depth] = nd;
                                        ++depth;
                                    }

                                    nd = n;
                                    pd = 0;

                                    chBoth = child[indexSort[(nnodes - 1) - nd]];
                                }
                            }
                            else if (!isCntrLeaf)
                            {
                                for (int pnt = cntrTree.range[k].first; pnt <= cntrTree.range[k].second; ++pnt)
                                {
                                    Point2D dz = cntrTree.object[cntrTree.mortonCodesIdx[pnt]] - ps;
                                    double dz2 = dz.length2();

                                    vParticles[pnt - cntrTree.range[k].first] += (gm / std::fmax(dz2, sgm2)) * dz;

                                    if (calcRadius)
                                    {
                                        double& m_3 = epsastParticles[pnt - cntrTree.range[k].first][2];
                                        double& m_2 = epsastParticles[pnt - cntrTree.range[k].first][1];
                                        double& m_1 = epsastParticles[pnt - cntrTree.range[k].first][0];

                                        if ((dz2 < m_3) && (dz2 > 0.0))
                                        {
                                            dst23 = std::fmin(dz2, m_2);
                                            m_3 = std::fmax(dz2, m_2);

                                            dst12 = std::fmin(dst23, m_1);
                                            m_2 = std::fmax(dst23, m_1);

                                            m_1 = dst12;
                                        }
                                    }
                                }
                            } //  БС                                                
                        }//while pd

                         // Оба потомка обработаны -> снимаем узел со стека
                        --depth;
                    }//while depth


                    if (isCntrLeaf)
                    {
                        vel[indexOfPoint] = IDPI * v.kcross();
                        if (calcRadius)
                            epsast[indexOfPoint] = 1.0 * sqrt((d_1 + d_2 + d_3) / 3.0);
                    }
                    else
                    {
                        for (int pnt = cntrTree.range[k].first; pnt <= cntrTree.range[k].second; ++pnt)
                        {
                            Point2D dz = cntrTree.object[cntrTree.mortonCodesIdx[pnt]] - p;
                            Point2D dzPow = dz;
                            v = Ek[0] + vParticles[pnt - cntrTree.range[k].first];
                            for (int q = 1; q < order; ++q)
                            {
                                v += ifac[q + 1] * multzA(Ek[q], dzPow);
                                dzPow = multz(dzPow, dz);
                            }
                            vel[cntrTree.mortonCodesIdx[pnt]] = IDPI * v.kcross();

                            if (calcRadius)
                            {
                                const auto& es = epsastParticles[pnt - cntrTree.range[k].first];
                                epsast[cntrTree.mortonCodesIdx[pnt]] = 1.0 * sqrt((es[0] + es[1] + es[2]) / 3.0);
                            }
                        }
                    }
                }//for k
            }
        }//pragma omp parallel

        float t2 = (float)omp_get_wtime();
        return (t2 - t1);
    }//DownwardTraversalVorticesToPoints(...)


  float CpuTreeInfo::DownwardTraversalVorticesToPanels(CpuTreeInfo& cntrTree, std::vector<double>& rhs, std::vector<double>& rhsLin, double theta, int order)
  {
      float t1 = (float)omp_get_wtime();

      using double4 = Point4D;
      using double2 = Point2D;
      using int2 = std::pair<int, int>;

      const int nbodies = (int)object.size(); //количество вихрей
      const int npoints = (int)cntrTree.object.size(); //количество панелей

      rhs.assign(npoints, 0.0);
      bool scheme = false;
      if (cntrTree.schemeType == scheme_T::linScheme)
          scheme = true;

      if(scheme)
          rhsLin.assign(npoints, 0.0); 
      else
          rhsLin.clear();

      if (nbodies <= 0 || npoints <= 0)
          return 0.0f;

//#pragma omp parallel
      {
          double itolsq = 1.0 / (theta * theta);
          const int nnodes = 2 * nbodies - 1; //количество узлов дерева вихрей (листья-вихри + внутренние узлы) 

          const int maxDepth = 32;

#pragma omp parallel for schedule(dynamic, 1)//Временно: обход по панелям как на GPU!!!
          for (int k = 0; k < npoints; ++k)
          {
              const int indexOfPoint = cntrTree.mortonCodesIdx[k];      //истинный индекс точки наблюдения	
              const double2 p = cntrTree.object[indexOfPoint];          //координаты центра панели
              const double4& pnl = cntrTree.gabForLeaves[indexOfPoint]; //начало и конец панели наблюдения
              const double2 beg{ pnl[0], pnl[1] };                      //отдельно начало
              const double2 end{ pnl[2], pnl[3] };                      //отдельно конец
              const double2 rPan = end - beg;                           //направляющий вектор панели
              const double dlen2 = rPan.length2();                      //длина панели в квадрате
              const double idlen = 1.0 / sqrt(dlen2);                   //обратная длина
              const double2 tau = idlen * rPan;                         //вектор касательной к панели, направленный от начала к концу
              double val = 0.0;                                         //обнуление результата (константная составляющая скорости)
              double vallin = 0.0;                                      //обнуление результата (линейная составляющая скорости)
              int posStack[maxDepth];
              int nodeStack[maxDepth];
              std::array<double2, orderAlignment> Eloc;                 //коэффициенты локальных разложений

              int depth = 0;
              posStack[0] = 0;
              nodeStack[0] = nnodes - 1;
              while (depth >= 0)
              {
                  int pd = posStack[depth];                             //0 или 1 --- какого потомка ячейки nd обходим (левого или правого)   
                  int nd = nodeStack[depth];                            //индекс родительской ячейки дерева, которую обходим, и которая находится на верхушке стека;  

                  int2 chBoth = child[indexSort[(nnodes - 1) - nd]];    //индексы обоих потомков узла nd
                  while (pd < 2)
                  {
                      const int chd = (pd == 0) ? chBoth.first : chBoth.second;
                      ++pd;
                      posStack[depth] = pd;
                      
                      const bool isVortex = (chd >= nbodies);           //признак того, что обрабатывамая ячейка --- лист (в правой половине дерева в порядке karrasOrder)

                      int n;                                            //для хранения индекса вихря или индекса ячейки
                      int srtT = -1;                                    //истинный индекс внутреннего узла дерева

                      double2 ps;                                       //координаты влияющего вихря, если обрабатываемая вершина --- лист, или центра влияющей ячейки, если обрабатывается внутренняя ячейка дерева
                      double gm = 0.0;                                  //циркуляция вихря, если обрабатываемая вершина --- лист, или 0-й мультипольный момент (суммарная циркуляция вихрей) если обрабатывается внутренняя ячейка дерева

                      double sumSide2 = 0.0;                            //габариты обрабатываемой ячейки

                      const Point2D* mom = nullptr;                     //указатель на набор мультипольных моментов влияющей ячейки
                      
                      if (isVortex)                                     //если лист
                      {
                          n = chd - nbodies;                            //номер вихря в Мортоновском порядке
                          const int vortexIndex = mortonCodesIdx[n];    //истинный номер вихря
                            
                          ps = object[vortexIndex];                     //координаты влияющего вихря
                          gm = gamma[vortexIndex];                      //циркуляция влияющего вихря

                          sumSide2 = 0.0;                               //листовая ячейка (вихрь) размера не имеет, т.к. является точкой
                      }//if (isVortex)
                      else
                      {
                          srtT = indexSortT[chd];                       //номер внутреннего узла в "развернутом порядке burtscherOrder" (когда корень --- 0-й), отвечающий узлу chd
                          n = (nnodes - 1) - srtT;                      //номер внутреннего узла в порядке burtscherOrder (когда корень --- последний)
                          ps = center[chd];                             //координаты центра внутреннего узла --- влияющей ячейки

                          const double4& gab = lowerupper[chd];                     //габарит влияющей ячейки
                          const double sumSide = gab[2] - gab[0] + gab[3] - gab[1]; //сумма габаритов (полупериметр)
                          sumSide2 = sumSide * sumSide;                             //квадрат суммы габаритов влияющей ячейки    
                      }//if (!isVortex)

                      const double2 dr = p - ps;                        //радиус-вектор из центра влияющей ячейки в точку наблюдения  
                      const double r2 = dr[0] * dr[0] + dr[1] * dr[1];  //квадрат модуля предыдущего

                      // Для вихря всегда считаем точное влияние. Для ячейки проверяем MAC.
                      if (isVortex || (sumSide2 + dlen2) * itolsq < r2)
                      {
                          if (isVortex)                     //если лист --- считаем напрямую
                          {
                              const double2 ss = ps - beg;
                              const double2 pp = ps - end;

                              const double alpha = atan2(pp[0] * ss[1] - pp[1] * ss[0], pp[0] * ss[0] + pp[1] * ss[1]);

                              if (r2 > 1e-20)
                                  val -= gm * alpha;

                              if (scheme) //если схема Т1
                              {
                                  const double txx = tau[0] * tau[0];
                                  const double txy = tau[0] * tau[1];
                                  const double tyy = tau[1] * tau[1];

                                  double2 u1;

                                  u1[0] = (pp[0] + ss[0]) * (txx - tyy) + 2.0 * (pp[1] + ss[1]) * txy;

                                  u1[1] = (pp[1] + ss[1]) * (tyy - txx) + 2.0 * (pp[0] + ss[0]) * txy;

                                  const double lambda = 0.5 * log(ss.length2() / pp.length2());

                                  const double tempVelLin = gm * (alpha * (u1[0] * tau[0] + u1[1] * tau[1]) +
                                                                 lambda * (-u1[1] * tau[0] + u1[0] * tau[1]));

                                  vallin -= 0.5 * idlen * tempVelLin;
                              }//if(linScheme)
                          }//if (isVortex)
                          else              //если не лист --- мультипольное приближение
                          {
                              mom = moms.data() + srtT * orderAlignment;
                              //std::vector<double2> Eloc(order, double2{ 0.0, 0.0 });
                              std::fill(Eloc.begin(), Eloc.begin() + order, double2{ 0.0, 0.0 }); //обнуляем коэффициенты локального разложения 

                              double2 thetaLoc = (1.0 / r2) * (p - ps);

                              // Коэффициенты локального разложения
                              for (int q = 0; q < order; ++q)
                              {
                                  for (int s = q; s >= 0; --s)
                                  {
                                      const double sign = ((2 * q - s) & 1) ? -1.0 : 1.0;

                                      Eloc[s] += sign * ifac[q - s + 1] * multzA(thetaLoc, mom[q - s]);
                                  }

                                  thetaLoc = ((q + 1.0) / r2) * multz(thetaLoc, p - ps);
                              }

                              // Интегрирование локального разложения по панели для T0 и T1
                              double2 v = 0.5 * Eloc[0];

                              double2 vL{ 0.0, 0.0 };

                              double2 kp = rPan;
                              double2 mulP = kp;

                              double2 taudL = (0.5 / dlen2) * rPan;

                              const double2 taudLc = taudL;

                              for (int kk = 1; kk < order; ++kk)
                              {
                                  mulP = multz(mulP, kp);

                                  taudL *= 0.5;
                                  
                                  if (!(kk & 1))              // Константная проекционная функция
                                    v += ifac[kk + 2] * multz(Eloc[kk], multzA(taudL, mulP));                                  
                                  else if(scheme)             // Линейная проекционная функция 
                                    vL += ifac[kk + 3] * multz(Eloc[kk], multzA(multz(taudL, taudLc), multz(mulP, (kk + 1.0) * rPan)));
                              }

                              val += 2.0 * (-v[1] * rPan[0] + v[0] * rPan[1]);
                              if(scheme)
                                  vallin += 2.0 * (-vL[1] * rPan[0] + vL[0] * rPan[1]);
                          }
                      }//MAC                
                      else //ячейка близкая --- спускаемся по дереву
                      {
                          if (depth + 1 < maxDepth)
                          {
                              if (pd == 1)          //если это была обработка левого потомка, и по нему пришлось идти вниз по дереву
                              {
                                  posStack[depth] = 1;
                                  nodeStack[depth] = nd;

                                  ++depth;
                              }

                              nd = n;
                              pd = 0;

                              chBoth = child[indexSort[(nnodes - 1) - nd]];
                          }
                      }
                  }//while (pd < 2)
                  --depth;          // Оба потомка обработаны -> снимаем узел со стека
              }//while (depth >= 0)

              const double cf = IDPI * idlen;

              rhs[indexOfPoint] = cf * val;
              if(scheme)
                  rhsLin[indexOfPoint] = cf * vallin;
          }//for npoints
      }

      float t2 = (float)omp_get_wtime();
      return (t2 - t1);
  }//DownwardTraversalVorticesToPanels(...)

}//namespace VM2D