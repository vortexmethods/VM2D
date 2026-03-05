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

        const unsigned int ki = mortonCodesKeyUnsort[i];
        const unsigned int kj = mortonCodesKeyUnsort[j];

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
            offsetOfPointInElement = 0;
            break;

        case object_T::point3:
            sizeOfElement = sizeof(Vortex2D);
            offsetOfPointInElement = 0;
            break;

        case object_T::panel:
            switch (treeType)
            {
            case tree_T::aux:
                sizeOfElement = sizeof(double) * 6;
                offsetOfPointInElement = 0;
                break;
            case tree_T::contr:
                sizeOfElement = sizeof(Point2D);
                offsetOfPointInElement = 0;
                break;
            case tree_T::vortex:
            case tree_T::source:
                sizeOfElement = sizeof(double) * 12;
                offsetOfPointInElement = 0;
                break;
            }
            break;
        }
        */
    };

    CpuTreeInfo::~CpuTreeInfo() {};

    void CpuTreeInfo::Update(const std::vector<Vortex2D>& vtx, double eps)
    {
        object.resize(vtx.size());
        gamma.resize(vtx.size());
        gabForLeaves.resize(vtx.size());
        mortonCodesKeyUnsort.resize(vtx.size());
        mortonCodesIdxUnsort.resize(vtx.size());

        if (treeType != tree_T::contr)
        {
            levelUnsort.resize(vtx.size() - 1);
            indexUnsort.resize(vtx.size() - 1);
            range.resize(vtx.size() - 1);
            parent.resize(vtx.size() - 1);
            child.resize(vtx.size() - 1);
            mass.resize(vtx.size() - 1);
        }

        for (size_t v = 0; v < vtx.size(); ++v)
        {
            Point2D r = vtx[v].r();
            object[v] = r;
            gamma[v] = vtx[v].g();
            gabForLeaves[v] = Point4D({ r[0] - eps, r[1] - eps, r[0] + eps, r[1] + eps });
        }
    }


    float CpuTreeInfo::Build()
    {
        float time = 0.0f;
        if (object.size() > 0)
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
            for (int bdy = 0; bdy < (int)object.size(); ++bdy)
            {
                Point2D rScaled = twoPowCodeLength * ((object[bdy] - rcen) * quadSideFactor + 0.5 * Point2D{ rbound, rbound });

                unsigned int xx = ExpandBits((unsigned int)rScaled[0]);
                unsigned int yy = ExpandBits((unsigned int)rScaled[1]);
                mortonCodesKeyUnsort[bdy] = yy | (xx << 1);
                mortonCodesIdxUnsort[bdy] = bdy;
            }

            printf("SORT morton codes!!!\n");
            exit(-1);


            //treeMortonInternalNodes
            if (treeType != tree_T::contr)
            {
                for (int i = 0; i < (int)object.size() - 1; ++i)
                {
                    int codei = mortonCodesKeyUnsort[i];

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
                    int childLeft = (Mmin == gammaPos) * (int)object.size() + left;
                    range[childLeft] = { Mmin, gammaPos };
                    parent[childLeft] = i;

                    //                -                         
                    int childRight = (Mmax == gammaPos + 1) * (int)object.size() + right;
                    range[childRight] = { gammaPos + 1, Mmax };
                    parent[childRight] = i;

                    child[i] = { childLeft, childRight };
                }
            }
            printf("SORT internal nodes!!!\n");
            exit(-1);
        }
        return time;
    }


    void CpuTreeInfo::summ12()
    {
        using double4 = Point4D;
        using double2 = Point2D;
        using int2 = std::pair<int, int>;

        int i, j, ch, flag;

        double4 lu[2];

        double2 mom0;
        double2 mom1;
        double2 mom2;
        double2 mom3;
        double2 mom4;
        double2 mom5;
        double2 mom6;
        double2 mom7;
        double2 mom8;
        double2 mom9;
        double2 mom10;
        double2 mom11;

        double2 cen, dr;

        int cm;
        int m[2];
        const int nnodes = 2 * (int)object.size() - 1;
        const int nbodies = (int)object.size();

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


            j = 0;
            flag = 0;
            // iterate over all cells assigned to thread
            while (flag == 0)
            {
                j = 2;
                int srt = indexSort[(nnodes - 1) - k];
                int2 chdPair = child[srt];

                int chdSorted[2];

                for (i = 0; i < 2; i++) {
                    int chd = i * chdPair.second + (1 - i) * chdPair.first;   // i==0 => .x;  i==1 => .y

                    ch = (chd >= nbodies) ? chd - nbodies : (chdSorted[i] = (nnodes - 1) - indexSortT[chd]);

                    if ((chd >= nbodies) || (mass[nnodes - 1 - ch] >= 0))
                        j--;
                }

                if (j == 0)
                {
                    // all children are ready
                    const int kch = ((nnodes - 1) - k) * orderAlignment;
                    cm = 0;

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
                        }
                        else
                        {
                            ch = (nnodes - 1) - chdSorted[i];
                            lu[i] = lowerupper[chd];
                        }
                    }//for i

                    lowerupper[srt] = double4{
                        ::fmin(lowerupper[0][0], lowerupper[1][0]), ::fmin(lowerupper[0][1], lowerupper[1][1]),
                        ::fmax(lowerupper[0][2], lowerupper[1][2]), ::fmax(lowerupper[0][3], lowerupper[1][3])
                    };



                    cen = center[srt] = double2{ 0.5 * (lowerupper[srt][0] + lowerupper[srt][2]),
                                                 0.5 * (lowerupper[srt][1] + lowerupper[srt][3]) };

                    const double2 zero = { 0.0, 0.0 };
                    double2 momh0 = zero;
                    double2 momh1 = zero;
                    double2 momh2 = zero;
                    double2 momh3 = zero;
                    double2 momh4 = zero;
                    double2 momh5 = zero;
                    double2 momh6 = zero;
                    double2 momh7 = zero;
                    double2 momh8 = zero;
                    double2 momh9 = zero;
                    double2 momh10 = zero;
                    double2 momh11 = zero;

                    for (i = 0; i < 2; i++)
                    {
                        //computation of ch = child[k*2+i]
                        const int chd = i * chdPair.second + (1 - i) * chdPair.first;
                        if (chd >= nbodies)
                        {
                            ch = chd - nbodies;
                            const int sortedBody = mortonCodesIdx[ch];

                            if (objectType == object_T::point3)
                            {
                                mom0 = double2{ gamma[sortedBody], 0.0 };
                                mom1 = mom2 = mom3 = mom4 = mom5 = mom6 = mom7 = mom8 = mom9 = mom10 = mom11 = double2{ 0.0, 0.0 };
                                double2 pos = object[sortedBody];
                                dr = pos - cen;
                                m[i] = 1;
                            } //objectType==0
/*
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
*/
                        }
                        else
                        {
                            const int srtT = indexSortT[chd];
                            ch = (nnodes - 1) - srtT;

                            const int nch = srtT * orderAlignment;

                            mom0 = double2{ moms[nch + 0][0], (double)0 };
                            mom1 = moms[nch + 1];
                            mom2 = moms[nch + 2];
                            mom3 = moms[nch + 3];
                            mom4 = moms[nch + 4];
                            mom5 = moms[nch + 5];
                            mom6 = moms[nch + 6];
                            mom7 = moms[nch + 7];
                            mom8 = moms[nch + 8];
                            mom9 = moms[nch + 9];
                            mom10 = moms[nch + 10];
                            mom11 = moms[nch + 11];

                            dr = center[chd] - cen;
                            m[i] = mass[nnodes - 1 - ch];
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

                    momh0[1] = 0;
                    moms[kch + 0] = momh0;
                    moms[kch + 1] = momh1;
                    moms[kch + 2] = momh2;
                    moms[kch + 3] = momh3;
                    moms[kch + 4] = momh4;
                    moms[kch + 5] = momh5;
                    moms[kch + 6] = momh6;
                    moms[kch + 7] = momh7;
                    moms[kch + 8] = momh8;
                    moms[kch + 9] = momh9;
                    moms[kch + 10] = momh10;
                    moms[kch + 11] = momh11;

                    cm += m[0] + m[1];

                    flag = 1;
                }

                //__threadfence();

                if (flag != 0)
                {
                    mass[nnodes - 1 - k] = cm;
                    //k += inc;
                    //flag = 0;
                }
            }//while flag==0
        }//for k
    }//summ12()


    float CpuTreeInfo::UpwardTraversal(int order)
    {
        float time = 0.0f;
        if (object.size() > 0)
        {
            //treeClearKernel
            mass.assign(object.size() - 1, -1);
            
            if (treeType != tree_T::contr && treeType != tree_T::aux)
                //treeSummarization;
                summ12();
            else if (treeType == tree_T::aux)
                //time += treeCalcAABBWrapper(*this);
            {
                printf("To do: treeCalcAABBWrapper\n");
                exit(-1);
            }
        }
        return time;
    }

}