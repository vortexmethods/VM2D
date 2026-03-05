/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: cpuTreeInfo.h                                                   |
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
\brief Заголовок класса дерева для реализации быстрых алгоритмов на CPU
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef CPUTREEINFO_H
#define CPUTREEINFO_H

#include "Vortex2D.h"
#include "Gpudefs.h"

namespace VM2D
{    


    ///*!    
    //\brief Структура, хранящая данные и указатели на массивы на GPU для оптимизации итерационного решения СЛАУ на GPU
    //
    //\author Марчевский Илья Константинович
    //\author Сокол Ксения Сергеевна
    //\author Рятина Евгения Павловна
    //\author Колганова Александра Олеговна
    //
    //\Version 1.14
    //\date 6 марта 2026 г.
    //*/
    //struct infoForMatVecMul
    //{
    //    std::vector<int> nClosePanels;    ///количество листьев дерева (панелей), находящихся в ближней зоне
    //    std::vector<int> nFarCells;       ///количество ячеек, находящихся в дальней зоне

    //    std::vector<int> closePrefixSum;  ///префиксная сумма для ближней зоны
    //    std::vector<int> farPrefixSum;    ///префиксная сумма для дальней зоны

    //    Point2D* i00D; ///указатель на GPU, хранящий влияние панелей из ближней зоны
    //    Point2D* i01D; ///аналог i00D, но используется для кусочно-линейной схемы
    //    Point2D* i10D; ///аналог i00D, но используется для кусочно-линейной схемы
    //    Point2D* i11D; ///аналог i00D, но используется для кусочно-линейной схемы

    //    int* nClosePanelsD;     ///количество листьев дерева (панелей), находящихся в ближней зоне
    //    int* nFarCellsD;        ///количество ячеек, находящихся в дальней зоне

    //    int* closeCellsIdxD;    ///индексы ячеек ближней зоны
    //    int* farCellsIdxD;      ///индексы ячеек дальней зоны

    //    int* closePrefixSumD;  ///префиксная сумма для ближней зоны
    //    int* farPrefixSumD;    ///префиксная сумма для дальней зоны
    //};


    /*!
    \brief Класс, определяющий структуры дерева для быстрых методов при их реализации на CPU

    \author Марчевский Илья Константинович
    \author Сокол Ксения Сергеевна
    \author Рятина Евгения Павловна
    \author Колганова Александра Олеговна

    \Version 1.14
    \date 6 марта 2026 г.
    */

    class CpuTreeInfo
    {
    public:
        //int nObject;        
        //int nNode;

        //bool duplicateObj;
        //double* objectD;
        std::vector<Point2D> object;
        std::vector<double> gamma;

        tree_T treeType;
        object_T objectType;         
        scheme_T schemeType;

        //int sizeOfElement;
        //int offsetOfPointInElement;

        std::vector<int> mass;	//массы (единица для точечного вихря, число вихрей для ячейки)

        Point2D maxr, minr;       //габаритный прямоугольник
        std::vector<Point2D> moms;//мультипольные моменты всех ячеек; хранятся в виде <mom_0x, mom_0y=0, mom_1x, mom_1y, ..., mom_px, mom_py>, <для второй ячейки> ...
        std::vector<Point2D> ED;  //к-ты локальных разложений

        //For Morton tree
        std::vector<int> mortonCodesKeyUnsort;
        std::vector<int> mortonCodesKey;

        std::vector<int> mortonCodesIdxUnsort; //0 1 2 3 ... nbodies-1		
        std::vector<int> mortonCodesIdx;

        std::vector<int> levelUnsort;
        std::vector<int> levelSort;

        std::vector<int> indexUnsort;   //0 1 2 3 ... nbodies-2
        std::vector<int> indexSort; 
        std::vector<int> indexSortT;

        std::vector<Point2D> center;       //Положения внутренних узлов в дерева Карраса
        std::vector<Point4D> lowerupper;   //Левый нижний и правый верхний углы ячейки
        std::vector<Point4D> gabForLeaves; //Левый нижний и правый верхний углы листа

        std::vector<int> parent;     //Номер ячейки-родителя
        std::vector<std::pair<int, int>> child;  //Потомки внутренних ячеек (в одной ячейке храним сразу два целых числа)
        std::vector<std::pair<int, int>> range;  //Диапазон частиц во внутренней ячейке (в одной ячейке храним сразу два целых числа)

        //void* sortObjectsBufferD;
        //int sortObjectsBufferSizeInBytes;

        //void* sortInternalCellsBufferD;
        //int sortInternalCellsBufferSizeInBytes;

        //const int nBlock;
        //int reservedMemorySizeItems;

        //infoForMatVecMul matVecMulInfo;

        CpuTreeInfo(tree_T treeType_, object_T objectType_, scheme_T schemeType_);
        ~CpuTreeInfo();

        void Update(const std::vector<Vortex2D>& vtx, double eps);
//        void UpdatePanelGeometry(int nObject_, double4* gabForLeavesD_);

//        void UpdatePanelFreeVortexIntensity(const double* dev_ptr_freeVortexSheet, const double* dev_ptr_freeVortexSheetLin);//функция нулит attached vortex sheet
//        void UpdatePanelAttachedVortexIntensity(const double* dev_ptr_attachedVortexSheet, const double* dev_ptr_attachedVortexSheetLin);//функция нулит free vortex sheet
//        void UpdatePanelAttachedSourceIntensity(const double* dev_ptr_attachedSourceSheet, const double* dev_ptr_attachedSourceSheetLin);
//        void UpdatePanelFreeAndAttachedVortexIntensity(const double* dev_ptr_freeVortexSheet, const double* dev_ptr_freeVortexSheetLin, const double* dev_ptr_attachedVortexSheet, const double* dev_ptr_attachedVortexSheetLin);

//        void MemoryAllocateForGMRES();
//        void MemoryFreeForGMRES();

//        void MemoryAllocate(int nCudaObject);
        
//        bool IsInitialized() const;

        int Delta(int i, int j) const;
        float Build();
        float UpwardTraversal(int order);
        void summ12();

//        float DownwardTraversalVorticesToPoints(CudaTreeInfo& cntrTree, Point2D* velD, double* epsastD, double eps2, double theta, int order, bool calcRadius);
//        float DownwardTraversalVorticesToPanels(CudaTreeInfo& cntrTree, double* rhsD, double* rhsLinD, double theta, int order);
//        float DownwardTraversalPanelsToPoints(CudaTreeInfo& cntrTree, Point2D* velD, double eps2, double theta, int order);
//        float DownwardTraversalGMRES(double* resD, double* resLinD, double theta, int order, int iter);

//        float I1I2CalculationWrapper(double minRd, double* __restrict I1D, Point2D* __restrict I2D, double* __restrict epsastD);
//        float I0I3CalculationWrapper(double minRd, float* __restrict I0D, Point2Df* __restrict I3D, double* __restrict epsastD, const double* __restrict meanEpsD, int nPan, double* __restrict panD, double* __restrict visstrD);
//        void RadixSortMortonCodes(int beginBit = 0, int endBit = 2 * codeLength);
//        void RadixSortInternalCells(int beginBit = 0, int endBit = 2 * codeLength);
    };

}


#endif