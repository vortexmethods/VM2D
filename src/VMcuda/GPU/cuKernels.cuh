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
| File name: cuKernels.cuh                                                    |
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
\brief Заголовки интерфейса с CUDA
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 29 августа 2023 г.
*/

#ifndef CUKERNELS_CUH_
#define CUKERNELS_CUH_


#include "types.cuh"

#ifdef __CUDACC__
#include <cuda.h>
#endif

#include <vector>

//#ifdef __CUDACC__
#include "Point2D.h"
#include "Vortex2D.h"
//#endif

#define THREADS1 32
#define THREADS2 512
#define THREADS3 32
#if __CUDA_ARCH__ >= 800
#define THREADS4 352
#else
#define THREADS4 192
#endif

#define THREADS5 1024
#define THREADS5W2W 1024
#define THREADS5S2W 1024

#define THREADS5rhs 64
#define THREADS5slae 32
#define THREADS5I1I2 32
#define THREADS5I0I3 32

#define FACTOR1 6
#define FACTOR2 2
#define FACTOR3 4
#define FACTOR4 4
#define FACTOR5 1
#define FACTOR5rhs 1
#define FACTOR5slae 1
#define FACTOR5I1I2 1
#define FACTOR5I0I3 1


namespace BHcu
{
    void CudaSelect(int dev);
    
    void setBlocks(int& blocks_);

    void* cudaNew(int n, size_t sizeType);

    void cudaDelete(void* cudaPtr, int code);

    void cudaCopyVecToDevice(void* hostPtr, void* cudaPtr, size_t n, size_t typeSize);

    void cudaCopyVecFromDevice(void* cudaPtr, void* hostPtr, size_t n, size_t typeSize);


    void CudaTest(const char* msg);

    void KernelsOptimization();

    /******************************************************************************/
    /*** initialize memory ********************************************************/
    /******************************************************************************/

    float cuInitializationKernel();


    /******************************************************************************/
    /*** compute center and radius ************************************************/
    /******************************************************************************/

    /// \brief Построение габаритного прямоугольника
    /// 
    /// \param[out] Mposl указатель на массив на device, куда записываются координаты центров внутренних ячеек, заполняется только нулевая ячейка(корень)
    /// \param[out] maxrl указатель на пару чисел на device, куда записываются координаты правого верхнего угла габаритного прямоугольника
    /// \param[out] minrl указатель на пару чисел на device, куда записываются координаты левого нижнего угла габаритного прямоугольника
    /// \param[in] nbodies количество вихрей
    /// \param[in] ptrl указатель на массив на device, где хранятся объекты 
    /// \param[in] sizeOfElement размер элемента в массиве ptrl в байтах
    /// \param[in] offsetOfPointInElement смещение Point2D в элементе массива ptrl в байтах
    /// 
    /// \return время исполнения
    float McuBoundingBoxKernelFree(
        realPoint* __restrict Mposl,
        realPoint* __restrict maxrl,
        realPoint* __restrict minrl,
        int nbodies,
        const void* __restrict ptrl,
        int sizeOfElement,
        int offsetOfPointInElement);
    

	float McuBoundingBoxKernel(
		CUDApointers ptr,
        int nbodiesd,
		/*const realVortex* __restrict vtxd*/
        const void* __restrict vtxd,
        int sizeOfElement,
        int offsetOfPointInElement);


	/******************************************************************************/
	/*** Morton codes *************************************************************/
	/******************************************************************************/

    /// \brief Вычисление кодов Мортона
    /// 
    /// \param[out] maxrl указатель на пару чисел на device, куда записываются координаты правого верхнего угла габаритного прямоугольника
    /// \param[out] minrl указатель на пару чисел на device, куда записываются координаты левого нижнего угла габаритного прямоугольника
    /// \param[out] MmortonCodesKeyUnsortl указатель на массив на device, куда записываются несортированные коды Мортона (в том же порядке, что и вихри в массиве vtxl)
    /// \param[out] MmortonCodesIdxUnsortl указатель на массив на device, куда записываются числа по порядку от 0 до nbodies-1
    /// \param[out] MmortonCodesKeyl указатель на массив на device, куда записываются отсортированные по возрастанию коды Мортона (в порядке обхода Z-кривой)
    /// \param[out] MmortonCodesIdxl указатель на массив на device, куда записываются правила перестановки кодов Мортона (получается "синхронной" сортировкой MmortonCodesKeyl) 
    /// \param[out] Mrangel указатель на массив на device, куда записываются диапазоны частиц (по индексам в отсортированном массиве MmortonCodesKeyl), содержащихся во внутренних ячейках (заполняется только для корня)
    /// \param[in] nbodies количество вихрей
    /// \param[in] vtxl указатель на массив на device, где хранятся вихри
    /// 
    /// \return время исполнения
    float McuMortonCodesKernelFree(
        const realPoint* __restrict maxrl,
        const realPoint* __restrict minrl,
        int* __restrict MmortonCodesKeyUnsortl,
        int* __restrict MmortonCodesIdxUnsortl,
        int* __restrict MmortonCodesKeyl,
        int* __restrict MmortonCodesIdxl,
        const intPair* __restrict Mrangel,

        int nbodiesd,
        const void* __restrict vtxd,
        int sizeOfElement,
        int offsetOfPointInElement, 
        bool sort);


	float McuMortonCodesKernel(
        CUDApointers ptr,
        int nbodiesd,
        const void* __restrict vtxd,
        int sizeOfElement,
        int offsetOfPointInElement);


    /******************************************************************************/
    /*** Morton Internal nodes build **********************************************/
    /******************************************************************************/
    
    /// \brief Определение топологии Мортоновского дерева 
    /// 
    /// \param[in] nbodies количество вихрей
    /// \param[in] MmortonCodesKeyl указатель на массив на device, куда записываются отсортированные по возрастанию коды Мортона (в порядке обхода Z-кривой)
    /// \param[out] Mparentl указатель на массив на device, куда записываются индексы родителей ячеек 
    /// \param[out] Mchildl указатель на массив на device, куда записываются пары индексов ячеек потомков 
    /// \param[out] Mrangel указатель на массив на device, куда записываются диапазоны частиц (по индексам в отсортированном массиве MmortonCodesKeyl), содержащихся во внутренних ячейках 
    /// 
    /// \return время исполнения
    float McuMortonInternalNodesKernel(
        CUDApointers ptr,
        int nbodiesd);

    /******************************************************************************/
    /*** Morton Internal nodes geometry calculation *******************************/
    /******************************************************************************/

    /// \brief Определение геометрических параметров внутренних ячеек Мортоновского дерева 
    /// 
    /// \param[in] nbodies количество вихрей
    /// \param[in] nnodes размер условного массива, хранящего все дерево (не менее, чем 2*nbodies)
    /// \param[in] MmortonCodesKeyl указатель на массив на device, куда записываются отсортированные по возрастанию коды Мортона (в порядке обхода Z-кривой)
    /// \param[out] Mposl указатель на массив на device, куда записываются координаты центров внутренних ячеек
    /// \param[in] Mrangel указатель на массив на device, где хранятся диапазоны частиц (по индексам в отсортированном массиве MmortonCodesKeyl), содержащихся во внутренних ячейках 
    /// \param[out] MlevelUnsortl указатель на массив на device, куда записываются уровни внутренних ячеек мортоновского дерева
    /// \param[out] MlevelSortl указатель на массив на device, куда записываются отсортированные уровни внутренних ячеек мортоновского дерева
    /// \param[out] MindexUnsortl указатель на массив на device, куда записываются числа по порядку от 0 до nbodies-2
    /// \param[out] MindexSortl указатель на массив на device, куда записываются правила перестановки MlevelSortl (получается "синхронной" сортировкой MlevelSortl) 
    /// \param[out] MindexSortTl указатель на массив на device, куда записываются правила обратной перестановки MlevelSortl
    /// 
    /// \return время исполнения
    float McuMortonInternalCellsGeometryKernel(
        CUDApointers ptr,
        int nbodiesd,
        int nnodesd);
   


    /******************************************************************************/
    /*** build tree ***************************************************************/
    /******************************************************************************/

    /// \brief Обнуление необходимых параметров 
    /// 
    /// \param[in] nnodes размер условного массива, хранящего все дерево (не менее, чем 2*nbodies)
    /// \param[in] nbodies количество вихрей
    /// \param[out] massl указатель на массив на device, где присваиваются -1 внутренним узлам  
    /// \param[out] momsl указатель на массив на device, где обнуляются мультипольные моменты всех ячеек
    /// 
    /// \return время исполнения
    float cuClearKernel2(
        CUDApointers ptr,
        const int order,
        int nnodesd, int nbodiesd);


    float cuAABBKernel2(
        CUDApointers ptr,
        const int nnodesd, const int nbodiesd,
        /*const realVortex* __restrict vtxd*/
        const void* __restrict ptrd,
        int sizeOfElement,
        int offsetOfPointInElement,
        bool bodiesZeroSize,
        const double* XYXY);




    /******************************************************************************/
    /*** compute center of mass ***************************************************/
    /******************************************************************************/

    /// \brief Вычисление мультипольных моментов 
    /// 
    /// \param[in] nnodes размер условного массива, хранящего все дерево (не менее, чем 2*nbodies)
    /// \param[in] nbodies количество вихрей
    /// \param[in] Mchildl указатель на массив на device, куда записываются пары индексов ячеек потомков 
    /// \param[out] massl указатель на массив на device, где присваиваются массы (единица для точечного вихря, число вихрей для ячейки)
    /// \param[out] momsl указатель на массив на device, куда записываются мультипольные моменты всех ячеек
    /// \param[in] vtxl указатель на массив на device, где хранятся вихри
    /// \param[in] MmortonCodesIdxl указатель на массив на device, где хранятся правила перестановки кодов Мортона 
    /// \param[in] Mposl указатель на массив на device, где хранятся координаты центров внутренних ячеек
    /// \param[in] MindexSortl указатель на массив на device, где хранятся правила перестановки массива уровней дерева (MlevelSortl)
    /// \param[in] MindexSortTl указатель на массив на device, где хранятся правила обратной перестановки массива уровней дерева (MlevelSortl)
    /// 
    /// \return время исполнения
    float cuSummarizationKernel2(
        CUDApointers ptr,
        const int order,
        const int nnodesd, const int nbodiesd,
        /*const realVortex* __restrict vtxd*/
        const double* __restrict vtxd,
        int objType);

    /******************************************************************************/
    /*** compute force ************************************************************/
    /******************************************************************************/

    /// \brief Вычисление скоростей
    /// 
    /// \param[in] nnodes размер условного массива, хранящего все дерево (не менее, чем 2*nbodies)
    /// \param[in] nbodies количество вихрей
    /// \param[in] itolsq параметр близости
    /// \param[in] epssq радиус вихря
    /// \param[in] Mchildl указатель на массив на device, где хранятся пары индексов ячеек потомков 
    /// \param[in] momsl указатель на массив на device, где хранятся мультипольные моменты всех ячеек
    /// \param[in] vtxl указатель на массив на device, где хранятся вихри
    /// \param[in] MmortonCodesIdxl указатель на массив на device, где хранятся правила перестановки кодов Мортона 
    /// \param[in] Mposl указатель на массив на device, где хранятся координаты центров внутренних ячеек
    /// \param[in] MindexSortl указатель на массив на device, где хранятся правила перестановки массива уровней дерева (MlevelSortl)
    /// \param[in] MindexSortTl указатель на массив на device, где хранятся правила обратной перестановки массива уровней дерева (MlevelSortl)
    /// \param[out] vell указатель на массив на device, куда записываются скорости вихрей
    /// \param[in] Msizel указатель на массив на device, где хранятся пары чисел - размеры внутренних ячеек по горизонтали и вертикали
    /// 
    /// \return время исполнения
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
        size_t nAfls, size_t* nVtxs, double** ptrVtxs);

    float cuForceCalculationKernelFromPanels2points(
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real itolsqd, real epssqd,
        const double* dev_ptr_r,
        const void* __restrict vtxd,
        int sizeOfElement,
        int offsetOfPointInElement,
        const int* __restrict MmortonCodesIdxl,
        int npointsd, const realVortex* pointsd,
        realPoint* __restrict veld,
        int schemeType);

    float cuI1I2CalculationKernel2(
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real epssqd,
        const realVortex* __restrict vtxd,
        real* __restrict I1d,
        realPoint* __restrict I2d,
        real* __restrict epsastd);

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
        real* __restrict visstr);

    float McuVerticesToControlPoints(int nTotPan, const double* dev_ptr_pt, double* pointsl);

    float McuVerticesToControlPointsVortexes(int nTotPan, const double* dev_ptr_pt, double* pointsl, const double* gaml);

    float McuVerticesAndSheetsToPanelPoints(int npnli, const double* dev_ptr_pt, 
        const double* dev_ptr_freeVortexSheet, const double* dev_ptr_freeVortexSheetLin,
        const double* dev_ptr_attachedVortexSheet, const double* dev_ptr_attachedVortexSheetLin,
        const double* dev_ptr_attachedSourceSheet, const double* dev_ptr_attachedSourceSheetLin,
        double* panelPoints, int schemeType);

    float cuRhsCalculationKernel(
        CUDApointers ptr,
        int order,
        int nnodesd, int nbodiesd,
        real itolsqd, 
        const realVortex* __restrict vtxd,
        const int* __restrict MmortonCodesIdxl,
        realPoint* __restrict El,
        int nTotPan, const real* dev_ptr_pt, const real* pointsd,
        real* __restrict veld,
        real* __restrict vellind
        );

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
        real* __restrict vellind);


    /******************************************************************************/
    /*** compute force (direct) ***************************************************/
    /******************************************************************************/
    float cuForceDirectCalculationKernel(
        const int nbodiesd,
        const real epssqd,
        const realVortex* __restrict vtxd,        
        volatile realPoint* __restrict veld);



}//namespace BHcu

#endif