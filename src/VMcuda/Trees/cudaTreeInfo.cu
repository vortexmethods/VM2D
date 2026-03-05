/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: cudaTreeInfo.cu                                                  |
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
\brief Реализация класса дерева для реализации быстрых алгоритмов на CUDA
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#include "cudaTreeInfo.h"
#include "treeKernels.cuh"
#include "templateKernels.cuh"
#include "Gpudefs.h"

//#include "cub/cub.cuh"   // or equivalently 
#include <cub/device/device_radix_sort.cuh>

namespace BHcu
{
    __global__
        void FillGabForPointsKernel(int n, const double* __restrict pntD, int sizeOfElement, double eps, double4* __restrict gabD)
    {
        size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < n)
        {
            double2 p;
            double3 vt;

            switch (sizeOfElement)
            {
            case 16:
                p = *(double2*)(&pntD[2 * i]);
                break;

            case 24:
                vt = *(double3*)(&pntD[3 * i]);
                p.x = vt.x;
                p.y = vt.y;
                break;

            default:
                printf("unknown size!\n");
            }

            gabD[i] = make_double4(p.x - eps, p.y - eps, p.x + eps, p.y + eps);
        }
    }

    __global__
        void FillObjectForControlPanelsKernel(int n, const double4* __restrict begEndD, double2* __restrict centresD)
    {
        size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < n)
        {
            double4 gab = begEndD[i];
            centresD[i] = make_double2(0.5 * (gab.x + gab.z), 0.5 * (gab.y + gab.w));
        }
    }


    __global__
        void FillObjectForInfluencingPanelsKernel(int n, const double4* __restrict begEndD, tree_T treeType, double* __restrict panelPoints)
    {
        const int stride = (treeType == tree_T::aux ? 6 : 12);

        size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i < n)
        {
            double4 gab = begEndD[i];
            panelPoints[stride * i + 0] = 0.5 * (gab.x + gab.z);
            panelPoints[stride * i + 1] = 0.5 * (gab.y + gab.w);

            panelPoints[stride * i + 2] = gab.x;
            panelPoints[stride * i + 3] = gab.y;
            panelPoints[stride * i + 4] = gab.z;
            panelPoints[stride * i + 5] = gab.w;     
        }
    }//FillObjectForInfluencingPanelsKernel(...)


    __global__
        void FillFreeVortexIntensityForPanelsKernel(int n, const double4* __restrict begEndD,
            const double* __restrict dev_ptr_freeVortexSheet, const double* __restrict dev_ptr_freeVortexSheetLin,
            tree_T treeType, double* __restrict panelPoints)

    {
        const int stride = (treeType == tree_T::aux ? 6 : 12);

        size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        double dx, dy, dl;
        if (i < n)
        {            
            double4 gab = begEndD[i];

            dx = gab.z - gab.x;
            dy = gab.w - gab.y;
            dl = sqrt(dx * dx + dy * dy);

            panelPoints[stride * i + 6] = dl * dev_ptr_freeVortexSheet[i];

            panelPoints[stride * i + 7] = 0.0;

            if (dev_ptr_freeVortexSheetLin != nullptr)
                panelPoints[stride * i + 9] = dl * dev_ptr_freeVortexSheetLin[i];
            else panelPoints[stride * i + 9] = 0.0;
            panelPoints[stride * i + 10] = 0.0;
        }
    }//FillFreeVortexIntensityForPanelsKernel(...)
    
    __global__
        void FillAttachedVortexIntensityForPanelsKernel(int n, const double4* __restrict begEndD,
            const double* __restrict dev_ptr_attachedVortexSheet, const double* __restrict dev_ptr_attachedVortexSheetLin,
            tree_T treeType, double* __restrict panelPoints)

    {
        const int stride = (treeType == tree_T::aux ? 6 : 12);

        size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        double dx, dy, dl;
        if (i < n)
        {
            double4 gab = begEndD[i];

            dx = gab.z - gab.x;
            dy = gab.w - gab.y;
            dl = sqrt(dx * dx + dy * dy);

            panelPoints[stride * i + 6] = 0.0; //free const vortex sheet

            if (dev_ptr_attachedVortexSheet != nullptr)
                panelPoints[stride * i + 7] = dl * dev_ptr_attachedVortexSheet[i];
            else
                panelPoints[stride * i + 7] = 0.0;

            panelPoints[stride * i + 9] = 0.0; //free lin vortex sheet

            if (dev_ptr_attachedVortexSheetLin != nullptr)
                panelPoints[stride * i + 10] = dl * dev_ptr_attachedVortexSheetLin[i];
            else
                panelPoints[stride * i + 10] = 0.0;
        }
    }//FillAttachedVortexIntensityForPanelsKernel(...)

    __global__
        void FillFreeAndAttachedVortexIntensityForPanelsKernel(int n, const double4* __restrict begEndD,
            const double* __restrict dev_ptr_freeVortexSheet, const double* __restrict dev_ptr_freeVortexSheetLin,
            const double* __restrict dev_ptr_attachedVortexSheet, const double* __restrict dev_ptr_attachedVortexSheetLin,
            tree_T treeType, double* __restrict panelPoints)

    {
        const int stride = (treeType == tree_T::aux ? 6 : 12);

        size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        double dx, dy, dl;
        if (i < n)
        {
            double4 gab = begEndD[i];

            dx = gab.z - gab.x;
            dy = gab.w - gab.y;
            dl = sqrt(dx * dx + dy * dy);

            panelPoints[stride * i + 6] = dl * dev_ptr_freeVortexSheet[i];

            if (dev_ptr_attachedVortexSheet != nullptr)
                panelPoints[stride * i + 7] = dl * dev_ptr_attachedVortexSheet[i];
            else
                panelPoints[stride * i + 7] = 0.0;

            if (dev_ptr_freeVortexSheetLin != nullptr)
                panelPoints[stride * i + 9] = dl * dev_ptr_freeVortexSheetLin[i];
            else panelPoints[stride * i + 9] = 0.0;

            if (dev_ptr_attachedVortexSheetLin != nullptr)
                panelPoints[stride * i + 10] = dl * dev_ptr_attachedVortexSheetLin[i];
            else
                panelPoints[stride * i + 10] = 0.0;
        }
    }//FillAttachedVortexIntensityForPanelsKernel(...)

    __global__
        void FillAttachedSourceIntensityForPanelsKernel(int n, const double4* __restrict begEndD,
            const double* __restrict dev_ptr_attachedSourceSheet, const double* __restrict dev_ptr_attachedSourceSheetLin,
            tree_T treeType, double* __restrict panelPoints)

    {
        const int stride = (treeType == tree_T::aux ? 6 : 12);

        size_t i = blockIdx.x * blockDim.x + threadIdx.x;
        double dx, dy, dl;
        if (i < n)
        {
            double4 gab = begEndD[i];

            dx = gab.z - gab.x;
            dy = gab.w - gab.y;
            dl = sqrt(dx * dx + dy * dy);
         
            if (dev_ptr_attachedSourceSheet != nullptr)
                panelPoints[stride * i + 8] = dl * dev_ptr_attachedSourceSheet[i];
            else
                panelPoints[stride * i + 8] = 0.0;

            if (dev_ptr_attachedSourceSheetLin != nullptr)
                panelPoints[stride * i + 11] = dl * dev_ptr_attachedSourceSheetLin[i];
            else
                panelPoints[stride * i + 11] = 0.0;
        }
    }//FillAttachedSourceIntensityForPanelsKernel(...)


    CudaTreeInfo::CudaTreeInfo(int nBlock_, tree_T treeType_, object_T objectType_, scheme_T schemeType_, bool duplicateObj_)
        :
        nObject(0),
        duplicateObj(duplicateObj_),
        objectD(nullptr),
        treeType(treeType_),
        objectType(objectType_),         
        schemeType(schemeType_), 
        sizeOfElement(0),
        offsetOfPointInElement(0),

        massD(nullptr),
        maxrD(nullptr), minrD(nullptr),

        momsD(nullptr),
        ED(nullptr),
        mortonCodesKeyUnsortD(nullptr),
        mortonCodesKeyD(nullptr),
        mortonCodesIdxUnsortD(nullptr),
        mortonCodesIdxD(nullptr),
        levelUnsortD(nullptr),
        levelSortD(nullptr),
        indexUnsortD(nullptr),
        indexSortD(nullptr),
        indexSortTD(nullptr),
        centerD(nullptr),
        lowerupperD(nullptr),
        gabForLeavesD(nullptr),
        parentD(nullptr),
        childD(nullptr),
        rangeD(nullptr),
        sortObjectsBufferD(nullptr),
        sortObjectsBufferSizeInBytes(0),
        sortInternalCellsBufferD(nullptr),
        sortInternalCellsBufferSizeInBytes(0),

        nBlock(nBlock_),
        reservedMemorySizeItems(0)
    {
        cudaMalloc(&maxrD, nBlock * FACTORgab * sizeof(Point2D));
        cudaMalloc(&minrD, nBlock * FACTORgab * sizeof(Point2D));
        
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
                sizeOfElement = sizeof(double)*6;
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
    };

    CudaTreeInfo::~CudaTreeInfo()
    {
        cudaFree(maxrD);
        cudaFree(minrD);

        cudaFree(massD);
        cudaFree(momsD);
        cudaFree(ED);
        cudaFree(mortonCodesKeyUnsortD);
        cudaFree(mortonCodesKeyD);
        cudaFree(mortonCodesIdxUnsortD);
        cudaFree(mortonCodesIdxD);
        cudaFree(levelUnsortD);
        cudaFree(levelSortD);
        cudaFree(indexUnsortD);
        cudaFree(indexSortD);
        cudaFree(indexSortTD);
        cudaFree(centerD);
        cudaFree(lowerupperD);
        cudaFree(parentD);
        cudaFree(childD);
        cudaFree(rangeD);
        cudaFree(sortObjectsBufferD);
        cudaFree(sortInternalCellsBufferD);
        if (duplicateObj)
        {
            cudaFree(objectD);
            cudaFree(gabForLeavesD);
        }
        else
        {
            switch (objectType)
            {
            case object_T::panel:
                cudaFree(objectD);
                break;

            case object_T::point2:
            case object_T::point3:
                cudaFree(gabForLeavesD);
                break;
            }
        }
    }

    void CudaTreeInfo::Update(int nObject_, double* objectD_, double eps)
    {
        if (nObject_ > reservedMemorySizeItems)
        {
            printf("Error in tree allocation!\n");
            exit(-235);
        }

        nObject = nObject_;
        
        if (nObject > 0)
        {
            if (!duplicateObj)
                objectD = objectD_;
            else
                cudaMemcpy(objectD, objectD_, nObject * sizeOfElement, cudaMemcpyDeviceToDevice);

            FillGabForPointsKernel << <(nObject + 31) / 32, 32 >> > (nObject, (double*)objectD, sizeOfElement, eps, (double4*)gabForLeavesD);

            CudaTestError("FillGabForPointsKernel launch failed");

        }//if (nObject > 0)

        nNode = nObject * 2;

        if (nNode < 1024 * nBlock)
            nNode = 1024 * nBlock;
        while ((nNode & (32 - 1)) != 0)  // 32 - это размер варпа
            nNode++;
        nNode--;
    }

    void CudaTreeInfo::UpdatePanelGeometry(int nObject_, double4* gabForLeavesD_)
    {
        if (nObject_ > reservedMemorySizeItems)
        {
            printf("Error in tree Panel allocation!\n");
            exit(-234);
        }

        nObject = nObject_;
        if (!duplicateObj)
            gabForLeavesD = (Point4D*)gabForLeavesD_;
        else
            cudaMemcpy(gabForLeavesD, gabForLeavesD_, nObject * sizeof(double4), cudaMemcpyDeviceToDevice);

        if(treeType == tree_T::contr)
            FillObjectForControlPanelsKernel << <(nObject + 31) / 32, 32 >> > (nObject, (double4*)gabForLeavesD, (double2*)objectD);
        else 
            FillObjectForInfluencingPanelsKernel << <(nObject + 31) / 32, 32 >> > (nObject, (double4*)gabForLeavesD, treeType, (double*)objectD);

        CudaTestError("FillObjectForPanelsKernel launch failed");

        nNode = nObject * 2;
        if (nNode < 1024 * nBlock)
            nNode = 1024 * nBlock;
        while ((nNode & (32 - 1)) != 0)  // 32 - это размер варпа
            nNode++;
        nNode--;
    }

    void CudaTreeInfo::UpdatePanelFreeVortexIntensity(const double*  dev_ptr_freeVortexSheet, const double*  dev_ptr_freeVortexSheetLin)
    {
        FillFreeVortexIntensityForPanelsKernel << <(nObject + 31) / 32, 32 >> > (nObject, (double4*)gabForLeavesD, dev_ptr_freeVortexSheet, dev_ptr_freeVortexSheetLin, treeType, (double*)objectD);

        CudaTestError("FillFreeVortexIntensityForPanelsKernel launch failed");
    } 

    void CudaTreeInfo::UpdatePanelAttachedVortexIntensity(const double* dev_ptr_attachedVortexSheet, const double* dev_ptr_attachedVortexSheetLin)
    {
        FillAttachedVortexIntensityForPanelsKernel << <(nObject + 31) / 32, 32 >> > (nObject, (double4*)gabForLeavesD, dev_ptr_attachedVortexSheet, dev_ptr_attachedVortexSheetLin, treeType, (double*)objectD);

        CudaTestError("FillAttachedVortexIntensityForPanelsKernel launch failed");
    }

    void CudaTreeInfo::UpdatePanelFreeAndAttachedVortexIntensity(const double* dev_ptr_freeVortexSheet, const double* dev_ptr_freeVortexSheetLin, const double* dev_ptr_attachedVortexSheet, const double* dev_ptr_attachedVortexSheetLin)
    {
        FillFreeAndAttachedVortexIntensityForPanelsKernel << <(nObject + 31) / 32, 32 >> > (nObject, (double4*)gabForLeavesD, dev_ptr_freeVortexSheet, dev_ptr_freeVortexSheetLin, dev_ptr_attachedVortexSheet, dev_ptr_attachedVortexSheetLin, treeType, (double*)objectD);

        CudaTestError("FillFreeAndAttachedVortexIntensityForPanelsKernel launch failed");

    }

    void CudaTreeInfo::UpdatePanelAttachedSourceIntensity(const double* dev_ptr_attachedSourceSheet, const double* dev_ptr_attachedSourceSheetLin)
    {
        FillAttachedSourceIntensityForPanelsKernel << <(nObject + 31) / 32, 32 >> > (nObject, (double4*)gabForLeavesD, dev_ptr_attachedSourceSheet, dev_ptr_attachedSourceSheetLin, treeType, (double*)objectD);

        CudaTestError("FillAttachedSourceIntensityForPanelsKernel launch failed");
    }

    void CudaTreeInfo::MemoryAllocate(int nCudaObject)
    {
        //реаллокация
        if (nCudaObject > reservedMemorySizeItems)
        {
            int nObjectUp = nCudaObject;
            int nNodeUp = nObjectUp * 2;
            if (nNodeUp < 1024 * nBlock)
                nNodeUp = 1024 * nBlock;
            while ((nNodeUp & (32 - 1)) != 0)  // 32 - это размер варпа
                nNodeUp++;
            nNodeUp--;

            if (reservedMemorySizeItems > 0)
            {
                cudaFree(massD);

                cudaFree(momsD);
                cudaFree(ED);
                cudaFree(mortonCodesKeyUnsortD);
                cudaFree(mortonCodesKeyD);
                cudaFree(mortonCodesIdxUnsortD);
                cudaFree(mortonCodesIdxD);
                cudaFree(levelUnsortD);
                cudaFree(levelSortD);
                cudaFree(indexUnsortD);
                cudaFree(indexSortD);
                cudaFree(indexSortTD);
                cudaFree(centerD);
                cudaFree(lowerupperD);
                cudaFree(parentD);
                cudaFree(childD);
                cudaFree(rangeD);
                cudaFree(sortObjectsBufferD);
                cudaFree(sortInternalCellsBufferD);
                
                switch (objectType)
                {
                case object_T::panel:
                    cudaFree(objectD);
                    if (duplicateObj)
                        cudaFree(gabForLeavesD);
                    break;

                case object_T::point2:
                case object_T::point3:
                    cudaFree(gabForLeavesD);
                    if (duplicateObj)
                        cudaFree(objectD);
                    break;
                }
            }           

            int nbodies = nObjectUp;
            int nnodes = nNodeUp;


            cudaMalloc(&massD, (nbodies - 1) * sizeof(int));

            if (treeType != tree_T::contr && treeType != tree_T::aux)
                cudaMalloc(&momsD, (nbodies - 1) * orderAlignment * sizeof(Point2D));
            else
                momsD = nullptr;

            ED = nullptr;

            ///For MortonTree
            cudaMalloc(&mortonCodesKeyUnsortD, nbodies * sizeof(int));
            cudaMalloc(&mortonCodesKeyD, nbodies * sizeof(int));
            cudaMalloc(&mortonCodesIdxUnsortD, nbodies * sizeof(int));
            cudaMalloc(&mortonCodesIdxD, nbodies * sizeof(int));

            cudaMalloc(&centerD, (nbodies - 1) * sizeof(Point2D));
            cudaMalloc(&lowerupperD, (nbodies - 1) * sizeof(Point4D));
            cudaMalloc(&gabForLeavesD, nbodies * sizeof(Point4D));

            cudaMalloc(&parentD, nnodes * sizeof(int));

            cudaMalloc(&childD, (nbodies - 1) * 2 * sizeof(int));

            cudaMalloc(&rangeD, nnodes * 2 * sizeof(int));

            cudaMalloc(&levelUnsortD, (nbodies - 1) * sizeof(int));
            cudaMalloc(&levelSortD, (nbodies - 1) * sizeof(int));
            cudaMalloc(&indexUnsortD, (nbodies - 1) * sizeof(int));
            cudaMalloc(&indexSortD, (nbodies - 1) * sizeof(int));
            cudaMalloc(&indexSortTD, (nbodies - 1) * sizeof(int));

            cudaMalloc(&sortObjectsBufferD, nbodies * 64 * sizeof(char));
            sortObjectsBufferSizeInBytes = nbodies * 64;
            cudaMalloc(&sortInternalCellsBufferD, (nbodies - 1) * 64 * sizeof(char));
            sortInternalCellsBufferSizeInBytes = (nbodies - 1) * 64;
            
            switch (objectType)
            {
            case object_T::panel:
                cudaMalloc(&objectD, nbodies * sizeOfElement);
                if(duplicateObj)
                    cudaMalloc(&gabForLeavesD, nbodies * sizeof(double4));
                break;

            case object_T::point2:
            case object_T::point3:
                cudaMalloc(&gabForLeavesD, nbodies * sizeof(double4));
                if (duplicateObj)
                    cudaMalloc(&objectD, nbodies * sizeOfElement);
                break;
            }

            reservedMemorySizeItems = nCudaObject;
        }//реаллокация
        CudaTestError("memoryAllocate launch failed");
    }//memoryAllocate(...)

    void CudaTreeInfo::RadixSortMortonCodes(int beginBit, int endBit)
    {
        size_t temp_storage_bytes = 0;

        cub::DeviceRadixSort::SortPairs(NULL, temp_storage_bytes, \
            mortonCodesKeyUnsortD, mortonCodesKeyD, \
            mortonCodesIdxUnsortD, mortonCodesIdxD, \
            nObject, beginBit, endBit);

        if (temp_storage_bytes > sortObjectsBufferSizeInBytes)
        {
            printf("REALLOC RadixSortMortonCodes: was reserved = %d, new num_items = %d, new reserver = %d (bytes)\n", sortObjectsBufferSizeInBytes, nObject, (int)temp_storage_bytes);

            cudaError_t err1 = cudaFree(sortObjectsBufferD);
            if (err1 != cudaSuccess)
                printf("%s (RadixSortMortonCodes)\n", cudaGetErrorString(err1));

            cudaMalloc(&sortObjectsBufferD, temp_storage_bytes);
            sortObjectsBufferSizeInBytes = (int)temp_storage_bytes;
            printf("REALLOCATION RadixSortMortonCodes is done!\n");
        }

        // Run sorting operation         
        cub::DeviceRadixSort::SortPairs(sortObjectsBufferD, temp_storage_bytes, \
            mortonCodesKeyUnsortD, mortonCodesKeyD, \
            mortonCodesIdxUnsortD, mortonCodesIdxD, \
            nObject, beginBit, endBit);
    }//RadixSortMortonCodes(...)

    void CudaTreeInfo::RadixSortInternalCells(int beginBit, int endBit)
    {
        size_t temp_storage_bytes = 0;

        cub::DeviceRadixSort::SortPairs(NULL, temp_storage_bytes, \
            levelUnsortD, levelSortD, \
            indexUnsortD, indexSortD, \
            nObject - 1, beginBit, endBit);

        if (temp_storage_bytes > sortInternalCellsBufferSizeInBytes)
        {
            printf("REALLOC RadixSortInternalCells: was reserved = %d, new num_items = %d, new reserver = %d (bytes)\n", sortInternalCellsBufferSizeInBytes, nObject - 1, (int)temp_storage_bytes);

            cudaError_t err1 = cudaFree(sortInternalCellsBufferD);
            if (err1 != cudaSuccess)
                printf("%s (RadixSortMortonCodes)\n", cudaGetErrorString(err1));

            cudaMalloc(&sortInternalCellsBufferD, temp_storage_bytes);
            sortInternalCellsBufferSizeInBytes = (int)temp_storage_bytes;
            printf("REALLOCATION RadixSortInternalCells is done!\n");
        }

        // Run sorting operation         
        cub::DeviceRadixSort::SortPairs(sortInternalCellsBufferD, temp_storage_bytes, \
            levelUnsortD, levelSortD, \
            indexUnsortD, indexSortD, \
            nObject - 1, beginBit, endBit);
    }//RadixSortInternalCells(...)

    bool CudaTreeInfo::IsInitialized() const
    {
        return objectD != nullptr;
    }//IsInitialized()

    float CudaTreeInfo::Build()
    {
        float time = 0.0f;
        if (nObject > 0)
        {
            time += treeBoundingBoxWrapper(*this);
            time += treeMortonCodesWrapper(*this, true);

            if (treeType != tree_T::contr)
                time += treeMortonInternalNodesWrapper(*this, true);
        }
        return time;
    }

    float CudaTreeInfo::UpwardTraversal(int order)
    {
        float time = 0.0f;
        if (nObject > 0)
        {
            time += treeClearKernelWrapper(*this, order);
            if (treeType != tree_T::contr && treeType != tree_T::aux)
                time += treeSummarizationWrapper(*this, order, true);
            else if (treeType == tree_T::aux)
                time += treeCalcAABBWrapper(*this);
        }
        return time;
    }

    float CudaTreeInfo::DownwardTraversalVorticesToPoints(CudaTreeInfo& cntrTree, Point2D* velD, double* epsastD, double eps2, double theta, int order, bool calcRadius)
    {
        float time = 0.0f;
        if (nObject > 0)
        {
            double itheta2 = 1.0 / (theta * theta);
            cudaEvent_t start, stop;
            float time;

            cudaEventCreate(&start);  cudaEventCreate(&stop);
            cudaEventRecord(start, 0);

            BHcu::treeVorticesToPointsCalculationKernel<12> << <nBlock * FACTORforces, THREADSforces >> > (
                nNode, nObject, itheta2, eps2, (int2*)childD, (double2*)momsD, \
                (double3*)objectD, mortonCodesIdxD, (double2*)centerD, indexSortD, \
                indexSortTD, (double4*)lowerupperD, \
                cntrTree.nObject, (double3*)cntrTree.objectD, cntrTree.mortonCodesIdxD, (double2*)velD, calcRadius, epsastD);

            cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

            CudaTestError("treeVorticesInfluenceCalculationKernel launch failed");

            cudaEventDestroy(start);  cudaEventDestroy(stop);
        }
        return time;    
    }   

    float CudaTreeInfo::DownwardTraversalVorticesToPanels(CudaTreeInfo& cntrTree, double* rhsD, double* rhsLinD, double theta, int order)
    {
        float time = 0.0f;
        if (nObject > 0)
        {
            double itheta2 = 1.0 / (theta * theta);
            cudaEvent_t start, stop;
            float time;

            cudaEventCreate(&start);  cudaEventCreate(&stop);
            cudaEventRecord(start, 0);

            int nTotPan = cntrTree.nObject;
            int blocks = nBlock;
            int th5rhs = nTotPan / blocks;
            if (th5rhs < order)
                th5rhs = order;
            else
                th5rhs = (int)std::ceil((double)nTotPan / blocks);

            int resultBlock = 1;
            if (th5rhs >= 64)  //Округляем вниз до предыдущей степени двойки, т.е.
            {                                           // (1024...2047) -> 512
                while (resultBlock <= th5rhs)           // ( 512...1023) -> 256
                    resultBlock <<= 1;                  // ( 256...511 ) -> 128
                th5rhs = min(1024, (resultBlock >> 2)); // ( 128...255 ) ->  64
            }                                           // (  64...127 ) ->  32        

            if ((th5rhs > 32) && (th5rhs < 64))
                th5rhs = 32;

            treeRhsCalculationKernel<12> << <blocks * FACTORrhs, th5rhs >> > (nNode, nObject, itheta2, (int2*)childD, \
                (double2*)momsD, (double3*)objectD, mortonCodesIdxD, (double2*)centerD, \
                indexSortD, indexSortTD, (double4*)lowerupperD, \
                cntrTree.nObject, (double4*)cntrTree.gabForLeavesD, (double2*)cntrTree.objectD, cntrTree.mortonCodesIdxD, rhsD, rhsLinD);

            cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);

            CudaTestError("treeRhsCalculationKernel launch failed");

            cudaEventDestroy(start);  cudaEventDestroy(stop);
        }
        return time;
    }

    float CudaTreeInfo::DownwardTraversalPanelsToPoints(CudaTreeInfo& cntrTree, Point2D* velD, double eps2, double theta,
        int order)
    {
        float time = 0.0f;
        if (nObject > 0)
        {
            double itheta2 = 1.0 / (theta * theta);
            time = treePanelsToPointsCalculationWrapper(*this, cntrTree, order, itheta2, eps2, velD);
        }
        return time;
    }

    void CudaTreeInfo::MemoryAllocateForGMRES()
    {
        cudaMalloc(&matVecMulInfo.nClosePanelsD, nObject * sizeof(int));
        cudaMalloc(&matVecMulInfo.nFarCellsD, nObject * sizeof(int));
        cudaMalloc(&matVecMulInfo.closePrefixSumD, (nObject +1)* sizeof(int));
        cudaMalloc(&matVecMulInfo.farPrefixSumD, (nObject + 1) * sizeof(int));
    }

    void CudaTreeInfo::MemoryFreeForGMRES()
    {
        cudaFree(matVecMulInfo.i00D);
        if (schemeType == scheme_T::linScheme)
        {
            cudaFree(matVecMulInfo.i01D);
            cudaFree(matVecMulInfo.i10D);
            cudaFree(matVecMulInfo.i11D);
        }
        cudaFree(matVecMulInfo.nClosePanelsD);
        cudaFree(matVecMulInfo.nFarCellsD);
        cudaFree(matVecMulInfo.closeCellsIdxD);
        cudaFree(matVecMulInfo.farCellsIdxD);
        cudaFree(matVecMulInfo.closePrefixSumD);
        cudaFree(matVecMulInfo.farPrefixSumD);
    }

    float CudaTreeInfo::DownwardTraversalGMRES(double* resD, double* resLinD, double theta, int order, int iter)
    {
        float time = 0.0f;
        double itheta2 = 1.0 / (theta * theta);

        if (iter == 0) 
        {
            time += treeMatrToVecNoCalculationWrapper(*this, itheta2);
            matVecMulInfo.nClosePanels.resize(nObject);
            matVecMulInfo.nFarCells.resize(nObject);

            cudaMemcpy(matVecMulInfo.nClosePanels.data(), matVecMulInfo.nClosePanelsD, nObject * sizeof(int), cudaMemcpyDeviceToHost);
            matVecMulInfo.closePrefixSum.resize(nObject + 1, 0);

            for (size_t i = 0; i < matVecMulInfo.nClosePanels.size(); ++i)
                matVecMulInfo.closePrefixSum[i + 1] = matVecMulInfo.closePrefixSum[i] + matVecMulInfo.nClosePanels[i];
            
            cudaMemcpy(matVecMulInfo.closePrefixSumD, matVecMulInfo.closePrefixSum.data(), (nObject + 1)* sizeof(int), cudaMemcpyHostToDevice);

            cudaMalloc(&matVecMulInfo.closeCellsIdxD, matVecMulInfo.closePrefixSum.back() * sizeof(int));

            cudaMalloc(&matVecMulInfo.i00D, matVecMulInfo.closePrefixSum.back() * sizeof(Point2D));

            if (schemeType == scheme_T::linScheme)
            {
                cudaMalloc(&matVecMulInfo.i01D, matVecMulInfo.closePrefixSum.back() * sizeof(Point2D));
                cudaMalloc(&matVecMulInfo.i10D, matVecMulInfo.closePrefixSum.back() * sizeof(Point2D));
                cudaMalloc(&matVecMulInfo.i11D, matVecMulInfo.closePrefixSum.back() * sizeof(Point2D));
            }

            cudaMemcpy(matVecMulInfo.nFarCells.data(), matVecMulInfo.nFarCellsD, nObject * sizeof(int), cudaMemcpyDeviceToHost);
            matVecMulInfo.farPrefixSum.resize(nObject + 1, 0);

            for (size_t i = 0; i < matVecMulInfo.nFarCells.size(); ++i)
                matVecMulInfo.farPrefixSum[i + 1] = matVecMulInfo.farPrefixSum[i] + matVecMulInfo.nFarCells[i];

            cudaMemcpy(matVecMulInfo.farPrefixSumD, matVecMulInfo.farPrefixSum.data(), (nObject + 1) * sizeof(int), cudaMemcpyHostToDevice);

            cudaMalloc(&matVecMulInfo.farCellsIdxD, matVecMulInfo.farPrefixSum.back() * sizeof(int));
        }//if (iter == 0)

        treeMatrToVecWrapper(*this, order, itheta2, resD, resLinD, iter);

        return time;

    }
    
    /******************************************************************************/
    /*** compute diffusive velo ***************************************************/
    /******************************************************************************/

    float CudaTreeInfo::I1I2CalculationWrapper(double minRd, double* __restrict I1D, Point2D* __restrict I2D, double* __restrict epsastD)
    {
        float time = 0;
        cudaEvent_t start, stop;
        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        treeI1I2CalculationKernel << <nBlock * FACTORI1I2, THREADSI1I2 >> > (
            nNode, nObject, minRd, (int2*)childD, (double2*)momsD,
            (double3*)objectD, mortonCodesIdxD,
            (double2*)centerD, indexSortD, indexSortTD,
            (double4*)lowerupperD,
            (double*)I1D,
            (double2*)I2D,
            epsastD);
        CudaTestError("treeI1I2CalculationKernel launch failed");

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);
        cudaEventDestroy(start);  cudaEventDestroy(stop);

        return time;
    }//I1I2CalculationWrapper(...)


    float CudaTreeInfo::I0I3CalculationWrapper(double minRd, float* __restrict I0D, Point2Df* __restrict I3D, double* __restrict epsastD, const double* __restrict meanEpsD, int nPan, double* __restrict panD, double* __restrict visstrD)
    {
        float time = 0;
        cudaEvent_t start, stop;
        cudaEventCreate(&start);  cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        ClearI0I3 << <nBlock * FACTORforces, THREADSforces >> > (nObject, I0D, (float2*)I3D);
        CudaTestError("ClearI0I3 launch failed");

        treeI0I3CalculationKernel << <nBlock * FACTORI0I3, THREADSI0I3 >> > (
            nNode, nObject, minRd, (int2*)childD, (double2*)momsD,
            (double3*)objectD, mortonCodesIdxD,
            (double2*)centerD, indexSortD, indexSortTD,
            (double4*)lowerupperD, (int*)rangeD,
            I0D, (float2*)I3D, epsastD, meanEpsD, nPan, panD, visstrD);
        CudaTestError("treeI0I3CalculationKernel failed");

        cudaEventRecord(stop, 0);  cudaEventSynchronize(stop);  cudaEventElapsedTime(&time, start, stop);
        cudaEventDestroy(start);  cudaEventDestroy(stop);
        return time;
    }//cuI0I3CalculationKernel2(...)

}