/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: types.cuh                                                        |
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
\brief Вспомогательная структура данных
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef TYPES_CUH
#define TYPES_CUH

#include "Point2D.h"
#include "Vortex2D.h"

namespace BHcu
{
    struct CUDApointers
    {

        int* massl;	//массы (единица для точечного вихря, число вихрей для ячейки)

        Point2D* maxrl, * minrl;      //габаритный прямоугольник
        Point2D* momsl;               //мультипольные моменты всех ячеек; хранятся в виде <mom_0x, mom_0y=0, mom_1x, mom_1y, ..., mom_px, mom_py>, <для второй ячейки> ...
        Point2D* El; //к-ты локальных разложений

        //For Morton tree
        int* MmortonCodesKeyUnsortl;
        int* MmortonCodesKeyl;

        int* MmortonCodesIdxUnsortl; //0 1 2 3 ... nbodies-1		
        int* MmortonCodesIdxl;

        int* MlevelUnsortl;
        int* MlevelSortl;

        int* MindexUnsortl;  //0 1 2 3 ... nbodies-2
        int* MindexSortl;
        int* MindexSortTl;

        Point2D* Mposl;  //Положения внутренних узлов в дерева Карраса
        Point4D* Mlowerupperl; //Левый нижний угол ячейки
                
        int* Mparentl;     //Номер ячейки-родителя
        long long* Mchildl;  //Потомки внутренних ячеек (в одной ячейке храним сразу два целых числа)
        long long* Mrangel;  //Диапазон частиц во внутренней ячейке (в одной ячейке храним сразу два целых числа)

        void* sortObjectsBuffer;
        int sortObjectsBufferSizeInBytes;
        
        void* sortInternalCellsBuffer;
        int sortInternalCellsBufferSizeInBytes;
    };
}


#endif