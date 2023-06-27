#pragma once 

#ifndef TYPES_CUH_
#define TYPES_CUH_

#define CALCinDOUBLE

#include "Point2D.h"
#include "Vortex2D.h"

namespace BHcu
{
    typedef double real;
    #define real2 double2
    #define real3 double3
    #define make_real2(x,y) make_double2(x,y)
    #define realmax max
    #define realmin min
    #define realPoint Point2D
    #define realVortex Vortex2D
    #define floatPoint Point2Df
    #define intPair long long


    ////////////////////////////////////////////////////////////////////////////////////////
    struct CUDApointersPoints
    {
        realPoint* maxrl, * minrl;       //габаритный прямоугольник

        int* MmortonCodesKeyUnsortl;
        int* MmortonCodesKeyl;

        int* MmortonCodesIdxUnsortl; //0 1 2 3 ... nbodies-1		
        int* MmortonCodesIdxl;
    };

    struct CUDApointers
    {
        /// ПРОКОММЕНТИРОВАННАЯ ЧАСТЬ
        //без буквы "l" на конце - на host,
        //c буквой "l" на конце - на device	

        int* massl;	//массы (единица для точечного вихря, число вихрей для ячейки)

        realPoint* maxrl, * minrl;       //габаритный прямоугольник
        realPoint* momsl;               //мультипольные моменты всех ячеек; хранятся в виде <mom_0x, mom_0y=0, mom_1x, mom_1y, ..., mom_px, mom_py>, <для второй ячейки> ...
        realPoint* El; //к-ты локальных разложений

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

        realPoint* Mposl;  //Положения внутренних узлов в дерева Карраса
        realPoint* Mlowerl; //Левый нижний угол ячейки
        realPoint* Mupperl; //Правый верхний угол ячейки
        
        int* Mparentl;     //Номер ячейки-родителя
        intPair* Mchildl;  //Потомки внутренних ячеек
        intPair* Mrangel;  //Диапазон частиц во внутренней ячейке
    };
}


#endif