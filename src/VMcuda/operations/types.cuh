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
        realPoint* maxrl, * minrl;       //���������� �������������

        int* MmortonCodesKeyUnsortl;
        int* MmortonCodesKeyl;

        int* MmortonCodesIdxUnsortl; //0 1 2 3 ... nbodies-1		
        int* MmortonCodesIdxl;
    };

    struct CUDApointers
    {
        /// ������������������� �����
        //��� ����� "l" �� ����� - �� host,
        //c ������ "l" �� ����� - �� device	

        int* massl;	//����� (������� ��� ��������� �����, ����� ������ ��� ������)

        realPoint* maxrl, * minrl;       //���������� �������������
        realPoint* momsl;               //������������� ������� ���� �����; �������� � ���� <mom_0x, mom_0y=0, mom_1x, mom_1y, ..., mom_px, mom_py>, <��� ������ ������> ...
        realPoint* El; //�-�� ��������� ����������

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

        realPoint* Mposl;  //��������� ���������� ����� � ������ �������
        realPoint* Mlowerl; //����� ������ ���� ������
        realPoint* Mupperl; //������ ������� ���� ������
        
        int* Mparentl;     //����� ������-��������
        intPair* Mchildl;  //������� ���������� �����
        intPair* Mrangel;  //�������� ������ �� ���������� ������
    };
}


#endif