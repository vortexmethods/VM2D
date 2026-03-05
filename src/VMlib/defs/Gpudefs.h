/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Gpudefs.h                                                        |
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
\brief Описание констант и параметров для взаимодействия с графическим ускорителем
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef GPUDEFS_H
#define GPUDEFS_H

/// вихрь Ламба
//#define LAMBVORTEX

/// разделенная механика
#define INITIAL

/// мост
//#define BRIDGE


/// \brief признак использования CUDA
/// 
/// Устанавливается автоматически средствами cmake
/// Если на компьютере найден CUDA Toolkit, то программа будет использовать графическую карту.
/// Чтобы этого не происходило, надо принудительно сделать раскомментировать undef USE_CUDA
// #define USE_CUDA
// #undef USE_CUDA

// далее -- технические "обертки" под CUDA

#define CUDA_ENABLE(a) a
#define CUDA_DISABLE(a)

#if defined(USE_CUDA)
#define CUDAdef CUDA_ENABLE
#else
#define CUDAdef CUDA_DISABLE
#endif

#define IFCUDA_(e, a) e(a)
#define IFCUDA(...) IFCUDA_(CUDAdef, __VA_ARGS__)


#define CUBLOCK (128)
#define INC_VORT_DEV (1024 * 100)

#define CU_I1I2
#define CU_RHS
#define CU_I0I3
#define CU_CONV_TOWAKE
#define CU_CONV_TOBOU
#define CU_CONVVIRT
#define CU_PAIRS
#define CU_VP

#define orderAlignment 16

#define codeLength 14

#define idpid (0.15915494309189533576888376337251)
#define valPi (3.1415926535897932384626433832795)
#define valPif (3.14159265f)

//for tree
#define FACTORgab 6
#define THREADSgab 32

#define FACTORupward 4
#define THREADSupward 32

#define FACTORforces 1
#define THREADSforces 1024

#define FACTORrhs 1
#define maxTHREADSrhs 1024

#define FACTORpanToPoint 1
#define THREADSpanToPoint 1024

#define FACTORnear 1
#define THREADSnear 32

#define FACTORsegintersect 1
#define THREADSsegintersect 32

#define THREADSslae 32
#define FACTORslae 1

#define maxTHREADSslae 1024

#define THREADSI1I2 32
#define THREADSI0I3 32

#define FACTORI1I2 1
#define FACTORI0I3 1




//Характеристики дерева
enum class tree_T
{
    contr, //контрольное
    aux,   //вспомогательное
    vortex, source   //влияющее
};

enum class object_T
{
    point2, point3, panel
};

enum class scheme_T
{
    constScheme, linScheme, noScheme
};

#endif
