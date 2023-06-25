/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.11   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2022/08/07     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: Gpudefs.h                                                        |
| Info: Source code of VMlib                                                  |
|                                                                             |
| This file is part of VMlib.                                                 |
| VMLib is free software: you can redistribute it and/or modify it            |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VMlib is distributed in the hope that it will be useful, but WITHOUT        |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VMlib.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Описание констант и параметров для взаимодействия с графическим ускорителем
\author Марчевский Илья Константинович
\version 1.11
\date 07 августа 2022 г.
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
#define INC_VORT_DEV (1024)

#define CU_I1I2
#define CU_RHS
#define CU_I0I3
#define CU_CONV_TOWAKE
#define CU_CONV_TOBOU
#define CU_CONVVIRT
#define CU_PAIRS
#define CU_VP


#endif
