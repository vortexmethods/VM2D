/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.3    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/09/26     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: gpudefs.h                                                        |
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
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/

#ifndef GPUDEFS_H
#define GPUDEFS_H

//#define USE_CUDA

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


#endif
