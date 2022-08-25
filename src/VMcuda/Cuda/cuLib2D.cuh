/*-------------------------------*- VMcuda -*----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.11   |
| ##  ## ### ### ##  ## ##  ##  |  VMcuda: VM2D/VM3D Library | 2022/08/07     |
| ##  ## ## # ##    ##  ##  ##  |  Open Source Code          *----------------*
|  ####  ##   ##   ##   ##  ##  |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: cuLib2D.cuh                                                      |
| Info: Source code of VMcuda                                                 |
|                                                                             |
| This file is part of VMcuda.                                                |
| VMcuda is free software: you can redistribute it and/or modify it           |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VMcuda is distributed in the hope that it will be useful, but WITHOUT       |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VMcuda.  If not, see <http://www.gnu.org/licenses/>.             |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Заголовочный файл с описанием функций библиотеки VMcuda для работы с CUDA
\author Марчевский Илья Константинович
\version 1.11
\date 07 августа 2022 г.
*/

#ifndef CUVELOCITYBIOTSAVART_CUH
#define CUVELOCITYBIOTSAVART_CUH

#include "Gpudefs.h"

#if defined(__CUDACC__) || defined(USE_CUDA)

#include "Vortex2D.h"

void cuDevice(int n);

void cuSetConstants(size_t pos_, size_t posR_, size_t posG_, int code = 0);
void cuSetAccelCoeff(double cft_, int code = 0);
void cuSetCollapseCoeff(double pos_, double refLength_, int code = 0);
void cuSetMaxGamma(double gam_, int code = 0);
void cuSetSchemeSwitcher(int schemeSwitcher_, int code);

void cuReserveDevMem(void*& ptr, size_t nBytes, int code = 0);
void cuClearWakeMem(size_t new_n, double* dev_ptr, int code = 0);
void cuCopyWakeToDev(size_t n, const Vortex2D* host_src, double* dev_ptr, int code = 0);
void cuCopyWakeToDevAsync(size_t n, const Vortex2D* host_src, double* dev_ptr, int code = 0);

void cuCopyFixedArray(void* dev_ptr, void* host_src, size_t nBytes, int code = 0);
void cuCopyFixedArrayPoint2D(double* dev_ptr, const Point2D* host_src, size_t npts, int code = 0);
void cuCopyFixedArrayPoint4D(double* dev_ptr, const Point2D* host_src, size_t npts, int code = 0);

void cuCopyMemFromDev(void* host_ptr, void* dev_ptr, size_t nBytes, int code = 0);
void cuDeleteFromDev(void* devPtr, int code = 0);

////////////////////////////////////////////////////////////////
void cuCalculateConvVeloWake(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, size_t nsr, double* sr, size_t nAfls, size_t* nVtxs, double** ptrVtxs, double* vel, double* rd, double eps2, bool calcVelo, bool calcRadius);
void cuCalculateConvVeloWakeFromVirtual(size_t myDisp, size_t myLen, double* pt, size_t npnl, double* r, double* freegamma, double* attgamma, double* attsource, double* vel, double eps2);

void cuCalculateDiffVeloWake(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, double* i1, double* i2, double* rd, double minRad);
void cuCalculateDiffVeloWakeMesh(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, int* mesh, double meshStep, double* i1, double* i2, double* rd);
void cuCalculateDiffVeloWakeFromPanels(size_t myDisp, size_t myLen, double* pt, size_t npnl, double* r, double* freegamma, double* i1, double* i2, double* rd, double minRad);

void cuCalculateSurfDiffVeloWake(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, double* i0, double* i3, double* rd, double* meanEps, double minRd, double* visstr);
void cuCalculateRhs(size_t myDisp, size_t myLen, size_t npt, double* pt, size_t nvt, double* vt, size_t nsr, double* sr, double eps2, double* rhs, double* rhsLin);
void cuCalculatePairs(size_t myDisp, size_t myLen, size_t npt, double* pt, int* mesh, int* nei, double meshStep, double epsCol2, int type);


//void cuTEST(const std::string& str);

void cuAlloc(void** ptr, size_t numBytes);
void cuDalloc(void* ptr);


#endif

#endif
