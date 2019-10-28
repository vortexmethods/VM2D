/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.6    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/10/28     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM    |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: cuLib2D.cuh                                                      |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM.                                                    |
| VM is free software: you can redistribute it and/or modify it               |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VM is distributed in the hope that it will be useful, but WITHOUT           |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM.  If not, see <http://www.gnu.org/licenses/>.                 |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Заголовочный файл с описанием функций библиотеки cuLib для работы с CUDA
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.6   
\date 28 октября 2019 г.
*/

#ifndef CUVELOCITYBIOTSAVART_CUH
#define CUVELOCITYBIOTSAVART_CUH

#include "Gpudefs.h"

#if defined(__CUDACC__) || defined(USE_CUDA)

#include "Vortex2D.h"

void cuDevice(int n);

void cuSetConstants(size_t pos_, size_t posR_, size_t posG_);
void cuSetAccelCoeff(double cft_);
void cuSetCollapseCoeff(double pos_, double refLength_);
void cuSetMaxGamma(double gam_);

void cuReserveDevMem(void*& ptr, size_t nBytes);
void cuClearWakeMem(size_t new_n, double* dev_ptr);
void cuCopyWakeToDev(size_t n, const Vortex2D* host_src, double* dev_ptr);
void cuCopyRsToDev(size_t n, const Point2D* host_src, double* dev_ptr);
void cuCopyFixedArray(void* dev_ptr, void* host_src, size_t nBytes);
void cuCopyMemFromDev(void* host_ptr, void* dev_ptr, size_t nBytes);
void cuDeleteFromDev(void* devPtr);

////////////////////////////////////////////////////////////////
void cuCalculateConvVeloWake(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, size_t nsr, double* sr, size_t nAfls, size_t* nVtxs, double** ptrVtxs, double* vel, double* rd, double minRd, double eps2, bool onlyRadius);
void cuCalculateConvVeloWakeFromVirtual(size_t myDisp, size_t myLen, double* pt, size_t npnl, double* r, double* freegamma, double* attgamma, double* attsource, double* vel, double eps2);

void cuCalculateDiffVeloWake(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, double* i1, double* i2, double* rd);
void cuCalculateDiffVeloWakeMesh(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, int* mesh, double meshStep, double* i1, double* i2, double* rd);
void cuCalculateDiffVeloWakeFromPanels(size_t myDisp, size_t myLen, double* pt, size_t npnl, double* r, double* freegamma, double* i1, double* i2, double* rd);

void cuCalculateSurfDiffVeloWake(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, double* i0, double* i3, double* rd);
void cuCalculateRhs(size_t myDisp, size_t myLen, size_t npt, double* pt, size_t nvt, double* vt, size_t nsr, double* sr, double* rhs);

void cuCalculatePairs(size_t myDisp, size_t myLen, size_t npt, double* pt, int* mesh, int* nei, double meshStep, double epsCol2, int type);


//void cuTEST(const std::string& str);

#endif

#endif
