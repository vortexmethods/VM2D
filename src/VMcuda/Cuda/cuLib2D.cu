/*-------------------------------*- VMcuda -*----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VMcuda: VM2D/VM3D Library | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  Open Source Code          *----------------*
|  ####  ##   ##   ##   ##  ##  |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2024 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: cuLib2D.cu                                                       |
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
\brief Файл с реализацией функций библиотеки VMcuda для работы с CUDA
\author Марчевский Илья Константинович
\Version 1.12
\date 14 января 2024 г.
*/

#include <iostream>
#include <algorithm>

#include "cuLib2D.cuh"

#include "cuda.h"
#include "Gpudefs.h"

#include <cusolverDn.h>


__device__ __constant__ size_t sizeVort;
__device__ __constant__ size_t posR;
__device__ __constant__ size_t posG;

__device__ __constant__ double accelCoeff;

__device__ __constant__ double maxGamma;
__device__ __constant__ double collapseRightBorder;
__device__ __constant__ double collapseScale;

__device__ __constant__ double iDPIminEpsAst2;

__device__ __constant__ int schemeSwitcher;
__device__ __constant__ int currentStep;


#define invdpi (0.15915494309189533576888376337251)
#define pi (3.1415926535897932384626433832795)


void cuAlloc(void** ptr, size_t numBytes)
{
	cudaHostAlloc(ptr, numBytes, cudaHostAllocDefault);
}

void cuDalloc(void* ptr)
{
	cudaError_t err1 = cudaFreeHost(ptr);

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (cuDalloc) " << std::endl;
}


__device__ double myAtomicAdd(double* address, double val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;

	do {
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +	__longlong_as_double(assumed)));

		// Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
	} while (assumed != old);
	return __longlong_as_double(old);
}





__device__ inline double myMax(double x, double y)
{
	return (x > y) ? x : y;
}

__device__ inline double myMin(double x, double y)
{
	return (x > y) ? y : x;
}

__device__ inline int sqr(int x)
{
	return x * x;
}

__device__ inline double sqr(double x)
{
	return x * x;
}


/// \brief Способ сглаживания скорости вихря (вихрь Рэнкина или вихрь Ламба)
__device__ inline double CUboundDenom(double r2, double eps2)
{
#ifndef LAMBVORTEX
	return myMax(r2, eps2);
#else
	if (r2 > eps2)
		return r2;
	else
		return (r2 < 1e-10) ? 1e-10 : r2  / (1.0 - exp(-6.0*r2 / eps2));
#endif
}



//Ниже - ядра (__global__), затем - "обычные" функции (__host__)

__global__ void CU_WakeToZero(size_t nvt, double* vt)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;	

	if (i < nvt)
	{
		vt[i * sizeVort + posR + 0] = 0.0;
		vt[i * sizeVort + posR + 1] = 0.0;
		vt[i * sizeVort + posG] = 0.0;
	}
}


__global__ void CU_calc_conv_epsast(
	size_t npt, double* pt,
	size_t nvt, double* vt,
	size_t nsr, double* sr,
	double eps2,
	double* vel, double* rad,
	size_t nAfls, size_t* nVtxs, double** ptrVtxs,
	bool calcVelo, bool calcRadius)
{	
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];
	__shared__ double shg[CUBLOCK];

	size_t i = blockIdx.x * blockDim.x + threadIdx.x;

	double ptx = pt[i*sizeVort + posR + 0];
	double pty = pt[i*sizeVort + posR + 1];

	double velx = 0.0;
	double vely = 0.0;

	double dx, dy, dr2;
	double izn;
		
#ifndef TESTONLYVELO
	double d_1 = 1e+5;
	double d_2 = 1e+5;
	double d_3 = 1e+5;
	double d_0 = 1e+5;
	double dst23, dst12, dst01;
#endif	

	//vortices
	for (size_t j = 0; j < nvt; j += CUBLOCK)
	{		
		shx[threadIdx.x] = vt[(j + threadIdx.x)*sizeVort + posR + 0];
		shy[threadIdx.x] = vt[(j + threadIdx.x)*sizeVort + posR + 1];
		shg[threadIdx.x] = vt[(j + threadIdx.x)*sizeVort + posG + 0];

		__syncthreads();
	
		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			if (j + q < nvt)
			{
				dx = ptx - shx[q];
				dy = pty - shy[q];
				dr2 = dx * dx + dy * dy;

				if (calcVelo)
				{
					izn = shg[q] / myMax(dr2, eps2);// / CUboundDenom(dr2, eps2); //Сглаживать надо!!!
										
					velx -= dy * izn;
					vely += dx * izn;
				}

#ifndef TESTONLYVELO
				if (calcRadius)
				{
					if (d_3 > dr2)
					{
						//if(i == 93)
						//printf("upd: from %d, dr2 = %f\n", (int)(j+q), dr2);

						dst23 = myMin(dr2, d_2);
						d_3 = myMax(dr2, d_2);

						dst12 = myMin(dst23, d_1);
						d_2 = myMax(dst23, d_1);

						dst01 = myMin(dst12, d_0);
						d_1 = myMax(dst12, d_0);
						d_0 = dst01;
					}
				}
#endif
			}
		}
		__syncthreads();
	}
	//if (i == 93)
	//	printf("A: d0,1,2,3: %f, %f, %f, %f\n", d_0, d_1, d_2, d_3);


	//sources
	if (calcVelo)
	{
		for (size_t j = 0; j < nsr; j += CUBLOCK)
		{
			shx[threadIdx.x] = sr[(j + threadIdx.x)*sizeVort + posR + 0];
			shy[threadIdx.x] = sr[(j + threadIdx.x)*sizeVort + posR + 1];
			shg[threadIdx.x] = sr[(j + threadIdx.x)*sizeVort + posG + 0] * accelCoeff;

			__syncthreads();

			for (size_t q = 0; q < CUBLOCK; ++q)
			{
				if (j + q < nsr)
				{
					dx = ptx - shx[q];
					dy = pty - shy[q];
					dr2 = dx * dx + dy * dy;

					izn = shg[q] / CUboundDenom(dr2, eps2); //Сглаживать надо!!!

					velx += dx * izn;
					vely += dy * izn;
				}
			}
			__syncthreads();
		}
	}

	if (calcRadius)
	{
#ifndef TESTONLYVELO
		/*
		size_t nTotVtxs = 0;
		for (size_t p = 0; p < nAfls; ++p)
			nTotVtxs += nVtxs[p];

		//for (size_t p = 0; p < nAfls; ++p)						
		size_t p = 0;
			//for (size_t j = 0; j < nVtxs[p]; j += CUBLOCK)
			for (size_t j = 0; j < nTotVtxs; j += CUBLOCK)
			{
				//if (((&ptrVtxs[p][0] + ((j + threadIdx.x) * sizeVort + posR + 1)) - &ptrVtxs[0][0]) > 4096 * 3)
				//{
				//	//printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %d\n", ((&ptrVtxs[p][0] + ((j + threadIdx.x) * sizeVort + posR + 1)) - &ptrVtxs[0][0]));
				//	printf("p = %d, nVtxs[p] = %d, j = %d\n", (int)p, (int)(nVtxs[p]), (int)j);
				//}

				shx[threadIdx.x] = ptrVtxs[p][(j + threadIdx.x)*sizeVort + posR + 0];
				shy[threadIdx.x] = ptrVtxs[p][(j + threadIdx.x)*sizeVort + posR + 1];

				__syncthreads();

				for (size_t q = 0; q < CUBLOCK; ++q)
				{
					if (j + q < nVtxs[p])
					{
						dx = ptx - shx[q];
						dy = pty - shy[q];
						dr2 = dx * dx + dy * dy;


						if (d_3 > dr2)
						{
							dst23 = myMin(dr2, d_2);
							d_3 = myMax(dr2, d_2);

							dst12 = myMin(dst23, d_1);
							d_2 = myMax(dst23, d_1);

							dst01 = myMin(dst12, d_0);
							d_1 = myMax(dst12, d_0);
							d_0 = dst01;
						}
					}
				}
				__syncthreads();
			}
			if (i == 0)
				printf("B: d0,1,2,3: %f, %f, %f, %f\n", d_0, d_1, d_2, d_3);
		*/
#endif
	}
	
	if (i < npt)
	{
		if (calcVelo)
		{
			vel[2 * i + 0] = velx * invdpi;
			vel[2 * i + 1] = vely * invdpi;
		}

		if (calcRadius)
		{
#ifndef TESTONLYVELO
			rad[i] =  1.0 * sqrt((d_1 + d_2 + d_3) * 0.3333333333333333);

			//printf("rad_%d = %f\n", (int)i, rad[i]);

			//rad[i] =  4.0 * sqrt((d_1 + d_2 + d_3) * 0.3333333333333333);
#endif
		}
	}	
}



__global__ void CU_calc_conv_From_Panels(
	size_t npt, double* pt,
	size_t npnl, double* r, 
	double* freegamma, double* freegammalin,
	double* attgamma, double* attgammalin,
	double* attsource, double* attsourcelin,
	double eps2,
	double* vel)
{
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];

	__shared__ double shdx[CUBLOCK];
	__shared__ double shdy[CUBLOCK];

	__shared__ double shlen[CUBLOCK];


	__shared__ double shfreegamma[CUBLOCK];
	__shared__ double shattgamma[CUBLOCK];
	__shared__ double shattsource[CUBLOCK];

	__shared__ double shfreegammaLin[CUBLOCK];
	__shared__ double shattgammaLin[CUBLOCK];
	__shared__ double shattsourceLin[CUBLOCK];

	size_t i = blockIdx.x * blockDim.x + threadIdx.x;

	double ptx = pt[i*sizeVort + posR + 0];
	double pty = pt[i*sizeVort + posR + 1];

	double velx = 0.0;
	double vely = 0.0;

	double sx, sy, px, py, s2, p2, alpha, lambda, taux, tauy, u1x, u1y, skos0x, skos0y, skos1x, skos1y;

	//if (i == 0)
	//	printf("Point[0]: (%f, %f)\n", ptx, pty);

	for (size_t j = 0; j < npnl; j += CUBLOCK)
	{
		shx[threadIdx.x] = r[(j + threadIdx.x) * 4 + 0];
		shy[threadIdx.x] = r[(j + threadIdx.x) * 4 + 1];
		
		shdx[threadIdx.x] = r[(j + threadIdx.x) * 4 + 2] - shx[threadIdx.x];
		shdy[threadIdx.x] = r[(j + threadIdx.x) * 4 + 3] - shy[threadIdx.x];

		shlen[threadIdx.x] = sqrt(shdx[threadIdx.x] * shdx[threadIdx.x] + shdy[threadIdx.x] * shdy[threadIdx.x]);

		shfreegamma[threadIdx.x] = freegamma[j + threadIdx.x];
		shattgamma[threadIdx.x] = attgamma[j + threadIdx.x];
		shattsource[threadIdx.x] = attsource[j + threadIdx.x];

		if (schemeSwitcher == 2)
		{
			shfreegammaLin[threadIdx.x] = freegammalin[j + threadIdx.x];
			shattgammaLin[threadIdx.x] = attgammalin[j + threadIdx.x];
			shattsourceLin[threadIdx.x] = attsourcelin[j + threadIdx.x];
		}

		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			if (j + q < npnl)
			{
				sx = ptx - shx[q];
				sy = pty - shy[q];

				px = sx - shdx[q];
				py = sy - shdy[q];

				alpha = atan2(px*sy - py*sx, px*sx + py*sy);
				
				s2 = sx * sx + sy * sy;
				p2 = px * px + py * py;

				if ((s2 > 1e-16) && (p2 > 1e-16)) 
					lambda = 0.5*log(s2 / p2);					
				else
					lambda = 0.0;

				taux = shdx[q] / shlen[q];
				tauy = shdy[q] / shlen[q];


				skos0x = alpha * tauy + lambda * taux;
				skos0y = -alpha * taux + lambda * tauy;

				if (schemeSwitcher == 2)
				{
				        u1x = 0.5 / shlen[q] * ((px + sx) * taux * taux \
					     + 2.0 * (py + sy) * taux * tauy - (px + sx) * tauy * tauy);

				        u1y = 0.5 / shlen[q] * (-(py + sy) * taux * taux \
					     + 2.0 * (px + sx) * taux * tauy + (py + sy) * tauy * tauy);

					skos1x = alpha * u1y + lambda * u1x - taux;
					skos1y = -alpha * u1x + lambda * u1y - tauy;
				}

				//if ((i == 0) && (j + q < 10))
				//	printf("j+q = %d, gam_c,gamatt_c,qatt_c = (%f, %f, %f), skos = (%f, %f)\n", (int)(j + q), shfreegamma[q], shattgamma[q], shattsource[q], skos0x, skos0y);

				velx += (shfreegamma[q] + shattgamma[q]) * (-skos0y) + shattsource[q] * skos0x;
				vely += (shfreegamma[q] + shattgamma[q]) * ( skos0x) + shattsource[q] * skos0y;

				if (schemeSwitcher == 2)
				{
					//if ((i == 0) && (j + q < 10))
					//	printf("j+q = %d, gam_l,gamatt_l,qatt_l = (%f, %f, %f), skos = (%f, %f)\n", (int)(j + q), shfreegammaLin[q], shattgammaLin[q], shattsourceLin[q], skos1x, skos1y);

					velx += (shfreegammaLin[q] + shattgammaLin[q]) * (-skos1y) + shattsourceLin[q] * skos1x;
					vely += (shfreegammaLin[q] + shattgammaLin[q]) * ( skos1x) + shattsourceLin[q] * skos1y;
				}
			}
		}
		__syncthreads();
	}

	if (i < npt)
	{
		vel[2 * i + 0] = velx * invdpi;
		vel[2 * i + 1] = vely * invdpi;
	}
}


__global__ void CU_calc_I1I2(
	size_t npt, double* pt,
	size_t nvt, double* vt,
	double* i1, double* i2,
	double* rd, double minRd)
{
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];
	__shared__ double shg[CUBLOCK];
	
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;

	double ptx = pt[i*sizeVort + posR + 0];
	double pty = pt[i*sizeVort + posR + 1];
	double rdi = myMax(rd[i], minRd);

	double val1 = 0.0;
	double val2x = 0.0;
	double val2y = 0.0;

	double dx, dy, dr;
	double expr, exprdivdr;
	
	double diffRadius = 8.0 * rdi;

	double left = ptx - diffRadius;
	double right = ptx + diffRadius;

	for (size_t j = 0; j < nvt; j += CUBLOCK)
	{
		shx[threadIdx.x] = vt[(j + threadIdx.x)*sizeVort + posR + 0];
		shy[threadIdx.x] = vt[(j + threadIdx.x)*sizeVort + posR + 1];
		shg[threadIdx.x] = vt[(j + threadIdx.x)*sizeVort + posG + 0];
		
		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			if (j + q < nvt)
			{
				if ((shx[q] < right) && (shx[q] > left))
				{
					dx = ptx - shx[q];
					dy = pty - shy[q];

					dr = sqrt(dx*dx + dy * dy);

					if ((dr < diffRadius) && (dr > 1e-10))
					{
						expr = shg[q] * exp(-dr / rdi);
						exprdivdr = expr / dr;
						val1 += expr;
						val2x += exprdivdr * dx;
						val2y += exprdivdr * dy;
					}//if (rij>1e-10)
				}
			}
		}
		__syncthreads();
	}

	//printf("thread = %d, ptx = %f, rd[i] = %f\n", (int)locI, ptx, rd[i]);

	if (i < npt)
	{
		i1[i] = val1;
		i2[2 * i + 0] = val2x;
		i2[2 * i + 1] = val2y;
	}

	//if ((currentStep == 1) && (i == 0))
	//	printf("fromWake:  i = %d, val1 = %f, val2.x = %f, val2.y = %f\n", (int)i, val1, val2x, val2y);
}


__global__ void CU_calc_I1I2mesh(
	size_t npt, double* pt,
	size_t nvt, double* vt,
	double* i1, double* i2,
	double* rd, int* dev_ptr_mesh)
{
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];
	__shared__ double shg[CUBLOCK];
	__shared__ int shmshx[CUBLOCK];
	__shared__ int shmshy[CUBLOCK];

	size_t i = blockIdx.x * blockDim.x + threadIdx.x;

	double ptx = pt[i*sizeVort + posR + 0];
	double pty = pt[i*sizeVort + posR + 1];
	double rdi = rd[i];

	double val1 = 0.0;
	double val2x = 0.0;
	double val2y = 0.0;

	double dx, dy, dr;
	double expr, exprdivdr;
	
	double diffRadius = 8.0*rdi;

	int imshx = dev_ptr_mesh[2 * i + 0];
	int imshy = dev_ptr_mesh[2 * i + 1];

	for (size_t j = 0; j < nvt; j += CUBLOCK)
	{
		shx[threadIdx.x] = vt[(j + threadIdx.x)*sizeVort + posR + 0];
		shy[threadIdx.x] = vt[(j + threadIdx.x)*sizeVort + posR + 1];
		shg[threadIdx.x] = vt[(j + threadIdx.x)*sizeVort + posG + 0];

		shmshx[threadIdx.x] = dev_ptr_mesh[(j + threadIdx.x)*2 + 0];
		shmshy[threadIdx.x] = dev_ptr_mesh[(j + threadIdx.x)*2 + 1];

		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			if (j + q < nvt)
			{
				if ((abs(imshx - shmshx[q]) < 25) && (abs(imshy - shmshy[q]) < 25))
				{
					dx = ptx - shx[q];
					dy = pty - shy[q];

					dr = sqrt(dx*dx + dy * dy);

					if ((dr < diffRadius) && (dr > 1e-10))
					{
						expr = shg[q] * exp(-dr / rdi);
						exprdivdr = expr / dr;
						val1 += expr;
						val2x += exprdivdr * dx;
						val2y += exprdivdr * dy;
					}//if (rij>1e-10)
				}
			}
		}
		__syncthreads();
	}

	if (i < npt) 
	{
		i1[i] = val1;
		i2[2 * i + 0] = val2x;
		i2[2 * i + 1] = val2y;
	}
}



__global__ void CU_calc_I1I2FromPanels(
	size_t npt, double* pt,
	size_t npnl, double* r, 
	double* freegamma, double* freegammalin,
	double* i1, double* i2,
	double* rd, double minRd)
{
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];
	__shared__ double shxp1[CUBLOCK];
	__shared__ double shyp1[CUBLOCK];

	__shared__ double shtaux[CUBLOCK];
	__shared__ double shtauy[CUBLOCK];

	__shared__ double shlen[CUBLOCK];
	__shared__ double shptG[CUBLOCK];

	size_t i = blockIdx.x * blockDim.x + threadIdx.x;

	double ptx, pty, rdi;

	if (i < npt)
	{
		ptx = pt[i * sizeVort + posR + 0];
		pty = pt[i * sizeVort + posR + 1];
		rdi = myMax(rd[i], minRd);
	}
	double val1 = 0.0;
	double val2x = 0.0;
	double val2y = 0.0;

	double x0, y0, mn;

	double dx, dy, dr;
	double expr, exprdivdr;

	double diffRadius = 8.0*rdi;
	double left = ptx - diffRadius;
	double right = ptx + diffRadius;

	const int nQuadPt = 3; 

	for (size_t j = 0; j < npnl; j += CUBLOCK)
	{
		if (j + threadIdx.x < npnl)
		{
			shx[threadIdx.x] = r[(j + threadIdx.x) * 4 + 0];
			shy[threadIdx.x] = r[(j + threadIdx.x) * 4 + 1];

			shxp1[threadIdx.x] = r[(j + threadIdx.x) * 4 + 2];
			shyp1[threadIdx.x] = r[(j + threadIdx.x) * 4 + 3];

			shtaux[threadIdx.x] = shxp1[threadIdx.x] - shx[threadIdx.x];
			shtauy[threadIdx.x] = shyp1[threadIdx.x] - shy[threadIdx.x];

			shlen[threadIdx.x] = sqrt(shtaux[threadIdx.x] * shtaux[threadIdx.x] + shtauy[threadIdx.x] * shtauy[threadIdx.x]);

			// todo: Учесть freegammalin
			shptG[threadIdx.x] = freegamma[j + threadIdx.x] * shlen[threadIdx.x] / nQuadPt;
		}

		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			if (j + q < npnl)
			{
				for (int s = 0; s < nQuadPt; ++s)
					{
					mn = (s + 0.5) / nQuadPt;
					x0 = shx[q] + shtaux[q] * mn;
					if ((x0 < right) && (x0 > left))
                    {
						y0 = shy[q] + shtauy[q] * mn;


						dx = ptx - x0;
						dy = pty - y0;

						dr = sqrt(dx*dx + dy * dy);

						if ((dr < diffRadius) && (dr > 1e-10))
						{
							expr = shptG[q] * exp(-dr / rdi);
							exprdivdr = expr / dr;
							val1 += expr;
							val2x += exprdivdr * dx;
							val2y += exprdivdr * dy;
						}//if (rij>1e-10)
					}
				}
			}
		}
		__syncthreads();
	}

	if (i < npt)
	{
		i1[i] = val1;
		i2[2 * i + 0] = val2x;
		i2[2 * i + 1] = val2y;
	}

	//if ((currentStep == 1) && (i == 0))
	//	printf("fromSheet: i = %d, val1 = %f, val2.x = %f, val2.y = %f\n", (int)i, val1, val2x, val2y);
}


__global__ void  CU_calc_I0I3(
	size_t npt, double* pt, //вихри
	size_t nvt, double* vt, //панели
	double* i0, double* i3,
	double* rd, double* meanEps,
	double minRd,
	double* visstr)
{
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;	

	double ptx = pt[i*sizeVort + posR + 0];
	double pty = pt[i*sizeVort + posR + 1];
	double ptg = pt[i*sizeVort + posG + 0];
	
	
	double rdi = myMax(rd[i], minRd);   // commented 30-05
	//double rdi = rd[i];

	double val0 = 0.0;
	double val3x = 0.0;
	double val3y = 0.0;

	double iDDomRad = 1.0 / rdi;

	double qx, qy, d/*, q*/;
	double begx, begy, endx, endy;
	double lenj, lenj_m;
	double taux, tauy;
	double s;
	double normx, normy;
	double v0x, v0y;
	double hx, hy;
	double xix, xiy, lxi;
	double expon;
	double mnx, mny;
	int new_n;
	double den;
	double xi_mx, xi_my, lxi_m;
	double mnog1;
	double vs;
	double meanepsj2;

	if (i < npt) //Здесь так делать можно, т.к. shared memory не используется
	{
		for (size_t j = 0; j < nvt; ++j)
		{
			begx = vt[j * 4 + 0];
			begy = vt[j * 4 + 1];
			endx = vt[j * 4 + 2];
			endy = vt[j * 4 + 3];
			
			vs = 0.0;

			qx = ptx - 0.5 * (begx + endx);
			qy = pty - 0.5 * (begy + endy);
			//q = sqrt(qx * qx + qy * qy);

			lenj = sqrt( (endx - begx)*(endx - begx) + (endy - begy)*(endy - begy) );

			taux = (endx - begx) / lenj;
			tauy = (endy - begy) / lenj;

			s = qx * taux + qy * tauy;

			normx = tauy;
			normy = -taux;
					   
			d = fabs(qx*normx + qy * normy);
			
			meanepsj2 = sqr(myMax(meanEps[j], minRd));    // commented 30-05
			//meanepsj2 = sqr(meanEps[j]);

			if ( (d < 50.0 * lenj) && (fabs(s) < 50.0 * lenj) )	//Почему зависит от длины панели???
			{
				v0x = taux * lenj;
				v0y = tauy * lenj;
								
				if ( (d > 5.0 * lenj) || (fabs(s) > 5.0 * lenj) )
				{
					xix = qx * iDDomRad;
					xiy = qy * iDDomRad;
					lxi = sqrt(xix * xix + xiy * xiy);

					expon = exp(-lxi) * lenj;
					mnx = normx * expon;
					mny = normy * expon;

					if (val0 != -pi * rdi)
					{
						val0 += (xix * mnx + xiy * mny) * (lxi + 1.0) / (lxi * lxi);
						val3x += mnx;
						val3y += mny;
					}

					vs = ptg * expon / (pi * meanepsj2);

				}					
				//else if ( (d >= 0 * lenj) || (fabs(s) >= 0 * lenj) )
				else if ( (d >= 0.01 * lenj) || (fabs(s) > 0.45  * lenj) )
				{
					//new_n = 100;
					//new_n = (int)(ceil(5.0 * lenj / q));

					den = (fabs(s) < 0.5 * lenj) ? d : (fabs(s) + d - 0.5 * lenj);
					new_n = (int)myMin((int)myMax(ceil(10.0 * lenj / den), 1.0), 20);

					hx = v0x / new_n;
					hy = v0y / new_n;

					for (int m = 0; m < new_n; ++m)
					{	
						//double midx = ptx - d * normx;
						//double midy = pty - d * normy;
						//double midintx = begx + hx * (m + 0.5);
						//double midinty = begx + hy * (m + 0.5);
						//int dopnew_n = 1;

						//if(((midintx - midx)* (midintx - midx) + (midinty - midy) * (midinty - midy)) < 100*d*d)
						//	dopnew_n = (int)myMin((int)myMax(ceil(1000.0 * lenj / den), 1.0), 2000);
						

						//for (int mm = 0; mm < dopnew_n; ++mm)
						//{
							xi_mx = (ptx - (begx + hx * (m + 0.5))) * iDDomRad;
							xi_my = (pty - (begy + hy * (m + 0.5))) * iDDomRad;

							lxi_m = sqrt(xi_mx * xi_mx + xi_my * xi_my);

							lenj_m = lenj / new_n;
							expon = exp(-lxi_m) * lenj_m;

							mnx = normx * expon;
							mny = normy * expon;

							if (val0 != -pi * rdi)
							{
								val0 += (xi_mx * mnx + xi_my * mny) * (lxi_m + 1.0) / (lxi_m * lxi_m);
								val3x += mnx;
								val3y += mny;
							}

							vs += expon;
						//}//for mm
					}//for m
					vs *= ptg / (pi * meanepsj2);

					/*
					if ((d <= 0.001 * lenj) && (fabs(s) < 0.45 * lenj))
					{
						double Rval0 = -pi * rdi;

						double Rmnog1 = 2.0 * rdi * (1.0 - exp(-lenj * 0.5 * iDDomRad) * cosh(fabs(s) * iDDomRad));
						double Rval3x = Rmnog1 * normx;
						double Rval3y = Rmnog1 * normy;
						double Rvs = Rmnog1 * ptg / (pi * meanepsj2);

						printf("v0: %f <-> %f, v3: (%f, %f) <-> (%f, %f), vs: %f <-> %f\n", val0, Rval0, val3x, val3y, Rval3x, Rval3y, vs, Rvs);
					}
					//*/

				}  				
				else
				{
					val0 = -pi * rdi;
										
					mnog1 = 2.0 * rdi * (1.0 - exp(-lenj * 0.5 * iDDomRad) * cosh(fabs(s) * iDDomRad));
					val3x = mnog1 * normx;
					val3y = mnog1 * normy;
					vs = mnog1 * ptg / (pi * meanepsj2);

				}
				
			}//if d<50 len 
			//if (i == 12)
			//	printf("j = %d, val3: %f, %f\n", (int)j, val3x, val3y);
#ifdef  __CUDA_ARCH__
#if __CUDA_ARCH__ <= 600
			myAtomicAdd(visstr + j, vs);
			//printf("MY: __CUDA_ARCH__ = %d\n", __CUDA_ARCH__);
#else
			atomicAdd(visstr + j, vs);
			//printf("CUDA: __CUDA_ARCH__ = %d\n", __CUDA_ARCH__);
#endif
#endif //  __CUDA_ARCH__			

		}//for j


	}

	if (i < npt)
	{
		i0[i] = val0;
		i3[2 * i + 0] = val3x;
		i3[2 * i + 1] = val3y;
	}
}


__global__ void CU_calc_RHS(	
	size_t npt, double* pt,
	size_t nvt, double* vt,
	size_t nsr, double* sr,
	double eps2,
	double* rhs,
	double* rhsLin)
{
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];
	__shared__ double shg[CUBLOCK];

	size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	   
	double begx = 0.0, begy = 0.0, endx = 0.0, endy = 0.0, dlen = 0.0, dix = 0.0, diy = 0.0, taux = 0.0, tauy = 0.0, u1x = 0.0, u1y = 0.0;
	if (i < npt)
	{
		begx = pt[i * 4 + 0];
		begy = pt[i * 4 + 1];
		endx = pt[i * 4 + 2];
		endy = pt[i * 4 + 3];

		dix = endx - begx;
		diy = endy - begy;

		dlen = sqrt(dix * dix + diy * diy);

		taux = dix / dlen;
		tauy = diy / dlen;
	}

	double val = 0.0;
	double valLin = 0.0;

	double sx, sy, px, py;
	double alpha, lambda, tempVel, tempVelLin; //из двух к-тов alpha и lambda в принципе можно для экономии пользоваться одной и той же переменной

	//vortices
	for (size_t j = 0; j < nvt; j += CUBLOCK)
	{
		shx[threadIdx.x] = vt[(j + threadIdx.x) * sizeVort + posR + 0];
		shy[threadIdx.x] = vt[(j + threadIdx.x) * sizeVort + posR + 1];
		shg[threadIdx.x] = vt[(j + threadIdx.x) * sizeVort + posG + 0];

		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			if (j + q < nvt)
			{
				if ((schemeSwitcher == 1) || (schemeSwitcher == 2))
				{
					sx = shx[q] - begx;
					sy = shy[q] - begy;

					px = shx[q] - endx;
					py = shy[q] - endy;

					alpha = atan2(px * sy - py * sx, px * sx + py * sy);

					tempVel = shg[q] * alpha;
					val -= tempVel;
					//printf("pnl = %d, val = %f\n", locI, val);
				}

				if (schemeSwitcher == 0)
					val += 0.5*(shg[q] * dlen / myMax(dlen*dlen, eps2)) * (taux * (sx+px) + tauy * (sy+py));					
								
				if (schemeSwitcher == 2)
				{
					u1x = 0.5 / dlen * ((px + sx) * taux * taux \
						+ 2.0 * (py + sy) * taux * tauy - (px + sx) * tauy * tauy);
					u1y = 0.5 / dlen * (-(py + sy) * taux * taux \
						+ 2.0 * (px + sx) * taux * tauy + (py + sy) * tauy * tauy);

					lambda = 0.5 * log((sx * sx + sy * sy) / (px * px + py * py));

					tempVelLin = shg[q] * (alpha * (u1x * taux + u1y * tauy) + lambda * (-u1y * taux + u1x * tauy));

					valLin -= tempVelLin;
				}
			}
		}
		__syncthreads();
	}

	//sources	
	for (size_t j = 0; j < nsr; j += CUBLOCK)
	{
		shx[threadIdx.x] = sr[(j + threadIdx.x) * sizeVort + posR + 0];
		shy[threadIdx.x] = sr[(j + threadIdx.x) * sizeVort + posR + 1];
		shg[threadIdx.x] = sr[(j + threadIdx.x) * sizeVort + posG + 0] * accelCoeff;

		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			if (j + q < nsr)
			{
				if ((schemeSwitcher == 1) || (schemeSwitcher == 2))
				{
					sx = shx[q] - begx;
					sy = shy[q] - begy;

					px = shx[q] - endx;
					py = shy[q] - endy;

					lambda = 0.5 * log((sx * sx + sy * sy) / (px * px + py * py));

					tempVel = shg[q] * lambda;
					val -= tempVel;
				}

				if (schemeSwitcher == 0)
				{
					val += 0.5*(shg[q] * dlen / myMax(dlen*dlen, eps2)) * (-tauy * (sx+px) + taux * (sy+py));
				}

				if (schemeSwitcher == 2)
				{
					u1x = 0.5 / dlen * ((px + sx) * taux * taux \
						+ 2.0 * (py + sy) * taux * tauy - (px + sx) * tauy * tauy);
					u1y = 0.5 / dlen * (-(px + sx) * taux * taux \
						+ 2.0 * (px + sx) * taux * tauy + (py + sy) * tauy * tauy);

					alpha = atan2(px * sy - py * sx, px * sx + py * sy);

					tempVelLin = shg[q] * (alpha * (u1y * taux - u1x * tauy) + lambda * (u1x * taux + u1y * tauy) - 1.0);

					valLin -= tempVelLin;
				}
			}
		}
		__syncthreads();
	}


	if (i < npt)
	{
		val *= invdpi / dlen;
		rhs[i] = val;

		if (schemeSwitcher == 2)
		{
			valLin *= invdpi / dlen;
			rhsLin[i] = valLin;
		}
	}
}


__global__ void CU_calc_mesh(
	size_t npt, double* pt,
	int* dev_ptr_mesh,
	double meshStep)
{
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;

	if (i < npt)
	{
		double ptx = pt[i * sizeVort + posR + 0];
		double pty = pt[i * sizeVort + posR + 1];

		dev_ptr_mesh[2 * i + 0] = floor(ptx / meshStep);
		dev_ptr_mesh[2 * i + 1] = floor(pty / meshStep);
	}
}


__global__ void CU_calc_nei(
	size_t npt, double* pt,
	int* dev_ptr_mesh, int* dev_ptr_nei,
	double epsCol2, int type)
{
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;
	size_t minI = blockIdx.x * blockDim.x;

	if (i < npt) // можно, т.к. не используем shared memory
	{
		int ix = dev_ptr_mesh[2 * i + 0];
		int iy = dev_ptr_mesh[2 * i + 1];

		int jx, jy;

		double ipx, ipy, jpx, jpy;
		ipx = pt[i*sizeVort + posR + 0];
		ipy = pt[i*sizeVort + posR + 1];

		double dx, dy, r2, r2test;
		dev_ptr_nei[i] = 0;

		double ig = pt[i*sizeVort + posG + 0];
		double jg;

				
		double cftmax = myMax(1.0, /* 2.0 * */ (ipx-collapseRightBorder) / collapseScale);
		//double cftmax = myMax(1.0, (ipx - 0.5) / 0.1);
		
		double cftmax2 = cftmax * cftmax;

		r2test = (type == 1) ? 4.0*epsCol2 * cftmax2 : epsCol2 * cftmax2;
		
		//if (cftmax > 1)
		//	r2test = sqr(0.005*collapseScale);

		int fracMesh = (int)(r2test / epsCol2);
   
		bool cond;

		for (size_t j = minI; j < npt; ++j)
		{
			jx = dev_ptr_mesh[2 * j + 0];
			jy = dev_ptr_mesh[2 * j + 1];

			if ((sqr(abs(ix - jx)) <= fracMesh) && (sqr(abs(iy - jy)) <= fracMesh) && (j > i))
			{
				jpx = pt[j*sizeVort + posR + 0];
				jpy = pt[j*sizeVort + posR + 1];

				dx = ipx - jpx;
				dy = ipy - jpy;

				r2 = dx*dx + dy*dy;
				jg = pt[j*sizeVort + posG + 0];

				cond = (r2 < r2test) && ((type == 1) ? ig*jg < 0 : (   (type == 0) ? ig*jg !=0 : ig*jg > 0  ) && (fabs(ig + jg) < cftmax2 * maxGamma) );
                if (cond)
				{
					dev_ptr_nei[i] = j;
					break;
				}
			}
		}
	}
}

void cuDevice(int n)
{
	cudaSetDevice(n);
}


//ниже - обычные (__host__) функции

int cuCalcBlocks(size_t new_n)
{
	size_t nBlocks = new_n / CUBLOCK;
	if (new_n % CUBLOCK)
		nBlocks++;
	return max((int)nBlocks, 1);
}

void cuSetConstants(size_t pos_, size_t posR_, size_t posG_, int code)
{
	cudaError_t err1 = cudaMemcpyToSymbol(sizeVort, &pos_,  sizeof(size_t));
	cudaError_t err2 = cudaMemcpyToSymbol(posR,     &posR_, sizeof(size_t));
	cudaError_t err3 = cudaMemcpyToSymbol(posG,     &posG_, sizeof(size_t));

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erSetConst01, code = " << code << ")" << std::endl;
	if (err2 != cudaSuccess)
		std::cout << cudaGetErrorString(err2) << " (erSetConst02, code =" << code << ")" << std::endl;
	if (err3 != cudaSuccess)
		std::cout << cudaGetErrorString(err3) << " (erSetConst03, code =" << code << ")" << std::endl;
}

void cuSetAccelCoeff(double cft_, int code)
{
	cudaError_t err1 = cudaMemcpyToSymbol(accelCoeff, &cft_, sizeof(double));

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erSetAccelCoeff01, code =" << code << ")" << std::endl;
}


void cuSetCollapseCoeff(double pos_, double refLength_, int code)
{
	cudaError_t err1 = cudaMemcpyToSymbol(collapseRightBorder, &pos_, sizeof(double));

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erSetCollapseCoeff01, code =" << code << ")" << std::endl;

	cudaError_t err2 = cudaMemcpyToSymbol(collapseScale, &refLength_, sizeof(double));

	if (err2 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erSetCollapseCoeff02, code =" << code << ")" << std::endl;
}


void cuSetMaxGamma(double gam_, int code)
{
	cudaError_t err1 = cudaMemcpyToSymbol(maxGamma, &gam_, sizeof(double));

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erSetMaxGamma01, code =" << code << ")" << std::endl;
}

void cuSetSchemeSwitcher(int schemeSwitcher_, int code)
{
	cudaError_t err1 = cudaMemcpyToSymbol(schemeSwitcher, &schemeSwitcher_, sizeof(int));

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erSchemeSwitcher01, code =" << code << ")" << std::endl;
}

void cuSetCurrentStep(int currentStep_, int code)
{
	cudaError_t err1 = cudaMemcpyToSymbol(currentStep, &currentStep_, sizeof(int));

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCurrentStep01, code =" << code << ")" << std::endl;
}



void cuReserveDevMem(void*& ptr, size_t nBytes, int code)
{
	cudaError_t err1 = cudaMalloc(&ptr, nBytes);
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erReserveDevMem01, code =" << code << ")" << std::endl;

}


void cuCheckError(int code)
{
	cudaDeviceSynchronize();
	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (cuCheckError) " << code << std::endl;
}


void cuClearWakeMem(size_t new_n, double* dev_ptr, int code)
{
	dim3 blocks(cuCalcBlocks(new_n)), threads(CUBLOCK);

	CU_WakeToZero<<<blocks, threads>>> (new_n, dev_ptr);

	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_WakeToZero01)" << std::endl;
}

void cuCopyWakeToDev(size_t n, const Vortex2D* host_src, double* dev_ptr, int code)
{
	size_t sizeOfVortex = sizeof(Vortex2D);
	cudaError_t err1 = cudaMemcpy(dev_ptr, host_src, sizeOfVortex * n, cudaMemcpyHostToDevice);
	if (err1 != cudaSuccess)
	{
		std::cout << cudaGetErrorString(err1) << " (erCopyWakeToDev01, code =" << code << ")" << std::endl;
	}

}

void cuCopyWakeToDevAsync(size_t n, const Vortex2D* host_src, double* dev_ptr, int code)
{
	size_t sizeOfVortex = sizeof(Vortex2D);
	cudaError_t err1 = cudaMemcpyAsync(dev_ptr, host_src, sizeOfVortex * n, cudaMemcpyHostToDevice);
	if (err1 != cudaSuccess)
	{
		std::cout << cudaGetErrorString(err1) << " (erCopyWakeToDevAsync01, code =" << code << ")" << std::endl;
	}

}


void cuCopyFixedArray(void* dev_ptr, void* host_src, size_t nBytes, int code)
{
	cudaError_t err1 = cudaMemcpy(dev_ptr, host_src, nBytes, cudaMemcpyHostToDevice);
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCopyFixedArray01, code =" << code << ")" << std::endl;
}

void cuCopyFixedArrayPoint2D(double* dev_ptr, const Point2D* host_src, size_t npts, int code)
{
	cudaError_t err1 = cudaMemcpy(dev_ptr, host_src, sizeof(double) * 2 * npts, cudaMemcpyHostToDevice);

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCopyFixedArrayPoint2D01, code = " << code << ")" << std::endl;
}

void cuCopyFixedArrayPoint4D(double* dev_ptr, const Point2D* host_src, size_t npts, int code)
{
	cudaError_t err1 = cudaMemcpy(dev_ptr, host_src, sizeof(double) * 4 * npts, cudaMemcpyHostToDevice);

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCopyFixedArrayPoint4D01, code = " << code << ")" << std::endl;
}

void cuCopyMemFromDev(void* host_ptr, void* dev_ptr, size_t nBytes, int code)
{
	cudaError_t err1 = cudaMemcpy(host_ptr, dev_ptr, nBytes, cudaMemcpyDeviceToHost);
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCopyMemFromDev01, code =" << code << ")" << std::endl;
}

void cuDeleteFromDev(void* devPtr, int code)
{
	cudaError_t err1 = cudaFree(devPtr);

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erDeleteFromDev01, code =" << code << ")" << std::endl;
}

/////////////////////////////////////////////////////////////
void cuCalculateConvVeloWake(size_t npt, double* pt, size_t nvt, double* vt, size_t nsr, double* sr, size_t nAfls, size_t* nVtxs, double** ptrVtxs, double* vel, double* rd, double eps2, bool calcVelo, bool calcRadius)
{	
	dim3 blocks(cuCalcBlocks(npt)), threads(CUBLOCK);


	/*
	cudaEvent_t start, stop;
	float gpu_time = 0.0f;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0); 
	*/

	CU_calc_conv_epsast <<< blocks, threads >>> (npt, pt, nvt, vt, nsr, sr, eps2, vel, rd, nAfls, nVtxs, ptrVtxs, calcVelo, calcRadius);

	/*
	cudaThreadSynchronize();
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&gpu_time, start, stop);
	printf("\nTime spent: %.5f\n", gpu_time/1000.0f);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	*/


	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_conv_epsast01)" << std::endl;

	cudaDeviceSynchronize();
}


void cuCalculateConvVeloWakeFromVirtual(size_t npt, double* pt, size_t npnl, double* r, double* freegamma, double* freegammalin, double* attgamma, double* attgammalin, double* attsource, double* attsourcelin, double* vel, double eps2)
{
	dim3 blocks(cuCalcBlocks(npt)), threads(CUBLOCK);
	CU_calc_conv_From_Panels<<<blocks, threads>>>(npt, pt, npnl, r, freegamma, freegammalin, attgamma, attgammalin, attsource, attsourcelin, eps2, vel);

	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_conv01)" << std::endl;

	cudaDeviceSynchronize();
}


void cuCalculateDiffVeloWake(size_t npt, double* pt, size_t nvt, double* vt, double* i1, double* i2, double* rd, double minRad)
{
	dim3 blocks(cuCalcBlocks(npt)), threads(CUBLOCK);
	CU_calc_I1I2 << < blocks, threads >> > (npt, pt, nvt, vt, i1, i2, rd, minRad);

	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_I1I201)" << std::endl;
}

void cuCalculateDiffVeloWakeMesh(size_t npt, double* pt, size_t nvt, double* vt, int* mesh, double meshStep, double* i1, double* i2, double* rd)
{
	dim3 blocks1(cuCalcBlocks(nvt)), blocks2(cuCalcBlocks(npt)), threads(CUBLOCK);
	CU_calc_mesh << < blocks1, threads >> > (nvt, vt, mesh, meshStep);
	
	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_mesh01)" << std::endl;

	CU_calc_I1I2mesh << < blocks2, threads >> > (npt, pt, nvt, vt, i1, i2, rd, mesh);

	cudaError_t err2 = cudaGetLastError();
	if (err2 != cudaSuccess)
		std::cout << cudaGetErrorString(err2) << " (erCU_calc_I1I2mesh01)" << std::endl;
}

void cuCalculateDiffVeloWakeFromPanels(size_t npt, double* pt, size_t npnl, double* r, double* freegamma, double* freegammalin, double* i1, double* i2, double* rd, double minRad)
{
	dim3 blocks(cuCalcBlocks(npt)), threads(CUBLOCK);
	CU_calc_I1I2FromPanels << < blocks, threads >> > (npt, pt, npnl, r, freegamma, freegammalin, i1, i2, rd, minRad);

	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_I1I201)" << std::endl;
}


void cuCalculateSurfDiffVeloWake(size_t npt, double* pt, size_t nvt, double* vt, double* i0, double* i3, double* rd, double* meanEps, double minRd, double* visstr)
{
	dim3 blocks(cuCalcBlocks(npt)), threads(CUBLOCK);
	CU_calc_I0I3<<<blocks, threads>>>(npt, pt, nvt, vt, i0, i3, rd, meanEps, minRd, visstr);
	
	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_I0I31)" << std::endl;
}


void cuCalculateRhs(size_t npt, double* pt, size_t nvt, double* vt, size_t nsr, double* sr, double eps2, double* rhs, double* rhsLin)
{
	dim3 blocks(cuCalcBlocks(npt)), threads(CUBLOCK);
	
	CU_calc_RHS << < blocks, threads >> > (npt, pt, nvt, vt, nsr, sr, eps2, rhs, rhsLin);
	
	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_RHS01)" << std::endl;
}



void cuCalculatePairs(size_t npt, double* pt, int* mesh, int* nei, double meshStep, double epsCol2, int type)
{
	dim3 blocksMesh(cuCalcBlocks(npt)), threadsMesh(CUBLOCK);
	CU_calc_mesh << < blocksMesh, threadsMesh >> > (npt, pt, mesh, meshStep);

	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_mesh01)" << std::endl;
	
	dim3 blocksNei(cuCalcBlocks(npt)), threadsNei(CUBLOCK);
	CU_calc_nei << < blocksNei, threadsNei >> > (npt, pt, mesh, nei, epsCol2, type);
	
	cudaError_t err2 = cudaGetLastError();
	if (err2 != cudaSuccess)
		std::cout << cudaGetErrorString(err2) << " (erCU_calc_nei01)" << std::endl;
}


/*
void cuTEST(const std::string& str)
{
	double* ppp;
	double mmm = 134;

	cudaMalloc(&ppp, 8);

	cudaError_t err1 = cudaMemcpy(ppp, &mmm, sizeof(double), cudaMemcpyHostToDevice);
	if (err1 != cudaSuccess)
		std::cout << str << ": " << cudaGetErrorString(err1) << " (CUDA_TEST_BREAK)" << std::endl;
	else
		std::cout << str << ": " << cudaGetErrorString(err1) << " (CUDA_TEST_PASSED)" << std::endl;
}
*/


void cuInverseMatrix(int n, double* matrPtr, double* invMatrPtr)
{
	cusolverDnHandle_t cusolverH = NULL;
	cudaStream_t stream = NULL;


	double *devMatr, *devInvMatr; // , * devSol;
	cudaMalloc(&devMatr, n * n * sizeof(double));
	cudaMalloc(&devInvMatr, n * n * sizeof(double));
	cudaError_t cudaStat1 = cudaMemcpy(devMatr, matrPtr, sizeof(double) * n * n, cudaMemcpyHostToDevice);
	if (cudaStat1 != cudaSuccess)
		std::cout << cudaGetErrorString(cudaStat1) << " (cuInverseMatrix01)" << std::endl;

	cudaError_t cudaStat2 = cudaMemcpy(devInvMatr, invMatrPtr, sizeof(double) * n * n, cudaMemcpyHostToDevice);
	if (cudaStat2 != cudaSuccess)
		std::cout << cudaGetErrorString(cudaStat2) << " (cuInverseMatrix02)" << std::endl;

	int* d_Ipiv = NULL; /* pivoting sequence */
	int* d_info = NULL; /* error info */

	int  lwork = 0;     /* size of workspace */
	double* d_work = NULL; /* device workspace for getrf */

	cudaError_t cudaStat3 = cudaMalloc(&d_Ipiv, sizeof(int) * n);
	if (cudaStat3 != cudaSuccess)
		std::cout << cudaGetErrorString(cudaStat3) << " (cuInverseMatrix03)" << std::endl;

	cudaError_t cudaStat4 = cudaMalloc(&d_info, sizeof(int));
	if (cudaStat4 != cudaSuccess)
		std::cout << cudaGetErrorString(cudaStat4) << " (cuInverseMatrix04)" << std::endl;

	/* step 1: create cusolver handle, bind a stream */
	cusolverStatus_t status = cusolverDnCreate(&cusolverH);
	if (CUSOLVER_STATUS_SUCCESS != status)
		std::cout << " (cuBlas01)" << std::endl;

	cudaError_t cudaStat5 = cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
	if (cudaStat5 != cudaSuccess)
		std::cout << cudaGetErrorString(cudaStat5) << " (cuInverseMatrix05)" << std::endl;

	status = cusolverDnSetStream(cusolverH, stream);
	if (CUSOLVER_STATUS_SUCCESS != status)
		std::cout << " (cuBlas02)" << std::endl;

	bool pivot_on = true;

	/* step 3: query working space of getrf */
	status = cusolverDnDgetrf_bufferSize(
		cusolverH,
		n,
		n,
		devMatr,
		n,
		&lwork);
	if (CUSOLVER_STATUS_SUCCESS != status)
		std::cout << " (cuBlas03)" << std::endl;

	cudaError_t cudaStat6 = cudaMalloc((void**)&d_work, sizeof(double) * lwork);
	if (cudaStat6 != cudaSuccess)
		std::cout << cudaGetErrorString(cudaStat6) << " (cuInverseMatrix06)" << std::endl;

	/* step 4: LU factorization */
	if (pivot_on) {
		status = cusolverDnDgetrf(
			cusolverH,
			n,
			n,
			devMatr,
			n,
			d_work,
			d_Ipiv,
			d_info);
	}
	else {
		status = cusolverDnDgetrf(
			cusolverH,
			n,
			n,
			devMatr,
			n,
			d_work,
			NULL,
			d_info);
	}
	cudaError_t cudaStat7 = cudaDeviceSynchronize();

	if (CUSOLVER_STATUS_SUCCESS != status)
		std::cout << " (cuBlas04)" << std::endl;
	
	if (cudaStat7 != cudaSuccess)
		std::cout << cudaGetErrorString(cudaStat7) << " (cuInverseMatrix07)" << std::endl;


	/* step 5: solve */
	if (pivot_on) {
		status = cusolverDnDgetrs(
			cusolverH,
			CUBLAS_OP_N,
			n,
			n, /* nrhs */
			devMatr,
			n,
			d_Ipiv,
			devInvMatr,
			n,
			d_info);
	}
	else {
		status = cusolverDnDgetrs(
			cusolverH,
			CUBLAS_OP_N,
			n,
			1, /* nrhs */
			devMatr,
			n,
			NULL,
			devInvMatr,
			n,
			d_info);
	}

	cudaError_t cudaStat8 = cudaDeviceSynchronize();

	if (CUSOLVER_STATUS_SUCCESS != status)
		std::cout << " (cuBlas05)" << std::endl;

	if (cudaStat8 != cudaSuccess)
		std::cout << cudaGetErrorString(cudaStat8) << " (cuInverseMatrix08)" << std::endl;

	cudaError_t cudaStat9 = cudaMemcpy(invMatrPtr, devInvMatr, sizeof(double) * n * n, cudaMemcpyDeviceToHost);
	if (cudaStat9 != cudaSuccess)
		std::cout << cudaGetErrorString(cudaStat9) << " (cuInverseMatrix09)" << std::endl;

	cusolverDnDestroy(cusolverH);
		
	cudaError_t cudafree1 = cudaFree(devMatr);
	if (cudafree1 != cudaSuccess)
		std::cout << cudaGetErrorString(cudafree1) << " (cuInverseMatrixFree1)" << std::endl;
	
	cudaError_t cudafree2 = cudaFree(devInvMatr);
	if (cudafree2 != cudaSuccess)
		std::cout << cudaGetErrorString(cudafree2) << " (cuInverseMatrixFree2)" << std::endl;
	
	cudaError_t cudafree3 = cudaFree(d_Ipiv);
	if (cudafree3 != cudaSuccess)
		std::cout << cudaGetErrorString(cudafree3) << " (cuInverseMatrixFree3)" << std::endl;
	
	cudaError_t cudafree4 = cudaFree(d_info);
	if (cudafree4 != cudaSuccess)
		std::cout << cudaGetErrorString(cudafree4) << " (cuInverseMatrixFree4)" << std::endl;
	
	cudaError_t cudafree5 = cudaFree(d_work);
	if (cudafree5 != cudaSuccess)
		std::cout << cudaGetErrorString(cudafree5) << " (cuInverseMatrixFree5)" << std::endl;
}