/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.2    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/06/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: cuLib.cu                                                         |
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
\brief Файл с реализацией функций библиотеки cuLib для работы с CUDA
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.2
\date 14 июня 2018 г.
*/

#include "cuLib.cuh"

//#include "cuda.h"
#include "gpudefs.h"


__device__ __constant__ size_t pos;
__device__ __constant__ size_t posR;
__device__ __constant__ size_t posG;

#define invdpi (0.15915494309189533576888376337251)
#define pi (3.1415926535897932384626433832795)


//Ниже - ядра (__global__), затем - "обычные" функции (__host__)

__global__ void CU_WakeToZero(size_t nvt, double* vt)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;

	vt[i*pos + posR + 0] = 0.0;
	vt[i*pos + posR + 1] = 0.0;
	vt[i*pos + posG + 0] = 0.0;
}


__global__ void CU_calc_conv_epsast(
	size_t disp, size_t len, double* pt,
	size_t nvt, double* vt,
	double eps2, double minRd,
	double* vel, double* rad,
	size_t nAfls, size_t* nPnls, double** ptrPnls)
{
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];
	__shared__ double shg[CUBLOCK];

	size_t locI = blockIdx.x * blockDim.x + threadIdx.x;
	size_t i = disp + locI;

	double ptx = pt[i*pos + posR + 0];
	double pty = pt[i*pos + posR + 1];

	double velx = 0.0;
	double vely = 0.0;

	double dx, dy, dr2;
	double izn;

	double d_1 = 1e+5;
	double d_2 = 1e+5;
	double d_3 = 1e+5;
	double d_0 = 1e+5;
	double dst23, dst12, dst01;

	for (size_t j = 0; j < nvt; j += CUBLOCK)
	{
		shx[threadIdx.x] = vt[(j + threadIdx.x)*pos + posR + 0];
		shy[threadIdx.x] = vt[(j + threadIdx.x)*pos + posR + 1];
		shg[threadIdx.x] = vt[(j + threadIdx.x)*pos + posG + 0];

		__syncthreads();
	
		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			dx = ptx - shx[q];
			dy = pty - shy[q];
			dr2 = dx*dx + dy*dy;
			izn = shg[q] / fmax(dr2, eps2);

			velx -= dy * izn;
			vely += dx * izn;

			if (d_3 > dr2) 
			{
				dst23 = fmin(dr2, d_2);
				d_3 = fmax(dr2, d_2);
				
				dst12 = fmin(dst23, d_1);
				d_2 = fmax(dst23, d_1);
				
				dst01 = fmin(dst12, d_0);
				d_1 = fmax(dst12, d_0);
				d_0 = dst01;
			}
		}
		__syncthreads();
	}

	for (size_t q = 0; q < nAfls; ++q)
	for (size_t j = 0; j < nPnls[q]; j += CUBLOCK)
	{
		shx[threadIdx.x] = ptrPnls[q][(j + threadIdx.x)*pos + posR + 0];
		shy[threadIdx.x] = ptrPnls[q][(j + threadIdx.x)*pos + posR + 1];
				
		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			dx = ptx - shx[q];
			dy = pty - shy[q];
			dr2 = dx*dx + dy*dy;

			if (d_3 > dr2) 
			{
				dst23 = fmin(dr2, d_2);
				d_3 = fmax(dr2, d_2);

				dst12 = fmin(dst23, d_1);
				d_2 = fmax(dst23, d_1);

				dst01 = fmin(dst12, d_0);
				d_1 = fmax(dst12, d_0);
				d_0 = dst01;
			}
		}
		__syncthreads();
	}

	vel[2 * locI + 0] = velx * invdpi;
	vel[2 * locI + 1] = vely * invdpi;
	rad[locI] = fmax( sqrt( (d_1 + d_2 + d_3) / 3.0), minRd);
}


__global__ void CU_calc_conv(
	size_t disp, size_t len, double* pt,
	size_t nvt, double* vt,
	double eps2, 
	double* vel)
{
/// \todo Не сделано влияние слоев, добавлен только учет виртуальных вихрей - корректно работает только для неподвижного профиля
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];
	__shared__ double shg[CUBLOCK];

	size_t locI = blockIdx.x * blockDim.x + threadIdx.x;
	size_t i = disp + locI;

	double ptx = pt[i*pos + posR + 0];
	double pty = pt[i*pos + posR + 1];

	double velx = 0.0;
	double vely = 0.0;

	double dx, dy, dr2;
	double izn;

	for (size_t j = 0; j < nvt; j += CUBLOCK)
	{
		shx[threadIdx.x] = vt[(j + threadIdx.x)*pos + posR + 0];
		shy[threadIdx.x] = vt[(j + threadIdx.x)*pos + posR + 1];
		shg[threadIdx.x] = vt[(j + threadIdx.x)*pos + posG + 0];

		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			dx = ptx - shx[q];
			dy = pty - shy[q];
			dr2 = dx*dx + dy*dy;
			izn = shg[q] / fmax(dr2, eps2);

			velx -= dy * izn;
			vely += dx * izn;		
		}
		__syncthreads();
	}	

	vel[2 * locI + 0] = velx * invdpi;
	vel[2 * locI + 1] = vely * invdpi;
}


__global__ void CU_calc_I1I2(
	size_t disp, size_t len, double* pt,
	size_t nvt, double* vt,
	double* i1, double* i2,
	double* rd)
{
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];
	__shared__ double shg[CUBLOCK];
	
	size_t locI = blockIdx.x * blockDim.x + threadIdx.x;
	size_t i = disp + locI;

	double ptx = pt[i*pos + posR + 0];
	double pty = pt[i*pos + posR + 1];
	double rdi = rd[i];

	double val1 = 0.0;
	double val2x = 0.0;
	double val2y = 0.0;

	double dx, dy, dr;
	double expr, exprdivdr;
	
	double diffRadius = 7.0*rdi;

	double left = ptx - diffRadius;
	double right = ptx + diffRadius;

	for (size_t j = 0; j < nvt; j += CUBLOCK)
	{
		shx[threadIdx.x] = vt[(j + threadIdx.x)*pos + posR + 0];
		shy[threadIdx.x] = vt[(j + threadIdx.x)*pos + posR + 1];
		shg[threadIdx.x] = vt[(j + threadIdx.x)*pos + posG + 0];
		
		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			if ((shx[q] < right) && (shx[q] > left))
			{
				dx = ptx - shx[q];
				dy = pty - shy[q];

				dr = sqrt(dx*dx + dy*dy);

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
		__syncthreads();
	}
	i1[locI] = val1;
	i2[2 * locI + 0] = val2x;
	i2[2 * locI + 1] = val2y;
}



__global__ void CU_calc_I1I2mesh(
	size_t disp, size_t len, double* pt,
	size_t nvt, double* vt,
	double* i1, double* i2,
	double* rd, int* dev_ptr_mesh)
{
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];
	__shared__ double shg[CUBLOCK];
	__shared__ int shmshx[CUBLOCK];
	__shared__ int shmshy[CUBLOCK];

	size_t locI = blockIdx.x * blockDim.x + threadIdx.x;
	size_t i = disp + locI;

	double ptx = pt[i*pos + posR + 0];
	double pty = pt[i*pos + posR + 1];
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
		shx[threadIdx.x] = vt[(j + threadIdx.x)*pos + posR + 0];
		shy[threadIdx.x] = vt[(j + threadIdx.x)*pos + posR + 1];
		shg[threadIdx.x] = vt[(j + threadIdx.x)*pos + posG + 0];
		shmshx[threadIdx.x] = dev_ptr_mesh[(j + threadIdx.x)*2 + 0];
		shmshy[threadIdx.x] = dev_ptr_mesh[(j + threadIdx.x)*2 + 1];

		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			if ((abs(imshx-shmshx[q])<25) && (abs(imshy - shmshy[q])<25))
			{
				dx = ptx - shx[q];
				dy = pty - shy[q];

				dr = sqrt(dx*dx + dy*dy);

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
		__syncthreads();
	}
	i1[locI] = val1;
	i2[2 * locI + 0] = val2x;
	i2[2 * locI + 1] = val2y;
}





__global__ void CU_calc_I0I3(
	size_t disp, size_t len, double* pt,
	size_t nvt, double* vt,
	double* i0, double* i3,
	double* rd)
{
	size_t locI = blockIdx.x * blockDim.x + threadIdx.x;
	size_t i = disp + locI;

	double ptx = pt[i*pos + posR + 0];
	double pty = pt[i*pos + posR + 1];
	double rdi = rd[i];

	double val0 = 0.0;
	double val3x = 0.0;
	double val3y = 0.0;

	double iDDomRad = 1.0 / rdi;

	double qx, qy, d;
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
	double xi_mx, xi_my, lxi_m;
	double mnog1;

	for (size_t j = 0; j < nvt - 1; ++j)
	{
		begx = vt[j * 2 + 0];
		begy = vt[j * 2 + 1];
		endx = vt[j * 2 + 2];
		endy = vt[j * 2 + 3];

		qx = ptx - 0.5 * (begx + endx);
		qy = pty - 0.5 * (begy + endy);

		lenj = sqrt((endx - begx)*(endx - begx) + (endy - begy)*(endy - begy));

		taux = (endx - begx) / lenj;
		tauy = (endy - begy) / lenj;

		s = qx * taux + qy * tauy;

		normx = tauy;
		normy = -taux;

		d = sqrt(qx*qx + qy*qy);

		if (d < 50.0 * lenj)	//Почему зависит от длины панели???
		{
			v0x = taux * lenj;
			v0y = tauy * lenj;

			if (d > 5.0 * lenj)
			{
				xix = qx * iDDomRad;
				xiy = qy * iDDomRad;
				lxi = sqrt(xix*xix + xiy*xiy);

				expon = exp(-lxi)*lenj;
				mnx = normx*expon;
				mny = normy*expon;

				if (val0 != -pi * rdi)
				{
					val0 += (xix * mnx + xiy * mny) * (lxi + 1.0) / (lxi*lxi);
					val3x += mnx;
					val3y += mny;
				}

				//viscousStress[j] += locPoints[i].g() * expon * iDPIepscol2;
			}
			else if ((d <= 5.0 * lenj) && (d >= 0.1 * lenj))
			{
				//vs = 0.0;
				//new_n = 100;
				new_n = (int)(ceil(5.0 * lenj / d));

				hx = v0x / new_n;
				hy = v0y / new_n;

				for (int m = 0; m < new_n; m++)
				{
					xi_mx = (ptx - (begx + hx * (m + 0.5))) * iDDomRad;
					xi_my = (pty - (begy + hy * (m + 0.5))) * iDDomRad;

					lxi_m = sqrt(xi_mx*xi_mx + xi_my*xi_my);

					lenj_m = lenj / new_n;
					expon = exp(-lxi_m)*lenj_m;
					

					mnx = normx * expon;
					mny = normy * expon;

					if (val0 != -pi * rdi)
					{
						val0 += (xi_mx*mnx + xi_my*mny) * (lxi_m + 1.0) / (lxi_m*lxi_m);
						val3x += mnx;
						val3y += mny;
					}
					
				}//for m
			}
			else if (d <= 0.1 * lenj)
			{
				val0 = -pi * rdi;
				
				if (fabs(s) > 0.5 * lenj)
				{
					mnog1 = 2.0 * rdi * (exp(-fabs(s)  * iDDomRad) * sinh(lenj * iDDomRad / 2.0));
					val3x = mnog1 * normx;
					val3y = mnog1 * normy;
				}
				else
				{
					mnog1 = 2.0 * rdi * (1.0 - exp(-lenj * iDDomRad / 2.0)*cosh(fabs(s) * iDDomRad));
					val3x = mnog1 * normx;
					val3y = mnog1 * normy;			
				}
				break;

			}
		}//if d<50 len 

	}//for j

	i0[locI] = val0;
	i3[2 * locI + 0] = val3x;
	i3[2 * locI + 1] = val3y;
}







__global__ void CU_calc_RHS(
	size_t disp, size_t len, double* pt,
	size_t nvt, double* vt,
	double* rhs)
{
	__shared__ double shx[CUBLOCK];
	__shared__ double shy[CUBLOCK];
	__shared__ double shg[CUBLOCK];
	
	size_t locI = blockIdx.x * blockDim.x + threadIdx.x;
	size_t i = disp + locI;

	double begx = pt[i*2 + 0];
	double begy = pt[i*2 + 1];
	double endx = pt[i*2 + 2];
	double endy = pt[i*2 + 3];

	double dlen = sqrt((endx-begx)*(endx-begx) + (endy-begy)*(endy-begy));

	double val = 0.0;

	double sx, sy, px, py;
	double alpha, tempVel;

	for (size_t j = 0; j < nvt; j += CUBLOCK)
	{
		shx[threadIdx.x] = vt[(j + threadIdx.x)*pos + posR + 0];
		shy[threadIdx.x] = vt[(j + threadIdx.x)*pos + posR + 1];
		shg[threadIdx.x] = vt[(j + threadIdx.x)*pos + posG + 0];
		
		__syncthreads();

		for (size_t q = 0; q < CUBLOCK; ++q)
		{
			sx = shx[q] - begx;
			sy = shy[q] - begy;

			px = shx[q] - endx;
			py = shy[q] - endy;

			alpha = atan2(px*sy - py*sx, px*sx + py*sy);
		
			tempVel = shg[q] * alpha;
			val -= tempVel;
		}
		__syncthreads();
	}

	val *= invdpi / dlen;
	rhs[locI] = val;
}


__global__ void CU_calc_mesh(
	size_t npt, double* pt,
	int* dev_ptr_mesh,
	double meshStep)
{
	size_t i = blockIdx.x * blockDim.x + threadIdx.x;

	double ptx = pt[i*pos + posR + 0];
	double pty = pt[i*pos + posR + 1];

	dev_ptr_mesh[2 * i + 0] = floor(ptx / meshStep);
	dev_ptr_mesh[2 * i + 1] = floor(pty / meshStep);
}


__global__ void CU_calc_nei(
	size_t disp, size_t len, size_t npt, double* pt,
	int* dev_ptr_mesh, int* dev_ptr_nei,
	double epsCol2, int type)
{
	size_t locI = blockIdx.x * blockDim.x + threadIdx.x;
	size_t i = disp + locI;

	int ix = dev_ptr_mesh[2 * i + 0];
	int iy = dev_ptr_mesh[2 * i + 1];

	int jx, jy;

	double ipx, ipy, jpx, jpy;
	ipx = pt[i*pos + posR + 0];
	ipy = pt[i*pos + posR + 1];

	double dx, dy, r2, r2test;
	dev_ptr_nei[locI] = 0;

	double ig = pt[i*pos + posG + 0];
	double jg;

	r2test = (type == 1) ? 4.0*epsCol2 * fmax(1.0, ipx*ipx) : epsCol2 * fmax(1.0, ipx*ipx);
	//r2test = epsCol2;

	bool cond;

	for (size_t j = 0; j < npt; ++j)
	{
		jx = dev_ptr_mesh[2 * j + 0];
		jy = dev_ptr_mesh[2 * j + 1];
		
		if ((abs(ix-jx) <= 1) && (abs(iy - jy) <= 1) && (j>i))
		{
			jpx = pt[j*pos + posR + 0];
			jpy = pt[j*pos + posR + 1];

			dx = ipx - jpx;
			dy = ipy - jpy;
			
			r2 = dx*dx + dy*dy;
			
			jg = pt[j*pos + posG + 0];

			cond = (r2 < r2test) && ((type == 1) ? ig*jg < 0 : (ig*jg > 0) && (fabs(ig + jg) < 0.03*0.25));
			//cond = (r2 < r2test);
			if (cond)
			{
				dev_ptr_nei[locI] = j;
				break;
			}
		}
	}	
}

//ниже - обычные (__host__) функции

int cuCalcBlocks(size_t new_n)
{
	size_t nBlocks = new_n / CUBLOCK;
	if (new_n%CUBLOCK)
		nBlocks++;
	return max((int)nBlocks, 1);
}

void cuSetConstants(size_t pos_, size_t posR_, size_t posG_)
{
	cudaError_t err1 = cudaMemcpyToSymbol(pos, &pos_, sizeof(size_t));
	cudaError_t err2 = cudaMemcpyToSymbol(posR, &posR_, sizeof(size_t));
	cudaError_t err3 = cudaMemcpyToSymbol(posG, &posG_, sizeof(size_t));

	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erSetConst01)" << std::endl;
	if (err2 != cudaSuccess)
		std::cout << cudaGetErrorString(err2) << " (erSetConst02)" << std::endl;
	if (err3 != cudaSuccess)
		std::cout << cudaGetErrorString(err3) << " (erSetConst03)" << std::endl;
}

void cuReserveDevMem(void*& ptr, size_t nBytes)
{
	cudaError_t err1 = cudaMalloc(&ptr, nBytes);
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erReserveDevMem01)" << std::endl;

}

void cuClearWakeMem(size_t new_n, double* dev_ptr)
{
	dim3 blocks(cuCalcBlocks(new_n)), threads(CUBLOCK);
	CU_WakeToZero << <blocks, threads >> > (new_n, dev_ptr);
}

void cuCopyWakeToDev(size_t n, const Vortex2D* host_src, double* dev_ptr)
{
	size_t sizeOfVortex = sizeof(Vortex2D);
	cudaError_t err1 = cudaMemcpy(dev_ptr, host_src, sizeOfVortex * n, cudaMemcpyHostToDevice);
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCopyWakeToDev01)" << std::endl;

}

void cuCopyRsToDev(size_t n, const Point2D* host_src, double* dev_ptr)
{
	cudaError_t err1 = cudaMemcpy(dev_ptr, host_src, sizeof(double)*2*n, cudaMemcpyHostToDevice);
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCopyRsToDev01)" << std::endl;

}

void cuCopyFixedArray(void* dev_ptr, void* host_src, size_t nBytes)
{
	cudaError_t err1 = cudaMemcpy(dev_ptr, host_src, nBytes, cudaMemcpyHostToDevice);
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCopyFixedArray01)" << std::endl;
}

void cuCopyMemFromDev(void* host_ptr, void* dev_ptr, size_t nBytes)
{
	cudaError_t err1 = cudaMemcpy(host_ptr, dev_ptr, nBytes, cudaMemcpyDeviceToHost);
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCopyMemFromDev01)" << std::endl;
}

void cuDeleteFromDev(void* devPtr)
{
	cudaError_t err1 = cudaFree(devPtr);
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erDeleteFromDev01)" << std::endl;
}

/////////////////////////////////////////////////////////////
void cuCalculateConvVeloWake(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, size_t nAfls, size_t* nPnls, double** ptrPnls, double* vel, double* rd, double minRd, double eps2)
{	
	dim3 blocks(cuCalcBlocks(myLen)), threads(CUBLOCK);
	CU_calc_conv_epsast << < blocks, threads >> > (myDisp, myLen, pt, nvt, vt, eps2, minRd, vel, rd, nAfls, nPnls, ptrPnls);
	
	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_conv_epsast01)" << std::endl;
}

void cuCalculateConvVeloWakeFromVirtual(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, double* vel, double eps2)
{
	dim3 blocks(cuCalcBlocks(myLen)), threads(CUBLOCK);
	CU_calc_conv << < blocks, threads >> > (myDisp, myLen, pt, nvt, vt, eps2, vel);

	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_conv01)" << std::endl;
}


void cuCalculateDiffVeloWake(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, double* i1, double* i2, double* rd)
{
	dim3 blocks(cuCalcBlocks(myLen)), threads(CUBLOCK);
	CU_calc_I1I2 << < blocks, threads >> > (myDisp, myLen, pt, nvt, vt, i1, i2, rd);

	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_I1I201)" << std::endl;
}

void cuCalculateDiffVeloWakeMesh(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, int* mesh, double meshStep, double* i1, double* i2, double* rd)
{
	dim3 blocks1(cuCalcBlocks(nvt)), blocks2(cuCalcBlocks(myLen)), threads(CUBLOCK);
	CU_calc_mesh << < blocks1, threads >> > (nvt, vt, mesh, meshStep);
	
	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_mesh01)" << std::endl;

	CU_calc_I1I2mesh << < blocks2, threads >> > (myDisp, myLen, pt, nvt, vt, i1, i2, rd, mesh);

	cudaError_t err2 = cudaGetLastError();
	if (err2 != cudaSuccess)
		std::cout << cudaGetErrorString(err2) << " (erCU_calc_I1I2mesh01)" << std::endl;
}


void cuCalculateSurfDiffVeloWake(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, double* i0, double* i3, double* rd)
{
	dim3 blocks(cuCalcBlocks(myLen)), threads(CUBLOCK);
	CU_calc_I0I3 << < blocks, threads >> > (myDisp, myLen, pt, nvt, vt, i0, i3, rd);

	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_I0I31)" << std::endl;
}






void cuCalculateRhs(size_t myDisp, size_t myLen, double* pt, size_t nvt, double* vt, double* rhs)
{
	dim3 blocks(cuCalcBlocks(myLen)), threads(CUBLOCK);
	CU_calc_RHS << < blocks, threads >> > (myDisp, myLen, pt, nvt, vt, rhs);

	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_RHS01)" << std::endl;
}



void cuCalculatePairs(size_t myDisp, size_t myLen, size_t npt, double* pt, int* mesh, int* nei, double meshStep, double epsCol2, int type)
{
	dim3 blocks(cuCalcBlocks(npt)), threads(CUBLOCK);
	CU_calc_mesh << < blocks, threads >> > (npt, pt, mesh, meshStep);

	cudaError_t err1 = cudaGetLastError();
	if (err1 != cudaSuccess)
		std::cout << cudaGetErrorString(err1) << " (erCU_calc_mesh01)" << std::endl;

	CU_calc_nei << < blocks, threads >> > (myDisp, myLen, npt, pt, mesh, nei, epsCol2, type);

	cudaError_t err2 = cudaGetLastError();
	if (err2 != cudaSuccess)
		std::cout << cudaGetErrorString(err2) << " (erCU_calc_nei01)" << std::endl;
}