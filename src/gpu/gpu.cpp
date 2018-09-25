/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: gpu.cpp                                                          |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.	                                      |
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
\brief Файл кода с описанием класса gpu
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/

#include "defs.h"
#include "gpu.h"
#include "cuLib.cuh"

#include <omp.h>


gpu::gpu(const std::vector<std::unique_ptr<Boundary>>& boundary_, const Wake& wake_, const Parallel& parallel_)
	: boundary(boundary_), wake(wake_), parallel(parallel_)
{
#if defined(__CUDACC__) || defined(USE_CUDA)
	
	
	cuSetConstants(sizeof(Vortex2D)/sizeof(double), Vortex2D::offsPos / sizeof(double), Vortex2D::offsGam / sizeof(double) );
	
	
	wake.devWakePtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(INC_VORT_DEV, wake.devNWake);
	
	dev_ptr_vel = ReserveDevMem<double, 2>(INC_VORT_DEV, n_CUDA_vel);
	dev_ptr_i0 = ReserveDevMem<double, 1>(INC_VORT_DEV, n_CUDA_i0);
	dev_ptr_i1 = ReserveDevMem<double, 1>(INC_VORT_DEV, n_CUDA_i1);
	dev_ptr_i2 = ReserveDevMem<double, 2>(INC_VORT_DEV, n_CUDA_i2);
	dev_ptr_i3 = ReserveDevMem<double, 2>(INC_VORT_DEV, n_CUDA_i3);
	dev_ptr_rad = ReserveDevMem<double, 1>(INC_VORT_DEV, n_CUDA_rad);

	dev_ptr_mesh = ReserveDevMem<int, 2>(INC_VORT_DEV, n_CUDA_mesh);
	dev_ptr_nei = ReserveDevMem<int, 1>(INC_VORT_DEV, n_CUDA_nei);

	vels.resize(n_CUDA_vel, { 0.0, 0.0 });
	rads.resize(n_CUDA_vel, 1.0e+5);

	i0.resize(n_CUDA_vel, 0.0);
	i1.resize(n_CUDA_vel, 0.0);
	i2.resize(n_CUDA_vel, { 0.0, 0.0 });
	i3.resize(n_CUDA_vel, { 0.0, 0.0 });
	nei.resize(n_CUDA_vel, 0);
	//mesh.resize(n_CUDA_vel, { 0, 0 });

	virtvels.resize(0);
	virtrads.resize(0);
	virti0.resize(0);
	virti1.resize(0);
	virti2.resize(0);
	virti3.resize(0);
	virtrhs.resize(0);

	wake.devNWake = 0;
	n_CUDA_afls = 0;
#endif
}


gpu::~gpu()
{
#if defined(__CUDACC__) || defined(USE_CUDA)
	cuDeleteFromDev(wake.devWakePtr);
	cuDeleteFromDev(dev_ptr_vel);
	cuDeleteFromDev(dev_ptr_rad);
	cuDeleteFromDev(dev_ptr_i0);
	cuDeleteFromDev(dev_ptr_i1);
	cuDeleteFromDev(dev_ptr_i2);
	cuDeleteFromDev(dev_ptr_i3);
	cuDeleteFromDev(dev_ptr_mesh);
	cuDeleteFromDev(dev_ptr_nei);

	for (size_t s = 0; s < n_CUDA_afls; ++s)
	{
		cuDeleteFromDev(host_ptr_ptr_pnl[s]);
		cuDeleteFromDev(host_ptr_ptr_r[s]);
		cuDeleteFromDev(host_ptr_ptr_i0[s]);
		cuDeleteFromDev(host_ptr_ptr_i1[s]);
		cuDeleteFromDev(host_ptr_ptr_i2[s]);
		cuDeleteFromDev(host_ptr_ptr_i3[s]);
		cuDeleteFromDev(host_ptr_ptr_rad[s]);
		cuDeleteFromDev(host_ptr_ptr_vel[s]);
		cuDeleteFromDev(host_ptr_ptr_rhs[s]);
	}

	if (n_CUDA_afls)
	{
		cuDeleteFromDev(dev_nPanels);
		cuDeleteFromDev(dev_ptr_ptr_pnl);
		cuDeleteFromDev(dev_ptr_ptr_r);
		cuDeleteFromDev(dev_ptr_ptr_i0);
		cuDeleteFromDev(dev_ptr_ptr_i1);
		cuDeleteFromDev(dev_ptr_ptr_i2);
		cuDeleteFromDev(dev_ptr_ptr_i3);
		cuDeleteFromDev(dev_ptr_ptr_rad);
		cuDeleteFromDev(dev_ptr_ptr_vel);
		cuDeleteFromDev(dev_ptr_ptr_rhs);
	}
#endif
}

#if defined(__CUDACC__) || defined(USE_CUDA)

void gpu::RefreshWake()
{
	const int& id = parallel.myidWork;

	if (wake.vtx.size() > 0)
	{
		//if (id == 0)
		{
			//size_t cuSize = cuGetCuSize();
			//size_t cpuSize = sizeof(Vortex2D);

			while (wake.vtx.size() > wake.devNWake)
			{
				size_t curLength = wake.devNWake;

				cuDeleteFromDev(wake.devWakePtr);
				cuDeleteFromDev(dev_ptr_vel);
				cuDeleteFromDev(dev_ptr_rad);
				cuDeleteFromDev(dev_ptr_i0);
				cuDeleteFromDev(dev_ptr_i1);
				cuDeleteFromDev(dev_ptr_i2);
				cuDeleteFromDev(dev_ptr_i3);
				cuDeleteFromDev(dev_ptr_mesh);
				cuDeleteFromDev(dev_ptr_nei);

				size_t sz = curLength + INC_VORT_DEV;

				wake.devWakePtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(sz, wake.devNWake);
				dev_ptr_vel = ReserveDevMem<double, 2>(sz, n_CUDA_vel);
				dev_ptr_rad = ReserveDevMem<double, 1>(sz, n_CUDA_rad);

				dev_ptr_i0 = ReserveDevMem<double, 1>(sz, n_CUDA_i0);
				dev_ptr_i1 = ReserveDevMem<double, 1>(sz, n_CUDA_i1);
				dev_ptr_i2 = ReserveDevMem<double, 2>(sz, n_CUDA_i2);
				dev_ptr_i3 = ReserveDevMem<double, 2>(sz, n_CUDA_i3);

				dev_ptr_mesh = ReserveDevMem<int, 2>(sz, n_CUDA_mesh);
				dev_ptr_nei = ReserveDevMem<int, 1>(sz, n_CUDA_nei);

				vels.resize(n_CUDA_vel, { 0.0, 0.0 });
				rads.resize(n_CUDA_vel, 1.0e+5);

				i0.resize(n_CUDA_vel, 0.0);
				i1.resize(n_CUDA_vel, 0.0);
				i2.resize(n_CUDA_vel, { 0.0, 0.0 });
				i3.resize(n_CUDA_vel, { 0.0, 0.0 });
				nei.resize(n_CUDA_vel, 0);

				//mesh.resize(n_CUDA_vel, { 0, 0 });

				std::cout << "resize: " << n_CUDA_vel << std::endl;
			}

			cuClearWakeMem(wake.devNWake, wake.devWakePtr);
			cuCopyWakeToDev(wake.vtx.size(), wake.vtx.data(), wake.devWakePtr);
		}
	}	
}


void gpu::RefreshAfls()
{
	const int& id = parallel.myidWork;
	if (boundary.size() > 0)
	{
		//if (id == 0)
		{
			if (boundary.size() > n_CUDA_afls)
			{
				host_ptr_ptr_pnl.resize(boundary.size(), NULL);
				host_ptr_ptr_r.resize(boundary.size(), NULL);
				host_ptr_ptr_i0.resize(boundary.size(), NULL);
				host_ptr_ptr_i1.resize(boundary.size(), NULL);
				host_ptr_ptr_i2.resize(boundary.size(), NULL);
				host_ptr_ptr_i3.resize(boundary.size(), NULL);
				host_ptr_ptr_rad.resize(boundary.size(), NULL);
				host_ptr_ptr_vel.resize(boundary.size(), NULL);
				host_ptr_ptr_rhs.resize(boundary.size(), NULL);

				n_CUDA_pnls.resize(boundary.size(), 0);
				n_CUDA_rs.resize(boundary.size(), 0);
				n_CUDA_i0s.resize(boundary.size(), 0);
				n_CUDA_i1s.resize(boundary.size(), 0);
				n_CUDA_i2s.resize(boundary.size(), 0);
				n_CUDA_i3s.resize(boundary.size(), 0);
				n_CUDA_rads.resize(boundary.size(), 0);
				n_CUDA_vels.resize(boundary.size(), 0);
				n_CUDA_rhss.resize(boundary.size(), 0);

				for (size_t s = 0; s < n_CUDA_afls; ++s)
				{
					cuDeleteFromDev(host_ptr_ptr_pnl[s]);
					cuDeleteFromDev(host_ptr_ptr_r[s]);
					cuDeleteFromDev(host_ptr_ptr_i0[s]);
					cuDeleteFromDev(host_ptr_ptr_i1[s]);
					cuDeleteFromDev(host_ptr_ptr_i2[s]);
					cuDeleteFromDev(host_ptr_ptr_i3[s]);
					cuDeleteFromDev(host_ptr_ptr_rad[s]);
					cuDeleteFromDev(host_ptr_ptr_vel[s]);
					cuDeleteFromDev(host_ptr_ptr_rhs[s]);
				}

				if (n_CUDA_afls)
				{
					cuDeleteFromDev(dev_nPanels);
					cuDeleteFromDev(dev_ptr_ptr_pnl);
					cuDeleteFromDev(dev_ptr_ptr_r);
					cuDeleteFromDev(dev_ptr_ptr_i0);
					cuDeleteFromDev(dev_ptr_ptr_i1);
					cuDeleteFromDev(dev_ptr_ptr_i2);
					cuDeleteFromDev(dev_ptr_ptr_i3);
					cuDeleteFromDev(dev_ptr_ptr_rad);
					cuDeleteFromDev(dev_ptr_ptr_vel);
					cuDeleteFromDev(dev_ptr_ptr_rhs);
				}

				std::vector<size_t> host_nPanels(0);

				virtvels.resize(boundary.size());
				virtrads.resize(boundary.size());
				virti0.resize(boundary.size());
				virti1.resize(boundary.size());
				virti2.resize(boundary.size());
				virti3.resize(boundary.size());
				virtrhs.resize(boundary.size());

				n_CUDA_afls = 0;
				for (size_t s = 0; s < boundary.size(); ++s)
				{
					const size_t& sz = boundary[s]->afl.np;//->virtualWake.size();
					host_ptr_ptr_pnl[s] = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(sz, n_CUDA_pnls[s]);
					host_ptr_ptr_r[s] = ReserveDevMem<double, 2>(sz+1, n_CUDA_rs[s]);
					host_ptr_ptr_i0[s] = ReserveDevMem<double, 1>(sz, n_CUDA_i1s[s]);
					host_ptr_ptr_i1[s] = ReserveDevMem<double, 1>(sz, n_CUDA_i1s[s]);
					host_ptr_ptr_i2[s] = ReserveDevMem<double, 2>(sz, n_CUDA_i2s[s]);
					host_ptr_ptr_i3[s] = ReserveDevMem<double, 2>(sz, n_CUDA_i2s[s]);
					host_ptr_ptr_rad[s] = ReserveDevMem<double, 1>(sz, n_CUDA_rads[s]);
					host_ptr_ptr_vel[s] = ReserveDevMem<double, 2>(sz, n_CUDA_vels[s]);
					host_ptr_ptr_rhs[s] = ReserveDevMem<double, 1>(sz, n_CUDA_rhss[s]);

					n_CUDA_afls++;
					host_nPanels.push_back(n_CUDA_pnls[s]);

					virtvels[s].resize(n_CUDA_vels[s], { 0.0, 0.0 });
					virtrads[s].resize(n_CUDA_vels[s], 1.0e+5);

					virti0[s].resize(n_CUDA_vels[s], 0.0);
					virti1[s].resize(n_CUDA_vels[s], 0.0);
					virti2[s].resize(n_CUDA_vels[s], { 0.0, 0.0 });
					virti3[s].resize(n_CUDA_vels[s], { 0.0, 0.0 });
					virtrhs[s].resize(n_CUDA_vels[s], 0.0);
				}

				dev_nPanels = ReserveDevMemAndCopyFixedArray(boundary.size(), host_nPanels.data());

				dev_ptr_ptr_pnl = ReserveDevMemAndCopyFixedArray(boundary.size(), host_ptr_ptr_pnl.data());
				dev_ptr_ptr_r = ReserveDevMemAndCopyFixedArray(boundary.size(), host_ptr_ptr_r.data());
				dev_ptr_ptr_i0 = ReserveDevMemAndCopyFixedArray(boundary.size(), host_ptr_ptr_i0.data());
				dev_ptr_ptr_i1 = ReserveDevMemAndCopyFixedArray(boundary.size(), host_ptr_ptr_i1.data());
				dev_ptr_ptr_i2 = ReserveDevMemAndCopyFixedArray(boundary.size(), host_ptr_ptr_i2.data());
				dev_ptr_ptr_i3 = ReserveDevMemAndCopyFixedArray(boundary.size(), host_ptr_ptr_i3.data());
				dev_ptr_ptr_rad = ReserveDevMemAndCopyFixedArray(boundary.size(), host_ptr_ptr_rad.data());
				dev_ptr_ptr_vel = ReserveDevMemAndCopyFixedArray(boundary.size(), host_ptr_ptr_vel.data());
				dev_ptr_ptr_rhs = ReserveDevMemAndCopyFixedArray(boundary.size(), host_ptr_ptr_rhs.data());
			}

			for (size_t s = 0; s < boundary.size(); ++s)
			{
				cuClearWakeMem(n_CUDA_pnls[s], host_ptr_ptr_pnl[s]);
				cuCopyWakeToDev(boundary[s]->virtualWake.size(), boundary[s]->virtualWake.data(), host_ptr_ptr_pnl[s]);
				cuCopyRsToDev(boundary[s]->afl.r.size(), boundary[s]->afl.r.data(), host_ptr_ptr_r[s]);
			}
		}
	}
}


//Вычисление конвективных скоростей и радиусов вихревых доменов в заданном наборе точек
//MPI+CUDA
void gpu::ExpCalcConvVeloToSetOfPoints
(
	const size_t npt, double* dev_ptr_pt,
	const size_t nvt, double* dev_ptr_vt,
	const size_t nbou, size_t* dev_nPanels, double** dev_ptr_ptr_pnl, 
	std::vector<Point2D>& Vel, std::vector<double>& Rad,
	std::vector<Point2D>& locvel, std::vector<double>& locrad,
	double* dev_ptr_vel, double* dev_ptr_rad, 
	double minRad, double eps2)
{
	const int& id = parallel.myidWork;
	parProp par = parallel.SplitMPI(npt, true);
	
	double tCUDASTART = 0.0, tCUDAEND = 0.0;

	tCUDASTART = omp_get_wtime();

	if (npt > 0)
	{
		cuCalculateConvVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, nbou, dev_nPanels, dev_ptr_ptr_pnl, dev_ptr_vel, dev_ptr_rad, minRad, eps2);

		CopyMemFromDev<double, 2>(par.myLen, dev_ptr_vel, (double*)&locvel[0]);
			
		CopyMemFromDev<double, 1>(par.myLen, dev_ptr_rad, &locrad[0]);

		std::vector<Point2D> newV;
		if (id == 0)
			newV.resize(Vel.size());

		MPI_Gatherv(locvel.data(), par.myLen, Point2D::mpiPoint2D, newV.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
		if (id == 0)
		for (size_t q = 0; q < Vel.size(); ++q)
			Vel[q] += newV[q];
		
		Rad.resize(par.totalLen);
		MPI_Allgatherv(locrad.data(), par.myLen, MPI_DOUBLE,          Rad.data(), par.len.data(), par.disp.data(), MPI_DOUBLE,          parallel.commWork);		

		cuCopyFixedArray(dev_ptr_rad, Rad.data(), sizeof(double) * Rad.size());
	}

	tCUDAEND = omp_get_wtime();

	std::cout << "CONV_GPU(" << parallel.myidWork << "): " << (tCUDAEND - tCUDASTART) << std::endl;
}//ExpCalcConvVeloToSetOfPoints(...)


//MPI+CUDA
void gpu::ExpGetConvVelocityToSetOfPointsFromVirtualVortexes(
	const size_t npt, double* dev_ptr_pt,	
	const size_t nvt, double* dev_ptr_vt,
	std::vector<Point2D>& Vel, std::vector<Point2D>& locvel, double* dev_ptr_vel,
	double eps2)
{
	const int& id = parallel.myidWork;
	parProp par = parallel.SplitMPI(npt, true);

	double tCUDASTART = 0.0, tCUDAEND = 0.0;

	tCUDASTART = omp_get_wtime();

	if (npt > 0)
	{
		cuCalculateConvVeloWakeFromVirtual(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, dev_ptr_vel, eps2);

		CopyMemFromDev<double, 2>(par.myLen, dev_ptr_vel, (double*)&locvel[0]);

		std::vector<Point2D> newV;
		if (id == 0)
			newV.resize(Vel.size());

		MPI_Gatherv(locvel.data(), par.myLen, Point2D::mpiPoint2D, newV.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
		
		if (id == 0)
		for (size_t q = 0; q < Vel.size(); ++q)
			Vel[q] += newV[q];
	}

	tCUDAEND = omp_get_wtime();

	std::cout << "CONV_VIRT_GPU: " << (tCUDAEND - tCUDASTART) << std::endl;
}


//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void gpu::ExpCalcDiffVeloI1I2ToSetOfPoints(
	const size_t npt, double* dev_ptr_pt, 
	double* dev_ptr_dr,
	const size_t nvt, double* dev_ptr_vt,
	std::vector<double>& I1, std::vector<Point2D>& I2,
	std::vector<double>& loci1, std::vector<Point2D>& loci2,
	double* dev_ptr_i1, double* dev_ptr_i2,
	bool useMesh)
{
	const int& id = parallel.myidWork;
	parProp par = parallel.SplitMPI(npt, true);


	double tCUDASTART = 0.0, tCUDAEND = 0.0, tCUDAENDCALC = 0.0;

	tCUDASTART = omp_get_wtime();

	if (nvt > 0)
	{
		if (!useMesh)
			cuCalculateDiffVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);
		else
			cuCalculateDiffVeloWakeMesh(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, dev_ptr_mesh, wake.param.epscol, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);

		CopyMemFromDev<double, 2>(par.myLen, dev_ptr_i2, (double*)&loci2[0]);
		CopyMemFromDev<double, 1>(par.myLen, dev_ptr_i1, &loci1[0]);
		
		std::vector<Point2D> newI2;
		std::vector<double> newI1;
		if (id == 0)
		{
			newI2.resize(I2.size());
			newI1.resize(I1.size());
		}
		
		MPI_Gatherv(loci2.data(), par.myLen, Point2D::mpiPoint2D, newI2.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
		MPI_Gatherv(loci1.data(), par.myLen, MPI_DOUBLE, newI1.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, parallel.commWork);

		tCUDAENDCALC = omp_get_wtime();

		if (id == 0)
		for (size_t q = 0; q < I2.size(); ++q)
		{
			I1[q] += newI1[q];
			I2[q] += newI2[q];
		}
	}
	tCUDAEND = omp_get_wtime();
	std::cout << "DIFF_GPU(" << id << ", " << par.myLen << "): " << (tCUDAEND - tCUDASTART) << " " << (tCUDAENDCALC - tCUDASTART) << std::endl;
}//ExpCalcDiffVeloI1I2ToSetOfPoints(...)



 //Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void gpu::ExpGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(
	const size_t npt, double* dev_ptr_pt, double* dev_ptr_rad,
	const size_t nr, double* dev_ptr_r,
	std::vector<double>& I0, std::vector<Point2D>& I3,
	std::vector<double>& loci0, std::vector<Point2D>& loci3,
	double* dev_ptr_i0, double* dev_ptr_i3)
{
	const int& id = parallel.myidWork;
	parProp par = parallel.SplitMPI(npt, true);

	double tCUDASTART = 0.0, tCUDAEND = 0.0;

	tCUDASTART = omp_get_wtime();

	if ((npt > 0) && (nr > 0))
	{
		cuCalculateSurfDiffVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nr, dev_ptr_r, dev_ptr_i0, dev_ptr_i3, dev_ptr_rad);
		CopyMemFromDev<double, 2>(par.myLen, dev_ptr_i3, (double*)&loci3[0]);
		CopyMemFromDev<double, 1>(par.myLen, dev_ptr_i0, &loci0[0]);

		std::vector<Point2D> newI3;
		std::vector<double> newI0;
		if (id == 0)
		{
			newI3.resize(I3.size());
			newI0.resize(I0.size());
		}

		MPI_Gatherv(loci3.data(), par.myLen, Point2D::mpiPoint2D, newI3.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
		MPI_Gatherv(loci0.data(), par.myLen, MPI_DOUBLE, newI0.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, parallel.commWork);

		if (id == 0)
		for (size_t q = 0; q < I3.size(); ++q)
		{
			I0[q] += newI0[q];
			I3[q] += newI3[q];
		}
	}

	tCUDAEND = omp_get_wtime();

	std::cout << "DIFF_SURF_GPU: " << (tCUDAEND - tCUDASTART) << std::endl;
}//ExpGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(...)


void gpu::ExpGetWakeInfluence(
	const size_t npt, double* dev_ptr_pt, 
	const size_t nvt, double* dev_ptr_vt, 
	std::vector<double>& rhs, std::vector<double>& locrhs, double* dev_ptr_rhs)
{
	const int& id = parallel.myidWork;
	parProp par = parallel.SplitMPI(npt, true);

	double tCUDASTART = 0.0, tCUDAEND = 0.0;

	tCUDASTART = omp_get_wtime();

	rhs.resize(npt, 0.0);

	if (nvt > 0)
	{
		cuCalculateRhs(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, dev_ptr_rhs);			
		CopyMemFromDev<double, 1>(par.myLen, dev_ptr_rhs, (double*)&locrhs[0]);			

		std::vector<double> newRhs;
		if (id == 0)
		{
			newRhs.resize(rhs.size());
		}

		MPI_Gatherv(locrhs.data(), par.myLen, MPI_DOUBLE, newRhs.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, parallel.commWork);

		if (id == 0)
		for (size_t q = 0; q < rhs.size(); ++q)
			rhs[q] = newRhs[q];
	}
	tCUDAEND = omp_get_wtime();
	std::cout << "RHS_GPU: " << (tCUDAEND - tCUDASTART) << std::endl;
}//ExpGetWakeInfluence(...)



void gpu::ExpGetPairs(int type, std::vector<int>& NEIB)
{
	size_t npt = wake.vtx.size();
	const int& id = parallel.myidWork;
	parProp par = parallel.SplitMPI(npt, true);

	double tCUDASTART = 0.0, tCUDAEND = 0.0;

	tCUDASTART = omp_get_wtime();

	NEIB.resize(npt, 0);


	if (npt > 0)
	{
		cuCalculatePairs(par.myDisp, par.myLen, npt, wake.devWakePtr, dev_ptr_mesh, dev_ptr_nei, wake.param.epscol, sqr(wake.param.epscol), type);

		CopyMemFromDev<int, 1>(par.myLen, dev_ptr_nei, &nei[0]);

		std::vector<int> newNei;

		newNei.resize(nei.size());

		MPI_Allgatherv(nei.data(), par.myLen, MPI_INT, newNei.data(), par.len.data(), par.disp.data(), MPI_INT,  parallel.commWork);

		for (size_t q = 0; q < NEIB.size(); ++q)
		{
			NEIB[q] = newNei[q];
			//std::cout << q << " " << NEIB[q] << std::endl;
		}
	}
	
	tCUDAEND = omp_get_wtime();

	std::cout << "GPU_Pairs: " << (tCUDAEND - tCUDASTART) << std::endl;
}

#endif