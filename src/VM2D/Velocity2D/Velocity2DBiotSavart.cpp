/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Velocity2DBiotSavart.cpp                                         |
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
\brief Файл кода с описанием класса VelocityBiotSavart
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#include "Velocity2DBiotSavart.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Wake2D.h"
#include "World2D.h"

#include "wrapper.h"

#include <algorithm>

using namespace VM2D;

/// Конструктор
VelocityBiotSavart::VelocityBiotSavart(const World2D& W_) :
	Velocity(W_)
{
};

/// Деструктор
VelocityBiotSavart::~VelocityBiotSavart()
{	
};

//Вычисление конвективных скоростей вихрей и виртуальных вихрей в вихревом следе, а также в точках wakeVP
void VelocityBiotSavart::CalcConvVelo()
{
	//const int checkStep = 3;


	double tt1 = omp_get_wtime();
	
	W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONV_TOWAKE))	
	//if (W.getCurrentStep() == checkStep)
	GPUWakeToWakeFAST(W.getWake(), wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake, true, true);
	//else
	//GPUCalcConvVeloToSetOfPointsFromWake(W.getWake(), wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake, true, true);	
#else
	CalcConvVeloToSetOfPointsFromWake(W.getWake(), wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake, true, true);
#endif
	



	double tt2 = omp_get_wtime();


/*
	if (W.getCurrentStep() == checkStep)
	{
		std::ofstream pos("checkPos.txt");
		pos.precision(17);
		pos << W.getWake().vtx.size() << std::endl;
		for (int i = 0; i < W.getWake().vtx.size(); ++i)
			pos  << W.getWake().vtx[i].r()[0] << " " << W.getWake().vtx[i].r()[1] << " " << W.getWake().vtx[i].g() << std::endl;
		pos.close();

		std::ofstream of("checkVelo.txt");
		of.precision(16);
		for (int i = 0; i < wakeVortexesParams.convVelo.size(); ++i)
			of << i << " " << wakeVortexesParams.convVelo[i][0] << " " << wakeVortexesParams.convVelo[i][1] << " " \
			<< wakeVortexesParams.epsastWake[i] << std::endl;
		of.close();
		exit(100500);
	}
*/

	//W.getInfo('t') << "Convective velocities time = " << int(((tt2 - tt1) * 1000)*10) / 10.0 << " ms, bodies = " << W.getWake().vtx.size() << std::endl;


	std::vector<Point2D> nullVector(0);

	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONV_TOBOU))	
		GPUCalcConvVeloToSetOfPointsFromWake(W.getBoundary(bou).virtualWake, nullVector, virtualVortexesParams[bou].epsastWake, false, true);
#else
		CalcConvVeloToSetOfPointsFromWake(W.getBoundary(bou).virtualWake, nullVector, virtualVortexesParams[bou].epsastWake, false, true);
#endif
	}
	
	//double tt3 = omp_get_wtime();

//вычисление конвективных скоростей по закону Био-Савара от виртуальных вихрей
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		//Влияние на пелену от виртуальных вихрей
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONVVIRT))
		W.getBoundary(bou).GPUCalcConvVelocityToSetOfPointsFromSheets(W.getWake(), wakeVortexesParams.convVelo);
#else
		W.getBoundary(bou).CalcConvVelocityToSetOfPointsFromSheets(W.getWake(), wakeVortexesParams.convVelo);
#endif
	}
		
	//double tt4 = omp_get_wtime();

	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		//Скороcти самих виртуальных вихрей
		W.getBoundary(bou).CalcConvVelocityAtVirtualVortexes(virtualVortexesParams[bou].convVelo);
	}

	//double tt5 = omp_get_wtime();

	//std::cout << tt2 - tt1 << " " << tt3 - tt2 << " " << tt4 - tt3 << " " << tt5 - tt4 << std::endl;

	W.getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();
	
	W.getTimestat().timeVP.first += omp_get_wtime();

	//вычисление скоростей в заданных точках только на соответствующем шаге по времени
	if ((W.getPassport().timeDiscretizationProperties.saveVPstep > 0) && (!(W.getCurrentStep() % W.getPassport().timeDiscretizationProperties.saveVPstep)))
		CalcVeloToWakeVP();	
	W.getTimestat().timeVP.second += omp_get_wtime();	
}//CalcConvVelo()


//Вычисление конвективных скоростей вихрей в точках wakeVP
void VelocityBiotSavart::CalcVeloToWakeVP()
{
	std::vector<Point2D> velConvWake;
	std::vector<std::vector<Point2D>> velConvBou;

	int addWSize = (int)W.getMeasureVP().getWakeVP().vtx.size();

	velConvWake.resize(addWSize, { 0.0, 0.0 });

	velConvBou.resize(W.getNumberOfBoundary());
	for (size_t i = 0; i < W.getNumberOfBoundary(); ++i)
		velConvBou[i].resize(addWSize, { 0.0, 0.0 });

#if (defined(USE_CUDA))
	W.getNonConstCuda().RefreshVP();
#endif

#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_VP))
	GPUCalcConvVeloToSetOfPointsFromWake(W.getMeasureVP().getWakeVP(), velConvWake, W.getNonConstMeasureVP().getNonConstDomainRadius(), true, false);
#else
	CalcConvVeloToSetOfPointsFromWake(W.getMeasureVP().getWakeVP(), velConvWake, W.getNonConstMeasureVP().getNonConstDomainRadius(), true, false);
#endif			   

	for (size_t i = 0; i < W.getNumberOfBoundary(); ++i)
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_VP))		
		W.getNonConstBoundary(i).GPUCalcConvVelocityToSetOfPointsFromSheets(W.getMeasureVP().getWakeVP(), velConvBou[i]);
#else
		W.getNonConstBoundary(i).CalcConvVelocityToSetOfPointsFromSheets(W.getMeasureVP().getWakeVP(), velConvBou[i]);
#endif			   

	std::vector<Point2D>& velocityRef = W.getNonConstMeasureVP().getNonConstVelocity();
	velocityRef.assign(addWSize, W.getPassport().physicalProperties.V0());

	for (int i = 0; i < addWSize; ++i)
		velocityRef[i] += velConvWake[i];

	for (size_t bou = 0; bou < velConvBou.size(); ++bou)
		for (int j = 0; j < addWSize; ++j)
			velocityRef[j] += velConvBou[bou][j];

}// CalcVeloToWakeVP()
 

inline void ModifyE2(double* ee2, double dst2)
{
	if (dst2 > 0)
	{
		if (dst2<ee2[0])
		{
			ee2[2] = ee2[1];
			ee2[1] = ee2[0];
			ee2[0] = dst2;
		}//if (dist2<ee2[0])
		else
		{
			if (dst2<ee2[1])
			{
				ee2[2] = ee2[1];
				ee2[1] = dst2;
			}// if (dist2<ee2[1])
			else
			if (dst2<ee2[2])
				ee2[2] = dst2;
		}//else
	}//if (dst2>0)
}

//Вычисление конвективных скоростей и радиусов вихревых доменов в заданном наборе точек от следа
void VelocityBiotSavart::CalcConvVeloToSetOfPointsFromWake(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius)
{
	std::vector<Point2D> selfVelo(pointsDb.vtx.size());
	domainRadius.resize(pointsDb.vtx.size());

	double cft = IDPI;

#pragma warning (push)
#pragma warning (disable: 4101)
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2eps, dst2;	
#pragma warning (pop)

	double eps2 = W.getPassport().wakeDiscretizationProperties.eps2;
	
	if (calcVelo)
	{
#pragma omp parallel for default(none) shared(selfVelo, cft, calcVelo, calcRadius, eps2, pointsDb, domainRadius) private(tempVel, velI, dst2, dst2eps) schedule(dynamic, DYN_SCHEDULE)
		for (int i = 0; i < pointsDb.vtx.size(); ++i)
		{
			double ee2[3] = { 10000.0, 10000.0, 10000.0 };

			velI.toZero();

			const Point2D& posI = pointsDb.vtx[i].r();

			for (size_t j = 0; j < W.getWake().vtx.size(); ++j)
			{
				const Point2D& posJ = W.getWake().vtx[j].r();

				dst2 = (posI-posJ).length2();

				//Модифицируем массив квадратов расстояний до ближайших вихрей из wake
#ifndef TESTONLYVELO
				if (calcRadius)
					VMlib::ModifyE2(ee2, dst2);
#endif // !TESTONLYVELO				

				const double& gamJ = W.getWake().vtx[j].g();

				tempVel.toZero();
				dst2eps = VMlib::boundDenom(dst2, eps2); //Сглаживать надо!!!				

				tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
				tempVel *= (gamJ / dst2eps);
				velI += tempVel;
			}
			
			for (size_t j = 0; j < W.getSource().vtx.size(); ++j)
			{
				const Point2D& posJ = W.getSource().vtx[j].r();
				const double& gamJ = W.getPassport().physicalProperties.accelCft() * W.getSource().vtx[j].g();

				tempVel.toZero();

				dst2 = dist2(posI, posJ);
				dst2eps = VMlib::boundDenom(dst2, W.getPassport().wakeDiscretizationProperties.eps2); //Сглаживать надо!!!

				tempVel = { posI[0] - posJ[0], posI[1] - posJ[1] };
				tempVel *= (gamJ / dst2eps);
				velI += tempVel;
			}
#ifndef TESTONLYVELO
			if (calcRadius)
			{
				for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
				{
					const auto& bou = W.getBoundary(s);
					//Модифицируем массив квадратов расстояний до ближайших вихрей из virtualWake					
					for (size_t j = 0; j < bou.virtualWake.vtx.size(); ++j)						
						ModifyE2(ee2, dist2(posI, bou.virtualWake.vtx[j].r()));			
				}
			}
#endif

			velI *= cft;
			selfVelo[i] = velI;

#ifndef TESTONLYVELO
			if (calcRadius)
				domainRadius[i] = 1.0 * sqrt((ee2[0] + ee2[1] + ee2[2]) / 3.0);
#endif
		}
	}
	else if (calcRadius)
	{
#pragma omp parallel for default(none) shared(selfVelo, cft, calcVelo, calcRadius, pointsDb, domainRadius) private(tempVel, velI, dst2) schedule(dynamic, DYN_SCHEDULE)
		for (int i = 0; i < pointsDb.vtx.size(); ++i)
		{
			double ee2[3] = { 10000.0, 10000.0, 10000.0 };

			velI.toZero();

			const Point2D& posI = pointsDb.vtx[i].r();

			for (size_t j = 0; j < W.getWake().vtx.size(); ++j)
			{
				const Point2D& posJ = W.getWake().vtx[j].r();

				dst2 = dist2(posI, posJ);

				//Модифицируем массив квадратов расстояний до ближайших вихрей из wake
				VMlib::ModifyE2(ee2, dst2);				
			}

			for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
			{
				for (size_t j = 0; j < W.getBoundary(s).virtualWake.vtx.size(); ++j)
				{
					const Point2D& posJ = W.getBoundary(s).virtualWake.vtx[j].r();
					dst2 = dist2(posI, posJ);

					//Модифицируем массив квадратов расстояний до ближайших вихрей из virtualWake
					ModifyE2(ee2, dst2);
					/*
					if (dst2 < ee2[0])
					{
						size_t pnl = W.getBoundary(s).virtualWake.aflPan[j].second;
						ee2[0] = ee2[1] = ee2[2] = 0.5 * W.getBoundary(s).afl.len[pnl] / (W.getBoundary(s).vortexBeginEnd[pnl].second - W.getBoundary(s).vortexBeginEnd[pnl].first);						
					}
					else
						ModifyE2(ee2, dst2);
					*/

				}
			}

			domainRadius[i] = 1.0 * sqrt((ee2[0] + ee2[1] + ee2[2]) / 3.0);
		}
	} //else


	if (calcVelo)
		for (size_t i = 0; i < velo.size(); ++i)
			velo[i] += selfVelo[i];
}//CalcConvVeloToSetOfPointsFromWake(...)


#if defined(USE_CUDA)
void VelocityBiotSavart::GPUWakeToWakeFAST(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius)
{
		size_t npt = pointsDb.vtx.size();
		double*& dev_ptr_pt = pointsDb.devVtxPtr;

		const size_t nvt = W.getWake().vtx.size();
		double*& dev_ptr_vt = W.getWake().devVtxPtr;
		
		const size_t nsr = W.getSource().vtx.size();
		double*& dev_ptr_sr = W.getSource().devVtxPtr;
		const size_t nbou = W.getNumberOfBoundary();

		//size_t* const& dev_nPanels = W.getCuda().dev_ptr_nPanels;
		size_t* const& dev_nVortices = W.getCuda().dev_ptr_nVortices;

		double** const& dev_ptr_ptr_vtx = W.getCuda().dev_ptr_ptr_vtx;

		std::vector<Point2D> Vel(npt);
		std::vector<double> Rad(npt);

		std::vector<Point2D> newV(npt);

		
		double*& dev_ptr_vel = pointsDb.devVelPtr;
		double*& dev_ptr_rad = pointsDb.devRadPtr;
		const double& eps2 = W.getPassport().wakeDiscretizationProperties.eps2;

		if (npt > 0)
		{
			//double t1 = omp_get_wtime();
			//cuCalculateConvVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, nsr, dev_ptr_sr, nbou, dev_nVortices, dev_ptr_ptr_vtx, dev_ptr_vel, dev_ptr_rad, eps2, calcVelo, calcRadius);
			
			double timings[7];

			wrapperInfluence((Vortex2D*)dev_ptr_pt, (Point2D*)dev_ptr_vel, dev_ptr_rad,
				W.getNonConstCuda().CUDAptrs, (int)npt, timings, sqrt(eps2), 1.30,
				W.getNonConstCuda().n_CUDA_bodies, (int)W.getNonConstCuda().n_CUDA_wake, 8,
				nbou, dev_nVortices, dev_ptr_ptr_vtx);

			//double t2 = omp_get_wtime();
			//std::cout << "KERNEL_TIME = " << t2 - t1 << std::endl;


			if (calcVelo)
			{
				//double tt1 = omp_get_wtime();

				W.getCuda().CopyMemFromDev<double, 2>(npt, dev_ptr_vel, (double*)newV.data(), 20);

				
/*
				////////////////////////////////////////////////////////
				size_t maxp = W.getWake().vtx.size();
				if (par.myLen >= maxp)
				{
					std::vector<Vortex2D> controlPoints_h(maxp);
					for (int i = 0; i < maxp; ++i)
					{
						controlPoints_h[i].r() = W.getWake().vtx[i].r();
						controlPoints_h[i].g() = 0.0;
					}
					
					realVortex* ptr_points;
					realPoint* ptr_velos;
					double* ptr_epsast;

					size_t maxpUp;
										
					ptr_points = (realVortex*)W.getNonConstCuda().ReserveDevMem<double, 3>(maxp, maxpUp);
					ptr_velos = (realPoint*)W.getNonConstCuda().ReserveDevMem<double, 2>(maxp, maxpUp);
					ptr_epsast = W.getNonConstCuda().ReserveDevMem<double, 1>(maxp, maxpUp);
										
					W.getNonConstCuda().CopyMemToDev<double, 3>(maxp, (double*)controlPoints_h.data(), (double*)ptr_points);

					double timingsToPoints[7];

					wrapperInfluenceToPoints((Vortex2D*)dev_ptr_pt, (Vortex2D*)ptr_points, (Point2D*)ptr_velos, ptr_epsast,
						W.getNonConstCuda().CUDAptrs, false, 
						(int)npt, (int)maxp, timingsToPoints, sqrt(eps2), 0.31,
						W.getNonConstCuda().n_CUDA_bodies, (int)W.getNonConstCuda().n_CUDA_wake, 14,
						nbou, dev_nVortices, dev_ptr_ptr_vtx);


					std::cout << "Time_base = " << timings[6] << ", Time_points = " << timingsToPoints[6] << std::endl;


					std::vector<Point2D> velos_h(maxp);
					W.getNonConstCuda().CopyMemFromDev<double, 2>(maxp, (double*)ptr_velos, (double*)&velos_h[0], 20);

					std::cout.precision(16);

					//for (int i = 0; i < maxp; ++i)
					//{
					//	std::cout << i << ": " << "vel = " << locvel[i] << ", new_vel = " << velos_h[i] << ", dv = " << locvel[i] - velos_h[i] << std::endl;
					//}

					W.getNonConstCuda().ReleaseDevMem(ptr_points, 800);
					W.getNonConstCuda().ReleaseDevMem(ptr_velos, 801);
					W.getNonConstCuda().ReleaseDevMem(ptr_epsast, 802);


					//if (maxp > 100)
					//	exit(-300);
				}
				*/
				///////////////////////////////////////////////////////////////////
				


				//double tt2 = omp_get_wtime();
				//std::cout << "COPY_TIME = " << tt2 - tt1 << std::endl;

				//double tt3 = omp_get_wtime();

				
				for (size_t q = 0; q < npt; ++q)
					Vel[q] = newV[q];

				for (size_t q = 0; q < npt; ++q)
					velo[q] += Vel[q];

				//double tt4 = omp_get_wtime();
				//std::cout << "GATHER_TIME = " << tt4 - tt3 << std::endl;

			}//if calcVelo

			//double tt5 = omp_get_wtime();

			if (calcRadius)
			{
				Rad.resize(npt);
				W.getCuda().CopyMemFromDev<double, 1>(npt, dev_ptr_rad, Rad.data(), 212);				

				//cuCopyFixedArray(dev_ptr_rad, Rad.data(), sizeof(double) * Rad.size());

				for (size_t q = 0; q < Rad.size(); ++q)
					domainRadius[q] = Rad[q];
					
				//std::ofstream filo("filo.txt");
				//for (int i = 0; i < npt; ++i)
				//	filo << i << " " << Rad[i] << std::endl;
				//	//filo << i << " " << locvel[i][0] << " " << locvel[i][1] << std::endl;
				//filo.close();

			}//if calcRadius

			//double tt6 = omp_get_wtime();
			//std::cout << "RADIUS_TIME = " << tt6 - tt5 << std::endl;

		}//if npt > 0
}
#endif


#if defined(USE_CUDA)
void VelocityBiotSavart::GPUCalcConvVeloToSetOfPointsFromWake(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius)
{	
	if ((&pointsDb == &W.getWake()) || (&pointsDb == &W.getBoundary(0).virtualWake) || (&pointsDb == &W.getMeasureVP().getWakeVP()))
	{		
		size_t npt = pointsDb.vtx.size();
		double*& dev_ptr_pt = pointsDb.devVtxPtr;

		if ((W.getNumberOfBoundary() >0) && (&pointsDb == &W.getBoundary(0).virtualWake))
		{
			for (size_t q = 1; q < W.getNumberOfBoundary(); ++q)
			npt += W.getBoundary(q).virtualWake.vtx.size();
		}
		
		const size_t nvt = W.getWake().vtx.size();
		double*& dev_ptr_vt = W.getWake().devVtxPtr;
		const size_t nsr = W.getSource().vtx.size();
		double*& dev_ptr_sr = W.getSource().devVtxPtr;
		const size_t nbou = W.getNumberOfBoundary();

		//size_t* const& dev_nPanels = W.getCuda().dev_ptr_nPanels;
		size_t* const& dev_nVortices = W.getCuda().dev_ptr_nVortices;

		double** const& dev_ptr_ptr_vtx = W.getCuda().dev_ptr_ptr_vtx;

		std::vector<Point2D> Vel(npt);
		std::vector<double> Rad(npt);

		std::vector<Point2D> newV(npt);
		
		double*& dev_ptr_vel = pointsDb.devVelPtr;
		double*& dev_ptr_rad = pointsDb.devRadPtr;
		const double& eps2 = W.getPassport().wakeDiscretizationProperties.eps2;
		
		if (npt > 0)
		{
			//double t1 = omp_get_wtime();
			cuCalculateConvVeloWake(npt, dev_ptr_pt, nvt, dev_ptr_vt, nsr, dev_ptr_sr, nbou, dev_nVortices, dev_ptr_ptr_vtx, dev_ptr_vel, dev_ptr_rad, eps2, calcVelo, calcRadius);
			//double t2 = omp_get_wtime();
			//std::cout << "KERNEL_TIME = " << t2 - t1 << std::endl;
			
			if (calcVelo)
			{
				//double tt1 = omp_get_wtime();
				
				W.getCuda().CopyMemFromDev<double, 2>(npt, dev_ptr_vel, (double*)newV.data(), 20);
			
				//double tt2 = omp_get_wtime();
				//std::cout << "COPY_TIME = " << tt2 - tt1 << std::endl;
				
				//double tt3 = omp_get_wtime();

				for (size_t q = 0; q < npt; ++q)
					Vel[q] = newV[q];

				if ((&pointsDb == &W.getWake()) || (&pointsDb == &W.getMeasureVP().getWakeVP()))
				{
					for (size_t q = 0; q < npt; ++q)
						velo[q] += Vel[q];
				}//if &pointsDb


				if ((W.getNumberOfBoundary() > 0) && (&pointsDb == &W.getBoundary(0).virtualWake))
				{
					//Сюда в принципе не должно быть попадания
					exit(123);
				}//if &pointsDb				

				//double tt4 = omp_get_wtime();
				//std::cout << "GATHER_TIME = " << tt4 - tt3 << std::endl;

			}//if calcVelo

			//double tt5 = omp_get_wtime();

			if (calcRadius)
			{
				Rad.resize(npt);
				W.getCuda().CopyMemFromDev<double, 1>(npt, dev_ptr_rad, Rad.data(), 211);

				//cuCopyFixedArray(dev_ptr_rad, Rad.data(), sizeof(double) * Rad.size());
				
				
				if ((&pointsDb == &W.getWake()) || (&pointsDb == &W.getMeasureVP().getWakeVP()))
				{
					for (size_t q = 0; q < Rad.size(); ++q)
						domainRadius[q] = Rad[q];
				}//if &pointsDb
				
				
				if ((W.getNumberOfBoundary() > 0) && (&pointsDb == &W.getBoundary(0).virtualWake))
				{
					size_t curGlobPnl = 0;
					for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s) //W.getNumberOfAirfoil()
					{
						size_t nv = W.getVelocity().virtualVortexesParams[s].epsastWake.size();
						for (size_t q = 0; q < nv; ++q)
							W.getNonConstVelocity().virtualVortexesParams[s].epsastWake[q] = Rad[curGlobPnl + q];
					
						curGlobPnl += nv;
					}//for s
				}//if &pointsDb
			}//if calcRadius
			
			//double tt6 = omp_get_wtime();
			//std::cout << "RADIUS_TIME = " << tt6 - tt5 << std::endl;

		}//if npt > 0
	}//if &pointsDb
}//GPUCalcConvVeloToSetOfPointsFromWake(...)
#endif



//Генерация вектора влияния вихревого следа на профиль
void VelocityBiotSavart::GetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeRhs) const
{
	size_t np = afl.getNumberOfPanels();
	size_t shDim = W.getBoundary(afl.numberInPassport).sheetDim;

	wakeRhs.resize(W.getBoundary(afl.numberInPassport).GetUnknownsSize());

	//локальные переменные для цикла	
	std::vector<double> velI(shDim, 0.0);

#pragma omp parallel for default(none) shared(shDim, afl, np, wakeRhs, IDPI) private(velI)
	for (int i = 0; i < np; ++i)
	{
		velI.assign(shDim, 0.0);

		if (W.getWake().vtx.size() > 0)
		{
			//Учет влияния следа
			afl.GetInfluenceFromVorticesToPanel(i, W.getWake().vtx.data(), W.getWake().vtx.size(), velI);
		}

		if (W.getSource().vtx.size() > 0)
		{
			//Учет влияния источников
			afl.GetInfluenceFromSourcesToPanel(i, W.getSource().vtx.data(), W.getSource().vtx.size(), velI);
		}

		for (size_t j = 0; j < shDim; ++j)
			velI[j] *= IDPI / afl.len[i];

		wakeRhs[i] = velI[0];

		if (shDim != 1)
			wakeRhs[np + i] = velI[1];
	}//for i
}//GetWakeInfluenceToRhs(...)



#if defined(USE_CUDA)
//Генерация вектора влияния вихревого следа на профиль
void VelocityBiotSavart::GPUGetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeVelo) const
{
	size_t shDim = W.getBoundary(afl.numberInPassport).sheetDim;

	const size_t& nvt = W.getWake().vtx.size();
	const size_t& nsr = W.getSource().vtx.size();
	const double eps2 = W.getPassport().wakeDiscretizationProperties.eps2;

	if (afl.numberInPassport == 0)
	{
		size_t nTotPan = 0;
		for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s)
			nTotPan += W.getAirfoil(s).getNumberOfPanels();

		double*& dev_ptr_pt = afl.devRPtr;
		double*& dev_ptr_vt = W.getWake().devVtxPtr;
		double*& dev_ptr_sr = W.getSource().devVtxPtr;
		double*& dev_ptr_rhs = afl.devRhsPtr;
		double*& dev_ptr_rhsLin = afl.devRhsLinPtr;

		std::vector<double> locrhs(nTotPan);
		std::vector<double> locrhsLin(nTotPan);



		if ((nvt > 0) || (nsr > 0))
		{
			cuCalculateRhs(nTotPan, dev_ptr_pt, nvt, dev_ptr_vt, nsr, dev_ptr_sr, eps2, dev_ptr_rhs, dev_ptr_rhsLin);

			std::vector<double> newRhs(nTotPan), newRhsLin(nTotPan);

			W.getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_rhs, newRhs.data(), 22);
			if (shDim != 1)
				W.getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_rhsLin, newRhsLin.data(), 22);


			//for (int qq = 0; qq < (int)newRhs.size(); ++qq)
			//	std::cout << "qq = " << qq << ", " << newRhs[qq] << std::endl;

			size_t curGlobPnl = 0;
			for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s)
			{
				std::vector<double>& tmpRhs = W.getNonConstAirfoil(s).tmpRhs;
				const size_t& np = W.getAirfoil(s).getNumberOfPanels();
				tmpRhs.resize(0);
				tmpRhs.insert(tmpRhs.end(), newRhs.begin() + curGlobPnl, newRhs.begin() + curGlobPnl + np);
				if (shDim != 1)
					tmpRhs.insert(tmpRhs.end(), newRhsLin.begin() + curGlobPnl, newRhsLin.begin() + curGlobPnl + np);

				curGlobPnl += np;
			}
		}
	}

	if ((nvt > 0) || (nsr > 0))
		wakeVelo = std::move(afl.tmpRhs);
	else
		wakeVelo.resize(afl.getNumberOfPanels() * (W.getPassport().numericalSchemes.boundaryCondition.second + 1), 0.0);

}//GPUGetWakeInfluenceToRhs(...)
#endif

#if defined(USE_CUDA)
//Генерация вектора влияния вихревого следа на профиль
void VelocityBiotSavart::GPUFASTGetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeVelo) const
{
	size_t shDim = W.getBoundary(afl.numberInPassport).sheetDim;

	const size_t& nvt = W.getWake().vtx.size();
	const size_t& nsr = W.getSource().vtx.size();


	if (afl.numberInPassport == 0)
	{
		size_t nTotPan = 0;
		for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s)
			nTotPan += W.getAirfoil(s).getNumberOfPanels();

		double*& dev_ptr_pt = afl.devRPtr;
		double*& dev_ptr_vt = W.getWake().devVtxPtr;
		double*& dev_ptr_sr = W.getSource().devVtxPtr;
		double*& dev_ptr_rhs = afl.devRhsPtr;
		double*& dev_ptr_rhsLin = afl.devRhsLinPtr;
		std::vector<double> locrhs(nTotPan);
		std::vector<double> locrhsLin(nTotPan);

		if ((nvt > 0) || (nsr > 0))
		{
			//cuCalculateRhs(par.myDisp, par.myLen, nTotPan, dev_ptr_pt, nvt, dev_ptr_vt, nsr, dev_ptr_sr, eps2, dev_ptr_rhs, dev_ptr_rhsLin);

			double timingsToRHS[7];

			wrapperInfluenceToRHS(
				(Vortex2D*)dev_ptr_vt,  //вихри в следе
				(double*)dev_ptr_pt,    //начала и концы панелей
				(double*)dev_ptr_rhs,   //куда сохранить результат 
				(double*)((shDim == 1) ? nullptr : dev_ptr_rhsLin),
				W.getNonConstCuda().CUDAptrs,  //указатели на дерево вихрей
				false,                  //признак перестроения дерева вихрей

				(int)nvt,               //число вихрей в следе
				(int)nTotPan,           //общее число панелей на всех профилях
				timingsToRHS,           //засечки времени
				//sqrt(eps2),             //eps
				1.2,                    //theta
				W.getNonConstCuda().n_CUDA_bodies,    //для следа
				(int)W.getNonConstCuda().n_CUDA_wake, //для следа
				8,                    //order
				W.getPassport().numericalSchemes.boundaryCondition.second
			);

			//printf("timings = %f: ( %f, %f, %f, %f, %f, %f )\n", timingsToRHS[6],
			//	timingsToRHS[0], timingsToRHS[1], timingsToRHS[2], timingsToRHS[3], timingsToRHS[4], timingsToRHS[5]);

			std::vector<double> newRhs(nTotPan), newRhsLin(nTotPan);

			W.getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_rhs, newRhs.data(), 22);
			if (shDim != 1)
				W.getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_rhsLin, newRhsLin.data(), 22);


			//for (int qq = 0; qq < (int)newRhs.size(); ++qq)
			//	std::cout << "qq = " << qq << ", " << newRhs[qq] << std::endl;

			size_t curGlobPnl = 0;
			for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s)
			{
				std::vector<double>& tmpRhs = W.getNonConstAirfoil(s).tmpRhs;
				const size_t& np = W.getAirfoil(s).getNumberOfPanels();
				tmpRhs.resize(0);
				tmpRhs.insert(tmpRhs.end(), newRhs.begin() + curGlobPnl, newRhs.begin() + curGlobPnl + np);
				if (shDim != 1)
					tmpRhs.insert(tmpRhs.end(), newRhsLin.begin() + curGlobPnl, newRhsLin.begin() + curGlobPnl + np);

				curGlobPnl += np;
			}

		}
	}

	if ((nvt > 0) || (nsr > 0))
		wakeVelo = std::move(afl.tmpRhs);
	else
		wakeVelo.resize(afl.getNumberOfPanels() * (W.getPassport().numericalSchemes.boundaryCondition.second + 1), 0.0);

}//GPUGetWakeInfluenceToRhsFAST(...)
#endif



void VelocityBiotSavart::FillRhs(Eigen::VectorXd& rhs) const
{
	Eigen::VectorXd locRhs;
	std::vector<double> lastRhs(W.getNumberOfBoundary());

	size_t currentRow = 0;

	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		//double tt0 = omp_get_wtime();

		const Airfoil& afl = W.getAirfoil(bou);
		size_t np = afl.getNumberOfPanels();

		size_t nVars;

		nVars = W.getBoundary(bou).GetUnknownsSize();
		locRhs.resize(nVars);


		std::vector<double> wakeRhs;

		double tt1 = omp_get_wtime();

#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_RHS))
		//GPUGetWakeInfluenceToRhs(afl, wakeRhs);		
		GPUFASTGetWakeInfluenceToRhs(afl, wakeRhs);
		//GetWakeInfluenceToRhs(afl, wakeRhs);

#else
		GetWakeInfluenceToRhs(afl, wakeRhs);
#endif			


		double tt2 = omp_get_wtime();
		//W.getInfo('t') << "Rhs time = " << int(((tt2 - tt1) * 1000) * 10) / 10.0 << " ms, bodies = " << W.getWake().vtx.size() << std::endl;


		std::vector<double> vInfRhs;
		afl.GetInfluenceFromVInfToPanel(vInfRhs);

		/*
		//if (W.currentStep == 200)
		{
			std::stringstream ss;
			ss << "rhsWake-" << W.currentStep;
			std::ofstream of(W.getPassport().dir + "dbg/" + ss.str());
			for (size_t i = 0; i < wakeRhs.size(); ++i)
				of << wakeRhs[i] << std::endl;
			of.close();
		}
		*/



#pragma omp parallel for \
	default(none) \
	shared(locRhs, afl, bou, wakeRhs, vInfRhs, np) \
	schedule(dynamic, DYN_SCHEDULE)
		for (int i = 0; i < (int)afl.getNumberOfPanels(); ++i)
		{
			locRhs(i) = -vInfRhs[i] - wakeRhs[i] + 0.25 * ((afl.getV(i) + afl.getV(i + 1)) & afl.tau[i]); //0.25 * (afl.getV(i) + afl.getV(i + 1))*afl.tau[i] - прямолинейные
			if (W.getBoundary(bou).sheetDim > 1)
				locRhs(np + i) = -vInfRhs[np + i] - wakeRhs[np + i];

			//влияние присоединенных слоев от самого себя и от других профилей				
			for (size_t q = 0; q < W.getNumberOfBoundary(); ++q)
			{
				const auto& sht = W.getBoundary(q).sheets;
				const auto& iq = W.getIQ(bou, q);

				const Airfoil& aflOther = W.getAirfoil(q);
				if (W.getMechanics(q).isDeform || W.getMechanics(q).isMoves)
				{
					for (size_t j = 0; j < aflOther.getNumberOfPanels(); ++j)
					{
						if ((i != j) || (bou != q))
						{
							locRhs(i) += -iq.first(i, j) * sht.attachedVortexSheet(j, 0);
							locRhs(i) += -iq.second(i, j) * sht.attachedSourceSheet(j, 0); // getIQ(bou, q).second(i, j) пока забито нулем для криволинейных
						}//if (i != j)
					}//for j
				}
			}//for q
		}//for i

		lastRhs[bou] = 0.0;

		for (size_t q = 0; q < afl.gammaThrough.size(); ++q)
			lastRhs[bou] += afl.gammaThrough[q];




		////////////////////////////////////////////////////////////////

		double dwcm; //wcm0, wcm1, wcmOld0, wcmOld1, dwc0, dwc1;
		const double currT = W.getPassport().timeDiscretizationProperties.currTime;
		const double dt = W.getPassport().timeDiscretizationProperties.dt;

		//wcm0 = W.getNonConstMechanics(bou).AngularVelocityOfAirfoil(currT);
		//wcmOld0 = W.getNonConstMechanics(bou).AngularVelocityOfAirfoil(currT - dt);
		//dwc0 = wcm0 - wcmOld0;

		dwcm = W.getNonConstMechanics(bou).AngularVelocityOfAirfoil(currT) - W.getNonConstMechanics(bou).AngularVelocityOfAirfoil(currT - dt);
		lastRhs[bou] += 2.0 * dwcm * W.getAirfoil(bou).area * (W.getAirfoil(bou).inverse ? 1.0 : -1.0);

		/*
		if (bou == 0)
		{
			lastRhs[bou] -= (2.0 * PI * 0.5) * dwc0 * 0.5;
		}

		if (bou == 1)
		{
			wcm1 = W.getNonConstMechanics(bou).AngularVelocityOfAirfoil(currT);
			wcmOld1 = W.getNonConstMechanics(bou).AngularVelocityOfAirfoil(currT - dt);
			dwc1 = wcm1 - wcmOld1;
			//std::cout << "dwc0,1 = " << dwc0 << " " << dwc1 << std::endl;
			lastRhs[bou] += (2.0 * PI * 1.0) * (dwc1) * 1.0;
		}
		*/
		////////////////////////////////////////////////////////////////

		//размазываем правую часть		
		for (size_t i = 0; i < nVars; ++i)
			rhs(i + currentRow) = locRhs(i);

		rhs(currentRow + nVars) = lastRhs[bou];

		currentRow += nVars + 1;
	}// for bou
}


void VelocityBiotSavart::CalcDiffVeloI1I2ToWakeFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	CalcDiffVeloI1I2ToSetOfPointsFromWake(pointsDb, domainRadius, vorticesDb, I1, I2);
}//CalcDiffVeloI1I2ToWakeFromWake(...)

void VelocityBiotSavart::CalcDiffVeloI1I2ToWakeFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	CalcDiffVeloI1I2ToSetOfPointsFromSheets(pointsDb, domainRadius, bnd, I1, I2);
}//CalcDiffVeloI1I2ToWakeFromSheets(...)