/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.9    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/07/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.9   
\date 22 июля 2020 г.
*/

#include "Velocity2DBiotSavart.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Tree2D.h"
#include "Wake2D.h"
#include "World2D.h"

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
	//double tt1 = omp_get_wtime();
	
	W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONV_TOWAKE))	
	GPUCalcConvVeloToSetOfPointsFromWake(W.getWake(), wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake, true, true);	
#else
	CalcConvVeloToSetOfPointsFromWake(W.getWake(), wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake, true, true);
#endif
		
	//double tt2 = omp_get_wtime();

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
	if ((W.getPassport().timeDiscretizationProperties.saveVP > 0) && (!(W.getCurrentStep() % W.getPassport().timeDiscretizationProperties.saveVP)))
		CalcVeloToWakeVP();	
	W.getTimestat().timeVP.second += omp_get_wtime();	
}//CalcConvVelo()

//Вычисление конвективных скоростей вихрей в точках wakeVP
void VelocityBiotSavart::CalcVeloToWakeVP()
{
	std::vector<Point2D> velConvWake;
	std::vector<std::vector<Point2D>> velConvBou;

	int addWSize = (int)W.getMeasureVP().getWakeVP().vtx.size();
	MPI_Bcast(&addWSize, 1, MPI_INT, 0, W.getParallel().commWork);

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

	if (W.getParallel().myidWork == 0)
	{
		std::vector<Point2D>& velocityRef = W.getNonConstMeasureVP().getNonConstVelocity();
		velocityRef.assign(addWSize, W.getPassport().physicalProperties.V0());
		
		for (int i = 0; i < addWSize; ++i)
			velocityRef[i] += velConvWake[i];

		for (size_t bou = 0; bou < velConvBou.size(); ++bou)
			for (int j = 0; j < addWSize; ++j)
				velocityRef[j] += velConvBou[bou][j];
	}
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
	std::vector<Point2D> selfVelo;
	
	const int& id = W.getParallel().myidWork;

	VMlib::parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size(), true);

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	//рассылка точек вычисления скоростей и радиусов доменов
	MPI_Scatterv(const_cast<std::vector<Vortex2D/*, VM2D::MyAlloc<VMlib::Vortex2D>*/>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		         locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	double cft = IDPI;

	std::vector<Point2D> locConvVelo;
	if (calcVelo)
		locConvVelo.resize(par.myLen);

	std::vector<double> locDomRadius;
	if (calcRadius)
		locDomRadius.resize(par.myLen);

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
#pragma omp parallel for default(none) shared(locConvVelo, locDomRadius, locPoints, cft, par, calcVelo, calcRadius, eps2) private(tempVel, velI, dst2, dst2eps) schedule(dynamic, DYN_SCHEDULE)
		for (int i = 0; i < par.myLen; ++i)
		{
			double ee2[3] = { 10000.0, 10000.0, 10000.0 };

			velI.toZero();

			const Point2D& posI = locPoints[i].r();

			for (size_t j = 0; j < W.getWake().vtx.size(); ++j)
			{
				const Point2D& posJ = W.getWake().vtx[j].r();

				dst2 = dist2(posI, posJ);

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
			locConvVelo[i] = velI;

#ifndef TESTONLYVELO
			if (calcRadius)
				locDomRadius[i] = sqrt((ee2[0] + ee2[1] + ee2[2]) / 3.0);
#endif
		}
	}
	else if (calcRadius)
	{
#pragma omp parallel for default(none) shared(locConvVelo, locDomRadius, locPoints, cft, par, calcVelo, calcRadius) private(tempVel, velI, dst2) schedule(dynamic, DYN_SCHEDULE)
		for (int i = 0; i < par.myLen; ++i)
		{
			double ee2[3] = { 10000.0, 10000.0, 10000.0 };

			velI.toZero();

			const Point2D& posI = locPoints[i].r();

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
				}
			}

			locDomRadius[i] = sqrt((ee2[0] + ee2[1] + ee2[2]) / 3.0);
		}
	} //else


	if (calcVelo)
	{
		if (id == 0) 
			selfVelo.resize(pointsDb.vtx.size());
		MPI_Gatherv(locConvVelo.data(), par.myLen, Point2D::mpiPoint2D, selfVelo.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
	}

#ifndef TESTONLYVELO
	if (calcRadius)
	{
		//if (id == 0)
		//это нужно всем, т.к. ниже стоит Allgatherv
		domainRadius.resize(par.totalLen);

		MPI_Allgatherv(locDomRadius.data(), par.myLen, MPI_DOUBLE, domainRadius.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, W.getParallel().commWork);
	}
#endif

	if ( (calcVelo) && (id == 0) )
		for (size_t i = 0; i < velo.size(); ++i)
			velo[i] += selfVelo[i];
}//CalcConvVeloToSetOfPointsFromWake(...)



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

		std::vector<Point2D> locvel(npt);
		std::vector<double> locrad(npt);

		double*& dev_ptr_vel = pointsDb.devVelPtr;
		double*& dev_ptr_rad = pointsDb.devRadPtr;
		const double& eps2 = W.getPassport().wakeDiscretizationProperties.eps2;

		const int& id = W.getParallel().myidWork;
		VMlib::parProp par = W.getParallel().SplitMPI(npt, true);
		
		if (npt > 0)
		{
			//double t1 = omp_get_wtime();
			cuCalculateConvVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, nsr, dev_ptr_sr, nbou, dev_nVortices, dev_ptr_ptr_vtx, dev_ptr_vel, dev_ptr_rad, eps2, calcVelo, calcRadius);
			//double t2 = omp_get_wtime();
			//std::cout << "KERNEL_TIME = " << t2 - t1 << std::endl;
			
			if (calcVelo)
			{
				//double tt1 = omp_get_wtime();
				
				W.getCuda().CopyMemFromDev<double, 2>(par.myLen, dev_ptr_vel, (double*)&locvel[0], 20);
			
				//double tt2 = omp_get_wtime();
				//std::cout << "COPY_TIME = " << tt2 - tt1 << std::endl;
				
				//double tt3 = omp_get_wtime();
				std::vector<Point2D> newV;
				if (id == 0)
					newV.resize(Vel.size());

				MPI_Gatherv(locvel.data(), par.myLen, Point2D::mpiPoint2D, newV.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
				if (id == 0)
					for (size_t q = 0; q < npt; ++q)
						Vel[q] = newV[q];
				
				if (id == 0)
					if ((&pointsDb == &W.getWake()) || (&pointsDb == &W.getMeasureVP().getWakeVP()))
					{
						for (size_t q = 0; q < npt; ++q)
							velo[q] += Vel[q];
					}//if &pointsDb

				if (id == 0)
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
				W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_rad, &locrad[0], 21);
				Rad.resize(par.totalLen);
				MPI_Allgatherv(locrad.data(), par.myLen, MPI_DOUBLE, Rad.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, W.getParallel().commWork);

				cuCopyFixedArray(dev_ptr_rad, Rad.data(), sizeof(double) * Rad.size());
				
				if (id == 0)
					if ((&pointsDb == &W.getWake()) || (&pointsDb == &W.getMeasureVP().getWakeVP()))
					{
						for (size_t q = 0; q < Rad.size(); ++q)
							domainRadius[q] = Rad[q];
					}//if &pointsDb
				
				if (id == 0)
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
	int id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(np);

	std::vector<double> locVeloWake, locVeloWakeLin;
	locVeloWake.resize(par.myLen);

	size_t shDim = W.getBoundary(afl.numberInPassport).sheetDim;

	if(shDim != 1)
		locVeloWakeLin.resize(par.myLen);
	
	//локальные переменные для цикла	
	std::vector<double> velI(shDim, 0.0);

#pragma omp parallel for default(none) shared(locVeloWake, par, shDim, afl, locVeloWakeLin, IDPI) private(velI)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI.assign(shDim, 0.0);

		if (W.getWake().vtx.size() > 0)
		{
			//Учет влияния следа
			afl.GetInfluenceFromVorticesToPanel(par.myDisp + i, &(*W.getWake().vtx.begin()), std::distance(W.getWake().vtx.begin(), W.getWake().vtx.end()), velI);
		}

		if (W.getSource().vtx.size() > 0)
		{
			//Учет влияния источников
			afl.GetInfluenceFromSourcesToPanel(par.myDisp + i, &(*W.getSource().vtx.begin()), std::distance(W.getSource().vtx.begin(), W.getSource().vtx.end()), velI);
		}
		
		for (size_t j = 0; j < shDim; ++j)
			velI[j] *= IDPI / afl.len[par.myDisp + i];
		
		locVeloWake[i] = velI[0];

		if (shDim != 1)
			locVeloWakeLin[i] = velI[1];
	}

	if (id == 0)
		wakeRhs.resize( W.getBoundary(afl.numberInPassport).GetUnknownsSize() );

	MPI_Gatherv(locVeloWake.data(), par.myLen, MPI_DOUBLE, wakeRhs.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

	if (shDim != 1)
		MPI_Gatherv(locVeloWakeLin.data(), par.myLen, MPI_DOUBLE, wakeRhs.data() + np, par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);
}//GetWakeInfluenceToRhs(...)



#if defined(USE_CUDA)
//Генерация вектора влияния вихревого следа на профиль
void VelocityBiotSavart::GPUGetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeVelo) const
{
	const int& id = W.getParallel().myidWork;	

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
		std::vector<double> locrhs(nTotPan);

		VMlib::parProp par = W.getParallel().SplitMPI(nTotPan, true);

		if ((nvt > 0) || (nsr > 0))
		{
			cuCalculateRhs(par.myDisp, par.myLen, nTotPan, dev_ptr_pt, nvt, dev_ptr_vt, nsr, dev_ptr_sr, dev_ptr_rhs);

			W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_rhs, (double*)&locrhs[0], 22);
				
			std::vector<double> newRhs;
			if (id == 0)
				newRhs.resize(nTotPan);

			MPI_Gatherv(locrhs.data(), par.myLen, MPI_DOUBLE, newRhs.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);			

			if (id == 0)
			{

				//for (int qq = 0; qq < (int)newRhs.size(); ++qq)
				//	std::cout << "qq = " << qq << ", " << newRhs[qq] << std::endl;

				size_t curGlobPnl = 0;
				for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s)
				{
					std::vector<double>& tmpRhs = W.getNonConstAirfoil(s).tmpRhs;
					const size_t& np = W.getAirfoil(s).getNumberOfPanels();
					tmpRhs.resize(0);
					tmpRhs.insert(tmpRhs.end(), newRhs.begin() + curGlobPnl, newRhs.begin() + curGlobPnl + np);
					curGlobPnl += np;
				}
			}
		}
	}


	if (id == 0)
	{
		if ((nvt > 0) || (nsr > 0))
			wakeVelo = std::move(afl.tmpRhs);
		else
			wakeVelo.resize(afl.getNumberOfPanels(), 0.0);		
	}
	
}//GPUGetWakeInfluenceToRhs(...)
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

		if (W.getParallel().myidWork == 0)
		{
			nVars = W.getBoundary(bou).GetUnknownsSize();

			locRhs.resize(nVars);
		}


		std::vector<double> wakeRhs;

			
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_RHS))
		GPUGetWakeInfluenceToRhs(afl, wakeRhs);
#else
		GetWakeInfluenceToRhs(afl, wakeRhs);
#endif				

		std::vector<double> vInfRhs;
		afl.GetInfluenceFromVInfToPanel(vInfRhs);
		
		if (W.getParallel().myidWork == 0)
		{
#pragma omp parallel for \
	default(none) \
	shared(locRhs, afl, bou, wakeRhs, vInfRhs, np) \
	schedule(dynamic, DYN_SCHEDULE)
			for (int i = 0; i < (int)afl.getNumberOfPanels(); ++i)
			{			
				/// \todo ВАЖНО!!!
				//Здесь учет скорости профиля и присоединенных слоев источников сделан только для прямолинейных панелей
				locRhs(i) = - vInfRhs[i] - wakeRhs[i] + 0.25 * ((afl.getV(i) + afl.getV(i + 1)) & afl.tau[i]); //0.25 * (afl.getV(i) + afl.getV(i + 1))*afl.tau[i] - прямолинейные
				if(W.getBoundary(bou).sheetDim > 1)
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
			
			//размазываем правую часть		
			for (size_t i = 0; i < nVars; ++i)
				rhs(i + currentRow) = locRhs(i);
			
			rhs(currentRow + nVars) = lastRhs[bou];

			currentRow += nVars + 1;

		}// if (parallel.myidWork == 0)

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