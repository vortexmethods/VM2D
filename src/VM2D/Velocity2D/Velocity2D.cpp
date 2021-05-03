/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Velocity2D.cpp                                                   |
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
\brief Файл кода с описанием класса Velocity
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
*/ 

#include "Velocity2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "StreamParser.h"
#include "Tree2D.h"
#include "Wake2D.h"
#include "WakeDataBase2D.h"
#include "World2D.h"
#include "Parallel.h"
#include "Passport2D.h"

#include "Velocity2DBarnesHut.h"
using namespace VM2D;

//Вычисление диффузионных скоростей вихрей и виртуальных вихрей в вихревом следе
void Velocity::CalcDiffVeloI1I2()
{	
	// !!! пелена на пелену
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))	
	if (W.getPassport().numericalSchemes.velocityComputation.second == 0)
		GPUCalcDiffVeloI1I2ToSetOfPointsFromWake(W.getWake(), wakeVortexesParams.epsastWake, W.getWake(), wakeVortexesParams.I1, wakeVortexesParams.I2);
	else
		CalcDiffVeloI1I2ToWakeFromWake(W.getWake(), wakeVortexesParams.epsastWake, W.getWake(), wakeVortexesParams.I1, wakeVortexesParams.I2);
#else
	CalcDiffVeloI1I2ToWakeFromWake(W.getWake(), wakeVortexesParams.epsastWake, W.getWake(), wakeVortexesParams.I1, wakeVortexesParams.I2);	
#endif
	
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		//не нужно, т.к. сделано выше перед началом вычисления скоростей
		//W.getNonConstBoundary(bou).virtualWake.WakeSynchronize();
		
		
	//виртуальные на границе на след
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))				
		if (W.getPassport().numericalSchemes.velocityComputation.second == 0)
			GPUCalcDiffVeloI1I2ToSetOfPointsFromSheets(W.getWake(), wakeVortexesParams.epsastWake, W.getBoundary(bou), wakeVortexesParams.I1, wakeVortexesParams.I2);
		else
			CalcDiffVeloI1I2ToWakeFromSheets(W.getWake(), wakeVortexesParams.epsastWake, W.getBoundary(bou), wakeVortexesParams.I1, wakeVortexesParams.I2);
#else			
		CalcDiffVeloI1I2ToWakeFromSheets(W.getWake(), wakeVortexesParams.epsastWake, W.getBoundary(bou), wakeVortexesParams.I1, wakeVortexesParams.I2);
#endif	

	// !!! след на виртуальные		
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))					
		GPUCalcDiffVeloI1I2ToSetOfPointsFromWake(W.getBoundary(bou).virtualWake, virtualVortexesParams[bou].epsastWake, W.getWake(), virtualVortexesParams[bou].I1, virtualVortexesParams[bou].I2);
#else
		CalcDiffVeloI1I2ToSetOfPointsFromWake(W.getBoundary(bou).virtualWake, virtualVortexesParams[bou].epsastWake, W.getWake(), virtualVortexesParams[bou].I1, virtualVortexesParams[bou].I2);
#endif
		for (size_t targetBou = 0; targetBou < W.getNumberOfBoundary(); ++targetBou)
		{
		// виртуальные на виртуальные
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))				
			GPUCalcDiffVeloI1I2ToSetOfPointsFromSheets(W.getBoundary(targetBou).virtualWake, virtualVortexesParams[targetBou].epsastWake, W.getBoundary(bou), virtualVortexesParams[targetBou].I1, virtualVortexesParams[targetBou].I2);
#else			
			CalcDiffVeloI1I2ToSetOfPointsFromSheets(W.getBoundary(targetBou).virtualWake, virtualVortexesParams[targetBou].epsastWake, W.getBoundary(bou), virtualVortexesParams[targetBou].I1, virtualVortexesParams[targetBou].I2);
#endif	
		}	
		
	} //for bou

	Point2D I2;
	double I1;
	
	for (size_t vt = 0; vt < wakeVortexesParams.diffVelo.size(); ++vt)
	{
		I2 = wakeVortexesParams.I2[vt];
		I1 = wakeVortexesParams.I1[vt];
		if (fabs(I1) < 1.e-8)
			wakeVortexesParams.diffVelo[vt] = { 0.0, 0.0 };
		else
			wakeVortexesParams.diffVelo[vt] = I2 * (1.0 / (I1 * std::max(wakeVortexesParams.epsastWake[vt], W.getPassport().wakeDiscretizationProperties.getMinEpsAst())));
	}	
	
	for (size_t targetBou = 0; targetBou < W.getNumberOfBoundary(); ++targetBou)
	for (size_t vt = 0; vt < virtualVortexesParams[targetBou].diffVelo.size(); ++vt)
	{
		I2 = virtualVortexesParams[targetBou].I2[vt];
		I1 = virtualVortexesParams[targetBou].I1[vt];

		if (fabs(I1) < 1.e-8)
			virtualVortexesParams[targetBou].diffVelo[vt] = { 0.0, 0.0 };
		else
			virtualVortexesParams[targetBou].diffVelo[vt] = I2 * (1.0 / (I1 * std::max(virtualVortexesParams[targetBou].epsastWake[vt], W.getPassport().wakeDiscretizationProperties.getMinEpsAst())));
	}
}//CalcDiffVeloI1I2()


void Velocity::CalcDiffVeloI0I3()
{
	for (size_t afl = 0; afl < W.getNumberOfAirfoil(); ++afl)
	{
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I0I3))		
		if (W.getPassport().numericalSchemes.velocityComputation.second == 0)
			W.getNonConstAirfoil(afl).GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(W.getWake(), wakeVortexesParams.epsastWake, wakeVortexesParams.I0, wakeVortexesParams.I3);
		else
			W.getNonConstAirfoil(afl).GetDiffVelocityI0I3ToWakeAndViscousStresses(W.getWake(), wakeVortexesParams.epsastWake, wakeVortexesParams.I0, wakeVortexesParams.I3);
#else
		W.getNonConstAirfoil(afl).GetDiffVelocityI0I3ToWakeAndViscousStresses(W.getWake(), wakeVortexesParams.epsastWake, wakeVortexesParams.I0, wakeVortexesParams.I3);		
#endif
	}	

	//Порядок циклов именно такой, т.к. CUDA оптимизирует вызов ядер, 
	//и на один "слой" виртуальных вихрей считается влияние сразу от всех профилей
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		for (size_t afl = 0; afl < W.getNumberOfAirfoil(); ++afl)
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I0I3))
			W.getNonConstAirfoil(afl).GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(W.getBoundary(bou).virtualWake, virtualVortexesParams[bou].epsastWake, virtualVortexesParams[bou].I0, virtualVortexesParams[bou].I3);
#else
			W.getNonConstAirfoil(afl).GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(W.getBoundary(bou).virtualWake, virtualVortexesParams[bou].epsastWake, virtualVortexesParams[bou].I0, virtualVortexesParams[bou].I3);
#endif			
	}
	//влияние поверхности
	Point2D I3;
	double I0;

	double domrad = 0.0;

	for (size_t vt = 0; vt < wakeVortexesParams.diffVelo.size(); ++vt)
	{
		domrad = std::max(wakeVortexesParams.epsastWake[vt], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());

		wakeVortexesParams.I0[vt] *= domrad;
		wakeVortexesParams.I0[vt] += DPI * sqr(domrad);

		I3 = wakeVortexesParams.I3[vt];
		I0 = wakeVortexesParams.I0[vt];

		if (fabs(I0) > 1.e-8)
			wakeVortexesParams.diffVelo[vt] += I3 * (1.0 / I0);
	}

	for (size_t targetBou = 0; targetBou < W.getNumberOfBoundary(); ++targetBou)
		for (size_t vt = 0; vt < virtualVortexesParams[targetBou].diffVelo.size(); ++vt)
		{
			domrad = std::max(virtualVortexesParams[targetBou].epsastWake[vt], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());

			virtualVortexesParams[targetBou].I0[vt] *= domrad;
			virtualVortexesParams[targetBou].I0[vt] += DPI * sqr(domrad);

			I3 = virtualVortexesParams[targetBou].I3[vt];
			I0 = virtualVortexesParams[targetBou].I0[vt];

			if (fabs(I0) > 1.e-8)
				virtualVortexesParams[targetBou].diffVelo[vt] += I3 * (1.0 / I0);
		}
}//CalcDiffVeloI0I3()

void Velocity::LimitDiffVelo(std::vector<Point2D>& diffVel)
{
	for (size_t i = 0; i < diffVel.size(); ++i)
	{
		Point2D& diffV = diffVel[i];

		diffV *= W.getPassport().physicalProperties.nu;

		if (diffV.length() > 1.5 * W.getPassport().physicalProperties.vRef)
			diffV.normalize(1.5 * W.getPassport().physicalProperties.vRef);
	}
}

// Вычисление диффузионных скоростей
void Velocity::CalcDiffVelo()
{
	W.getTimestat().timeCalcVortexDiffVelo.first += omp_get_wtime();
	if (W.getPassport().physicalProperties.nu > 0.0)
	{		
		CalcDiffVeloI1I2();		

		CalcDiffVeloI0I3();
				
		//контроль застрелов диффузионной скорости
		if (W.getParallel().myidWork == 0)
		{
			LimitDiffVelo(wakeVortexesParams.diffVelo);

			///omp
			for (size_t bou = 0; bou < virtualVortexesParams.size(); ++bou)
				LimitDiffVelo(virtualVortexesParams[bou].diffVelo);
		}
		
		for (size_t afl = 0; afl < W.getNumberOfAirfoil(); ++afl)
			for (size_t i = 0; i < W.getAirfoil(afl).viscousStress.size(); ++i)
				W.getNonConstAirfoil(afl).viscousStress[i] *= W.getPassport().physicalProperties.nu;

		SaveVisStress();

	}
	W.getTimestat().timeCalcVortexDiffVelo.second += omp_get_wtime();
}// CalcDiffVelo()


//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void Velocity::CalcDiffVeloI1I2ToSetOfPointsFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	double tCPUSTART, tCPUEND;

	tCPUSTART = omp_get_wtime();

	std::vector<double> selfI1;
	std::vector<Point2D> selfI2;

	const int& id = W.getParallel().myidWork;

	VMlib::parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D/*, VM2D::MyAlloc<VMlib::Vortex2D>*/>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	std::vector<double> locI1(par.myLen, 0.0);
	std::vector<Point2D> locI2(par.myLen, { 0.0, 0.0 });

#pragma warning (push)
#pragma warning (disable: 4101)
	//Локальные переменные для цикла
	Point2D Rij;
	double rij, expr;
	double diffRadius, domRad;
	double left;
	double right;
	double posJx;
#pragma warning (pop)

#pragma omp parallel for default(none) shared(locI1, locI2, domainRadius, locPoints, vorticesDb, par) private(Rij, rij, expr, diffRadius, domRad, left, right, posJx)
	for (int i = 0; i < par.myLen; ++i)
	{
		const Vortex2D& vtxI = locPoints[i];

		domRad = std::max(domainRadius[i + par.myDisp], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());

		/// \todo Понять природу магической константы 8.0 и синхронизировать с GPU
		diffRadius = 8.0 * domRad;

		left = vtxI.r()[0] - diffRadius;
		right = vtxI.r()[0] + diffRadius;

		for (size_t j = 0; j < vorticesDb.vtx.size(); ++j)
		{
			const Vortex2D& vtxJ = vorticesDb.vtx[j];
			posJx = vtxJ.r()[0];

			if ((left < posJx) && (posJx < right))
			{
				Rij = vtxI.r() - vtxJ.r();
				rij = Rij.length();
				if (rij < diffRadius && rij > 1.e-10)
				{
					expr = exp(-rij / domRad);
					locI2[i] += (vtxJ.g()* expr / rij) * Rij;
					locI1[i] += vtxJ.g()*expr;
				}
			}//if (rij>1e-6)
		}//for j
	} // for r

	if (id == 0)
	{
		selfI1.resize(pointsDb.vtx.size(), 0.0);
		selfI2.resize(pointsDb.vtx.size(), { 0.0, 0.0 });
	}

	MPI_Gatherv(locI1.data(), par.myLen, MPI_DOUBLE, selfI1.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

	MPI_Gatherv(locI2.data(), par.myLen, Point2D::mpiPoint2D, selfI2.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);

	if (id == 0)
		for (size_t i = 0; i < I1.size(); ++i)
		{
			I1[i] += selfI1[i];
			I2[i] += selfI2[i];
		}

	tCPUEND = omp_get_wtime();
	//W.getInfo('t') << "DIFF_CPU: " << tCPUEND - tCPUSTART << std::endl;
}//CalcDiffVeloI1I2ToSetOfPointsFromWake(...)


//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void Velocity::CalcDiffVeloI1I2ToSetOfPointsFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	double tCPUSTART, tCPUEND;

	tCPUSTART = omp_get_wtime();

	std::vector<double> selfI1;
	std::vector<Point2D> selfI2;

	const int& id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size());


	//синхронизация свободного вихревого слоя
	bnd.sheets.FreeSheetSynchronize();

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D/*, VM2D::MyAlloc<VMlib::Vortex2D>*/>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	std::vector<double> locI1(par.myLen, 0.0);
	std::vector<Point2D> locI2(par.myLen, { 0.0, 0.0 });

#pragma warning (push)
#pragma warning (disable: 4101)
	//Локальные переменные для цикла
	Point2D Rij;
	double rij, expr;
	double diffRadius;
	double left;
	double right;
	double posJx;
	double domRad;
#pragma warning (pop)


#pragma omp parallel for default(none) shared(locI1, locI2, domainRadius, locPoints, bnd, par, std::cout) private(Rij, rij, expr, domRad, diffRadius, left, right, posJx)
	for (int i = 0; i < par.myLen; ++i)
	{
		const Vortex2D& vtxI = locPoints[i];

		domRad = std::max(domainRadius[i + par.myDisp], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());

		/// \todo Понять природу магической константы 8.0 и синхронизировать с GPU
		diffRadius = 8.0 * domRad;

		left = vtxI.r()[0] - diffRadius;
		right = vtxI.r()[0] + diffRadius;

		for (size_t j = 0; j < bnd.afl.getNumberOfPanels(); ++j)
		{
			/// \todo Сделать переменной и синхронизировать с GPU
			const int nQuadPt = 3;

			/// \todo Учитываем пока только нулевой момент решения
			const double ptG = bnd.sheets.freeVortexSheet(j, 0) * bnd.afl.len[j] / nQuadPt;

			for (int q = 0; q < nQuadPt; ++q)
			{
				const Point2D& ptJ = bnd.afl.getR(j) + bnd.afl.tau[j] * (q + 0.5) * bnd.afl.len[j] * (1.0 / nQuadPt);  // vorticesDb.vtx[j];
				posJx = ptJ[0];

				if ((left < posJx) && (posJx < right))
				{
					Rij = vtxI.r() - ptJ;
					rij = Rij.length();
					if (rij < diffRadius && rij > 1.e-10)
					{
						expr = exp(-rij / domRad);
						locI2[i] += (ptG * expr / rij) * Rij;
						locI1[i] += ptG * expr;
					}
				}//if (rij>1e-6)
			}
		}//for j
	} // for r

	if (id == 0)
	{
		selfI1.resize(pointsDb.vtx.size(), 0.0);
		selfI2.resize(pointsDb.vtx.size(), { 0.0, 0.0 });
	}

	MPI_Gatherv(locI1.data(), par.myLen, MPI_DOUBLE, selfI1.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);
	MPI_Gatherv(locI2.data(), par.myLen, Point2D::mpiPoint2D, selfI2.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);


	if (id == 0)
		for (size_t i = 0; i < I1.size(); ++i)
		{
			I1[i] += selfI1[i];
			I2[i] += selfI2[i];
		}

	tCPUEND = omp_get_wtime();
	//W.getInfo('t') << "DIFF_CPU: " << tCPUEND - tCPUSTART << std::endl;
}//CalcDiffVeloI1I2ToSetOfPointsFromSheets(...)



#if defined (USE_CUDA)
//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void Velocity::GPUCalcDiffVeloI1I2ToSetOfPointsFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh)
{
	if ( (&pointsDb == &W.getWake()) || (&pointsDb == &W.getBoundary(0).virtualWake) )
	{
		size_t npt = pointsDb.vtx.size();

		if ((W.getNumberOfBoundary() > 0) && (&pointsDb == &W.getBoundary(0).virtualWake))
		{
			for (size_t q = 1; q < W.getNumberOfBoundary(); ++q)
				npt += W.getBoundary(q).virtualWake.vtx.size();
		}

		double*& dev_ptr_pt = pointsDb.devVtxPtr;
		double*& dev_ptr_dr = pointsDb.devRadPtr;

		const size_t nvt = vorticesDb.vtx.size();
		double*& dev_ptr_vt = vorticesDb.devVtxPtr;

		std::vector<double> loci1(npt);
		std::vector<Point2D> loci2(npt);
				
		double*& dev_ptr_i1 = pointsDb.devI1Ptr;
		double*& dev_ptr_i2 = pointsDb.devI2Ptr;
		double minRad = W.getPassport().wakeDiscretizationProperties.getMinEpsAst();

		const int& id = W.getParallel().myidWork;
		VMlib::parProp par = W.getParallel().SplitMPI(npt, true);

		if ((nvt > 0) && (npt > 0))
		{
			//СЕТКА
			if (!useMesh)
				cuCalculateDiffVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr, minRad);
			else
				cuCalculateDiffVeloWakeMesh(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, W.getWake().devMeshPtr, W.getPassport().wakeDiscretizationProperties.epscol, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);

			W.getCuda().CopyMemFromDev<double, 2>(par.myLen, dev_ptr_i2, (double*)&loci2[0], 10);
			W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_i1, &loci1[0], 11);

			std::vector<Point2D> newI2;
			std::vector<double> newI1;
			if (id == 0)
			{
				newI2.resize(npt); // I2.size());
				newI1.resize(npt); // I1.size());
			}

			MPI_Gatherv(loci2.data(), par.myLen, Point2D::mpiPoint2D, newI2.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
			MPI_Gatherv(loci1.data(), par.myLen, MPI_DOUBLE, newI1.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

			if (&pointsDb == &W.getWake())
			{
				if (id == 0)
					for (size_t q = 0; q < I2.size(); ++q)
					{
						I1[q] += newI1[q];
						I2[q] += newI2[q];
					}
			}

			if ((W.getNumberOfBoundary() > 0) && (&pointsDb == &W.getBoundary(0).virtualWake))
			{
				size_t curGlobPnl = 0;
				if (id == 0)
					for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s) //W.getNumberOfAirfoil()
					{
						size_t nv = W.getBoundary(s).virtualWake.vtx.size();
						for (size_t q = 0; q < nv; ++q)
						{
							W.getNonConstVelocity().virtualVortexesParams[s].I1[q] += newI1[curGlobPnl + q];
							W.getNonConstVelocity().virtualVortexesParams[s].I2[q] += newI2[curGlobPnl + q];
						}
						curGlobPnl += nv;
					}
			}
		}	
	}	
}//GPUCalcDiffVeloI1I2ToSetOfPointsFromWake(...)


//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void Velocity::GPUCalcDiffVeloI1I2ToSetOfPointsFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bou, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh)
{
	if ((&pointsDb == &W.getWake()) || (&pointsDb == &W.getBoundary(0).virtualWake))
	{	
		if (bou.afl.numberInPassport == 0)
		{
			size_t npt = pointsDb.vtx.size();
			double*& dev_ptr_pt = pointsDb.devVtxPtr;
			double*& dev_ptr_dr = pointsDb.devRadPtr;

			if ((W.getNumberOfBoundary() > 0) && (&pointsDb == &W.getBoundary(0).virtualWake))
			{
				for (size_t q = 1; q < W.getNumberOfBoundary(); ++q)
					npt += W.getBoundary(q).virtualWake.vtx.size();
			}

			size_t npnl = bou.afl.getNumberOfPanels(); //vorticesDb.vtx.size();
			double*& dev_ptr_r = bou.afl.devRPtr;
			double*& dev_ptr_freeVortexSheet = bou.afl.devFreeVortexSheetPtr;

			for (size_t q = 1; q < W.getNumberOfBoundary(); ++q)
				npnl += W.getBoundary(q).afl.getNumberOfPanels();

			std::vector<double> loci1(npt);
			std::vector<Point2D> loci2(npt);

			double*& dev_ptr_i1 = pointsDb.devI1Ptr;
			double*& dev_ptr_i2 = pointsDb.devI2Ptr;
			double minRad = W.getPassport().wakeDiscretizationProperties.getMinEpsAst();

			const int& id = W.getParallel().myidWork;
			VMlib::parProp par = W.getParallel().SplitMPI(npt, true);

			if ((npnl > 0) && (npt > 0))
			{
				/// \todo Реализовать версию с сеткой
				//if (!useMesh)
				cuCalculateDiffVeloWakeFromPanels(par.myDisp, par.myLen, dev_ptr_pt, npnl, dev_ptr_r, dev_ptr_freeVortexSheet, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr, minRad);
				//else
				//	cuCalculateDiffVeloWakeMesh(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, W.getWake().devMeshPtr, W.getPassport().wakeDiscretizationProperties.epscol, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);

				W.getCuda().CopyMemFromDev<double, 2>(par.myLen, dev_ptr_i2, (double*)&loci2[0], 12);
				W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_i1, &loci1[0], 13);

				std::vector<Point2D> newI2;
				std::vector<double> newI1;
				if (id == 0)
				{
					newI2.resize(npt);
					newI1.resize(npt);
				}

				MPI_Gatherv(loci2.data(), par.myLen, Point2D::mpiPoint2D, newI2.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
				MPI_Gatherv(loci1.data(), par.myLen, MPI_DOUBLE, newI1.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

				if (&pointsDb == &W.getWake())
				{
					if (id == 0)
						for (size_t q = 0; q < I2.size(); ++q)
						{
							I1[q] += newI1[q];
							I2[q] += newI2[q];
						}
				}

				if ((W.getNumberOfBoundary() > 0) && (&pointsDb == &W.getBoundary(0).virtualWake))
				{
					size_t curGlobPnl = 0;
					if (id == 0)
						for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s) //W.getNumberOfAirfoil()
						{
							size_t nv = W.getBoundary(s).virtualWake.vtx.size();
							for (size_t q = 0; q < nv; ++q)
							{
								W.getNonConstVelocity().virtualVortexesParams[s].I1[q] += newI1[curGlobPnl + q];
								W.getNonConstVelocity().virtualVortexesParams[s].I2[q] += newI2[curGlobPnl + q];
							}
							curGlobPnl += nv;
						}
				}
			}
		}
	}
	
}//GPUCalcDiffVeloI1I2ToSetOfPointsFromSheets(...)
#endif

// Очистка старых массивов под хранение скоростей, выделение новой памяти и обнуление
void Velocity::ResizeAndZero()
{
	if (W.getParallel().myidWork == 0)
	{
		W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
		wakeVortexesParams.convVelo.clear();
		wakeVortexesParams.convVelo.resize(W.getWake().vtx.size(), { 0.0, 0.0 });
		W.getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();

		W.getTimestat().timeCalcVortexDiffVelo.first += omp_get_wtime();
		wakeVortexesParams.I0.clear();
		wakeVortexesParams.I0.resize(W.getWake().vtx.size(), 0.0);

		wakeVortexesParams.I1.clear();
		wakeVortexesParams.I1.resize(W.getWake().vtx.size(), 0.0);

		wakeVortexesParams.I2.clear();
		wakeVortexesParams.I2.resize(W.getWake().vtx.size(), { 0.0, 0.0 });

		wakeVortexesParams.I3.clear();
		wakeVortexesParams.I3.resize(W.getWake().vtx.size(), { 0.0, 0.0 });

		wakeVortexesParams.diffVelo.clear();
		wakeVortexesParams.diffVelo.resize(W.getWake().vtx.size(), { 0.0, 0.0 });

		wakeVortexesParams.epsastWake.clear();
		wakeVortexesParams.epsastWake.resize(W.getWake().vtx.size(), 0.0);
		W.getTimestat().timeCalcVortexDiffVelo.second += omp_get_wtime();

		//Создаем массивы под виртуальные вихри
		W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
		virtualVortexesParams.clear();
		virtualVortexesParams.resize(W.getNumberOfBoundary());
		W.getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();

		for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
		{
			W.getTimestat().timeCalcVortexConvVelo.first += omp_get_wtime();
			virtualVortexesParams[bou].convVelo.clear();
			virtualVortexesParams[bou].convVelo.resize(W.getBoundary(bou).virtualWake.vtx.size(), { 0.0, 0.0 });
			W.getTimestat().timeCalcVortexConvVelo.second += omp_get_wtime();

			W.getTimestat().timeCalcVortexDiffVelo.first += omp_get_wtime();
			virtualVortexesParams[bou].I0.clear();
			virtualVortexesParams[bou].I0.resize(W.getBoundary(bou).virtualWake.vtx.size(), 0.0);

			virtualVortexesParams[bou].I1.clear();
			virtualVortexesParams[bou].I1.resize(W.getBoundary(bou).virtualWake.vtx.size(), 0.0);

			virtualVortexesParams[bou].I2.clear();
			virtualVortexesParams[bou].I2.resize(W.getBoundary(bou).virtualWake.vtx.size(), { 0.0, 0.0 });

			virtualVortexesParams[bou].I3.clear();
			virtualVortexesParams[bou].I3.resize(W.getBoundary(bou).virtualWake.vtx.size(), { 0.0, 0.0 });

			virtualVortexesParams[bou].diffVelo.clear();
			virtualVortexesParams[bou].diffVelo.resize(W.getBoundary(bou).virtualWake.vtx.size(), { 0.0, 0.0 });

			virtualVortexesParams[bou].epsastWake.clear();
			virtualVortexesParams[bou].epsastWake.resize(W.getBoundary(bou).virtualWake.vtx.size(), 0.0);

			W.getTimestat().timeCalcVortexDiffVelo.second += omp_get_wtime();
		}
	}
}//ResizeAndZero()


void Velocity::SaveVisStress()
{
	if ((W.getPassport().timeDiscretizationProperties.saveVTK > 0) && (W.ifDivisible(W.getPassport().timeDiscretizationProperties.saveVisStress)))
		//		if ((W.getPassport().timeDiscretizationProperties.saveVTK > 0) && (W.ifDivisible(10)) && (W.getNumberOfAirfoil() > 0))
	{
		for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
		{
			if (q == 0)
				VMlib::CreateDirectory(W.getPassport().dir, "visStress");

			std::stringstream ss;
			ss << "VisStress_" << q << "-";
			std::string fname = VMlib::fileNameStep(ss.str(), W.getPassport().timeDiscretizationProperties.nameLength, W.getCurrentStep(), "txt");
			std::ofstream outfile;
			outfile.open(W.getPassport().dir + "visStress/" + fname);

			outfile << W.getAirfoil(q).viscousStress.size() << std::endl; //Сохранение числа вихрей в пелене

			for (size_t i = 0; i < W.getAirfoil(q).viscousStress.size(); ++i)
			{
				const Point2D& r = 0.5 * (W.getAirfoil(q).getR(i + 1) + W.getAirfoil(q).getR(i));
				double gi = W.getAirfoil(q).viscousStress[i];
				outfile << static_cast<int>(i) << " " << r[0] << " " << r[1] << " " << gi << std::endl;
			}//for i	
			outfile.close();
		}
	}
}//SaveVisStress()