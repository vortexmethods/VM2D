/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
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
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/ 

#include "Velocity2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "StreamParser.h"
#include "Wake2D.h"
#include "World2D.h"

#include "BarnesHut.h"

using namespace VM2D;

//Вычисление диффузионных скоростей вихрей и виртуальных вихрей в вихревом следе
void Velocity::CalcDiffVeloI1I2()
{	
	// !!! пелена на пелену
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))	
	if (W.getPassport().numericalSchemes.velocityComputation.second == 0)
		GPUCalcDiffVeloI1I2ToSetOfPointsFromWake(W.getWake(), wakeVortexesParams.epsastWake, W.getWake(), wakeVortexesParams.I1, wakeVortexesParams.I2);		
	else
		GPUDiffVeloFAST(wakeVortexesParams.epsastWake, wakeVortexesParams.I1, wakeVortexesParams.I2);		
#else
	CalcDiffVeloI1I2ToWakeFromWake(W.getWake(), wakeVortexesParams.epsastWake, W.getWake(), wakeVortexesParams.I1, wakeVortexesParams.I2);	
#endif

	 
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		//не нужно, т.к. сделано выше перед началом вычисления скоростей
		//W.getNonConstBoundary(bou).virtualWake.WakeSynchronize();
		
		
	//виртуальные на границе на след
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))				
		GPUCalcDiffVeloI1I2ToSetOfPointsFromSheets(W.getWake(), wakeVortexesParams.epsastWake, W.getBoundary(bou), wakeVortexesParams.I1, wakeVortexesParams.I2);
#else			
		CalcDiffVeloI1I2ToWakeFromSheets(W.getWake(), wakeVortexesParams.epsastWake, W.getBoundary(bou), wakeVortexesParams.I1, wakeVortexesParams.I2);
#endif			
	} //for bou
	
	
	double tDivstart, tDivfinish;

	Point2D I2;
	double I1;
	
	tDivstart = omp_get_wtime();

#pragma omp parallel for private(I1, I2)
	for (int vt = 0; vt < (int)wakeVortexesParams.diffVelo.size(); ++vt)
	{
		I2 = wakeVortexesParams.I2[vt];
		I1 = wakeVortexesParams.I1[vt];

		if (fabs(I1) < 1.e-8)
			wakeVortexesParams.diffVelo[vt] = { 0.0, 0.0 };
		else
			wakeVortexesParams.diffVelo[vt] = I2 * (1.0 / (I1 * std::max(wakeVortexesParams.epsastWake[vt], W.getPassport().wakeDiscretizationProperties.getMinEpsAst())));
	}	

	tDivfinish = omp_get_wtime();

}//CalcDiffVeloI1I2()


void Velocity::CalcDiffVeloI0I3()
{
	for (size_t afl = 0; afl < W.getNumberOfAirfoil(); ++afl)
	{
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I0I3))		
		W.getNonConstAirfoil(afl).GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(wakeVortexesParams.epsastWake, wakeVortexesParams.I0, wakeVortexesParams.I3);
#else
		W.getNonConstAirfoil(afl).GetDiffVelocityI0I3ToWakeAndViscousStresses(W.getWake(), wakeVortexesParams.epsastWake, wakeVortexesParams.I0, wakeVortexesParams.I3);		
#endif
	}	


	//влияние поверхности
	Point2D I3;
	double I0;

	double domrad = 0.0;

#pragma omp parallel for private(I0, I3, domrad)
	for (int vt = 0; vt < (int)wakeVortexesParams.diffVelo.size(); ++vt)
	{
		domrad = std::max(wakeVortexesParams.epsastWake[vt], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());

		wakeVortexesParams.I0[vt] *= domrad;
		wakeVortexesParams.I0[vt] += DPI * sqr(domrad);
			

		I3 = wakeVortexesParams.I3[vt];
		I0 = wakeVortexesParams.I0[vt];
		

		if (fabs(I0) > 1.e-8)
			wakeVortexesParams.diffVelo[vt] += I3 * (1.0 / I0);		
	}
}//CalcDiffVeloI0I3()


void Velocity::LimitDiffVelo(std::vector<Point2D>& diffVel)
{
	for (size_t i = 0; i < diffVel.size(); ++i)
	{
		Point2D& diffV = diffVel[i];		

		if (diffV.length() > 1.5 * W.getPassport().physicalProperties.vRef)
			diffV.normalize(1.5 * W.getPassport().physicalProperties.vRef);
	}
}

// Вычисление диффузионных скоростей
void Velocity::CalcDiffVelo()
{
	W.getTimers().start("DiffVel");

	double t12start, t12finish;
	double t03start, t03finish;
	double tOtherstart, tOtherfinish;

	if (W.getPassport().physicalProperties.nu > 0.0)
	{
		t12start = omp_get_wtime();
		CalcDiffVeloI1I2();
		t12finish = omp_get_wtime();

		t03start = omp_get_wtime();
		CalcDiffVeloI0I3();
		t03finish = omp_get_wtime();
  
		tOtherstart = omp_get_wtime();
		for (Point2D& diffV : wakeVortexesParams.diffVelo)
			diffV *= W.getPassport().physicalProperties.nu;

		//контроль застрелов диффузионной скорости		
		LimitDiffVelo(wakeVortexesParams.diffVelo);

		for (size_t afl = 0; afl < W.getNumberOfAirfoil(); ++afl)
			for (size_t i = 0; i < W.getAirfoil(afl).viscousStress.size(); ++i)
				W.getNonConstAirfoil(afl).viscousStress[i] *= W.getPassport().physicalProperties.nu;

		//SaveVisStress();
		

		//Заполнение структуры данных виртуальных вихрей
		size_t nVirtVortices = 0;
		for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
			nVirtVortices += W.getBoundary(bou).virtualWake.vtx.size();

		size_t counter = wakeVortexesParams.convVelo.size() - nVirtVortices;

		for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
			for (size_t v = 0; v < W.getBoundary(bou).virtualWake.vtx.size(); ++v)
			{			
				virtualVortexesParams[bou].diffVelo[v] = wakeVortexesParams.diffVelo[counter];
				virtualVortexesParams[bou].I0[v] = wakeVortexesParams.I0[counter];
				virtualVortexesParams[bou].I1[v] = wakeVortexesParams.I1[counter];
				virtualVortexesParams[bou].I2[v] = wakeVortexesParams.I2[counter];
				virtualVortexesParams[bou].I3[v] = wakeVortexesParams.I3[counter];				
				++counter;
			}

		tOtherfinish = omp_get_wtime();		
	}
	W.getTimers().stop("DiffVel");

}// CalcDiffVelo()


//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void Velocity::CalcDiffVeloI1I2ToSetOfPointsFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	double tCPUSTART, tCPUEND;

	tCPUSTART = omp_get_wtime();

	std::vector<double> selfI1(pointsDb.vtx.size(), 0.0);
	std::vector<Point2D> selfI2(pointsDb.vtx.size(), { 0.0, 0.0 });	

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

#pragma omp parallel for default(none) shared(selfI1, selfI2, domainRadius, vorticesDb, pointsDb) private(Rij, rij, expr, diffRadius, domRad, left, right, posJx)
	for (int i = 0; i < pointsDb.vtx.size(); ++i)
	{
		const Vortex2D& vtxI = pointsDb.vtx[i];

		domRad = std::max(domainRadius[i], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());

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
					selfI2[i] += (vtxJ.g()* expr / rij) * Rij;
					selfI1[i] += vtxJ.g()*expr;
				}
			}//if (rij>1e-6)
		}//for j
	} // for r


	for (size_t i = 0; i < I1.size(); ++i)
	{
		I1[i] += selfI1[i];
		I2[i] += selfI2[i];
	}

	tCPUEND = omp_get_wtime();
}//CalcDiffVeloI1I2ToSetOfPointsFromWake(...)


//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void Velocity::CalcDiffVeloI1I2ToSetOfPointsFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	double tCPUSTART, tCPUEND;

	tCPUSTART = omp_get_wtime();

	std::vector<double> selfI1(pointsDb.vtx.size(), 0.0);
	std::vector<Point2D> selfI2(pointsDb.vtx.size(), { 0.0, 0.0 });

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


#pragma omp parallel for default(none) shared(selfI1, selfI2, domainRadius, bnd, pointsDb, std::cout) private(Rij, rij, expr, domRad, diffRadius, left, right, posJx)
	for (int i = 0; i < pointsDb.vtx.size(); ++i)
	{
		const Vortex2D& vtxI = pointsDb.vtx[i];

		domRad = std::max(domainRadius[i], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());

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
						selfI2[i] += (ptG * expr / rij) * Rij;
						selfI1[i] += ptG * expr;
					}
				}//if (rij>1e-6)
			}
		}//for j
	} // for r

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

		std::vector<Point2D> newI2(npt); // I2.size());
		std::vector<double> newI1(npt); // I1.size());
				
		double*& dev_ptr_i1 = pointsDb.devI1Ptr;
		double*& dev_ptr_i2 = pointsDb.devI2Ptr;
		double minRad = W.getPassport().wakeDiscretizationProperties.getMinEpsAst();

		if ((nvt > 0) && (npt > 0))
		{
			//СЕТКА
			if (!useMesh)
				cuCalculateDiffVeloWake(npt, dev_ptr_pt, nvt, dev_ptr_vt, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr, minRad);
			else
				cuCalculateDiffVeloWakeMesh(npt, dev_ptr_pt, nvt, dev_ptr_vt, W.getWake().devMeshPtr, W.getPassport().wakeDiscretizationProperties.epscol, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);

			W.getCuda().CopyMemFromDev<double, 2>(npt, dev_ptr_i2, (double*)newI2.data(), 10);
			W.getCuda().CopyMemFromDev<double, 1>(npt, dev_ptr_i1, newI1.data(), 11);

			if (&pointsDb == &W.getWake())
			{
				for (size_t q = 0; q < I2.size(); ++q)
				{
					I1[q] += newI1[q];
					I2[q] += newI2[q];
				}
			}

			if ((W.getNumberOfBoundary() > 0) && (&pointsDb == &W.getBoundary(0).virtualWake))
			{
				size_t curGlobPnl = 0;

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
			double*& dev_ptr_freeVortexSheetLin = bou.afl.devFreeVortexSheetLinPtr;

			for (size_t q = 1; q < W.getNumberOfBoundary(); ++q)
				npnl += W.getBoundary(q).afl.getNumberOfPanels();

			std::vector<Point2D> newI2(npt);
			std::vector<double> newI1(npt);

			double*& dev_ptr_i1 = pointsDb.devI1Ptr;
			double*& dev_ptr_i2 = pointsDb.devI2Ptr;
			double minRad = W.getPassport().wakeDiscretizationProperties.getMinEpsAst();

			if ((npnl > 0) && (npt > 0))
			{
				/// \todo Реализовать версию с сеткой
				//if (!useMesh)
				cuCalculateDiffVeloWakeFromPanels(npt, dev_ptr_pt, npnl, dev_ptr_r, dev_ptr_freeVortexSheet, dev_ptr_freeVortexSheetLin, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr, minRad);
				//else
				//	cuCalculateDiffVeloWakeMesh(npt, dev_ptr_pt, nvt, dev_ptr_vt, W.getWake().devMeshPtr, W.getPassport().wakeDiscretizationProperties.epscol, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);

				W.getCuda().CopyMemFromDev<double, 2>(npt, dev_ptr_i2, (double*)newI2.data(), 12);
				W.getCuda().CopyMemFromDev<double, 1>(npt, dev_ptr_i1, newI1.data(), 13);
								
				if (&pointsDb == &W.getWake())
				{
					for (size_t q = 0; q < I2.size(); ++q)
					{
						I1[q] += newI1[q];
						I2[q] += newI2[q];
					}
				}

				if ((W.getNumberOfBoundary() > 0) && (&pointsDb == &W.getBoundary(0).virtualWake))
				{
					size_t curGlobPnl = 0;

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

#if defined(USE_CUDA)
void Velocity::GPUDiffVeloFAST(const std::vector<double>& domainRadius, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	size_t nvt = W.getWake().vtx.size();

	std::vector<Point2D> newI2(nvt); 
	std::vector<double> newI1(nvt); 

	double*& dev_ptr_i1 = W.getWake().devI1Ptr;
	double*& dev_ptr_i2 = W.getWake().devI2Ptr;

	size_t* const& dev_nVortices = W.getCuda().dev_ptr_nVortices;

	if (nvt > 0)
	{
		W.getNonConstCuda().RefreshWake(4);
		auto& treeWake = *W.getCuda().inflTreeWake;
		treeWake.MemoryAllocate((int)W.getCuda().n_CUDA_wake);
		treeWake.Update((int)W.getWake().vtx.size(), W.getWake().devVtxPtr, W.getPassport().wakeDiscretizationProperties.eps);
		treeWake.Build();
		treeWake.UpwardTraversal(W.getPassport().numericalSchemes.nbodyMultipoleOrder);

		double time = treeWake.I1I2CalculationWrapper(W.getPassport().wakeDiscretizationProperties.getMinEpsAst(), dev_ptr_i1, (Point2D*)dev_ptr_i2, W.getWake().devRadPtr);

		W.getCuda().CopyMemFromDev<double, 2>(nvt, dev_ptr_i2, (double*)newI2.data(), 10);
		W.getCuda().CopyMemFromDev<double, 1>(nvt, dev_ptr_i1, newI1.data(), 11);

		for (size_t q = 0; q < I2.size(); ++q)
		{
			I1[q] += newI1[q];
			I2[q] += newI2[q];
		}
	}
}
#endif

// Очистка старых массивов под хранение скоростей, выделение новой памяти и обнуление
void Velocity::ResizeAndZero()
{
	W.getTimers().start("ConvVel");	
	wakeVortexesParams.convVelo.clear();
	wakeVortexesParams.convVelo.resize(W.getWake().vtx.size(), { 0.0, 0.0 });
	W.getTimers().stop("ConvVel");

	W.getTimers().start("DiffVel");
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
	W.getTimers().stop("DiffVel");

	//Создаем массивы под виртуальные вихри
	W.getTimers().start("ConvVel");
	virtualVortexesParams.clear();
	virtualVortexesParams.resize(W.getNumberOfBoundary());
	W.getTimers().stop("ConvVel");

	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		W.getTimers().start("ConvVel");
		virtualVortexesParams[bou].convVelo.clear();
		virtualVortexesParams[bou].convVelo.resize(W.getBoundary(bou).virtualWake.vtx.size(), { 0.0, 0.0 });
		W.getTimers().stop("ConvVel");

		W.getTimers().start("DiffVel");
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

		W.getTimers().stop("DiffVel");
	}
}//ResizeAndZero()


void Velocity::SaveVisStress()
{
	if ((W.getPassport().timeDiscretizationProperties.saveVtxStep > 0) && (W.ifDivisible(W.getPassport().timeDiscretizationProperties.saveVisStress)))
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



//Генерация вектора влияния вихревого следа на профиль
void Velocity::GetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeRhs) const
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
void Velocity::GPUGetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeVelo) const
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
		wakeVelo.resize(afl.getNumberOfPanels() * (W.getPassport().numericalSchemes.boundaryCondition.second), 0.0);

}//GPUGetWakeInfluenceToRhs(...)
#endif

#if defined(USE_CUDA)
//Генерация вектора влияния вихревого следа на профиль
void Velocity::GPUFASTGetWakeInfluenceToRhs(const Airfoil & afl, std::vector<double>&wakeVelo) const
{
	size_t shDim = W.getBoundary(afl.numberInPassport).sheetDim;

	const size_t& nvt = W.getWake().vtx.size();
	const size_t& nsr = W.getSource().vtx.size();

	auto& inflTree = *W.getCuda().inflTreeWake;
	const int& order = W.getPassport().numericalSchemes.nbodyMultipoleOrder;

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
			//double timingsToRHS[7];

			double* linPtr = (double*)((shDim == 1) ? nullptr : dev_ptr_rhsLin);

			auto& inflTree = *W.getCuda().inflTreeWake;
			auto& cntrTree = *W.getCuda().cntrTreePnl;

			inflTree.DownwardTraversalVorticesToPanels(cntrTree, (double*)dev_ptr_rhs, 
				linPtr, multipoleTheta, multipoleOrder);


			std::vector<double> newRhs(nTotPan), newRhsLin(nTotPan);

			W.getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_rhs, newRhs.data(), 22);
			if (shDim != 1)
				W.getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_rhsLin, newRhsLin.data(), 22);

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
		wakeVelo.resize(afl.getNumberOfPanels() * (W.getPassport().numericalSchemes.boundaryCondition.second), 0.0);

}//GPUGetWakeInfluenceToRhsFAST(...)
#endif



void Velocity::FillRhs(Eigen::VectorXd& rhsReord) const
{
	Eigen::VectorXd locRhs;
	std::vector<double> lastRhs(W.getNumberOfBoundary());

	size_t currentRow = 0;
	size_t currentSkosRow = 0;

	size_t nAllVars = 0;
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
		nAllVars += W.getBoundary(bou).GetUnknownsSize();

	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
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
#else
		GetWakeInfluenceToRhs(afl, wakeRhs);
#endif	
		std::vector<double> vInfRhs;
		afl.GetInfluenceFromVInfToPanel(vInfRhs);

#pragma omp parallel for \
	default(none) \
	shared(locRhs, afl, bou, wakeRhs, vInfRhs, np) \
	schedule(dynamic, DYN_SCHEDULE)
		for (int i = 0; i < (int)afl.getNumberOfPanels(); ++i)
		{
			locRhs(i) = -vInfRhs[i] - wakeRhs[i] + 0.25 * ((afl.getV(i) + afl.getV(i + 1)) & afl.tau[i]); //0.25 * (afl.getV(i) + afl.getV(i + 1))*afl.tau[i] - прямолинейные
			if (W.getBoundary(bou).sheetDim > 1)
				locRhs(np + i) = -vInfRhs[np + i] - wakeRhs[np + i];

			//!!!!!!!!!!!влияние присоединенных слоев от самого себя и от других профилей				
			//Если включен быстрый метод и есть подвижные тела, то пока что все равно вычисляем IQ
			if (W.isAnyMovableOrDeformable())
			{
				for (size_t q = 0; q < W.getNumberOfBoundary(); ++q)
				{
					const auto& sht = W.getBoundary(q).sheets;
					const auto& iq = W.getIQ(bou, q);

					const Airfoil& aflOther = W.getAirfoil(q);
					if (W.getMechanics(q).isMoves)
					{
						for (size_t j = 0; j < aflOther.getNumberOfPanels(); ++j)
						{
							if ((i != j) || (bou != q))
							{
								locRhs(i) += -iq.first(i, j) * sht.attachedVortexSheet(j, 0);
								locRhs(i) += -iq.second(i, j) * sht.attachedSourceSheet(j, 0);
							}//if (i != j)

						}//for j
					}
				}//for q
			}
		}//for i

		lastRhs[bou] = 0.0;

#pragma omp for
		for (int q = 0; q < (int)afl.gammaThrough.size(); ++q)
			lastRhs[bou] += afl.gammaThrough[q];

		////////////////////////////////////////////////////////////////

		const double currT = W.getCurrentTime();
		const double dt = W.getPassport().timeDiscretizationProperties.dt;

		lastRhs[bou] += (W.getMechanics(bou).circulation - W.getMechanics(bou).circulationOld) * (W.getAirfoil(bou).inverse ? 1.0 : -1.0);

		////////////////////////////////////////////////////////////////

		//размазываем правую часть		
		for (size_t i = 0; i < nVars; ++i)
			rhsReord(i + currentSkosRow) = locRhs(i);

		rhsReord(nAllVars + bou) = lastRhs[bou];

		currentRow += nVars + 1;
		currentSkosRow += nVars;
	}// for bou
}




//Вычисление конвективных скоростей вихрей и виртуальных вихрей в вихревом следе, а также в точках wakeVP
void Velocity::CalcConvVelo()
{
	W.getTimers().start("ConvVel");

	//Влияние следа на след
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONV_TOWAKE))		
	//Вызов виртуальной функции
	std::unique_ptr<BHcu::CudaTreeInfo>& cntrTree = W.getNonConstCuda().cntrTreeWake;
	if (W.getPassport().numericalSchemes.velocityComputation.second == 1)
	{		
		cntrTree->MemoryAllocate((int)W.getCuda().n_CUDA_wake);
		cntrTree->Update((int)W.getWake().vtx.size(), W.getWake().devVtxPtr, W.getPassport().wakeDiscretizationProperties.eps);
		cntrTree->Build();
	}
	GPUCalcConvVeloToSetOfPointsFromWake(W.getNonConstCuda().cntrTreeWake, W.getWake(), wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake, true, true);
	
#else
	//Вызов виртуальной функции
	CalcConvVeloToSetOfPointsFromWake(W.getWake(), wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake, true, true);
#endif

	std::vector<Point2D> nullVector(0);
	//вычисление конвективных скоростей по закону Био-Савара от виртуальных вихрей (точнее, от слоев)
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		//Влияние на след от панелей
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONVVIRT))	

		switch (W.getPassport().numericalSchemes.velocityComputation.second)
		{
		case 0:
			W.getBoundary(bou).GPUCalcConvVelocityToSetOfPointsFromSheets(W.getWake(), wakeVortexesParams.convVelo);
			break;
		case 1:
			if (bou == 0)
				GPUCalcConvVelocityToSetOfPointsFromSheets(W.getNonConstCuda().cntrTreeWake, W.getWake(), wakeVortexesParams.convVelo);
			break;
		}
		//		

#else
		W.getBoundary(bou).CalcConvVelocityToSetOfPointsFromSheets(W.getWake(), wakeVortexesParams.convVelo);
#endif
	}

	//Скорости только что рожденных вихрей
	{
		size_t nVirtVortices = 0;
		for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
			nVirtVortices += W.getBoundary(bou).virtualWake.vtx.size();

		size_t counter = wakeVortexesParams.convVelo.size() - nVirtVortices;

		for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
			for (size_t v = 0; v < W.getBoundary(bou).virtualWake.vtx.size(); ++v)
			{
				wakeVortexesParams.convVelo[counter] = W.getBoundary(bou).afl.getV(W.getBoundary(bou).virtualWake.aflPan[v].second) \
					+ W.getBoundary(bou).virtualWake.vecHalfGamma[v] \
					- W.getV0();
				++counter;
			}
	}

	//Заполнение структур данных виртуальных вихрей
	{
		size_t nVirtVortices = 0;
		for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
			nVirtVortices += W.getBoundary(bou).virtualWake.vtx.size();

		size_t counter = wakeVortexesParams.convVelo.size() - nVirtVortices;

		for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
			for (size_t v = 0; v < W.getBoundary(bou).virtualWake.vtx.size(); ++v)
			{
				virtualVortexesParams[bou].convVelo[v] = wakeVortexesParams.convVelo[counter];
				virtualVortexesParams[bou].epsastWake[v] = wakeVortexesParams.epsastWake[counter];
				++counter;
			}
	}


	W.getTimers().stop("ConvVel");
	W.getTimers().start("VelPres");

	//вычисление скоростей в заданных точках только на соответствующем шаге по времени
	if ((W.getPassport().timeDiscretizationProperties.saveVPstep > 0) && (!(W.getCurrentStep() % W.getPassport().timeDiscretizationProperties.saveVPstep)))
		CalcVeloToWakeVP();

	W.getTimers().stop("VelPres");
}//CalcConvVelo()



//Вычисление конвективных скоростей вихрей в точках wakeVP
void Velocity::CalcVeloToWakeVP()
{
	std::vector<Point2D> velConvWake;
	std::vector<std::vector<Point2D>> velConvBou;

	int addWSize = (int)W.getMeasureVP().getWakeVP().vtx.size();

	velConvWake.resize(addWSize, { 0.0, 0.0 });

	velConvBou.resize(W.getNumberOfBoundary());
	for (size_t i = 0; i < W.getNumberOfBoundary(); ++i)
		velConvBou[i].resize(addWSize, { 0.0, 0.0 });

	std::vector<Point2D>& velocityRef = W.getNonConstMeasureVP().getNonConstVelocity();
	velocityRef.assign(addWSize, W.getV0());


#if (defined(USE_CUDA))
	W.getNonConstCuda().RefreshVP();
#endif

#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_VP))
	std::unique_ptr<BHcu::CudaTreeInfo>& cntrTree = W.getNonConstCuda().cntrTreeVP;
	if (W.getPassport().numericalSchemes.velocityComputation.second == 1)
	{
		if (W.getCurrentStep() == 0 || W.isAnyMovableOrDeformable())
		{
			cntrTree->MemoryAllocate((int)W.getCuda().n_CUDA_velVP);
			cntrTree->Update((int)W.getMeasureVP().getWakeVP().vtx.size(), W.getMeasureVP().getWakeVP().devVtxPtr, 0.0);
			cntrTree->Build();
		}
	}
	GPUCalcConvVeloToSetOfPointsFromWake(cntrTree, W.getMeasureVP().getWakeVP(), velocityRef, W.getNonConstMeasureVP().getNonConstDomainRadius(), true, false);
#else
	CalcConvVeloToSetOfPointsFromWake(W.getMeasureVP().getWakeVP(), velConvWake, W.getNonConstMeasureVP().getNonConstDomainRadius(), true, false);
#endif			   

	for (size_t i = 0; i < W.getNumberOfBoundary(); ++i)
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_VP))			
		switch (W.getPassport().numericalSchemes.velocityComputation.second)
		{
		case 0:
			W.getNonConstBoundary(i).GPUCalcConvVelocityToSetOfPointsFromSheets(W.getMeasureVP().getWakeVP(), velConvBou[i]);

			for (int i = 0; i < addWSize; ++i)
				velocityRef[i] += velConvWake[i];

			for (size_t bou = 0; bou < velConvBou.size(); ++bou)
				for (int j = 0; j < addWSize; ++j)
					velocityRef[j] += velConvBou[bou][j];
			break;
		case 1:
			if (i == 0)
				GPUCalcConvVelocityToSetOfPointsFromSheets(cntrTree, W.getMeasureVP().getWakeVP(), velocityRef);
			break;
		}
#else
		W.getNonConstBoundary(i).CalcConvVelocityToSetOfPointsFromSheets(W.getMeasureVP().getWakeVP(), velConvBou[i]);
#endif
}// CalcVeloToWakeVP()


void Velocity::CalcDiffVeloI1I2ToWakeFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	CalcDiffVeloI1I2ToSetOfPointsFromWake(pointsDb, domainRadius, vorticesDb, I1, I2);
}//CalcDiffVeloI1I2ToWakeFromWake(...)

void Velocity::CalcDiffVeloI1I2ToWakeFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	CalcDiffVeloI1I2ToSetOfPointsFromSheets(pointsDb, domainRadius, bnd, I1, I2);
}//CalcDiffVeloI1I2ToWakeFromSheets(...)