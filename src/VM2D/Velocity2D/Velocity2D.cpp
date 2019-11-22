/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
| VM is distributed in the hope that it will be useful, but WITHOUT           |
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
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7   
\date 22 ноября 2019 г.
*/ 

#include "Velocity2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "StreamParser.h"
#include "Wake2D.h"
#include "WakeDataBase2D.h"
#include "World2D.h"
#include "Parallel.h"
#include "Passport2D.h"

using namespace VM2D;

//Вычисление диффузионных скоростей вихрей и виртуальных вихрей в вихревом следе
void Velocity::CalcDiffVelo()
{	
	double tTemp[4];

	tTemp[0] = omp_get_wtime();

	//W.getInfo('t') << "CalcDiffVelo: "
	//	<< virtualVortexesParams[0].I2[0] << " " << virtualVortexesParams[0].I1[0] << std::endl;


	//пелена на пелену
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))	
	GPUCalcDiffVeloI1I2ToSetOfPoints(W.getWake(), wakeVortexesParams.epsastWake, W.getWake(), wakeVortexesParams.I1, wakeVortexesParams.I2);
#else
	CalcDiffVeloI1I2ToSetOfPoints(W.getWake(), wakeVortexesParams.epsastWake, W.getWake(), wakeVortexesParams.I1, wakeVortexesParams.I2);
#endif


	tTemp[1] = omp_get_wtime();
		
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		//не нужно, т.к. сделано выше перед началом вычисления скоростей
		//W.getNonConstBoundary(bou).virtualWake.WakeSynchronize();
		
		//виртуальные на границе на след
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))				
		GPUCalcDiffVeloI1I2ToSetOfPointsFromPanels(W.getWake(), wakeVortexesParams.epsastWake, W.getBoundary(bou), wakeVortexesParams.I1, wakeVortexesParams.I2);
#else			
		CalcDiffVeloI1I2ToSetOfPointsFromPanels(W.getWake(), wakeVortexesParams.epsastWake, W.getBoundary(bou), wakeVortexesParams.I1, wakeVortexesParams.I2);
#endif
		


		// след на виртуальные		
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))					
		GPUCalcDiffVeloI1I2ToSetOfPoints(W.getBoundary(bou).virtualWake, virtualVortexesParams[bou].epsastWake, W.getWake(), virtualVortexesParams[bou].I1, virtualVortexesParams[bou].I2);
#else
		CalcDiffVeloI1I2ToSetOfPoints(W.getBoundary(bou).virtualWake, virtualVortexesParams[bou].epsastWake, W.getWake(), virtualVortexesParams[bou].I1, virtualVortexesParams[bou].I2);
#endif

		
		for (size_t targetBou = 0; targetBou < W.getNumberOfBoundary(); ++targetBou)
		{
		// виртуальные на виртуальные
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))				
			GPUCalcDiffVeloI1I2ToSetOfPointsFromPanels(W.getBoundary(targetBou).virtualWake, virtualVortexesParams[targetBou].epsastWake, W.getBoundary(bou), virtualVortexesParams[targetBou].I1, virtualVortexesParams[targetBou].I2);
#else			
			CalcDiffVeloI1I2ToSetOfPointsFromPanels(W.getBoundary(targetBou).virtualWake, virtualVortexesParams[targetBou].epsastWake, W.getBoundary(bou), virtualVortexesParams[targetBou].I1, virtualVortexesParams[targetBou].I2);
#endif	
		}

		//W.getInfo('t') << "CalcDiffVelo: "
		//	<< virtualVortexesParams[bou].I2[0] << " " << virtualVortexesParams[bou].I1[0] << std::endl;

		//exit(-123);
	} //for bou

tTemp[2] = omp_get_wtime();

	Point2D I2;
	double I1;
	
	for (size_t vt = 0; vt < wakeVortexesParams.diffVelo.size(); ++vt)
	{
		I2 = wakeVortexesParams.I2[vt];
		I1 = wakeVortexesParams.I1[vt];
		if (fabs(I1) < 1.e-8)
			wakeVortexesParams.diffVelo[vt] = { 0.0, 0.0 };
		else
			wakeVortexesParams.diffVelo[vt] = I2 * (1.0 / (I1 * wakeVortexesParams.epsastWake[vt]));
	}
	
	
	
	for (size_t targetBou = 0; targetBou < W.getNumberOfBoundary(); ++targetBou)
	for (size_t vt = 0; vt < virtualVortexesParams[targetBou].diffVelo.size(); ++vt)
	{
		I2 = virtualVortexesParams[targetBou].I2[vt];
		I1 = virtualVortexesParams[targetBou].I1[vt];

		if (fabs(I1) < 1.e-8)
			virtualVortexesParams[targetBou].diffVelo[vt] = { 0.0, 0.0 };
		else
			virtualVortexesParams[targetBou].diffVelo[vt] = I2 * (1.0 / (I1 * virtualVortexesParams[targetBou].epsastWake[vt]));
	}

	

	tTemp[3] = omp_get_wtime();

#pragma warning (push)
#pragma warning (disable: 4010)
	//W.getInfo('t') << "Times_diff:" << \
		tTemp[1] - tTemp[0] << " " << \
		tTemp[2] - tTemp[1] << " " << \
		tTemp[3] - tTemp[2] << " " << \
		std::endl;
#pragma warning (pop)
	
}//CalcDiffVelo()


//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void Velocity::CalcDiffVeloI1I2ToSetOfPoints(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	double tCPUSTART, tCPUEND;

	tCPUSTART = omp_get_wtime();

	std::vector<double> selfI1;
	std::vector<Point2D> selfI2;

	const int& id = W.getParallel().myidWork;

	VMlib::parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
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
#pragma warning (pop)

#pragma omp parallel for default(none) shared(locI1, locI2, domainRadius, locPoints, vorticesDb, par) private(Rij, rij, expr, diffRadius, left, right, posJx)
	for (int i = 0; i < par.myLen; ++i)
	{
		const Vortex2D& vtxI = locPoints[i];

		/// \todo Понять природу магической константы 8.0 и синхронизировать с GPU

		diffRadius = 8.0 * domainRadius[i + par.myDisp];

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
					expr = exp(-rij / domainRadius[i + par.myDisp]);
					locI2[i] += (vtxJ.g()* expr / rij) * Rij;
					locI1[i] += vtxJ.g()*expr;
				}
			}//if (rij>1e-6)
		}//for j
	} // for r

	//std::ostringstream sss;
	//sss << "velo_";
	//std::ofstream veloFile(sss.str());
	//for (int i = 0; i < domainRadius.size(); ++i)
	//	veloFile << domainRadius[i] << std::endl;
	//veloFile.close();

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
}//CalcDiffVeloI1I2ToSetOfPoints(...)



//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void Velocity::CalcDiffVeloI1I2ToSetOfPointsFromPanels(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2)
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

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
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
#pragma warning (pop)


#pragma omp parallel for default(none) shared(locI1, locI2, domainRadius, locPoints, bnd, par, std::cout) private(Rij, rij, expr, diffRadius, left, right, posJx)
	for (int i = 0; i < par.myLen; ++i)
	{
		const Vortex2D& vtxI = locPoints[i];

		/// \todo Понять природу магической константы 8.0 и синхронизировать с GPU
		diffRadius = 8.0 * domainRadius[i + par.myDisp];

		left = vtxI.r()[0] - diffRadius;
		right = vtxI.r()[0] + diffRadius;

		for (size_t j = 0; j < bnd.sheets.getSheetSize(); ++j)
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
						expr = exp(-rij / domainRadius[i + par.myDisp]);
						locI2[i] += (ptG * expr / rij) * Rij;
						locI1[i] += ptG * expr;
					}
				}//if (rij>1e-6)
			}
		}//for j
	} // for r

	//std::ostringstream sss;
	//sss << "velo_";
	//std::ofstream veloFile(sss.str());
	//for (int i = 0; i < domainRadius.size(); ++i)
	//	veloFile << domainRadius[i] << std::endl;
	//veloFile.close();

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
}//CalcDiffVeloI1I2ToSetOfPointsFromPanels(...)


#if defined (USE_CUDA)
//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void Velocity::GPUCalcDiffVeloI1I2ToSetOfPoints(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh)
{
	const size_t npt = pointsDb.vtx.size();
	double*& dev_ptr_pt = pointsDb.devVtxPtr;
	double*& dev_ptr_dr = pointsDb.devRadPtr;

	const size_t nvt = vorticesDb.vtx.size();
	double*& dev_ptr_vt = vorticesDb.devVtxPtr;

	std::vector<double>& loci1 = pointsDb.tmpI1;
	std::vector<Point2D>& loci2 = pointsDb.tmpI2;
	double*& dev_ptr_i1 = pointsDb.devI1Ptr;
	double*& dev_ptr_i2 = pointsDb.devI2Ptr;

	const int& id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(npt, true);


	double tCUDASTART = 0.0, tCUDAEND = 0.0, tCUDAENDCALC = 0.0;

	tCUDASTART = omp_get_wtime();

	if ((nvt > 0) && (npt > 0))
	{
		if (!useMesh)
			cuCalculateDiffVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);
		else
			cuCalculateDiffVeloWakeMesh(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, W.getWake().devMeshPtr, W.getPassport().wakeDiscretizationProperties.epscol, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);

		W.getCuda().CopyMemFromDev<double, 2>(par.myLen, dev_ptr_i2, (double*)&loci2[0]);
		W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_i1, &loci1[0]);

		std::vector<Point2D> newI2;
		std::vector<double> newI1;
		if (id == 0)
		{
			newI2.resize(I2.size());
			newI1.resize(I1.size());
		}

		MPI_Gatherv(loci2.data(), par.myLen, Point2D::mpiPoint2D, newI2.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
		MPI_Gatherv(loci1.data(), par.myLen, MPI_DOUBLE, newI1.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

		tCUDAENDCALC = omp_get_wtime();

		if (id == 0)
			for (size_t q = 0; q < I2.size(); ++q)
			{
				I1[q] += newI1[q];
				I2[q] += newI2[q];
			}
	}
	tCUDAEND = omp_get_wtime();
	//W.getInfo('t') << "DIFF_GPU(" << id << ", " << par.myLen << "): " << (tCUDAEND - tCUDASTART) << " " << (tCUDAENDCALC - tCUDASTART) << std::endl;
}//GPUCalcDiffVeloI1I2ToSetOfPoints(...)



//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void Velocity::GPUCalcDiffVeloI1I2ToSetOfPointsFromPanels(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bou, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh)
{
	const size_t npt = pointsDb.vtx.size();
	double*& dev_ptr_pt = pointsDb.devVtxPtr;
	double*& dev_ptr_dr = pointsDb.devRadPtr;

	const size_t npnl = bou.afl.getNumberOfPanels(); //vorticesDb.vtx.size();
	double*& dev_ptr_r = bou.afl.devRPtr;
	double*& dev_ptr_freeVortexSheet = bou.afl.devFreeVortexSheetPtr;
	//double*& dev_ptr_attachedVortexSheet = bou.afl.devAttachedVortexSheetPtr;
	//double*& dev_ptr_attachedSourceSheet = bou.afl.devAttachedSourceSheetPtr;

	std::vector<double>& loci1 = pointsDb.tmpI1;
	std::vector<Point2D>& loci2 = pointsDb.tmpI2;
	double*& dev_ptr_i1 = pointsDb.devI1Ptr;
	double*& dev_ptr_i2 = pointsDb.devI2Ptr;

	const int& id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(npt, true);


	double tCUDASTART = 0.0, tCUDAEND = 0.0, tCUDAENDCALC = 0.0;

	tCUDASTART = omp_get_wtime();

	if ((npnl > 0) && (npt > 0))
	{
		/// \todo Реализовать версию с сеткой
		//if (!useMesh)
		
		
		/*
		{

		//	std::cout << "sizeof(Vortex2D) = " << sizeof(Vortex2D) << std::endl;
		//	std::cout << "sizeof(Point2D) = " << sizeof(Point2D) << std::endl;		
			std::ofstream fo("gpuData.txt");
			fo << "pt: " << par.myDisp << " " << par.myLen << std::endl;
			fo << "--------" << std::endl;

			double* TPT = new double[3 * par.myLen];
			W.getCuda().CopyMemFromDev<double, 3>(par.myLen, dev_ptr_pt, TPT);
			for (int i = 0; i < par.myLen; ++i)
			{
				fo << i << " " << TPT[3 * i] << " " << TPT[3 * i + 1] << std::endl;
			}
			delete[] TPT;


			fo << "pnl: " << npnl << std::endl;
			fo << "--------" << std::endl;

			double* TPNL = new double[2 * npnl];
			W.getCuda().CopyMemFromDev<double, 2>(npnl, dev_ptr_r, TPNL);
			for (int i = 0; i < npnl; ++i)
			{
				fo << i << " " << TPNL[2 * i] << " " << TPNL[2 * i + 1] << std::endl;
			}
			delete[] TPNL;

			fo << "gam: " << npnl << std::endl;
			fo << "--------" << std::endl;

			double* TGM = new double[npnl];
			W.getCuda().CopyMemFromDev<double, 1>(npnl, dev_ptr_freeVortexSheet, TGM);
			for (int i = 0; i < npnl; ++i)
			{
				fo << i << " " << TGM[i] << std::endl;
			}
			delete[] TGM;


			fo << "rd: " << par.myLen << std::endl;
			fo << "--------" << std::endl;

			double* TRD = new double[par.myLen];
			W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_dr, TRD);
			for (int i = 0; i < par.myLen; ++i)
			{
				fo << i << " " << TRD[i] << std::endl;
			}
			delete[] TRD;


			fo.close();

		//	//exit(-122);
		}		
		*/

		cuCalculateDiffVeloWakeFromPanels(par.myDisp, par.myLen, dev_ptr_pt, npnl, dev_ptr_r, dev_ptr_freeVortexSheet, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);
		//
		
		/*
		{

		//	std::cout << "sizeof(Vortex2D) = " << sizeof(Vortex2D) << std::endl;
		//	std::cout << "sizeof(Point2D) = " << sizeof(Point2D) << std::endl;

			std::ofstream fo("gpuDataOUT.txt");

			double* TI2 = new double[2 * par.myLen];
			W.getCuda().CopyMemFromDev<double, 2>(par.myLen, dev_ptr_i2, TI2);
			for (int i = 0; i < par.myLen; ++i)
			{
				fo << i << " " << TI2[2 * i] << " " << TI2[2 * i + 1] << std::endl;
			}
			delete[] TI2;

			fo << "-----------" << std::endl;

			double* TI1 = new double[par.myLen];
			W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_i1, TI1);
			for (int i = 0; i < par.myLen; ++i)
			{
				fo << i << " " << TI1[i] << std::endl;
			}
			delete[] TI1;
			


			fo.close();

			//exit(-122);
		}
		*/
		
		
		
		//else
		//	cuCalculateDiffVeloWakeMesh(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, W.getWake().devMeshPtr, W.getPassport().wakeDiscretizationProperties.epscol, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);

		W.getCuda().CopyMemFromDev<double, 2>(par.myLen, dev_ptr_i2, (double*)&loci2[0]);
		W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_i1, &loci1[0]);

		std::vector<Point2D> newI2;
		std::vector<double> newI1;
		if (id == 0)
		{
			newI2.resize(I2.size());
			newI1.resize(I1.size());
		}

		MPI_Gatherv(loci2.data(), par.myLen, Point2D::mpiPoint2D, newI2.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
		MPI_Gatherv(loci1.data(), par.myLen, MPI_DOUBLE, newI1.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

		tCUDAENDCALC = omp_get_wtime();

		if (id == 0)
			for (size_t q = 0; q < I2.size(); ++q)
			{
				I1[q] += newI1[q];
				I2[q] += newI2[q];
			}
	}
	tCUDAEND = omp_get_wtime();
	//W.getInfo('t') << "DIFF_GPU(" << id << ", " << par.myLen << "): " << (tCUDAEND - tCUDASTART) << " " << (tCUDAENDCALC - tCUDASTART) << std::endl;
}//GPUCalcDiffVeloI1I2ToSetOfPoints(...)
#endif