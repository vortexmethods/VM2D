/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: VelocityBiotSavart.cpp                                           |
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
\version 1.4
\date 16 октября 2018 г.
*/


#include "VelocityBiotSavart.h"
#include "World2D.h"


/// Конструктор
VelocityBiotSavart::VelocityBiotSavart(const World2D& W_) :
	Velocity(W_)
{
};

/// Деструктор
VelocityBiotSavart::~VelocityBiotSavart()
{	
};

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

//Вычисление конвективных скоростей и радиусов вихревых доменов в заданном наборе точек
void VelocityBiotSavart::CalcConvVeloToSetOfPoints(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius)
{
	double tCPUSTART, tCPUEND;
	
	tCPUSTART = omp_get_wtime();
	std::vector<Point2D> selfVelo;
	
	const int& id = W.getParallel().myidWork;

	parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size(), true);

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	//Заполнение "своей" части массива скоростей
	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		         locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	double cft = IDPI;

	std::vector<Point2D> locConvVelo;
	locConvVelo.resize(par.myLen);

	std::vector<double> locDomRadius;
	locDomRadius.resize(par.myLen);

#pragma warning (push)
#pragma warning (disable: 4101)
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2eps, dst2;	
#pragma warning (pop)
	
#pragma omp parallel for default(none) shared(locConvVelo, locDomRadius, locPoints, cft, par) private(tempVel, velI, dst2, dst2eps) schedule(dynamic, DYN_SCHEDULE)
	for (int i = 0; i < par.myLen; ++i)
	{
		double ee2[3] = { 10000.0, 10000.0, 10000.0 };
		
		/*
		if (wake.vtx.size() == 0)
		{
			ee2[0] = 0.0;
			ee2[1] = 0.0;
			ee2[2] = 0.0;
		}
		*/
		
		velI.toZero();
		
		const Point2D& posI = locPoints[i].r();

		for (size_t j = 0; j < W.getWake().vtx.size(); ++j)
		{
			const Point2D& posJ = W.getWake().vtx[j].r();
			const double& gamJ = W.getWake().vtx[j].g();

			tempVel.toZero();

			dst2 = dist2(posI, posJ);

			//Модифицируем массив квадратов расстояний до ближайших вихрей из wake
			ModifyE2(ee2, dst2);

			dst2eps = std::max(dst2, W.getPassport().wakeDiscretizationProperties.eps2);
			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2eps);
			velI += tempVel;
		}
		
		for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
		{
			for (size_t j = 0; j < W.getBoundary(s).virtualWake.vtx.size(); ++j)
			{
				const Point2D& posJ = W.getBoundary(s).virtualWake.vtx[j].r();
				const double& gamJ = W.getBoundary(s).virtualWake.vtx[j].g();

				dst2 = dist2(posI, posJ);

				//Модифицируем массив квадратов расстояний до ближайших вихрей из virtualWake
				ModifyE2(ee2, dst2);
			}
		}

		velI *= cft;
		locConvVelo[i] = velI;

		locDomRadius[i] = std::max(sqrt((ee2[0] + ee2[1] + ee2[2]) / 3.0), 2.0*W.getPassport().wakeDiscretizationProperties.epscol);
	}


	if (id == 0)
		selfVelo.resize(pointsDb.vtx.size());


	//if (id == 0)
	//это нужно всем, т.к. ниже стоит Allgatherv
	domainRadius.resize(par.totalLen);

	MPI_Gatherv(locConvVelo.data(), par.myLen, Point2D::mpiPoint2D, selfVelo.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
	
	//parallel.BCastAllLenDisp();
	MPI_Allgatherv(locDomRadius.data(), par.myLen, MPI_DOUBLE, domainRadius.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, W.getParallel().commWork);

	if (id == 0)
		for (size_t i = 0; i < velo.size(); ++i)
			velo[i] += selfVelo[i];
	
	tCPUEND = omp_get_wtime();
	//W.getInfo('t') << "CONV_CPU: " << tCPUEND - tCPUSTART << std::endl;

}//CalcConvVeloToSetOfPoints(...)

#if defined(USE_CUDA)
void VelocityBiotSavart::GPUCalcConvVeloToSetOfPoints(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius)
{	
	const size_t npt = pointsDb.vtx.size();
	double*& dev_ptr_pt = pointsDb.devVtxPtr;

	const size_t nvt = W.getWake().vtx.size();
	double*& dev_ptr_vt = W.getWake().devVtxPtr;
	const size_t nbou = W.getNumberOfBoundary();
	
	size_t* const& dev_nPanels = W.getCuda().dev_ptr_nPanels;
	size_t* const& dev_nVortices = W.getCuda().dev_ptr_nVortices;


	double** const & dev_ptr_ptr_vtx = W.getCuda().dev_ptr_ptr_vtx;
	
	std::vector<Point2D>& Vel = velo;
	std::vector<double>& Rad = domainRadius;
	std::vector<Point2D>& locvel = pointsDb.tmpVel;
	std::vector<double>& locrad = pointsDb.tmpRad;
	double*& dev_ptr_vel = pointsDb.devVelPtr;
	double*& dev_ptr_rad = pointsDb.devRadPtr;
	double minRad = 2.0*W.getPassport().wakeDiscretizationProperties.epscol;
	const double& eps2 = W.getPassport().wakeDiscretizationProperties.eps2;
	

	const int& id = W.getParallel().myidWork;
	parProp par = W.getParallel().SplitMPI(npt, true);

	double tCUDASTART = 0.0, tCUDAEND = 0.0;

	tCUDASTART = omp_get_wtime();

	if (npt > 0)
	{
		cuCalculateConvVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nvt, dev_ptr_vt, nbou, dev_nVortices, dev_ptr_ptr_vtx, dev_ptr_vel, dev_ptr_rad, minRad, eps2);

		W.getCuda().CopyMemFromDev<double, 2>(par.myLen, dev_ptr_vel, (double*)&locvel[0]);
		W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_rad, &locrad[0]);

		std::vector<Point2D> newV;
		if (id == 0)
			newV.resize(Vel.size());

		MPI_Gatherv(locvel.data(), par.myLen, Point2D::mpiPoint2D, newV.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
		if (id == 0)
			for (size_t q = 0; q < Vel.size(); ++q)
				Vel[q] += newV[q];

		Rad.resize(par.totalLen);
		MPI_Allgatherv(locrad.data(), par.myLen, MPI_DOUBLE, Rad.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, W.getParallel().commWork);

		cuCopyFixedArray(dev_ptr_rad, Rad.data(), sizeof(double) * Rad.size());
	}

	tCUDAEND = omp_get_wtime();

	//W.getInfo('t') << "CONV_GPU(" << parallel.myidWork << "): " << (tCUDAEND - tCUDASTART) << std::endl;
}//GPUCalcConvVeloToSetOfPoints(...)
#endif




//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void VelocityBiotSavart::CalcDiffVeloI1I2ToSetOfPoints(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	double tCPUSTART, tCPUEND;

	tCPUSTART = omp_get_wtime();

	std::vector<double> selfI1;
	std::vector<Point2D> selfI2;

	const int& id = W.getParallel().myidWork;

	parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	std::vector<double> locI1(par.myLen, 0.0);
	std::vector<Point2D> locI2(par.myLen, {0.0, 0.0});

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

	MPI_Gatherv(locI1.data(), par.myLen, MPI_DOUBLE,          selfI1.data(), par.len.data(), par.disp.data(), MPI_DOUBLE,          0, W.getParallel().commWork);
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
void VelocityBiotSavart::CalcDiffVeloI1I2ToSetOfPointsFromPanels(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	double tCPUSTART, tCPUEND;

	tCPUSTART = omp_get_wtime();

	std::vector<double> selfI1;
	std::vector<Point2D> selfI2;

	const int& id = W.getParallel().myidWork;
	parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size());

	
	/// \todo Временно сделана довольно извращенная, но синхронизация свободного вихревого слоя
	std::vector<double> freeVortexSheetGamma;
	int sz;

	if (id == 0)
	{
		sz = (int)(bnd.sheets.freeVortexSheet.size());
		freeVortexSheetGamma.resize(sz, 0.0);
		for (size_t j = 0; j < bnd.sheets.freeVortexSheet.size(); ++j)
			freeVortexSheetGamma[j] = bnd.sheets.freeVortexSheet[j][0];
	}
	MPI_Bcast(&sz, 1, MPI_INT, 0, W.getParallel().commWork);
	if (id != 0)
		freeVortexSheetGamma.resize(sz, 0.0);
	//std::cout << "id = " << id << ", size = " << freeVortexSheetGamma.size() << std::endl;
	MPI_Bcast(freeVortexSheetGamma.data(), (int)sz, MPI_DOUBLE, 0, W.getParallel().commWork);


	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	std::vector<double> locI1(par.myLen, 0.0);
	std::vector<Point2D> locI2(par.myLen, {0.0, 0.0});

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


#pragma omp parallel for default(none) shared(locI1, locI2, domainRadius, locPoints, bnd, par, std::cout, freeVortexSheetGamma) private(Rij, rij, expr, diffRadius, left, right, posJx)
	for (int i = 0; i < par.myLen; ++i)
	{
		const Vortex2D& vtxI = locPoints[i];

		/// \todo Понять природу магической константы 8.0 и синхронизировать с GPU
		diffRadius = 8.0 * domainRadius[i + par.myDisp];

		left = vtxI.r()[0] - diffRadius;
		right = vtxI.r()[0] + diffRadius;

		//std::cout << bnd.sheets.freeVortexSheet[4][0] << std::endl;

		for (size_t j = 0; j < bnd.sheets.freeVortexSheet.size(); ++j)
		{			
			/// \todo Сделать переменной и синхронизировать с GPU
			const int nQuadPt = 1;
			
			const double ptG = freeVortexSheetGamma[j] * bnd.afl.len[j] / nQuadPt;

			for (int q = 0; q < nQuadPt; ++q)
			{
				const Point2D& ptJ = bnd.CC[j] + bnd.afl.tau[j] * (q + 0.5) * bnd.afl.len[j] * (1.0 / nQuadPt);  // vorticesDb.vtx[j];
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
void VelocityBiotSavart::GPUCalcDiffVeloI1I2ToSetOfPoints(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh)
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
	parProp par = W.getParallel().SplitMPI(npt, true);


	double tCUDASTART = 0.0, tCUDAEND = 0.0, tCUDAENDCALC = 0.0;

	tCUDASTART = omp_get_wtime();

	if ( (nvt > 0) && (npt > 0) )
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
//vo VelocityBiotSavart::GPUCalcDiffVeloI1I2ToSetOfPoints          (const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh)
void VelocityBiotSavart::GPUCalcDiffVeloI1I2ToSetOfPointsFromPanels(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bou,            std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh)
{
	const size_t npt = pointsDb.vtx.size();
	double*& dev_ptr_pt = pointsDb.devVtxPtr;
	double*& dev_ptr_dr = pointsDb.devRadPtr;

	const size_t npnl = bou.afl.np; //vorticesDb.vtx.size();
	double*& dev_ptr_r = bou.afl.devRPtr;
	double*& dev_ptr_freeVortexSheet = bou.afl.devFreeVortexSheetPtr;
	//double*& dev_ptr_attachedVortexSheet = bou.afl.devAttachedVortexSheetPtr;
	//double*& dev_ptr_attachedSourceSheet = bou.afl.devAttachedSourceSheetPtr;

	std::vector<double>& loci1 = pointsDb.tmpI1;
	std::vector<Point2D>& loci2 = pointsDb.tmpI2;
	double*& dev_ptr_i1 = pointsDb.devI1Ptr;
	double*& dev_ptr_i2 = pointsDb.devI2Ptr;

	const int& id = W.getParallel().myidWork;
	parProp par = W.getParallel().SplitMPI(npt, true);


	double tCUDASTART = 0.0, tCUDAEND = 0.0, tCUDAENDCALC = 0.0;

	tCUDASTART = omp_get_wtime();

	if ( (npnl > 0) && (npt > 0) )
	{
		/// \todo Реализовать версию с сеткой
		//if (!useMesh)
			cuCalculateDiffVeloWakeFromPanels(par.myDisp, par.myLen, dev_ptr_pt, npnl, dev_ptr_r, dev_ptr_freeVortexSheet, dev_ptr_i1, dev_ptr_i2, dev_ptr_dr);
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