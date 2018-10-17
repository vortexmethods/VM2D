/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MeasureVP.cpp                                                    |
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
\brief Файл кода с описанием класса MeasureVP
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#if defined(_WIN32)
#include <direct.h>
#endif

#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "MeasureVP.h"
#include "World2D.h"
#include "Airfoil.h"
#include "Velocity.h"
#include "Boundary.h"
#include "Preprocessor.h"
#include "WakeDataBase.h"

#include <string.h>



//Чтение точек, в которых нужно посчитать давление и скорость
void MeasureVP::ReadPointsFromFile(const std::string& dir)
{
	std::string filename = dir + "pointsVP";
	std::ifstream VPFile;

	if (fileExistTest(filename, W.getInfo()))
	{
		std::stringstream VPFile(Preprocessor(filename).resultString);

		StreamParser VPParser(W.getInfo(), "velocity & pressure parser", VPFile);

		VPParser.get("n", nPts);
		VPParser.get("points", initialPoints);
	}
}//ReadPointsFromFile(...)

//Расчет и сохранение в файл поля скоростей и давления
void MeasureVP::CalcSaveVP(const std::string& dir, size_t step, timePeriod& time)
{
	time.first = omp_get_wtime();
	std::vector<Point2D> velConvWake;
	std::vector<std::vector<Point2D>> velConvBou;

	if (W.getParallel().myidWork == 0)
	{
		additionalWake.vtx.clear();
		additionalWake.vtx.resize(0);

		Vortex2D addvtx;

		velocity.clear();
		velocity.resize(0);

		pressure.clear();
		pressure.resize(0);

		domainRadius.clear();
		domainRadius.resize(0);

		//	std::vector<WakeDataBase> additionalWakeSource;	//дополнительная база данных вихрей на местах источников для расчета скоростей
		//	additionalWakeSource.resize(W.getNumberOfBoundary());

		//	std::vector<std::vector<Point2D>> additionalWakeSourceConvVelo;
		//	std::vector<std::vector<Point2D>> additionalWakeSourceDiffVelo;

		//определяем точки, которые не находятся внутри профилей
		bool inside;
		for (size_t i = 0; i < initialPoints.size(); i++)
		{
			inside = false;
			for (size_t j = 0; j < W.getNumberOfAirfoil(); j++)
			{
				if (W.getAirfoil(j).isInsideGabarits(initialPoints[i]) && W.getAirfoil(j).IsPointInAirfoil(initialPoints[i]))
				{
					inside = true;
					break;
				}
			}
			if (!inside)
			{
				addvtx.r() = initialPoints[i];
				additionalWake.vtx.push_back(addvtx);
			}
		}

		velConvWake.resize(additionalWake.vtx.size(), { 0.0, 0.0 });
	}

	int addWSize = (int)additionalWake.vtx.size();
	MPI_Bcast(&addWSize, 1, MPI_INT, 0, W.getParallel().commWork);

	velConvBou.resize(W.getNumberOfBoundary());
	for (size_t i = 0; i < W.getNumberOfBoundary(); i++)
		velConvBou[i].resize(additionalWake.vtx.size(), { 0.0, 0.0 });

	/*#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONV_TOWAKE))
		{
		VelocityClass.GPUCalcConvVeloToSetOfPoints(additionalWake, velConvWake, domainRadius);
		for (size_t i = 0; i < W.getNumberOfBoundary(); i++)
		{
		Boundary& BoundaryClass = W.getNonConstBoundary(i);
		BoundaryClass.GPUGetConvVelocityToSetOfPointsFromVirtualVortexes(additionalWake, velConvBou[i]);
		}
		}
		#else
		*/
	{		
	
		W.getNonConstVelocity().CalcConvVeloToSetOfPoints(additionalWake, velConvWake, domainRadius);
			
		for (size_t i = 0; i < W.getNumberOfBoundary(); i++)
		{			
			W.getNonConstBoundary(i).GetConvVelocityToSetOfPointsFromVirtualVortexes(additionalWake, velConvBou[i]);
		}
	}
//#endif


	if (W.getParallel().myidWork == 0)
	{
		velocity.assign(addWSize, W.getPassport().physicalProperties.V0());
		pressure.resize(addWSize);

		for (size_t i = 0; i < additionalWake.vtx.size(); i++)
			velocity[i] += velConvWake[i];

		for (size_t bou = 0; bou < velConvBou.size(); bou++)
		for (size_t j = 0; j < additionalWake.vtx.size(); j++)
			velocity[j] += velConvBou[bou][j];
	}

	int id = W.getParallel().myidWork;
	parProp par = W.getParallel().SplitMPI(addWSize);

	/// \todo Для успокоения добавить WakeSyncronize(), но после обучения его проверке hash

	std::vector<Point2D> locVelocity;
	locVelocity.resize(par.myLen);
	MPI_Scatterv(velocity.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, locVelocity.data(), par.myLen, Point2D::mpiPoint2D, 0, W.getParallel().commWork);

	std::vector<double> locPressure;
	locPressure.resize(par.myLen, 0.0);

	WakeDataBase locAdditionalWake;
	locAdditionalWake.vtx.resize(par.myLen);
	MPI_Scatterv(additionalWake.vtx.data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, locAdditionalWake.vtx.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	
	W.getNonConstVelocity().wakeVortexesParams.convVelo.resize(W.getWake().vtx.size());
	MPI_Bcast(W.getNonConstVelocity().wakeVortexesParams.convVelo.data(), (int)(W.getVelocity().wakeVortexesParams.convVelo.size()), Point2D::mpiPoint2D, 0, W.getParallel().commWork);

	
	/// \todo Пока передаем только средние значения свободного слоя на панелях
	std::vector<std::vector<double>> gamAverPan;
	gamAverPan.resize(W.getNumberOfBoundary());

	if (W.getParallel().myidWork == 0)
	for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
	{
		gamAverPan[bou].resize(W.getBoundary(bou).sheets.freeVortexSheet.size());
		for (int j = 0; j < W.getBoundary(bou).sheets.freeVortexSheet.size(); ++j)
			gamAverPan[bou][j] = W.getBoundary(bou).sheets.freeVortexSheet[j][0];
	}
	
	for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
	{
		int sz = (int)(W.getBoundary(bou).sheets.freeVortexSheet.size());
		MPI_Bcast(&sz, 1, MPI_INT, 0, W.getParallel().commWork);
		if (id > 0)
			gamAverPan[bou].resize(sz);
		MPI_Bcast(gamAverPan[bou].data(), sz, MPI_DOUBLE, 0, W.getParallel().commWork);
	}

	for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
	{
		int sz = (int)(W.getAirfoil(bou).gammaThrough.size());
		MPI_Bcast(&sz, 1, MPI_INT, 0, W.getParallel().commWork);
		MPI_Bcast(W.getNonConstAirfoil(bou).gammaThrough.data(), sz, MPI_DOUBLE, 0, W.getParallel().commWork);
	}

	for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
	{
		int sz = (int)(W.getBoundary(bou).virtualWake.vecHalfGamma.size());
		MPI_Bcast(&sz, 1, MPI_INT, 0, W.getParallel().commWork);
		if (id > 0)
			W.getNonConstBoundary(bou).virtualWake.vecHalfGamma.resize(sz);
		MPI_Bcast(W.getNonConstBoundary(bou).virtualWake.vecHalfGamma.data(), sz, Point2D::mpiPoint2D, 0, W.getParallel().commWork);
	}
	
	Point2D dri;
	Point2D Vi;
	Point2D vi;


	const Point2D& V0 = W.getPassport().physicalProperties.V0();
	const double& eps2 = W.getPassport().wakeDiscretizationProperties.eps2;
	const double& dt = W.getPassport().timeDiscretizationProperties.dt;


	double P0 = /*1*/ 0.5 * V0 * V0;
	

	Point2D cPan;

#pragma warning (push)
#pragma warning (disable: 4101)
	double alpha;
	double dst2eps;
#pragma warning (pop)

#pragma omp parallel for default(none) private(alpha, dri, Vi, vi, dst2eps, cPan) shared(P0, dt, eps2, V0, gamAverPan, id, locAdditionalWake, locPressure, locVelocity, par, std::cout) 
	for (int locI = 0; locI < par.myLen; ++locI)
	{
		//int i = locI + par.myDisp;
		locPressure[locI] = P0;
		locPressure[locI] -=  0.5 * (locVelocity[locI] * locVelocity[locI]); //2

		const Point2D& pt = locAdditionalWake.vtx[locI].r();

		for (int j = 0; j < W.getWake().vtx.size(); j++)
		{
			dri = pt - W.getWake().vtx[j].r();
			Vi = W.getVelocity().wakeVortexesParams.convVelo[j] + V0;

			dst2eps = std::max(dri.length2(), eps2);

			vi = IDPI * W.getWake().vtx[j].g() / dst2eps * dri.kcross();
			locPressure[locI] += vi * Vi; //3
		}

		for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
		for (int j = 0; j < W.getBoundary(bou).virtualWake.vtx.size(); j++)
		{
			dri = pt - W.getBoundary(bou).virtualWake.vtx[j].r();
			dst2eps = std::max(dri.length2(), eps2);

			Vi = W.getBoundary(bou).virtualWake.vecHalfGamma[j];
			vi = IDPI * W.getBoundary(bou).virtualWake.vtx[j].g() / dst2eps * dri.kcross();

			locPressure[locI] +=  vi * Vi; //4
		}

		//for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
		//for (int j = 0; j < W.getAirfoil(bou).np; j++)
		//{
		//	x = (W.getAirfoil(bou).rcm - additionalWake.vtx[i].r()) * (0.5 * (W.getAirfoil(bou).r[j] + W.getAirfoil(bou).r[j + 1]) - additionalWake.vtx[i].r());
		//	y = (W.getAirfoil(bou).rcm - additionalWake.vtx[i].r()) ^ (0.5 * (W.getAirfoil(bou).r[j] + W.getAirfoil(bou).r[j + 1]) - additionalWake.vtx[i].r());
		//	alpha = atan2(y, x);
		//	pressure[i] += IDPI * alpha *\
				//		( (W.getBoundary(bou).sheets.attachedVortexSheet[j][0] - W.getBoundary(bou).oldSheets.attachedVortexSheet[j][0]\
				//		+ W.getBoundary(bou).sheets.freeVortexSheet[j][0]) * W.getAirfoil(bou).len[j] 
				//		 - W.getAirfoil(bou).gammaThrough[j] //5 + 6 + 7
				//		)
		//		/ W.getPassport().timeDiscretizationProperties.dt;
		//}
			
		for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
		{
			const Point2D& rcm = W.getAirfoil(bou).rcm;

			for (int j = 0; j < W.getAirfoil(bou).np; j++)
			{
				cPan = 0.5 * (W.getAirfoil(bou).r[j] + W.getAirfoil(bou).r[j + 1]);
				alpha = atan2((cPan - pt) ^ (rcm - pt),   (cPan - pt)*(rcm - pt));
					
				locPressure[locI] += IDPI * alpha *	gamAverPan[bou][j] * W.getAirfoil(bou).len[j] / dt;

				locPressure[locI] -= IDPI * alpha * W.getAirfoil(bou).gammaThrough[j] / dt;
			}
		}//for bou
	}//for i

	MPI_Gatherv(locPressure.data(), par.myLen, MPI_DOUBLE, pressure.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);


	if (W.getParallel().myidWork == 0)
	{
		std::string fname = "velocityPressure";
		if (step < 10) fname += "0";
		if (step < 100) fname += "0";
		if (step < 1000) fname += "0";
		if (step < 10000) fname += "0";

		std::ostringstream ss;
		ss << step;
		fname += ss.str();
		fname += ".vtk";

		std::ofstream outfile;

#if defined(_WIN32)
		_mkdir((dir + "velocityPressure").c_str());
#else
		mkdir((dir + "velocityPressure").c_str(), S_IRWXU | S_IRGRP | S_IROTH);
#endif

		outfile.open(dir + "velocityPressure/" + fname);

		outfile << "# vtk DataFile Version 2.0" << std::endl;
		outfile << (dir + "velocityPressure/" + fname).c_str() << std::endl;
		outfile << "ASCII" << std::endl;
		outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
		outfile << "POINTS " << additionalWake.vtx.size() << " float" << std::endl;

		for (size_t i = 0; i < additionalWake.vtx.size(); i++)
		{
			double xi = (additionalWake.vtx[i].r())[0];
			double yi = (additionalWake.vtx[i].r())[1];
			outfile << xi << " " << yi << " " << "0.0" << std::endl;
		}//for i

		outfile << "CELLS " << additionalWake.vtx.size() << " " << 2 * additionalWake.vtx.size() << std::endl;
		for (size_t i = 0; i < additionalWake.vtx.size(); ++i)
			outfile << "1 " << i << std::endl;

		outfile << "CELL_TYPES " << additionalWake.vtx.size() << std::endl;
		for (size_t i = 0; i < additionalWake.vtx.size(); ++i)
			outfile << "1" << std::endl;

		outfile << std::endl;
		outfile << "POINT_DATA " << additionalWake.vtx.size() << std::endl;

		outfile << "VECTORS V float" << std::endl;
		for (size_t i = 0; i < additionalWake.vtx.size(); i++)
		{
			outfile << velocity[i][0] << " " << velocity[i][1] << " 0.0" << std::endl;
		}//for i

		outfile << std::endl;

		outfile << "SCALARS P float 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;

		for (size_t i = 0; i < additionalWake.vtx.size(); i++)
		{
			outfile << pressure[i] << std::endl;
		}//for i

		outfile.close();
	}


	time.second = omp_get_wtime();

}//CalcSaveVP()


