/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.3    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/09/26     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MeasureVelocityPressure.cpp                                      |
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
\brief Файл кода с описанием класса MeasureVelocityPressure
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/

#if !defined(__linux__)
#include <direct.h>
#endif

#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "MeasureVelocityPressure.h"
#include "World2D.h"
#include "Airfoil.h"
#include "Velocity.h"
#include "Boundary.h"
#include "Preprocessor.h"
#include "WakeDataBase.h"

#include <string.h>



//Чтение точек, в которых нужно посчитать давление и скорость
void MeasureVelocityPressure::ReadPointsFromFile(const std::string& dir)
{
	std::string filename = dir + "pointsVelocityPressure";
	std::ifstream velocityPressureFile;

	std::ostream *Pinfo, *Perr;
	if (W.getParallel().myidWork == 0)
	{
		Pinfo = defaults::defaultPinfo;
		Perr = defaults::defaultPerr;
	}
	else
	{
		Pinfo = nullptr;
		Perr = nullptr;
	}

	if (fileExistTest(filename, Pinfo, Perr, "velocityPressure"))
	{
		std::stringstream velocityPressureFile(Preprocessor(filename).resultString);

		StreamParser velocityPressureParser(velocityPressureFile);

		velocityPressureParser.get("n", nPts);

		velocityPressureParser.get("points", initialPoints);
	}
}//ReadPointsFromFile(...)

//Расчет и сохранение в файл поля скоростей и давления
void MeasureVelocityPressure::CalcSaveVelocityPressure(const std::string& dir, size_t step, timePeriod& time)
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

		if (W.getParallel().myidWork > 0)
			additionalWake.vtx.resize(addWSize);

		MPI_Bcast(additionalWake.vtx.data(), addWSize, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

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

		velocity.assign(additionalWake.vtx.size(), W.getPassport().physicalProperties.V0());
		pressure.resize(additionalWake.vtx.size());

		for (size_t i = 0; i < additionalWake.vtx.size(); i++)
			velocity[i] += velConvWake[i];

		for (size_t bou = 0; bou < velConvBou.size(); bou++)
		for (size_t j = 0; j < additionalWake.vtx.size(); j++)
			velocity[j] += velConvBou[bou][j];

		Point2D dri;
		Point2D Vi;
		Point2D vi;

		double alpha;

		double P0 = /*1*/ 0.5 * W.getPassport().physicalProperties.V0() * W.getPassport().physicalProperties.V0();
		double dst2eps;
				
		for (size_t i = 0; i < additionalWake.vtx.size(); i++)
		{
			pressure[i] = P0;
			pressure[i] -= /*2*/ 0.5 * velocity[i] * velocity[i];

			const Point2D& pt = additionalWake.vtx[i].r();

			for (int j = 0; j < W.getWake().vtx.size(); j++)
			{
				dri = pt - W.getWake().vtx[j].r();
				Vi = W.getVelocity().wakeVortexesParams.convVelo[j] + \
					W.getPassport().physicalProperties.V0();

				dst2eps = std::max(dri.length2(), W.getPassport().wakeDiscretizationProperties.eps2);

				vi = IDPI * W.getWake().vtx[j].g() / dst2eps * dri.kcross();
				pressure[i] += /*3*/ vi * Vi;
			}

			for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
			for (int j = 0; j < W.getBoundary(bou).virtualWake.vtx.size(); j++)
			{
				dri = pt - W.getBoundary(bou).virtualWake.vtx[j].r();

				//			Vi = 2.0 * W.getBoundary(bou).virtualWake.vtx[j].g() / (W.getAirfoil(bou).len[i]);
				dst2eps = std::max(dri.length2(), W.getPassport().wakeDiscretizationProperties.eps2);

				Vi = W.getBoundary(bou).virtualWakeVelocity[j];
				vi = IDPI * W.getBoundary(bou).virtualWake.vtx[j].g() / dst2eps * dri.kcross();

				pressure[i] += /*4*/ vi * Vi;
			}

			//for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
			//for (int j = 0; j < W.getAirfoil(bou).np; j++)
			//{
			//	x = (W.getAirfoil(bou).rcm - additionalWake.vtx[i].r()) * (0.5 * (W.getAirfoil(bou).r[j] + W.getAirfoil(bou).r[j + 1]) - additionalWake.vtx[i].r());
			//	y = (W.getAirfoil(bou).rcm - additionalWake.vtx[i].r()) ^ (0.5 * (W.getAirfoil(bou).r[j] + W.getAirfoil(bou).r[j + 1]) - additionalWake.vtx[i].r());
			//	alpha = atan2(y, x);
			//	pressure[i] += IDPI * alpha *\
					//		(/*5 + 6*/ (W.getBoundary(bou).sheets.attachedVortexSheet[j][0] - W.getBoundary(bou).oldSheets.attachedVortexSheet[j][0]\
					//		+ W.getBoundary(bou).sheets.freeVortexSheet[j][0]) * W.getAirfoil(bou).len[j]
			//		/*7*/ - W.getAirfoil(bou).gammaThrough[j]
			//		)
			//		/ W.getPassport().timeDiscretizationProperties.dt;
			//}
			
			for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
			{
				const Point2D& rcm = W.getAirfoil(bou).rcm;

				for (int j = 0; j < W.getAirfoil(bou).np; j++)
				{
					const Point2D& cPan = 0.5 * (W.getAirfoil(bou).r[j] + W.getAirfoil(bou).r[j + 1]);
					alpha = atan2((cPan - pt) ^ (rcm - pt),   (cPan - pt)*(rcm - pt));
					pressure[i] += IDPI * alpha *\
						(W.getBoundary(bou).sheets.freeVortexSheet[j][0]) * W.getAirfoil(bou).len[j] / W.getPassport().timeDiscretizationProperties.dt;

					pressure[i] += -IDPI * alpha * W.getAirfoil(bou).gammaThrough[j] / W.getPassport().timeDiscretizationProperties.dt;
				}
			}//for bou
		}//for i
				
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

#if !defined(__linux__)
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

		/*
		outfile << "SCALARS Vx float 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;

		for (size_t i = 0; i < additionalWake.vtx.size(); i++)
		{
			outfile << velocity[i][0] << std::endl;
		}//for i

		outfile << std::endl;

		outfile << "SCALARS Vy float 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;

		for (size_t i = 0; i < additionalWake.vtx.size(); i++)
		{
			outfile << velocity[i][1] << std::endl;
		}//for i
		*/
		
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

}//CalcSaveVelocityPressure()


