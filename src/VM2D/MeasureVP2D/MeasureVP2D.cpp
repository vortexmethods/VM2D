/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.6    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/10/28     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MeasureVP2D.cpp                                                  |
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
\brief Файл кода с описанием класса MeasureVP
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.6   
\date 28 октября 2019 г.
*/

#if defined(_WIN32)
#include <direct.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <string.h>

#include "MeasureVP2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "Boundary2DConstLayerAver.h"

#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

//Конструктор
MeasureVP::MeasureVP(const World2D& W_)
	: W(W_)
{
	additionalWake.reset(new WakeDataBase(W_));
};//MeasureVP(...)

//Чтение точек, в которых нужно посчитать давление и скорость
void MeasureVP::ReadPointsFromFile(const std::string& dir)
{
	std::string filename = dir + "pointsVP";
	std::ifstream VPFile;

	if (fileExistTest(filename, W.getInfo()))
	{
		std::stringstream VPFile(VMlib::Preprocessor(filename).resultString);

		VMlib::StreamParser VPParser(W.getInfo(), "velocity & pressure parser", VPFile);
		
		VPParser.get("points", initialPoints);
		VPParser.get("history", historyPoints);
				
		if (historyPoints.size() > 0)
		{
			VMlib::CreateDirectory(dir, "velPres");
			std::string VPFileNameList = W.getPassport().dir + "velPres/listPoints";
			std::ofstream VPFileList(VPFileNameList.c_str());

			VMlib::PrintLogoToTextFile(VPFileList, VPFileNameList.c_str(), "List of points, where velocity and pressure are measured in csv-files");

			VMlib::PrintHeaderToTextFile(VPFileList, "No.     FileName      X     Y");

			for (int q = 0; q < historyPoints.size(); ++q)
			{
				VPFileList << '\n' << q << '\t' << VMlib::fileNameStep("VP-atPoint-", 2, q, "csv") << '\t' << historyPoints[q][0] << '\t' << historyPoints[q][1];

				std::string VPFileNameCsv;
				VPFileNameCsv = W.getPassport().dir + "velPres/" + VMlib::fileNameStep("VP-atPoint-", 2, q, "csv");

				std::ofstream VPFileCsv(VPFileNameCsv.c_str());

				VPFileCsv << "t,Vx,Vy,p" << std::endl;
				VPFileCsv.close();
				VPFileCsv.clear();
			}

			VPFileList.close();
			VPFileList.clear();
		}
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
		additionalWake->vtx.clear();
		additionalWake->vtx.resize(0);

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


		//отключаем проверку того, что точки расположены вне профиля
		/*
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
		*/

		//вычисляем скорости и давления во всех точках, а не только в тех, которые вне профиля
		for (size_t i = 0; i < initialPoints.size(); i++)
		{
			addvtx.r() = initialPoints[i];
			additionalWake->vtx.push_back(addvtx);
		}

		for (size_t i = 0; i < historyPoints.size(); i++)
		{
			addvtx.r() = historyPoints[i];
			additionalWake->vtx.push_back(addvtx);
		}

		velConvWake.resize(additionalWake->vtx.size(), { 0.0, 0.0 });
	}


	int addWSize = (int)additionalWake->vtx.size();
	MPI_Bcast(&addWSize, 1, MPI_INT, 0, W.getParallel().commWork);

	velConvBou.resize(W.getNumberOfBoundary());
	for (size_t i = 0; i < W.getNumberOfBoundary(); i++)
		velConvBou[i].resize(additionalWake->vtx.size(), { 0.0, 0.0 });



/*
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONV_TOWAKE))
	{
		W.getNonConstVelocity().GPUCalcConvVeloToSetOfPoints(*additionalWake, velConvWake, domainRadius);
		for (size_t i = 0; i < W.getNumberOfBoundary(); i++)
		{			
			W.getNonConstBoundary(i).GPUGetConvVelocityToSetOfPointsFromVirtualVortexes(*additionalWake, velConvBou[i]);
		}
	}
#else
*/
	{	
		/// \todo Сделать GPU-версию вычисления скорости и давления
		W.getNonConstVelocity().CalcConvVeloToSetOfPoints(*additionalWake, velConvWake, domainRadius);

		for (size_t i = 0; i < W.getNumberOfBoundary(); i++)
		{			
			W.getNonConstBoundary(i).GetConvVelocityToSetOfPointsFromVirtualVortexes(*additionalWake, velConvBou[i]);
		}
	}
//#endif


	if (W.getParallel().myidWork == 0)
	{
		velocity.assign(addWSize, W.getPassport().physicalProperties.V0());
		pressure.resize(addWSize);

		for (size_t i = 0; i < additionalWake->vtx.size(); i++)
			velocity[i] += velConvWake[i];

		for (size_t bou = 0; bou < velConvBou.size(); bou++)
		for (size_t j = 0; j < additionalWake->vtx.size(); j++)
			velocity[j] += velConvBou[bou][j];
	}

	int id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(addWSize);

	/// \todo Для успокоения добавить WakeSyncronize(), но после обучения его проверке hash

	std::vector<Point2D> locVelocity;
	locVelocity.resize(par.myLen);
	MPI_Scatterv(velocity.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, locVelocity.data(), par.myLen, Point2D::mpiPoint2D, 0, W.getParallel().commWork);

	std::vector<double> locPressure;
	locPressure.resize(par.myLen, 0.0);

	WakeDataBase locAdditionalWake(W);
	locAdditionalWake.vtx.resize(par.myLen);
	MPI_Scatterv(additionalWake->vtx.data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, locAdditionalWake.vtx.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	
	W.getNonConstVelocity().wakeVortexesParams.convVelo.resize(W.getWake().vtx.size());
	MPI_Bcast(W.getNonConstVelocity().wakeVortexesParams.convVelo.data(), (int)(W.getVelocity().wakeVortexesParams.convVelo.size()), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
		
	for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
	{
		W.getNonConstBoundary(bou).sheets.FreeSheetSynchronize();	
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

#pragma omp parallel for default(none) private(alpha, dri, Vi, vi, dst2eps, cPan) shared(P0, dt, eps2, V0, /*gamAverPan,*/ id, locAdditionalWake, locPressure, locVelocity, par, std::cout) 
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

			dst2eps = VMlib::boundDenom(dri.length2(), eps2);
		
			vi = IDPI * W.getWake().vtx[j].g() / dst2eps * dri.kcross();
			locPressure[locI] += vi * Vi; //3
		}

		for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
		for (int j = 0; j < W.getBoundary(bou).virtualWake.vtx.size(); j++)
		{
			dri = pt - W.getBoundary(bou).virtualWake.vtx[j].r();
			
			dst2eps = VMlib::boundDenom(dri.length2(), eps2);

			Vi = W.getBoundary(bou).virtualWake.vecHalfGamma[j];
			vi = IDPI * W.getBoundary(bou).virtualWake.vtx[j].g() / dst2eps * dri.kcross();

			locPressure[locI] +=  vi * Vi; //4
		}
			
		for (int bou = 0; bou < W.getNumberOfBoundary(); bou++)
		{
			const Point2D& rcm = W.getAirfoil(bou).rcm;

			for (int j = 0; j < W.getAirfoil(bou).getNumberOfPanels(); j++)
			{
				cPan = 0.5 * (W.getAirfoil(bou).getR(j) + W.getAirfoil(bou).getR(j + 1));
				alpha = atan2((cPan - pt) ^ (rcm - pt),   (cPan - pt)*(rcm - pt));
					
				/// \todo Пока используем только средние значения свободного слоя на панелях
				locPressure[locI] += IDPI * alpha *	W.getBoundary(bou).sheets.freeVortexSheet(j, 0) * W.getAirfoil(bou).len[j] / dt;

				locPressure[locI] -= IDPI * alpha * W.getAirfoil(bou).gammaThrough[j] / dt;
			}
		}//for bou
	}//for i

	MPI_Gatherv(locPressure.data(), par.myLen, MPI_DOUBLE, pressure.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);


	if (W.getParallel().myidWork == 0)
	{
		std::string fname = VMlib::fileNameStep("VelPres", W.getPassport().timeDiscretizationProperties.nameLength, step, "vtk");
		
		std::ofstream outfile;

		VMlib::CreateDirectory(dir, "velPres");
		outfile.open(dir + "velPres/" + fname);

		outfile << "# vtk DataFile Version 2.0" << std::endl;
		outfile << "VM2D VTK result: " << (dir + "velPres/" + fname).c_str() << " saved " << VMlib::CurrentDataTime() << std::endl;
		outfile << "ASCII" << std::endl;
		outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
		outfile << "POINTS " << additionalWake->vtx.size() << " float" << std::endl;

		for (size_t i = 0; i < additionalWake->vtx.size(); i++)
		{
			double xi = (additionalWake->vtx[i].r())[0];
			double yi = (additionalWake->vtx[i].r())[1];
			outfile << xi << " " << yi << " " << "0.0" << std::endl;
		}//for i

		outfile << "CELLS " << additionalWake->vtx.size() << " " << 2 * additionalWake->vtx.size() << std::endl;
		for (size_t i = 0; i < additionalWake->vtx.size(); ++i)
			outfile << "1 " << i << std::endl;

		outfile << "CELL_TYPES " << additionalWake->vtx.size() << std::endl;
		for (size_t i = 0; i < additionalWake->vtx.size(); ++i)
			outfile << "1" << std::endl;

		outfile << std::endl;
		outfile << "POINT_DATA " << additionalWake->vtx.size() << std::endl;

		outfile << "VECTORS V float" << std::endl;
		for (size_t i = 0; i < additionalWake->vtx.size(); i++)
		{
			outfile << velocity[i][0] << " " << velocity[i][1] << " 0.0" << std::endl;
		}//for i

		outfile << std::endl;

		outfile << "SCALARS P float 1" << std::endl;
		outfile << "LOOKUP_TABLE default" << std::endl;

		for (size_t i = 0; i < additionalWake->vtx.size(); i++)
		{
			outfile << pressure[i] << std::endl;
		}//for i

		outfile.close();

		//Вывод в csv-файлы
		for (int q = 0; q < historyPoints.size(); ++q)
		{
			std::string VPFileNameCsv;
			VPFileNameCsv = W.getPassport().dir + "velPres/" + VMlib::fileNameStep("VP-atPoint-", 2, q, "csv");

			std::ofstream VPFileCsv(VPFileNameCsv.c_str(), std::ios::app);

			VPFileCsv << W.getPassport().physicalProperties.getCurrTime() << ","
				<< velocity[initialPoints.size() + q][0] << ","
				<< velocity[initialPoints.size() + q][1] << ","
				<< pressure[initialPoints.size() + q] << std::endl;
			VPFileCsv.close();
		}

	}


	time.second = omp_get_wtime();

}//CalcSaveVP()


