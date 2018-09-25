/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
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
\version 1.1
\date 2 апреля 2018 г.
*/


#include "VelocityBiotSavart.h"

/// Конструктор
VelocityBiotSavart::VelocityBiotSavart(const Parallel& parallel_, gpu& cuda_, const Wake& wake_, const std::vector<std::unique_ptr<Boundary>>& boundary_) :
	Velocity(parallel_, cuda_, wake_, boundary_)
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
void VelocityBiotSavart::CalcConvVeloToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo, std::vector<double>& domainRadius)
{
	double tCPUSTART, tCPUEND;
	
	tCPUSTART = omp_get_wtime();
	std::vector<Point2D> selfVelo;
	
	const int& id = parallel.myidWork;

	parProp par = parallel.SplitMPI(points.size(), true);

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	//Заполнение "своей" части массива скоростей
	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(points).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		         locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, parallel.commWork);

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
	
#pragma omp parallel for default(none) shared(locConvVelo, locDomRadius, locPoints, cft, par) private(tempVel, velI, dst2, dst2eps)
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

		for (size_t j = 0; j < wake.vtx.size(); ++j)
		{
			const Point2D& posJ = wake.vtx[j].r();
			const double& gamJ = wake.vtx[j].g();

			tempVel.toZero();

			dst2 = dist2(posI, posJ);

			//Модифицируем массив квадратов расстояний до ближайших вихрей из wake
			ModifyE2(ee2, dst2);

			dst2eps = std::max(dst2, wake.param.eps2);
			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2eps);
			velI += tempVel;
		}
		
		for (size_t s = 0; s < boundary.size(); ++s)
		{
			for (size_t j = 0; j < boundary[s]->virtualWake.size(); ++j)
			{
				const Point2D& posJ = boundary[s]->virtualWake[j].r();
				const double& gamJ = boundary[s]->virtualWake[j].g();

				dst2 = dist2(posI, posJ);

				//Модифицируем массив квадратов расстояний до ближайших вихрей из virtualWake
				ModifyE2(ee2, dst2);
			}
		}

		velI *= cft;
		locConvVelo[i] = velI;

		locDomRadius[i] = std::max(sqrt((ee2[0] + ee2[1] + ee2[2]) / 3.0), 2.0*wake.param.epscol);
	}


	if (id == 0)
		selfVelo.resize(points.size());


	//if (id == 0)
	//это нужно всем, т.к. ниже стоит Allgatherv
	domainRadius.resize(par.totalLen);

	MPI_Gatherv(locConvVelo.data(), par.myLen, Point2D::mpiPoint2D, selfVelo.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
	


	//parallel.BCastAllLenDisp();
	MPI_Allgatherv(locDomRadius.data(), par.myLen, MPI_DOUBLE, domainRadius.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, parallel.commWork);

	if (id == 0)
		for (size_t i = 0; i < velo.size(); ++i)
			velo[i] += selfVelo[i];
	
	tCPUEND = omp_get_wtime();
	//std::cout << "CONV_CPU: " << tCPUEND - tCPUSTART << std::endl;

}//CalcConvVeloToSetOfPoints(...)


//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
void VelocityBiotSavart::CalcDiffVeloI1I2ToSetOfPoints(const std::vector<Vortex2D>& points, const std::vector<double>& domainRadius, const std::vector<Vortex2D>& vortices, std::vector<double>& I1, std::vector<Point2D>& I2)
{
	double tCPUSTART, tCPUEND;

	tCPUSTART = omp_get_wtime();

	std::vector<double> selfI1;
	std::vector<Point2D> selfI2;

	const int& id = parallel.myidWork;

	parProp par = parallel.SplitMPI(points.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(points).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, parallel.commWork);

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

	
#pragma omp parallel for default(none) shared(locI1, locI2, domainRadius, locPoints, vortices, par) private(Rij, rij, expr, diffRadius, left, right, posJx)
	for (int i = 0; i < par.myLen; ++i)
	{
		const Vortex2D& vtxI = locPoints[i];

		diffRadius = 7.0 * domainRadius[i + par.myDisp];

		left = vtxI.r()[0] - diffRadius;
		right = vtxI.r()[0] + diffRadius;
		
		for (size_t j = 0; j < vortices.size(); ++j)
		{
			const Vortex2D& vtxJ = vortices[j];
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
		selfI1.resize(points.size(), 0.0);
		selfI2.resize(points.size(), { 0.0, 0.0 });
	}

	MPI_Gatherv(locI1.data(), par.myLen, MPI_DOUBLE,          selfI1.data(), par.len.data(), par.disp.data(), MPI_DOUBLE,          0, parallel.commWork);
	MPI_Gatherv(locI2.data(), par.myLen, Point2D::mpiPoint2D, selfI2.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);


	if (id == 0)
	for (size_t i = 0; i < I1.size(); ++i)
	{
		I1[i] += selfI1[i];
		I2[i] += selfI2[i];
	}

	tCPUEND = omp_get_wtime();
	//std::cout << "DIFF_CPU: " << tCPUEND - tCPUSTART << std::endl;
}//CalcDiffVeloI1I2ToSetOfPoints(...)