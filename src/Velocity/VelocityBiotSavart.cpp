/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2017/12/01     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina       |
*-----------------------------------------------------------------------------*
| File name: VelocityBiotSavart.cpp                                           |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.	                                      |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.	                                                          |
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
\version 1.0
\date 1 декабря 2017 г.
*/


#include "VelocityBiotSavart.h"


//Вычисление конвективных скоростей в заданном наборе точек
void VelocityBiotSavart::CalcConvVeloToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo, std::vector<double>& domainRadius)
{
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

	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2eps, dst2;	

	
#pragma omp parallel for default(none) shared(locConvVelo, locDomRadius, locPoints, cft, par) private(tempVel, velI, dst2, dst2eps)
	for (int i = 0; i < par.myLen; ++i)
	{
		double ee2[3] = { 10000.0, 10000.0, 10000.0 };
		if (wake.vtx.size() == 0)
		{
			ee2[0] = 0.0;
			ee2[1] = 0.0;
			ee2[2] = 0.0;
		}
		
		velI.toZero();
		
		const Point2D& posI = locPoints[i].r();

		for (size_t j = 0; j < wake.vtx.size(); ++j)
		{
			const Point2D& posJ = wake.vtx[j].r();
			const double& gamJ = wake.vtx[j].g();

			tempVel.toZero();

			dst2 = dist2(posI, posJ);

			//Вычисление радиуса домена
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


			dst2eps = std::max(dst2, wake.param.eps2);
			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2eps);
			velI += tempVel;
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
}//CalcConvVeloToSetOfPoints(...)


//Вычисление диффузионных скоростей в заданном наборе точек
void VelocityBiotSavart::CalcDiffVeloToSetOfPoints(const std::vector<Vortex2D>& points, const std::vector<double>& domainRadius, const std::vector<Vortex2D>& vortices, std::vector<Point2D>& velo)
{
	std::vector<Point2D> selfVelo;

	const int& id = parallel.myidWork;

	parProp par = parallel.SplitMPI(points.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(points).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, parallel.commWork);

	std::vector<Point2D> locDiffVelo;
	locDiffVelo.resize(par.myLen);

	//Локальные переменные для цикла
	Point2D velI;
	Point2D I2;
	Point2D Rij;
	double I1, rij, expr;



#pragma omp parallel for default(none) shared(locDiffVelo, domainRadius, locPoints, vortices, par) private(velI, I1, I2, Rij, rij, expr)
	for (int i = 0; i < par.myLen; ++i)
	{
		const Vortex2D& vtxI = locPoints[i];
		velI.toZero();

		I2.toZero();
		I1 = 0.0;

		for (size_t j = 0; j < vortices.size(); ++j)
		{
			const Vortex2D& vtxJ = vortices[j];

			Rij = vtxI.r() - vtxJ.r();
			rij = Rij.length();

			if (rij > 1e-10)
			{
				expr = exp(-rij / domainRadius[i + par.myDisp]);
				I2 += (vtxJ.g()* expr / rij) * Rij;
				I1 += vtxJ.g()*expr;
			}//if (rij>1e-6)
		}//for j

		if (fabs(I1) > 1e-10)
			velI = I2* (1.0 / (I1 * domainRadius[i + par.myDisp]));

		locDiffVelo[i] = velI;
	} // for r

	//std::ostringstream sss;
	//sss << "velo_";
	//std::ofstream veloFile(sss.str());
	//for (int i = 0; i < domainRadius.size(); ++i)
	//	veloFile << domainRadius[i] << std::endl;
	//veloFile.close();

	if (id == 0)
		selfVelo.resize(points.size(), { 0.0, 0.0 });

	MPI_Gatherv(locDiffVelo.data(), par.myLen, Point2D::mpiPoint2D, selfVelo.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);


	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
}//CalcDiffVeloToSetOfPoints(...)