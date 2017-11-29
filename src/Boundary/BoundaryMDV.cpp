/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2017/12/01     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina       |
*-----------------------------------------------------------------------------*
| File name: BoundaryMDV.cpp                                                  |
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
\brief Файл кода с описанием класса BoundaryMDV
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include "BoundaryMDV.h"


//Конструктор
BoundaryMDV::BoundaryMDV(const Passport& passport_, const Airfoil& afl_, const std::vector<std::unique_ptr<Boundary>>& allBoundary_, const Wake& wake_, const Parallel& parallel_)
: Boundary(passport_, afl_, allBoundary_, 1, wake_, parallel_)
{
	size_t np = afl.np;
	const std::vector<Point2D>& CC = afl.r;

	//задаем контрольные точки (точки коллокации)
	for (size_t i = 0; i < np; ++i)
		KK.push_back(0.5*(CC[i] + CC[i + 1]));
}//BoundaryMDV(...)


//Вычисляет скос на точку R от точки X
Point2D BoundaryMDV::Skos(const Point2D& R, const Point2D& X)
{
	double dst2 = dist2(R, X);
	Point2D res = { -(R[1] - X[1]), (R[0] - X[0]) };
	res *= IDPI / dst2;
	return res;
}//Skos(...)


//Генерация блока матрицы
void BoundaryMDV::FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol)
{
	size_t np = afl.np;

	for (size_t i = 0; i < np; ++i)
	{
		lactCol(i) = afl.len[i];
		lastLine(i) = 1.0;
	}

	for (size_t i = 0; i < np; ++i)
	for (size_t j = 0; j < np; ++j)
	{
		matr(i, j) = Skos(KK[i], CC[j]) * afl.nrm[i];
	}

}//FillMatrixSelf(...)


//Генерация вектора влияния вихревого следа на профиль
void BoundaryMDV::GetWakeInfluence(std::vector<double>& wakeVelo) const
{
	size_t np = afl.np;
	int id = parallel.myidWork;

	parProp par = parallel.SplitMPI(np);

	std::vector<double> locVeloWake;
	locVeloWake.resize(par.myLen);

	//локальные переменные для цикла
	double velI = 0.0;
	double tempVel = 0.0;
	double dst2 = 0.0;

#pragma omp parallel for default(none) shared(locVeloWake) private(velI, tempVel, dst2, par)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI = 0.0;

		const Point2D& posI = KK[par.myDisp + i];
		const Point2D& tau = afl.tau[par.myDisp + i];

		for (size_t j = 0; j < wake.vtx.size(); ++j)
		{
			const Point2D& posJ = wake.vtx[j].r();
			const double& gamJ = wake.vtx[j].g();

			dst2 = std::max(dist2(posI, posJ), wake.param.eps2); //Сглаживать надо!!!
			tempVel = tau * (posI - posJ);
			tempVel *= (gamJ / dst2);
			velI -= tempVel;
		}

		velI *= IDPI;
		locVeloWake[i] = velI;
	}

	if (id == 0)
		wakeVelo.resize(np);
		

	MPI_Gatherv(locVeloWake.data(), par.myLen, MPI_DOUBLE, wakeVelo.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, parallel.commWork);
}//GetWakeInfluence(...)


//Вычисление скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле
void BoundaryMDV::GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const
{
	std::vector<Point2D> selfVelo;

	size_t np = afl.np;

	int id = parallel.myidWork;

	parProp par = parallel.SplitMPI(points.size());

	std::vector<Point2D> locVelo;
	locVelo.resize(par.myLen);

	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2 = 0.0;

#pragma omp parallel for default(none) shared(locVelo, points, par) private(velI, tempVel, dst2)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI.toZero();

		const Point2D& posI = points[par.myDisp + i].r();

		for (size_t j = 0; j < afl.np; ++j)
		{
			const Point2D& posJ = KK[j];
			double gamJ = sheets.freeVortexSheet[j][0] * afl.len[j];

			dst2 = std::max(dist2(posI, posJ), wake.param.eps2); //Сглаживать надо!!!
			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2);
			velI += tempVel;
		}

		velI *= IDPI;
		locVelo[i] = velI;
	}

	if (id == 0)
		selfVelo.resize(points.size());

	MPI_Gatherv(locVelo.data(), par.myLen, Point2D::mpiPoint2D, selfVelo.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);

	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
}//GetVelocityToSetOfPoints(...)


//Заполнение правой части
void BoundaryMDV::FillRhs(const Point2D& V0, Eigen::VectorXd& rhs, double* lastRhs)
{
	size_t np = afl.np;
	int id = parallel.myidWork;

	std::vector<double> wakeVelo;

	GetWakeInfluence(wakeVelo);

	if (id == 0)
	{
		for (size_t i = 0; i < np; ++i)
			rhs(i) = -(V0*afl.nrm[i]) - wakeVelo[i];
	}

	*lastRhs = 0.0;

	for (size_t q = 0; q < afl.gammaThrough.size(); ++q)
	//for each (double g in afl.gammaThrough)
	{
		*lastRhs += afl.gammaThrough[q];
		//*lastRhs += g;
	}
}//FillRhs(...)


//Возврат размерности вектора решения 
int BoundaryMDV::GetUnknownsSize() const
{
	return afl.np;
}//GetUnknownsSize()


//Пересчет решения на интенсивность вихревого слоя и на рождаемые вихри на конкретном профиле
void BoundaryMDV::SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol)
{
	for (size_t j = 0; j < (afl.np)-1; ++j)
		sheets.freeVortexSheet[j][0] = 0.5*( sol(j) + sol(j+1) ) / afl.len[j];
	
	sheets.freeVortexSheet[(afl.np) - 1][0] = 0.5*(sol( (afl.np) - 1 ) + sol(0)) / afl.len[(afl.np) - 1];

	//"закольцовываем"
	sheets.freeVortexSheet[afl.np][0] = 0.5*( sol(0) + sol(1) ) / afl.len[0];


	Vortex2D virtVort;
	Point2D midNorm;

	/// \todo delta в паспорт
	double delta = passport.wakeDiscretizationProperties.delta;
	//double delta = 0.5;

	
	//Сбрасываем с начала 1-й панели:
	midNorm = (afl.nrm[0] + afl.nrm[afl.np - 1]).unit(delta);
	virtVort.r() = afl.r[0] + midNorm;
	virtVort.g() = sol(0);
	virtualWake.push_back(virtVort);

	for (size_t i = 1; i < afl.np; ++i)
	{
		midNorm = (afl.nrm[i] + afl.nrm[i - 1]).unit(delta);
		virtVort.r() = afl.r[i] + midNorm;

		/// \todo сделать формирование присоединенных вихрей и источников
		virtVort.g() = sol(i);

		virtualWake.push_back(virtVort);
	}
}//SolutionToFreeVortexSheetAndVirtualVortex(...)



