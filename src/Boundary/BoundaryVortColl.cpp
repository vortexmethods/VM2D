/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.3    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/09/26     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: BoundaryVortColl.cpp                                             |
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
\brief Файл кода с описанием класса BoundaryVortColl
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/

#include "BoundaryVortColl.h"
#include "World2D.h"

//Конструктор
BoundaryVortColl::BoundaryVortColl(const World2D& W_, size_t numberInPassport_)
: Boundary(W_, numberInPassport_, 1)
{
	size_t np = afl.np;
	const std::vector<Point2D>& CC = afl.r;

	//задаем контрольные точки (точки коллокации)
	for (size_t i = 0; i < np; ++i)
		KK.push_back(0.5*(CC[i] + CC[i + 1]));
}//BoundaryVortColl(...)


//Вычисляет скос на точку R от точки X
Point2D BoundaryVortColl::Skos(const Point2D& R, const Point2D& X)
{
	double dst2 = dist2(R, X);
	Point2D res = { -(R[1] - X[1]), (R[0] - X[0]) };
	res *= IDPI / dst2;
	return res;
}//Skos(...)


//Генерация блока матрицы
void BoundaryVortColl::FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol)
{
	size_t np = afl.np;
	
	for (size_t i = 0; i < np; ++i)
	{
		lactCol(i) = 1.0;
		lastLine(i) = 1.0;
	}

	for (size_t i = 0; i < np; ++i)
	for (size_t j = 0; j < np; ++j)
	{
		matr(i, j) = Skos(KK[i], KK[j]) * afl.tau[i];
	}

	for (size_t i = 0; i < np; ++i)
	{
		//(afl.tau[i] ^ afl.nrm[i]) для учета внешней нормали
		matr(i, i) = (0.5 * (afl.tau[i] ^ afl.nrm[i]) )/ afl.len[i];
	}	
}//FillMatrixSelf(...)


//Генерация вектора влияния вихревого следа на профиль
void BoundaryVortColl::GetWakeInfluence(std::vector<double>& wakeVelo) const
{
	size_t np = afl.np;
	int id = W.getParallel().myidWork;

	parProp par = W.getParallel().SplitMPI(np);

	std::vector<double> locVeloWake;
	locVeloWake.resize(par.myLen);

	//локальные переменные для цикла
	double velI = 0.0;
	double tempVel = 0.0;
	double dst2 = 0.0;

#pragma omp parallel for default(none) shared(locVeloWake, par) private(velI, tempVel, dst2)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI = 0.0;

		const Point2D& posI = KK[par.myDisp + i];
		const Point2D& nrm = afl.nrm[par.myDisp + i];

		Point2D posJ;
		double gamJ;

		for (size_t j = 0; j < W.getWake().vtx.size(); ++j)
		{
			posJ = W.getWake().vtx[j].r();
			gamJ = W.getWake().vtx[j].g();

			dst2 = std::max(dist2(posI, posJ), W.getPassport().wakeDiscretizationProperties.eps2); //Сглаживание
			tempVel = nrm * (posI - posJ);
			tempVel *= (gamJ / dst2);
			velI += tempVel;
		}

		velI *= IDPI;
		locVeloWake[i] = velI;
	}

	if (id == 0)
		wakeVelo.resize(np); 

	MPI_Gatherv(locVeloWake.data(), par.myLen, MPI_DOUBLE, wakeVelo.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);
}//GetWakeInfluence(...)


//Генерация вектора влияния вихревого следа на профиль
#if defined(USE_CUDA)
void BoundaryVortColl::GPUGetWakeInfluence(std::vector<double>& wakeVelo) const
{
	std::cout << "GPU-computing is not implemented now! CPU-computing is active" << std::endl;
	GetWakeInfluence(wakeVelo);
}
#endif

//Вычисление скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле
void BoundaryVortColl::GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const
{
	std::vector<Point2D> selfVelo;

	size_t np = afl.np;
		
	int id = W.getParallel().myidWork;

	parProp par = W.getParallel().SplitMPI(points.size());

	std::vector<Point2D> locVelo;
	locVelo.resize(par.myLen);

	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2 = 0.0;
	double cft = IDPI;
	
#pragma omp parallel for default(none) shared(locVelo, cft, points, par) private(velI, tempVel, dst2)
	for (int i = 0; i < par.myLen; ++i)
	{	
		velI.toZero(); 

		const Point2D& posI = points[par.myDisp + i].r();

		Point2D posJ;
		double gamJ;
		
		for (size_t j = 0; j < afl.np; ++j)
		{
			posJ = KK[j];
			gamJ = sheets.freeVortexSheet[j][0] * afl.len[j];

			dst2 = std::max(dist2(posI, posJ), W.getPassport().wakeDiscretizationProperties.eps2); //Сглаживать надо!!!
			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2);
			velI += tempVel;
		}

		velI *= cft;
		locVelo[i] = velI;
	}

	if (id == 0)
	{
		selfVelo.resize(points.size());
	}

	MPI_Gatherv(locVelo.data(), par.myLen, Point2D::mpiPoint2D, selfVelo.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);

	if (id == 0)
		for (size_t i = 0; i < velo.size(); ++i)
			velo[i] += selfVelo[i];
}//GetVelocityToSetOfPoints(...)

//Заполнение в правой части влияния набегающего потока, следа и присоединенных слоев, действующих от самого себя
void BoundaryVortColl::FillRhs(const Point2D& V0, Eigen::VectorXd& rhs, double* lastRhs, bool move, bool deform)
{
	size_t np = afl.np;
	int id = W.getParallel().myidWork;

	std::vector<double> wakeVelo;

	GetWakeInfluence(wakeVelo);

	Point2D r, xi, dr;
	double drdr;

	if (id == 0)
	{
		for (size_t i = 0; i < np; ++i)
		{
			rhs(i) = -(V0*afl.tau[i] - afl.v[i] * afl.tau[i]) - wakeVelo[i];
			for (size_t j = 0; j < np; ++j)
			if (i != j)
			{
				r = 0.5 * (afl.r[i] + afl.r[i + 1]);
				xi = 0.5 * (afl.r[j] + afl.r[j + 1]);
				dr = r - xi;
				drdr = dr.length2();

				rhs(i) += -afl.tau[i] * (IDPI * sheets.attachedVortexSheet[j][0] * afl.len[j] * Skos(r, xi));
				rhs(i) += -afl.tau[i] * (IDPI * sheets.attachedSourceSheet[j][0] * afl.len[j] / drdr * dr);
			}//if (i != j)
			/// \todo 0.5 или 1.0???
			//	(afl.tau[i] ^ afl.nrm[i]) чтобы учесть внешнюю нормаль
			rhs(i) += -0.5 * sheets.attachedVortexSheet[i][0] * (afl.tau[i] ^ afl.nrm[i]);
		}//for i
	}// if (id == 0)

	*lastRhs = 0.0;
	
	for (size_t q = 0; q < afl.gammaThrough.size(); ++q)
		//for each (double g in afl.gammaThrough)
	{
		*lastRhs += afl.gammaThrough[q];
		//*lastRhs += g;
	}

}//FillRhs(...)

//Заполнение в правой части влияния присоединенных слоев, действующих на один профиль от другого
void BoundaryVortColl::FillRhsFromOther(const Airfoil& otherAirfoil, Eigen::VectorXd& rhs)
{
	Point2D r, xi, dr;
	double drdr;

	for (size_t i = 0; i < afl.np; ++i)
	{
		for (size_t j = 0; j < otherAirfoil.np; j++)
		{
			r = 0.5 * (afl.r[i] + afl.r[i + 1]);
			xi = 0.5 * (otherAirfoil.r[j] + otherAirfoil.r[j + 1]);
			dr = r - xi;
			drdr = dr.length2();

			rhs(i) += -afl.tau[i] * (IDPI * sheets.attachedVortexSheet[j][0] * otherAirfoil.len[j] * Skos(r, xi));
			rhs(i) += -afl.tau[i] * (IDPI * sheets.attachedSourceSheet[j][0] * otherAirfoil.len[j] / drdr * dr);
		}//for j
	}//for i
}//FillRhsFromOther(...)


//Возврат размерности вектора решения 
size_t BoundaryVortColl::GetUnknownsSize() const
{
	return afl.np;
}//GetUnknownsSize()


//Пересчет решения на интенсивность вихревого слоя и на рождаемые вихри на конкретном профиле
void BoundaryVortColl::SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol)
{
	for (size_t j = 0; j < afl.np; ++j)
		sheets.freeVortexSheet[j][0] = sol(j) / afl.len[j];

	//"закольцовываем"
	sheets.freeVortexSheet[afl.np][0] = sol(0) / afl.len[0];
	
	
	Vortex2D virtVort;
	Point2D midNorm;

	double delta = W.getPassport().wakeDiscretizationProperties.delta;
	//double delta = 0.5;

	
	//Сбрасываем с начала 1-й панели:
	midNorm = (afl.nrm[0] + afl.nrm[afl.np - 1]).unit(delta);
	virtVort.r() = afl.r[0] + midNorm;
	virtVort.g() = 0.5 * (sheets.freeVortexSheet[0][0] * afl.len[0] + sheets.freeVortexSheet[afl.np - 1][0] * afl.len[afl.np - 1]);
	virtualWake.vtx.push_back(virtVort);

	for (size_t i = 1; i < afl.np; ++i)
	{
		midNorm = (afl.nrm[i] + afl.nrm[i - 1]).unit(delta);
		virtVort.r() = afl.r[i] + midNorm;

		/// \todo сделать формирование присоединенных вихрей и источников
		virtVort.g() = 0.5 * (sheets.freeVortexSheet[i][0] * afl.len[i] + sheets.freeVortexSheet[i-1][0] * afl.len[i-1]);

		virtualWake.vtx.push_back(virtVort);
	}
}//SolutionToFreeVortexSheetAndVirtualVortex(...)

void BoundaryVortColl::ComputeAttachedSheetsIntensity()
{
	for (size_t i = 0; i < sheets.attachedVortexSheet.size(); ++i)
	{
		sheets.attachedVortexSheet[i][0] = afl.v[i] * afl.tau[i];
		sheets.attachedSourceSheet[i][0] = afl.v[i] * afl.nrm[i];
	}
}//ComputeAttachedSheetsIntensity()