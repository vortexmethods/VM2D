/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2017/12/01     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina       |
*-----------------------------------------------------------------------------*
| File name: BoundaryConstLayerAver.cpp                                       |
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
\brief Файл кода с описанием класса BoundaryConstLayerAver
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#include "BoundaryConstLayerAver.h"


//Возврат размерности вектора решения 
int BoundaryConstLayerAver::GetUnknownsSize() const
{
	return afl.np;
}//GetUnknownsSize()


//Пересчет решения на интенсивность вихревого слоя
void BoundaryConstLayerAver::SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol)
{
	//TODO: Пока нет слоев, сразу сбрасываются вихри

	Vortex2D virtVort;
	Point2D midNorm;

	double delta = passport.wakeDiscretizationProperties.delta;

	//Сбрасываем с начала 1-й панели:
	//TODO: Убрали среднюю нормаль
	//midNorm = (afl.nrm[0] + afl.nrm[afl.np - 1]).unit(delta);
	//midNorm = afl.nrm[0] * delta;
	//virtVort.r() = afl.r[0] + midNorm;

	//TODO: Переделала сброс на сброс с конца панели
	//virtVort.g() = 0.5*(sol(0)*afl.len[0] + sol(afl.np - 1)*afl.len[afl.np - 1]);
	//virtualWake.push_back(virtVort);

	//for (size_t i = 1; i < afl.np; ++i)
	//{
	//	//TODO: Убрали среднюю нормаль
	//	//midNorm = (afl.nrm[i] + afl.nrm[i - 1]).unit(delta);
	//	midNorm = afl.nrm[i] * delta;
	//	virtVort.r() = afl.r[i] + midNorm;

	//	/// \todo сделать формирование присоединенных вихрей и источников
	//	virtVort.g() = 0.5*(sol(i)*afl.len[i] + sol(i - 1)*afl.len[i - 1]);;

	//	virtualWake.push_back(virtVort);
	//}

	
	for (size_t i = 0; i < afl.np - 1; ++i)
	{
		//POLARA
		//Убрали среднюю нормаль
		//midNorm = (afl.nrm[i] + afl.nrm[i - 1]).unit(delta);
		midNorm = afl.nrm[i] * delta;
		virtVort.r() = afl.r[i + 1] + midNorm;

		/// \todo сделать формирование присоединенных вихрей и источников
		virtVort.g() = 0.5*(sol(i)*afl.len[i] + sol(i + 1)*afl.len[i + 1]);;

		virtualWake.push_back(virtVort);
	}

	midNorm = afl.nrm[afl.np - 1] * delta;
	virtVort.r() = afl.r[0] + midNorm;
	virtVort.g() = 0.5*(sol(0)*afl.len[0] + sol(afl.np - 1)*afl.len[afl.np - 1]);
	virtualWake.push_back(virtVort);


	for (size_t j = 0; j < afl.np; ++j)
		sheets.freeVortexSheet[j][0] = sol(j);

	//"закольцовываем"
	sheets.freeVortexSheet.push_back(sheets.freeVortexSheet[0]);

/*
	Vortex2D virtVort;
	Point2D midNorm;

	/// \todo delta в паспорт
	double delta = passport.wakeDiscretizationProperties.delta;
	//double delta = 0.5;


	//Сбрасываем с начала 1-й панели:
	midNorm = (afl.nrm[0] + afl.nrm[afl.np - 1]).unit(delta);
	virtVort.r() = afl.r[0] + midNorm;
	virtVort.g() = 0.5*(sol(0)*afl.len[0] + sol(afl.np - 1)*afl.len[afl.np - 1]);
	virtualWake.push_back(virtVort);

	for (size_t i = 1; i < afl.np; ++i)
	{
		midNorm = (afl.nrm[i] + afl.nrm[i - 1]).unit(delta);
		virtVort.r() = afl.r[i] + midNorm;

		/// \todo сделать формирование присоединенных вихрей и источников
		virtVort.g() = 0.5*(sol(i)*afl.len[i] + sol(i - 1)*afl.len[i - 1]);;

		virtualWake.push_back(virtVort);
	}*/
	
}//SolutionToFreeVortexSheetAndVirtualVortex(...)


//Генерация блока матрицы
void BoundaryConstLayerAver::FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol)
{
	size_t np = afl.np;

	//Panel vectors
	std::vector<Point2D> dd;
	for (size_t i = 0; i < np; ++i)
		dd.push_back(afl.tau[i] * afl.len[i]);

	for (size_t i = 0; i < np; ++i)
	{
		lactCol(i) = 1.0;
		lastLine(i) = afl.len[i];
	}

	//auxillary scalars
	numvector<double, 3> alpha, lambda;

	//auxillary vectors
	Point2D p1, s1, p2, s2, i00;
	numvector<Point2D, 3> v;

	for (size_t i = 0; i < np; ++i)
	for (size_t j = 0; j < np; ++j)
	{
		if (i != j)
		{
			const Point2D& di = dd[i];
			const Point2D& dj = dd[j];

			const Point2D& taui = afl.tau[i];
			const Point2D& tauj = afl.tau[j];

			p1 = CC[i + 1] - CC[j + 1];
			s1 = CC[i + 1] - CC[j];
			p2 = CC[i] - CC[j + 1];
			s2 = CC[i] - CC[j];

			alpha = { \
				afl.isAfter(j, i) ? 0.0 : Alpha(s2, s1), \
				Alpha(s2, p1), \
				afl.isAfter(i, j) ? 0.0 : Alpha(p1, p2) \
			};

			lambda = { \
				afl.isAfter(j, i) ? 0.0 : Lambda(s2, s1), \
				Lambda(s2, p1), \
				afl.isAfter(i, j) ? 0.0 : Lambda(p1, p2) \
			};

			v = { Omega(s1, taui, tauj), -Omega(di, taui, tauj), Omega(p2, taui, tauj) };

			i00 = IDPI / afl.len[i] * ((alpha[0] * v[0] + alpha[1] * v[1] + alpha[2] * v[2]) + (lambda[0] * v[0] + lambda[1] * v[1] + lambda[2] * v[2]).kcross());

			//i00 = IDPI / afl.len[i] * (-(alpha[0] * v[0] + alpha[1] * v[1] + alpha[2] * v[2]).kcross() + (lambda[0] * v[0] + lambda[1] * v[1] + lambda[2] * v[2]));


			matr(i, j) = i00 * afl.tau[i];
		}
	}

	for (size_t i = 0; i < np; ++i)
	{
		matr(i, i) = -0.5;
	}	
}//FillMatrixSelf(...)




//Генерация блока матрицы влияния от другого профиля того же типа
void BoundaryConstLayerAver::FillMatrixFromOther(const Boundary& otherBoundary, Eigen::MatrixXd& matr)
{
	size_t np = afl.np;
	size_t npOther = otherBoundary.afl.np;

	//Panel vectors
	std::vector<Point2D> dd;
	for (size_t i = 0; i < np; ++i)
		dd.push_back(afl.tau[i] * afl.len[i]);

	std::vector<Point2D> ddOther;
	for (size_t j = 0; j < npOther; ++j)
		ddOther.push_back(otherBoundary.afl.tau[j] * otherBoundary.afl.len[j]);

	//auxillary scalars
	numvector<double, 3> alpha, lambda;

	//auxillary vectors
	Point2D p1, s1, p2, s2, i00;
	numvector<Point2D, 3> v;

	for (size_t i = 0; i < np; ++i)
	for (size_t j = 0; j < npOther; ++j)
	{
		const Point2D& di = dd[i];
		const Point2D& dj = ddOther[j];

		const Point2D& taui = afl.tau[i];
		const Point2D& tauj = otherBoundary.afl.tau[j];

		p1 = CC[i + 1] - otherBoundary.CC[j + 1];
		s1 = CC[i + 1] - otherBoundary.CC[j];
		p2 = CC[i] - otherBoundary.CC[j + 1];
		s2 = CC[i] - otherBoundary.CC[j];

		alpha = { \
			Alpha(s2, s1), \
			Alpha(s2, p1), \
			Alpha(p1, p2) \
		};

		lambda = { \
			Lambda(s2, s1), \
			Lambda(s2, p1), \
			Lambda(p1, p2) \
		};

		v = { Omega(s1, taui, tauj), -Omega(di, taui, tauj), Omega(p2, taui, tauj) };

		i00 = IDPI / afl.len[i] * ((alpha[0] * v[0] + alpha[1] * v[1] + alpha[2] * v[2]) + (lambda[0] * v[0] + lambda[1] * v[1] + lambda[2] * v[2]).kcross());

		//i00 = IDPI / afl.len[i] * (-(alpha[0] * v[0] + alpha[1] * v[1] + alpha[2] * v[2]).kcross() + (lambda[0] * v[0] + lambda[1] * v[1] + lambda[2] * v[2]));


		matr(i, j) = i00 * afl.tau[i];
	}

}//FillMatrixFromOther(...)



//Генерация вектора влияния вихревого следа на профиль
void BoundaryConstLayerAver::GetWakeInfluence(std::vector<double>& wakeVelo) const
{
	size_t np = afl.np;
	int id = parallel.myidWork;
	parProp par = parallel.SplitMPI(np);

	std::vector<double> locVeloWake;
	locVeloWake.resize(par.myLen);

	//локальные переменные для цикла
	double velI = 0.0;
	double tempVel = 0.0;

#pragma omp parallel for default(none) shared(locVeloWake, par) private(velI, tempVel)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI = 0.0;

		const Point2D& posI0 = CC[par.myDisp + i];
		const Point2D& posI1 = CC[par.myDisp + i + 1];
		const Point2D& tau = afl.tau[par.myDisp + i];

		for (size_t j = 0; j < wake.vtx.size(); ++j)
		{
			const Point2D& posJ = wake.vtx[j].r();
			const double& gamJ = wake.vtx[j].g();

			Point2D s = posJ - posI0;
			Point2D p = posJ - posI1;

			double alpha = Alpha(p, s);
			double lambda = Lambda(p, s);

			tempVel = gamJ * alpha;
			velI -= tempVel;
		}

		velI *= IDPI / afl.len[par.myDisp + i];
		locVeloWake[i] = velI;
	}

	if (id == 0)
		wakeVelo.resize(np);

	MPI_Gatherv(locVeloWake.data(), par.myLen, MPI_DOUBLE, wakeVelo.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, parallel.commWork);
}//GetWakeInfluence(...)

//Переделанная в соответствии с дисс. Моревой
//void BoundaryConstLayerAver::GetWakeInfluence(std::vector<double>& wakeVelo) const
//{
//	size_t np = afl.np;
//	int id = parallel.myidWork;
//
//	parProp par = parallel.SplitMPI(np);
//
//	std::vector<double> locVeloWake;
//	locVeloWake.resize(par.myLen);
//
//	//локальные переменные для цикла
//	double velI = 0.0;
//	double tempVel = 0.0;
//
//#pragma omp parallel for default(none) shared(locVeloWake, id, par) private(velI, tempVel)
//	for (int i = 0; i < par.myLen; ++i)
//	{
//		velI = 0.0;
//
//		const Point2D& posI0 = CC[par.myDisp + i];
//		const Point2D& posI1 = CC[par.myDisp + i + 1];
//		const Point2D& tau = afl.tau[par.myDisp + i];
//		Point2D d = CC[par.myDisp + i + 1] - CC[par.myDisp + i];
//
//		for (size_t j = 0; j < wake.vtx.size(); ++j)
//		{
//			const Point2D& posJ = wake.vtx[j].r();
//			const double& gamJ = wake.vtx[j].g();
//
//			Point2D s = -posJ + posI1;
//			Point2D s0 = -posJ + posI0;
//
//			double z0 = cross3(d, s0);
//
//			double alpha = atan2(s0 * d, z0) - atan2(s * d, z0);
//
//			double lambda = Lambda(s0, s);
//
//			tempVel = tau * (alpha * d + lambda * d.kcross());
//			tempVel *= gamJ;
//			velI -= tempVel;
//		}
//
//		velI *= IDPI / (afl.len[par.myDisp + i] * afl.len[par.myDisp + i]);
//		locVeloWake[i] = velI;
//	}
//
//	if (id == 0)
//		wakeVelo.resize(np);
//
//
//	MPI_Gatherv(locVeloWake.data(), par.myLen, MPI_DOUBLE, wakeVelo.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, parallel.commWork);
//}//GetWakeInfluence(...)

//Вычисление скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле
void BoundaryConstLayerAver::GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const
{
	std::vector<Point2D> selfVelo;

	size_t np = afl.np;

	int id = parallel.myidWork;

	parProp par = parallel.SplitMPI(points.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(points).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, parallel.commWork);

	std::vector<Point2D> locVelo;
	locVelo.resize(par.myLen);

	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;

#pragma omp parallel for default(none) shared(locVelo, locPoints, par) private(velI, tempVel)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI.toZero();

		const Point2D& posI = locPoints[i].r();

		for (size_t j = 0; j < afl.np; ++j)
		{
			const Point2D& posJ0 = afl.r[j];
			const Point2D& posJ1 = afl.r[j + 1];
			const Point2D& tau = afl.tau[j];

			double gamJ = sheets.freeVortexSheet[j][0];
			
			Point2D s = posI - posJ0;
			Point2D p = posI - posJ1;

			double alpha = Alpha(p, s);
			
			double lambda = Lambda(p, s);
			
			tempVel = alpha* tau + lambda * tau.kcross();
			tempVel *= gamJ;
			velI += tempVel;
		} //for j

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


//Вычисление скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле как от виртуальных вихрей
void BoundaryConstLayerAver::GetConvVelocityToSetOfPointsFromVirtualVortexes(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const
{
	std::vector<Point2D> selfVelo;

	size_t np = afl.np;

	int id = parallel.myidWork;

	parProp par = parallel.SplitMPI(points.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(points).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, parallel.commWork);

	double cft = IDPI;

	std::vector<Point2D> locConvVelo;
	locConvVelo.resize(par.myLen);
	
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2eps, dst2;


#pragma omp parallel for default(none) shared(locConvVelo, locPoints, cft, par) private(velI, tempVel, dst2, dst2eps)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI.toZero();

		const Point2D& posI = locPoints[i].r();

		for (size_t j = 0; j < virtualWake.size(); ++j)
		{
			const Point2D& posJ = virtualWake[j].r();
			const double& gamJ = virtualWake[j].g();

			tempVel.toZero();

			dst2 = dist2(posI, posJ);

			dst2eps = std::max(dst2, wake.param.eps2);
			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2eps);
			velI += tempVel;
		}

		velI *= cft;
		locConvVelo[i] = velI;
	}

	if (id == 0)
		selfVelo.resize(points.size());

	MPI_Gatherv(locConvVelo.data(), par.myLen, Point2D::mpiPoint2D, selfVelo.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);

	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
}//GetVelocityToSetOfPointsFromVirtualVortexes(...)

//Заполнение в правой части влияния набегающего потока, следа и присоединенных слоев, действующих от самого себя
void BoundaryConstLayerAver::FillRhs(const Point2D& V0, Eigen::VectorXd& rhs, double* lastRhs, bool move, bool deform)
{
	std::vector<double> wakeVelo;
	
	GetWakeInfluence(wakeVelo);

	if (parallel.myidWork == 0)
	for (size_t i = 0; i < afl.np; ++i)
	{
		rhs(i) = -afl.tau[i] * (V0 -afl.v[i] ) - wakeVelo[i];


		//влияние присоединенных слоев от самого себя
		if (move || deform)
		for (size_t j = 0; j < afl.np; j++)
		{
			if (i != j)
			{
				Point2D di = afl.r[i + 1] - afl.r[i];
				Point2D dj = afl.r[j + 1] - afl.r[j];
				Point2D s1 = afl.r[i + 1] - afl.r[j];
				Point2D s2 = afl.r[i] - afl.r[j];
				Point2D p1 = afl.r[i + 1] - afl.r[j + 1];
				Point2D p2 = afl.r[i] - afl.r[j + 1];

				double a1 = Alpha(s2, s1);
				double a2 = Alpha(s2, p1);
				double a3 = Alpha(p1, p2);

				double lambda1 = Lambda(s1, s2);
				double lambda2 = Lambda(p1, s2);
				double lambda3 = Lambda(p2, p1);

				Point2D v1 = Omega(s1, afl.tau[i], afl.tau[j]);
				Point2D v2 = -Omega(di, afl.tau[i], afl.tau[j]);
				Point2D v3 = Omega(p2, afl.tau[i], afl.tau[j]);

				if ((i == j + 1) || ((i == 0) && (j == afl.np-1)))
				{
					a3 = 0.0; lambda3 = 0.0;
				}
				else if ((j == i + 1) || ((j == 0) && (i == afl.np - 1)))
				{
					a1 = 0.0; lambda1 = 0.0;
				}
				
				rhs(i) += -IDPI / afl.len[i] * sheets.attachedVortexSheet[j][0] * afl.tau[i] * (((-(a1 * v1 + a2 * v2 + a3 * v3).kcross()) + lambda1 * v1 + lambda2 * v2 + lambda3 * v3).kcross());
				rhs(i) += -IDPI / afl.len[i] * sheets.attachedSourceSheet[j][0] * afl.tau[i] * (((-(a1 * v1 + a2 * v2 + a3 * v3).kcross()) + lambda1 * v1 + lambda2 * v2 + lambda3 * v3));
				
			}//if (i != j)
		}//for j
		rhs(i) += -0.5 * sheets.attachedVortexSheet[i][0];

	}//for i
	
	if (parallel.myidWork == 0)
	{
		*lastRhs = 0.0;
		
		for (size_t q = 0; q < afl.gammaThrough.size(); ++q)
			//for each (double g in afl.gammaThrough)
		{
			*lastRhs += afl.gammaThrough[q];
			//*lastRhs += g;
		}
	}
}//FillRhs(...)

//Заполнение в правой части влияния присоединенных слоев, действующих на один профиль от другого
void BoundaryConstLayerAver::FillRhsFromOther(const Airfoil& otherAirfoil, Eigen::VectorXd& rhs)
{
	for (size_t i = 0; i < afl.np; ++i)
	{
		for (size_t j = 0; j < otherAirfoil.np; j++)
		{
			Point2D di = afl.r[i + 1] - afl.r[i];
			Point2D dj = otherAirfoil.r[j + 1] - otherAirfoil.r[j];
			Point2D s1 = afl.r[i + 1] - otherAirfoil.r[j];
			Point2D s2 = afl.r[i] - otherAirfoil.r[j];
			Point2D p1 = afl.r[i + 1] - otherAirfoil.r[j + 1];
			Point2D p2 = afl.r[i] - otherAirfoil.r[j + 1];

			double a1 = Alpha(s2, s1);
			double a2 = Alpha(s2, p1);
			double a3 = Alpha(p1, p2);

			double lambda1 = Lambda(s1, s2);
			double lambda2 = Lambda(p1, s2);
			double lambda3 = Lambda(p2, p1);

			Point2D v1 = Omega(s1, afl.tau[i], otherAirfoil.tau[j]);
			Point2D v2 = -Omega(di, afl.tau[i], otherAirfoil.tau[j]);
			Point2D v3 = Omega(p2, afl.tau[i], otherAirfoil.tau[j]);

			rhs(i) += -IDPI / afl.len[i] * sheets.attachedVortexSheet[j][0] * afl.tau[i] * (((-(a1 * v1 + a2 * v2 + a3 * v3).kcross()) + lambda1 * v1 + lambda2 * v2 + lambda3 * v3).kcross());
			rhs(i) += -IDPI / afl.len[i] * sheets.attachedSourceSheet[j][0] * afl.tau[i] * (((-(a1 * v1 + a2 * v2 + a3 * v3).kcross()) + lambda1 * v1 + lambda2 * v2 + lambda3 * v3));
		}//for j
	}//for i
}//FillRhsFromOther(...)


//Вычисление интенсивностей присоединенного вихревого слоя и присоединенного слоя источников
void BoundaryConstLayerAver::ComputeAttachedSheetsIntensity()
{
	for (size_t i = 0; i < sheets.attachedVortexSheet.size(); ++i)
	{
		sheets.attachedVortexSheet[i][0] = afl.v[i] * afl.tau[i];
		sheets.attachedSourceSheet[i][0] = afl.v[i] * afl.nrm[i];
	}
}//ComputeAttachedSheetsIntensity()


