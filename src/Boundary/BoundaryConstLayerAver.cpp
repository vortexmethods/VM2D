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

	parallel.SplitMPI(np);

	std::vector<double> locVeloWake;
	locVeloWake.resize(parallel.myLen);

	//локальные переменные для цикла
	double velI = 0.0;
	double tempVel = 0.0;

#pragma omp parallel for default(none) shared(locVeloWake, id) private(velI, tempVel)
	for (int i = 0; i < parallel.myLen; ++i)
	{
		velI = 0.0;

		const Point2D& posI0 = CC[parallel.myDisp + i];
		const Point2D& posI1 = CC[parallel.myDisp + i + 1];
		const Point2D& tau = afl.tau[parallel.myDisp + i];

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

		velI *= IDPI / afl.len[parallel.myDisp + i];
		locVeloWake[i] = velI;
	}

	if (id == 0)
		wakeVelo.resize(np);

	MPI_Gatherv(locVeloWake.data(), parallel.myLen, MPI_DOUBLE, wakeVelo.data(), parallel.len.data(), parallel.disp.data(), MPI_DOUBLE, 0, parallel.commWork);
}//GetWakeInfluence(...)

//Переделанная в соответствии с дисс. Моревой
//void BoundaryConstLayerAver::GetWakeInfluence(std::vector<double>& wakeVelo) const
//{
//	size_t np = afl.np;
//	int id = parallel.myidWork;
//
//	parallel.SplitMPI(np);
//
//	std::vector<double> locVeloWake;
//	locVeloWake.resize(parallel.len[id]);
//
//	//локальные переменные для цикла
//	double velI = 0.0;
//	double tempVel = 0.0;
//
//#pragma omp parallel for default(none) shared(locVeloWake, id) private(velI, tempVel)
//	for (int i = 0; i < parallel.len[id]; ++i)
//	{
//		velI = 0.0;
//
//		const Point2D& posI0 = CC[parallel.disp[id] + i];
//		const Point2D& posI1 = CC[parallel.disp[id] + i + 1];
//		const Point2D& tau = afl.tau[parallel.disp[id] + i];
//		Point2D d = CC[parallel.disp[id] + i + 1] - CC[parallel.disp[id] + i];
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
//		velI *= IDPI / (afl.len[parallel.disp[id] + i] * afl.len[parallel.disp[id] + i]);
//		locVeloWake[i] = velI;
//	}
//
//	if (id == 0)
//		wakeVelo.resize(np);
//
//
//	MPI_Gatherv(locVeloWake.data(), parallel.len[id], MPI_DOUBLE, wakeVelo.data(), parallel.len.data(), parallel.disp.data(), MPI_DOUBLE, 0, parallel.commWork);
//}//GetWakeInfluence(...)

//Вычисление скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле
void BoundaryConstLayerAver::GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const
{
	std::vector<Point2D> selfVelo;

	size_t np = afl.np;

	int id = parallel.myidWork;

	parallel.SplitMPI(points.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(parallel.myLen);

	MPI_Scatterv(points.data(), parallel.len.data(), parallel.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), parallel.myLen, Vortex2D::mpiVortex2D, 0, parallel.commWork);

	std::vector<Point2D> locVelo;
	locVelo.resize(parallel.myLen);

	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;

#pragma omp parallel for default(none) shared(locVelo, id, locPoints) private(velI, tempVel)
	for (int i = 0; i < parallel.myLen; ++i)
	{
		velI = { 0.0, 0.0 };

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

	MPI_Gatherv(locVelo.data(), parallel.myLen, Point2D::mpiPoint2D, selfVelo.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);

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

	parallel.SplitMPI(points.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(parallel.myLen);

	MPI_Scatterv(points.data(), parallel.len.data(), parallel.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), parallel.myLen, Vortex2D::mpiVortex2D, 0, parallel.commWork);

	double cft = IDPI;

	std::vector<Point2D> locConvVelo;
	locConvVelo.resize(parallel.myLen);
	
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2eps, dst2;


#pragma omp parallel for default(none) shared(locConvVelo, locPoints, id, cft) private(velI, tempVel, dst2, dst2eps)
	for (int i = 0; i < parallel.myLen; ++i)
	{
		velI = { 0.0, 0.0 };

		const Point2D& posI = locPoints[i].r();

		for (size_t j = 0; j < virtualWake.size(); ++j)
		{
			const Point2D& posJ = virtualWake[j].r();
			const double& gamJ = virtualWake[j].g();

			tempVel = { 0.0, 0.0 };

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

	MPI_Gatherv(locConvVelo.data(), parallel.myLen, Point2D::mpiPoint2D, selfVelo.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);

	parallel.BCastAllLenDisp();

	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
}//GetVelocityToSetOfPointsFromVirtualVortexes(...)

//Заполнение правой части
void BoundaryConstLayerAver::FillRhs(const Point2D& V0, Eigen::VectorXd& rhs, double* lastRhs)
{
	std::vector<double> wakeVelo;
	GetWakeInfluence(wakeVelo);

	if (parallel.myidWork == 0)
	for (size_t i = 0; i < afl.np; ++i)
		rhs(i) = -(V0*afl.tau[i]) - wakeVelo[i];
	
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


