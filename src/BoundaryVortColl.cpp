/*!
\file
\brief Файл кода с описанием класса BoundaryVortColl
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "BoundaryVortColl.h"


//Конструктор
BoundaryVortColl::BoundaryVortColl(const std::unique_ptr<Airfoil>& afl_, const Wake& wake_, const Parallel& parallel_)
	: Boundary(afl_, 1, wake_, parallel_)
{
	size_t np = afl->np;
	const std::vector<Point2D>& CC = afl->r;

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
	size_t np = afl->np;
	
	//формируем векторы панелей
	std::vector<Point2D> dd;
	for (size_t i = 0; i < np; ++i)
		dd.push_back(afl->tau[i] * afl->len[i]);

	for (size_t i = 0; i < np; ++i)
	{
		lactCol(i) = 1.0;
		lastLine(i) = 1.0;
	}

	for (size_t i = 0; i < np; ++i)
	for (size_t j = 0; j < np; ++j)
	{
		matr(i, j) = Skos(KK[i], CC[j]) * afl->tau[i];
	}

	for (size_t i = 0; i < np; ++i)
	{
		matr(i, i) = -0.5 / afl->len[i];
	}	
}//FillMatrixSelf(...)


//Генерация вектора влияния вихревого следа на профиль
void BoundaryVortColl::GetWakeInfluence(std::vector<double>& wakeVelo) const
{
	size_t np = afl->np;
	int id = parallel.myidWork;

	parallel.SplitMPI(np);

	std::vector<double> locVeloWake;
	locVeloWake.resize(parallel.len[id]);

	//локальные переменные для цикла
	double velI = 0.0;
	double tempVel = 0.0;
	double dst2 = 0.0;

#pragma omp parallel for default(none) shared(locVeloWake, id) private(velI, tempVel, dst2)
	for (int i = 0; i < parallel.len[id]; ++i)
	{
		velI = 0.0;

		const Point2D& posI = KK[parallel.disp[id] + i];
		const Point2D& nrm = afl->nrm[parallel.disp[id] + i];

		for (int j = 0; j < wake.nv; ++j)
		{
			const Point2D& posJ = wake.vtx[j].r();
			const double& gamJ = wake.vtx[j].g();

			dst2 = std::max(dist2(posI, posJ), 1e-10); //Сглаживать не надо!!!
			tempVel = nrm * (posI - posJ);
			tempVel *= (gamJ / dst2);
			velI += tempVel;
		}

		velI *= IDPI;
		locVeloWake[i] = velI;
	}

	if (id == 0)
		wakeVelo.resize(np);

	MPI_Gatherv(locVeloWake.data(), parallel.len[id], MPI_DOUBLE, wakeVelo.data(), parallel.len.data(), parallel.disp.data(), MPI_DOUBLE, 0, parallel.commWork);
}//GetWakeInfluence(...)


//Вычисление скоростей вихрей в вихревом следе, вызываемых наличием завихренности и источников на профиле
void BoundaryVortColl::GetWakeVelocity(std::vector<Point2D>& wakeVelo, double dt) const
{
	size_t np = afl->np;
	double eps2 = sqr(wake.param.eps);
	
	int id = parallel.myidWork;

	parallel.SplitMPI(wake.nv);

	std::vector<Point2D> locVeloWake;
	locVeloWake.resize(parallel.len[id]);

	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2 = 0.0;
	double cft = IDPI*dt;
	
#pragma omp parallel for default(none) shared(locVeloWake, id, eps2, cft) private(velI, tempVel, dst2)
	for (int i = 0; i < parallel.len[id]; ++i)
	{	
		velI = { 0.0, 0.0 }; 

		const Point2D& posI = wake.vtx[parallel.disp[id] + i].r();
		
		for (size_t j = 0; j < afl->np; ++j)
		{
			const Point2D& posJ = KK[j];
			double gamJ = sheets.freeVortexSheet[j][0] * afl->len[j];

			dst2 = std::max(dist2(posI, posJ), eps2); //Сглаживать надо!!!
			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2);
			velI += tempVel;
		}

		velI *= cft;
		locVeloWake[i] = velI;
	}

	if (id == 0)
		wakeVelo.resize(wake.nv);

	MPI_Gatherv(locVeloWake.data(), parallel.len[id], Point2D::mpiPoint2D, wakeVelo.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
}//GetWakeVelocity(...)


//Заполнение правой части
void BoundaryVortColl::FillRhs(const Point2D& V0, Eigen::VectorXd& rhs)
{
	size_t np = afl->np;
	int id = parallel.myidWork;

	std::vector<double> wakeVelo;

	GetWakeInfluence(wakeVelo);

	if (id == 0)
	{
		for (size_t i = 0; i < np; ++i)
			rhs(i) = -(V0*afl->tau[i]) - wakeVelo[i];
	}
}//FillRhs(...)


//Возврат размерности вектора решения 
int BoundaryVortColl::GetUnknownsSize() const
{
	return afl->np;
}//GetUnknownsSize()


//Пересчет решения на интенсивность вихревого слоя
void BoundaryVortColl::SolutionToFreeVortexSheet(const Eigen::VectorXd& sol)
{
	for (size_t j = 0; j < afl->np; ++j)
		sheets.freeVortexSheet[j][0] = sol(j) / afl->len[j];
}//SolutionToFreeVortexSheet(...)