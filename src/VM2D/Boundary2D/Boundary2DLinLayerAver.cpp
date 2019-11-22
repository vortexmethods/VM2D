/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Boundary2DLinLayerAver.cpp                                       |
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
\brief Файл кода с описанием класса BoundaryLinLayerAver
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.5
\date 20 февраля 2019 г.
*/



#include "Boundary2DLinLayerAver.h"

#include "Airfoil2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;


//Пересчет решения на интенсивность вихревого слоя //Lin
void BoundaryLinLayerAver::SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol)
{
	Vortex2D virtVort;
	Point2D midNorm;

	double delta = W.getPassport().wakeDiscretizationProperties.delta;

	int nVortPerPan = W.getPassport().wakeDiscretizationProperties.minVortexPerPanel;

	//Очистка и резервирование памяти
	virtualWake.vecHalfGamma.clear();
	virtualWake.vecHalfGamma.reserve(afl.getNumberOfPanels() * nVortPerPan);

	//Очистка и резервирование памяти
	virtualWake.aflPan.clear();
	virtualWake.aflPan.reserve(afl.getNumberOfPanels() * nVortPerPan);

	//Очистка и резервирование памяти
	virtualWake.vtx.reserve(afl.getNumberOfPanels() * nVortPerPan);

	double maxG = W.getPassport().wakeDiscretizationProperties.maxGamma;

	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		midNorm = afl.nrm[i] * delta;

		//длина участков, на которые разбивается панель для сброса вихрей
		//определяется через наибольшее значение решения на профиле, т.е. в крайней точке
		double lenSmallPan = maxG / (sol(i) + 0.5 * fabs(sol(afl.getNumberOfPanels() + i)) * afl.len[i]);

		size_t NEWnVortPerPan = (size_t)std::max((int)std::ceil(afl.len[i] / lenSmallPan), nVortPerPan);

		Point2D dr = 1.0 / NEWnVortPerPan * (afl.getR(i + 1) - afl.getR(i));

		for (int j = 0; j < NEWnVortPerPan; j++)
		{
			virtVort.r() = afl.getR(i) + dr * (j * 1.0 + 0.5) + midNorm;
			virtVort.g() = sol(i) + sol(afl.getNumberOfPanels() + i) * (-0.5 * afl.len[i] + lenSmallPan * (j * 1.0 + 0.5));


			virtualWake.vtx.push_back(virtVort);

			// пока как в константной схеме //todolin
			virtualWake.vecHalfGamma.push_back(0.5 * sol(i)  * afl.tau[i]);
			virtualWake.aflPan.push_back({ numberInPassport, i });
		}
	}


	for (size_t j = 0; j < afl.getNumberOfPanels(); ++j)
	{
		sheets.freeVortexSheet(j, 0) = sol(j);
		sheets.freeVortexSheet(j, 1) = sol(afl.getNumberOfPanels() + j);
	}

}//SolutionToFreeVortexSheetAndVirtualVortex(...)



//Генерация блока матрицы //Lin
void BoundaryLinLayerAver::FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol)
{
	size_t np = afl.getNumberOfPanels();

	std::vector<double> res(4, 0.0);

	for (size_t i = 0; i < np; ++i)
	{
		lactCol(i) = 1.0;
		lastLine(i) = afl.len[i];
	}

	for (size_t i = 0; i < np; ++i)
	for (size_t j = 0; j < np; ++j)
	{
		res = afl.getA(2, i, afl, j);
		matr(i, j) = res[0];
		matr(i, np + j) = res[1];
		matr(np + i, j) = res[2];
		matr(np + i, np + j) = res[3];
	}
}//FillMatrixSelf(...)

//Генерация блока матрицы //Lin
void BoundaryLinLayerAver::FillIQSelf(std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ)
{
	afl.calcIQ(1, afl, IQ);
}//FillIQSelf(...)


//Генерация блока матрицы влияния от другого профиля того же типа // Lin
void BoundaryLinLayerAver::FillMatrixFromOther(const Boundary& otherBoundary, Eigen::MatrixXd& matr)
{
	size_t np = afl.getNumberOfPanels();
	size_t npOther = otherBoundary.afl.getNumberOfPanels();

	std::vector<double> res(4, 0.0);

	for (size_t i = 0; i < np; ++i)
	for (size_t j = 0; j < npOther; ++j)
	{
		res = afl.getA(2, i, otherBoundary.afl, j);

		matr(i, j) = res[0];
		matr(i, np + j) = res[1];
		matr(np + i, j) = res[2];
		matr(np + i, np + j) = res[3];
	}
}//FillMatrixFromOther(...)

void BoundaryLinLayerAver::FillIQFromOther(const Boundary& otherBoundary, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ)
{
	afl.calcIQ(1, otherBoundary.afl, IQ);
}//FillIQFromOther(...)


//Вычисление скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле //Lin
void BoundaryLinLayerAver::GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const
{
	std::vector<Point2D> selfVelo;

	size_t np = afl.getNumberOfPanels();

	int id = W.getParallel().myidWork;

	VMlib::parProp par = W.getParallel().SplitMPI(points.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(points).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	std::vector<Point2D> locVelo;
	locVelo.resize(par.myLen);
	  
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel0, tempVel1;

#pragma omp parallel for default(none) shared(locVelo, locPoints, par) private(velI, tempVel0, tempVel1)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI.toZero();

		const Point2D& posI = locPoints[i].r();

		for (size_t j = 0; j < afl.getNumberOfPanels(); ++j)
		{
			const Point2D& posJ0 = afl.getR(j);
			const Point2D& posJ1 = afl.getR(j + 1);
			const Point2D& tau = afl.tau[j];

			double gamJ0 = sheets.freeVortexSheet(j, 0);
			double gamJ1 = sheets.freeVortexSheet(j, 1);

			Point2D s = posI - posJ0;
			Point2D p = posI - posJ1;
			Point2D dj = posJ1 - posJ0;
			Point2D u1 = 0.5 / dj.length() * VMlib::Omega(p + s, tau, tau);

			double alpha = VMlib::Alpha(p, s);

			double lambda = VMlib::Lambda(p, s);

			tempVel0 = (alpha* tau + lambda * tau.kcross()) * gamJ0;
			tempVel1 = (-alpha * u1.kcross() + lambda * u1 - IDPI * tau).kcross() * gamJ1;

			velI += tempVel0 + tempVel1;
		} //for j

		velI *= IDPI;
		locVelo[i] = velI;
	}

	if (id == 0)
		selfVelo.resize(points.size());

	MPI_Gatherv(locVelo.data(), par.myLen, Point2D::mpiPoint2D, selfVelo.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);

	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
}//GetVelocityToSetOfPoints(...)


//Вычисление скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле как от виртуальных вихрей
void BoundaryLinLayerAver::GetConvVelocityToSetOfPointsFromVirtualVortexes(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const
{
	std::vector<Point2D> selfVelo;

	size_t np = afl.getNumberOfPanels();

	int id = W.getParallel().myidWork;

	VMlib::parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size());



	//синхронизация свободного вихревого слоя
	std::vector<double> freeVortexSheetGamma;
	int sz;

	if (id == 0)
	{
		sz = (int)sheets.getSheetSize();
		freeVortexSheetGamma.resize(sz, 0.0);
		for (size_t j = 0; j < sheets.getSheetSize(); ++j)
			freeVortexSheetGamma[j] = sheets.freeVortexSheet(j, 0);
	}
	MPI_Bcast(&sz, 1, MPI_INT, 0, W.getParallel().commWork);
	if (id != 0)
		freeVortexSheetGamma.resize(sz, 0.0);

	MPI_Bcast(freeVortexSheetGamma.data(), (int)sz, MPI_DOUBLE, 0, W.getParallel().commWork);
	//свободный вихревой слой синхронизирован



	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	double cft = IDPI;

	std::vector<Point2D> locConvVelo;
	locConvVelo.resize(par.myLen);

#pragma warning (push)
#pragma warning (disable: 4101)
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	//double dst2eps, dst2;
#pragma warning (pop)

#pragma omp parallel for default(none) shared(locConvVelo, locPoints, cft, par, freeVortexSheetGamma, std::cout) private(velI, tempVel/*, dst2, dst2eps*/)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI.toZero();

		const Point2D& posI = locPoints[i].r();

		/// \todo Тут надо разобраться, как должно быть...
		/// \todo сделать  if(move || deform)
		for (size_t j = 0; j < sheets.getSheetSize(); j++)
		{
			Point2D dj = afl.getR(j + 1) - afl.getR(j);
			Point2D tauj = dj.unit();

			Point2D s = posI - afl.getR(j);
			Point2D p = posI - afl.getR(j + 1);

			double a = VMlib::Alpha(p, s);

			double lambda = VMlib::Lambda(p, s);

			velI += (freeVortexSheetGamma[j] * (-a * tauj.kcross() + lambda * tauj)).kcross();
			velI += (sheets.attachedVortexSheet(j, 0) * (-a * tauj.kcross() + lambda * tauj)).kcross();
			velI += (sheets.attachedSourceSheet(j, 0) * (-a * tauj.kcross() + lambda * tauj));
		}//for j

		velI *= cft;
		locConvVelo[i] = velI;
	}//for i

	if (id == 0)
		selfVelo.resize(pointsDb.vtx.size());

	MPI_Gatherv(locConvVelo.data(), par.myLen, Point2D::mpiPoint2D, selfVelo.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);

	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
}//GetVelocityToSetOfPointsFromVirtualVortexes(...)

#if defined(USE_CUDA)
void BoundaryLinLayerAver::GPUGetConvVelocityToSetOfPointsFromVirtualVortexes(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const
{
	const size_t npt = pointsDb.vtx.size();
	double*& dev_ptr_pt = pointsDb.devVtxPtr;
	const size_t npnl = afl.getNumberOfPanels(); //virtualWake.vtx.size();

	double*& dev_ptr_r = afl.devRPtr;
	double*& dev_ptr_freeVortexSheet = afl.devFreeVortexSheetPtr;
	double*& dev_ptr_attachedVortexSheet = afl.devAttachedVortexSheetPtr;
	double*& dev_ptr_attachedSourceSheet = afl.devAttachedSourceSheetPtr;

	std::vector<Point2D>& Vel = velo;
	std::vector<Point2D>& locvel = pointsDb.tmpVel;
	double*& dev_ptr_vel = pointsDb.devVelPtr;
	double eps2 = W.getPassport().wakeDiscretizationProperties.eps2;

	const int& id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(npt, true);

	//Явная синхронизация слоев не нужна, т.к. она выполняется в Gpu::RefreshAfls() 

	double tCUDASTART = 0.0, tCUDAEND = 0.0;

	tCUDASTART = omp_get_wtime();

	if (npt > 0)
	{
		cuCalculateConvVeloWakeFromVirtual(par.myDisp, par.myLen, dev_ptr_pt, npnl, dev_ptr_r, dev_ptr_freeVortexSheet, dev_ptr_attachedVortexSheet, dev_ptr_attachedSourceSheet, dev_ptr_vel, eps2);

		W.getCuda().CopyMemFromDev<double, 2>(par.myLen, dev_ptr_vel, (double*)&locvel[0]);

		std::vector<Point2D> newV;
		if (id == 0)
			newV.resize(Vel.size());

		MPI_Gatherv(locvel.data(), par.myLen, Point2D::mpiPoint2D, newV.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);

		if (id == 0)
		for (size_t q = 0; q < Vel.size(); ++q)
			Vel[q] += newV[q];
	}

	//tCUDAEND = omp_get_wtime();

	//W.getInfo('t') << "CONV_VIRT_GPU: " << (tCUDAEND - tCUDASTART) << std::endl;
}
//GPUGetVelocityToSetOfPointsFromVirtualVortexes(...)
#endif


void BoundaryLinLayerAver::GetConvVelocityAtVirtualVortexes(std::vector<Point2D>& velo) const
{
	const int& id = W.getParallel().myidWork;
	std::vector<Point2D>& Vel = velo;

	//Скорости виртуальных вихрей
	if (id == 0)
	{
#pragma omp parallel for default(none) shared(Vel)		
		for (int i = 0; i < Vel.size(); ++i)
			Vel[i] = afl.getV(virtualWake.aflPan[i].second) \
			+ virtualWake.vecHalfGamma[i] \
			- W.getPassport().physicalProperties.V0();
	}
}
//GetConvVelocityAtVirtualVortexes(...)



//Заполнение в правой части влияния присоединенных слоев, действующих на один профиль от другого //Lin
//!!!!!! Вызывать после того, как будет вызвана соответствующая FillMatrixFromOther
void BoundaryLinLayerAver::FillRhsFromOther(const Airfoil& otherAirfoil, Eigen::VectorXd& rhs, size_t currentRow, size_t currentCol)
{
	//Переделать с учетом линейных и сделать так, чтобы еще раз не пересчитывались все интегралы в функции getInfAttFromOther0
	std::vector<double> attOtherVelo;

	afl.getInfAttFromOther1(attOtherVelo, otherAirfoil, currentRow, currentCol);

	for (size_t i = 0; i < GetUnknownsSize(); ++i)
	{
		rhs(i) += attOtherVelo[i];
	}

}//FillRhsFromOther(...)


//Вычисление интенсивностей присоединенного вихревого слоя и присоединенного слоя источников
void BoundaryLinLayerAver::ComputeAttachedSheetsIntensity()
{
	for (size_t i = 0; i < sheets.getSheetSize(); ++i)
	{
		oldSheets.attachedVortexSheet(i, 0) = sheets.attachedVortexSheet(i, 0);
		oldSheets.attachedSourceSheet(i, 0) = sheets.attachedVortexSheet(i, 0);
	}

	for (size_t i = 0; i < sheets.getSheetSize(); ++i)
	{
		sheets.attachedVortexSheet(i, 0) = afl.getV(i) * afl.tau[i];
		sheets.attachedSourceSheet(i, 0) = afl.getV(i) * afl.nrm[i];
	}
}//ComputeAttachedSheetsIntensity()


void BoundaryLinLayerAver::GetInfluenceFromVorticesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const
{
	double& velI = wakeRhs[0];
	double& velILin = wakeRhs[1];

	const Point2D& posI0 = afl.getR(panel);
	const Point2D& posI1 = afl.getR(panel + 1);
	Point2D di = posI1 - posI0;
	const Point2D& taui = afl.tau[panel];

	for (size_t it = 0; it != count; ++it)
	{
		const Vortex2D& vt = ptr[it];
		const Point2D& posJ = vt.r();
		const double& gamJ = vt.g();

		Point2D s = posJ - posI0;
		Point2D p = posJ - posI1;

		Point2D u1 = 0.5 / di.length() * VMlib::Omega(p + s, taui, taui);

		double alpha = VMlib::Alpha(p, s);
		double lambda = VMlib::Lambda(p, s);

		velI -= gamJ * alpha;

		velILin -= gamJ * (alpha *u1 * taui + lambda * u1 * (-taui.kcross()));
	}
}//GetInfluenceFromVorticesToRectPanel(...)


//Вычисляет влияния части подряд идущих источников в области течения на прямолинейную панель для правой части
void BoundaryLinLayerAver::GetInfluenceFromSourcesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const
{
	double& velI = wakeRhs[0];
	double& velILin = wakeRhs[1];

	const Point2D& posI0 = afl.getR(panel);
	const Point2D& posI1 = afl.getR(panel + 1);
	Point2D di = posI1 - posI0;
	const Point2D& taui = afl.tau[panel];

	for (size_t it = 0; it != count; ++it)
	{
		const Vortex2D& vt = ptr[it];
		const Point2D& posJ = vt.r();
		const double& gamJ = vt.g();

		Point2D s = posJ - posI0;
		Point2D p = posJ - posI1;

		Point2D u1 = 0.5 / di.length() * VMlib::Omega(p + s, taui, taui);

		double alpha = VMlib::Alpha(p, s);
		double lambda = VMlib::Lambda(p, s);

		velI -= gamJ * lambda;

		velILin -= gamJ * (alpha *u1 *  (taui.kcross()) + lambda * u1 * taui - 1.0);
	}

}//GetInfluenceFromSourcesToRectPanel(...)