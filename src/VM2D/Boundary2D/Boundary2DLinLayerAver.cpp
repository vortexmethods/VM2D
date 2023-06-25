/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.11   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2022/08/07     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
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
\brief Файл кода с описанием класса BoundaryLinLayerAver
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.11
\date 07 августа 2022 г.
*/



#include "Boundary2DLinLayerAver.h"

#include "Airfoil2D.h"
#include "Airfoil2DCurv.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Tree2D.h"
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
	virtualWake.vtx.clear();
	virtualWake.vtx.reserve(afl.getNumberOfPanels() * nVortPerPan);

	//Очистка и резервирование памяти
	vortexBeginEnd.clear();
	vortexBeginEnd.reserve(afl.getNumberOfPanels());

	double maxG = W.getPassport().wakeDiscretizationProperties.maxGamma;

	std::pair<int, int> pair;

	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		midNorm = afl.nrm[i] * delta;

		//число участков, на которые разбивается панель для сброса вихрей
		//определяется через наибольшее значение решения на профиле, т.е. в крайней точке
		 

		pair.first = (int)virtualWake.vtx.size();

		double a = sol(i) - 0.5 * sol(afl.getNumberOfPanels() + i);
		double b = sol(i) + 0.5 * sol(afl.getNumberOfPanels() + i);

		if (fabs(a - b) < 1e-10)
		{
			size_t NEWnVortPerPan = (size_t)std::max((int)std::ceil(fabs(sol(i)*afl.len[i]) / maxG), nVortPerPan);
			Point2D dr = 1.0 / NEWnVortPerPan * (afl.getR(i + 1) - afl.getR(i));

			for (size_t j = 0; j < NEWnVortPerPan; ++j)
			{
				virtVort.r() = afl.getR(i) + dr * (j * 1.0 + 0.5) + midNorm;
				virtVort.g() = sol(i) * afl.len[i] / NEWnVortPerPan;
				virtualWake.vtx.push_back(virtVort);

				virtualWake.vecHalfGamma.push_back(0.5 * sol(i)  * afl.tau[i]);
				virtualWake.aflPan.push_back({ numberInPassport, i });
			}
		}
		else if (a * b >= 0.0)
		{
			size_t NEWnVortPerPan = (size_t)std::max((int)std::ceil(fabs(0.5*(a+b)*afl.len[i]) / maxG), nVortPerPan);
			std::vector<double> ds(NEWnVortPerPan+1, 0.0);
			
			for (size_t k = 1; k <= NEWnVortPerPan; ++k)
			{
				if ((a > 0) || (b > 0))
					ds[k] = afl.len[i] / (a - b) * (a - sqrt((b*b * k + a * a*(NEWnVortPerPan - k)) / NEWnVortPerPan));
				else
					ds[k] = afl.len[i] / (a - b) * (a + sqrt((b*b * k + a * a*(NEWnVortPerPan - k)) / NEWnVortPerPan));
			}

			for (size_t j = 0; j < NEWnVortPerPan; ++j)
			{
				virtVort.r() = afl.getR(i) + 0.5*(ds[j]+ds[j+1])*afl.tau[i] + midNorm;
				virtVort.g() = 0.5 * (a+b) * afl.len[i] / NEWnVortPerPan;
				virtualWake.vtx.push_back(virtVort);

				virtualWake.vecHalfGamma.push_back(0.5 * (a + 0.5*(ds[j] + ds[j + 1]) * (b-a)/ afl.len[i])  * afl.tau[i]);
				virtualWake.aflPan.push_back({ numberInPassport, i });
			}
		}
		else
		{
			double sast = -a * afl.len[i] / (b - a);

			//from 0 to sast
			{
				size_t NEWnVortPerPan = (size_t)std::max((int)std::ceil(fabs(0.5*(a + 0.0)*sast) / maxG), (int)std::ceil(nVortPerPan * sast / afl.len[i]));
				std::vector<double> ds(NEWnVortPerPan + 1, 0.0);

				for (size_t k = 1; k <= NEWnVortPerPan; ++k)
				{
					if (a > 0)
						ds[k] = sast / (a - 0.0) * (a - sqrt((a * a*(NEWnVortPerPan - k)) / NEWnVortPerPan));
					else
						ds[k] = sast / (a - 0.0) * (a + sqrt((a * a*(NEWnVortPerPan - k)) / NEWnVortPerPan));
				}

				for (size_t j = 0; j < NEWnVortPerPan; ++j)
				{
					virtVort.r() = afl.getR(i) + 0.5*(ds[j] + ds[j + 1])*afl.tau[i] + midNorm;
					virtVort.g() = 0.5 * (a + 0.0) * sast / NEWnVortPerPan;
					virtualWake.vtx.push_back(virtVort);

					virtualWake.vecHalfGamma.push_back(0.5 * (a + 0.5*(ds[j] + ds[j + 1]) * (0.0 - a) / sast)  * afl.tau[i]);
					virtualWake.aflPan.push_back({ numberInPassport, i });
				}
			}

			double sastast = afl.len[i] - sast;

			//from sats to len
			{
				size_t NEWnVortPerPan = (size_t)std::max((int)std::ceil(fabs(0.5*(0.0 + b)*sastast) / maxG), (int)std::ceil(nVortPerPan * sastast / afl.len[i]));
				std::vector<double> ds(NEWnVortPerPan + 1, 0.0);

				for (size_t k = 1; k <= NEWnVortPerPan; ++k)
				{
					if (b > 0)
						ds[k] = sastast / (0.0 - b) * (0.0 - sqrt((b*b * k) / NEWnVortPerPan));
					else
						ds[k] = sastast / (0.0 - b) * (0.0 + sqrt((b*b * k) / NEWnVortPerPan));
				}

				for (size_t j = 0; j < NEWnVortPerPan; ++j)
				{
					virtVort.r() = afl.getR(i) + sast * afl.tau[i] + 0.5*(ds[j] + ds[j + 1])*afl.tau[i] + midNorm;
					virtVort.g() = 0.5 * (0.0 + b) * sastast / NEWnVortPerPan;
					virtualWake.vtx.push_back(virtVort);

					virtualWake.vecHalfGamma.push_back(0.5 * (0.0 + 0.5*(ds[j] + ds[j + 1]) * (b - 0.0) / sastast)  * afl.tau[i]);
					virtualWake.aflPan.push_back({ numberInPassport, i });
				}

			}
		}

		pair.second = (int)virtualWake.vtx.size();
		vortexBeginEnd.push_back(pair);
	}


	for (size_t j = 0; j < afl.getNumberOfPanels(); ++j)
	{
		sheets.freeVortexSheet(j, 0) = sol(j);
		sheets.freeVortexSheet(j, 1) = sol(afl.getNumberOfPanels() + j);
	}

}//SolutionToFreeVortexSheetAndVirtualVortex(...)



//Генерация блока матрицы //Lin
void BoundaryLinLayerAver::FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lastCol)
{
	size_t np = afl.getNumberOfPanels();

	std::vector<double> res(4, 0.0);

	for (size_t i = 0; i < np; ++i)
	{
		lastCol(i) = 1.0;
		lastLine(i) = afl.len[i];
		lastCol(np + i) = 0.0;
		lastLine(np + i) = 0.0;
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
	afl.calcIQ(2, afl, IQ);
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
	afl.calcIQ(2, otherBoundary.afl, IQ);
}//FillIQFromOther(...)


//Вычисление скоростей в наборе точек, вызываемых наличием слоев вихрей и источников на профиле
void BoundaryLinLayerAver::CalcConvVelocityToSetOfPointsFromSheets(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const
{
	std::vector<Point2D> selfVelo;

	//size_t np = afl.getNumberOfPanels();

	int id = W.getParallel().myidWork;

	VMlib::parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size());


	sheets.FreeSheetSynchronize();

	////синхронизация свободного вихревого слоя
	//std::vector<double> freeVortexSheetGamma;
	//int sz;

	//if (id == 0)
	//{
	//	sz = (int)sheets.getSheetSize();
	//	freeVortexSheetGamma.resize(sz, 0.0);
	//	for (size_t j = 0; j < afl.getNumberOfPanels(); ++j)
	//	{
	//		freeVortexSheetGamma[j] = sheets.freeVortexSheet(j, 0);
	//		freeVortexSheetGamma[afl.getNumberOfPanels() + j] = sheets.freeVortexSheet(j, 1);
	//	}
	//}
	//MPI_Bcast(&sz, 1, MPI_INT, 0, W.getParallel().commWork);
	//if (id != 0)
	//	freeVortexSheetGamma.resize(sz, 0.0);

	//MPI_Bcast(freeVortexSheetGamma.data(), (int)sz, MPI_DOUBLE, 0, W.getParallel().commWork);
	////свободный вихревой слой синхронизирован



	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D/*, VM2D::MyAlloc<VMlib::Vortex2D>*/>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	std::vector<Point2D> locConvVelo;
	locConvVelo.resize(par.myLen);

#pragma warning (push)
#pragma warning (disable: 4101)
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	//double dst2eps, dst2;
#pragma warning (pop)

#pragma omp parallel for default(none) shared(locConvVelo, locPoints, par, std::cout) private(velI, tempVel/*, dst2, dst2eps*/)
	for (int i = 0; i < par.myLen; ++i)
	{
		velI.toZero();

		const Point2D& posI = locPoints[i].r();

		/// \todo Тут надо разобраться, как должно быть...
		/// \todo сделать  if(move || deform)
		for (size_t j = 0; j < afl.getNumberOfPanels(); ++j)
		{
			Point2D dj = afl.getR(j + 1) - afl.getR(j);
			Point2D tauj = dj.unit();

			Point2D s = posI - afl.getR(j);
			Point2D p = posI - afl.getR(j + 1);

			double a = VMlib::Alpha(p, s);

			double lambda = VMlib::Lambda(p, s);

			Point2D u1 = 0.5 / dj.length() * VMlib::Omega(p + s, tauj, tauj) ;

			Point2D skos0 = -a * tauj.kcross() + lambda * tauj;
			Point2D skos1 = -a * u1.kcross() + lambda * u1 - tauj;

			/// \todo почему не sheets.freeVortexSheet(j, 0)?  
			velI += sheets.freeVortexSheet(j, 0) * skos0.kcross() + sheets.freeVortexSheet(j, 1) * skos1.kcross();
			velI += sheets.attachedVortexSheet(j, 0) * skos0.kcross() + sheets.attachedVortexSheet(j, 1) * skos1.kcross();
			velI += sheets.attachedSourceSheet(j, 0) * skos0 + sheets.attachedSourceSheet(j, 1) * skos1;
		}//for j

		velI *= IDPI;
		locConvVelo[i] = velI;
	}//for i

	if (id == 0)
		selfVelo.resize(pointsDb.vtx.size());

	MPI_Gatherv(locConvVelo.data(), par.myLen, Point2D::mpiPoint2D, selfVelo.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);

	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
}//CalcConvVelocityToSetOfPointsFromSheets(...)

#if defined(USE_CUDA)
void BoundaryLinLayerAver::GPUCalcConvVelocityToSetOfPointsFromSheets(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const
{
	const size_t npt = pointsDb.vtx.size();
	double*& dev_ptr_pt = pointsDb.devVtxPtr;
	const size_t npnl = afl.getNumberOfPanels(); //virtualWake.vtx.size();

	double*& dev_ptr_r = afl.devRPtr;
	double*& dev_ptr_freeVortexSheet = afl.devFreeVortexSheetPtr;
	double*& dev_ptr_attachedVortexSheet = afl.devAttachedVortexSheetPtr;
	double*& dev_ptr_attachedSourceSheet = afl.devAttachedSourceSheetPtr;

	std::vector<Point2D>& Vel = velo;
	std::vector<Point2D> locvel(npt);
	double*& dev_ptr_vel = pointsDb.devVelPtr;
	double eps2 = W.getPassport().wakeDiscretizationProperties.eps2;

	const int& id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(npt, true);

	//Явная синхронизация слоев не нужна, т.к. она выполняется в Gpu::RefreshAfls() 
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
//GPUCalcConvVelocityToSetOfPointsFromSheets(...)
#endif


void BoundaryLinLayerAver::CalcConvVelocityAtVirtualVortexes(std::vector<Point2D>& velo) const
{
	const int& id = W.getParallel().myidWork;
	std::vector<Point2D>& Vel = velo;

	//Скорости виртуальных вихрей
	if (id == 0)
	{
#pragma omp parallel for default(none) shared(Vel)		
		for (int i = 0; i < (int)Vel.size(); ++i)
			Vel[i] = afl.getV(virtualWake.aflPan[i].second) \
			+ virtualWake.vecHalfGamma[i] \
			- W.getPassport().physicalProperties.V0();
	}
}
//CalcConvVelocityAtVirtualVortexes(...)


//Вычисление интенсивностей присоединенного вихревого слоя и присоединенного слоя источников
void BoundaryLinLayerAver::ComputeAttachedSheetsIntensity()
{
	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		oldSheets.attachedVortexSheet(i, 0) = sheets.attachedVortexSheet(i, 0);
		oldSheets.attachedVortexSheet(i, 1) = sheets.attachedVortexSheet(i, 1);
		oldSheets.attachedSourceSheet(i, 0) = sheets.attachedSourceSheet(i, 0);
		oldSheets.attachedSourceSheet(i, 1) = sheets.attachedSourceSheet(i, 1);
	}

	const Airfoil* oldAfl = (W.getCurrentStep() == 0) ? &afl : &W.getOldAirfoil(numberInPassport);

	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		sheets.attachedVortexSheet(i, 0) = 0.5 * (afl.getV(i + 1) + afl.getV(i)) & afl.tau[i];
		sheets.attachedVortexSheet(i, 1) = (afl.len[i] - oldAfl->len[i]) / W.getPassport().timeDiscretizationProperties.dt;			
			
		sheets.attachedSourceSheet(i, 0) = 0.5 * (afl.getV(i + 1) + afl.getV(i)) & afl.nrm[i];
		sheets.attachedSourceSheet(i, 1) = -afl.len[i] * (afl.phiAfl - oldAfl->phiAfl) / W.getPassport().timeDiscretizationProperties.dt;
                
                 
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

		velILin -= gamJ * (alpha * (u1 & taui) + lambda * (u1 & (-taui.kcross())));
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

		velILin -= gamJ * (alpha * (u1 &  taui.kcross()) + lambda * (u1 & taui) - 1.0);
	}

}//GetInfluenceFromSourcesToRectPanel(...)

//Вычисление влияния слоя источников конкретной прямолинейной панели на вихрь в области течения
void BoundaryLinLayerAver::GetInfluenceFromSourceSheetAtRectPanelToVortex(size_t panel, const Vortex2D& ptr, Point2D& vel) const
{
	vel.toZero();

	const Point2D& posI = ptr.r();

	Point2D dj = afl.getR(panel + 1) - afl.getR(panel);
	Point2D tauj = dj.unit();

	Point2D s = posI - afl.getR(panel);
	Point2D p = posI - afl.getR(panel + 1);

	Point2D u1 = 0.5 / dj.length()* VMlib::Omega(p + s, tauj, tauj);

	double a = VMlib::Alpha(p, s);

	double lambda;
	if ((s.length2() > 1e-16) && (p.length2() > 1e-16))
		lambda = VMlib::Lambda(p, s);
	else
		lambda = 0.0;

	vel += sheets.attachedSourceSheet(panel, 0) * (-a * u1.kcross() + lambda * u1 - tauj);
	vel *= IDPI;
}// GetInfluenceFromSourceSheetAtRectPanelToVortex(...)

//Вычисление влияния вихревых слоев (свободный + присоединенный) конкретной прямолинейной панели на вихрь в области течения
void BoundaryLinLayerAver::GetInfluenceFromVortexSheetAtRectPanelToVortex(size_t panel, const Vortex2D& ptr, Point2D& vel) const
{
	vel.toZero();

	const Point2D& posI = ptr.r();

	Point2D dj = afl.getR(panel + 1) - afl.getR(panel);
	Point2D tauj = dj.unit();

	Point2D s = posI - afl.getR(panel);
	Point2D p = posI - afl.getR(panel + 1);

	Point2D u1 = 0.5 / dj.length()* VMlib::Omega(p + s, tauj, tauj);

	double lambda;
	if ((s.length2() > 1e-16) && (p.length2() > 1e-16))
		lambda = VMlib::Lambda(p, s);
	else
		lambda = 0.0;

	vel += sheets.freeVortexSheet(panel, 0) * (lambda * u1.kcross() - tauj.kcross());
	vel += sheets.attachedVortexSheet(panel, 0) * (lambda * u1.kcross() - tauj.kcross());
	vel *= IDPI;
}// GetInfluenceFromVortexSheetAtRectPanelToVortex(...)

//Вычисляет влияния набегающего потока на прямолинейную панель для правой части
void BoundaryLinLayerAver::GetInfluenceFromVInfToRectPanel(std::vector<double>& vInfRhs) const
{
	size_t np = afl.getNumberOfPanels();
	//int id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(np);


	vInfRhs.resize(2 * np);

#pragma omp parallel for default(none) shared(vInfRhs, par, np)
	for (int i = 0; i < par.myLen; ++i)
	{
		vInfRhs[par.myDisp + i] = afl.tau[par.myDisp + i] & W.getPassport().physicalProperties.V0();
		vInfRhs[par.myDisp + np + i] =0;
	}
}// GetInfluenceFromVInfToRectPanel(...)


//Вычисляет влияния набегающего потока на криволинейную панель для правой части
void BoundaryLinLayerAver::GetInfluenceFromVInfToCurvPanel(std::vector<double>& vInfRhs) const
{
	size_t np = afl.getNumberOfPanels();
	//int id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(np);


	vInfRhs.resize(2*np);

#pragma omp parallel for default(none) shared(vInfRhs, par, np)
	for (int i = 0; i < par.myLen; ++i)
	{
		
		vInfRhs[par.myDisp + i] = - (afl.tau[par.myDisp + i] & W.getPassport().physicalProperties.V0()) + 1.0 / 24.0 * ((afl.nrm[par.myDisp + i] & (W.getPassport().physicalProperties.V0())) * ((AirfoilCurv)afl).getKc(par.myDisp + i) + \
			+ (afl.tau[par.myDisp + i] & W.getPassport().physicalProperties.V0()) * ((AirfoilCurv)afl).getDkc(par.myDisp + i) * ((AirfoilCurv)afl).getKc(par.myDisp + i)) * afl.len[par.myDisp + i] * afl.len[par.myDisp + i];
	
		vInfRhs[par.myDisp + np + i] = 1.0 / 12.0 * (afl.nrm[par.myDisp + i] & W.getPassport().physicalProperties.V0()) * ((AirfoilCurv)afl).getKc(par.myDisp + i) * afl.len[par.myDisp + i];

	}
}// GetInfluenceFromVorticesToCurvPanel(...)
