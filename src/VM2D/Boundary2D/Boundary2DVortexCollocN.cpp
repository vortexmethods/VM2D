/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Boundary2DVortexCollocN.cpp                                      |
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
\brief Файл кода с описанием класса BoundaryVortexCollocN
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/



#include "Boundary2DVortexCollocN.h"

#include "Airfoil2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

//Конструктор
BoundaryVortexCollocN::BoundaryVortexCollocN(const World2D& W_, size_t numberInPassport_) :
	Boundary(W_, numberInPassport_, 1) 
	{
		c.reserve(afl.getNumberOfPanels());
		for (size_t i=0; i < afl.getNumberOfPanels(); ++i)
			c.push_back(0.5*(afl.getR(i)+afl.getR(i+1)));
	};


//Пересчет решения на интенсивность вихревого слоя
void BoundaryVortexCollocN::SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol)
{
	Vortex2D virtVort;
	Point2D midNorm;

	size_t np = afl.getNumberOfPanels();

	double delta = W.getPassport().wakeDiscretizationProperties.delta;
	
	int nVortPerPan = W.getPassport().wakeDiscretizationProperties.minVortexPerPanel;

	//Очистка и резервирование памяти
	virtualWake.vecHalfGamma.clear();
	virtualWake.vecHalfGamma.reserve(np * nVortPerPan);

	//Очистка и резервирование памяти
	virtualWake.aflPan.clear();
	virtualWake.aflPan.reserve(np * nVortPerPan);

	//Резервирование памяти
	virtualWake.vtx.clear();
	virtualWake.vtx.reserve(np * nVortPerPan);

	//Очистка и резервирование памяти
	vortexBeginEnd.clear();
	vortexBeginEnd.reserve(np);

	double maxG = W.getPassport().wakeDiscretizationProperties.maxGamma;

	std::pair<int, int> pair; 

	for (size_t i = 0; i < np; ++i)
	{
		midNorm = afl.nrm[i] * delta;
		size_t iNext = (i < np-1) ? i + 1 : 0;
		size_t iPrev = (i > 0) ? i - 1 : np - 1;

		double Gamma = (sol(i) / (afl.len[iPrev] + afl.len[i]) + sol(iNext) / (afl.len[iNext] + afl.len[i])) * afl.len[i];
		
		size_t NEWnVortPerPan = std::max(static_cast<int>(std::ceil(fabs(Gamma) / W.getPassport().wakeDiscretizationProperties.maxGamma)), nVortPerPan);
		
		pair.first = (int)virtualWake.vtx.size();


		Point2D dr = 1.0 / NEWnVortPerPan * (afl.getR(i + 1) - afl.getR(i));
		virtVort.g() = Gamma / NEWnVortPerPan;
		
		for (size_t j = 0; j < NEWnVortPerPan; ++j)
		{
			virtVort.r() = afl.getR(i) + (0.5 + j) * dr + delta * afl.nrm[i]; 
			
			virtualWake.vtx.push_back(virtVort);

			virtualWake.vecHalfGamma.push_back(0.5 * Gamma / afl.len[i] * afl.tau[i]); 
			virtualWake.aflPan.push_back({ numberInPassport, i });
		}

		pair.second = (int)virtualWake.vtx.size();
		vortexBeginEnd.push_back(pair);

		sheets.freeVortexSheet(i, 0) = Gamma / afl.len[i]; 		
	}
	

}//SolutionToFreeVortexSheetAndVirtualVortex(...)


//Генерация блока матрицы
void BoundaryVortexCollocN::FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol)
{
	size_t np = afl.getNumberOfPanels();
		
	for (size_t i = 0; i < np; ++i)
	{
		lactCol(i) = 1.0;
		lastLine(i) = 1.0;
	}

	
#pragma omp parallel for shared(matr)
	for (int i = 0; i < np; ++i)
		for (size_t j = 0; j < np; ++j)
		{						
			matr(i, j) = -(afl.tau[i] & (c[i]-afl.getR(j)) ) * IDPI / (c[i]-afl.getR(j)).length2();
		}
}//FillMatrixSelf(...)

void BoundaryVortexCollocN::FillIQSelf(std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ)
{
	//afl.calcIQ(1, afl, IQ);
}//FillIQSelf(...)

//Генерация блока матрицы влияния от другого профиля того же типа
void BoundaryVortexCollocN::FillMatrixFromOther(const Boundary& otherBoundary, Eigen::MatrixXd& matr)
{
#pragma omp parallel for shared(matr)
	for (int i = 0; i < afl.getNumberOfPanels(); ++i)
	for (size_t j = 0; j < otherBoundary.afl.getNumberOfPanels(); ++j)		
	{
		
		matr(i, j) = -(afl.tau[i] & (c[i]-otherBoundary.afl.getR(j)) ) * IDPI / (c[i]-otherBoundary.afl.getR(j)).length2();
	}		
}//FillMatrixFromOther(...)


void BoundaryVortexCollocN::FillIQFromOther(const Boundary& otherBoundary, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ)
{
	//afl.calcIQ(1, otherBoundary.afl, IQ);
}//FillIQFromOther(...)


//Вычисление скоростей в наборе точек, вызываемых наличием слоев вихрей и источников на профиле
void BoundaryVortexCollocN::CalcConvVelocityToSetOfPointsFromSheets(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const
{	
	std::vector<Point2D> selfVelo(pointsDb.vtx.size());	

	double cft = IDPI;	 

#pragma warning (push)
#pragma warning (disable: 4101)
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	//double dst2eps, dst2;
#pragma warning (pop)

#pragma omp parallel for default(none) shared(selfVelo, pointsDb, cft, std::cout) private(velI, tempVel)
	for (int i = 0; i < pointsDb.vtx.size(); ++i)
	{
		velI.toZero();

		const Point2D& posI = pointsDb.vtx[i].r();
		
		/// \todo Тут надо разобраться, как должно быть...
		/// \todo сделать  if(move || deform)
		for (size_t j = 0; j < sheets.getSheetSize(); ++j)
		{
			Point2D dj = afl.getR(j + 1) - afl.getR(j);
			Point2D tauj = dj.unit();

			Point2D s = posI - afl.getR(j);
			Point2D p = posI - afl.getR(j + 1);

			double a = VMlib::Alpha(p, s);

			double lambda;
			if ( (s.length2() > 1e-16) && (p.length2() > 1e-16) )
				lambda = VMlib::Lambda(p, s);
			else
				lambda = 0.0;

			Point2D skos = -a * tauj.kcross() + lambda * tauj;

			velI += sheets.freeVortexSheet(j, 0) * skos.kcross();
			velI += sheets.attachedVortexSheet(j, 0) * skos.kcross();
			velI += sheets.attachedSourceSheet(j, 0) * skos;
		}//for j
		
		velI *= cft;
		selfVelo[i] = velI;

	}//for i

	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
}//CalcConvVelocityToSetOfPointsFromSheets(...)


#if defined(USE_CUDA)
void BoundaryVortexCollocN::GPUCalcConvVelocityToSetOfPointsFromSheets(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const
{	
	if (afl.numberInPassport == 0)
	{
		const size_t npt = pointsDb.vtx.size();
		double*& dev_ptr_pt = pointsDb.devVtxPtr;

		size_t npnl = afl.getNumberOfPanels();
		for (size_t q = 1; q < W.getNumberOfAirfoil(); ++q)
			npnl += W.getAirfoil(q).getNumberOfPanels();

		double*& dev_ptr_r = afl.devRPtr;
		double*& dev_ptr_freeVortexSheet = afl.devFreeVortexSheetPtr;
		double*& dev_ptr_attachedVortexSheet = afl.devAttachedVortexSheetPtr;
		double*& dev_ptr_attachedSourceSheet = afl.devAttachedSourceSheetPtr;

		double*& dev_ptr_freeVortexSheetLin = afl.devFreeVortexSheetLinPtr;
		double*& dev_ptr_attachedVortexSheetLin = afl.devAttachedVortexSheetLinPtr;
		double*& dev_ptr_attachedSourceSheetLin = afl.devAttachedSourceSheetLinPtr;

		std::vector<Point2D>& Vel = velo;		
		std::vector<Point2D> newV(npt);
		double*& dev_ptr_vel = pointsDb.devVelPtr;
		double eps2 = W.getPassport().wakeDiscretizationProperties.eps2;

		//Явная синхронизация слоев не нужна, т.к. она выполняется в Gpu::RefreshAfls() 
		if (npt > 0)
		{
			//double tt1 = omp_get_wtime();
			cuCalculateConvVeloWakeFromVirtual(npt, dev_ptr_pt, npnl, dev_ptr_r, \
				dev_ptr_freeVortexSheet, dev_ptr_freeVortexSheetLin, \
				dev_ptr_attachedVortexSheet, dev_ptr_attachedVortexSheetLin, \
				dev_ptr_attachedSourceSheet, dev_ptr_attachedSourceSheetLin, \
				dev_ptr_vel, eps2);			
			//double tt2 = omp_get_wtime();
			//std::cout << "SHEET: " << tt2 - tt1 << std::endl;

			W.getCuda().CopyMemFromDev<double, 2>(npt, dev_ptr_vel, (double*)newV.data());


			for (size_t q = 0; q < Vel.size(); ++q)
				Vel[q] += newV[q];
		}
	}
}
//GPUCalcConvVelocityToSetOfPointsFromSheets(...)
#endif


//Вычисление интенсивностей присоединенного вихревого слоя и присоединенного слоя источников
void BoundaryVortexCollocN::ComputeAttachedSheetsIntensity()
{

	for (size_t i = 0; i < sheets.getSheetSize(); ++i)
	{
		oldSheets.attachedVortexSheet(i, 0) = sheets.attachedVortexSheet(i, 0);
		oldSheets.attachedSourceSheet(i, 0) = sheets.attachedSourceSheet(i, 0);
	}

	for (size_t i = 0; i < sheets.getSheetSize(); ++i)
	{
		sheets.attachedVortexSheet(i, 0) = 0.5 * (afl.getV(i) + afl.getV(i + 1)) & afl.tau[i];
		sheets.attachedSourceSheet(i, 0) = 0.5 * (afl.getV(i) + afl.getV(i + 1)) & afl.nrm[i];
	}	
}//ComputeAttachedSheetsIntensity()


//Вычисляет влияния части подряд идущих вихрей из вихревого следа на прямолинейную панель для правой части
void BoundaryVortexCollocN::GetInfluenceFromVorticesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const
{
	double& velI = wakeRhs[0];

	const Point2D& posI = c[panel];
	
	for (size_t it = 0; it != count; ++it)
	{
		const Vortex2D& vt = ptr[it];

		const Point2D& posJ = vt.r();
		const double& gamJ = vt.g();

		Point2D d = posI - posJ;
		//std::cout << "d = " << d << std::endl;

		velI -= (gamJ * afl.len[panel] / std::max(d.length2(), W.getPassport().wakeDiscretizationProperties.eps2)) * (afl.tau[panel] & d);
	}
}// GetInfluenceFromVorticesToRectPanel(...)



//Вычисляет влияния части подряд идущих источников в области течения на прямолинейную панель для правой части
void BoundaryVortexCollocN::GetInfluenceFromSourcesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const
{
	double& velI = wakeRhs[0];

	const Point2D& posI = 0.5 * (afl.getR(panel) + afl.getR(panel + 1));

	for (size_t it = 0; it != count; ++it)
	{
		const Vortex2D& vt = ptr[it];

		const Point2D& posJ = vt.r();
		const double& gamJ = vt.g();

		Point2D d = posI - posJ;

		velI += (gamJ * afl.len[panel] / std::max(d.length2(), W.getPassport().wakeDiscretizationProperties.eps2)) * (afl.nrm[panel] & d);
	}	
}// GetInfluenceFromSourcesToRectPanel(...)


//Вычисление влияния слоя источников конкретной прямолинейной панели на вихрь в области течения
void BoundaryVortexCollocN::GetInfluenceFromSourceSheetAtRectPanelToVortex(size_t panel, const Vortex2D& ptr, Point2D& vel) const
{		
	vel.toZero();

	const Point2D& posI = ptr.r();

	Point2D dj = afl.getR(panel + 1) - afl.getR(panel);
	Point2D tauj = dj.unit();

	Point2D s = posI - afl.getR(panel);
	Point2D p = posI - afl.getR(panel + 1);

	double a = VMlib::Alpha(p, s);

	double lambda;
	if ((s.length2() > 1e-16) && (p.length2() > 1e-16))
		lambda = VMlib::Lambda(p, s);
	else
		lambda = 0.0;

	vel += sheets.attachedSourceSheet(panel, 0) * (-a * tauj.kcross() + lambda * tauj);
}// GetInfluenceFromSourceSheetAtRectPanelToVortex(...)

//Вычисление влияния вихревых слоев (свободный + присоединенный) конкретной прямолинейной панели на вихрь в области течения
void BoundaryVortexCollocN::GetInfluenceFromVortexSheetAtRectPanelToVortex(size_t panel, const Vortex2D& ptr, Point2D& vel) const
{
	vel.toZero();

	const Point2D& posI = ptr.r();

	Point2D dj = afl.getR(panel + 1) - afl.getR(panel);
	Point2D tauj = dj.unit();

	Point2D s = posI - afl.getR(panel);
	Point2D p = posI - afl.getR(panel + 1);
	double a = VMlib::Alpha(p, s);

	double lambda;
	if ((s.length2() > 1e-16) && (p.length2() > 1e-16))
		lambda = VMlib::Lambda(p, s);
	else
		lambda = 0.0;

	Point2D skos = -a * tauj.kcross() + lambda * tauj;

	vel += sheets.freeVortexSheet(panel, 0) * skos.kcross();
	vel += sheets.attachedVortexSheet(panel, 0) * skos.kcross();

}// GetInfluenceFromVortexSheetAtRectPanelToVortex(...)


//Вычисляет влияния набегающего потока на прямолинейную панель для правой части
void BoundaryVortexCollocN::GetInfluenceFromVInfToRectPanel(std::vector<double>& vInfRhs) const
{
	size_t np = afl.getNumberOfPanels();

#pragma omp parallel for default(none) shared(vInfRhs, np)
	for (int i = 0; i < np; ++i)	
		vInfRhs[i] = afl.nrm[i] & W.getPassport().physicalProperties.V0();
	
}// GetInfluenceFromVInfToRectPanel(...)

