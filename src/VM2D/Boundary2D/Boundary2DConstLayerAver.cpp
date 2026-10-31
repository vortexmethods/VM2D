/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Boundary2DConstLayerAver.cpp                                     |
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
\brief Файл кода с описанием класса BoundaryConstLayerAver
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.14
\date 6 марта 2026 г.
*/



#include "Boundary2DConstLayerAver.h"

#include "Airfoil2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"
#include "Gmres2D.h"

#include <complex>

using namespace VM2D;


//Пересчет решения на интенсивность вихревого слоя
void BoundaryConstLayerAver::SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol)
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

		size_t NEWnVortPerPan = (size_t)std::max((int)std::ceil(fabs(sol(i) * afl.len[i]) / maxG), nVortPerPan);
		
		pair.first = (int)virtualWake.vtx.size();


		Point2D dr = 1.0 / NEWnVortPerPan * (afl.getR(i + 1) - afl.getR(i));

		for (size_t j = 0; j < NEWnVortPerPan; ++j)
		{
			virtVort.r() = afl.getR(i) + dr * (j * 1.0 + 0.5) + midNorm;
			virtVort.g() = sol(i) * afl.len[i] / NEWnVortPerPan;
			virtualWake.vtx.push_back(virtVort);

			virtualWake.vecHalfGamma.push_back(0.5 * sol(i)  * afl.tau[i]);
			virtualWake.aflPan.push_back({ numberInPassport, i });
		}

		pair.second = (int)virtualWake.vtx.size();
		vortexBeginEnd.push_back(pair);
	}
	

	for (size_t j = 0; j < np; ++j)
		sheets.freeVortexSheet(j, 0) = sol(j);
		
}//SolutionToFreeVortexSheetAndVirtualVortex(...)


//Генерация блока матрицы
void BoundaryConstLayerAver::FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol)
{
	size_t np = afl.getNumberOfPanels();
		
	for (size_t i = 0; i < np; ++i)
	{
		lactCol(i) = 1.0;
		lastLine(i) = afl.len[i];
	}

	for (size_t i = 0; i < np; ++i)
	for (size_t j = 0; j < np; ++j)
		matr(i, j) = afl.getA(1, i, afl, j)[0];

}//FillMatrixSelf(...)

void BoundaryConstLayerAver::FillIQSelf(std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ)
{
	afl.calcIQ(1, afl, IQ);
}//FillIQSelf(...)

//Генерация блока матрицы влияния от другого профиля того же типа
void BoundaryConstLayerAver::FillMatrixFromOther(const Boundary& otherBoundary, Eigen::MatrixXd& matr)
{
	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	for (size_t j = 0; j < otherBoundary.afl.getNumberOfPanels(); ++j)
		matr(i, j) = afl.getA(1, i, otherBoundary.afl, j)[0];
}//FillMatrixFromOther(...)


void BoundaryConstLayerAver::FillIQFromOther(const Boundary& otherBoundary, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ)
{
	afl.calcIQ(1, otherBoundary.afl, IQ);
}//FillIQFromOther(...)



/*
//Вычисление скоростей в наборе точек, вызываемых наличием слоев вихрей и источников на профиле
void BoundaryConstLayerAver::CalcConvVelocityToSetOfPointsFromSheets(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const
{	
	std::vector<Point2D> selfVelo(pointsDb.vtx.size());
	
	double cft = IDPI;

#pragma warning (push)
#pragma warning (disable: 4101)
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
#pragma warning (pop)

#pragma omp parallel for default(none) shared(selfVelo, cft, pointsDb, std::cout) private(velI, tempVel)
	for (int i = 0; i < pointsDb.vtx.size(); ++i)
	{
		/// \todo сделать вызов функции GetInfluenceFromVortexSheetAtRectPanelToVortex

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
*/

//KFU
//Вклад одной панели в скорость (без множителя 1/(2 pi))
static inline Point2D ConstLayerPanelVelo(const Point2D& posI, const Point2D& rj, const Point2D& rj1, double gFree, double gAtt, double qAtt)
{
	Point2D dj = rj1 - rj;
	Point2D tauj = dj.unit();

	Point2D s = posI - rj;
	Point2D p = posI - rj1;

	double a = VMlib::Alpha(p, s);

	double lambda;
	if ((s.length2() > 1e-16) && (p.length2() > 1e-16))
		lambda = VMlib::Lambda(p, s);
	else
		lambda = 0.0;

	Point2D skos = -a * tauj.kcross() + lambda * tauj;

	Point2D velI = gFree * skos.kcross();
	velI += gAtt * skos.kcross();
	velI += qAtt * skos;
	return velI;
}

//KFU
namespace
{
	/// Порядок мультипольных разложений для слоев на профиле
	const int sheetMpOrder = 20;

	/// Разложение используется для точек, удаленных от центра более чем на sheetMpRatio радиусов группы панелей;
	/// погрешность порядка sheetMpRatio^-(sheetMpOrder+1) ~ 5e-9
	const double sheetMpRatio = 2.5;

	/// Число панелей в группе
	const size_t sheetMpGroupSize = 8;

	/// Мультипольное разложение слоев на группе панелей [jBegin, jEnd):
	/// u - i v = 1/(2 pi) * sum_m mom[m] / (z - c)^(m+1)
	struct SheetMultipole
	{
		std::complex<double> c;
		double r2Far; //квадрат радиуса, вне которого применяется разложение
		size_t jBegin, jEnd;
		std::complex<double> mom[sheetMpOrder + 1];

		std::complex<double> Eval(const std::complex<double>& z) const
		{
			const std::complex<double> t = 1.0 / (z - c);
			std::complex<double> w = mom[sheetMpOrder];
			for (int m = sheetMpOrder - 1; m >= 0; --m)
				w = w * t + mom[m];
			return w * t;
		}
	};

	/// Вычисление моментов: для панели от A до B (относительно центра) с плотностью (q - i gamma)
	/// интеграл (zeta - c)^m ds = L * S_m / (m+1), где S_m = sum_{k=0..m} A^(m-k) B^k
	template<typename GetR, typename GetDens>
	SheetMultipole MakeSheetMultipole(size_t jBegin, size_t jEnd, GetR getR, GetDens getDens)
	{
		SheetMultipole mp;
		mp.jBegin = jBegin;
		mp.jEnd = jEnd;

		Point2D cnt = { 0.0, 0.0 };
		for (size_t j = jBegin; j <= jEnd; ++j)
			cnt += getR(j);
		cnt *= 1.0 / (jEnd - jBegin + 1);
		mp.c = { cnt[0], cnt[1] };

		double r2 = 0.0;
		for (size_t j = jBegin; j <= jEnd; ++j)
			r2 = std::max(r2, (getR(j) - cnt).length2());
		mp.r2Far = sheetMpRatio * sheetMpRatio * r2;

		for (int m = 0; m <= sheetMpOrder; ++m)
			mp.mom[m] = 0.0;

		for (size_t j = jBegin; j < jEnd; ++j)
		{
			const Point2D& ra = getR(j);
			const Point2D& rb = getR(j + 1);
			const std::complex<double> A(ra[0] - cnt[0], ra[1] - cnt[1]);
			const std::complex<double> B(rb[0] - cnt[0], rb[1] - cnt[1]);
			const std::complex<double> dens = getDens(j) * (rb - ra).length();

			std::complex<double> S = 1.0, Bm = 1.0;
			for (int m = 0; m <= sheetMpOrder; ++m)
			{
				mp.mom[m] += dens * S * (1.0 / (m + 1));
				Bm *= B;
				S = A * S + Bm;
			}
		}
		return mp;
	}
}

//KFU
//Вычисление скоростей в наборе точек, вызываемых наличием слоев вихрей и источников на профиле
void BoundaryConstLayerAver::CalcConvVelocityToSetOfPointsFromSheets(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const
{	
	std::vector<Point2D> selfVelo(pointsDb.vtx.size());
	
	double cft = IDPI;

	const size_t nPnl = sheets.getSheetSize();

	//В быстром режиме влияние удаленных групп панелей учитывается через мультипольные разложения
	const bool useMultipoles = (W.getPassport().numericalSchemes.velocityComputation.second == 1) && (nPnl > 0);

	auto getR = [this](size_t j) -> const Point2D& { return afl.getR(j); };
	auto getDens = [this](size_t j) -> std::complex<double> 
	{ 
		return std::complex<double>(sheets.attachedSourceSheet(j, 0), -(sheets.freeVortexSheet(j, 0) + sheets.attachedVortexSheet(j, 0)));
	};

	SheetMultipole mpAll;
	std::vector<SheetMultipole> mpGroups;
	if (useMultipoles)
	{
		mpAll = MakeSheetMultipole(0, nPnl, getR, getDens);
		for (size_t jb = 0; jb < nPnl; jb += sheetMpGroupSize)
			mpGroups.push_back(MakeSheetMultipole(jb, std::min(jb + sheetMpGroupSize, nPnl), getR, getDens));
	}

#pragma omp parallel for default(none) shared(selfVelo, cft, pointsDb, mpAll, mpGroups, useMultipoles, nPnl) schedule(dynamic, 256)
	for (int i = 0; i < (int)pointsDb.vtx.size(); ++i)
	{
		/// \todo сделать вызов функции GetInfluenceFromVortexSheetAtRectPanelToVortex
		Point2D velI = { 0.0, 0.0 };

		const Point2D& posI = pointsDb.vtx[i].r();
		const std::complex<double> z(posI[0], posI[1]);

		/// \todo Тут надо разобраться, как должно быть...
		/// \todo сделать  if(move || deform)
		if (!useMultipoles)
		{
			for (size_t j = 0; j < nPnl; ++j)
				velI += ConstLayerPanelVelo(posI, afl.getR(j), afl.getR(j + 1), sheets.freeVortexSheet(j, 0), sheets.attachedVortexSheet(j, 0), sheets.attachedSourceSheet(j, 0));
		}
		else if (std::norm(z - mpAll.c) > mpAll.r2Far)
		{
			const std::complex<double> w = mpAll.Eval(z);
			velI = { w.real(), -w.imag() };
		}
		else
		{
			std::complex<double> w = 0.0;
			for (const SheetMultipole& mp : mpGroups)
			{
				if (std::norm(z - mp.c) > mp.r2Far)
					w += mp.Eval(z);
				else
					for (size_t j = mp.jBegin; j < mp.jEnd; ++j)
						velI += ConstLayerPanelVelo(posI, afl.getR(j), afl.getR(j + 1), sheets.freeVortexSheet(j, 0), sheets.attachedVortexSheet(j, 0), sheets.attachedSourceSheet(j, 0));
			}
			velI += Point2D({ w.real(), -w.imag() });
		}
		
		velI *= cft;
		selfVelo[i] = velI;
	}//for i
		
#ifdef VM_CHECK_FAST
	if (useMultipoles && (W.getCurrentStep() % 25 == 0))
	{
		double err = 0.0, mx = 0.0;
		for (int i = 0; i < (int)pointsDb.vtx.size(); ++i)
		{
			Point2D v = { 0.0, 0.0 };
			for (size_t j = 0; j < nPnl; ++j)
				v += ConstLayerPanelVelo(pointsDb.vtx[i].r(), afl.getR(j), afl.getR(j + 1), sheets.freeVortexSheet(j, 0), sheets.attachedVortexSheet(j, 0), sheets.attachedSourceSheet(j, 0));
			v *= cft;
			err = std::max(err, (v - selfVelo[i]).length());
			mx = std::max(mx, v.length());
		}
		std::cout << "CHECK_FAST step " << W.getCurrentStep() << " sheets->points (n=" << pointsDb.vtx.size() << "): maxErr/maxVel = " << err / mx << std::endl;
	}
#endif

	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
}//CalcConvVelocityToSetOfPointsFromSheets(...)


#if defined(USE_CUDA)
void BoundaryConstLayerAver::GPUCalcConvVelocityToSetOfPointsFromSheets(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const
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
		double*& dev_ptr_vel = pointsDb.devVelPtr;
		

		if (npt > 0)
		{
			cuCalculateConvVeloWakeFromVirtual(npt, dev_ptr_pt, npnl, dev_ptr_r, \
				dev_ptr_freeVortexSheet, dev_ptr_freeVortexSheetLin, \
				dev_ptr_attachedVortexSheet, dev_ptr_attachedVortexSheetLin, \
				dev_ptr_attachedSourceSheet, dev_ptr_attachedSourceSheetLin, \
				dev_ptr_vel);
	
			std::vector<Point2D> newV(Vel.size());
			W.getCuda().CopyMemFromDev<double, 2>(npt, dev_ptr_vel, (double*)newV.data());

			for (size_t q = 0; q < Vel.size(); ++q)
				Vel[q] += newV[q];

		/*	VMlib::CreateUserDirectory(W.getPassport().dir, "dbg");
			std::ostringstream sss;
			sss << "velWakeFromBou";
			sss << W.currentStep;
			std::ofstream prmtFile(W.getPassport().dir + "dbg/" + sss.str());
			prmtFile << "i x y g B2Wx B2Wy" << std::endl;
			for (size_t i = 0; i < npt; ++i)
				prmtFile << i << " " \
				<< pointsDb.vtx[i].r()[0] << " " << pointsDb.vtx[i].r()[1] << " " \
				<< pointsDb.vtx[i].g() << " " \
				<< newV[i][0] << " " << newV[i][1] << " " \
				<< std::endl;
			prmtFile.close();*/
		}
	}
}
//GPUCalcConvVelocityToSetOfPointsFromSheets(...)
#endif



//Вычисление интенсивностей присоединенного вихревого слоя и присоединенного слоя источников
void BoundaryConstLayerAver::ComputeAttachedSheetsIntensity()
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
void BoundaryConstLayerAver::GetInfluenceFromVorticesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const
{
	double& velI = wakeRhs[0];

	const Point2D& posI0 = afl.getR(panel);
	const Point2D& posI1 = afl.getR(panel + 1);


	for (size_t it = 0; it != count; ++it)
	{
		const Vortex2D& vt = ptr[it];

		const Point2D& posJ = vt.r();
		const double& gamJ = vt.g();

		Point2D s = posJ - posI0;
		Point2D p = posJ - posI1;

		double alpha = VMlib::Alpha(p, s);

		velI -= gamJ * alpha;
	}
}// GetInfluenceFromVorticesToRectPanel(...)



//Вычисляет влияния части подряд идущих источников в области течения на прямолинейную панель для правой части
void BoundaryConstLayerAver::GetInfluenceFromSourcesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const
{
	double& velI = wakeRhs[0];

	const Point2D& posI0 = afl.getR(panel);
	const Point2D& posI1 = afl.getR(panel + 1);

	for (size_t it = 0; it != count; ++it)
	{
		const Vortex2D& vt = ptr[it];
		const Point2D& posJ = vt.r();
		const double& gamJ = vt.g();

		Point2D s = posJ - posI0;
		Point2D p = posJ - posI1;

		double lambda = VMlib::Lambda(p, s);

		velI -= gamJ * lambda;
	}
}// GetInfluenceFromSourcesToRectPanel(...)


//Вычисление влияния слоя источников конкретной прямолинейной панели на вихрь в области течения
void BoundaryConstLayerAver::GetInfluenceFromSourceSheetAtRectPanelToVortex(size_t panel, const Vortex2D& ptr, Point2D& vel) const
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
void BoundaryConstLayerAver::GetInfluenceFromVortexSheetAtRectPanelToVortex(size_t panel, const Vortex2D& ptr, Point2D& vel) const
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
void BoundaryConstLayerAver::GetInfluenceFromVInfToRectPanel(std::vector<double>& vInfRhs) const
{
	size_t np = afl.getNumberOfPanels();	
	vInfRhs.resize(np);

#pragma omp parallel for default(none) shared(vInfRhs, np)
	for (int i = 0; i < np; ++i)
		vInfRhs[i] = afl.tau[i] & W.getV0();	

}// GetInfluenceFromVInfToRectPanel(...)