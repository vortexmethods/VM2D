/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Airfoil2DRect.cpp                                                |
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
\brief Файл кода с описанием класса AirfoilRect
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
*/

#include "Airfoil2DCurv.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "nummatrix.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Tree2D.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"
#include <algorithm>

using namespace VM2D;

//Считывание профиля из файла
void AirfoilCurv::ReadFromFile(const std::string& dir) //загрузка профиля из файла, его поворот и масштабирование
{
	const AirfoilParams& param = W.getPassport().airfoilParams[numberInPassport];
	std::string filename = dir + param.fileAirfoil;
	std::ifstream airfoilFile;

	if (fileExistTest(filename, W.getInfo()))
	{
		std::stringstream airfoilFile(VMlib::Preprocessor(filename).resultString);

		VMlib::StreamParser airfoilParser(W.getInfo(), "airfoil file parser", airfoilFile);

		m = 0.0; //TODO
		J = 0.0; //TODO
		rcm = { 0.0, 0.0 }; //TODO
		inverse = param.inverse;

		airfoilParser.get("r", r_);
		if (inverse)
			std::reverse(r_.begin(), r_.end());

		airfoilParser.get("rc", rc_);
		if (inverse)
			std::reverse(rc_.begin(), rc_.end());

		airfoilParser.get("L", len);
		if (inverse)
			std::reverse(len.begin(), len.end());

		airfoilParser.get("k", k_);
		if (inverse)
			std::reverse(k_.begin(), k_.end());


		airfoilParser.get("dk", dk_);
		if (inverse)
			std::reverse(dk_.begin(), dk_.end());

		airfoilParser.get("ddk", ddk_);
		if (inverse)
			std::reverse(ddk_.begin(), ddk_.end());

		airfoilParser.get("kc", kc_);
		if (inverse)
			std::reverse(kc_.begin(), kc_.end());

		airfoilParser.get("dkc", dkc_);
		if (inverse)
			std::reverse(dkc_.begin(), dkc_.end());


		airfoilParser.get("ddkc", ddkc_);
		if (inverse)
			std::reverse(ddkc_.begin(), ddkc_.end());

		airfoilParser.get("tauc", tau);
		if (inverse)
			std::reverse(tau.begin(), tau.end());

		v_.resize(0);
		for (size_t q = 0; q < r_.size(); ++q)
			v_.push_back({ 0.0, 0.0 });


		gammaThrough.clear();
		gammaThrough.resize(r_.size(), 0.0);

		Move(param.basePoint);
		Scale(param.scale);
		Rotate(-param.angle);
	}
}//ReadFromFile(...)


//Перемещение профиля 
void AirfoilCurv::Move(const Point2D& dr)	//перемещение профиля как единого целого на вектор dr
{
	for (size_t i = 0; i < r_.size(); ++i)
	{
		r_[i] += dr;
		rc_[i] += dr;
	}	

	rcm += dr;
	GetGabarits();
}//Move(...)


//Поворот профиля 
void AirfoilCurv::Rotate(double alpha)	//поворот профиля на угол alpha вокруг центра масс
{
	phiAfl += alpha;
	nummatrix<double, 2, 2> rotMatrix = { { cos(alpha), -sin(alpha) }, { sin(alpha), cos(alpha) } };

	for (size_t i = 0; i < r_.size(); ++i)
	{
		r_[i] = rcm + (rotMatrix & (r_[i] - rcm));
		rc_[i] = rcm + (rotMatrix & (rc_[i] - rcm));
		tau[i] = (rotMatrix & tau[i]);
	}	


	CalcNrm();

	GetGabarits();
}//Rotate(...)


//Масштабирование профиля
void AirfoilCurv::Scale(double factor)	//масштабирование профиля на коэффициент factor относительно центра масс
{
	for (size_t i = 0; i < r_.size(); ++i)
	{
		r_[i] = rcm + factor*(r_[i] - rcm);
		rc_[i] = rcm + factor*(rc_[i] - rcm);
		len[i] *= factor;
	}

	GetGabarits();
}//Scale(...)


// Вычисление нормалей, касательных и длин панелей по текущему положению вершин
void AirfoilCurv::CalcNrm()
{
	if (nrm.size() != r_.size())
	{
		nrm.resize(r_.size());
	}


	for (size_t i = 0; i < r_.size(); ++i)
	{
		nrm[i] = { tau[i][1], -tau[i][0] };
	}
}//CalcNrmTauLen()

//Вычисляет габаритный прямоугольник профиля
void AirfoilCurv::GetGabarits(double gap)	//определение габаритного прямоугольника
{
	lowLeft = { 1E+10, 1E+10 };
	upRight = { -1E+10, -1E+10 };

	for (size_t i = 0; i < r_.size(); ++i)
	{
		lowLeft[0] = std::min(lowLeft[0], r_[i][0]);
		lowLeft[1] = std::min(lowLeft[1], r_[i][1]);

		upRight[0] = std::max(upRight[0], r_[i][0]);
		upRight[1] = std::max(upRight[1], r_[i][1]);
	}

	for (size_t i = 0; i < r_.size(); ++i)
	{
		lowLeft[0] = std::min(lowLeft[0], rc_[i][0]);
		lowLeft[1] = std::min(lowLeft[1], rc_[i][1]);

		upRight[0] = std::max(upRight[0], rc_[i][0]);
		upRight[1] = std::max(upRight[1], rc_[i][1]);
	}

	Point2D size = upRight - lowLeft;
	lowLeft -= size*gap;
	upRight += size*gap;
}//GetGabarits(...)


//Переделать для криволинейных панелей!!! Пока скопирована из AirfoilRect
//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек, обусловленных геометрией профиля, и вычисление вязкого трения
void AirfoilCurv::GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3)
{
	std::vector<double> selfI0;
	std::vector<Point2D> selfI3;

	if (&pointsDb == &(W.getWake()))
	{
		viscousStress.clear();
		viscousStress.resize(r_.size(), 0.0);
	}

	////-------------------------------------------------------------------------------------------
	// Вычисление средних значений eps на панелях (с одной панели может рождаться несколько вихрей)
	std::vector<double> midEpsOverPanel;
	double midEps;
	const Boundary& bnd = W.getBoundary(numberInPassport);
	auto vtxBegin = bnd.vortexBeginEnd[0].first;
	auto virtVortParams = W.getVelocity().virtualVortexesParams[numberInPassport];

	for (size_t i = 0; i < r_.size() - 1; ++i)
	{
		midEps = 0.0;
		for (auto it = bnd.vortexBeginEnd[i].first; it != bnd.vortexBeginEnd[i].second; ++it)
		{
			auto numVtx = it - vtxBegin;
			midEps += virtVortParams.epsastWake[numVtx];
		}
		midEps /= (bnd.vortexBeginEnd[i].second - bnd.vortexBeginEnd[i].first);
		midEpsOverPanel.push_back(midEps);
	}
	//----------------------------------------------------------------------------------------------

	std::vector<double> locViscousStress;
	locViscousStress.resize(r_.size(), 0.0);


	std::vector<double> currViscousStress;
	currViscousStress.resize(r_.size(), 0.0);


	const int& id = W.getParallel().myidWork;

	VMlib::parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D/*, VM2D::MyAlloc<VMlib::Vortex2D>*/>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	std::vector<double> locI0(par.myLen, 0.0);
	std::vector<Point2D> locI3(par.myLen, { 0.0, 0.0 });


#pragma warning (push)
#pragma warning (disable: 4101)
	//Локальные переменные для цикла
	Point2D q, xi, xi_m, v0;
	double lxi, lxi_m, lenj_m;

	Point2D mn;
	int new_n;
	Point2D h;

	Point2D vec;
	double s, d;

	double vs;
	double iDDomRad, domRad, expon;
#pragma warning (pop)

#pragma omp parallel for \
	default(none) \
	shared(locI0, locI3, domainRadius, locPoints, id, par, locViscousStress, midEpsOverPanel, PI) \
	private(xi, xi_m, lxi, lxi_m, lenj_m, v0, q, new_n, mn, h, d, s, vec, vs, expon, domRad, iDDomRad) schedule(dynamic, DYN_SCHEDULE)
	for (int i = 0; i < par.myLen; ++i)
	{
		domRad = std::max(domainRadius[i + par.myDisp], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());
		iDDomRad = 1.0 / domRad;

		for (size_t j = 0; j < r_.size(); ++j)
		{
			vs = 0.0;
			q = locPoints[i].r() - 0.5 * (getR(j) + getR(j + 1));
			vec = tau[j];

			s = q & vec;
			d = q.length();

			if (d < 50.0 * len[j])	//Почему зависит от длины панели???
			{
				v0 = vec * len[j];

				if (d > 5.0 * len[j])
				{
					xi = q * iDDomRad;
					lxi = xi.length();

					expon = exp(-lxi) * len[j];
					mn = nrm[j] * expon;

					if (locI0[i] != -PI * domRad)
					{
						locI0[i] += (xi & mn) * (lxi + 1.0) / (lxi*lxi);
						locI3[i] += mn;
					}

					vs = locPoints[i].g() * expon / (PI * midEpsOverPanel[j] * midEpsOverPanel[j]);
				}
				else if ((d <= 5.0 * len[j]) && (d >= 0.1 * len[j]))
				{
					vs = 0.0;
					new_n = static_cast<int>(ceil(5.0 * len[j] / d));
					//new_n = 100;
					h = v0 * (1.0 / new_n);

					for (int m = 0; m < new_n; ++m)
					{
						xi_m = (locPoints[i].r() - (getR(j) + h * (m + 0.5))) * iDDomRad;
						lxi_m = xi_m.length();

						lenj_m = len[j] / new_n;
						expon = exp(-lxi_m) * lenj_m;

						mn = nrm[j] * expon;
						if (locI0[i] != -PI * domRad)
						{
							locI0[i] += (xi_m & mn) * (lxi_m + 1.0) / (lxi_m*lxi_m);
							locI3[i] += mn;
						}
						vs += expon;
					}//for m
					vs *= locPoints[i].g() / (PI * midEpsOverPanel[j] * midEpsOverPanel[j]);
				}
				else if (d <= 0.1 * len[j])
				{
					if (locI0[i] != -PI * domRad)
					{
						locI0[i] = -PI * domRad;

						if (fabs(s) > 0.5 * len[j])
						{
							locI3[i] = 2.0 * nrm[j] * domRad * (exp(-fabs(s)  * iDDomRad) * sinh(len[j] * iDDomRad / 2.0));
							//viscousStress[j] += 2.0 * locPoints[i].g()*(exp(-fabs(s)  * iDDomRad) * sinh(len[j] * iDDomRad / 2.0))  * iDDomRad / (PI * len[j]);
						}
						else
						{
							locI3[i] = 2.0 * nrm[j] * domRad * (1.0 - exp(-len[j] * iDDomRad / 2.0)*cosh(fabs(s) * iDDomRad));
							//viscousStress[j] += 2.0 * locPoints[i].g()* (1.0 - exp(-len[j] * iDDomRad / 2.0)*cosh(fabs(s)  * iDDomRad))  * iDDomRad / (PI * len[j]);
						}
					}

					vs = 0.0;
					new_n = static_cast<int>(ceil(5.0 * len[j] / d));
					//new_n = 100;
					h = v0 * (1.0 / new_n);

					for (int m = 0; m < new_n; ++m)
					{
						xi_m = (locPoints[i].r() - (getR(j) + h * (m + 0.5))) * iDDomRad;
						lxi_m = xi_m.length();

						lenj_m = len[j] / new_n;
						expon = exp(-lxi_m) * lenj_m;
						vs += expon;
					}//for m
					vs *= locPoints[i].g() / (PI * midEpsOverPanel[j] * midEpsOverPanel[j]);

					//break;

				}
			}//if d<50 len 

#pragma omp atomic
			locViscousStress[j] += vs;
		}//for j
	}//for i

	if (id == 0)
	{
		selfI0.resize(pointsDb.vtx.size(), 0.0);
		selfI3.resize(pointsDb.vtx.size(), { 0.0, 0.0 });
	}

	MPI_Gatherv(locI0.data(), par.myLen, MPI_DOUBLE, selfI0.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);
	MPI_Gatherv(locI3.data(), par.myLen, Point2D::mpiPoint2D, selfI3.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);


	MPI_Reduce(locViscousStress.data(), currViscousStress.data(), (int)getNumberOfPanels(), MPI_DOUBLE, MPI_SUM, 0, W.getParallel().commWork);

	if (id == 0)
		for (size_t i = 0; i < viscousStress.size(); ++i)
		{
			viscousStress[i] += currViscousStress[i];
		}


	if (id == 0)
		for (size_t i = 0; i < I0.size(); ++i)
		{
			domRad = std::max(domainRadius[i + par.myDisp], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());

			if (I0[i] != -PI * domRad)
			{
				if (selfI0[i] == -PI * domRad)
				{
					I0[i] = selfI0[i];
					I3[i] = selfI3[i];
				}
				else
				{
					I0[i] += selfI0[i];
					I3[i] += selfI3[i];
				}
			}
		}
}; //GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(...)

//Переделать для криволинейных панелей!!! Пока скопирована из AirfoilRect
#if defined(USE_CUDA)
void AirfoilCurv::GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3)
{
	const size_t npt = pointsDb.vtx.size();
	double*& dev_ptr_pt = pointsDb.devVtxPtr;
	double*& dev_ptr_rad = pointsDb.devRadPtr;
	double*& dev_ptr_meanEps = devMeanEpsOverPanelPtr;
	const size_t nr = r_.size();
	double*& dev_ptr_r = devRPtr;
	std::vector<double> loci0(npt);
	std::vector<Point2D> loci3(npt);
	double*& dev_ptr_i0 = pointsDb.devI0Ptr;
	double*& dev_ptr_i3 = pointsDb.devI3Ptr;
	double minRad = W.getPassport().wakeDiscretizationProperties.getMinEpsAst();

	double*& dev_ptr_visstr = devViscousStressesPtr;
	std::vector<double> locvisstr(nr);



	const int& id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(npt, true);

	double tCUDASTART = 0.0, tCUDAEND = 0.0;

	tCUDASTART = omp_get_wtime();


	if (&pointsDb == &(W.getWake()))
	{
		viscousStress.clear();
		viscousStress.resize(r_.size(), 0.0);
	}

	std::vector<double> zeroVec(r_.size(), 0.0);
	cuCopyFixedArray(devViscousStressesPtr, zeroVec.data(), zeroVec.size() * sizeof(double));


	W.getCuda().CopyMemFromDev<double, 1>(nr, dev_ptr_visstr, &locvisstr[0]);


	if ((npt > 0) && (nr > 0))
	{
		cuCalculateSurfDiffVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nr, dev_ptr_r, dev_ptr_i0, dev_ptr_i3, dev_ptr_rad, dev_ptr_meanEps, minRad, dev_ptr_visstr);

		W.getCuda().CopyMemFromDev<double, 2>(par.myLen, dev_ptr_i3, (double*)&loci3[0]);
		W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_i0, &loci0[0]);
		W.getCuda().CopyMemFromDev<double, 1>(nr, dev_ptr_visstr, &locvisstr[0]);

		std::vector<Point2D> newI3;
		std::vector<double> newI0;
		std::vector<double> newViscousStress;

		if (id == 0)
		{
			newI3.resize(I3.size());
			newI0.resize(I0.size());
			newViscousStress.resize(nr);
		}

		MPI_Gatherv(loci3.data(), par.myLen, Point2D::mpiPoint2D, newI3.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
		MPI_Gatherv(loci0.data(), par.myLen, MPI_DOUBLE, newI0.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

		MPI_Reduce(locvisstr.data(), newViscousStress.data(), (int)nr, MPI_DOUBLE, MPI_SUM, 0, W.getParallel().commWork);


		if (id == 0)
			for (size_t i = 0; i < viscousStress.size(); ++i)
			{
				viscousStress[i] += newViscousStress[i];
			}

		if (id == 0)
			for (size_t q = 0; q < I3.size(); ++q)
			{
				if (newI0[q] == -PI * std::max(domainRadius[q], minRad))
				{
					I0[q] = newI0[q];
					I3[q] = newI3[q];
				}
				else if (I0[q] != -PI * std::max(domainRadius[q], minRad))
				{
					I0[q] += newI0[q];
					I3[q] += newI3[q];
				}
			}
	}

	tCUDAEND = omp_get_wtime();

	//W.info('t') << "DIFF_SURF_GPU: " << (tCUDAEND - tCUDASTART) << std::endl;


}//GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses
#endif

//Переделать для криволинейных панелей!!! Пока скопирована из AirfoilRect
bool AirfoilCurv::IsPointInAirfoil(const Point2D& point) const
{
	double angle = 0.0;

	Point2D v1, v2;

	for (size_t i = 0; i < r_.size(); ++i)
	{
		v1 = getR(i) - point;
		v2 = getR(i + 1) - point;
		angle += atan2(v1^v2, v1 & v2);
	}

	if (fabs(angle) < 0.1)
		return false;
	else
		return true;
}//IsPointInAirfoil(...)


//Вычисление коэффициентов матрицы A для расчета влияния панели на панель
std::vector<double> AirfoilCurv::getA(size_t p, size_t i, const Airfoil& airfoil, size_t j) const
{
	std::vector<double> res(p*p, 0.0);

	bool self = (&airfoil == this);

	if ((i == j) && self)
	{
		res[0] = 0.5 * (tau[i] ^ nrm[i]) + W.getIQ(numberInPassport, airfoil.numberInPassport).first(i, i);
		if (p == 1)
			return res;

		res[3] = (1.0 / 24.0) * (tau[i] ^ nrm[i]) + W.getIQ(numberInPassport, airfoil.numberInPassport).first(getNumberOfPanels() + i, getNumberOfPanels() + i);
		if (p == 2)
			return res;

		if (p > 2)
			throw (-42);

	}//if i==j

	res[0] = W.getIQ(numberInPassport, airfoil.numberInPassport).first(i, j);

	if (p == 1)
		return res;

	res[1] = W.getIQ(numberInPassport, airfoil.numberInPassport).first(i, airfoil.getNumberOfPanels() + j);

	res[2] = W.getIQ(numberInPassport, airfoil.numberInPassport).first(getNumberOfPanels() + i, j);

	res[3] = W.getIQ(numberInPassport, airfoil.numberInPassport).first(getNumberOfPanels() + i, airfoil.getNumberOfPanels() + j);

	if (p == 2)
		return res;

	if (p > 2)
		throw (-42);

	//Хотя сюда никогда и не попадем
	return res;


}//getA(...)


//Вычисление коэффициентов матрицы A для расчета влияния профиля самого на себя
void AirfoilCurv::calcIQ(size_t p, const Airfoil& otherAirfoil, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& matrPair) const
{
	bool self = (&otherAirfoil == this);

	size_t aflNSelf = numberInPassport;
	size_t aflNOther = otherAirfoil.numberInPassport;

	//auxillary vectors
	double alpha, beta;
	double a00, a01, a02, a10, a11, a20;

	size_t npI = getNumberOfPanels();
	size_t npJ = otherAirfoil.getNumberOfPanels();
#pragma omp parallel for \
	default(none) \
	shared(otherAirfoil, self, aflNSelf, aflNOther, matrPair, p, npI, npJ, PI) \
	private(alpha, beta, a00, a01, a02, a10, a11, a20) schedule(dynamic, DYN_SCHEDULE)
	for (int i = 0; i < (int)npI; ++i)
		for (int j = 0; j < (int)npJ; ++j)
		{		
			if ((i == j) && self)
			{

				matrPair.first(i, j) =  (36.0 * len[i] * kc_[i] + len[i] * len[i] * ddkc_[i]) / (144.0 * PI);
				matrPair.second(i, j) = 0;

				if (p == 2)
				{
					matrPair.first(i, npJ + j) = (len[i] * len[i] * dkc_[i]) / (144.0 * PI);//a10
					matrPair.first(npI + i, j) = (len[i] * len[i] * dkc_[i]) / (72.0 * PI);//a01

					matrPair.second(i, npJ + j) = 0;
					matrPair.second(npI + i, j) = 0;

					matrPair.first(npI + i, npJ + j) =  (len[i] * len[i] * len[i] * ddkc_[i]) / (3456.0 * PI);
					matrPair.second(npI + i, npJ + j) = 0;
				}
				if (p > 2)
				{
					matrPair.first(i, 2 * npJ + j) = (len[i] * len[i] * len[i] * ddkc_[i]) / (2160.0 * PI);//20
					matrPair.first(2 * npI + i, j) = (len[i] * len[i] * len[i] * ddkc_[i]) / (720.0 * PI);//02

					matrPair.second(i, 2 * npJ + j) = 0;//20
					matrPair.second(2 * npI + i, j) = 0;//02
				}

			}//if i==j
			else
			{

			const Point2D& taui = tau[i];
			const Point2D& tauj = otherAirfoil.tau[j];

			Point2D dij = getRc(i) - getRc(j);
			double d = dij.length();

			if (isAfter(i, j))
			{
				if (self)
					a00 = len[j] / (288.0 * PI) * (72.0 * k_[i] - 12.0 * (len[j] - 2.0 * len[i]) * dk_[i] + \
						(6.0 * len[i] * len[i] - 3.0 * len[i] * len[j] + 2.0 * len[i] * len[j] + 2.0 * len[j] * len[j]) * ddk_[i]);
				else a00 = 0;
			}
			if (isAfter(j, i))
			{
				if (self)
					a00 = len[j] / (288.0 * PI) * (72.0 * k_[j] + 12.0 * (len[j] - 2.0 * len[i]) * dk_[j] + \
						(6.0 * len[i] * len[i] - 3.0 * len[i] * len[j] + 2.0 * len[i] * len[j] + 2.0 * len[j] * len[j]) * ddk_[j]);
				else a00 = 0;
			}
			else
			{
				alpha = VMlib::Alpha(dij, taui);
				beta = VMlib::Alpha(dij, tauj);

				a00 = len[j] / (48.0 * PI * d * d * d) * \
					(2.0 * (len[j] * len[j] * sin(alpha + 2.0 * beta) + 12.0 * d * d * sin(alpha) + len[i] * len[i] * sin(3.0 * alpha))\
						+ d * (len[j] * len[j] * kc_[j] * cos(alpha + beta) + len[i] * len[i] * (d * dkc_[i] * cos(alpha) - kc_[i] * (3.0 * cos(2.0 * alpha) + d * kc_[i] * sin(alpha)))));
			}

			matrPair.first(i, j) = a00;
			matrPair.second(i, j)=a00;


			if (p == 2)
			{

				if (isAfter(i, j))
				{
					if (self)
						a01 = (len[j] * len[j]) / (576.0 * PI) * (4.0 * dk_[i] - (len[j] - len[i]) * ddk_[i]);
					else a01 = 0;
				}
				else if (isAfter(j, i))
				{
					if (self)
						a01 = (len[j] * len[j]) / (576.0 * PI) * (4.0 * dk_[j] + (len[j] - len[i]) * ddk_[j]);
					else a01 = 0;
				}
			else
			{
				alpha = VMlib::Alpha(dij, tau[i]);
				beta = VMlib::Alpha(dij, tau[j]);

				a01= (len[j] * len[j]) / (24.0 * PI * d * d) * sin(alpha + beta);
			}

				matrPair.first(i, npJ + j) = a01;
				matrPair.second(i, npJ + j) = a01;

				if (isAfter(i, j))
				{
					if (self)
						a10 = (len[j] * len[j]) / (576.0 * PI) * (8.0 * dk_[i] - (len[j] - 3.0 * len[i]) * ddk_[i]);
					else a10 = 0;
				}
				else if (isAfter(j, i))
				{
					if (self)
						a10 = (len[j] * len[j]) / (576.0 * PI) * (8.0 * dk_[j] - (len[j] - 3.0 * len[i]) * ddk_[j]);
					else a10 = 0;
				}
			else
			{
				alpha = VMlib::Alpha(dij, tau[i]);
				beta = VMlib::Alpha(dij, tau[j]);

				a10= (len[i] * len[j]) / (24.0 * PI * d * d) * (cos(alpha) * (d + k_[i] - 2.0 * sin(alpha)));
			}

				matrPair.first(npI + i, j) =a10;
				matrPair.second(npI + i, j) =a10;


				if (isAfter(i, j))
				{
					if (self)
						a11 = (len[i] * len[j] * len[j] * ddkc_[i]) / (3456.0 * PI);
					else a11 = 0;
				}
				else if (isAfter(j, i))
				{
					if (self)
						a11 = (len[i] * len[j] * len[j] * ddkc_[j]) / (3456.0 * PI);
					else a11 = 0;
				}
			else
			{
				alpha = VMlib::Alpha(dij, tau[i]);
				beta = VMlib::Alpha(dij, tau[j]);

				a11= (len[i] * len[j] * len[j]) / (288.0 * PI * d * d * d) * (d * k_[i] * cos(alpha + beta) - 2.0 * sin(2.0 * alpha + beta));
			}


				matrPair.first(npI + i, npJ + j) =a11;
				matrPair.second(npI + i, npJ + j) =a11;
			}


			if (p > 2)
			{

				if (isAfter(i, j))
				{
					if (self)
						a02 = (len[j] * len[j] * len[j]) / (2160.0 * PI) * ddk_[i];
					else a02 = 0;
				}
				else if (isAfter(j, i))
				{
					if (self)
						a02 = (len[j] * len[j] * len[j]) / (2160.0 * PI) * ddk_[j];
					else a02 = 0;
				}
				else
				{
					alpha = VMlib::Alpha(dij, tau[i]);
					beta = VMlib::Alpha(dij, tau[j]);

					a02 = (len[j] * len[j] * len[j]) / (180.0 * PI * d * d * d) * (d * k_[j] * cos(alpha + beta) + 2.0 * sin(alpha + 2.0 * beta));
				}


				matrPair.first(i, 2 * npJ + j) = a02;
				matrPair.second(i, 2 * npJ + j) = a02;

				if (isAfter(i, j))
				{
					if (self)
						a20 = (len[i] * len[i] * len[j] * ddk_[i]) / (720.0 * PI);
					else a20 = 0;
				}
				else if (isAfter(j, i))
				{
					if (self)
						a20 = (len[i] * len[i] * len[j] * ddk_[j]) / (720.0 * PI);
					else a20 = 0;
				}
				else
				{
					alpha = VMlib::Alpha(dij, tau[i]);
					beta = VMlib::Alpha(dij, tau[j]);

					a20 = (len[i] * len[i] * len[j]) / (720.0 * PI * d * d * d) * (2.0 * sin(3.0 * alpha) - d * (k_[i] * (3.0 * cos(2.0 * alpha) + d * k_[i] * sin(alpha)) - d * dk_[i] * cos(alpha)));
				}

				matrPair.first(2 * npI + i, j) = a20;
				matrPair.second(2 * npI + i, j) = a20;
			}
			}//else (i == j)
		}//for j
}//getIQ(...)


void AirfoilCurv::GetInfluenceFromVorticesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const
{
	W.getBoundary(numberInPassport).GetInfluenceFromVorticesToCurvPanel(panel, ptr, count, panelRhs);
}//GetInfluenceFromVorticesToPanel(...)

void AirfoilCurv::GetInfluenceFromVInfToPanel(std::vector<double>& vInfRhs) const
{
	W.getBoundary(numberInPassport).GetInfluenceFromVInfToCurvPanel(vInfRhs);
}//GetInfluenceFromVInfToPanel(...)