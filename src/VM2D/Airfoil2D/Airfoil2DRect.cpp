/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
\brief Файл кода с описанием класса AirfoilRect
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7   
\date 22 ноября 2019 г.
*/

#include "Airfoil2DRect.h"

#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

//Считывание профиля из файла
void AirfoilRect::ReadFromFile(const std::string& dir) //загрузка профиля из файла, его поворот и масштабирование
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


		v_.resize(0);
		for (size_t q = 0; q < r_.size(); ++q)
			v_.push_back({ 0.0, 0.0 });		

		Move(param.basePoint);
		Scale(param.scale);
		Rotate(param.angle);
		//в конце Rotate нормали, касательные и длины вычисляются сами

		gammaThrough.clear();
		gammaThrough.resize(r_.size(), 0.0);
	}
}//ReadFromFile(...)

//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек, обусловленных геометрией профиля, и вычисление вязкого трения
void AirfoilRect::GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3)
{
	std::vector<double> selfI0;
	std::vector<Point2D> selfI3;

	if (&pointsDb == &(W.getWake()))
	{
		viscousStress.clear();
		viscousStress.resize(r_.size(), 0.0);
	}
	
	std::vector<double> locViscousStress;
	locViscousStress.resize(r_.size(), 0.0);

	
	std::vector<double> currViscousStress;
	currViscousStress.resize(r_.size(), 0.0);


	const int& id = W.getParallel().myidWork;

	VMlib::parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
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
	double iDDomRad, expon;
#pragma warning (pop)

	double iDPIeps2 = 1.0 / (PI * sqr(2.0*W.getPassport().wakeDiscretizationProperties.epscol));

#pragma omp parallel for \
	default(none) \
	shared(locI0, locI3, domainRadius, locPoints, id, iDPIeps2, par, locViscousStress) \
	private(xi, xi_m, lxi, lxi_m, lenj_m, v0, q, new_n, mn, h, d, s, vec, vs, expon, iDDomRad) schedule(dynamic, DYN_SCHEDULE)
	for (int i = 0; i < par.myLen; ++i)
	{
		iDDomRad = 1.0 / domainRadius[i + par.myDisp];

		for (size_t j = 0; j < r_.size() - 1; j++)
		{
			vs = 0.0;
			q = locPoints[i].r() - 0.5 * (r_[j] + r_[j + 1]);
			vec = tau[j];

			s = q * vec;
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

					if (locI0[i] != -PI * domainRadius[i + par.myDisp])
					{
						locI0[i] += xi * mn * (lxi + 1.0) / (lxi*lxi);
						locI3[i] += mn;
					}

					vs = locPoints[i].g() * expon * iDPIeps2;
				}
				else if ((d <= 5.0 * len[j]) && (d >= 0.1 * len[j]))
				{
					vs = 0.0;
					new_n = static_cast<int>(ceil(5.0 * len[j] / d));
					//new_n = 100;
					h = v0 * (1.0 / new_n);

					for (int m = 0; m < new_n; m++)
					{
						xi_m = (locPoints[i].r() - (r_[j] + h * (m + 0.5))) * iDDomRad;
						lxi_m = xi_m.length();
					
						lenj_m = len[j] / new_n;
						expon = exp(-lxi_m) * lenj_m;

						mn = nrm[j] * expon;
						if (locI0[i] != -PI * domainRadius[i + par.myDisp])
						{
							locI0[i] += xi_m*  mn * (lxi_m + 1.0) / (lxi_m*lxi_m);
							locI3[i] += mn;
						}
						vs += expon;
					}//for m
					vs *= locPoints[i].g() * iDPIeps2;
				}
				else if (d <= 0.1 * len[j])
				{
					if (locI0[i] != -PI * domainRadius[i + par.myDisp])
					{
						locI0[i] = -PI * domainRadius[i + par.myDisp];
						
						if (fabs(s) > 0.5 * len[j])
						{
							locI3[i] = 2.0 * nrm[j] * domainRadius[i + par.myDisp] * (exp(-fabs(s)  * iDDomRad) * sinh(len[j] * iDDomRad / 2.0));
							//viscousStress[j] += 2.0 * locPoints[i].g()*(exp(-fabs(s)  * iDDomRad) * sinh(len[j] * iDDomRad / 2.0))  * iDDomRad / (PI * len[j]);
						}
						else
						{
							locI3[i] = 2.0 * nrm[j] * domainRadius[i + par.myDisp] * (1.0 - exp(-len[j] * iDDomRad / 2.0)*cosh(fabs(s) * iDDomRad));
							//viscousStress[j] += 2.0 * locPoints[i].g()* (1.0 - exp(-len[j] * iDDomRad / 2.0)*cosh(fabs(s)  * iDDomRad))  * iDDomRad / (PI * len[j]);
						}
					}

					vs = 0.0;
					new_n = static_cast<int>(ceil(5.0 * len[j] / d));
					//new_n = 100;
					h = v0 * (1.0 / new_n);

					for (int m = 0; m < new_n; m++)
					{
						xi_m = (locPoints[i].r() - (r_[j] + h * (m + 0.5))) * iDDomRad;
						lxi_m = xi_m.length();

						lenj_m = len[j] / new_n;
						expon = exp(-lxi_m) * lenj_m;													
						vs += expon;
					}//for m
					vs *= locPoints[i].g() * iDPIeps2;

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
		if (I0[i] != -PI * domainRadius[i])
		{
			if (selfI0[i] == -PI * domainRadius[i])
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

#if defined(USE_CUDA)
void AirfoilRect::GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3)
{
	const size_t npt = pointsDb.vtx.size();
	double*& dev_ptr_pt = pointsDb.devVtxPtr;
	double*& dev_ptr_rad = pointsDb.devRadPtr;
	const size_t nr = r_.size();
	double*& dev_ptr_r = devRPtr; 
	std::vector<double>& loci0 = pointsDb.tmpI0;
	std::vector<Point2D>& loci3 = pointsDb.tmpI3;
	double*& dev_ptr_i0 = pointsDb.devI0Ptr;
	double*& dev_ptr_i3 = pointsDb.devI3Ptr;

	double*& dev_ptr_visstr = devViscousStressesPtr;
	std::vector<double>& locvisstr = tmpViscousStresses;



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
		cuCalculateSurfDiffVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nr, dev_ptr_r, dev_ptr_i0, dev_ptr_i3, dev_ptr_rad, dev_ptr_visstr);
			
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

		MPI_Reduce(locvisstr.data(), newViscousStress.data(), nr, MPI_DOUBLE, MPI_SUM, 0, W.getParallel().commWork);


		if (id == 0)
			for (size_t i = 0; i < viscousStress.size(); ++i)
			{
				viscousStress[i] += newViscousStress[i];
			}

		if (id == 0)
			for (size_t q = 0; q < I3.size(); ++q)
			{
				if (I0[q] != -PI * domainRadius[q])
				{
					I0[q] = newI0[q];
					I3[q] = newI3[q];
				}
				else
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

bool AirfoilRect::IsPointInAirfoil(const Point2D& point) const
{
	double angle = 0.0;

	Point2D v1, v2;

	for (size_t i = 0; i < r_.size(); i++)
	{
		v1 = getR(i) - point;
		v2 = getR(i+1) - point;
		angle += atan2(v1^v2, v1 * v2);
	}

	if (fabs(angle) < 0.1)
		return false;
	else 
		return true;
}//IsPointInAirfoil(...)

//Перемещение профиля 
void AirfoilRect::Move(const Point2D& dr)	//перемещение профиля как единого целого на вектор dr
{
	for (size_t i = 0; i < r_.size(); i++)
		r_[i] += dr;
	rcm += dr;
	CalcNrmTauLen();
	GetGabarits();
}//Move(...)


//Поворот профиля 
void AirfoilRect::Rotate(double alpha)	//поворот профиля на угол alpha вокруг центра масс
{
	phiAfl += alpha;
	numvector<numvector<double, 2>, 2> rotMatrix = { { cos(alpha), sin(alpha) }, { -sin(alpha), cos(alpha) } };

	for (size_t i = 0; i < r_.size(); i++)
		r_[i] = rcm + dot(rotMatrix, r_[i] - rcm);

	CalcNrmTauLen();
	GetGabarits();
}//Rotate(...)


//Масштабирование профиля
void AirfoilRect::Scale(double factor)	//масштабирование профиля на коэффициент factor относительно центра масс
{
	for (size_t i = 0; i < r_.size(); i++)
		r_[i] = rcm + factor*(r_[i] - rcm);

	CalcNrmTauLen(); //строго говоря, меняются при этом только длины, нормали и касательные - без изменения
	GetGabarits();
}//Scale(...)

// Вычисление нормалей, касательных и длин панелей по текущему положению вершин
void AirfoilRect::CalcNrmTauLen()
{
	if (nrm.size() != r_.size())
	{
		nrm.resize(r_.size());
		tau.resize(r_.size());
		len.resize(r_.size());
	}

	Point2D rpan;

	for (size_t i = 0; i < r_.size(); ++i)
	{
		rpan = (getR(i + 1) - getR(i));
		len[i] = rpan.length();
		tau[i] = rpan.unit();

		nrm[i] = { tau[i][1], -tau[i][0] };
	}

	//закольцовываем
	nrm.push_back(nrm[0]);
	tau.push_back(tau[0]);
	len.push_back(len[0]);
}//CalcNrmTauLen()

//Вычисляет габаритный прямоугольник профиля
void AirfoilRect::GetGabarits(double gap)	//определение габаритного прямоугольника
{
	lowLeft = { 1E+10, 1E+10 };
	upRight = { -1E+10, -1E+10 };

	for (size_t i = 0; i < r_.size(); i++)
	{
		lowLeft[0] = std::min(lowLeft[0], r_[i][0]);
		lowLeft[1] = std::min(lowLeft[1], r_[i][1]);

		upRight[0] = std::max(upRight[0], r_[i][0]);
		upRight[1] = std::max(upRight[1], r_[i][1]);
	}

	Point2D size = upRight - lowLeft;
	lowLeft -= size*gap;
	upRight += size*gap;
}//GetGabarits(...)



//Вычисление коэффициентов матрицы A для расчета влияния профиля самого на себя
std::vector<double> AirfoilRect::getA(size_t p, size_t i, const Airfoil& otherAirfoil, size_t j) const
{
	std::vector<double> res(p*p, 0.0);

	bool self = (&otherAirfoil == this);

	if ( (i == j) && self )
	{
		res[0] = 0.5 * (tau[i] ^ nrm[i]);
		if (p == 1)
			return res;

		res[3] = (1.0 / 12.0) * res[0];
		if (p == 2)
			return res;

		if (p > 2)
			throw (-42);

	}//if i==j
		
	const Point2D& taui = tau[i];

	res[0] = W.getIQ(numberInPassport, otherAirfoil.numberInPassport).first(i, j);

	if (p == 1)
		return res;

	res[1] = W.getIQ(numberInPassport, otherAirfoil.numberInPassport).first(i, otherAirfoil.getNumberOfPanels() + j);

	res[2] = W.getIQ(numberInPassport, otherAirfoil.numberInPassport).first(getNumberOfPanels() + i, j);

	res[3] = W.getIQ(numberInPassport, otherAirfoil.numberInPassport).first(getNumberOfPanels() + i, otherAirfoil.getNumberOfPanels() + j);

	if (p == 2)
		return res;

	if (p > 2)
		throw (-42);
	
	//Хотя сюда никогда и не попадем
	return res;


}//getA(...)

//Вычисление коэффициентов матрицы A для расчета влияния профиля самого на себя
void AirfoilRect::calcIQ(size_t p, const Airfoil& otherAirfoil, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& matrPair) const
{
	bool self = (&otherAirfoil == this);

	size_t aflNSelf = numberInPassport;
	size_t aflNOther = otherAirfoil.numberInPassport;

	size_t npI;
	size_t npJ;

	numvector<double, 3> alpha, lambda;

	//auxillary vectors
	Point2D p1, s1, p2, s2, di, dj, i00, i01, i10, i11;
	numvector<Point2D, 3> v00, v11;
	numvector<Point2D, 2> v01, v10;

#pragma omp parallel for \
	default(none) \
	shared(otherAirfoil, self, aflNSelf, aflNOther, matrPair, p) \
	private(npI, npJ, alpha, lambda, p1, s1, p2, s2, di, dj, i00, i01, i10, i11, v00, v11, v01, v10) schedule(dynamic, DYN_SCHEDULE)
	for(int i = 0; i < getNumberOfPanels(); ++i)
	for (int j = 0; j < otherAirfoil.getNumberOfPanels(); ++j)
	{
		npI = getNumberOfPanels();
		npJ = otherAirfoil.getNumberOfPanels();

		if ((i == j) && self)
		{
			
			matrPair.first(i, j) = 0.0;
			matrPair.second(i, j) = 0.0;

			if (p == 2)
			{
				matrPair.first(i, npJ + j) = -0.5 * IDPI * len[i] * tau[i] & nrm[i];
				matrPair.first(npI + i, j) = -0.5 * IDPI * len[i] * tau[i] & nrm[i];

				matrPair.second(i, npJ + j) = -0.5 * IDPI * len[i] * tau[i] & tau[i];
				matrPair.second(npI + i, j) = -0.5 * IDPI * len[i] * tau[i] & tau[i];

				matrPair.first(npI + i, npJ + j) = 0.0;
				matrPair.second(npI + i, npJ + j) = 0.0;

			}
				if (p == 2)
					throw (-42);
		}//if i==j


		const Point2D& taui = tau[i];
		const Point2D& tauj = otherAirfoil.tau[j];

		p1 = getR(i + 1) - otherAirfoil.getR(j + 1);
		s1 = getR(i + 1) - otherAirfoil.getR(j);
		p2 = getR(i) - otherAirfoil.getR(j + 1);
		s2 = getR(i) - otherAirfoil.getR(j);
		di = getR(i + 1) - getR(i);
		dj = otherAirfoil.getR(j + 1) - otherAirfoil.getR(j);

		alpha = { \
			(self && isAfter(j, i)) ? 0.0 : VMlib::Alpha(s2, s1), \
			VMlib::Alpha(s2, p1), \
			(self && isAfter(i, j)) ? 0.0 : VMlib::Alpha(p1, p2) \
		};

		lambda = { \
			(self && isAfter(j, i)) ? 0.0 : VMlib::Lambda(s2, s1), \
			VMlib::Lambda(s2, p1), \
			(self && isAfter(i, j)) ? 0.0 : VMlib::Lambda(p1, p2) \
		};

		v00 = {
			VMlib::Omega(s1, taui, tauj),
			-VMlib::Omega(di, taui, tauj),
			VMlib::Omega(p2, taui, tauj)
		};

		i00 = IDPI / len[i] * (-(alpha[0] * v00[0] + alpha[1] * v00[1] + alpha[2] * v00[2]).kcross() \
			+ (lambda[0] * v00[0] + lambda[1] * v00[1] + lambda[2] * v00[2]));

		matrPair.first(i, j) = i00 & nrm[i];
		matrPair.second(i, j) = i00 & tau[i];


		if (p == 2)
		{

			v01 = {
				0.5 / (dj.length()) * ((p1 + s1) * tauj * VMlib::Omega(s1, taui, tauj) - s1.length2() * taui),
				-0.5 * di.length() / dj.length() * VMlib::Omega(s1 + p2, tauj, tauj)
			};

			i01 = IDPI / len[i] * (-((alpha[0] + alpha[2]) * v01[0] + (alpha[1] + alpha[2]) * v01[1]).kcross() \
				+ ((lambda[0] + lambda[2]) * v01[0] + (lambda[1] + lambda[2]) * v01[1]) - 0.5 * di.length() * tau[j]);

			matrPair.first(i, npJ + j) = i01 & nrm[i];
			matrPair.second(i, npJ + j) = i01 & tau[i];

			v10 = {
				0.5 / di.length() * ((s1 + s2) * taui * VMlib::Omega(s1, taui, tauj) - s1.length2() * tauj),
				0.5 * dj.length() / di.length() * VMlib::Omega(s1 + p2, taui, taui)
			};

			i10 = IDPI / len[i] * (-((alpha[0] + alpha[2]) * v10[0] + alpha[2] * v10[1]).kcross() \
				+ ((lambda[0] + lambda[2]) * v10[0] + lambda[2] * v10[1]) + 0.5 * dj.length() * tau[i]);

			matrPair.first(npI + i, j) = i10 & nrm[i];
			matrPair.second(npI + i, j) = i10 & tau[i];


			v11 = {
				1.0 / (12.0 * di.length() * dj.length()) * (2.0 * (s1 * VMlib::Omega(s1 - 3.0 * p2, taui, tauj)) * VMlib::Omega(s1, taui, tauj) - s1.length2() * (s1 - 3.0 * p2)) - 0.25 * VMlib::Omega(s1, taui, tauj),
				-di.length() / (12.0 * dj.length()) * VMlib::Omega(di, tauj, tauj),
				-dj.length() / (12.0 * di.length()) * VMlib::Omega(dj, taui, taui)
			};

			i11 = IDPI / len[i] * (-((alpha[0] + alpha[2]) * v11[0] + (alpha[1] + alpha[2]) * v11[1] + alpha[2] * v11[2]).kcross()\
				+ (lambda[0] + lambda[2]) * v11[0] + (lambda[1] + lambda[2]) * v11[1] + lambda[2] * v11[2] \
				+ 1.0 / 12.0 * (dj.length() * taui + di.length() * tauj - 2.0 * VMlib::Omega(s1, taui, tauj)));

			matrPair.first(npI + i, npJ + j) = i11 & nrm[i];
			matrPair.second(npI + i, npJ + j) = i11 & tau[i];
		}


		if (p > 2)
			throw (-42);
	}


}//getIQ(...)



//Вычисление влияния присоединенных слоев от другого профиля (константные базисные функции)
void AirfoilRect::getInfAttFromOther0(std::vector<double>& attOtherVelo, const Airfoil& otherAirfoil, size_t currentRow, size_t currentCol) const
{
	attOtherVelo.resize(getNumberOfPanels(), 0.0);
	for (size_t i = 0; i < getNumberOfPanels(); ++i)
	{
		for (size_t j = 0; j < otherAirfoil.getNumberOfPanels(); j++)
		{
			Point2D di = getR(i + 1) - getR(i);
			Point2D dj = otherAirfoil.getR(j + 1) - otherAirfoil.getR(j);
			Point2D s1 = getR(i + 1) - otherAirfoil.getR(j);
			Point2D s2 = getR(i) - otherAirfoil.getR(j);
			Point2D p1 = getR(i + 1) - otherAirfoil.getR(j + 1);
			Point2D p2 = getR(i) - otherAirfoil.getR(j + 1);

			double a1 = VMlib::Alpha(s2, s1);
			double a2 = VMlib::Alpha(s2, p1);
			double a3 = VMlib::Alpha(p1, p2);

			double lambda1 = VMlib::Lambda(s1, s2);
			double lambda2 = VMlib::Lambda(p1, s2);
			double lambda3 = VMlib::Lambda(p2, p1);

			Point2D v1 = VMlib::Omega(s1, tau[i], otherAirfoil.tau[j]);
			Point2D v2 = -VMlib::Omega(di, tau[i], otherAirfoil.tau[j]);
			Point2D v3 = VMlib::Omega(p2, tau[i], otherAirfoil.tau[j]);

			attOtherVelo[i] += -IDPI / len[i] * W.getBoundary(W.getNumberOfAirfoil()).sheets.attachedVortexSheet(j, 0) * tau[i] * (((-(a1 * v1 + a2 * v2 + a3 * v3).kcross()) + lambda1 * v1 + lambda2 * v2 + lambda3 * v3).kcross());
			attOtherVelo[i] += -IDPI / len[i] * W.getBoundary(W.getNumberOfAirfoil()).sheets.attachedSourceSheet(j, 0) * tau[i] * (((-(a1 * v1 + a2 * v2 + a3 * v3).kcross()) + lambda1 * v1 + lambda2 * v2 + lambda3 * v3));
		}//for j
	}//for i
}//getInfAttFromOther0(...)

//Вычисление влияния присоединенных слоев от другого профиля (константные и линейные базисные функции)
void AirfoilRect::getInfAttFromOther1(std::vector<double>& attOtherVelo, const Airfoil& otherAirfoil, size_t currentRow, size_t currentCol) const
{

}//getInfAttFromOther1(...)


void AirfoilRect::GetInfluenceFromVorticesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const
{
	W.getBoundary(numberInPassport).GetInfluenceFromVorticesToRectPanel(panel, ptr, count, panelRhs);
}//GetInfluenceFromVorticesToPanel(...)


//Вычисление влияния части подряд идущих источников из области течения на панель для правой части
void AirfoilRect::GetInfluenceFromSourcesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const
{
	W.getBoundary(numberInPassport).GetInfluenceFromSourcesToRectPanel(panel, ptr, count, panelRhs);
}//GetInfluenceFromSourcesToPanel(...)