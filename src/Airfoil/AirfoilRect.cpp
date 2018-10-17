/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: AirfoilRect.cpp                                                  |
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
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#include "AirfoilRect.h"
#include "Preprocessor.h"
#include "World2D.h"


//Считывание профиля из файла
void AirfoilRect::ReadFromFile(const std::string& dir) //загрузка профиля из файла, его поворот и масштабирование
{
	const AirfoilParams& param = W.getPassport().airfoilParams[numberInPassport];
	std::string filename = dir + param.fileAirfoil;
	std::ifstream airfoilFile;
	
	if (fileExistTest(filename, W.getInfo()))
	{
		std::stringstream airfoilFile(Preprocessor(filename).resultString);
		
		StreamParser airfoilParser(W.getInfo(), "airfoil file parser", airfoilFile);

		m = 0.0; //TODO
		J = 0.0; //TODO
		rcm = { 0.0, 0.0 }; //TODO
		
		airfoilParser.get("np", np);

		airfoilParser.get("r", r);

		//замыкаем --- в конец приписываем первую точку (для простоты)
		r.push_back(r[0]);

		v.resize(0);
		for (size_t q = 0; q < r.size(); ++q)
			v.push_back({ 0.0, 0.0 });

		Move(param.basePoint);
		Scale(param.scale);
		Rotate(param.angle);
		//в конце Rotate нормали, касательные и длины вычисляются сами

		gammaThrough.clear();
		gammaThrough.resize(np, 0.0);
	}
}//ReadFromFile(...)

//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек, обусловленных геометрией профиля, и вычисление вязкого трения
void AirfoilRect::GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3)
{
	std::vector<double> selfI0;
	std::vector<Point2D> selfI3;

	viscousStress.clear();
	viscousStress.resize(np, 0.0);

	const int& id = W.getParallel().myidWork;

	parProp par = W.getParallel().SplitMPI(pointsDb.vtx.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(par.myLen);

	MPI_Scatterv(const_cast<std::vector<Vortex2D>&>(pointsDb.vtx).data(), par.len.data(), par.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), par.myLen, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);

	std::vector<double> locI0(par.myLen, 0.0);
	std::vector<Point2D> locI3(par.myLen, {0.0, 0.0});
	

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

	double iDPIepscol2 = 1.0 / (PI * sqr(W.getPassport().wakeDiscretizationProperties.epscol));

#pragma omp parallel for \
	default(none) \
	shared(locI0, locI3, domainRadius, locPoints, id, iDPIepscol2, par) \
	private(xi, xi_m, lxi, lxi_m, lenj_m, v0, q, new_n, mn, h, d, s, vec, vs, expon, iDDomRad) schedule(dynamic, DYN_SCHEDULE)
	for (int i = 0; i < par.myLen; ++i)
	{
		iDDomRad = 1.0 / domainRadius[i+par.myDisp];

		for (size_t j = 0; j < r.size() - 1; j++)
		{
			q = locPoints[i].r() - 0.5 * (r[j] + r[j + 1]);
			vec = tau[j];

			s = q * vec;

			//POLARA
			//Сделала расстояние до центра, а не до конца панели
			/*if (fabs(s) > 0.5 * len[j])
				d = sqrt(sqr(fabs(q*vec) - 0.5*len[j]) + sqr(q*nrm[j]));
			else
				d = fabs(q*nrm[j]);*/
			d = q.length();

			if (d < 50.0 * len[j])	//Почему зависит от длины панели???
			{
				v0 = vec * len[j];

				if (d > 5.0 * len[j])
				{
					xi = q * iDDomRad;
					lxi = xi.length();
				
					expon = exp(-lxi);
					mn = len[j] * nrm[j] * expon;
					
					if (locI0[i] != -PI * domainRadius[i + par.myDisp])
					{
						locI0[i] += xi * mn * (lxi + 1.0) / (lxi*lxi);
						locI3[i] += mn;
					}

					viscousStress[j] += locPoints[i].g() * expon * iDPIepscol2;
				}
				else if ((d <= 5.0 * len[j]) && (d >= 0.1 * len[j]))
				{
					vs = 0.0;
					new_n = static_cast<int>(ceil(5.0 * len[j] / d));
					//new_n = 100;
					h = v0 * (1.0 / new_n);

					for (int m = 0; m < new_n; m++)
					{
						xi_m = (locPoints[i].r() - (r[j] + h * (m + 0.5))) * iDDomRad;
						lxi_m = xi_m.length();

						expon = exp(-lxi_m);
						lenj_m = len[j] / new_n;
						mn = lenj_m *nrm[j] * expon;
						if (locI0[i] != -PI * domainRadius[i + par.myDisp])
						{
							locI0[i] += xi_m*  mn * (lxi_m + 1.0) / (lxi_m*lxi_m);
							locI3[i] += mn;
						}
						vs += locPoints[i].g() * expon;
					}//for m
					viscousStress[j] += vs / new_n * iDPIepscol2;
				}
				else if (d <= 0.1 * len[j])
				{
					locI0[i] = -PI * domainRadius[i + par.myDisp];
					//locI3[i] = 2.0 * nrm[j] * domainRadius[i + par.myDisp] * (1.0 - exp(-0.5 * len[j] / domainRadius[i + par.myDisp]));
					//break;
					
					if (fabs(s) > 0.5 * len[j])
					{						
						locI3[i] = 2.0 * nrm[j] * domainRadius[i + par.myDisp] * (exp(-fabs(s)  * iDDomRad) * sinh(len[j] * iDDomRad / 2.0));
						//	viscousStress[j] += locPoints[i].g() * (mn.kcross() * tau[j]) / (PI * sqr(epscol));
						//	viscousStress[j] += locPoints[i].g() * (2.0 * nrm[j] * domainRadius[i+parallel.myDisp] * (1.0 - exp(-0.5 * len[j] / domainRadius[i+parallel.myDisp]))).kcross() * tau[j] / (PI * sqr(epscol) *len[j]);
						//viscousStress[j] += locPoints[i].g()*(exp(-fabs(s)  * iDDomRad) * sinh(len[j] * iDDomRad / 2.0))  * iDDomRad / len[j]; //locPoints[i].g() *(2.0 *  epscol* (1.0 - exp(-0.5 * len[j] / epscol))) / (PI * sqr(epscol) *len[j]);
						viscousStress[j] += 2.0 * locPoints[i].g()*(exp(-fabs(s)  * iDDomRad) * sinh(len[j] * iDDomRad / 2.0))  * iDDomRad / (PI * len[j]); 
					}
					else
					{
						locI3[i] = 2.0 * nrm[j] * domainRadius[i + par.myDisp] * (1.0 - exp(-len[j] * iDDomRad / 2.0)*cosh(fabs(s) * iDDomRad));
						//viscousStress[j] += points[i].g()* (1.0 - exp(-len[j] * iDDomRad / 2.0)*cosh(fabs(s)  * iDDomRad))  * iDDomRad / len[j]; // points[i].g()*(2.0 *  epscol* (1.0 - exp(-0.5 * len[j] / epscol))) / (PI * sqr(epscol) *len[j]);
						viscousStress[j] += 2.0 * locPoints[i].g()* (1.0 - exp(-len[j] * iDDomRad / 2.0)*cosh(fabs(s)  * iDDomRad))  * iDDomRad / (PI * len[j]);
					}
					break;

				}
			}//if d<50 len 
		}//for j
	}//for i

	if (id == 0)
	{
		selfI0.resize(pointsDb.vtx.size(), 0.0);
		selfI3.resize(pointsDb.vtx.size(), {0.0, 0.0});
	}

	MPI_Gatherv(locI0.data(), par.myLen, MPI_DOUBLE,          selfI0.data(), par.len.data(), par.disp.data(), MPI_DOUBLE,          0, W.getParallel().commWork);
	MPI_Gatherv(locI3.data(), par.myLen, Point2D::mpiPoint2D, selfI3.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
	
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
	const size_t nr = r.size();
	double*& dev_ptr_r = devRPtr; 
	std::vector<double>& loci0 = pointsDb.tmpI0;
	std::vector<Point2D>& loci3 = pointsDb.tmpI3;
	double*& dev_ptr_i0 = pointsDb.devI0Ptr;
	double*& dev_ptr_i3 = pointsDb.devI3Ptr;

	const int& id = W.getParallel().myidWork;
	parProp par = W.getParallel().SplitMPI(npt, true);

	double tCUDASTART = 0.0, tCUDAEND = 0.0;

	tCUDASTART = omp_get_wtime();

	if ((npt > 0) && (nr > 0))
	{
		cuCalculateSurfDiffVeloWake(par.myDisp, par.myLen, dev_ptr_pt, nr, dev_ptr_r, dev_ptr_i0, dev_ptr_i3, dev_ptr_rad);
			
		W.getCuda().CopyMemFromDev<double, 2>(par.myLen, dev_ptr_i3, (double*)&loci3[0]);
		W.getCuda().CopyMemFromDev<double, 1>(par.myLen, dev_ptr_i0, &loci0[0]);

		std::vector<Point2D> newI3;
		std::vector<double> newI0;
		
		if (id == 0)
		{
			newI3.resize(I3.size());
			newI0.resize(I0.size());
		}

		MPI_Gatherv(loci3.data(), par.myLen, Point2D::mpiPoint2D, newI3.data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, 0, W.getParallel().commWork);
		MPI_Gatherv(loci0.data(), par.myLen, MPI_DOUBLE, newI0.data(), par.len.data(), par.disp.data(), MPI_DOUBLE, 0, W.getParallel().commWork);

		if (id == 0)
			for (size_t q = 0; q < I3.size(); ++q)
			{
				I0[q] += newI0[q];
				I3[q] += newI3[q];
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

	for (size_t i = 0; i < r.size() - 1; i++)
	{
		v1 = r[i] - point;
		v2 = r[i+1] - point;
		angle += atan2(v1^v2, v1 * v2);
	}

	if (fabs(angle) < 0.1)
		return false;
	else return true;
}//IsPointInAirfoil(...)