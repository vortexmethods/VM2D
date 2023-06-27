/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Airfoil2D.cpp                                                    |
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
\brief Файл кода с описанием класса Airfoil
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#include "Airfoil2D.h"

#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "nummatrix.h"
#include "Passport2D.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "WakeDataBase2D.h"
#include "World2D.h"
#include "Wake2D.h"

#include "spline/spline.h"

using namespace VM2D;

// Конструктор
Airfoil::Airfoil(const World2D& W_, const size_t numberInPassport_)
	: W(W_), numberInPassport(numberInPassport_)
{ }


//Проверка, идет ли вершина i следом за вершиной j
bool Airfoil::isAfter(size_t i, size_t j) const
{
	return ((i == j + 1) || (i == 0 && j == r_.size() - 1));
}//isAfter(...)


//Определяет, находится ли точка с радиус-вектором r внутри габаритного прямоугольника профиля
bool Airfoil::isInsideGabarits(const Point2D& r) const
{
	return (r[0] <= upRight[0] && (r[0] >= lowLeft[0] && r[1] >= lowLeft[1] && r[1] <= upRight[1]));
}//isInsideGabarits(...)


//Определяет, находится ли точка с радиус-вектором r вне габаритного прямоугольника профиля
bool Airfoil::isOutsideGabarits(const Point2D& r) const
{
	return (r[0] > upRight[0] || (r[0] < lowLeft[0] || r[1] < lowLeft[1] || r[1] > upRight[1]));
}//isOutsideGabarits(...)

//Вычисление средних значений eps на панелях
void Airfoil::calcMeanEpsOverPanel()
{
	meanEpsOverPanel.clear();
	meanEpsOverPanel.resize(getNumberOfPanels());


	double midEps;

	const Boundary& bnd = W.getBoundary(numberInPassport);
	VortexesParams virtVortParams = W.getVelocity().virtualVortexesParams[numberInPassport];

	for (size_t i = 0; i < getNumberOfPanels(); ++i)
	{
		midEps = 0.0;
		for (int j = bnd.vortexBeginEnd[i].first; j < bnd.vortexBeginEnd[i].second; ++j)
			midEps += virtVortParams.epsastWake[j];
			//midEps += std::max(virtVortParams.epsastWake[j], 0.5 * len[i] / (bnd.vortexBeginEnd[i].second - bnd.vortexBeginEnd[i].first));

		midEps /= (bnd.vortexBeginEnd[i].second - bnd.vortexBeginEnd[i].first);
		meanEpsOverPanel[i] = midEps;
	}//for i
}//calcMeanEpsOverPanel()

//Считывание профиля из файла
void Airfoil::ReadFromFile(const std::string& dir) //загрузка профиля из файла, его поворот и масштабирование
{
	const AirfoilParams& param = W.getPassport().airfoilParams[numberInPassport];
	std::string filename = dir + param.fileAirfoil;
	std::ifstream airfoilFile;

	if (fileExistTest(filename, W.getInfo(), { "txt", "TXT" }))
	{
		std::stringstream airfoilFile(VMlib::Preprocessor(filename).resultString);

		VMlib::StreamParser airfoilParser(W.getInfo(), "airfoil file parser", airfoilFile);

		m = 0.0; //TODO
		J = 0.0; //TODO
		rcm = { 0.0, 0.0 }; //TODO
		phiAfl = 0.0;

		inverse = param.inverse;



		//*

		////////////////////////////////////////////////////////////////////////
		//Проверяем отсутствие "задвоенных" подряд идущих точек
		//airfoilParser.get("r", r_);
		std::vector<Point2D> rFromFile;
		airfoilParser.get("r", rFromFile);

		if (rFromFile.size() > 0)
		{
			if (param.requiredNPanels <= rFromFile.size()) {

				r_.reserve(rFromFile.size());
				r_.push_back(rFromFile[0]);
				for (size_t i = 1; i < rFromFile.size(); ++i)
					if ((rFromFile[i] - rFromFile[i - 1]).length2() > 1e-12)
						r_.push_back(rFromFile[i]);

				//Если первая совпадает с последней, то убираем последнюю
				if ((r_.back() - r_.front()).length2() < 1e-12)
					r_.resize(r_.size() - 1);
			}
			else
			{
				W.getInfo('e') << "Airfoil shape is given explicitely, it can not be automatically split!" << std::endl;
				exit(200);
			}
		}
		else
		{
			std::vector<GeomPoint> geomFromFile;

			airfoilParser.get("geometry", geomFromFile);

			if ((geomFromFile.front().r - geomFromFile.back().r).length2() < 1e-12)
				geomFromFile.resize(geomFromFile.size() - 1);

			size_t reqN = param.requiredNPanels;

			std::vector<double> L(geomFromFile.size());
			double totalLength = 0.0;
			std::vector<size_t> nc = {};
			std::vector<size_t> ni = {};

			for (size_t i = 0; i < geomFromFile.size(); ++i)
			{
				const Point2D& p1 = geomFromFile[i].r;
				const Point2D& p2 = geomFromFile[(i + 1) % geomFromFile.size()].r;

				L[i] = (p2 - p1).length();
				totalLength += L[i];
			}

			/*
			Point2D p1, p2, p3, cc, dd;
			double ie = 0;
			int numOfDivPnt = 12;
			double alpha;
			double rscale = 0.1;
			std::vector<Point2D> pnt(numOfDivPnt);
			std::vector<Point2D> newGeom(numOfDivPnt * geomFromFile.size());

			for (size_t i = 0; i < geomFromFile.size(); ++i)
			{
				p1 = (i == 0 ? geomFromFile[(int)geomFromFile.size() - 1].r : geomFromFile[i - 1].r);
				p2 = geomFromFile[i].r;
				p3 = geomFromFile[(i + 1) % geomFromFile.size()].r;

				Point2D dc = (p1 - p2).unit() + (p3 - p2).unit();
				ie = dc.dist2To(dc.proj(p1 - p2));
				if (ie < 1e-3) ie = 0;

				alpha = (PI - acos(((p3 - p2) & (p1 - p2)) / (sqrt(((p3 - p2) & (p3 - p2)) * ((p1 - p2) & (p1 - p2)))))) / (numOfDivPnt - 1);

				cc = p2 + (rscale * ie) * dc;

				dd = p2 + (cc - p2).proj(p1 - p2) - cc;
				for (int j = 0; j < numOfDivPnt; ++j)
					newGeom[i * numOfDivPnt + j] = cc + dd.rotated(j * alpha);
			}
			*/

			//Ищем первый сплайн
			size_t i = 0;
			std::vector<size_t> splineStart, splineFinish;

			while (i < geomFromFile.size())
			{
				while (i < geomFromFile.size() && !(geomFromFile[i].type == "c"))
					++i;

				if (i < geomFromFile.size())
				{
					splineStart.push_back(i);
					++i;
					//ищем второй конец
					while (i < geomFromFile.size() && !(geomFromFile[i].type == "c"))
						++i;

					splineFinish.push_back(i);
				}
			}

			if ((splineStart.size() > 0) && (splineStart[0] != 0))
				splineFinish.back() = splineStart[0];

			if (reqN == 0)
			{
				r_.reserve(geomFromFile.size());
				for (size_t i = 0; i < geomFromFile.size(); ++i)
					r_.push_back(geomFromFile[i].r);
			}
			else
			{
				double hRef = totalLength / reqN;

				std::vector<int> nPanels;

				bool cyclic = (splineStart.size() == 0);

				if (cyclic)
				{
					splineStart.push_back(0);
					splineFinish.push_back(geomFromFile.size());
				}

				for (size_t s = 0; s < splineStart.size(); ++s)
				{
					//Сплайн				
					std::vector<double> T, X, Y;

					size_t splineLegs = splineFinish[s] - splineStart[s];
					if (splineFinish[s] <= splineStart[s])
						splineLegs += geomFromFile.size();

					T.reserve(splineLegs + 1);
					X.reserve(T.size());
					Y.reserve(T.size());

					double tbuf = 0.0;

					for (size_t i = splineStart[s]; i <= splineStart[s] + splineLegs; ++i)
					{
						const Point2D& pt = geomFromFile[i % geomFromFile.size()].r;

						T.push_back(tbuf);
						X.push_back(pt[0]);
						Y.push_back(pt[1]);

						if ((i == splineStart[s]) && (splineFinish[s] - splineStart[s] == 1))
						{
							double halfL = (i < geomFromFile.size()) ? 0.5 * L[i] : 0.5 * (geomFromFile.front().r - geomFromFile.back().r).length();
							tbuf += halfL;
							T.push_back(tbuf);
							Point2D nextPt = geomFromFile[(i + 1) % geomFromFile.size()].r;

							X.push_back(0.5 * (pt[0] + nextPt[0]));
							Y.push_back(0.5 * (pt[1] + nextPt[1]));

							tbuf += halfL;
						}
						else
							tbuf += (i < geomFromFile.size()) ? L[i] : (geomFromFile.front().r - geomFromFile.back().r).length();
					}

					tk::spline s1, s2;
					if (cyclic)
					{
						s1.set_boundary(tk::spline::cyclic);
						s2.set_boundary(tk::spline::cyclic);
					}
					else
					{
						s1.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::second_deriv, 0.0);
						s2.set_boundary(tk::spline::second_deriv, 0.0, tk::spline::second_deriv, 0.0);
					}

					s1.set_points(T, X);
					s2.set_points(T, Y);

					int NumberOfPanelsSpline = (int)ceil(T.back() / hRef);
					double hSpline = T.back() / NumberOfPanelsSpline;

					for (int i = 0; i < NumberOfPanelsSpline; ++i)
						r_.push_back({ s1(hSpline * i), s2(hSpline * i) });

					nPanels.push_back(NumberOfPanelsSpline);

				}
			}
		}//else geometry

		//std::ofstream of("airfoil.txt");		
		//for (size_t i = 0; i < r_.size(); ++i)
		//	of << r_[i] << "," << std::endl;
		//of.close();

		if (r_.size() == 0)
		{
			W.getInfo('e') << "No points on airfoil contour!" << std::endl;
			exit(200);
		}

		if (inverse)
			std::reverse(r_.begin(), r_.end());

		v_.resize(0);
		for (size_t q = 0; q < r_.size(); ++q)
			v_.push_back({ 0.0, 0.0 });


		int nPossibleWays;
		int defaultNPossibleWays = 0;
		airfoilParser.get("nPossibleWays", nPossibleWays, &defaultNPossibleWays, false);

		possibleWays.resize(nPossibleWays);
		for (int q = 0; q < nPossibleWays; ++q)
			airfoilParser.get("possibleWay" + std::to_string(q + 1), possibleWays[q]);






		Move(param.basePoint);
		Scale(param.scale);

		double rotationAngle = param.angle;
		if (W.getPassport().geographicalAngles)
			rotationAngle = -param.angle - 0.5 * PI;

		Rotate(-rotationAngle);
		//в конце Rotate нормали, касательные и длины вычисляются сами

		gammaThrough.clear();
		gammaThrough.resize(r_.size(), 0.0);

		//Вычисляем площадь
		area = 0.0;
		for (size_t q = 0; q < r_.size(); ++q)
		{
			const Point2D& cntq = 0.5 * (getR(q) + getR(q + 1));
			const Point2D& drq = getR(q + 1) - getR(q);
			area += (cntq[0] * drq[1] - cntq[1] * drq[0]);
		}
		area *= 0.5;
		area = fabs(area);

		//Тестируем на "освещенность"
		wayToVertex.resize(getNumberOfPanels(), -1);

		for (size_t i = 0; i < getNumberOfPanels(); ++i)
		{
			const Point2D& dest = 0.5 * (getR(i) + getR(i + 1));

			bool flag = false;
			int wayNumber = -1;
			do
			{
				flag = false;
				++wayNumber;
				size_t j = 0;
				for (; j < getNumberOfPanels(); ++j)
				{
					const Point2D& aflRj = getR(j);
					const Point2D& aflRj1 = getR(j + 1);

					auto check = [aflRj, aflRj1](const Point2D& start, const Point2D& finish)
					{
						return ((((aflRj - start) ^ (finish - start)) * ((aflRj1 - start) ^ (finish - start)) <= 0) && \
							(((start - aflRj) ^ (aflRj1 - aflRj)) * ((finish - aflRj) ^ (aflRj1 - aflRj)) <= 0));
					};

					if (i == j)
						continue;

					for (size_t q = 0; (!flag) && (q < ((wayNumber == 0) ? 1 : possibleWays[wayNumber - 1].size() + 1)); ++q)
					{
						Point2D start = ((q == 0) ? rcm : possibleWays[wayNumber - 1][q - 1]);
						Point2D finish = (wayNumber == 0 || q == possibleWays[wayNumber - 1].size()) ? dest : possibleWays[wayNumber - 1][q];

						flag = (flag || check(start, finish));
					}

					if (flag)
						break;
				}//for j

				if (j == getNumberOfPanels())
				{
					wayToVertex[i] = wayNumber;
					break;
				}

			} while ((wayToVertex[i] == -1) && (wayNumber < possibleWays.size()));

			if (wayToVertex[i] == -1)
			{
				std::cout << "Way to vertex inside airfoil not found" << std::endl;
				//std::cout << "!!!" << std::endl;
				//std::cout << "q = " << q << ", path[q] = " << path[q] << std::endl;
				//std::cout << "i = " << i << std::endl;
				//std::cout << "rcm = " << rcm << std::endl;
				//std::cout << "dest = " << 0.5 * (getR(i) + getR(i + 1)) << std::endl;
				//exit(767);
			}

		}//for i
	}
}//ReadFromFile(...)

//Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек, обусловленных геометрией профиля, и вычисление вязкого трения
void Airfoil::GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3)
{
	std::vector<double> selfI0;
	std::vector<Point2D> selfI3;

	std::vector<double> currViscousStress;
	currViscousStress.resize(r_.size(), 0.0);


	selfI0.resize(pointsDb.vtx.size(), 0.0);
	selfI3.resize(pointsDb.vtx.size(), { 0.0, 0.0 });

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
	shared(domainRadius, pointsDb, selfI0, selfI3, currViscousStress, PI) \
	private(xi, xi_m, lxi, lxi_m, lenj_m, v0, q, new_n, mn, h, d, s, vec, vs, expon, domRad, iDDomRad) schedule(dynamic, DYN_SCHEDULE)
	for (int i = 0; i < pointsDb.vtx.size(); ++i)
	{
		const Vortex2D& vtxI = pointsDb.vtx[i];

		domRad = std::max(domainRadius[i], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());
		iDDomRad = 1.0 / domRad;

		for (size_t j = 0; j < r_.size(); ++j)
		{
			vs = 0.0;
			q = vtxI.r() - 0.5 * (getR(j) + getR(j + 1));
			vec = tau[j];

			s = q & vec;
			d = fabs(q & nrm[j]);

			if ((d < 50.0 * len[j]) && (fabs(s) < 50.0 * len[j]))
			{
				v0 = vec * len[j];

				if ((d > 5.0 * len[j]) || (fabs(s) > 5.0 * len[j]))
				{
					xi = q * iDDomRad;
					lxi = xi.length();

					expon = exp(-lxi) * len[j];
					mn = nrm[j] * expon;

					if (selfI0[i] != -PI * domRad)
					{
						selfI0[i] += (xi & mn) * (lxi + 1.0) / (lxi * lxi);
						selfI3[i] += mn;
					}

					vs = vtxI.g() * expon / (PI * sqr(meanEpsOverPanel[j]));
				}
				else if ((d >= 0.1 * len[j]) || (fabs(s) > 0.5 * len[j]))
				{
					vs = 0.0;
					double den = (fabs(s) < 0.5 * len[j]) ? d : (fabs(s) + d - 0.5 * len[j]);

					new_n = std::max(1, static_cast<int>(ceil(5.0 * len[j] / den)));
					//new_n = 100;
					h = v0 * (1.0 / new_n);

					for (int m = 0; m < new_n; ++m)
					{
						xi_m = (vtxI.r() - (getR(j) + h * (m + 0.5))) * iDDomRad;
						lxi_m = xi_m.length();

						lenj_m = len[j] / new_n;
						expon = exp(-lxi_m) * lenj_m;

						mn = nrm[j] * expon;
						if (selfI0[i] != -PI * domRad)
						{
							selfI0[i] += (xi_m & mn) * (lxi_m + 1.0) / (lxi_m * lxi_m);
							selfI3[i] += mn;
						}
						vs += expon;
					}//for m
					vs *= vtxI.g() / (PI * sqr(meanEpsOverPanel[j]));
				}
				else
				{

					selfI0[i] = -PI * domRad;
					double mnog = 2.0 * domRad * (1.0 - exp(-len[j] * iDDomRad / 2.0) * cosh(fabs(s) * iDDomRad));
					selfI3[i] = nrm[j] * mnog;
					vs = mnog * vtxI.g() / (PI * sqr(meanEpsOverPanel[j]));
				}
			}//if d<50 len 

#pragma omp atomic
			currViscousStress[j] += vs;
		}//for j
	}//for i


	if (&pointsDb == &(W.getWake()))
		for (size_t i = 0; i < viscousStress.size(); ++i)
		{
			viscousStress[i] += currViscousStress[i];
		}

	for (size_t i = 0; i < I0.size(); ++i)
	{
		domRad = std::max(domainRadius[i], W.getPassport().wakeDiscretizationProperties.getMinEpsAst());

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


void Airfoil::GetDiffVelocityI0I3ToWakeAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3)
{
	if (&pointsDb == &(W.getWake()))
	{
		viscousStress.clear();
		viscousStress.resize(r_.size(), 0.0);
	}

	GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(pointsDb, domainRadius, I0, I3);
}//GetDiffVelocityI0I3ToWakeAndViscousStresses(...)


#if defined(USE_CUDA)
void Airfoil::GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3)
{
	std::vector<double> newViscousStress;

	//Обнуление вязких напряжений
	if (&pointsDb == &(W.getWake()))
	{
		viscousStress.clear();
		viscousStress.resize(r_.size(), 0.0);
	}

	//CUDA-ядро вызывается 1 раз и учитывает влияние сразу всех профилей
	if ((numberInPassport == 0) && ((&pointsDb == &(W.getWake())) || (&pointsDb == &(W.getBoundary(0).virtualWake))))
	{
		size_t npt = pointsDb.vtx.size();

		if (&pointsDb == &(W.getBoundary(0).virtualWake))
		{
			for (size_t q = 1; q < W.getNumberOfAirfoil(); ++q)
				npt += W.getBoundary(q).virtualWake.vtx.size();
		}

		double*& dev_ptr_pt = pointsDb.devVtxPtr;
		double*& dev_ptr_rad = pointsDb.devRadPtr;
		double*& dev_ptr_meanEps = devMeanEpsOverPanelPtr;
		const size_t nr = r_.size();
		double*& dev_ptr_r = devRPtr;

		double*& dev_ptr_i0 = pointsDb.devI0Ptr;
		double*& dev_ptr_i3 = pointsDb.devI3Ptr;

		float*& dev_ptr_i0f = pointsDb.devI0fPtr;
		float*& dev_ptr_i3f = pointsDb.devI3fPtr;

		double minRad = W.getPassport().wakeDiscretizationProperties.getMinEpsAst();

		size_t nTotPan = 0;
		for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s)
			nTotPan += W.getAirfoil(s).getNumberOfPanels();


		std::vector<Point2D> newI3(npt); //(I3.size());
		std::vector<double> newI0(npt); //(I0.size());
		std::vector<Point2Df> newI3f(npt); //(I3.size());
		std::vector<float> newI0f(npt); //(I0.size());

		newViscousStress.resize(nTotPan);


		if ((npt > 0) && (nr > 0))
		{
			double*& dev_ptr_visstr = devViscousStressesPtr;

			std::vector<double> locvisstr(nTotPan);

			std::vector<double> zeroVec(nTotPan, 0.0);
			cuCopyFixedArray(dev_ptr_visstr, zeroVec.data(), nTotPan * sizeof(double), 101);



			/*
			cuCalculateSurfDiffVeloWake(npt, dev_ptr_pt, nTotPan, dev_ptr_r, dev_ptr_i0, dev_ptr_i3, dev_ptr_rad, dev_ptr_meanEps, minRad, dev_ptr_visstr);
			W.getCuda().CopyMemFromDev<double, 2>(npt, dev_ptr_i3, (double*)newI3.data());
			W.getCuda().CopyMemFromDev<double, 1>(npt, dev_ptr_i0, newI0.data());
			W.getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_visstr, newViscousStress.data());
			if (&pointsDb == &(W.getWake()))
			{
				for (size_t q = 0; q < I3.size(); ++q)
				{
					I0[q] += newI0[q];
					I3[q] += newI3[q];
				}
			}
			//*/

			//FAST 31-05
			//*
			double timings[7];
			BHcu::wrapperDiffusiveVeloI0I3((Vortex2D*)dev_ptr_pt, dev_ptr_i0f, (Point2Df*)dev_ptr_i3f, dev_ptr_rad, dev_ptr_r, W.getNonConstCuda().CUDAptrs, true, (int)npt, nTotPan, dev_ptr_visstr, timings, dev_ptr_meanEps, minRad, \
				W.getNonConstCuda().n_CUDA_bodies, (int)W.getNonConstCuda().n_CUDA_wake, 8);
			W.getCuda().CopyMemFromDev<float, 2>(npt, dev_ptr_i3f, (float*)newI3f.data());
			W.getCuda().CopyMemFromDev<float, 1>(npt, dev_ptr_i0f, newI0f.data());
			W.getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_visstr, newViscousStress.data());

			if (&pointsDb == &(W.getWake()))
			{
				for (size_t q = 0; q < I3.size(); ++q)
				{
					I0[q] += newI0f[q];
					I3[q] += newI3f[q];
				}
			}
			//*/

			size_t curGlobPnl = 0;
			for (size_t s = 0; s < W.getNumberOfAirfoil(); ++s)
			{
				std::vector<double>& tmpVisStress = W.getNonConstAirfoil(s).tmpViscousStresses;
				const size_t& np = W.getAirfoil(s).getNumberOfPanels();
				tmpVisStress.resize(0);
				tmpVisStress.insert(tmpVisStress.end(), newViscousStress.begin() + curGlobPnl, newViscousStress.begin() + curGlobPnl + np);
				curGlobPnl += np;
			}

		}
	}//if numberInPassport==0


	if ((&pointsDb == &(W.getWake())) || (&pointsDb == &(W.getBoundary(0).virtualWake)))
		for (size_t i = 0; i < viscousStress.size(); ++i)
			viscousStress[i] += tmpViscousStresses[i];

}//GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses
#endif

bool Airfoil::IsPointInAirfoil(const Point2D& point) const
{
	double sumAngle = 0.0;

	Point2D v1, v2;

	for (size_t i = 0; i < r_.size(); ++i)
	{
		v1 = getR(i) - point;
		v2 = getR(i + 1) - point;
		sumAngle += atan2(v1 ^ v2, v1 & v2);
	}

	if (fabs(sumAngle) < 0.1)
		return false;
	else
		return true;
}//IsPointInAirfoil(...)

//Вычисляет габаритный прямоугольник профиля
void Airfoil::GetGabarits(double gap)	//определение габаритного прямоугольника
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

	Point2D size = upRight - lowLeft;
	lowLeft -= size * gap;
	upRight += size * gap;
}//GetGabarits(...)

// Вычисление нормалей, касательных и длин панелей по текущему положению вершин
void Airfoil::CalcNrmTauLen()
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

//Перемещение профиля 
void Airfoil::Move(const Point2D& dr)	//перемещение профиля как единого целого на вектор dr
{
	for (size_t i = 0; i < r_.size(); ++i)
		r_[i] += dr;
	rcm += dr;

	for (size_t q = 0; q < possibleWays.size(); ++q)
		for (Point2D& pts : possibleWays[q])
			pts += dr;

	CalcNrmTauLen();
	GetGabarits();
}//Move(...)

//Поворот профиля 
void Airfoil::Rotate(double alpha)	//поворот профиля на угол alpha вокруг центра масс
{
	phiAfl += alpha;
	nummatrix<double, 2, 2> rotMatrix = { { cos(alpha), -sin(alpha) }, { sin(alpha), cos(alpha) } };

	for (size_t i = 0; i < r_.size(); ++i)
		r_[i] = rcm + (rotMatrix & (r_[i] - rcm));

	for (size_t q = 0; q < possibleWays.size(); ++q)
		for (Point2D& pts : possibleWays[q])
		{
			Point2D oldPts = pts;
			pts = rcm + (rotMatrix & (oldPts - rcm));
		}

	CalcNrmTauLen();
	GetGabarits();
}//Rotate(...)


//Масштабирование профиля
void Airfoil::Scale(const Point2D& factor)	//масштабирование профиля на коэффициент factor относительно центра масс
{
	for (size_t i = 0; i < r_.size(); ++i)
		r_[i] = rcm + Point2D{ factor[0] * (r_[i] - rcm)[0], factor[1] * (r_[i] - rcm)[1] };

	for (size_t q = 0; q < possibleWays.size(); ++q)
		for (Point2D& pts : possibleWays[q])
		{
			Point2D oldPts = pts;
			pts = rcm + Point2D{ factor[0] * (oldPts - rcm)[0], factor[1] * (oldPts - rcm)[1] };
		}

	CalcNrmTauLen();
	GetGabarits();
}//Scale(...)

//Вычисление коэффициентов матрицы A для расчета влияния панели на панель
std::vector<double> Airfoil::getA(size_t p, size_t i, const Airfoil& airfoil, size_t j) const
{
	std::vector<double> res(p * p, 0.0);

	if ((i == j) && (&airfoil == this))
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

	const auto& miq = W.getIQ(numberInPassport, airfoil.numberInPassport);

	res[0] = miq.first(i, j);

	if (p == 1)
		return res;

	res[1] = miq.first(i, airfoil.getNumberOfPanels() + j);
	res[2] = miq.first(getNumberOfPanels() + i, j);
	res[3] = miq.first(getNumberOfPanels() + i, airfoil.getNumberOfPanels() + j);

	if (p == 2)
		return res;

	if (p > 2)
		throw (-42);

	//Хотя сюда никогда и не попадем
	return res;
}//getA(...)

//Вычисление коэффициентов матрицы A для расчета влияния профиля самого на себя
void Airfoil::calcIQ(size_t p, const Airfoil& otherAirfoil, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& matrPair) const
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
	shared(otherAirfoil, self, aflNSelf, aflNOther, matrPair, p, IDPI, IQPI) \
	private(npI, npJ, alpha, lambda, p1, s1, p2, s2, di, dj, i00, i01, i10, i11, v00, v11, v01, v10) schedule(dynamic, DYN_SCHEDULE)
	for (int i = 0; i < (int)getNumberOfPanels(); ++i)
		for (int j = 0; j < (int)otherAirfoil.getNumberOfPanels(); ++j)
		{
			npI = getNumberOfPanels();
			npJ = otherAirfoil.getNumberOfPanels();

			if ((i == j) && self)
			{

				matrPair.first(i, j) = 0.0;
				matrPair.second(i, j) = 0.0;

				if (p == 2)
				{
					matrPair.first(i, npJ + j) = 0.0;
					matrPair.first(npI + i, j) = 0.0;

					matrPair.second(i, npJ + j) = -IQPI;
					matrPair.second(npI + i, j) = IQPI;

					matrPair.first(npI + i, npJ + j) = 0.0;
					matrPair.second(npI + i, npJ + j) = 0.0;

				}
				if (p > 2)
					throw (-42);
			}//if i==j
			else
			{
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
				matrPair.second(i, j) = i00 & taui;


				if (p == 2)
				{
					v01 = {
						0.5 / (dj.length()) * (((p1 + s1) & tauj) * VMlib::Omega(s1, taui, tauj) - s1.length2() * taui),
						-0.5 * di.length() / dj.length() * VMlib::Omega(s1 + p2, tauj, tauj)
					};

					i01 = IDPI / len[i] * (-((alpha[0] + alpha[2]) * v01[0] + (alpha[1] + alpha[2]) * v01[1]).kcross() \
						+ ((lambda[0] + lambda[2]) * v01[0] + (lambda[1] + lambda[2]) * v01[1]) - 0.5 * di.length() * tauj);

					matrPair.first(i, npJ + j) = i01 & nrm[i];
					matrPair.second(i, npJ + j) = i01 & taui;


					v10 = {
						-0.5 / di.length() * (((s1 + s2) & taui) * VMlib::Omega(s1, taui, tauj) - s1.length2() * tauj),
						0.5 * dj.length() / di.length() * VMlib::Omega(s1 + p2, taui, taui)
					};

					i10 = IDPI / len[i] * (-((alpha[0] + alpha[2]) * v10[0] + alpha[2] * v10[1]).kcross() \
						+ ((lambda[0] + lambda[2]) * v10[0] + lambda[2] * v10[1]) + 0.5 * dj.length() * taui);

					matrPair.first(npI + i, j) = i10 & nrm[i];
					matrPair.second(npI + i, j) = i10 & taui;


					v11 = {
						1.0 / (12.0 * di.length() * dj.length()) * (2.0 * (s1 & VMlib::Omega(s1 - 3.0 * p2, taui, tauj)) * VMlib::Omega(s1, taui, tauj) - s1.length2() * (s1 - 3.0 * p2)) - 0.25 * VMlib::Omega(s1, taui, tauj),
						-di.length() / (12.0 * dj.length()) * VMlib::Omega(di, tauj, tauj),
						-dj.length() / (12.0 * di.length()) * VMlib::Omega(dj, taui, taui)
					};

					i11 = IDPI / len[i] * (-((alpha[0] + alpha[2]) * v11[0] + (alpha[1] + alpha[2]) * v11[1] + alpha[2] * v11[2]).kcross()\
						+ (lambda[0] + lambda[2]) * v11[0] + (lambda[1] + lambda[2]) * v11[1] + lambda[2] * v11[2] \
						+ 1.0 / 12.0 * (dj.length() * taui + di.length() * tauj - 2.0 * VMlib::Omega(s1, taui, tauj)));

					matrPair.first(npI + i, npJ + j) = i11 & nrm[i];
					matrPair.second(npI + i, npJ + j) = i11 & taui;
				}


				if (p > 2)
					throw (-42);
			}//else(i == j)

		}//for(...)
}//getIQ(...)

void Airfoil::GetInfluenceFromVorticesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const
{
	W.getBoundary(numberInPassport).GetInfluenceFromVorticesToRectPanel(panel, ptr, count, panelRhs);
}//GetInfluenceFromVorticesToPanel(...)


//Вычисление влияния части подряд идущих источников из области течения на панель для правой части
void Airfoil::GetInfluenceFromSourcesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const
{
	W.getBoundary(numberInPassport).GetInfluenceFromSourcesToRectPanel(panel, ptr, count, panelRhs);
}//GetInfluenceFromSourcesToPanel(...)

//Вычисление влияния слоя источников конкретной прямолинейной панели на вихрь в области течения
void Airfoil::GetInfluenceFromSourceSheetToVortex(size_t panel, const Vortex2D& vtx, Point2D& vel) const
{
	W.getBoundary(numberInPassport).GetInfluenceFromSourceSheetAtRectPanelToVortex(panel, vtx, vel);
}//GetInfluenceFromSourceSheetToVortex(...)

//Вычисление влияния вихревых слоев (свободный + присоединенный) конкретной прямолинейной панели на вихрь в области течения
void Airfoil::GetInfluenceFromVortexSheetToVortex(size_t panel, const Vortex2D& vtx, Point2D& vel) const
{
	W.getBoundary(numberInPassport).GetInfluenceFromVortexSheetAtRectPanelToVortex(panel, vtx, vel);
}//GetInfluenceFromVortexSheetToVortex(...)

void Airfoil::GetInfluenceFromVInfToPanel(std::vector<double>& vInfRhs) const
{
	W.getBoundary(numberInPassport).GetInfluenceFromVInfToRectPanel(vInfRhs);
}//GetInfluenceFromVInfToPanel(...)
