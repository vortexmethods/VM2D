/*!
\file
\brief Файл кода с описанием класса AirfoilRect
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "AirfoilRect.h"


//Считывание профиля из файла
void AirfoilRect::ReadFromFile(const std::string& dir) //загрузка профиля из файла, его поворот и масштабирование
{
	const AirfoilParams& param = passport.airfoilParams[numberInPassport];
	std::string filename = dir + param.fileAirfoil;
	std::ifstream airfoilFile;

	std::ostream *Pinfo, *Perr;
	if (parallel.myidWork == 0)
	{
		Pinfo = defaults::defaultPinfo;
		Perr = defaults::defaultPerr;
	}
	else
	{
		Pinfo = nullptr;
		Perr = nullptr;
	}

	if (fileExistTest(filename, Pinfo, Perr, "airfoil"))
	{
		std::stringstream airfoilFile(Preprocessor(filename).resultString);

		StreamParser airfoilParser(airfoilFile);

		rcm = param.basePoint;
		m = 0.0; //TODO
		J = 0.0; //TODO
		
		airfoilParser.get("np", np);

		airfoilParser.get("r", r);

		//замыкаем --- в конец приписываем первую точку (для простоты)
		r.push_back(r[0]);

		v.resize(0);
		for (size_t q = 0; q < r.size(); ++q)
			v.push_back({ 0.0, 0.0 });

		Move(rcm);
		Scale(param.scale);
		Rotate(param.angle);
		//в конце Rotate нормали, касательные и длины вычисляются сами

		gammaThrough.clear();
		gammaThrough.resize(np, 0.0);
	}
}//ReadFromFile(...)

//Вычисление диффузионных скоростей в наборе точек, обусловленных геометрией профиля, и вычисление вязкого трения
void AirfoilRect::GetDiffVelocityToSetOfPointsAndViscousStresses(const std::vector<Vortex2D>& points, std::vector<double>& domainRadius, std::vector<Point2D>& velo)
{
	std::vector<Point2D> selfVelo;

	viscousStress.clear();
	viscousStress.resize(np, 0.0);

	const int& id = parallel.myidWork;

	parallel.SplitMPI(points.size());

	std::vector<Vortex2D> locPoints;
	locPoints.resize(parallel.myLen);

	MPI_Scatterv(points.data(), parallel.len.data(), parallel.disp.data(), Vortex2D::mpiVortex2D, \
		locPoints.data(), parallel.myLen, Vortex2D::mpiVortex2D, 0, parallel.commWork);

	std::vector<Point2D> locDiffVelo;
	locDiffVelo.resize(parallel.myLen);

	//Локальные переменные для цикла
	double I0;
	Point2D I3;

	Point2D q, xi, xi_m, v0;
	double lxi, lxi_m, lenj_m;

	Point2D mn;
	int new_n;
	Point2D h;

	Point2D vec;
	double s, d;

	double vs, i0;
	double iDDomRad, expon;

	double iDPIepscol2 = 1.0 / (PI * sqr(passport.wakeDiscretizationProperties.epscol));

	#pragma omp parallel for \
		default(none) \
		shared(locDiffVelo, domainRadius, locPoints, id, iDPIepscol2) \
		private(I0, i0, I3, xi, xi_m, lxi, lxi_m, lenj_m, v0, q, new_n, mn, h, d, s, vec, vs, expon, iDDomRad) schedule(dynamic, 1)
	for (int i = 0; i < parallel.myLen; ++i)
	{
		I0 = 0.0;
		I3 = { 0.0, 0.0 };

		i0 = 0.0;

		iDDomRad = 1.0 / domainRadius[i+parallel.myDisp];

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
					I0 += xi * mn * (lxi + 1.0) / xi.length2();
					I3 += mn;
					viscousStress[j] += locPoints[i].g() * expon * iDPIepscol2;
				}
				else if ((d <= 5.0 * len[j]) && (d >= 0.1 * len[j]))
				{
					vs = 0.0;
					new_n = (int)(ceil(10.0 * len[j] / d)) + 1;
					h = v0 * (1.0 / new_n);

					for (int m = 0; m < new_n; m++)
					{
						xi_m = (locPoints[i].r() - (r[j] + h * (m + 0.5))) * iDDomRad;
						lxi_m = xi_m.length();

						expon = exp(-lxi_m);
						lenj_m = len[j] / new_n;
						mn = lenj_m *nrm[j] * expon;
						I0 += xi_m*  mn * (lxi_m + 1.0) / xi_m.length2();
						I3 += mn;
						vs += locPoints[i].g() * expon;
					}//for m
					viscousStress[j] += vs / new_n * iDPIepscol2;
				}
				else if (d <= 0.1 * len[j])
				{
					i0 = PI * sqr(domainRadius[i+parallel.myDisp]);
					if (fabs(s) > 0.5 * len[j])
					{						
						I3 += 2.0 * nrm[j] * domainRadius[i + parallel.myDisp] * (exp(-fabs(s)  * iDDomRad) * sinh(len[j] * iDDomRad / 2.0));
						//	viscousStress[j] += locPoints[i].g() * (mn.kcross() * tau[j]) / (PI * sqr(epscol));
						//	viscousStress[j] += locPoints[i].g() * (2.0 * nrm[j] * domainRadius[i+parallel.myDisp] * (1.0 - exp(-0.5 * len[j] / domainRadius[i+parallel.myDisp]))).kcross() * tau[j] / (PI * sqr(epscol) *len[j]);
						//viscousStress[j] += locPoints[i].g()*(exp(-fabs(s)  * iDDomRad) * sinh(len[j] * iDDomRad / 2.0))  * iDDomRad / len[j]; //locPoints[i].g() *(2.0 *  epscol* (1.0 - exp(-0.5 * len[j] / epscol))) / (PI * sqr(epscol) *len[j]);
						viscousStress[j] += 2.0 * locPoints[i].g()*(exp(-fabs(s)  * iDDomRad) * sinh(len[j] * iDDomRad / 2.0))  * iDDomRad / (PI * len[j]); 
					}
					else
					{
						I3 += 2.0 * nrm[j] * domainRadius[i + parallel.myDisp] * (1.0 - exp(-len[j] * iDDomRad / 2.0)*cosh(fabs(s) * iDDomRad));
						//viscousStress[j] += points[i].g()* (1.0 - exp(-len[j] * iDDomRad / 2.0)*cosh(fabs(s)  * iDDomRad))  * iDDomRad / len[j]; // points[i].g()*(2.0 *  epscol* (1.0 - exp(-0.5 * len[j] / epscol))) / (PI * sqr(epscol) *len[j]);
						viscousStress[j] += 2.0 * locPoints[i].g()* (1.0 - exp(-len[j] * iDDomRad / 2.0)*cosh(fabs(s)  * iDDomRad))  * iDDomRad / (PI * len[j]);
					}
				}
			}//if d<50 len 
		}//for j

		if (i0 != 0.0) 
			I0 = i0;
		else
		{
			I0 *= domainRadius[i + parallel.myDisp];
			I0 += 2.0 * PI * sqr(domainRadius[i + parallel.myDisp]);
		}
		locDiffVelo[i] = I3 * (1.0 / I0);
	}//for k

	if (id == 0)
		selfVelo.resize(points.size());

	MPI_Gatherv(locDiffVelo.data(), parallel.myLen, Point2D::mpiPoint2D, selfVelo.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
	
	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];

};