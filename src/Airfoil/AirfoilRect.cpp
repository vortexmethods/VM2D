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
void AirfoilRect::ReadFromFile(const std::string& dir, const AirfoilParams& param) //загрузка профиля из файла, его поворот и масштабирование
{
	std::string filename = dir + param.fileAirfoil;
	std::ifstream airfoilFile;
	if (fileExistTest(filename, *(defaults::defaultPinfo), *(defaults::defaultPerr), "airfoil"))
	{
		airfoilFile.open(filename);
		StreamParser airfoilParser(airfoilFile);
		airfoilFile.close();

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
void AirfoilRect::GetDiffVelocityToSetOfPointsAndViscousStresses(const std::vector<Vortex2D>& points, std::vector<double>& domainRadius, std::vector<Point2D>& velo, double epscol)
{
	std::vector<Point2D> selfVelo;

	viscousStress.clear();
	viscousStress.resize(np, 0.0);

	const int& id = parallel.myidWork;

	parallel.SplitMPI(points.size());

	std::vector<Point2D> locDiffVelo;
	locDiffVelo.resize(parallel.len[id]);

	//Локальные переменные для цикла
	double I0;
	Point2D I3;

	Point2D q, xi, xi_m, v0;
	double lxi, lq, lxi_m, lenj_m;

	Point2D mn;
	int new_n;
	Point2D h;

	Point2D vec;	
	double s, d;


//#pragma omp parallel for \
	default(none) \
	shared(locDiffVelo, domainRadius, points, id, epscol) \
	private(I0, I3, xi, xi_m, lxi, lxi_m, lenj_m, v0, q, new_n, mn, lq, h, d, s, vec)
	for (int k = 0; k < parallel.len[id]; ++k)
	{
		int i = k + parallel.disp[id];

		I0 = 0.0;
		I3 = { 0.0, 0.0 };

		//TODO пока профиль строго 1
		for (size_t j = 0; j < r.size() - 1; j++)
		{
			q = points[i].r() - 0.5 * (r[j] + r[j + 1]);
			vec = tau[j];

			s = q * vec;

			if (fabs(s) > 0.5 * len[j])
				d = sqrt(sqr(fabs(q*vec) - 0.5*len[j]) + sqr(q*nrm[j]));
			else
				d = fabs(q*nrm[j]);


			xi = q * (1.0 / domainRadius[i]);

			lxi = xi.length();
			lq = q.length();

			v0 = vec * len[j];

			if (d < 50.0 * len[j])
			{
				if (d > 5.0 * len[j])
				{
					mn = nrm[j] * exp(-lxi);
					I0 += xi * len[j] * mn * (lxi + 1.0) / xi.length2();
					I3 += mn * len[j];
					viscousStress[j] += points[i].g() * (mn.kcross() * tau[j]) / (PI * sqr(epscol));
				}
				else if ((d <= 5.0 * len[j]) && (d >= 0.001 * len[j]))
				{
					new_n = (int)(ceil(10.0 * len[j] / d)) + 1;
					h = v0 * (1.0 / new_n);

					for (int m = 0; m < new_n; m++)
					{
						xi_m = (points[i].r() - (r[j] + h * (m + 0.5))) * (1.0 / domainRadius[i]);
						lxi_m = xi_m.length();
						lenj_m = len[j] / new_n;
						mn = nrm[j] * exp(-lxi_m);
						I0 += xi_m* lenj_m * mn * (lxi_m + 1.0) / xi_m.length2();
						I3 += mn * lenj_m;
						viscousStress[j] += points[i].g() * (mn.kcross() * tau[j]) / (PI * sqr(epscol));
					}//for m
				}
				else if (d <= 0.001 * len[j])
				{
					I0 += PI * sqr(domainRadius[i]);
					I3 += 2.0 * nrm[j] * domainRadius[i] * (1.0 - exp(-0.5 * len[j] / domainRadius[i]));
				//	viscousStress[j] += points[i].g() * (mn.kcross() * tau[j]) / (PI * sqr(epscol));
				//	viscousStress[j] += points[i].g() * (2.0 * nrm[j] * domainRadius[i] * (1.0 - exp(-0.5 * len[j] / domainRadius[i]))).kcross() * tau[j] / (PI * sqr(epscol) *len[j]);
					viscousStress[j] += points[i].g() *(2.0 *  epscol* (1.0 - exp(-0.5 * len[j] / epscol))) / (PI * sqr(epscol) *len[j]);
				}
			}//if d<50 len 
		}//for j
	
		I0 *= domainRadius[i];
		I0 += 2.0 * PI * sqr(domainRadius[i]);

		locDiffVelo[k] = I3 * (1.0 / I0);
	}//for k
	
	if (id == 0)
		selfVelo.resize(points.size());

	MPI_Gatherv(locDiffVelo.data(), parallel.len[id], Point2D::mpiPoint2D, selfVelo.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];

};