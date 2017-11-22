/*!
\file
\brief Файл кода с описанием класса VelocityBiotSavart
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "VelocityBiotSavart.h"


//Вычисление конвективных скоростей в заданном наборе точек
void VelocityBiotSavart::CalcConvVeloToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo, std::vector<double>& domainRadius)
{
	std::vector<Point2D> selfVelo;
	
	const int& id = parallel.myidWork;

	parallel.SplitMPI(points.size());

	//Заполнение "своей" части массива скоростей
	double cft = IDPI;

	std::vector<Point2D> locConvVelo;
	locConvVelo.resize(parallel.len[id]);

	std::vector<double> locDomRadius;
	locDomRadius.resize(parallel.len[id]);

	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2eps, dst2;
	
	
#pragma omp parallel for default(none) shared(locConvVelo, locDomRadius, points, id, cft) private(velI, tempVel, dst2, dst2eps)
	for (int i = 0; i < parallel.len[id]; ++i)
	{
		double ee2[3] = { 10000.0, 10000.0, 10000.0 };
		if (wake.vtx.size() == 0)
		{
			ee2[0] = 0.0;
			ee2[1] = 0.0;
			ee2[2] = 0.0;
		}
		
		velI = { 0.0, 0.0 };
		
		const Point2D& posI = points[parallel.disp[id] + i].r();

		for (size_t j = 0; j < wake.vtx.size(); ++j)
		{
			const Point2D& posJ = wake.vtx[j].r();
			const double& gamJ = wake.vtx[j].g();

			tempVel = { 0.0, 0.0 };

			dst2 = dist2(posI, posJ);

			//Вычисление радиуса домена
			if (dst2 > 0)
			{
				if (dst2<ee2[0])
				{
					ee2[2] = ee2[1];
					ee2[1] = ee2[0];
					ee2[0] = dst2;
				}//if (dist2<ee2[0])
				else
				{
					if (dst2<ee2[1])
					{
						ee2[2] = ee2[1];
						ee2[1] = dst2;
					}// if (dist2<ee2[1])
					else
					if (dst2<ee2[2])
						ee2[2] = dst2;
				}//else
			}//if (dst2>0)


			dst2eps = std::max(dst2, wake.param.eps2);
			tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
			tempVel *= (gamJ / dst2eps);
			velI += tempVel;
		}

		velI *= cft;
		locConvVelo[i] = velI;

		locDomRadius[i] = std::max(sqrt((ee2[0] + ee2[1] + ee2[2]) / 3.0), 2.0*wake.param.epscol);
	}

	if (id == 0)
		selfVelo.resize(points.size());
	

	domainRadius.resize(points.size());

	MPI_Gatherv(locConvVelo.data(), parallel.len[id], Point2D::mpiPoint2D, selfVelo.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
	MPI_Allgatherv(locDomRadius.data(), parallel.len[id], MPI_DOUBLE, domainRadius.data(), parallel.len.data(), parallel.disp.data(), MPI_DOUBLE, parallel.commWork);

	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
	
}//CalcConvVeloToSetOfPoints(...)


//Вычисление диффузионных скоростей в заданном наборе точек
void VelocityBiotSavart::CalcDiffVeloToSetOfPoints(const std::vector<Vortex2D>& points, const std::vector<double>& domainRadius, const std::vector<Vortex2D>& vortices, std::vector<Point2D>& velo)
{

	std::vector<Point2D> selfVelo;

	const int& id = parallel.myidWork;

	parallel.SplitMPI(points.size());

	std::vector<Point2D> locDiffVelo;
	locDiffVelo.resize(parallel.len[id]);

	//Локальные переменные для цикла
	Point2D velI;
	Point2D I2;
	Point2D Rij;
	double I1, rij, expr;



#pragma omp parallel for default(none) shared(locDiffVelo, domainRadius, points, vortices, id) private(velI, I1, I2, Rij, rij, expr)
	for (int r = 0; r < parallel.len[id]; ++r)
	{
		int i = r + parallel.disp[id];
		const Vortex2D& vtxI = points[i];
		velI = { 0.0, 0.0 };

		I2 = { 0.0, 0.0 };
		I1 = 0.0;

		for (size_t j = 0; j < vortices.size(); ++j)
		{
			const Vortex2D& vtxJ = vortices[j];

			Rij = vtxI.r() - vtxJ.r();
			rij = Rij.length();

			if (rij > 1e-10)
			{
				expr = exp(-rij / domainRadius[i]);
				I2 += (vtxJ.g()* expr / rij) * Rij;
				I1 += vtxJ.g()*expr;
			}//if (rij>1e-6)
		}//for j

		if (fabs(I1) > 1e-10)
			velI = I2* (1.0 / (I1 * domainRadius[i]));

		locDiffVelo[r] = velI;

	} // for r

	//std::ostringstream sss;
	//sss << "velo_";
	//std::ofstream veloFile(sss.str());
	//for (int i = 0; i < domainRadius.size(); ++i)
	//	veloFile << domainRadius[i] << std::endl;
	//veloFile.close();

	if (id == 0)
		selfVelo.resize(points.size(), { 0.0, 0.0 });

	MPI_Gatherv(locDiffVelo.data(), parallel.len[id], Point2D::mpiPoint2D, selfVelo.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);


	if (id == 0)
	for (size_t i = 0; i < velo.size(); ++i)
		velo[i] += selfVelo[i];
}//CalcDiffVeloToSetOfPoints(...)