/*!
\file
\brief Файл кода с описанием класса VelocityDirect
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "VelocityDirect.h"


//Вычисление конвективных скоростей вихрей в вихревом следе
void VelocityDirect::CalcConvVelo(double dt)
{
	const int& id = parallel.myidWork;

	parallel.SplitMPI(wake.nv);

	//Заполнение "своей" части массива convVelo
	double cft = IDPI*dt;

	std::vector<Point2D> locConvVelo;
	locConvVelo.resize(parallel.len[id]);

	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2;

#pragma omp parallel for default(none) shared(locConvVelo, id, cft) private(velI, tempVel, dst2)
	for (int i = 0; i < parallel.len[id]; ++i)
	{
		velI = { 0.0, 0.0 };
		
		const Point2D& posI = wake.vtx[parallel.disp[id] + i].r();

		for (int j = 0; j < wake.nv; ++j)
		{
			const Point2D& posJ = wake.vtx[j].r();
			const double& gamJ = wake.vtx[j].g();

			tempVel = { 0.0, 0.0 };

			dst2 = std::max(dist2(posI, posJ), wake.param.eps);
			tempVel = { -(posI[1] - posJ[1]), (posI[0] - posJ[0]) };
			tempVel *= (gamJ / dst2);
			velI += tempVel;
		}

		velI *= cft;
		locConvVelo[i] = velI;
	}

	if (id == 0)
		convVelo.resize(wake.nv);

	MPI_Gatherv(locConvVelo.data(), parallel.len[id], Point2D::mpiPoint2D, convVelo.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, 0, parallel.commWork);
}//CalcConvVelo(...)

