/*!
\file
\brief Файл кода с описанием класса Wake
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#include "Wake.h"


//Считывание вихревого следа из файла 
void Wake::ReadFromFile(const std::string& dir)
{
	std::string filename = dir + param.fileWake;
	std::ifstream wakeFile;
	if (fileExistTest(filename, *(defaults::defaultPinfo), *(defaults::defaultPerr), "wake"))
	{
		wakeFile.open(filename);
		std::vector<std::string> vekKeys;
		StreamParser wakeParser(wakeFile);
				
		wakeParser.get("nv", nv);

		wakeParser.get("vtx", vtx);

		wakeFile.close();
	}
}//ReadFromFile(...)


//Сохранение вихревого следа в файл
void Wake::SaveKadr(const std::string& dir, int step) const
{
	std::string fname = dir+"/Kadr";
	if (step < 10) fname += "0";
	if (step < 100) fname += "0";
	if (step < 1000) fname += "0";
	if (step < 10000) fname += "0";
	
	std::ostringstream ss;
	ss << step;
	fname += ss.str();
	fname += ".txt";

	std::ofstream outfile;

	outfile.open(fname);
	outfile << nv << std::endl; //Сохранение числа вихрей в пелене

	for (int i = 0; i < nv; i++)
	{
		Point2D rrr = vtx[i].r();
		
		double xi = (vtx[i].r())[0];
		double yi = (vtx[i].r())[1];
		double gi = vtx[i].g();

		outfile << (int)(i) << " " << (double)(param.eps) << " " << xi << " " << yi << " " << "0.0" << " " << "0.0" << " " << "0.0" << " " << gi << std::endl;
		//нули пишутся для совместимости с трехмерной программой и обработчиком ConMDV	
	}//for i	
	outfile.close();
}//SaveKadr(...)


//MPI-синхронизация вихревого следа
void Wake::WakeSinchronize()
{
	int nV;
	if (parallel.myidWork == 0)
		nV = vtx.size();
	MPI_Bcast(&nV, 1, MPI_INT, 0, parallel.commWork);

	if (parallel.myidWork > 0)
		vtx.resize(nV);

	MPI_Bcast(vtx.data(), nV, Vortex2D::mpiVortex2D, 0, parallel.commWork);
}//WakeSinchronize()


//Проверка пересечения вихрями следа профиля при перемещении
std::vector<double> Wake::Inside(const std::vector<Point2D>& newPos, const Airfoil& afl)
{
	const double porog_r = 1e-12;
	int id = parallel.myidWork;

	parallel.SplitMPI(nv);

	std::vector<Point2D> locNewPos;
	locNewPos.resize(parallel.len[id]);
	
	std::vector<double> gamma;
	gamma.resize(afl.np, 0);

	std::vector<double> locGamma;
	locGamma.resize(afl.np, 0);

	MPI_Scatterv(newPos.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, locNewPos.data(), parallel.len[id], Point2D::mpiPoint2D, 0, parallel.commWork);

#pragma omp parallel for default(none) shared(locGamma, locNewPos, afl, id)
	for (int locI = 0; locI < parallel.len[id]; ++locI)
	{
		int i = parallel.disp[id] + locI;

		double minDist = 1.0E+10; //расстояние до пробиваемой панели
		int minN = -1;            //номер пробиваемой панели

		//проверка габ. прямоугольника
		if (afl.isOutsideGabarits(locNewPos[locI]) && afl.isOutsideGabarits(vtx[i].r()))
			continue;

		//если внутри габ. прямоугольника - проводим контроль
		bool hit = false;

		//Определение прямой: Ax+By+D=0 - перемещение вихря
		double A = locNewPos[locI][1] - vtx[i].r()[1];
		double B = vtx[i].r()[0] - locNewPos[locI][0];
		double D = vtx[i].r()[1] * locNewPos[locI][0] - vtx[i].r()[0] * locNewPos[locI][1];

		double A1, B1, D1;

		//Проверка на пересечение
		double r0 = 0, r1 = 0, r2 = 0, r3 = 0;
		for (size_t j = 0; j < afl.np; ++j)
		{
			r0 = A*afl.r[j][0] + B*afl.r[j][1] + D;
			r1 = A*afl.r[j + 1][0] + B*afl.r[j + 1][1] + D;

			if (fabs(r0) < porog_r) r0 = 0.0;
			if (fabs(r1)<porog_r) r1 = 0.0;

			if (r0*r1>0)
				continue;

			//Определение прямой:A1x+B1y+D1=0 - панель
			A1 = afl.r[j + 1][1] - afl.r[j][1];
			B1 = afl.r[j][0] - afl.r[j + 1][0];
			D1 = afl.r[j][1] * afl.r[j + 1][0] - afl.r[j][0] * afl.r[j + 1][1];

			r2 = A1*vtx[i].r()[0] + B1*vtx[i].r()[1] + D1;
			r3 = A1*locNewPos[locI][0] + B1*locNewPos[locI][1] + D1;

			if (fabs(r2) < porog_r) r2 = 0.0;
			if (fabs(r3) < porog_r) r3 = 0.0;

			if (r2*r3 > 0)
				continue;

			hit = true;// пробила!
			double d2 = (vtx[i].r()[0] - (B*D1 - D*B1) / (A*B1 - B*A1))*(vtx[i].r()[0] - (B*D1 - D*B1) / (A*B1 - B*A1)) + \
				        (vtx[i].r()[1] - (A1*D - D1*A) / (A*B1 - B*A1))*(vtx[i].r()[1] - (A1*D - D1*A) / (A*B1 - B*A1));
			if (d2 < minDist)
			{
				minDist = d2;
				minN = j;
			}//if d2
		}//for j

		//если есть протыкание
		if (hit)
		{
			locGamma[minN] = vtx[i].g();
			vtx[i].g() = 0.0;
			//cout << "i = " << i << " hit! " << endl;
		}
	}//for locI

	MPI_Reduce(locGamma.data(), gamma.data(), afl.np, MPI_DOUBLE, MPI_SUM, 0, parallel.commWork);

	return gamma;
}//Inside(...)