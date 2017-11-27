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

#if !defined(__linux__)
#include <direct.h>
#endif

#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

//Считывание вихревого следа из файла 
void Wake::ReadFromFile(const std::string& dir)
{
	std::string filename = dir + param.fileWake;
	std::ifstream wakeFile;

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

	if (fileExistTest(filename, Pinfo, Perr, "wake"))
	{
		wakeFile.open(filename);
		std::vector<std::string> vekKeys;
		StreamParser wakeParser(wakeFile);
		
		int nv;
		wakeParser.get("nv", nv);

		wakeParser.get("vtx", vtx);

		wakeFile.close();
	}
}//ReadFromFile(...)


//Сохранение вихревого следа в файл
void Wake::SaveKadr(const std::string& dir, int step, timePeriod& time) const
{
	time.first = omp_get_wtime();

	std::string fname = "Kadr";
	if (step < 10) fname += "0";
	if (step < 100) fname += "0";
	if (step < 1000) fname += "0";
	if (step < 10000) fname += "0";
	
	std::ostringstream ss;
	ss << step;
	fname += ss.str();
	fname += ".txt";

	std::ofstream outfile;


	int numberNonZero = 0;
	
	for (size_t q = 0; q < vtx.size(); ++q)
	{
		if (vtx[q].g() != 0.0)
			numberNonZero++;
	}

#if !defined(__linux__)
	_mkdir((dir + "snapshots").c_str());
#else
	mkdir((dir + "snapshots").c_str(), S_IRWXU | S_IRGRP | S_IROTH);
#endif

	outfile.open(dir + "snapshots/" + fname);
//	PrintLogoToTextFile(outfile, dir + "snapshots/" + fname, "Positions and circulations of vortices in the wake");

//	PrintHeaderToTextFile(outfile, "Number of vortices");
	outfile << vtx.size() << std::endl; //Сохранение числа вихрей в пелене
//	outfile << std::endl << "// " << numberNonZero << std::endl; //Сохранение числа вихрей в пелене

//	PrintHeaderToTextFile(outfile, "x_i     y_i     G_i");

	for (size_t i = 0; i < vtx.size(); i++)
	{
		Point2D rrr = vtx[i].r();
		
		double xi = (vtx[i].r())[0];
		double yi = (vtx[i].r())[1];
		double gi = vtx[i].g();

		if (gi != 0.0)
		{
			outfile << (int)(i) << " " << (double)(param.eps) << " " << xi << " " << yi << " " << "0.0" << " " << "0.0" << " " << "0.0" << " " << gi << std::endl;
			//нули пишутся для совместимости с трехмерной программой и обработчиком ConMDV	
			//outfile << std::endl << xi << " " << yi << " " << gi;
		}
	}//for i	
	outfile.close();

	time.second = omp_get_wtime();
}//SaveKadr(...)


//MPI-синхронизация вихревого следа
void Wake::WakeSynchronize()
{
	int nV;
	if (parallel.myidWork == 0)
		nV = vtx.size();
	MPI_Bcast(&nV, 1, MPI_INT, 0, parallel.commWork);

	if (parallel.myidWork > 0)
		vtx.resize(nV);

	MPI_Bcast(vtx.data(), nV, Vortex2D::mpiVortex2D, 0, parallel.commWork);
}//WakeSinchronize()


bool Wake::MoveInside(const Point2D& newPos, const Point2D& oldPos, const Airfoil& afl, int& panThrough)
{
	const double porog_r = 1e-12;
	
	double minDist = 1.0E+10; //расстояние до пробиваемой панели
	panThrough = -1;

	//проверка габ. прямоугольника
	if (afl.isOutsideGabarits(newPos) && afl.isOutsideGabarits(oldPos))
		return false;

	//если внутри габ. прямоугольника - проводим контроль
	bool hit = false;

	//Определение прямой: Ax+By+D=0 - перемещение вихря
	double A = newPos[1] - oldPos[1];
	double B = oldPos[0] - newPos[0];
	double D = oldPos[1] * newPos[0] - oldPos[0] * newPos[1];

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

		r2 = A1*oldPos[0] + B1*oldPos[1] + D1;
		r3 = A1*newPos[0] + B1*newPos[1] + D1;

		if (fabs(r2) < porog_r) r2 = 0.0;
		if (fabs(r3) < porog_r) r3 = 0.0;

		if (r2*r3 > 0)
			continue;

		hit = true;// пробила!
		double d2 = (oldPos[0] - (B*D1 - D*B1) / (A*B1 - B*A1))*(oldPos[0] - (B*D1 - D*B1) / (A*B1 - B*A1)) + \
			(oldPos[1] - (A1*D - D1*A) / (A*B1 - B*A1))*(oldPos[1] - (A1*D - D1*A) / (A*B1 - B*A1));
		
		if (d2 < minDist)
		{
			minDist = d2;
			panThrough = j;
		}//if d2
	}//for j

	return hit;
}

//Проверка пересечения вихрями следа профиля при перемещении
void Wake::Inside(const std::vector<Point2D>& newPos, Airfoil& afl)
{
	int id = parallel.myidWork;

	WakeSynchronize();

	parallel.SplitMPI(vtx.size());

	std::vector<Point2D> locNewPos;
	locNewPos.resize(parallel.myLen);
	
	std::vector<double> gamma;
	gamma.resize(afl.np, 0.0);

	std::vector<double> locGamma;
	locGamma.resize(afl.np, 0.0);

	std::vector<int> through;
	if (parallel.myidWork == 0)
		through.resize(parallel.totalLen);

	std::vector<int> locThrough;
	locThrough.resize(parallel.myLen, 0);
	
	MPI_Scatterv(newPos.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, locNewPos.data(), parallel.myLen, Point2D::mpiPoint2D, 0, parallel.commWork);

#pragma omp parallel for default(none) shared(locGamma, locNewPos, locThrough, afl, id, std::cout) 
	for (int locI = 0; locI < parallel.myLen; ++locI)
	{
		int i = parallel.myDisp + locI;
		int minN;
		if (MoveInside(locNewPos[locI], vtx[i].r(), afl, minN))
		{
			//std::cout << "Inside " << i << std::endl;
#pragma omp atomic
			locGamma[minN] += vtx[i].g();
			
			//vtx[i].g() = 0.0;
			locThrough[locI] = 1;

			//cout << "i = " << i << " hit! " << endl;
		}
	}//for locI

	MPI_Reduce(locGamma.data(), gamma.data(), afl.np, MPI_DOUBLE, MPI_SUM, 0, parallel.commWork);

	MPI_Gatherv(locThrough.data(), parallel.myLen, MPI_INT, through.data(), parallel.len.data(), parallel.disp.data(), MPI_INT, 0, parallel.commWork);

	//std::ostringstream sss;
	//sss << "through_";
	//std::ofstream throughFile(sss.str());
	//for (size_t i = 0; i < gamma.size(); ++i)
	//	throughFile << gamma[i] << std::endl;
	//throughFile.close();

	/// \todo Только нулевой процессор или все?
	if (parallel.myidWork == 0)
	{
		afl.gammaThrough = gamma;

		for (size_t q = 0; q < through.size(); ++q)
		if (through[q])
			vtx[q].g() = 0.0;
	}
}//Inside(...)


//Поиск ближайшего соседа
void Wake::GetPairs(int type)
{
	int id = parallel.myidWork;
	parallel.SplitMPI(vtx.size());

	/// \todo Временно для профиля из 1000 панелей
	const double max_g = 0.006;		//максимальная циркуляция вихря, получаемого на первом шаге расчета
	const double coeff_max_g = 0.5; // коэффициент, определяющий максимально возможную циркуляцию вихря при коллапсе

	std::vector<int> locNeighb; //локальный массив соседей (для данного процессора)
	locNeighb.resize(parallel.myLen);

	Point2D Ri, Rk;

#pragma omp parallel for default(none) shared(type, locNeighb, id)
	for (int locI = 0; locI < parallel.myLen; ++locI)
	{
		size_t s = locI + parallel.myDisp;
		const Vortex2D& vtxI = vtx[s];
		
		locNeighb[locI] = 0;//по умолчанию

		bool found = false;

		double r2, r2test;

		while ( (!found) && ( s < vtx.size() -1 ) )
		{
			s++;
			const Vortex2D& vtxK = vtx[s];

			r2 = dist2(vtxI.r(), vtxK.r());

			/// \todo Линейное увеличение радиуса коллапса, нужно сделать более универсальный алгоритм

			r2test = sqr( param.epscol * std::max(1.0, vtxI.r()[0]) );

			if (type == 1)
				r2test *= 4.0; //Увеличение радиуса коллапса в 2 раза для коллапса вихрей разных знаков			

			if (r2<r2test)
			{
				switch (type)
				{
					case 0: 
						found = ( fabs(vtxI.g()*vtxK.g()) != 0.0) && (fabs(vtxI.g() + vtxK.g()) < max_g*coeff_max_g);
						break;
					case 1:
						found = (vtxI.g()*vtxK.g() < 0.0);
						break;
					case 2:
						found = (vtxI.g()*vtxK.g() > 0.0) && (fabs(vtxI.g() + vtxK.g()) < max_g*coeff_max_g);
						break;
				}
			}//if r2<r2_test
		}//while

		if (found)
			locNeighb[locI] = s;
	}//for locI

	if (id == 0)
		neighb.resize(vtx.size());

	MPI_Gatherv(locNeighb.data(), parallel.myLen, MPI_INT, neighb.data(), parallel.len.data(), parallel.disp.data(), MPI_INT, 0, parallel.commWork);

}//GetPairs(...)


// Коллапс вихрей
int Wake::Collaps(int type, int times)
{
	int nHlop = 0; //общее число убитых вихрей

	//int loc_hlop = 0; //

	std::vector<bool> flag;	//как только вихрь сколлапсировался  flag = 1

	for (int z = 0; z < times; ++z)
	{
		GetPairs(type);
		if (parallel.myidWork == 0)
		{
			//loc_hlop = 0;//число схлопнутых вихрей

			flag.clear();
			flag.resize(vtx.size(), false);

			double sumAbsGam, iws;
			Point2D newPos;

			for (size_t vt = 0; vt < vtx.size() - 1; ++vt)
			{
				Vortex2D& vtxI = vtx[vt];

				if (!flag[vt])
				{
					int ssd = neighb[vt];
					if ((ssd != 0) && (!flag[ssd]))
					{
						Vortex2D& vtxK = vtx[ssd];

						flag[ssd] = true;

						Vortex2D sumVtx;
						sumVtx.g() = vtxI.g() + vtxK.g();

						switch (type)
						{
						case 0:
						case 2:
							sumAbsGam = fabs(vtxI.g()) + fabs(vtxK.g());

							iws = sumAbsGam > 1e-10 ? 1.0 / sumAbsGam : 1.0;

							sumVtx.r() = (vtxI.r()*fabs(vtxI.g()) + vtxK.r()*fabs(vtxK.g())) * iws;
							break;

						case 1:
							sumVtx.r() = (fabs(vtxI.g()) > fabs(vtxK.g())) ? vtxI.r() : vtxK.r();
							break;
						}

						bool fl_hit = true;
						//double Ch1[2];
						int hitpan = -1;

						for(size_t afl = 0; afl < airfoils.size(); afl++)
						{
							//проверим, не оказался ли новый вихрь внутри контура
							if (MoveInside(sumVtx.r(), vtxI.r(), *airfoils[afl], hitpan) || MoveInside(sumVtx.r(), vtxK.r(), *airfoils[afl], hitpan))
								fl_hit = false;
						}//for

						if (fl_hit)
						{
							vtxI = sumVtx;
							vtxK.g() = 0.0;
							nHlop++;
						}//if (fl_hit)

					}//if ((ssd!=0)&&(!flag[ssd]))
				}//if !flag[vt] 
			}//for vt
		}

		WakeSynchronize();
	}//for z

	return nHlop;
}//Collaps(...)


int Wake::KillFar()
{
	int nFar = 0;
	double distKill2 = sqr(param.distKill);
	///TODO Пока профиль 1, расстояние от его центра
	Point2D zerovec = { 0.0, 0.0 };
#pragma omp parallel for default(none) shared(distKill2,zerovec) reduction(+:nFar)
	for (int i = 0; i <static_cast<int>(vtx.size()); ++i)
	{

	//	!!!!!!!!!!!!
	//	if (dist2(vtx[i].r(), (*airfoils[0]).rcm) > distKill2)
		if (dist2(vtx[i].r(), zerovec) > distKill2)
		{
			vtx[i].g() = 0.0;
			nFar++; 
		}
	}

	WakeSynchronize();
	return nFar;
}


int Wake::RemoveZero()
{
	const double porog_g = 1e-12;

	std::vector<Vortex2D> newWake;

	newWake.reserve(vtx.size());

	for (size_t q = 0; q < vtx.size(); ++q)
	if (fabs(vtx[q].g()) > porog_g)
		newWake.push_back(vtx[q]);

	int delta = vtx.size() - newWake.size();

	newWake.swap(vtx);

	WakeSynchronize();

	return delta;
}


//Реструктуризация вихревого следа
void Wake::Restruct(timePeriod& time)
{
	time.first = omp_get_wtime();

	WakeSynchronize();

	Collaps(1, 1);
	Collaps(2, 1);
	KillFar();
	RemoveZero();

	time.second = omp_get_wtime();
}//Restruct()