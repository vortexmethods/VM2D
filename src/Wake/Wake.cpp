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
		
		int nv;
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


	int numberNonZero = 0;
	for each (Vortex2D ivtx in vtx)
	{
		if (ivtx.g() != 0.0)
			numberNonZero++;
	}

	outfile.open(fname);
	//outfile << vtx.size() << std::endl; //Сохранение числа вихрей в пелене
	outfile << numberNonZero << std::endl; //Сохранение числа вихрей в пелене

	for (size_t i = 0; i < vtx.size(); i++)
	{
		Point2D rrr = vtx[i].r();
		
		double xi = (vtx[i].r())[0];
		double yi = (vtx[i].r())[1];
		double gi = vtx[i].g();

		if (gi != 0.0)
			outfile << (int)(i) << " " << (double)(param.eps) << " " << xi << " " << yi << " " << "0.0" << " " << "0.0" << " " << "0.0" << " " << gi << std::endl;
		//нули пишутся для совместимости с трехмерной программой и обработчиком ConMDV	
	}//for i	
	outfile.close();
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

	parallel.SplitMPI(vtx.size());

	std::vector<Point2D> locNewPos;
	locNewPos.resize(parallel.len[id]);
	
	std::vector<double> gamma;
	gamma.resize(afl.np, 0);

	std::vector<double> locGamma;
	locGamma.resize(afl.np, 0);

	MPI_Scatterv(newPos.data(), parallel.len.data(), parallel.disp.data(), Point2D::mpiPoint2D, locNewPos.data(), parallel.len[id], Point2D::mpiPoint2D, 0, parallel.commWork);

#pragma omp parallel for default(none) shared(locGamma, locNewPos, afl, id, std::cout) 
	for (int locI = 0; locI < parallel.len[id]; ++locI)
	{
		int i = parallel.disp[id] + locI;
		int minN;
		if (MoveInside(locNewPos[locI], vtx[i].r(), afl, minN))
		{
			//std::cout << "Inside " << i << std::endl;
#pragma omp atomic
			locGamma[minN] += vtx[i].g();
			
			vtx[i].g() = 0.0;
			//cout << "i = " << i << " hit! " << endl;
		}
	}//for locI

	MPI_Reduce(locGamma.data(), gamma.data(), afl.np, MPI_DOUBLE, MPI_SUM, 0, parallel.commWork);

	//std::ostringstream sss;
	//sss << "through_";
	//sss << ".txt";
	//std::ofstream throughFile(sss.str());
	//for (size_t i = 0; i < gamma.size(); ++i)
	//	throughFile << gamma[i] << std::endl;
	//throughFile.close();

	/// \todo Только нулевой процессор или все?
	if (parallel.myidWork == 0)
		afl.gammaThrough = gamma;
}//Inside(...)


//Поиск ближайшего соседа
void Wake::GetPairs(int type)
{
	int id = parallel.myidWork;
	parallel.SplitMPI(vtx.size());

	/// \todo Временно для профилей из 60 панелей
	const double max_g = 0.1;		//максимальная циркуляция вихря, получаемого на первом шаге расчета
	const double coeff_max_g = 0.2; // коэффициент, определяющий максимально возможную циркуляцию вихря при коллапсе

	std::vector<int> locNeighb; //локальный массив соседей (для данного процессора)
	locNeighb.resize(parallel.len[id]);

	Point2D Ri, Rk;

#pragma omp parallel for default(none) shared(type, locNeighb, id)
	for (int locI = 0; locI < parallel.len[id]; ++locI)
	{
		size_t s = locI + parallel.disp[id];
		const Vortex2D& vtxI = vtx[s];
		
		locNeighb[locI] = 0;//по умолчанию

		bool found = false;

		double r2, r2test;

		while ( (!found) && ( s < vtx.size() -1 ) )
		{
			s++;
			const Vortex2D& vtxK = vtx[s];

			r2 = dist2(vtxI.r(), vtxK.r());

			/// \todo Линейное увеличение радиуса коллапса
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

	MPI_Gatherv(locNeighb.data(), parallel.len[id], MPI_INT, neighb.data(), parallel.len.data(), parallel.disp.data(), MPI_INT, 0, parallel.commWork);

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

						///TODO Профиль 1
						if (airfoils.size() > 0)
						{
							//проверим, не оказался ли новый вихрь внутри контура
							if (MoveInside(sumVtx.r(), vtxI.r(), *airfoils[0], hitpan) || MoveInside(sumVtx.r(), vtxK.r(), *airfoils[0], hitpan))
								fl_hit = false;
						}//if airfoils.size()>0

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
	
#pragma omp parallel for default(none) shared(distKill2) reduction(+:nFar)
	for (int i = 0; i < static_cast<int>(vtx.size()); ++i)
	{
		if (dist2(vtx[i].r(), (*airfoils[0]).rcm) > distKill2)
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

	for each (Vortex2D vt in vtx)
	{
		if (fabs(vt.g()) > porog_g)
			newWake.push_back(vt);
	}

	int delta = vtx.size() - newWake.size();

	newWake.swap(vtx);

	WakeSynchronize();

	return delta;
}


//Реструктуризация вихревого следа
void Wake::Restruct()
{
	Collaps(1, 1);
	Collaps(2, 1);
	KillFar();
	RemoveZero();
}//Restruct()