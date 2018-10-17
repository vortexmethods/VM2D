/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Wake.cpp                                                         |
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
\brief Файл кода с описанием класса Wake
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#if defined(_WIN32)
#include <direct.h>
#endif

#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "Preprocessor.h"
#include "Wake.h"
#include "World2D.h"

//Считывание вихревого следа из файла 
void Wake::ReadFromFile(const std::string& dir)
{
	std::string filename = dir + W.getPassport().wakeDiscretizationProperties.fileWake;
	std::ifstream wakeFile;
	
	if (fileExistTest(filename, W.getInfo()))
	{
		std::stringstream wakeFile(Preprocessor(filename).resultString);
		
		LogStream XXX;
		StreamParser wakeParser(XXX, "vortex wake file parser", wakeFile);
		
		int nv;
		wakeParser.get("nv", nv);
		wakeParser.get("vtx", vtx);
	}
}//ReadFromFile(...)


//Сохранение вихревого следа в файл
void Wake::SaveKadr(const std::vector<Vortex2D>& outVtx, const std::string& dir, size_t step, timePeriod& time) const
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


	size_t numberNonZero = 0;
	
	for (size_t q = 0; q < outVtx.size(); ++q)
	{
		if (outVtx[q].g() != 0.0)
			numberNonZero++;
	}

#if defined(_WIN32)
	_mkdir((dir + "snapshots").c_str());
#else
	mkdir((dir + "snapshots").c_str(), S_IRWXU | S_IRGRP | S_IROTH);
#endif

	outfile.open(dir + "snapshots/" + fname);
//	PrintLogoToTextFile(outfile, dir + "snapshots/" + fname, "Positions and circulations of vortices in the wake");

//	PrintHeaderToTextFile(outfile, "Number of vortices");
	outfile << numberNonZero << std::endl; //Сохранение числа вихрей в пелене
//	outfile << std::endl << "// " << numberNonZero << std::endl; //Сохранение числа вихрей в пелене

//	PrintHeaderToTextFile(outfile, "x_i     y_i     G_i");

	for (size_t i = 0; i < outVtx.size(); i++)
	{
		const Point2D& r = outVtx[i].r();
		double gi = outVtx[i].g();

		if (gi != 0.0)
		{
			outfile << static_cast<int>(i) << " " << W.getPassport().wakeDiscretizationProperties.eps << " " << r[0] << " " << r[1] << " " << "0.0" << " " << "0.0" << " " << "0.0" << " " << gi << std::endl;
			//нули пишутся для совместимости с трехмерной программой и обработчиком ConMDV	
			//outfile << std::endl << xi << " " << yi << " " << gi;
		}
	}//for i	
	outfile.close();

	time.second = omp_get_wtime();
}//SaveKadr(...)



void Wake::SaveKadrVtk(const std::vector<Vortex2D>& outVtx, const std::string& dir, size_t step, timePeriod& time) const
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
	fname += ".vtk";

	std::ofstream outfile;


	size_t numberNonZero = 0;
	
	for (size_t q = 0; q < outVtx.size(); ++q)
	{
		//if (outVtx[q].g() != 0.0)
			numberNonZero++;
	}

#if defined(_WIN32)
	_mkdir((dir + "snapshots").c_str());
#else
	mkdir((dir + "snapshots").c_str(), S_IRWXU | S_IRGRP | S_IROTH);
#endif

	outfile.open(dir + "snapshots/" + fname);
	
	outfile << "# vtk DataFile Version 2.0" << std::endl;
	outfile << (dir + "snapshots/" + fname).c_str() << std::endl;
	outfile << "ASCII" << std::endl;
	outfile << "DATASET UNSTRUCTURED_GRID" << std::endl;
	outfile << "POINTS " << numberNonZero << " float" << std::endl;
	
	
	for (size_t i = 0; i < outVtx.size(); i++)
	{
		Point2D rrr = outVtx[i].r();
		
		double xi = (outVtx[i].r())[0];
		double yi = (outVtx[i].r())[1];
		double gi =outVtx[i].g();

		//if (gi != 0.0)
		{
			outfile << xi << " " << yi << " " << "0.0"  << std::endl;
		}
	}//for i
	
	outfile << "CELLS " << numberNonZero << " " << 2*numberNonZero << std::endl;
	for (size_t i=0; i<numberNonZero; ++i)
	    outfile << "1 " << i << std::endl;
	    
	outfile << "CELL_TYPES " << numberNonZero << std::endl;
	for (size_t i=0; i<numberNonZero; ++i)
	    outfile << "1" << std::endl;
	
	outfile << std::endl;
	outfile << "POINT_DATA " << numberNonZero << std::endl;
	outfile << "SCALARS Gamma float 1" << std::endl;
	outfile << "LOOKUP_TABLE default" << std::endl;
	
	for (size_t i = 0; i < outVtx.size(); i++)
	{
	    outfile << outVtx[i].g() << std::endl;
	}//for i
	
	
	outfile.close();

	time.second = omp_get_wtime();
}//SaveKadrVtk(...)



//MPI-синхронизация вихревого следа
void Wake::WakeSynchronize()
{
	int nV;
	if (W.getParallel().myidWork == 0)
		nV = (int)vtx.size();
	MPI_Bcast(&nV, 1, MPI_INT, 0, W.getParallel().commWork);

	if (W.getParallel().myidWork > 0)
		vtx.resize(nV);

	MPI_Bcast(vtx.data(), nV, Vortex2D::mpiVortex2D, 0, W.getParallel().commWork);
}//WakeSinchronize()



bool Wake::MoveInside(const Point2D& newPos, const Point2D& oldPos, const Airfoil& afl, size_t& panThrough)
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
}//MoveInside(...)

bool Wake::MoveInsideMovingBoundary(const Point2D& newPos, const Point2D& oldPos, const Airfoil& oldAfl, const Airfoil& afl, size_t& panThrough)
{
	panThrough = -1;

	/// \todo сравнить производительности двух inside-ов

	//проверка габ. прямоугольника
	if (afl.isOutsideGabarits(newPos) && afl.isOutsideGabarits(oldPos))
		return false;

	bool hit = false;

	double angle = 0;
	double cs, sn;

	double dist2 = 1000000000.0;

	for (int i = 0; i<afl.np; i++)
	{
		Point2D v1, v2, vv;
		v1 = afl.r[i] - newPos;
						
		v2 = afl.r[i + 1] - newPos;
						
		vv = afl.r[i + 1] - afl.r[i];

		double dst = v1.length2() + v2.length2();
		if (dst < dist2)
		{
			dist2 = dst;
			panThrough = i;
		}

		cs = v1 * v2;
		sn = v1 ^ v2;
		
		angle += atan2(sn, cs);
	}//for i

	hit = ((angle > 3.00) || (angle < -3.00));

	return hit;
}//MoveInsideMovingBoundary(...)



//Проверка пересечения вихрями следа профиля при перемещении
void Wake::Inside(const std::vector<Point2D>& newPos, Airfoil& afl, bool isMoves, const Airfoil& oldAfl)
{
	int id = W.getParallel().myidWork;

	WakeSynchronize();

	parProp par = W.getParallel().SplitMPI(vtx.size());

	std::vector<Point2D> locNewPos;
	locNewPos.resize(par.myLen);
	
	std::vector<double> gamma;
	gamma.resize(afl.np, 0.0);

	std::vector<double> locGamma;
	locGamma.resize(afl.np, 0.0);

	std::vector<int> through;
	if (W.getParallel().myidWork == 0)
		through.resize(par.totalLen);

	std::vector<int> locThrough;
	locThrough.resize(par.myLen, -1);

	MPI_Scatterv(const_cast<std::vector<Point2D> &>(newPos).data(), par.len.data(), par.disp.data(), Point2D::mpiPoint2D, locNewPos.data(), par.myLen, Point2D::mpiPoint2D, 0, W.getParallel().commWork);

#pragma omp parallel for default(none) shared(locGamma, locNewPos, locThrough, afl, oldAfl, par, isMoves) 
	for (int locI = 0; locI < par.myLen; ++locI)
	{
		int i = par.myDisp + locI;
		size_t minN;

		bool crit = isMoves ? MoveInsideMovingBoundary(locNewPos[locI], vtx[i].r(), oldAfl, afl, minN) : MoveInside(locNewPos[locI], vtx[i].r(), afl, minN);

		if (crit)
			locThrough[static_cast<size_t>(locI)] = static_cast<int>(minN);
	}//for locI	

	MPI_Gatherv(locThrough.data(), par.myLen, MPI_INT, through.data(), par.len.data(), par.disp.data(), MPI_INT, 0, W.getParallel().commWork);

	//std::ostringstream sss;
	//sss << "through_";
	//std::ofstream throughFile(sss.str());
	//for (size_t i = 0; i < gamma.size(); ++i)
	//	throughFile << gamma[i] << std::endl;
	//throughFile.close();

	/// \todo Только нулевой процессор или все?
	if (W.getParallel().myidWork == 0)
	{
		for (size_t q = 0; q < through.size(); ++q)
		if (through[q] > -1)
		{
			locGamma[static_cast<size_t>(through[q])] += vtx[q].g();
			vtx[q].g() = 0.0;
		}
	}

	MPI_Reduce(locGamma.data(), gamma.data(), (int)afl.np, MPI_DOUBLE, MPI_SUM, 0, W.getParallel().commWork);

	if (W.getParallel().myidWork == 0)
		afl.gammaThrough = gamma;	

}//Inside(...)



//Поиск ближайшего соседа
void Wake::GetPairs(int type)
{
	int id = W.getParallel().myidWork;
	parProp par = W.getParallel().SplitMPI(vtx.size());

	/// \todo Временно для профиля из 1000 панелей
	const double max_g = 0.03;		//максимальная циркуляция вихря, получаемого на первом шаге расчета
	const double coeff_max_g = 0.25; // коэффициент, определяющий максимально возможную циркуляцию вихря при коллапсе

	std::vector<int> locNeighb; //локальный массив соседей (для данного процессора)
	locNeighb.resize(par.myLen);

	Point2D Ri, Rk;

#pragma omp parallel for default(none) shared(type, locNeighb, par) schedule(dynamic, DYN_SCHEDULE)
	for (int locI = 0; locI < par.myLen; ++locI)
	{
		int s = locI + par.myDisp;
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

			r2test = sqr( W.getPassport().wakeDiscretizationProperties.epscol * std::max(1.0, vtxI.r()[0]) );

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

	MPI_Gatherv(locNeighb.data(), par.myLen, MPI_INT, neighb.data(), par.len.data(), par.disp.data(), MPI_INT, 0, W.getParallel().commWork);

}//GetPairs(...)


#if defined(USE_CUDA)
void Wake::GPUGetPairs(int type)
{
	size_t npt = vtx.size();
	const int& id = W.getParallel().myidWork;
	parProp par = W.getParallel().SplitMPI(npt, true);

	double tCUDASTART = 0.0, tCUDAEND = 0.0;
	
	tCUDASTART = omp_get_wtime();
	tmpNei.resize(npt, 0);
	neighb.resize(npt, 0);
	
	if (npt > 0)
	{
		cuCalculatePairs(par.myDisp, par.myLen, npt, devVtxPtr, devMeshPtr, devNeiPtr, 2.0*W.getPassport().wakeDiscretizationProperties.epscol, sqr(W.getPassport().wakeDiscretizationProperties.epscol), type);
		
		W.getCuda().CopyMemFromDev<int, 1>(par.myLen, devNeiPtr, &tmpNei[0]);

		std::vector<int> newNei;

		newNei.resize(neighb.size());

		MPI_Allgatherv(tmpNei.data(), par.myLen, MPI_INT, newNei.data(), par.len.data(), par.disp.data(), MPI_INT, W.getParallel().commWork);

		for (size_t q = 0; q < neighb.size(); ++q)
		{
			neighb[q] = newNei[q];			
			//std::cout << q << " " << vtx[q].r() << " " << vtx[q].g() << " " << neighb[q] << " " << vtx[neighb[q]].r() << " " << vtx[neighb[q]].g() << std::endl;		
		}
	}

	tCUDAEND = omp_get_wtime();

	//W.getInfo('t') << "GPU_Pairs: " << (tCUDAEND - tCUDASTART) << std::endl;
}
#endif


// Коллапс вихрей
int Wake::Collaps(int type, int times)
{
	int nHlop = 0; //общее число убитых вихрей

	//int loc_hlop = 0; //

	std::vector<bool> flag;	//как только вихрь сколлапсировался  flag = 1
	for (int z = 0; z < times; ++z)
	{
	
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_PAIRS))	
		//W.getInfo('t') << "nvtx (" << parallel.myidWork << ") = " << vtx.size() << std::endl;
		const_cast<Gpu&>(W.getCuda()).RefreshWake();
		GPUGetPairs(type);
#else
		GetPairs(type);
#endif		
		
		if (W.getParallel().myidWork == 0)
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
						size_t hitpan = -1;

						for(size_t afl = 0; afl < W.getNumberOfAirfoil(); afl++)
						{
							//проверим, не оказался ли новый вихрь внутри контура
							if (MoveInside(sumVtx.r(), vtxI.r(), W.getAirfoil(afl), hitpan) || MoveInside(sumVtx.r(), vtxK.r(), W.getAirfoil(afl), hitpan))
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


int Wake::RemoveFar()
{
	int nFar = 0;
	double distFar2 = sqr(W.getPassport().wakeDiscretizationProperties.distFar);
	/// \todo Пока профиль 1, расстояние от его центра
	Point2D zerovec = { 0.0, 0.0 };
#pragma omp parallel for default(none) shared(distFar2, zerovec) reduction(+:nFar)
	for (int i = 0; i <static_cast<int>(vtx.size()); ++i)
	{

	//	!!!!!!!!!!!!
	//	if (dist2(vtx[i].r(), (*airfoils[0]).rcm) > distFar2)
		if (dist2(vtx[i].r(), zerovec) > distFar2)
		{
			vtx[i].g() = 0.0;
			nFar++; 
		}
	}

	WakeSynchronize();
	return nFar;
}//RemoveFar()


size_t Wake::RemoveZero()
{
	const double porog_g = 1e-15;

	std::vector<Vortex2D> newWake;

	newWake.reserve(vtx.size());

	for (size_t q = 0; q < vtx.size(); ++q)
	if (fabs(vtx[q].g()) > porog_g)
		newWake.push_back(vtx[q]);

	size_t delta = vtx.size() - newWake.size();

	newWake.swap(vtx);

	WakeSynchronize();

	return delta;
}//RemoveZero()


//Реструктуризация вихревого следа
void Wake::Restruct(timePeriod& time)
{
	time.first = omp_get_wtime();

	WakeSynchronize();

	Collaps(1, 1);
	Collaps(2, 1);

	RemoveFar();
	RemoveZero();

	time.second = omp_get_wtime();
}//Restruct(...)


//bool Wake::isPointInsideAirfoil(const Point2D& pos, const Airfoil& afl)
//{
//	Point2D posFar = { 100.0, 100.0 };
//
//	Point2D r1, r2;
//
//	int nIntersec = 0;
//
//	for (int i = 0; i < afl.np; ++i)
//	{
//		r1 = afl.r[i];
//		r2 = afl.r[i+1];
//
//		if ((area(pos, posFar, r1) * area(pos, posFar, r2) <= -1.e-10) && (area(r1, r2, pos) * area(r1, r2, posFar) <= -1.e-10))
//			nIntersec++;
//	}
//	return nIntersec % 2;
//}
//
//double Wake::area(Point2D p1, Point2D p2, Point2D p3)
//{
//	return (p2[0] - p1[0]) * (p3[1] - p1[1]) - (p2[1] - p1[1]) * (p3[0] - p1[0]);
//}
