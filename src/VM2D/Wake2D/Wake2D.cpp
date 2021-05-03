/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Wake2D.cpp                                                       |
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
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
*/

#if defined(_WIN32)
 #include <direct.h>
#endif

#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

#include "Wake2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Tree2D.h"
#include "Velocity2D.h"
#include "World2D.h"
#include <algorithm>

#include "Velocity2DBarnesHut.h"

using namespace VM2D;



bool Wake::MoveInside(const Point2D& newPos, const Point2D& oldPos, const Airfoil& afl, size_t& panThrough)
{
	//const double porog_r = 1e-12;
	
	//double minDist = 1.0E+10; //расстояние до пробиваемой панели
	panThrough = -1;

	//проверка габ. прямоугольника
	if (afl.isOutsideGabarits(newPos) && afl.isOutsideGabarits(oldPos))
		return false;

	//если внутри габ. прямоугольника - проводим контроль
	bool hit = false;

	//std::pair<double, double> gabPointX = std::minmax(newPos[0], oldPos[0]);
	//std::pair<double, double> gabPointY = std::minmax(newPos[1], oldPos[1]);

	//std::pair<double, double> gabPanelX, gabPanelY;

	//Проверка на пересечение
	//double r0 = 0, r1 = 0, r2 = 0, r3 = 0;

	//int cnt = 0;
	for (size_t j = 0; j < afl.getNumberOfPanels(); ++j)
	{
		const Point2D& aflRj = afl.getR(j);
		const Point2D& aflRj1 = afl.getR(j + 1);

//		gabPanelX = std::minmax(aflRj[0], aflRj1[0]);
//		gabPanelY = std::minmax(aflRj[1], aflRj1[1]);

//		if ((gabPointX.second >= gabPanelX.first) && (gabPanelX.second >= gabPointX.first) && (gabPointY.second >= gabPanelY.first) && (gabPanelY.second >= gabPointY.first))
//		{
			if ((((aflRj - oldPos) ^ (newPos - oldPos)) * ((aflRj1 - oldPos) ^ (newPos - oldPos)) <= 0) && \
				(((oldPos - aflRj) ^ (aflRj1 - aflRj)) * ((newPos - aflRj) ^ (aflRj1 - aflRj)) <= 0))
			{
				hit = true;
				panThrough = j;
				break;
			}
//		}
//		else { ++cnt; }
	}//for j

//	std::cout << cnt << std::endl;

	return hit;
}//MoveInside(...)


//bool Wake::MoveInside(const Point2D& newPos, const Point2D& oldPos, const Airfoil& afl, size_t& panThrough)
//{
//	const double porog_r = 1e-12;
//
//	double minDist = 1.0E+10; //расстояние до пробиваемой панели
//	panThrough = -1;
//
//	//проверка габ. прямоугольника
//	if (afl.isOutsideGabarits(newPos) && afl.isOutsideGabarits(oldPos))
//		return false;
//
//	//если внутри габ. прямоугольника - проводим контроль
//	bool hit = false;
//
//	//Определение прямой: Ax+By+D=0 - перемещение вихря
//	double A = newPos[1] - oldPos[1];
//	double B = oldPos[0] - newPos[0];
//	double D = oldPos[1] * newPos[0] - oldPos[0] * newPos[1];
//
//	double A1, B1, D1;
//
//	std::pair<double, double> gabPointX = std::minmax(newPos[0], oldPos[0]);
//	std::pair<double, double> gabPointY = std::minmax(newPos[1], oldPos[1]);
//
//	std::pair<double, double> gabPanelX, gabPanelY;
//
//	//Проверка на пересечение
//	double r0 = 0, r1 = 0, r2 = 0, r3 = 0;
//
//	for (size_t j = 0; j < afl.getNumberOfPanels(); ++j)
//	{
//		r0 = A * afl.getR(j)[0] + B * afl.getR(j)[1] + D;
//		r1 = A * afl.getR(j + 1)[0] + B * afl.getR(j + 1)[1] + D;
//
//		if (fabs(r0) < porog_r) r0 = 0.0;
//		if (fabs(r1) < porog_r) r1 = 0.0;
//
//		if (r0*r1 > 0)
//			continue;
//
//		//Определение прямой:A1x+B1y+D1=0 - панель
//		A1 = afl.getR(j + 1)[1] - afl.getR(j)[1];
//		B1 = afl.getR(j)[0] - afl.getR(j + 1)[0];
//		D1 = afl.getR(j)[1] * afl.getR(j + 1)[0] - afl.getR(j)[0] * afl.getR(j + 1)[1];
//
//		r2 = A1 * oldPos[0] + B1 * oldPos[1] + D1;
//		r3 = A1 * newPos[0] + B1 * newPos[1] + D1;
//
//		if (fabs(r2) < porog_r) r2 = 0.0;
//		if (fabs(r3) < porog_r) r3 = 0.0;
//
//		if (r2*r3 > 0)
//			continue;
//
//		hit = true;// пробила!
//		double d2 = (oldPos[0] - (B*D1 - D * B1) / (A*B1 - B * A1))*(oldPos[0] - (B*D1 - D * B1) / (A*B1 - B * A1)) + \
//			(oldPos[1] - (A1*D - D1 * A) / (A*B1 - B * A1))*(oldPos[1] - (A1*D - D1 * A) / (A*B1 - B * A1));
//
//		if (d2 < minDist)
//		{
//			minDist = d2;
//			panThrough = j;
//		}//if d2
//	}//for j
//
//
//	return hit;
//}//MoveInside(...)

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

	for (size_t i = 0; i < afl.getNumberOfPanels(); ++i)
	{
		Point2D v1, v2, vv;
		v1 = afl.getR(i) - newPos;
						
		v2 = afl.getR(i + 1) - newPos;
						
		vv = afl.getR(i + 1) - afl.getR(i);

		double dst = v1.length2() + v2.length2();
		if (dst < dist2)
		{
			dist2 = dst;
			panThrough = i;
		}

		cs = v1 & v2;
		sn = v1 ^ v2;
		
		angle += atan2(sn, cs);
	}//for i

	hit = ((angle > 3.14) || (angle < -3.14));

	return hit;
}//MoveInsideMovingBoundary(...)

//bool Wake::MoveInsideMovingBoundary(const Point2D& newPos, const Point2D& oldPos, const Airfoil& oldAfl, const Airfoil& afl, size_t& panThrough)
//{
//	panThrough = -1;
//
//	/// \todo сравнить производительности двух inside-ов
//
//	//проверка габ. прямоугольника
//	if (afl.isOutsideGabarits(newPos) && afl.isOutsideGabarits(oldPos))
//		return false;
//
//	bool hit = false;
//
//	double angle = 0;
//	double cs, sn;
//
//	double dist2 = 1000000000.0;
//
//	for (int i = 0; i < afl.getNumberOfPanels(); ++i)
//	{
//		Point2D v1, v2, vv;
//		v1 = afl.getR(i) - newPos;
//
//		v2 = afl.getR(i + 1) - newPos;
//
//		vv = afl.getR(i + 1) - afl.getR(i);
//
//		double dst = v1.length2() + v2.length2();
//		if (dst < dist2)
//		{
//			dist2 = dst;
//			panThrough = i;
//		}
//
//		cs = v1 * v2;
//		sn = v1 ^ v2;
//
//		angle += atan2(sn, cs);
//	}//for i
//
//	hit = ((angle > 3.00) || (angle < -3.00));
//
//	return hit;
//}//MoveInsideMovingBoundary(...)



//Проверка пересечения вихрями следа профиля при перемещении
void Wake::Inside(const std::vector<Point2D>& newPos, Airfoil& afl, bool isMoves, const Airfoil& oldAfl)
{
	//int id = W.getParallel().myidWork;

	WakeSynchronize();

	VMlib::parProp par = W.getParallel().SplitMPI(vtx.size());

	std::vector<Point2D> locNewPos;
	locNewPos.resize(par.myLen);
	
	std::vector<double> gamma;
	gamma.resize(afl.getNumberOfPanels(), 0.0);

	std::vector<double> locGamma;
	locGamma.resize(afl.getNumberOfPanels(), 0.0);

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


	//std::stringstream sss;
	//sss << "through_" << W.currentStep;
	//std::ofstream of(W.getPassport().dir + "dbg/" + sss.str());
	//for (size_t i = 0; i < gamma.size(); ++i)
	//	of << gamma[i] << std::endl;
	//of.close();

	if (W.getParallel().myidWork == 0)
	{
		for (size_t q = 0; q < through.size(); ++q)
		if (through[q] > -1)
		{
			locGamma[static_cast<size_t>(through[q])] += vtx[q].g();
			vtx[q].g() = 0.0;
		}
	}

	MPI_Reduce(locGamma.data(), gamma.data(), (int)afl.getNumberOfPanels(), MPI_DOUBLE, MPI_SUM, 0, W.getParallel().commWork);

	if (W.getParallel().myidWork == 0)
		afl.gammaThrough = gamma;	

}//Inside(...)

//Поиск ближайшего соседа
void Wake::GetPairs(int type)
{
	switch (W.getPassport().numericalSchemes.velocityComputation.second)
	{
	case 0:
		GetPairsBS(type);
		break;
	case 1:
		GetPairsBH(type);
		break;
	}
}

//Поиск ближайшего соседа
void Wake::GetPairsBS(int type)
{
	int id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(vtx.size());

	std::vector<int> locNeighb; //локальный массив соседей (для данного процессора)
	locNeighb.resize(par.myLen);

	Point2D Ri, Rk;

	double maxG = W.getPassport().wakeDiscretizationProperties.maxGamma;

#pragma omp parallel for default(none) shared(type, locNeighb, par, maxG) schedule(dynamic, DYN_SCHEDULE)
	for (int locI = 0; locI < par.myLen; ++locI)
	{
		int s = locI + par.myDisp;
		const Vortex2D& vtxI = vtx[s];
		
		locNeighb[locI] = 0;//по умолчанию

		bool found = false;

		double r2, r2test;

		const double& cSP = collapseScaleParameter;
		const double& cRBP = collapseRightBorderParameter;

		while ( (!found) && ( s + 1 < (int)vtx.size() ) )
		{
			s++;
			const Vortex2D& vtxK = vtx[s];

			r2 = dist2(vtxI.r(), vtxK.r());

			//линейное увеличение радиуса коллапса				    
			double mnog = std::max(1.0, /* 2.0 * */ (vtxI.r()[0] - cRBP) / cSP);

			r2test = sqr( W.getPassport().wakeDiscretizationProperties.epscol * mnog );

			if (type == 1)
				r2test *= 4.0; //Увеличение радиуса коллапса в 2 раза для коллапса вихрей разных знаков			

			//if (mnog > 1.0)    
			//	r2test = sqr(0.005*cSP); 			

			if (r2 < r2test)
			{
				switch (type)
				{
					case 0:
						found = ( fabs(vtxI.g()*vtxK.g()) != 0.0) && (fabs(vtxI.g() + vtxK.g()) < sqr(mnog) * maxG);
						break;
					case 1:
						found = (vtxI.g()*vtxK.g() < 0.0);
						break;
					case 2:
						found = (vtxI.g()*vtxK.g() > 0.0) && (fabs(vtxI.g() + vtxK.g()) < sqr(mnog) * maxG);
						break;
				}
			}//if r2 < r2_test
		}//while

		if (found)
			locNeighb[locI] = s;
	}//for locI

	if (id == 0)
		neighb.resize(vtx.size());

	MPI_Gatherv(locNeighb.data(), par.myLen, MPI_INT, neighb.data(), par.len.data(), par.disp.data(), MPI_INT, 0, W.getParallel().commWork);

}//GetPairsBS(...)


void Wake::GetPairsBH(int type)
{
	/*
	int id = W.getParallel().myidWork;

	if (id == 0)
	{
		W.BuildTree(PointType::wake);

		auto& tree = *(W.treeWake);

		neighb.resize(vtx.size());
		double maxG = W.getPassport().wakeDiscretizationProperties.maxGamma;

		int s;
#pragma omp parallel for private(s) shared(maxG) schedule(dynamic, 16)
		for (int i = 0; i < (int)tree.lowCells.size(); ++i)
		{
			tree.lowCells[i]->CalcABCDandCloseCellsToLowLevel(tree.rootCell, false);
			for (auto it = tree.lowCells[i]->itBegin; it != tree.lowCells[i]->itEnd; ++it)
			{
				const Vortex2D& vtxI = it.getVtx();

				bool found = false;

				double r2, r2test;

				const double& cSP = collapseScaleParameter;
				const double& cRBP = collapseRightBorderParameter;

				for(size_t k = 0; k < tree.lowCells[i]->closeCells.size() && (!found); ++k)
					for (auto itCl = tree.lowCells[i]->closeCells[k]->itBegin; itCl < tree.lowCells[i]->closeCells[k]->itEnd && (!found); ++itCl)
					{
						s = itCl.getNum();
						if (s > it.getNum())
						{
							const Vortex2D& vtxK = itCl.getVtx();

							r2 = dist2(vtxI.r(), vtxK.r());

							//линейное увеличение радиуса коллапса
							double mnog = std::max(1.0,  (vtxI.r()[0] - cRBP) / cSP);
							                          //тут м.б. умножение на 2, см. выше

							r2test = sqr(W.getPassport().wakeDiscretizationProperties.epscol * mnog);

							if (type == 1)
								r2test *= 4.0; //Увеличение радиуса коллапса в 2 раза для коллапса вихрей разных знаков			

							//if (mnog > 1.0)     
							//	r2test = sqr(0.005*cSP);

							if (r2 < r2test)
							{
								switch (type)
								{
								case 0:
									found = (fabs(vtxI.g() * vtxK.g()) != 0.0) && (fabs(vtxI.g() + vtxK.g()) < sqr(mnog) * maxG);
									break;
								case 1:
									found = (vtxI.g() * vtxK.g() < 0.0);
									break;
								case 2:
									found = (vtxI.g() * vtxK.g() > 0.0) && (fabs(vtxI.g() + vtxK.g()) < sqr(mnog) * maxG);
									break;
								}
							}//if r2 < r2_test
						}
					}
				if (found)
					neighb[it.getNum()] = s;
				else neighb[it.getNum()] = 0;
			}//for locI
		}
	}// if (id == 0)
	*/
}//GetPairsBH(...)


#if defined(USE_CUDA)
void Wake::GPUGetPairs(int type)
{
	size_t npt = vtx.size();
	//const int& id = W.getParallel().myidWork;
	VMlib::parProp par = W.getParallel().SplitMPI(npt, true);

	double tCUDASTART = 0.0, tCUDAEND = 0.0;
	
	tCUDASTART += omp_get_wtime();
	std::vector<int> tnei(npt, 0);
	neighb.resize(npt, 0);
	
	if (npt > 0)
	{
		cuCalculatePairs(par.myDisp, par.myLen, npt, devVtxPtr, devMeshPtr, devNeiPtr, 2.0*W.getPassport().wakeDiscretizationProperties.epscol, sqr(W.getPassport().wakeDiscretizationProperties.epscol), type);
		
		W.getCuda().CopyMemFromDev<int, 1>(par.myLen, devNeiPtr, &tnei[0], 30);

		std::vector<int> newNei;

		newNei.resize(neighb.size());

		MPI_Allgatherv(tnei.data(), par.myLen, MPI_INT, newNei.data(), par.len.data(), par.disp.data(), MPI_INT, W.getParallel().commWork);

		for (size_t q = 0; q < neighb.size(); ++q)
		{
			neighb[q] = newNei[q];			
			//std::cout << q << " " << vtx[q].r() << " " << vtx[q].g() << " " << neighb[q] << " " << vtx[neighb[q]].r() << " " << vtx[neighb[q]].g() << std::endl;		
		}
	}

	tCUDAEND += omp_get_wtime();

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

		if (W.getPassport().numericalSchemes.velocityComputation.second == 0)
		{
			const_cast<Gpu&>(W.getCuda()).RefreshWake(3);
			GPUGetPairs(type);
		}
		else
			GetPairs(type);
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

			for (size_t vt = 0; vt + 1 < vtx.size(); ++vt)
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

						for(size_t afl = 0; afl < W.getNumberOfAirfoil(); ++afl)
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
	/// \todo Пока профиль 1, расстояние от его центра; сделать от самого правого профиля
	Point2D zerovec = { 0.0, 0.0 };
#pragma omp parallel for default(none) shared(distFar2, zerovec) reduction(+:nFar)
	for (int i = 0; i <static_cast<int>(vtx.size()); ++i)
	{
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

	std::vector<Vortex2D/*, VM2D::MyAlloc<VMlib::Vortex2D>*/> newWake;

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
void Wake::Restruct()
{
	W.getTimestat().timeRestruct.first += omp_get_wtime();


	// Определение параметров, отвечающих за увеличение радиуса коллапса
	std::vector<double> rightBorder, horizSpan;
	rightBorder.reserve(W.getNumberOfAirfoil());
	horizSpan.reserve(W.getNumberOfAirfoil());

	for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
	{
		
		rightBorder.emplace_back(W.getAirfoil(q).upRight[0]);
		horizSpan.emplace_back(W.getAirfoil(q).upRight[0] - W.getAirfoil(q).lowLeft[0]);
	}

	if (W.getNumberOfAirfoil() > 0)
	{
		W.getNonConstWake().collapseRightBorderParameter = *std::max_element(rightBorder.begin(), rightBorder.end());
		W.getNonConstWake().collapseScaleParameter = *std::max_element(horizSpan.begin(), horizSpan.end());
	}
	else
	{
		W.getNonConstWake().collapseRightBorderParameter = 0.0;
		W.getNonConstWake().collapseScaleParameter = 1.0;
	}
#if defined(__CUDACC__) || defined(USE_CUDA)
	W.getNonConstCuda().setCollapseCoeff(W.getWake().collapseRightBorderParameter, W.getWake().collapseScaleParameter);
#endif 
	

	WakeSynchronize();

	Collaps(1, 1);
	Collaps(2, 1);

	RemoveFar();
	RemoveZero();

	W.getTimestat().timeRestruct.second += omp_get_wtime();
}//Restruct()


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
