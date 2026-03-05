/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
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
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#if defined(_WIN32)
 #include <direct.h>
#endif

#include "Wake2D.h"

#include "knnCPU.h"
#include "cuKnn.cuh"

#include "treeKernels.cuh"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "World2D.h"

using namespace VM2D;

bool Wake::MoveInside(const Point2D& newPos, const Point2D& oldPos, const Airfoil& afl, size_t& panThrough) const
{
	//const double porog_r = 1e-12;
	
	//double minDist = 1.0E+10; //расстояние до пробиваемой панели
	panThrough = -1;

	//проверка габ. прямоугольника
	
	if (afl.isOutsideGabarits(newPos) && afl.isOutsideGabarits(oldPos))
	{ 
		++W.gabb;
		return false;
	}

	//если внутри габ. прямоугольника - проводим контроль
	bool hit = false;

	
	for (size_t j = 0; j < afl.getNumberOfPanels(); ++j)
	{
		const Point2D& aflRj = afl.getR(j);
		const Point2D& aflRj1 = afl.getR(j + 1);

		++W.checkPan;

		if ((((aflRj - oldPos) ^ (newPos - oldPos)) * ((aflRj1 - oldPos) ^ (newPos - oldPos)) <= 0) && \
			(((oldPos - aflRj) ^ (aflRj1 - aflRj)) * ((newPos - aflRj) ^ (aflRj1 - aflRj)) <= 0))
		{
			hit = true;
			panThrough = j;						
			break;
		}
	}//for j


	
	if (hit)
		++W.check01;
	else
		++W.check02;

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


bool Wake::MoveInsideMovingBoundary(const Point2D& newPos, const Point2D& oldPos, const AirfoilGeometry& oldAfl, const Airfoil& afl, size_t& panThrough) const
{
	panThrough = -1;

	/// \todo сравнить производительности двух inside-ов

	//проверка габ. прямоугольника
	if (!afl.inverse && afl.isOutsideGabarits(newPos) && afl.isOutsideGabarits(oldPos))
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
	
	if (afl.inverse)
		hit = !hit;

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
//		cs = (v1 & v2);
//		sn = (v1 ^ v2);
//
//		angle += atan2(sn, cs);
//	}//for i
//
//	hit = ((angle > 3.00) || (angle < -3.00));
//
//	return hit;
//}//MoveInsideMovingBoundary(...)



//Проверка пересечения вихрями следа профиля при перемещении
void Wake::Inside(const std::vector<Point2D>& newPos, Airfoil& afl, bool isMoves, const AirfoilGeometry& oldAfl)
{
	std::vector<double> gamma;
	gamma.resize(afl.getNumberOfPanels(), 0.0);

	std::vector<int> through;
	through.resize(vtx.size(), -1);

//#pragma omp parallel for default(none) shared(afl, oldAfl, isMoves, through, newPos) 
	for (int i = 0; i < (int)vtx.size(); ++i)
	{		
		size_t minN;

		bool crit = isMoves ? MoveInsideMovingBoundary(newPos[i], vtx[i].r(), oldAfl, afl, minN) : MoveInside(newPos[i], vtx[i].r(), afl, minN);
		
		if (crit)
			through[i] = (int)minN;
	}//for i	


	//std::stringstream sss;
	//sss << "through_" << W.currentStep;
	//std::ofstream of(W.getPassport().dir + "dbg/" + sss.str());
	//for (size_t i = 0; i < gamma.size(); ++i)
	//	of << gamma[i] << std::endl;
	//of.close();

	//std::stringstream sss;
	//sss << "through_" << W.currentStep;
	//std::ofstream of(W.getPassport().dir + "dbg/" + sss.str());
	//for (size_t i = 0; i < gamma.size(); ++i)
	//	of << through[i] << std::endl;
	//of.close();

	for (size_t q = 0; q < through.size(); ++q)
	if (through[q] > -1)
	{
		gamma[through[q]] += vtx[q].g();
		vtx[q].g() = 0.0;
	}

	afl.gammaThrough = gamma;	
}//Inside(...)

//Поиск ближайшего соседа
void Wake::GetPairs(int type)
{
	GetPairsBS(type);
}

//Поиск ближайшего соседа
void Wake::GetPairsBS(int type)
{
	neighb.resize(vtx.size(), 0);
	
	Point2D Ri, Rk;

	double maxG = W.getPassport().wakeDiscretizationProperties.maxGamma;

#pragma omp parallel for default(none) shared(type, maxG) schedule(dynamic, DYN_SCHEDULE)
	for (int i = 0; i < vtx.size(); ++i)
	{
		neighb[i] = 0;
		int s = i;
		const Vortex2D& vtxI = vtx[i];	

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

			if (r2 < r2test)
			{
				switch (type)
				{
					case 0:
						found = ( vtxI.g()*vtxK.g() != 0.0) && (fabs(vtxI.g() + vtxK.g()) < sqr(mnog) * maxG);
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
			neighb[i] = s;
	}//for locI	
}//GetPairsBS(...)


//Поиск ближайшего соседа
void Wake::GetPairsClosestNeib(int type)
{
	neighb.resize(vtx.size(), 0);

	Point2D Ri, Rk;

	double maxG = W.getPassport().wakeDiscretizationProperties.maxGamma;

#pragma omp parallel for default(none) shared(type, maxG) schedule(dynamic, DYN_SCHEDULE)
	for (int i = 0; i < vtx.size(); ++i)
	{
		neighb[i] = 0;
		int s = i;
		const Vortex2D& vtxI = vtx[i];

		bool found = false;

		double r2, r2test = 1e+10;

		const double& cSP = collapseScaleParameter;
		const double& cRBP = collapseRightBorderParameter;

		while (/*(!found) &&*/ (s + 1 < (int)vtx.size()))
		{
			s++;
			const Vortex2D& vtxK = vtx[s];

			r2 = dist2(vtxI.r(), vtxK.r());

			if (r2 < r2test)
			{
				r2test = r2;
				neighb[i] = s;
			}//if r2 < r2_test
		}//while		
			
	}//for locI	
}//GetPairsClosestNeib(...)



#if defined(USE_CUDA)
void Wake::GPUGetPairs(int type)
{
	size_t npt = vtx.size();

	double tCUDASTART = 0.0, tCUDAEND = 0.0;
	
	tCUDASTART += omp_get_wtime();
	std::vector<int> tnei(npt, 0);
	neighb.resize(npt, 0);
	
	if (npt > 0)
	{
		cuCalculatePairs(npt, devVtxPtr, devMeshPtr, devNeiPtr, 2.0*W.getPassport().wakeDiscretizationProperties.epscol, sqr(W.getPassport().wakeDiscretizationProperties.epscol), type);
		W.getCuda().CopyMemFromDev<int, 1>(npt, devNeiPtr, neighb.data(), 30);		
	}

	tCUDAEND += omp_get_wtime();

	//W.getInfo('t') << "GPU_Pairs: " << (tCUDAEND - tCUDASTART) << std::endl;
}
#endif

#if defined(USE_CUDA)
void Wake::GPUGetPairsClosestNeib(int type)
{
	size_t npt = vtx.size();

	double tCUDASTART = 0.0, tCUDAEND = 0.0;

	tCUDASTART += omp_get_wtime();
	std::vector<int> tnei(npt, 0);
	neighb.resize(npt, 0);

	if (npt > 0)
	{
		cuCalculatePairsClosestNeib(npt, devVtxPtr, devMeshPtr, devNeiPtr, 2.0 * W.getPassport().wakeDiscretizationProperties.epscol, sqr(W.getPassport().wakeDiscretizationProperties.epscol), type);

		W.getCuda().CopyMemFromDev<int, 1>(npt, devNeiPtr, neighb.data(), 30);
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
		W.timerMerging.reset();
		W.timerMerging.start();
		//if (W.getPassport().numericalSchemes.velocityComputation.second == 0)
		{
			const_cast<Gpu&>(W.getCuda()).RefreshWake(3);
			
			double NeibStart = omp_get_wtime();
			GPUGetPairs(type);
			double NeibFinish = omp_get_wtime();
			std::cout << "GPU_direct_closest_neib = " << NeibFinish - NeibStart << std::endl;
		}
		W.timerMerging.stop();

		//else
		//{
		//	double NeibStart = omp_get_wtime();
		//	GetPairs(type);
		//	//GetPairsClosestNeib(type);
		//	double NeibFinish = omp_get_wtime();
		//	std::cout << "CPU_direct_closest_neib = " << NeibFinish - NeibStart << std::endl;					
		//}
#else
		GetPairs(type);
#endif			

		
		//std::ofstream neibFile(W.getPassport().dir + "neib" + std::to_string(W.currentStep));
		//for (size_t q = 0; q < neighb.size(); ++q)
		//	neibFile << q << " " << neighb[q] << std::endl;
		//neibFile.close();
		

		//loc_hlop = 0;//число схлопнутых вихрей
//*
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

				if ((ssd < 0) || (ssd >= flag.size()))
					std::cout << "ssd = " << ssd << ", flag.size() = " << flag.size() << ", vtx.size() = " << vtx.size() << std::endl;

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

						sumVtx.r() = (vtxI.r() * fabs(vtxI.g()) + vtxK.r() * fabs(vtxK.g())) * iws;
						break;

					case 1:
						sumVtx.r() = (fabs(vtxI.g()) > fabs(vtxK.g())) ? vtxI.r() : vtxK.r();
						break;
					}

					bool fl_hit = true;
					//double Ch1[2];
					size_t hitpan = -1;

					for (size_t afl = 0; afl < W.getNumberOfAirfoil(); ++afl)
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
//*/
	}//for z

	return nHlop;
}//Collaps(...)

// Коллапс вихрей
int Wake::CollapsNew(int type, int times)
{
	int nHlop = 0; //общее число убитых вихрей

	const double& cSP = collapseScaleParameter;
	const double& cRBP = collapseRightBorderParameter;

	double maxG = W.getPassport().wakeDiscretizationProperties.maxGamma;

	//int loc_hlop = 0; // neigh

	std::vector<bool> flag;	//как только вихрь сколлапсировался  flag = 1
	
	for (int z = 0; z < times; ++z)
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

				for (int s = 0; s < knbForRestruct; ++s)
				{
					//int ssd = neighb[vt];
					int ssd = neighbNew[vt * knbForRestruct + s];
					if (ssd == 0)
						continue;

					Vortex2D& vtxK = vtx[ssd];

					if (/*(ssd != 0) &&*/ (!flag[ssd])) //ssd==0 always true
					{
						

						Vortex2D sumVtx;
						sumVtx.g() = vtxI.g() + vtxK.g();

						switch (type)
						{
						case 0:
						case 2:
							sumAbsGam = fabs(vtxI.g()) + fabs(vtxK.g());
							iws = sumAbsGam > 1e-10 ? 1.0 / sumAbsGam : 1.0;
							sumVtx.r() = (vtxI.r() * fabs(vtxI.g()) + vtxK.r() * fabs(vtxK.g())) * iws;
							break;

						case 1:
							sumVtx.r() = (fabs(vtxI.g()) > fabs(vtxK.g())) ? vtxI.r() : vtxK.r();
							break;
						}

						bool fl_hit = true;
						size_t hitpan = -1;

						for (size_t afl = 0; afl < W.getNumberOfAirfoil(); ++afl)
						{
							//проверим, не оказался ли новый вихрь внутри контура
							if (MoveInside(sumVtx.r(), vtxI.r(), W.getAirfoil(afl), hitpan) || MoveInside(sumVtx.r(), vtxK.r(), W.getAirfoil(afl), hitpan))
							{
								//std::cout << "HIT" << std::endl;
								fl_hit = false;
							}
						}//for

						if (fl_hit)
						{
							vtxI = sumVtx;
							vtxK.g() = 0.0;
							nHlop++;
							flag[vt] = true;
							flag[ssd] = true;
						}//if (fl_hit)

					}//if ((ssd!=0)&&(!flag[ssd]))

					if (flag[vt])
						break;
				}// for s
			}//if !flag[vt] 
		}//for vt

	}//for z

	return nHlop;
}//CollapsNew(...)

// Коллапс вихрей
int Wake::CollapsNewFast(int type, int times, std::vector<Vortex2D>& ri, std::vector<Vortex2D>& rj, std::vector<Point2D>& rnew, std::vector<std::pair<int, int>>& rindex)
{
	int nHlop = 0; //общее число убитых вихрей

	const double& cSP = collapseScaleParameter;
	const double& cRBP = collapseRightBorderParameter;

	double maxG = W.getPassport().wakeDiscretizationProperties.maxGamma;

	//int loc_hlop = 0; // neigh

	std::vector<bool> flag;	//как только вихрь сколлапсировался  flag = 1

	for (int z = 0; z < times; ++z)
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

				for (int s = 0; s < knbForRestruct; ++s)
				{
					//int ssd = neighb[vt];
					int ssd = neighbNew[vt * knbForRestruct + s];
					if (ssd == 0)
						continue;

					Vortex2D& vtxK = vtx[ssd];

					if (/*(ssd != 0) &&*/ (!flag[ssd])) //ssd==0 always true
					{


						Vortex2D sumVtx;
						sumVtx.g() = vtxI.g() + vtxK.g();

						switch (type)
						{
						case 0:
						case 2:
							sumAbsGam = fabs(vtxI.g()) + fabs(vtxK.g());
							iws = sumAbsGam > 1e-10 ? 1.0 / sumAbsGam : 1.0;
							sumVtx.r() = (vtxI.r() * fabs(vtxI.g()) + vtxK.r() * fabs(vtxK.g())) * iws;
							break;

						case 1:
							sumVtx.r() = (fabs(vtxI.g()) > fabs(vtxK.g())) ? vtxI.r() : vtxK.r();
							break;
						}

						bool fl_hit = true;
						size_t hitpan = -1;

						//for (size_t afl = 0; afl < W.getNumberOfAirfoil(); ++afl) //кол-во профилей
						//{
							ri.push_back(Vortex2D(vtxI.r(), vtxI.g()));
							rj.push_back(Vortex2D(vtxK.r(), vtxK.g()));
							rindex.push_back({ (int)vt, ssd });
							rnew.push_back(sumVtx.r());							
							nHlop++;
							flag[vt] = true;
							flag[ssd] = true;
							//проверим, не оказался ли новый вихрь внутри контура
							//if (MoveInside(sumVtx.r(), vtxI.r(), W.getAirfoil(afl), hitpan) || MoveInside(sumVtx.r(), vtxK.r(), W.getAirfoil(afl), hitpan))
							//{
							//	//std::cout << "HIT" << std::endl;
							//	fl_hit = false;
							//}
						//}//for

						//if (fl_hit)
						//{
						//	vtxI = sumVtx;
						//	vtxK.g() = 0.0;
						//	//nHlop++;
						//	flag[vt] = true;
						//	flag[ssd] = true;
						//}//if (fl_hit)

					}//if ((ssd!=0)&&(!flag[ssd]))

					if (flag[vt])
						break;
				}// for s
			}//if !flag[vt] 
		}//for vt

	}//for z

	return nHlop;
}//CollapsNew(...)



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

	return delta;
}//RemoveZero()




//Реструктуризация вихревого следа
void Wake::Restruct()
{
	W.getTimers().start("Restr");	

	if (W.getPassport().wakeDiscretizationProperties.epscol > 0)
	{
		double timePreKnn = -omp_get_wtime();

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
		
		timePreKnn += omp_get_wtime();

		//std::cout << "timePreKnn = " << timePreKnn * 1000 << " ms" << std::endl;
		

		/////////////////////////////////////////////////////////////

		
		if (W.getPassport().numericalSchemes.velocityComputation.second == 0)
		{
			//CPU/GPU - прямой алгоритм
			//double timeCollaps = -omp_get_wtime();
			Collaps(0, 1);
			//Collaps(1, 1);
			//Collaps(2, 1);
			//timeCollaps += omp_get_wtime();
			//std::cout << "timeCollaps = " << timeCollaps * 1000 << " ms" << std::endl;
		}
		else
		{
			//быстрый алгоритм
			for (int collapsStep = 0; collapsStep < 1; ++collapsStep)
			{
				//ttB = -omp_get_wtime();

				double timeResize = -omp_get_wtime();
				neighbNew.resize(vtx.size() * knbForRestruct);
				timeResize += omp_get_wtime();

				//std::cout << "timeResize = " << timeResize * 1000 << " ms" << std::endl;

				const double& cSP = collapseScaleParameter;
				const double& cRBP = collapseRightBorderParameter;
				const double& maxG = W.getPassport().wakeDiscretizationProperties.maxGamma;
				const double& epsCol = W.getPassport().wakeDiscretizationProperties.epscol;

#ifndef USE_CUDA
				//CPU
				std::vector<std::vector<std::pair<double, size_t>>> initdist(vtx.size());
				for (auto& d : initdist)
					d.resize(2 * knbForRestruct, { -1.0, -1 });
				double timeKnn = -omp_get_wtime();
				WakekNNnewForCollaps(vtx, knbForRestruct, initdist, cSP, cRBP, maxG, epsCol, collapsStep);//CPU
				timeKnn += omp_get_wtime();
				//std::cout << "Time_Knn_CPU = " << timeKnn * 1000 << " ms" << std::endl;
#else			
				double timeAlloc = -omp_get_wtime();
				std::vector<std::pair<double, size_t>> initdistcuda(knbForRestruct * vtx.size());  //CUDA
				timeAlloc += omp_get_wtime();
				//std::cout << "timeAlloc = " << timeAlloc * 1000 << " ms" << std::endl;

				double timeKnn = -omp_get_wtime();
				W.getNonConstCuda().RefreshWake(5);

				W.timerMerging.reset();
				W.timerMerging.start();
				//Построение дерева для вихрей
				BHcu::CudaTreeInfo knnTree(W.getCuda().blocks, tree_T::contr, object_T::point3, scheme_T::noScheme, false);
				knnTree.MemoryAllocate((int)W.getCuda().n_CUDA_wake);
				knnTree.Update((int)W.getWake().vtx.size(), W.getWake().devVtxPtr, W.getPassport().wakeDiscretizationProperties.eps);
				knnTree.Build();

				Point2D minr, maxr;

				cudaMemcpy(&minr, knnTree.minrD, sizeof(Point2D), cudaMemcpyDeviceToHost);
				cudaMemcpy(&maxr, knnTree.maxrD, sizeof(Point2D), cudaMemcpyDeviceToHost);


				kNNcuda<knbForRestruct>(minr, maxr, W.getCuda().blocks, vtx, initdistcuda, vecForKnn, cSP, cRBP, maxG, epsCol, collapsStep);                             //CUDA
				W.timerMerging.stop();
				timeKnn += omp_get_wtime();
				//std::cout << "Time_Knn_GPU = " << timeKnn * 1000 << " ms" << std::endl;
#endif			

				double timeCopy = -omp_get_wtime();
#pragma omp parallel for
				for (int i = 0; i < vtx.size(); ++i)
				{
#ifndef USE_CUDA
					for (int j = 0; j < knbForRestruct; ++j)
						neighbNew[i * knbForRestruct + j] = (int)initdist[i][j].second;
#else
					for (int j = 0; j < knbForRestruct; ++j)
						neighbNew[i * knbForRestruct + j] = (int)initdistcuda[i * knbForRestruct + j].second;
#endif
				}
				timeCopy += omp_get_wtime();
				//std::cout << "timeCopy = " << timeCopy * 1000 << " ms" << std::endl;

				//
				//std::ofstream initdistFile(W.getPassport().dir + "initdist-GPU" + std::to_string(W.currentStep));
				//for (int i = 0; i < vtx.size(); ++i)
				//{
				//	initdistFile << i;
				//	for (int k = 0; k < knb; ++k)
				//		initdistFile << " " << neighbNew[i * (knb)+k] << " " << (vtx[i].r() - vtx[neighbNew[i * (knb)+k]].r()).length();
				//	initdistFile << std::endl;
				//	
				//	//for (int k = 0; k < knb; ++k)
				//		//initdistFile << " " << neighbNew[i * (knb)+k] << " " << vtx[neighbNew[i * (knb)+k]].g() << " " << (vtx[neighbNew[i * (knb)+k]].r() - vtx[i].r()).length() << "; ";					
				//	//initdistFile << std::endl;
				//}
				//initdistFile.close();
				//


				double timeReserve = -omp_get_wtime();
				std::vector<Vortex2D> ri, rj;
				std::vector<Point2D> rnew;
				std::vector<std::pair<int, int>> rindex;
				ri.reserve(vtx.size());
				rj.reserve(vtx.size());
				rnew.reserve(vtx.size());

				rindex.reserve(vtx.size());
				timeReserve += omp_get_wtime();
				//std::cout << "timeReserve = " << timeReserve * 1000 << " ms" << std::endl;


#ifdef	USE_CUDA
			//поиск пар для объединения			
				int nHlop = CollapsNewFast(0, 1, ri, rj, rnew, rindex); //тут заполняются ri, rj, rnew, rindex

				//формирование траекторий объединямых вихрей
				std::vector<std::pair<Point2D, Point2D>> segments(2 * ri.size());
				for (size_t i = 0; i < ri.size(); ++i)
				{
					segments[2 * i + 0].first = ri[i];
					segments[2 * i + 0].second = rnew[i];

					segments[2 * i + 1].first = rj[i];
					segments[2 * i + 1].second = rnew[i];
				}

				//контроль протыкания
				std::vector<int> hit(segments.size());
				if (ri.size() > 0)
				{
					double* devSegments_ptr;

					cudaMalloc(&devSegments_ptr, segments.size() * sizeof(double) * 4);
					cudaMemcpy(devSegments_ptr, segments.data(), segments.size() * sizeof(double) * 4, cudaMemcpyHostToDevice);

					W.timerInside.start();
					auto& cntrTreeSeg = *W.getCuda().cntrTreeSegment;
					cntrTreeSeg.MemoryAllocate((int)W.getCuda().n_CUDA_wake);
					cntrTreeSeg.UpdatePanelGeometry((int)segments.size(), (double4*)devSegments_ptr);
					cntrTreeSeg.Build();

					BHcu::treePanelsSegmentsIntersectionCalculationWrapper(*W.getNonConstCuda().auxTreePnl, cntrTreeSeg,
						W.getWake().devNearestPanelPtr);

					W.timerInside.stop();

					cudaMemcpy(hit.data(), W.getWake().devNearestPanelPtr, segments.size() * sizeof(int), cudaMemcpyDeviceToHost);
					cudaFree(devSegments_ptr);
				}//if ri.size() > 0

				for (size_t i = 0; i < ri.size(); ++i)
				{
					if (hit[2 * i + 0] == -1 && hit[2 * i + 1] == -1)
					{
						vtx[rindex[i].first].r() = rnew[i];
						vtx[rindex[i].first].g() += vtx[rindex[i].second].g();

						vtx[rindex[i].second].g() = 0.0;
					}
				}

#else
				double timeCollaps;
				//CollapsNew(0, 1);
				timeCollaps = -omp_get_wtime();
				int nHlop = CollapsNewFast(0, 1, ri, rj, rnew, rindex);

				timeCollaps += omp_get_wtime();

				size_t hitA, hitB;
				bool ifInside;
#pragma omp parallel for private (hitA, hitB, ifInside)
				for (int i = 0; i < (int)ri.size(); ++i)
				{
					ifInside = false;
					for (size_t q = 0; q < W.getNumberOfAirfoil(); ++q)
					{
						MoveInside(rnew[i], ri[i], W.getAirfoil(q), hitA);
						MoveInside(rnew[i], rj[i], W.getAirfoil(q), hitB);
						if ((hitA != size_t(-1)) || (hitB != size_t(-1)))
							ifInside = true;
					}
					if (!ifInside)
					{
						vtx[rindex[i].first].r() = rnew[i];
						vtx[rindex[i].first].g() += vtx[rindex[i].second].g();

						vtx[rindex[i].second].g() = 0.0;

					}
				}
#endif
			}
		}
	
	}
	double timeRemove = -omp_get_wtime();
	RemoveFar();
	RemoveZero();
	timeRemove += omp_get_wtime();
	//std::cout << "timeRemove = " << timeRemove * 1000 << std::endl;

	W.getTimers().stop("Restr");
//*
//	std::cout << "Pre-knn time = " << timePreKnn * 1000.0 << " ms" << std::endl;
//	std::cout << "Pre-knn resize time = " << timeResize * 1000.0 << " ms" << std::endl;
//	std::cout << "Time_Knn_GPU = " << timeKnn * 1000 << " ms" << std::endl;
//	std::cout << "Post-knn copy time = " << timeCopy * 1000.0 << " ms" << std::endl;
//	std::cout << "Reserve_time = " << timeReserve * 1000 << " ms" << std::endl;
//	std::cout << "Collaps_time = " << timeCollaps * 1000 << " ms" << std::endl;
//	std::cout << "Copy2_time = " << timeCopy2 * 1000 << " ms" << std::endl;
//	std::cout << "Ray_time = " << timeRay * 1000 << " ms" << std::endl;
//	std::cout << "Collaps2_time = " << timeCollaps2 * 1000 << " ms" << std::endl;
//	std::cout << "Remove_time = " << timeRemove * 1000 << " ms" << std::endl;
	

//*/
}//Restruct()