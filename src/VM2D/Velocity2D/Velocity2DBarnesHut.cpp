/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Velocity2DBarnesHut.cpp                                          |
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
\brief Файл кода с описанием класса VelocityBarnesHut
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#include "Velocity2DBarnesHut.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "StreamParser.h"
#include "Wake2D.h"
#include "World2D.h"
#include "Gmres2D.h"

using namespace VM2D;

/// Конструктор
VelocityBarnesHut::VelocityBarnesHut(const World2D& W_) :
	Velocity(W_)
{
};

/// Деструктор
VelocityBarnesHut::~VelocityBarnesHut()
{	
};


double VelocityBarnesHut::CalcConvVeloToSetOfPointsFromWake(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius)
{
	int order = W.getPassport().numericalSchemes.nbodyMultipoleOrder;
	double theta = W.getPassport().numericalSchemes.nbodyTheta;

	VMlib::vmTimer timerA;
	timerA.start();
	W.getInflTreeWake().DownwardTraversalVorticesToPoints(W.getCntrTreeWake(), velo, domainRadius, theta, order, calcRadius);
	timerA.stop();

	//std::ofstream treeTimeFile;
	//if (W.getCurrentStep() == 0)
	//{
	//	treeTimeFile.open(W.getPassport().dir + "/dbg/treeTime.csv");
	//	treeTimeFile << "step,time,N,tBUI,tUPW,tDNW\n";
	//}
	//else
	//	treeTimeFile.open(W.getPassport().dir + "/dbg/treeTime.csv", std::ios::app);
	//
	//treeTimeFile << W.getCurrentStep() << ',' << W.getCurrentTime() << ',' << W.getInflTreeWake().object.size() << ',' << W.getCntrTreeWake().object.size() << ',' << timerA.duration() << '\n';
	//
	//treeTimeFile.close();
	
#ifdef USE_CUDA	
	if (calcVelo)
		W.getCuda().CopyMemToDev<double, 2>(velo.size(), (double*)velo.data(), pointsDb.devVelPtr);
	if (calcRadius)
		W.getCuda().CopyMemToDev<double, 1>(domainRadius.size(), domainRadius.data(), pointsDb.devRadPtr);
#endif

	return timerA.duration();
}


inline void ModifyE2(double* ee2, double dst2)
{
	if (dst2 > 0)
	{
		if (dst2 < ee2[0])
		{
			ee2[2] = ee2[1];
			ee2[1] = ee2[0];
			ee2[0] = dst2;
		}//if (dist2<ee2[0])
		else
		{
			if (dst2 < ee2[1])
			{
				ee2[2] = ee2[1];
				ee2[1] = dst2;
			}// if (dist2<ee2[1])
			else
				if (dst2 < ee2[2])
					ee2[2] = dst2;
		}//else
	}//if (dst2>0)
}





//Вычисление конвективных скоростей и радиусов вихревых доменов в заданном наборе точек от следа
void VelocityBarnesHut::CalcConvVPVeloToSetOfPointsFromWake(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius)
{
	std::vector<Point2D> selfVelo(pointsDb.vtx.size());
	domainRadius.resize(pointsDb.vtx.size());

	double cft = IDPI;

#pragma warning (push)
#pragma warning (disable: 4101)
	//Локальные переменные для цикла
	Point2D velI;
	Point2D tempVel;
	double dst2eps, dst2;
#pragma warning (pop)

	if (calcVelo)
	{
#pragma omp parallel for default(none) shared(selfVelo, cft, calcVelo, calcRadius, pointsDb, domainRadius) private(tempVel, velI, dst2, dst2eps) schedule(dynamic, DYN_SCHEDULE)
		for (int i = 0; i < pointsDb.vtx.size(); ++i)
		{
			double ee2[3] = { 10000.0, 10000.0, 10000.0 };

			velI.toZero();

			const Point2D& posI = pointsDb.vtx[i].r();

			for (size_t j = 0; j < W.getWake().vtx.size(); ++j)
			{
				const Point2D& posJ = W.getWake().vtx[j].r();

				dst2 = (posI - posJ).length2();

//				//Модифицируем массив квадратов расстояний до ближайших вихрей из wake
//#ifndef TESTONLYVELO
//				if (calcRadius)
//					VMlib::ModifyE2(ee2, dst2);
//#endif //!TESTONLYVELO				

				const double& gamJ = W.getWake().vtx[j].g();

				tempVel.toZero();
				dst2eps = VMlib::boundDenom(dst2, sqr(W.getWake().vtx[j].sigma())); //Сглаживать надо!!!				

				tempVel = { -posI[1] + posJ[1], posI[0] - posJ[0] };
				tempVel *= (gamJ / dst2eps);
				velI += tempVel;
			}

			for (size_t j = 0; j < W.getSource().vtx.size(); ++j)
			{
				const Point2D& posJ = W.getSource().vtx[j].r();
				const double& gamJ = W.getPassport().physicalProperties.accelCft(W.getCurrentTime()) * W.getSource().vtx[j].g();

				tempVel.toZero();

				dst2 = dist2(posI, posJ);
				dst2eps = VMlib::boundDenom(dst2, sqr(W.getSource().vtx[j].sigma())); //Сглаживать надо!!!

				tempVel = { posI[0] - posJ[0], posI[1] - posJ[1] };
				tempVel *= (gamJ / dst2eps);
				velI += tempVel;
			}
//#ifndef TESTONLYVELO
//			if (calcRadius)
//			{
//				for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
//				{
//					const auto& bou = W.getBoundary(s);
//					//Модифицируем массив квадратов расстояний до ближайших вихрей из virtualWake					
//					for (size_t j = 0; j < bou.virtualWake.vtx.size(); ++j)
//						ModifyE2(ee2, dist2(posI, bou.virtualWake.vtx[j].r()));
//				}
//			}
//#endif

			velI *= cft;
			selfVelo[i] = velI;

//#ifndef TESTONLYVELO
//			if (calcRadius)
//				domainRadius[i] = 1.0 * sqrt((ee2[0] + ee2[1] + ee2[2]) / 3.0);
//#endif
		}
	}
//	else if (calcRadius)
//	{
//#pragma omp parallel for default(none) shared(selfVelo, cft, calcVelo, calcRadius, pointsDb, domainRadius) private(tempVel, velI, dst2) schedule(dynamic, DYN_SCHEDULE)
//		for (int i = 0; i < pointsDb.vtx.size(); ++i)
//		{
//			double ee2[3] = { 10000.0, 10000.0, 10000.0 };
//
//			velI.toZero();
//
//			const Point2D& posI = pointsDb.vtx[i].r();
//
//			for (size_t j = 0; j < W.getWake().vtx.size(); ++j)
//			{
//				const Point2D& posJ = W.getWake().vtx[j].r();
//
//				dst2 = dist2(posI, posJ);
//
//				//Модифицируем массив квадратов расстояний до ближайших вихрей из wake
//				VMlib::ModifyE2(ee2, dst2);
//			}
//
//			for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
//			{
//				for (size_t j = 0; j < W.getBoundary(s).virtualWake.vtx.size(); ++j)
//				{
//					const Point2D& posJ = W.getBoundary(s).virtualWake.vtx[j].r();
//					dst2 = dist2(posI, posJ);
//
//					//Модифицируем массив квадратов расстояний до ближайших вихрей из virtualWake
//					ModifyE2(ee2, dst2);
//				}
//			}
//
//			domainRadius[i] = 1.0 * sqrt((ee2[0] + ee2[1] + ee2[2]) / 3.0);
//		}
//	} //else


	if (calcVelo)
		for (size_t i = 0; i < velo.size(); ++i)
			velo[i] += selfVelo[i];

#ifdef USE_CUDA	
	if (calcVelo)
		W.getCuda().CopyMemToDev<double, 2>(velo.size(), (double*)velo.data(), pointsDb.devVelPtr);
//	if (calcRadius)
//		W.getCuda().CopyMemToDev<double, 1>(domainRadius.size(), domainRadius.data(), pointsDb.devRadPtr);
#endif

}//CalcConvVPVeloToSetOfPointsFromWake(...)




#if defined(USE_CUDA)
double VelocityBarnesHut::GPUCalcConvVeloToSetOfPointsFromWake(std::unique_ptr<BHcu::CudaTreeInfo>& cntrTree, const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius)
{	
	const double& theta = W.getPassport().numericalSchemes.nbodyTheta;
	const int& order = W.getPassport().numericalSchemes.nbodyMultipoleOrder;
	
	int npt = cntrTree->nObject;

	float tBLD = 0.0f, tUPW = 0.0f, tDNV;

	auto& inflTree = *W.getCuda().inflTreeWake;		


	VMlib::vmTimer timerA;
	timerA.start();
	if (cntrTree->nObject > 1024 * 80)
		tDNV = inflTree.DownwardTraversalVorticesToPoints<1024>(*cntrTree, (Point2D*)pointsDb.devVelPtr, pointsDb.devRadPtr, theta, order, calcRadius);
	else
		tDNV = inflTree.DownwardTraversalVorticesToPoints<32>(*cntrTree, (Point2D*)pointsDb.devVelPtr, pointsDb.devRadPtr, theta, order, calcRadius);
	timerA.stop();

	//std::cout << "nvt = " << npt << std::endl;
	//std::cout << "tBLD = " << tBLD << std::endl;
	//std::cout << "tUPW = " << tUPW << std::endl;
	//std::cout << "tDNV = " << tDNV << std::endl;

	std::vector<Point2D> Vel;
	if (calcVelo)
	{
		Vel.resize(npt);
		W.getCuda().CopyMemFromDev<double, 2>(npt, (double*)pointsDb.devVelPtr, (double*)Vel.data(), 20);

		for (size_t q = 0; q < npt; ++q)
			velo[q] += Vel[q];
	}

	if (calcRadius)
		W.getCuda().CopyMemFromDev<double, 1>(npt, pointsDb.devRadPtr, domainRadius.data(), 212);

	return timerA.duration();
}


void VelocityBarnesHut::GPUCalcConvVelocityToSetOfPointsFromSheets(std::unique_ptr<BHcu::CudaTreeInfo>& cntrTree, const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const
{
	double*& velD = pointsDb.devVelPtr;

	const double& theta = W.getPassport().numericalSchemes.nbodyTheta;
	const int& order = W.getPassport().numericalSchemes.nbodyMultipoleOrder;

	std::vector<Point2D> newV(velo.size());

	int npt = cntrTree->nObject;

	if (npt > 0)
	{	
		int nTotPan = (int)W.getAirfoil(0).getNumberOfPanels();
		for (size_t q = 1; q < W.getNumberOfAirfoil(); ++q)
			nTotPan += (int)W.getAirfoil(q).getNumberOfPanels();
		auto& afl = W.getAirfoil(0);

		//Влияние вихревых слоев
		W.getCuda().inflTreePnlVortex->DownwardTraversalPanelsToPoints(*cntrTree, (Point2D*)velD, theta, order);
		W.getCuda().CopyMemFromDev<double, 2>(npt, velD, (double*)newV.data());

		for (size_t q = 0; q < velo.size(); ++q)
			velo[q] += newV[q];

		if (W.isAnyMovableOrDeformable())
		{
			//Влияние слоев источников
			W.getCuda().inflTreePnlSource->DownwardTraversalPanelsToPoints(*cntrTree, (Point2D*)velD, theta, order);
			W.getCuda().CopyMemFromDev<double, 2>(npt, velD, (double*)newV.data());

			for (size_t q = 0; q < velo.size(); ++q)
				velo[q] += newV[q];
		}
	}
}//GPUCalcConvVelocityToSetOfPointsFromSheets(...)
#endif





 



