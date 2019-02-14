/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.5    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/02/20     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Velocity2D.cpp                                                   |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VM is distributed in the hope that it will be useful, but WITHOUT           |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Файл кода с описанием класса Velocity
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.5   
\date 20 февраля 2019 г.
*/ 

#include "Velocity2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "StreamParser.h"
#include "Wake2D.h"
#include "WakeDataBase2D.h"
#include "World2D.h"

using namespace VM2D;

//Вычисление конвективных скоростей вихрей и виртуальных вихрей в вихревом следе
void Velocity::CalcConvVelo()
{
	/*
	if (parallel.myidWork == 0)
	{
		std::ostringstream ss1;
		ss1 << "wakeFile";
		std::ofstream wakefile(ss1.str());
		for (int i = 0; i < wake.vtx.size(); i++)
			wakefile << wake.vtx[i].r()[0] << " " << wake.vtx[i].r()[1] << " " << wake.vtx[i].g() << std::endl;
		wakefile.close();
	}
	*/
	
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONV_TOWAKE))
	GPUCalcConvVeloToSetOfPoints(W.getWake(), wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake);
#else
	CalcConvVeloToSetOfPoints(W.getWake(), wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake);
#endif

	std::vector<Point2D> nullVector(0);

	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONV_TOBOU))	
		GPUCalcConvVeloToSetOfPoints(W.getBoundary(bou).virtualWake, nullVector, virtualVortexesParams[bou].epsastWake, true);
#else
		CalcConvVeloToSetOfPoints(W.getBoundary(bou).virtualWake, nullVector, virtualVortexesParams[bou].epsastWake, true);
#endif	
	}
}//CalcConvVelo()


//Вычисление диффузионных скоростей вихрей и виртуальных вихрей в вихревом следе
void Velocity::CalcDiffVelo()
{	
	double tTemp[4];

	tTemp[0] = omp_get_wtime();

	//пелена на пелену
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))	
	GPUCalcDiffVeloI1I2ToSetOfPoints(W.getWake(), wakeVortexesParams.epsastWake, W.getWake(), wakeVortexesParams.I1, wakeVortexesParams.I2);
#else
	CalcDiffVeloI1I2ToSetOfPoints(W.getWake(), wakeVortexesParams.epsastWake, W.getWake(), wakeVortexesParams.I1, wakeVortexesParams.I2);
#endif


	tTemp[1] = omp_get_wtime();
		
	for (size_t bou = 0; bou < W.getNumberOfBoundary(); ++bou)
	{
		//не нужно, т.к. сделано выше перед началом вычисления скоростей
		//W.getNonConstBoundary(bou).virtualWake.WakeSynchronize();
		
		//виртуальные на границе на след
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))				
		GPUCalcDiffVeloI1I2ToSetOfPointsFromPanels(W.getWake(), wakeVortexesParams.epsastWake, W.getBoundary(bou), wakeVortexesParams.I1, wakeVortexesParams.I2);
#else			
		CalcDiffVeloI1I2ToSetOfPointsFromPanels(W.getWake(), wakeVortexesParams.epsastWake, W.getBoundary(bou), wakeVortexesParams.I1, wakeVortexesParams.I2);
#endif
		
		// след на виртуальные		
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))					
		GPUCalcDiffVeloI1I2ToSetOfPoints(W.getBoundary(bou).virtualWake, virtualVortexesParams[bou].epsastWake, W.getWake(), virtualVortexesParams[bou].I1, virtualVortexesParams[bou].I2);
#else
		CalcDiffVeloI1I2ToSetOfPoints(W.getBoundary(bou).virtualWake, virtualVortexesParams[bou].epsastWake, W.getWake(), virtualVortexesParams[bou].I1, virtualVortexesParams[bou].I2);
#endif

		for (size_t targetBou = 0; targetBou < W.getNumberOfBoundary(); ++targetBou)
		{
		// виртуальные на виртуальные
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))				
			GPUCalcDiffVeloI1I2ToSetOfPointsFromPanels(W.getBoundary(targetBou).virtualWake, virtualVortexesParams[targetBou].epsastWake, W.getBoundary(bou), virtualVortexesParams[targetBou].I1, virtualVortexesParams[targetBou].I2);
#else			
			CalcDiffVeloI1I2ToSetOfPointsFromPanels(W.getBoundary(targetBou).virtualWake, virtualVortexesParams[targetBou].epsastWake, W.getBoundary(bou), virtualVortexesParams[targetBou].I1, virtualVortexesParams[targetBou].I2);
#endif	
		}
	} //for bou

tTemp[2] = omp_get_wtime();

	Point2D I2;
	double I1;

	for (size_t vt = 0; vt < wakeVortexesParams.diffVelo.size(); ++vt)
	{
		I2 = wakeVortexesParams.I2[vt];
		I1 = wakeVortexesParams.I1[vt];
		if (fabs(I1) < 1.e-8)
			wakeVortexesParams.diffVelo[vt] = { 0.0, 0.0 };
		else
			wakeVortexesParams.diffVelo[vt] = I2 * (1.0 / (I1 * wakeVortexesParams.epsastWake[vt]));
	}
		
	for (size_t targetBou = 0; targetBou < W.getNumberOfBoundary(); ++targetBou)
	for (size_t vt = 0; vt < virtualVortexesParams[targetBou].diffVelo.size(); ++vt)
	{
		I2 = virtualVortexesParams[targetBou].I2[vt];
		I1 = virtualVortexesParams[targetBou].I1[vt];

		if (fabs(I1) < 1.e-8)
			virtualVortexesParams[targetBou].diffVelo[vt] = { 0.0, 0.0 };
		else
			virtualVortexesParams[targetBou].diffVelo[vt] = I2 * (1.0 / (I1 * virtualVortexesParams[targetBou].epsastWake[vt]));
	}

	tTemp[3] = omp_get_wtime();

#pragma warning (push)
#pragma warning (disable: 4010)
	//W.getInfo('t') << "Times_diff:" << \
		tTemp[1] - tTemp[0] << " " << \
		tTemp[2] - tTemp[1] << " " << \
		tTemp[3] - tTemp[2] << " " << \
		std::endl;
#pragma warning (pop)
	
}//CalcDiffVelo()
