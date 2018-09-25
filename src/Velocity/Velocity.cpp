/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Velocity.cpp                                                     |
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
\brief Файл кода с описанием класса Velocity
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/ 


#include "Velocity.h"

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
	cuda.ExpCalcConvVeloToSetOfPoints(
		wake.vtx.size(), 
		wake.devWakePtr, 
		wake.vtx.size(), 
		wake.devWakePtr, 
		boundary.size(), 
		cuda.dev_nPanels, 
		cuda.dev_ptr_ptr_pnl, 
		wakeVortexesParams.convVelo, 
		wakeVortexesParams.epsastWake, 
		cuda.vels, 
		cuda.rads, 
		cuda.dev_ptr_vel, 
		cuda.dev_ptr_rad, 
		2.0*wake.param.epscol, wake.param.eps2);
#else
	CalcConvVeloToSetOfPoints(wake.vtx, wakeVortexesParams.convVelo, wakeVortexesParams.epsastWake);
#endif

	/*
	if (parallel.myidWork == 0)
	{
		std::ostringstream ss;
		ss << "bouVeloNew";
		std::ofstream bouVeloNewFile(ss.str());
		bouVeloNewFile << wakeVortexesParams.convVelo << std::endl;
		bouVeloNewFile.close();
	}
	*/
	
	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_CONV_TOBOU))	
		cuda.ExpCalcConvVeloToSetOfPoints(boundary[bou]->virtualWake.size(), cuda.host_ptr_ptr_pnl[bou], wake.vtx.size(), wake.devWakePtr, boundary.size(), cuda.dev_nPanels, cuda.dev_ptr_ptr_pnl, virtualVortexesParams[bou].convVelo, virtualVortexesParams[bou].epsastWake, cuda.virtvels[bou], cuda.virtrads[bou], cuda.host_ptr_ptr_vel[bou], cuda.host_ptr_ptr_rad[bou], 2.0*wake.param.epscol, wake.param.eps2);
#else
		CalcConvVeloToSetOfPoints(boundary[bou]->virtualWake, virtualVortexesParams[bou].convVelo, virtualVortexesParams[bou].epsastWake);
#endif	
		/*std::ostringstream sss;
		sss << "virtVelo_";
		std::ofstream virtVeloFile(sss.str());
		for (int i = 0; i < convVirtualVelo[0].size(); ++i)
			virtVeloFile << convVirtualVelo[0][i] << std::endl;
		virtVeloFile.close();*/
	}

}//CalcConvVelo()


//Вычисление диффузионных скоростей вихрей и виртуальных вихрей в вихревом следе
void Velocity::CalcDiffVelo()
{	
	double tTemp[4];

	tTemp[0] = omp_get_wtime();

	/// \todo Сделать влияние присоединенных слоев

	//пелена на пелену
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))	
	cuda.ExpCalcDiffVeloI1I2ToSetOfPoints(wake.vtx.size(), wake.devWakePtr, cuda.dev_ptr_rad, wake.vtx.size(), wake.devWakePtr, wakeVortexesParams.I1, wakeVortexesParams.I2, cuda.i1, cuda.i2, cuda.dev_ptr_i1, cuda.dev_ptr_i2, true);
#else
	CalcDiffVeloI1I2ToSetOfPoints(wake.vtx, wakeVortexesParams.epsastWake, wake.vtx, wakeVortexesParams.I1, wakeVortexesParams.I2);
#endif

	tTemp[1] = omp_get_wtime();

	/// \todo Не заменить ли ее на VirtualWakeSynchronize?
	for (size_t bou = 0; bou < boundary.size(); ++bou)
	{
		int lenVirtWake = (int)(boundary[bou]->virtualWake.size());
		MPI_Bcast(&lenVirtWake, 1, MPI_INT, 0, parallel.commWork);
		if (parallel.myidWork != 0)
			boundary[bou]->virtualWake.resize(lenVirtWake);
		MPI_Bcast(boundary[bou]->virtualWake.data(), lenVirtWake, Vortex2D::mpiVortex2D, 0, parallel.commWork);

		//виртуальные на границе на след
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))		
		cuda.ExpCalcDiffVeloI1I2ToSetOfPoints(wake.vtx.size(), wake.devWakePtr, cuda.dev_ptr_rad, boundary[bou]->virtualWake.size(), cuda.host_ptr_ptr_pnl[bou], wakeVortexesParams.I1, wakeVortexesParams.I2, cuda.i1, cuda.i2, cuda.dev_ptr_i1, cuda.dev_ptr_i2, false);
#else		
		CalcDiffVeloI1I2ToSetOfPoints(wake.vtx, wakeVortexesParams.epsastWake, boundary[bou]->virtualWake, wakeVortexesParams.I1, wakeVortexesParams.I2);
#endif

		// след на виртуальные		
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))			
		cuda.ExpCalcDiffVeloI1I2ToSetOfPoints(boundary[bou]->virtualWake.size(), cuda.host_ptr_ptr_pnl[bou], cuda.host_ptr_ptr_rad[bou], wake.vtx.size(), wake.devWakePtr, virtualVortexesParams[bou].I1, virtualVortexesParams[bou].I2, cuda.virti1[bou], cuda.virti2[bou], cuda.host_ptr_ptr_i1[bou], cuda.host_ptr_ptr_i2[bou], false);
#else
		CalcDiffVeloI1I2ToSetOfPoints(boundary[bou]->virtualWake, virtualVortexesParams[bou].epsastWake, wake.vtx, virtualVortexesParams[bou].I1, virtualVortexesParams[bou].I2);
#endif

		for (size_t targetBou = 0; targetBou < boundary.size(); ++targetBou)
		{
		// виртуальные на виртуальные
#if (defined(__CUDACC__) || defined(USE_CUDA)) && (defined(CU_I1I2))	
			cuda.ExpCalcDiffVeloI1I2ToSetOfPoints(boundary[targetBou]->virtualWake.size(), cuda.host_ptr_ptr_pnl[targetBou], cuda.host_ptr_ptr_rad[targetBou], boundary[bou]->virtualWake.size(), cuda.host_ptr_ptr_pnl[bou], virtualVortexesParams[targetBou].I1, virtualVortexesParams[targetBou].I2, cuda.virti1[targetBou], cuda.virti2[targetBou], cuda.host_ptr_ptr_i1[targetBou], cuda.host_ptr_ptr_i2[targetBou], false);
#else			
			CalcDiffVeloI1I2ToSetOfPoints(boundary[targetBou]->virtualWake, virtualVortexesParams[targetBou].epsastWake, boundary[bou]->virtualWake, virtualVortexesParams[targetBou].I1, virtualVortexesParams[targetBou].I2);
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
		
	for (size_t targetBou = 0; targetBou < boundary.size(); ++targetBou)
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

	//std::cout << "Times_diff:" << \
		tTemp[1] - tTemp[0] << " " << \
		tTemp[2] - tTemp[1] << " " << \
		tTemp[3] - tTemp[2] << " " << \
		std::endl;
	
}//CalcDiffVelo()
