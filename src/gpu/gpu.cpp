/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.3    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/09/26     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: gpu.cpp                                                          |
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
\brief Файл кода с описанием класса gpu
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/

#include "gpu.h"
#include "World2D.h"

gpu::gpu(const World2D& W_)
	: W(W_)
{
#if defined(__CUDACC__) || defined(USE_CUDA)
	
	
	cuSetConstants(sizeof(Vortex2D)/sizeof(double), Vortex2D::offsPos / sizeof(double), Vortex2D::offsGam / sizeof(double) );
	
	
	W.getWake().devWakePtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(INC_VORT_DEV, W.getWake().devNWake);
	
	W.getWake().devVelsPtr = ReserveDevMem<double, 2>(INC_VORT_DEV, n_CUDA_vel);
	W.getWake().devI0Ptr = ReserveDevMem<double, 1>(INC_VORT_DEV, n_CUDA_i0);
	W.getWake().devI1Ptr = ReserveDevMem<double, 1>(INC_VORT_DEV, n_CUDA_i1);
	W.getWake().devI2Ptr = ReserveDevMem<double, 2>(INC_VORT_DEV, n_CUDA_i2);
	W.getWake().devI3Ptr = ReserveDevMem<double, 2>(INC_VORT_DEV, n_CUDA_i3);
	W.getWake().devRadsPtr = ReserveDevMem<double, 1>(INC_VORT_DEV, n_CUDA_rad);

	W.getWake().devMeshPtr = ReserveDevMem<int, 2>(INC_VORT_DEV, n_CUDA_mesh);
	W.getWake().devNeiPtr = ReserveDevMem<int, 1>(INC_VORT_DEV, n_CUDA_nei);

	W.getWake().tmpVels.resize(n_CUDA_vel, { 0.0, 0.0 });
	W.getWake().tmpRads.resize(n_CUDA_vel, 1.0e+5);

	W.getWake().tmpI0.resize(n_CUDA_vel, 0.0);
	W.getWake().tmpI1.resize(n_CUDA_vel, 0.0);
	W.getWake().tmpI2.resize(n_CUDA_vel, { 0.0, 0.0 });
	W.getWake().tmpI3.resize(n_CUDA_vel, { 0.0, 0.0 });
	W.getWake().tmpNei.resize(n_CUDA_vel, 0);
	//mesh.resize(n_CUDA_vel, { 0, 0 });
		

	W.getWake().devNWake = 0;
	n_CUDA_afls = 0;
#endif
}


gpu::~gpu()
{
#if defined(__CUDACC__) || defined(USE_CUDA)
	cuDeleteFromDev(W.getWake().devWakePtr);
	cuDeleteFromDev(W.getWake().devVelsPtr);
	cuDeleteFromDev(W.getWake().devRadsPtr);
	cuDeleteFromDev(W.getWake().devI0Ptr);
	cuDeleteFromDev(W.getWake().devI1Ptr);
	cuDeleteFromDev(W.getWake().devI2Ptr);
	cuDeleteFromDev(W.getWake().devI3Ptr);
	cuDeleteFromDev(W.getWake().devMeshPtr);
	cuDeleteFromDev(W.getWake().devNeiPtr);

	for (size_t s = 0; s < n_CUDA_afls; ++s)
	{
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devWakePtr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devI0Ptr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devI1Ptr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devI2Ptr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devI3Ptr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devRadsPtr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devVelsPtr);
		
		cuDeleteFromDev(W.getBoundary(s).afl.devR);
		cuDeleteFromDev(W.getBoundary(s).afl.devRhs);
	}

	if (n_CUDA_afls)
	{
		cuDeleteFromDev(dev_nPanels);
		cuDeleteFromDev(dev_nVortices);

		cuDeleteFromDev(dev_ptr_ptr_vtx);		
		cuDeleteFromDev(dev_ptr_ptr_i0);
		cuDeleteFromDev(dev_ptr_ptr_i1);
		cuDeleteFromDev(dev_ptr_ptr_i2);
		cuDeleteFromDev(dev_ptr_ptr_i3);
		cuDeleteFromDev(dev_ptr_ptr_rad);
		cuDeleteFromDev(dev_ptr_ptr_vel);
		
		cuDeleteFromDev(dev_ptr_ptr_r);
		cuDeleteFromDev(dev_ptr_ptr_rhs);
	}
#endif
}

#if defined(__CUDACC__) || defined(USE_CUDA)

void gpu::RefreshWake() 
{
	const int& id = W.getParallel().myidWork;

	if (W.getWake().vtx.size() > 0)
	{
		//if (id == 0)
		{
			//size_t cuSize = cuGetCuSize();
			//size_t cpuSize = sizeof(Vortex2D);

			while (W.getWake().vtx.size() > W.getWake().devNWake)
			{
				size_t curLength = W.getWake().devNWake;

				cuDeleteFromDev(W.getWake().devWakePtr);
				cuDeleteFromDev(W.getWake().devVelsPtr);
				cuDeleteFromDev(W.getWake().devRadsPtr);
				cuDeleteFromDev(W.getWake().devI0Ptr);
				cuDeleteFromDev(W.getWake().devI1Ptr);
				cuDeleteFromDev(W.getWake().devI2Ptr);
				cuDeleteFromDev(W.getWake().devI3Ptr);
				cuDeleteFromDev(W.getWake().devMeshPtr);
				cuDeleteFromDev(W.getWake().devNeiPtr);

				size_t sz = curLength + INC_VORT_DEV;

				W.getWake().devWakePtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(sz, W.getWake().devNWake);
				W.getWake().devVelsPtr = ReserveDevMem<double, 2>(sz, n_CUDA_vel);
				W.getWake().devRadsPtr = ReserveDevMem<double, 1>(sz, n_CUDA_rad);

				W.getWake().devI0Ptr = ReserveDevMem<double, 1>(sz, n_CUDA_i0);
				W.getWake().devI1Ptr = ReserveDevMem<double, 1>(sz, n_CUDA_i1);
				W.getWake().devI2Ptr = ReserveDevMem<double, 2>(sz, n_CUDA_i2);
				W.getWake().devI3Ptr = ReserveDevMem<double, 2>(sz, n_CUDA_i3);

				W.getWake().devMeshPtr = ReserveDevMem<int, 2>(sz, n_CUDA_mesh);
				W.getWake().devNeiPtr = ReserveDevMem<int, 1>(sz, n_CUDA_nei);

				W.getWake().tmpVels.resize(n_CUDA_vel, { 0.0, 0.0 });
				W.getWake().tmpRads.resize(n_CUDA_vel, 1.0e+5);

				W.getWake().tmpI0.resize(n_CUDA_vel, 0.0);
				W.getWake().tmpI1.resize(n_CUDA_vel, 0.0);
				W.getWake().tmpI2.resize(n_CUDA_vel, { 0.0, 0.0 });
				W.getWake().tmpI3.resize(n_CUDA_vel, { 0.0, 0.0 });
				W.getWake().tmpNei.resize(n_CUDA_vel, 0);

				//mesh.resize(n_CUDA_vel, { 0, 0 });

				std::cout << "resize: " << n_CUDA_vel << std::endl;
			}

			cuClearWakeMem(W.getWake().devNWake, W.getWake().devWakePtr);
			
			cuCopyWakeToDev(W.getWake().vtx.size(), W.getWake().vtx.data(), W.getWake().devWakePtr);			
		}
	}	
}


void gpu::RefreshAfls() 
{
	const int& id = W.getParallel().myidWork;
	if (W.getNumberOfBoundary() > 0)
	{
		//if (id == 0)
		{
			if (W.getNumberOfBoundary() > n_CUDA_afls)
			{	
				n_CUDA_vtxs.resize(W.getNumberOfBoundary(), 0);
				n_CUDA_i0s.resize(W.getNumberOfBoundary(), 0);
				n_CUDA_i1s.resize(W.getNumberOfBoundary(), 0);
				n_CUDA_i2s.resize(W.getNumberOfBoundary(), 0);
				n_CUDA_i3s.resize(W.getNumberOfBoundary(), 0);
				n_CUDA_rads.resize(W.getNumberOfBoundary(), 0);
				n_CUDA_vels.resize(W.getNumberOfBoundary(), 0);
				
				n_CUDA_rs.resize(W.getNumberOfBoundary(), 0);
				n_CUDA_rhss.resize(W.getNumberOfBoundary(), 0);

				for (size_t s = 0; s < n_CUDA_afls; ++s)
				{
					cuDeleteFromDev(W.getBoundary(s).virtualWake.devWakePtr);
					cuDeleteFromDev(W.getBoundary(s).virtualWake.devI0Ptr);
					cuDeleteFromDev(W.getBoundary(s).virtualWake.devI1Ptr);
					cuDeleteFromDev(W.getBoundary(s).virtualWake.devI2Ptr);
					cuDeleteFromDev(W.getBoundary(s).virtualWake.devI3Ptr);
					cuDeleteFromDev(W.getBoundary(s).virtualWake.devRadsPtr);
					cuDeleteFromDev(W.getBoundary(s).virtualWake.devVelsPtr);
					
					cuDeleteFromDev(W.getBoundary(s).afl.devR);
					cuDeleteFromDev(W.getBoundary(s).afl.devRhs);
				}

				if (n_CUDA_afls)
				{
					cuDeleteFromDev(dev_nPanels);
					cuDeleteFromDev(dev_nVortices);

					cuDeleteFromDev(dev_ptr_ptr_vtx);
					cuDeleteFromDev(dev_ptr_ptr_i0);
					cuDeleteFromDev(dev_ptr_ptr_i1);
					cuDeleteFromDev(dev_ptr_ptr_i2);
					cuDeleteFromDev(dev_ptr_ptr_i3);
					cuDeleteFromDev(dev_ptr_ptr_rad);
					cuDeleteFromDev(dev_ptr_ptr_vel);
					
					cuDeleteFromDev(dev_ptr_ptr_r);
					cuDeleteFromDev(dev_ptr_ptr_rhs);
				}

				std::vector<size_t> host_nPanels(0);
				std::vector<size_t> host_nVortices(0);
				
				n_CUDA_afls = 0;

				for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
				{
					const size_t& szPan = W.getBoundary(s).afl.np;
					const size_t& szVtx = std::max(W.getBoundary(s).virtualWake.vtx.size(), W.getPassport().wakeDiscretizationProperties.vortexPerPanel * W.getBoundary(s).afl.np);

					//std::cout << "szVtx = " << szVtx << std::endl;

					W.getBoundary(s).virtualWake.devWakePtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(szVtx, n_CUDA_vtxs[s]);
					W.getBoundary(s).virtualWake.devI0Ptr = ReserveDevMem<double, 1>(szVtx, n_CUDA_i1s[s]);
					W.getBoundary(s).virtualWake.devI1Ptr = ReserveDevMem<double, 1>(szVtx, n_CUDA_i1s[s]);
					W.getBoundary(s).virtualWake.devI2Ptr = ReserveDevMem<double, 2>(szVtx, n_CUDA_i2s[s]);
					W.getBoundary(s).virtualWake.devI3Ptr = ReserveDevMem<double, 2>(szVtx, n_CUDA_i2s[s]);
					W.getBoundary(s).virtualWake.devRadsPtr = ReserveDevMem<double, 1>(szVtx, n_CUDA_rads[s]);
					W.getBoundary(s).virtualWake.devVelsPtr = ReserveDevMem<double, 2>(szVtx, n_CUDA_vels[s]);
					
					W.getBoundary(s).afl.devR = ReserveDevMem<double, 2>(szPan + 1, n_CUDA_rs[s]);
					W.getBoundary(s).afl.devRhs = ReserveDevMem<double, 1>(szPan, n_CUDA_rhss[s]);

					n_CUDA_afls++;
					host_nVortices.push_back(n_CUDA_vtxs[s]);
					host_nPanels.push_back(n_CUDA_rhss[s]);
					

					W.getBoundary(s).virtualWake.tmpVels.resize(n_CUDA_vtxs[s], { 0.0, 0.0 });
					W.getBoundary(s).virtualWake.tmpRads.resize(n_CUDA_vtxs[s], 1.0e+5);

					W.getBoundary(s).virtualWake.tmpI0.resize(n_CUDA_vtxs[s], 0.0);
					W.getBoundary(s).virtualWake.tmpI1.resize(n_CUDA_vtxs[s], 0.0);
					W.getBoundary(s).virtualWake.tmpI2.resize(n_CUDA_vtxs[s], { 0.0, 0.0 });
					W.getBoundary(s).virtualWake.tmpI3.resize(n_CUDA_vtxs[s], { 0.0, 0.0 });
					
					W.getBoundary(s).afl.tmpRhs.resize(n_CUDA_rhss[s], 0.0);
				}

				dev_nPanels = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_nPanels.data());
				dev_nVortices = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_nVortices.data());

				std::vector<double*> host_ptr_ptr_vtx_tmp;
				std::vector<double*> host_ptr_ptr_vel_tmp;
				std::vector<double*> host_ptr_ptr_rad_tmp;
				std::vector<double*> host_ptr_ptr_i0_tmp;
				std::vector<double*> host_ptr_ptr_i1_tmp;
				std::vector<double*> host_ptr_ptr_i2_tmp;
				std::vector<double*> host_ptr_ptr_i3_tmp;

				std::vector<double*> host_ptr_ptr_r_tmp;
				std::vector<double*> host_ptr_ptr_rhs_tmp;

				for (size_t q = 0; q < W.getNumberOfBoundary(); ++q)
				{
					host_ptr_ptr_vtx_tmp.push_back(W.getBoundary(q).virtualWake.devWakePtr);
					host_ptr_ptr_vel_tmp.push_back(W.getBoundary(q).virtualWake.devVelsPtr);
					host_ptr_ptr_rad_tmp.push_back(W.getBoundary(q).virtualWake.devRadsPtr);
					host_ptr_ptr_i0_tmp.push_back(W.getBoundary(q).virtualWake.devI0Ptr);
					host_ptr_ptr_i1_tmp.push_back(W.getBoundary(q).virtualWake.devI1Ptr);
					host_ptr_ptr_i2_tmp.push_back(W.getBoundary(q).virtualWake.devI2Ptr);
					host_ptr_ptr_i3_tmp.push_back(W.getBoundary(q).virtualWake.devI3Ptr);
					
					host_ptr_ptr_r_tmp.push_back(W.getBoundary(q).afl.devR);
					host_ptr_ptr_rhs_tmp.push_back(W.getBoundary(q).afl.devRhs);
				}

				dev_ptr_ptr_vtx = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_ptr_vtx_tmp.data());
				dev_ptr_ptr_rad = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_ptr_rad_tmp.data());
				dev_ptr_ptr_vel = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_ptr_vel_tmp.data());
				dev_ptr_ptr_i0 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_ptr_i0_tmp.data());
				dev_ptr_ptr_i1 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_ptr_i1_tmp.data());
				dev_ptr_ptr_i2 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_ptr_i2_tmp.data());
				dev_ptr_ptr_i3 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_ptr_i3_tmp.data());
				
				dev_ptr_ptr_r = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_ptr_r_tmp.data());
				dev_ptr_ptr_rhs = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_ptr_rhs_tmp.data());
			}

			for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
			{
				cuClearWakeMem(n_CUDA_vtxs[s], W.getBoundary(s).virtualWake.devWakePtr);
				cuCopyWakeToDev(W.getBoundary(s).virtualWake.vtx.size(), W.getBoundary(s).virtualWake.vtx.data(), W.getBoundary(s).virtualWake.devWakePtr);
				cuCopyRsToDev(W.getBoundary(s).afl.r.size(), W.getBoundary(s).afl.r.data(), W.getBoundary(s).afl.devR);
			}
		}
	}
}


#endif