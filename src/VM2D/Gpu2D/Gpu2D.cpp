/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.11   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2022/08/07     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Gpu2D.cpp                                                        |
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
\brief Файл кода с описанием класса Gpu
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.11
\date 07 августа 2022 г.
*/

#include "Gpu2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Tree2D.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;


//int Gpu::nReserve = 0;

Gpu::Gpu(const World2D& W_)
	: W(W_)
{
#if defined(__CUDACC__) || defined(USE_CUDA)


// Откомментировать следующую строку, если запускается счет на кластере, каждый узел которого  
// имеет несколько видеокарт, при этом хочется одновременно решать несколько задач --- каждую на своей видеокарте ---
// каждая задача по своему номеру, деленному по модулю числа видеокарт на узле будет привязана к своей видеокарте;
// MPI-распаралеливание внутри каждой задачи НЕ ПРОИЗВОДИТЬ;
// на каждый узел при этом отправлять СТОЛЬКО MPI-нитей, СКОЛЬКО ТАМ ВИДЕОКАРТ;
// число задач НАСТОЯТЕЛЬНО РЕКОМЕНДУЕТСЯ ВЫБИРАТЬ ТОЧНО РАВНЫМ СУММАРНОМУ ЧИСЛУ ВИДЕОКАРТ,
// т.е. чтобы все задачи стартовали сразу же.
//
// Uncomment the following string if the program runs on the computer cluster with several graphic cards on every node
// and you want to solve several tasks simultaneously --- EVERY TASK ON ITS OWN GRAPHIC CARD;
// every task will be associated with separate graphic card;
// DO NOT use MPI-parallelization within one task;
// send THE SAME AMOUNT OF MPI-THREADS for the node as THE NUMBER OF GRAPHIC CARDS on this node;
// IT IS STRONGLY RECOMMENDED TO CHOOSE THE NUMBER OF TASKS EXACTLY EQUAL TO TOTAL VIDEO CARDs NUMBERS,
// i.e. to start all the tasks simultaneously.
										
//cuDevice(W.getPassport().problemNumber % 4); //The index of the used video card will be equal to the task number
                                             // in the task list (to modulo 4 --- number of graphic cards on each node) 


	
// Откомментировать следующую строку, если запускается счет на кластере, каждый узел которого  
// имеет несколько видеокарт, при этом хочется каждую задачу решать при помощи нескольких видеокарт ---
// каждый MPI-процесс, решающий текущую задачу, будет привязан к своей видеокарте;
// число порождаемых MPI-нитей (на каждую задачу) должно было не больше числа видеокарт;
// на каждый узел при этом отправлять ТОЛЬКО ОДНУ задачу;
// число задач может быть любым, они будут исполняться по очереди
//
//
// Uncomment the following string if the program runs on the computer cluster with several graphic cards on every node
// and you want to solve EVERY TASK USING SEVERAL GRAPHIC CARDS ---
// every MPI-thread (solving the current task) will be associated with its own graphic card;
// the number of MPI-threads (for every task) have to be not much than number of graphic cards;
// send only one task to every node;
// number of problems can be arbitrary, they will be executed in turn.

// cuDevice(W.getParallel().myidWork);          //The index of the used video card will be equal to the number of MPI-thread


	   
	cuSetConstants(sizeof(Vortex2D)/sizeof(double), Vortex2D::offsPos / sizeof(double), Vortex2D::offsGam / sizeof(double) );
	
	n_CUDA_wake = 0;
	n_CUDA_source = 0;
	n_CUDA_afls = 0;

	n_CUDA_velVP = 0;
#endif
}


Gpu::~Gpu()
{
#if defined(__CUDACC__) || defined(USE_CUDA)
	ReleaseDevMem(W.getWake().devVtxPtr, 1);
	ReleaseDevMem(W.getWake().devVelPtr, 2);
	ReleaseDevMem(W.getWake().devRadPtr, 3);
	ReleaseDevMem(W.getWake().devI0Ptr, 4);
	ReleaseDevMem(W.getWake().devI1Ptr, 5);
	ReleaseDevMem(W.getWake().devI2Ptr, 6);
	ReleaseDevMem(W.getWake().devI3Ptr, 7);

	ReleaseDevMem(W.getWake().devMeshPtr, 8);
	ReleaseDevMem(W.getWake().devNeiPtr, 9);

	if (W.getSource().vtx.size() > 0)
		ReleaseDevMem(W.getSource().devVtxPtr, 10);

	for (size_t s = 0; s < 1/*n_CUDA_afls*/; ++s)
	{
		ReleaseDevMem(W.getBoundary(s).virtualWake.devVtxPtr, 11);
		ReleaseDevMem(W.getBoundary(s).virtualWake.devVelPtr, 12);
		ReleaseDevMem(W.getBoundary(s).virtualWake.devRadPtr, 13);
		ReleaseDevMem(W.getBoundary(s).virtualWake.devI0Ptr, 14);
		ReleaseDevMem(W.getBoundary(s).virtualWake.devI1Ptr, 15);
		ReleaseDevMem(W.getBoundary(s).virtualWake.devI2Ptr, 16);
		ReleaseDevMem(W.getBoundary(s).virtualWake.devI3Ptr, 17);

		ReleaseDevMem(W.getBoundary(s).afl.devRPtr, 18);
		ReleaseDevMem(W.getBoundary(s).afl.devRhsPtr, 19);
		ReleaseDevMem(W.getBoundary(s).afl.devRhsLinPtr, 191);

		ReleaseDevMem(W.getBoundary(s).afl.devFreeVortexSheetPtr, 20);
		ReleaseDevMem(W.getBoundary(s).afl.devAttachedVortexSheetPtr, 21);
		ReleaseDevMem(W.getBoundary(s).afl.devAttachedSourceSheetPtr, 22);

		ReleaseDevMem(W.getBoundary(s).afl.devMeanEpsOverPanelPtr, 23);
		ReleaseDevMem(W.getBoundary(s).afl.devViscousStressesPtr, 24);
	}

	if (n_CUDA_afls)
	{		
		ReleaseDevMem(dev_ptr_nPanels, 25);
		ReleaseDevMem(dev_ptr_nVortices, 26);

		ReleaseDevMem(dev_ptr_ptr_vtx, 27);
		ReleaseDevMem(dev_ptr_ptr_vel, 28);
		ReleaseDevMem(dev_ptr_ptr_rad, 29);
		ReleaseDevMem(dev_ptr_ptr_i0, 30);
		ReleaseDevMem(dev_ptr_ptr_i1, 31);
		ReleaseDevMem(dev_ptr_ptr_i2, 32);
		ReleaseDevMem(dev_ptr_ptr_i3, 33);
				
		ReleaseDevMem(dev_ptr_ptr_r, 34);
		ReleaseDevMem(dev_ptr_ptr_rhs, 35);

		ReleaseDevMem(dev_ptr_ptr_freeVortexSheet, 36);
		ReleaseDevMem(dev_ptr_ptr_attachedVortexSheet, 37);
		ReleaseDevMem(dev_ptr_ptr_attachedSourceSheet, 38);

		ReleaseDevMem(dev_ptr_ptr_meanEpsOverPanel, 39);

		ReleaseDevMem(dev_ptr_ptr_viscousStresses, 40);
	}
		
	if (W.getMeasureVP().getWakeVP().vtx.size())
	{
		ReleaseDevMem(W.getMeasureVP().getWakeVP().devVtxPtr, 41);
		ReleaseDevMem(W.getMeasureVP().getWakeVP().devVelPtr, 42);
		ReleaseDevMem(W.getMeasureVP().getWakeVP().devRadPtr, 43);
	}
#endif
}

#if defined(__CUDACC__) || defined(USE_CUDA)


//Обновление состояния следа wake
void Gpu::RefreshWake(int code) 
{
	if (W.getWake().vtx.size() > 0)
	{		
		//Если зарезервировано меньше, чем вихрей в пелене
		if (W.getWake().vtx.size() > n_CUDA_wake) 
		{
			size_t curLength = n_CUDA_wake;

			//Освобождаем всю память на видеокарте
			if (curLength > 0)
			{				
				ReleaseDevMem(W.getWake().devVtxPtr, 44);
				ReleaseDevMem(W.getWake().devVelPtr, 45);
				ReleaseDevMem(W.getWake().devRadPtr, 46);
				ReleaseDevMem(W.getWake().devI0Ptr, 47);
				ReleaseDevMem(W.getWake().devI1Ptr, 48);
				ReleaseDevMem(W.getWake().devI2Ptr, 49);
				ReleaseDevMem(W.getWake().devI3Ptr, 50);

				ReleaseDevMem(W.getWake().devMeshPtr, 51);
				ReleaseDevMem(W.getWake().devNeiPtr, 52);
			}

			size_t sz = curLength;
			while (W.getWake().vtx.size() > sz)
				sz += INC_VORT_DEV;

			//Резервируем новое количество памяти
			W.getWake().devVtxPtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(sz, n_CUDA_wake);

			W.getWake().devVelPtr = ReserveDevMem<double, 2>(sz, n_CUDA_wake);
			W.getWake().devRadPtr = ReserveDevMem<double, 1>(sz, n_CUDA_wake);

			W.getWake().devI0Ptr = ReserveDevMem<double, 1>(sz, n_CUDA_wake);
			W.getWake().devI1Ptr = ReserveDevMem<double, 1>(sz, n_CUDA_wake);
			W.getWake().devI2Ptr = ReserveDevMem<double, 2>(sz, n_CUDA_wake);
			W.getWake().devI3Ptr = ReserveDevMem<double, 2>(sz, n_CUDA_wake);

			W.getWake().devMeshPtr = ReserveDevMem<int, 2>(sz, n_CUDA_wake);
			W.getWake().devNeiPtr = ReserveDevMem<int, 1>(sz, n_CUDA_wake);
			
			//mesh.resize(n_CUDA_vel, { 0, 0 });

			W.getInfo('i') << "CUDA memory resize: " << curLength << " -> " << n_CUDA_wake << " vortices" << std::endl;
		}// if (W.getWake().vtx.size() > n_CUDA_wake)

		//Обнуляем память, выделенную выше для хранения следа		
		cuClearWakeMem(W.getWake().vtx.size(), W.getWake().devVtxPtr);

		//Копирование следа на видеокарту
		//double t1 = omp_get_wtime();
		//cuCopyWakeToDevAsync(W.getWake().vtx.size(), W.getWake().vtx.data(), W.getWake().devVtxPtr, 1);
		cuCopyWakeToDev(W.getWake().vtx.size(), W.getWake().vtx.data(), W.getWake().devVtxPtr, 1);
		//double t2 = omp_get_wtime();
		//std::cout << "CopyTime = " << t2 - t1 << std::endl;
	}

	if (W.getSource().vtx.size() > 0)
	{
		//Если зарезервировано меньше, чем источников в пелене
		if (W.getSource().vtx.size() > n_CUDA_source) 
		{
			size_t curLength = n_CUDA_source; 
			
			//Освобождаем всю память на ведеокарте
			if (curLength > 0)				
				ReleaseDevMem(W.getSource().devVtxPtr, 53);
			
			size_t sz = curLength;
			while (W.getSource().vtx.size() > sz)
				sz += INC_VORT_DEV;

			//Резервируем новое количество памяти			
			W.getSource().devVtxPtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(sz, n_CUDA_source);


			W.getInfo('i') << "CUDA memory resize: " << curLength << " -> " << sz << " sources" << std::endl;
		}// if (W.getSource().vtx.size() > N_CUDA_source)


		//Обнуляем память, выделенную выше для хранения следа
		cuClearWakeMem(W.getSource().vtx.size(), W.getSource().devVtxPtr);

		//Копирование следа на видеокарту		
		cuCopyWakeToDev(W.getSource().vtx.size(), W.getSource().vtx.data(), W.getSource().devVtxPtr, 2);
	
		
	}
}



//Обновление состояния сетки для вычисления VP
void Gpu::RefreshVP(int code)
{
	if (W.getMeasureVP().getWakeVP().vtx.size() > 0)
	{
		//Если зарезервировано меньше, чем вихрей в пелене
		if (W.getMeasureVP().getWakeVP().vtx.size() > n_CUDA_velVP)
		{
			size_t curLength = n_CUDA_velVP;

			//Освобождаем всю память на видеокарте
			if (curLength > 0)
			{
				ReleaseDevMem(W.getMeasureVP().getWakeVP().devVtxPtr, 54);
				ReleaseDevMem(W.getMeasureVP().getWakeVP().devVelPtr, 55);
				ReleaseDevMem(W.getMeasureVP().getWakeVP().devRadPtr, 56);
			}

			size_t sz = curLength;
			while (W.getMeasureVP().getWakeVP().vtx.size() > sz)
				sz += INC_VORT_DEV;

			//Резервируем новое количество памяти
			W.getMeasureVP().getWakeVP().devVtxPtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(sz, n_CUDA_velVP);
			W.getMeasureVP().getWakeVP().devVelPtr = ReserveDevMem<double, 2>(sz, n_CUDA_velVP);
			W.getMeasureVP().getWakeVP().devRadPtr = ReserveDevMem<double, 1>(sz, n_CUDA_velVP);

			W.getInfo('i') << "CUDA memory resize: " << curLength << " -> " << n_CUDA_velVP << " points_VP" << std::endl;
		}// if (W.getWake().vtx.size() > N_CUDA_velvp)

		//Обнуляем память, выделенную выше для хранения следа
		cuClearWakeMem(W.getMeasureVP().getWakeVP().vtx.size(), W.getMeasureVP().getWakeVP().devVtxPtr);

		//Копирование следа на видеокарту
		cuCopyWakeToDev(W.getMeasureVP().getWakeVP().vtx.size(), W.getMeasureVP().getWakeVP().vtx.data(), W.getMeasureVP().getWakeVP().devVtxPtr, 3);
	}
}



//Обновление состояния всех профилей и слоев на них
void Gpu::RefreshAfls(int code)
{
	//std::cout << "RefreshAfls (code = " << code << ") start" << std::endl;

	if (W.getNumberOfBoundary() > 0)
	{
		if (W.getNumberOfBoundary() > n_CUDA_afls)
		{
			//Обнуление массива на хосте, содержащего количество панелей на профилях
			n_CUDA_panel.resize(W.getNumberOfBoundary(), 0);

			if (n_CUDA_afls)
				for (size_t s = 0; s < 1/*n_CUDA_afls*/; ++s)
				{
					ReleaseDevMem(W.getBoundary(s).afl.devRPtr, 57);
					ReleaseDevMem(W.getBoundary(s).afl.devRhsPtr, 58);
					ReleaseDevMem(W.getBoundary(s).afl.devRhsLinPtr, 581);
					ReleaseDevMem(W.getBoundary(s).afl.devMeanEpsOverPanelPtr, 59);

					ReleaseDevMem(W.getBoundary(s).afl.devFreeVortexSheetPtr, 60);
					ReleaseDevMem(W.getBoundary(s).afl.devAttachedVortexSheetPtr, 61);
					ReleaseDevMem(W.getBoundary(s).afl.devAttachedSourceSheetPtr, 62);

					ReleaseDevMem(W.getBoundary(s).afl.devViscousStressesPtr, 63);
				}

			if (n_CUDA_afls)
			{
				ReleaseDevMem(dev_ptr_nPanels, 64);
				ReleaseDevMem(dev_ptr_ptr_r, 65);
				ReleaseDevMem(dev_ptr_ptr_rhs, 66);

				ReleaseDevMem(dev_ptr_ptr_freeVortexSheet, 67);
				ReleaseDevMem(dev_ptr_ptr_attachedVortexSheet, 68);
				ReleaseDevMem(dev_ptr_ptr_attachedSourceSheet, 69);

				ReleaseDevMem(dev_ptr_nVortices, 70);
				ReleaseDevMem(dev_ptr_ptr_vtx, 71);
				ReleaseDevMem(dev_ptr_ptr_rad, 72);
				ReleaseDevMem(dev_ptr_ptr_vel, 73);
				ReleaseDevMem(dev_ptr_ptr_i0, 74);
				ReleaseDevMem(dev_ptr_ptr_i1, 75);
				ReleaseDevMem(dev_ptr_ptr_i2, 76);
				ReleaseDevMem(dev_ptr_ptr_i3, 77);
			}

			//Временный массив для агрегации числа панелей (без округления вверх до блока) на профилях для их последующей отправки в dev_ptr_nPanels
			std::vector<size_t> host_nPanels(0);
			
			size_t totnPanels = 0, totnPanelsUP;
			for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
			{
				const size_t& nps = W.getBoundary(s).afl.getNumberOfPanels();
				host_nPanels.push_back(nps);
				totnPanels += nps;
			}
			dev_ptr_nPanels = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_nPanels.data());

			const int totNVars = (int)totnPanels * (W.getPassport().numericalSchemes.boundaryCondition.second + 1);
			
			W.getBoundary(0).afl.devRPtr = ReserveDevMem<double, 4>(totnPanels, totnPanelsUP);
			W.getBoundary(0).afl.devRhsPtr = ReserveDevMem<double, 1>(totnPanels, totnPanelsUP);
			W.getBoundary(0).afl.devRhsLinPtr = ReserveDevMem<double, 1>(totnPanels, totnPanelsUP);
			W.getBoundary(0).afl.devMeanEpsOverPanelPtr = ReserveDevMem<double, 1>(totnPanels, totnPanelsUP);//MeanEpsOverPanel

			W.getBoundary(0).afl.devFreeVortexSheetPtr = ReserveDevMem<double, 1>(totNVars, totnPanelsUP);//devFreeVortexSheetPtr
			W.getBoundary(0).afl.devAttachedVortexSheetPtr = ReserveDevMem<double, 1>(totNVars, totnPanelsUP);//devAttachedVortexSheetPtr
			W.getBoundary(0).afl.devAttachedSourceSheetPtr = ReserveDevMem<double, 1>(totNVars, totnPanelsUP);//devAttachedSourceSheetPtr

			W.getBoundary(0).afl.devViscousStressesPtr = ReserveDevMem<double, 1>(totnPanels, totnPanelsUP);//ViscousStress

			for (size_t s = 1; s < W.getNumberOfBoundary(); ++s)
			{
				W.getBoundary(s).afl.devRPtr = W.getBoundary(s - 1).afl.devRPtr + 4 * host_nPanels[s - 1];
				W.getBoundary(s).afl.devRhsPtr = W.getBoundary(s - 1).afl.devRhsPtr + 1 * host_nPanels[s - 1];
				W.getBoundary(s).afl.devRhsLinPtr = W.getBoundary(s - 1).afl.devRhsLinPtr + 1 * host_nPanels[s - 1];
				W.getBoundary(s).afl.devMeanEpsOverPanelPtr = W.getBoundary(s - 1).afl.devMeanEpsOverPanelPtr + 1 * host_nPanels[s - 1];

				W.getBoundary(s).afl.devFreeVortexSheetPtr = W.getBoundary(s - 1).afl.devFreeVortexSheetPtr + 1 * host_nPanels[s - 1];
				W.getBoundary(s).afl.devAttachedVortexSheetPtr = W.getBoundary(s - 1).afl.devAttachedVortexSheetPtr + 1 * host_nPanels[s - 1];
				W.getBoundary(s).afl.devAttachedSourceSheetPtr = W.getBoundary(s - 1).afl.devAttachedSourceSheetPtr + 1 * host_nPanels[s - 1];

				W.getBoundary(s).afl.devViscousStressesPtr = W.getBoundary(s - 1).afl.devViscousStressesPtr + 1 * host_nPanels[s - 1];
			}

			for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
			{
				W.getBoundary(s).afl.tmpRhs.resize(totNVars, 0.0);
				W.getBoundary(s).afl.tmpViscousStresses.resize(W.getBoundary(s).afl.getNumberOfPanels(), 0.0);
			}// for s

			//Временные массивы для агрегации указателей на видеокарте
			std::vector<double*> host_ptr_r;
			std::vector<double*> host_ptr_rhs;
			std::vector<double*> host_ptr_rhsLin;

			std::vector<double*> host_ptr_freeVortexSheet;
			std::vector<double*> host_ptr_attachedVortexSheet;
			std::vector<double*> host_ptr_attachedSourceSheet;

			std::vector<double*> host_ptr_meanEpsOverPanel;

			std::vector<double*> host_ptr_viscousStresses;

			for (size_t q = 0; q < W.getNumberOfBoundary(); ++q)
			{
				host_ptr_r.push_back(W.getBoundary(q).afl.devRPtr);
				host_ptr_rhs.push_back(W.getBoundary(q).afl.devRhsPtr);
				host_ptr_rhsLin.push_back(W.getBoundary(q).afl.devRhsLinPtr);

				host_ptr_freeVortexSheet.push_back(W.getBoundary(q).afl.devFreeVortexSheetPtr);
				host_ptr_attachedVortexSheet.push_back(W.getBoundary(q).afl.devAttachedVortexSheetPtr);
				host_ptr_attachedSourceSheet.push_back(W.getBoundary(q).afl.devAttachedSourceSheetPtr);

				host_ptr_meanEpsOverPanel.push_back(W.getBoundary(q).afl.devMeanEpsOverPanelPtr);

				host_ptr_viscousStresses.push_back(W.getBoundary(q).afl.devViscousStressesPtr);
			}

			dev_ptr_ptr_r = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_r.data());
			dev_ptr_ptr_rhs = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_rhs.data());

			dev_ptr_ptr_freeVortexSheet = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_freeVortexSheet.data());
			dev_ptr_ptr_attachedVortexSheet = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_attachedVortexSheet.data());
			dev_ptr_ptr_attachedSourceSheet = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_attachedSourceSheet.data());

			dev_ptr_ptr_viscousStresses = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_viscousStresses.data());

			dev_ptr_ptr_meanEpsOverPanel = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_meanEpsOverPanel.data());


			std::vector<double*> zeroPtrVec(W.getNumberOfBoundary(), nullptr);			
			dev_ptr_ptr_vtx = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), zeroPtrVec.data());
			dev_ptr_ptr_rad = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), zeroPtrVec.data());
			dev_ptr_ptr_vel = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), zeroPtrVec.data());
			dev_ptr_ptr_i0 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), zeroPtrVec.data());
			dev_ptr_ptr_i1 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), zeroPtrVec.data());
			dev_ptr_ptr_i2 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), zeroPtrVec.data());
			dev_ptr_ptr_i3 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), zeroPtrVec.data());

			std::vector<size_t> zeroVec(W.getNumberOfBoundary(), 0);
			dev_ptr_nVortices = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), zeroVec.data());

		}//if (W.getNumberOfBoundary() > n_CUDA_afls)


		for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
		{
			size_t np = W.getBoundary(s).afl.getNumberOfPanels();
			size_t nv = W.getBoundary(s).GetUnknownsSize();

			//Копирование вершин профиля на видеокарту			 
			std::vector<double> rbegend(4 * np);
			for (size_t q = 0; q < np; ++q)
			{
				rbegend[4 * q + 0] = W.getBoundary(s).afl.getR(q)[0];
				rbegend[4 * q + 1] = W.getBoundary(s).afl.getR(q)[1];
				rbegend[4 * q + 2] = W.getBoundary(s).afl.getR(q+1)[0];
				rbegend[4 * q + 3] = W.getBoundary(s).afl.getR(q+1)[1];
			}
	
			cuCopyFixedArrayPoint4D(W.getBoundary(s).afl.devRPtr, (Point2D*)rbegend.data(), np, code);
			
			//Копирование слоев на видеокарту			
			std::vector<double> host_freeVortexSheet(nv);
			std::vector<double> host_attachedVortexSheet(nv);
			std::vector<double> host_attachedSourceSheet(nv);

			const Sheet& sh = W.getBoundary(s).sheets;

			int sch = W.getPassport().numericalSchemes.boundaryCondition.second;

			for (size_t p = 0; p < np; ++p)
			{
				host_attachedVortexSheet[p] = sh.attachedVortexSheet(p, 0);
				host_attachedSourceSheet[p] = sh.attachedSourceSheet(p, 0);

				if (W.getParallel().myidWork == 0)
					host_freeVortexSheet[p] = sh.freeVortexSheet(p, 0);
			}//for p

			if(sch == 1)
			for (size_t p = 0; p < np; ++p)
			{
				host_attachedVortexSheet[np + p] = sh.attachedVortexSheet(p, 1);
				host_attachedSourceSheet[np + p] = sh.attachedSourceSheet(p, 1);

				if (W.getParallel().myidWork == 0)
					host_freeVortexSheet[np + p] = sh.freeVortexSheet(p, 1);
			}//for p

			MPI_Bcast(host_freeVortexSheet.data(), (int)host_freeVortexSheet.size(), MPI_DOUBLE, 0, W.getParallel().commWork);

			cuCopyFixedArray(W.getBoundary(s).afl.devFreeVortexSheetPtr, host_freeVortexSheet.data(), sizeof(double) * host_freeVortexSheet.size());
			cuCopyFixedArray(W.getBoundary(s).afl.devAttachedVortexSheetPtr, host_attachedVortexSheet.data(), sizeof(double)* host_attachedVortexSheet.size());
			cuCopyFixedArray(W.getBoundary(s).afl.devAttachedSourceSheetPtr, host_attachedSourceSheet.data(), sizeof(double)* host_attachedSourceSheet.size());
		}
	}	
	n_CUDA_afls = W.getNumberOfBoundary();
}//RefreshAfls()


//Обновление состояния "виртуальных следов" - только что рожденных вихрей на профилях
void Gpu::RefreshVirtualWakes(int code)
{	
	if (n_CUDA_virtWake.size() == 0)
	{
		n_CUDA_virtWake.resize(W.getNumberOfBoundary(), 0);
		n_CUDA_totalVirtWake = 0;
	}

	size_t totnVirt = 0, totnPan = 0;
	for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
	{
		totnVirt += W.getBoundary(s).virtualWake.vtx.size();
		totnPan += W.getBoundary(s).afl.getNumberOfPanels();
	}

	const size_t& szVtx = std::max(totnVirt, W.getPassport().wakeDiscretizationProperties.minVortexPerPanel * totnPan);
	
	//Чистим все массивы
	if (szVtx > n_CUDA_totalVirtWake)
	{
		if (n_CUDA_totalVirtWake > 0)
			for (size_t s = 0; s < 1/*W.getNumberOfBoundary()*/; ++s)
			{
				ReleaseDevMem(W.getBoundary(s).virtualWake.devVtxPtr, 70);
				ReleaseDevMem(W.getBoundary(s).virtualWake.devVelPtr, 71);
				ReleaseDevMem(W.getBoundary(s).virtualWake.devRadPtr, 72);
				ReleaseDevMem(W.getBoundary(s).virtualWake.devI0Ptr, 73);
				ReleaseDevMem(W.getBoundary(s).virtualWake.devI1Ptr, 74);
				ReleaseDevMem(W.getBoundary(s).virtualWake.devI2Ptr, 75);
				ReleaseDevMem(W.getBoundary(s).virtualWake.devI3Ptr, 76);
			}

		if (n_CUDA_afls)
		{
			W.getBoundary(0).virtualWake.devVtxPtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(szVtx, n_CUDA_totalVirtWake);
			W.getBoundary(0).virtualWake.devVelPtr = ReserveDevMem<double, 2>(szVtx, n_CUDA_totalVirtWake);
			W.getBoundary(0).virtualWake.devRadPtr = ReserveDevMem<double, 1>(szVtx, n_CUDA_totalVirtWake);
			W.getBoundary(0).virtualWake.devI0Ptr = ReserveDevMem<double, 1>(szVtx, n_CUDA_totalVirtWake);
			W.getBoundary(0).virtualWake.devI1Ptr = ReserveDevMem<double, 1>(szVtx, n_CUDA_totalVirtWake);
			W.getBoundary(0).virtualWake.devI2Ptr = ReserveDevMem<double, 2>(szVtx, n_CUDA_totalVirtWake);
			W.getBoundary(0).virtualWake.devI3Ptr = ReserveDevMem<double, 2>(szVtx, n_CUDA_totalVirtWake);
		}
	}

	for (size_t s = 1; s < W.getNumberOfBoundary(); ++s)
	{
		size_t nprev = W.getBoundary(s - 1).virtualWake.vtx.size();
		W.getBoundary(s).virtualWake.devVtxPtr = W.getBoundary(s - 1).virtualWake.devVtxPtr + (sizeof(Vortex2D) / sizeof(double)) * nprev;
		W.getBoundary(s).virtualWake.devVelPtr = W.getBoundary(s - 1).virtualWake.devVelPtr + 2 * nprev;
		W.getBoundary(s).virtualWake.devRadPtr = W.getBoundary(s - 1).virtualWake.devRadPtr + 1 * nprev;
		W.getBoundary(s).virtualWake.devI0Ptr = W.getBoundary(s - 1).virtualWake.devI0Ptr + 1 * nprev;
		W.getBoundary(s).virtualWake.devI1Ptr = W.getBoundary(s - 1).virtualWake.devI1Ptr + 1 * nprev;
		W.getBoundary(s).virtualWake.devI2Ptr = W.getBoundary(s - 1).virtualWake.devI2Ptr + 2 * nprev;
		W.getBoundary(s).virtualWake.devI3Ptr = W.getBoundary(s - 1).virtualWake.devI3Ptr + 2 * nprev;
	}
	
	std::vector<size_t> host_nVortices(0);
	host_nVortices.reserve(W.getNumberOfBoundary());

	for (size_t q = 0; q < W.getNumberOfBoundary(); ++q)
		host_nVortices.push_back(W.getBoundary(q).virtualWake.vtx.size());

	cuCopyFixedArray(dev_ptr_nVortices, host_nVortices.data(), W.getNumberOfBoundary()*sizeof(size_t));

	for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
	{
		//Временные массивы для агрегации указателей на видеокарте
		std::vector<double*> host_ptr_vtx;
		std::vector<double*> host_ptr_vel;
		std::vector<double*> host_ptr_rad;
		std::vector<double*> host_ptr_i0;
		std::vector<double*> host_ptr_i1;
		std::vector<double*> host_ptr_i2;
		std::vector<double*> host_ptr_i3;

		for (size_t q = 0; q < W.getNumberOfBoundary(); ++q)
		{
			host_ptr_vtx.push_back(W.getBoundary(q).virtualWake.devVtxPtr);
			host_ptr_vel.push_back(W.getBoundary(q).virtualWake.devVelPtr);
			host_ptr_rad.push_back(W.getBoundary(q).virtualWake.devRadPtr);
			host_ptr_i0.push_back(W.getBoundary(q).virtualWake.devI0Ptr);
			host_ptr_i1.push_back(W.getBoundary(q).virtualWake.devI1Ptr);
			host_ptr_i2.push_back(W.getBoundary(q).virtualWake.devI2Ptr);
			host_ptr_i3.push_back(W.getBoundary(q).virtualWake.devI3Ptr);
		}

		size_t nBytes = W.getNumberOfBoundary() * sizeof(double*);
		cuCopyFixedArray(dev_ptr_ptr_vtx, host_ptr_vtx.data(), nBytes);
		cuCopyFixedArray(dev_ptr_ptr_rad, host_ptr_rad.data(), nBytes);
		cuCopyFixedArray(dev_ptr_ptr_vel, host_ptr_vel.data(), nBytes);
		cuCopyFixedArray(dev_ptr_ptr_i0, host_ptr_i0.data(), nBytes);
		cuCopyFixedArray(dev_ptr_ptr_i1, host_ptr_i1.data(), nBytes);
		cuCopyFixedArray(dev_ptr_ptr_i2, host_ptr_i2.data(), nBytes);
		cuCopyFixedArray(dev_ptr_ptr_i3, host_ptr_i3.data(), nBytes);
	}//if (n_CUDA_virtWake[s] < W.getBoundary(s).virtualWake.vtx.size())

	//Обнуляем память, выделенную выше для хранения виртуального следа
	if (W.getNumberOfBoundary() > 0)
	{
		cuClearWakeMem(szVtx, W.getBoundary(0).virtualWake.devVtxPtr);
		for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
			cuCopyWakeToDev(W.getBoundary(s).virtualWake.vtx.size(), W.getBoundary(s).virtualWake.vtx.data(), W.getBoundary(s).virtualWake.devVtxPtr, 4);
	}
}

#endif