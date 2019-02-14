/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.5    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/02/20     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
\brief Файл кода с описанием класса Gpu
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.5   
\date 20 февраля 2019 г.
*/

#include "Gpu2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Parallel.h"
#include "Passport2D.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "Wake2D.h"
#include "World2D.h"

using namespace VM2D;

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
										
// cuDevice(W.getPassport().problemNumber % 4); //Индекс используемой видеокарты будет равен номеру задачи в перечне задач (по модулю 4)


	
// Откомментировать следующую строку, если запускается счет на кластере, каждый узел которого  
// имеет несколько видеокарт, при этом хочется каждую задачу решать при помощи нескольких видеокарт ---
// каждый MPI-процесс, решающий текущую задачу, будет привязан к своей видеокарте;
// число порождаемых MPI-нитей (на каждую задачу) должно было не больше числа видеокарт;
// на каждый узел при этом отправлять ТОЛЬКО ОДНУ задачу;
// число задач может быть любым, они будут исполняться по очереди

// cuDevice(W.getParallel().myidWork);          //Индекс используемой видеокарты будет равен номеру MPI-ной нити




	cuSetConstants(sizeof(Vortex2D)/sizeof(double), Vortex2D::offsPos / sizeof(double), Vortex2D::offsGam / sizeof(double) );
	
/*
	W.getWake().devVtxPtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(INC_VORT_DEV, W.getWake().devNWake);
	
	W.getWake().devVelPtr = ReserveDevMem<double, 2>(INC_VORT_DEV, n_CUDA_wake);
	W.getWake().devRadPtr = ReserveDevMem<double, 1>(INC_VORT_DEV, n_CUDA_wake);
	W.getWake().devI0Ptr = ReserveDevMem<double, 1>(INC_VORT_DEV, n_CUDA_wake);
	W.getWake().devI1Ptr = ReserveDevMem<double, 1>(INC_VORT_DEV, n_CUDA_wake);
	W.getWake().devI2Ptr = ReserveDevMem<double, 2>(INC_VORT_DEV, n_CUDA_wake);
	W.getWake().devI3Ptr = ReserveDevMem<double, 2>(INC_VORT_DEV, n_CUDA_wake);
	

	W.getWake().devMeshPtr = ReserveDevMem<int, 2>(INC_VORT_DEV, n_CUDA_wake);
	W.getWake().devNeiPtr = ReserveDevMem<int, 1>(INC_VORT_DEV, n_CUDA_wake);

	W.getWake().tmpVel.resize(n_CUDA_wake, { 0.0, 0.0 });
	W.getWake().tmpRad.resize(n_CUDA_wake, 1.0e+5);

	W.getWake().tmpI0.resize(n_CUDA_wake, 0.0);
	W.getWake().tmpI1.resize(n_CUDA_wake, 0.0);
	W.getWake().tmpI2.resize(n_CUDA_wake, { 0.0, 0.0 });
	W.getWake().tmpI3.resize(n_CUDA_wake, { 0.0, 0.0 });
	W.getWake().tmpNei.resize(n_CUDA_wake, 0);
	//mesh.resize(n_CUDA_vel, { 0, 0 });
*/		
	n_CUDA_wake = 0;
	n_CUDA_afls = 0;

#endif
}


Gpu::~Gpu()
{
#if defined(__CUDACC__) || defined(USE_CUDA)
	cuDeleteFromDev(W.getWake().devVtxPtr);
	cuDeleteFromDev(W.getWake().devVelPtr);
	cuDeleteFromDev(W.getWake().devRadPtr);
	cuDeleteFromDev(W.getWake().devI0Ptr);
	cuDeleteFromDev(W.getWake().devI1Ptr);
	cuDeleteFromDev(W.getWake().devI2Ptr);
	cuDeleteFromDev(W.getWake().devI3Ptr);
	cuDeleteFromDev(W.getWake().devMeshPtr);
	cuDeleteFromDev(W.getWake().devNeiPtr);

	if (W.getSource().vtx.size() > 0)
		cuDeleteFromDev(W.getSource().devVtxPtr);

	for (size_t s = 0; s < n_CUDA_afls; ++s)
	{
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devVtxPtr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devVelPtr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devRadPtr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devI0Ptr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devI1Ptr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devI2Ptr);
		cuDeleteFromDev(W.getBoundary(s).virtualWake.devI3Ptr);
		
		cuDeleteFromDev(W.getBoundary(s).afl.devRPtr);
		cuDeleteFromDev(W.getBoundary(s).afl.devRhsPtr);

		cuDeleteFromDev(W.getBoundary(s).afl.devFreeVortexSheetPtr);
		cuDeleteFromDev(W.getBoundary(s).afl.devAttachedVortexSheetPtr);
		cuDeleteFromDev(W.getBoundary(s).afl.devAttachedSourceSheetPtr);
	}

	if (n_CUDA_afls)
	{
		cuDeleteFromDev(dev_ptr_nPanels);
		cuDeleteFromDev(dev_ptr_nVortices);

		cuDeleteFromDev(dev_ptr_ptr_vtx);		
		cuDeleteFromDev(dev_ptr_ptr_vel);
		cuDeleteFromDev(dev_ptr_ptr_rad);
		cuDeleteFromDev(dev_ptr_ptr_i0);
		cuDeleteFromDev(dev_ptr_ptr_i1);
		cuDeleteFromDev(dev_ptr_ptr_i2);
		cuDeleteFromDev(dev_ptr_ptr_i3);
				
		cuDeleteFromDev(dev_ptr_ptr_r);
		cuDeleteFromDev(dev_ptr_ptr_rhs);

		cuDeleteFromDev(dev_ptr_ptr_freeVortexSheet);
		cuDeleteFromDev(dev_ptr_ptr_attachedVortexSheet);
		cuDeleteFromDev(dev_ptr_ptr_attachedSourceSheet);
	}
#endif
}

#if defined(__CUDACC__) || defined(USE_CUDA)


//Обновление состояния следа wake
void Gpu::RefreshWake() 
{
	const int& id = W.getParallel().myidWork;

	if (W.getWake().vtx.size() > 0)
	{		
		//Если зарезервировано меньше, чем вихрей в пелене
		if (W.getWake().vtx.size() > W.getWake().devNWake)
		{
			size_t curLength = W.getWake().devNWake;

			//Освобождаем всю память на видеокарте
			if (curLength > 0)
			{
				cuDeleteFromDev(W.getWake().devVtxPtr);
				cuDeleteFromDev(W.getWake().devVelPtr);
				cuDeleteFromDev(W.getWake().devRadPtr);
				cuDeleteFromDev(W.getWake().devI0Ptr);
				cuDeleteFromDev(W.getWake().devI1Ptr);
				cuDeleteFromDev(W.getWake().devI2Ptr);
				cuDeleteFromDev(W.getWake().devI3Ptr);
				cuDeleteFromDev(W.getWake().devMeshPtr);
				cuDeleteFromDev(W.getWake().devNeiPtr);
			}

			size_t sz = curLength;
			while (W.getWake().vtx.size() > sz)
				sz += INC_VORT_DEV;

			//Резервируем новое количество памяти
			W.getWake().devVtxPtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(sz, W.getWake().devNWake);

			W.getWake().devVelPtr = ReserveDevMem<double, 2>(sz, n_CUDA_wake);
			W.getWake().devRadPtr = ReserveDevMem<double, 1>(sz, n_CUDA_wake);

			W.getWake().devI0Ptr = ReserveDevMem<double, 1>(sz, n_CUDA_wake);
			W.getWake().devI1Ptr = ReserveDevMem<double, 1>(sz, n_CUDA_wake);
			W.getWake().devI2Ptr = ReserveDevMem<double, 2>(sz, n_CUDA_wake);
			W.getWake().devI3Ptr = ReserveDevMem<double, 2>(sz, n_CUDA_wake);

			W.getWake().devMeshPtr = ReserveDevMem<int, 2>(sz, n_CUDA_wake);
			W.getWake().devNeiPtr = ReserveDevMem<int, 1>(sz, n_CUDA_wake);

			//Подготовка массивов на хосте для последующей переправки на видеокарту
			W.getWake().tmpVel.resize(n_CUDA_wake, { 0.0, 0.0 });
			W.getWake().tmpRad.resize(n_CUDA_wake, 1.0e+5);

			W.getWake().tmpI0.resize(n_CUDA_wake, 0.0);
			W.getWake().tmpI1.resize(n_CUDA_wake, 0.0);
			W.getWake().tmpI2.resize(n_CUDA_wake, { 0.0, 0.0 });
			W.getWake().tmpI3.resize(n_CUDA_wake, { 0.0, 0.0 });
			W.getWake().tmpNei.resize(n_CUDA_wake, 0);

			//mesh.resize(n_CUDA_vel, { 0, 0 });

			W.getInfo('i') << "CUDA memory resize: " << curLength << " -> " << n_CUDA_wake << " vortices" << std::endl;
		}// if (W.getWake().vtx.size() > W.getWake().devNWake)

		//Обнуляем память, выделенную выше для хранения следа
		cuClearWakeMem(W.getWake().devNWake, W.getWake().devVtxPtr);

		//Копирование следа на видеокарту
		cuCopyWakeToDev(W.getWake().vtx.size(), W.getWake().vtx.data(), W.getWake().devVtxPtr);
	}



	if (W.getSource().vtx.size() > 0)
	{

		//Если зарезервировано меньше, чем источников в пелене
		if (W.getSource().vtx.size() > W.getSource().devNWake)
		{
			size_t curLength = W.getSource().devNWake;
			
			//Освобождаем всю память на ведеокарте
			if (curLength > 0)
				cuDeleteFromDev(W.getSource().devVtxPtr);
			
			size_t sz = curLength;
			while (W.getSource().vtx.size() > sz)
				sz += INC_VORT_DEV;

			//Резервируем новое количество памяти
			W.getSource().devVtxPtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(sz, W.getSource().devNWake);


			W.getInfo('i') << "CUDA memory resize: " << curLength << " -> " << sz << " sources" << std::endl;
		}// if (W.getSource().vtx.size() > W.getSource().devNWake)


		 //Обнуляем память, выделенную выше для хранения следа
		cuClearWakeMem(W.getSource().devNWake, W.getSource().devVtxPtr);

		//Копирование следа на видеокарту
		cuCopyWakeToDev(W.getSource().vtx.size(), W.getSource().vtx.data(), W.getSource().devVtxPtr);
	}


}


//Обновление состояния всех профилей и всего, что с ними тесно связано
void Gpu::RefreshAfls() 
{
	const int& id = W.getParallel().myidWork;
	if (W.getNumberOfBoundary() > 0)
	{
		if (W.getNumberOfBoundary() > n_CUDA_afls)
		{
			//Обнуление массива на хосте, содержащего количество вихрей в виртуальном следе на профилях
			n_CUDA_virtWake.resize(W.getNumberOfBoundary(), 0);
			
			n_CUDA_r.resize(W.getNumberOfBoundary(), 0);
			
			n_CUDA_panel.resize(W.getNumberOfBoundary(), 0);

			for (size_t s = 0; s < n_CUDA_afls; ++s)
			{
				//Чистим все массивы
				cuDeleteFromDev(W.getBoundary(s).virtualWake.devVtxPtr);
				cuDeleteFromDev(W.getBoundary(s).virtualWake.devVelPtr);
				cuDeleteFromDev(W.getBoundary(s).virtualWake.devRadPtr);
				cuDeleteFromDev(W.getBoundary(s).virtualWake.devI0Ptr);
				cuDeleteFromDev(W.getBoundary(s).virtualWake.devI1Ptr);
				cuDeleteFromDev(W.getBoundary(s).virtualWake.devI2Ptr);
				cuDeleteFromDev(W.getBoundary(s).virtualWake.devI3Ptr);
				
				cuDeleteFromDev(W.getBoundary(s).afl.devRPtr);
				cuDeleteFromDev(W.getBoundary(s).afl.devRhsPtr);

				cuDeleteFromDev(W.getBoundary(s).afl.devFreeVortexSheetPtr);
				cuDeleteFromDev(W.getBoundary(s).afl.devAttachedVortexSheetPtr);
				cuDeleteFromDev(W.getBoundary(s).afl.devAttachedSourceSheetPtr);
			}

			if (n_CUDA_afls)
			{
				cuDeleteFromDev(dev_ptr_nPanels);
				cuDeleteFromDev(dev_ptr_nVortices);

				cuDeleteFromDev(dev_ptr_ptr_vtx);
				cuDeleteFromDev(dev_ptr_ptr_vel);
				cuDeleteFromDev(dev_ptr_ptr_rad);
				cuDeleteFromDev(dev_ptr_ptr_i0);
				cuDeleteFromDev(dev_ptr_ptr_i1);
				cuDeleteFromDev(dev_ptr_ptr_i2);
				cuDeleteFromDev(dev_ptr_ptr_i3);

				cuDeleteFromDev(dev_ptr_ptr_r);
				cuDeleteFromDev(dev_ptr_ptr_rhs);

				cuDeleteFromDev(dev_ptr_ptr_freeVortexSheet);
				cuDeleteFromDev(dev_ptr_ptr_attachedVortexSheet);
				cuDeleteFromDev(dev_ptr_ptr_attachedSourceSheet);
			}

			//Временные массивы для агрегации числа панелей и числа виртуальных вихрей (округленные вверх до блока) на профилях для их последующей отправки в dev_ptr_nPanels и dev_ptr_nVortices
			std::vector<size_t> host_nPanels(0);
			std::vector<size_t> host_nVortices(0);

			n_CUDA_afls = 0;

			for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
			{
				const size_t& szPan = W.getBoundary(s).afl.np;
				const size_t& szVtx = std::max(W.getBoundary(s).virtualWake.vtx.size(),   W.getPassport().wakeDiscretizationProperties.vortexPerPanel * W.getBoundary(s).afl.np);
				
				//W.getInfo('t') << "szVtx = " << szVtx << std::endl;

				W.getBoundary(s).virtualWake.devVtxPtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(szVtx, n_CUDA_virtWake[s]);
				W.getBoundary(s).virtualWake.devVelPtr = ReserveDevMem<double, 2>(szVtx, n_CUDA_virtWake[s]);
				W.getBoundary(s).virtualWake.devRadPtr = ReserveDevMem<double, 1>(szVtx, n_CUDA_virtWake[s]);
				W.getBoundary(s).virtualWake.devI0Ptr = ReserveDevMem<double, 1>(szVtx, n_CUDA_virtWake[s]);
				W.getBoundary(s).virtualWake.devI1Ptr = ReserveDevMem<double, 1>(szVtx, n_CUDA_virtWake[s]);
				W.getBoundary(s).virtualWake.devI2Ptr = ReserveDevMem<double, 2>(szVtx, n_CUDA_virtWake[s]);
				W.getBoundary(s).virtualWake.devI3Ptr = ReserveDevMem<double, 2>(szVtx, n_CUDA_virtWake[s]);

				W.getBoundary(s).virtualWake.devNWake = n_CUDA_virtWake[s];


				W.getBoundary(s).afl.devRPtr = ReserveDevMem<double, 2>(szPan + 1, n_CUDA_r[s]);
				W.getBoundary(s).afl.devRhsPtr = ReserveDevMem<double, 1>(szPan, n_CUDA_panel[s]);

				W.getBoundary(s).afl.devFreeVortexSheetPtr = ReserveDevMem<double, 1>(szPan, n_CUDA_panel[s]);
				W.getBoundary(s).afl.devAttachedVortexSheetPtr = ReserveDevMem<double, 1>(szPan, n_CUDA_panel[s]);
				W.getBoundary(s).afl.devAttachedSourceSheetPtr = ReserveDevMem<double, 1>(szPan, n_CUDA_panel[s]);


				n_CUDA_afls++;
				host_nVortices.push_back(n_CUDA_virtWake[s]);
				host_nPanels.push_back(n_CUDA_panel[s]);

				W.getBoundary(s).virtualWake.tmpVel.resize(n_CUDA_virtWake[s], { 0.0, 0.0 });
				W.getBoundary(s).virtualWake.tmpRad.resize(n_CUDA_virtWake[s], 1.0e+5);

				W.getBoundary(s).virtualWake.tmpI0.resize(n_CUDA_virtWake[s], 0.0);
				W.getBoundary(s).virtualWake.tmpI1.resize(n_CUDA_virtWake[s], 0.0);
				W.getBoundary(s).virtualWake.tmpI2.resize(n_CUDA_virtWake[s], { 0.0, 0.0 });
				W.getBoundary(s).virtualWake.tmpI3.resize(n_CUDA_virtWake[s], { 0.0, 0.0 });

				W.getBoundary(s).afl.tmpRhs.resize(n_CUDA_panel[s], 0.0);
			}// for s

			dev_ptr_nPanels = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_nPanels.data());
			dev_ptr_nVortices = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_nVortices.data());

			//Временные массивы для агрегации указателей на видеокарте
			std::vector<double*> host_ptr_vtx;
			std::vector<double*> host_ptr_vel;
			std::vector<double*> host_ptr_rad;
			std::vector<double*> host_ptr_i0;
			std::vector<double*> host_ptr_i1;
			std::vector<double*> host_ptr_i2;
			std::vector<double*> host_ptr_i3;

			std::vector<double*> host_ptr_r;
			std::vector<double*> host_ptr_rhs;

			std::vector<double*> host_ptr_freeVortexSheet;
			std::vector<double*> host_ptr_attachedVortexSheet;
			std::vector<double*> host_ptr_attachedSourceSheet;

			for (size_t q = 0; q < W.getNumberOfBoundary(); ++q)
			{
				host_ptr_vtx.push_back(W.getBoundary(q).virtualWake.devVtxPtr);
				host_ptr_vel.push_back(W.getBoundary(q).virtualWake.devVelPtr);
				host_ptr_rad.push_back(W.getBoundary(q).virtualWake.devRadPtr);
				host_ptr_i0.push_back(W.getBoundary(q).virtualWake.devI0Ptr);
				host_ptr_i1.push_back(W.getBoundary(q).virtualWake.devI1Ptr);
				host_ptr_i2.push_back(W.getBoundary(q).virtualWake.devI2Ptr);
				host_ptr_i3.push_back(W.getBoundary(q).virtualWake.devI3Ptr);

				host_ptr_r.push_back(W.getBoundary(q).afl.devRPtr);
				host_ptr_rhs.push_back(W.getBoundary(q).afl.devRhsPtr);

				host_ptr_freeVortexSheet.push_back(W.getBoundary(q).afl.devFreeVortexSheetPtr);
				host_ptr_attachedVortexSheet.push_back(W.getBoundary(q).afl.devAttachedVortexSheetPtr);
				host_ptr_attachedSourceSheet.push_back(W.getBoundary(q).afl.devAttachedSourceSheetPtr);

			}

			dev_ptr_ptr_vtx = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_vtx.data());
			dev_ptr_ptr_rad = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_rad.data());
			dev_ptr_ptr_vel = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_vel.data());
			dev_ptr_ptr_i0 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_i0.data());
			dev_ptr_ptr_i1 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_i1.data());
			dev_ptr_ptr_i2 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_i2.data());
			dev_ptr_ptr_i3 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_i3.data());

			dev_ptr_ptr_r = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_r.data());
			dev_ptr_ptr_rhs = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_rhs.data());

			dev_ptr_ptr_freeVortexSheet = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_freeVortexSheet.data());
			dev_ptr_ptr_attachedVortexSheet = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_attachedVortexSheet.data());
			dev_ptr_ptr_attachedSourceSheet = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_attachedSourceSheet.data());
		}//if (W.getNumberOfBoundary() > n_CUDA_afls)
		

		for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
		{
			
			/*
			//Обнуляем память, выделенную выше для хранения виртуального следа
			cuClearWakeMem(n_CUDA_virtWake[s], W.getBoundary(s).virtualWake.devVtxPtr);

			//Копирование виртуального следа на видеокарту
			cuCopyWakeToDev(W.getBoundary(s).virtualWake.vtx.size(), W.getBoundary(s).virtualWake.vtx.data(), W.getBoundary(s).virtualWake.devVtxPtr);
			*/

			//Копирование вершин профиля на видеокарту
			cuCopyRsToDev(W.getBoundary(s).afl.r.size(), W.getBoundary(s).afl.r.data(), W.getBoundary(s).afl.devRPtr);

			//Копирование слоев на видеокарту			
			std::vector<double> host_freeVortexSheet(W.getBoundary(s).afl.np);
			std::vector<double> host_attachedVortexSheet(W.getBoundary(s).afl.np);
			std::vector<double> host_attachedSourceSheet(W.getBoundary(s).afl.np);

			const Sheet& sh = W.getBoundary(s).sheets;
						
			for (size_t p = 0; p < W.getBoundary(s).afl.np; ++p)
			{				
				host_attachedVortexSheet[p] = sh.attachedVortexSheet[p][0];
				host_attachedSourceSheet[p] = sh.attachedSourceSheet[p][0];
				
				if (W.getParallel().myidWork == 0)
					host_freeVortexSheet[p] = sh.freeVortexSheet[p][0];
			}//for p
			
			MPI_Bcast(host_freeVortexSheet.data(), (int)host_freeVortexSheet.size(), MPI_DOUBLE, 0, W.getParallel().commWork);
			
			cuCopyFixedArray(W.getBoundary(s).afl.devFreeVortexSheetPtr, host_freeVortexSheet.data(), sizeof(double)*host_freeVortexSheet.size());
			cuCopyFixedArray(W.getBoundary(s).afl.devAttachedVortexSheetPtr, host_attachedVortexSheet.data(), sizeof(double)*host_attachedVortexSheet.size());
			cuCopyFixedArray(W.getBoundary(s).afl.devAttachedSourceSheetPtr, host_attachedSourceSheet.data(), sizeof(double)*host_attachedSourceSheet.size());
		}
	}
}//RefreshAfls()


//Обновление состояния всех профилей и всего, что с ними тесно связано
void Gpu::RefreshVirtualWakes()
{
	const int& id = W.getParallel().myidWork;
	for(size_t s = 0; s < W.getNumberOfBoundary(); ++s)
	{
		if (n_CUDA_virtWake[s] < W.getBoundary(s).virtualWake.vtx.size())
		{			
			//std::cout << "s = " << s << "n_CUDA_virtWake[s] = " << n_CUDA_virtWake[s] << "W.getBoundary(s).virtualWake.vtx.size() = " << W.getBoundary(s).virtualWake.vtx.size() << std::endl;
			
			//Чистим все массивы
			cuDeleteFromDev(W.getBoundary(s).virtualWake.devVtxPtr);
			cuDeleteFromDev(W.getBoundary(s).virtualWake.devVelPtr);
			cuDeleteFromDev(W.getBoundary(s).virtualWake.devRadPtr);
			cuDeleteFromDev(W.getBoundary(s).virtualWake.devI0Ptr);
			cuDeleteFromDev(W.getBoundary(s).virtualWake.devI1Ptr);
			cuDeleteFromDev(W.getBoundary(s).virtualWake.devI2Ptr);
			cuDeleteFromDev(W.getBoundary(s).virtualWake.devI3Ptr);

			const size_t& szPan = W.getBoundary(s).afl.np;
			const size_t& szVtx = std::max(W.getBoundary(s).virtualWake.vtx.size(), W.getPassport().wakeDiscretizationProperties.vortexPerPanel * W.getBoundary(s).afl.np);

			//W.getInfo('t') << "szVtx = " << szVtx << std::endl;

			W.getBoundary(s).virtualWake.devVtxPtr = ReserveDevMem<double, sizeof(Vortex2D) / sizeof(double)>(szVtx, n_CUDA_virtWake[s]);
			W.getBoundary(s).virtualWake.devVelPtr = ReserveDevMem<double, 2>(szVtx, n_CUDA_virtWake[s]);
			W.getBoundary(s).virtualWake.devRadPtr = ReserveDevMem<double, 1>(szVtx, n_CUDA_virtWake[s]);
			W.getBoundary(s).virtualWake.devI0Ptr = ReserveDevMem<double, 1>(szVtx, n_CUDA_virtWake[s]);
			W.getBoundary(s).virtualWake.devI1Ptr = ReserveDevMem<double, 1>(szVtx, n_CUDA_virtWake[s]);
			W.getBoundary(s).virtualWake.devI2Ptr = ReserveDevMem<double, 2>(szVtx, n_CUDA_virtWake[s]);
			W.getBoundary(s).virtualWake.devI3Ptr = ReserveDevMem<double, 2>(szVtx, n_CUDA_virtWake[s]);

			W.getInfo('i') << "CUDA virtual wake["<< s << "] resize: " << W.getBoundary(s).virtualWake.devNWake << " -> " << n_CUDA_virtWake[s] << " vortices" << std::endl;
			W.getBoundary(s).virtualWake.devNWake = n_CUDA_virtWake[s];
			

			W.getBoundary(s).virtualWake.tmpVel.resize(n_CUDA_virtWake[s], { 0.0, 0.0 });
			W.getBoundary(s).virtualWake.tmpRad.resize(n_CUDA_virtWake[s], 1.0e+5);

			W.getBoundary(s).virtualWake.tmpI0.resize(n_CUDA_virtWake[s], 0.0);
			W.getBoundary(s).virtualWake.tmpI1.resize(n_CUDA_virtWake[s], 0.0);
			W.getBoundary(s).virtualWake.tmpI2.resize(n_CUDA_virtWake[s], { 0.0, 0.0 });
			W.getBoundary(s).virtualWake.tmpI3.resize(n_CUDA_virtWake[s], { 0.0, 0.0 });
			
			std::vector<size_t> host_nVortices(0);
			host_nVortices.reserve(W.getNumberOfBoundary());
			for (size_t q = 0; q < W.getNumberOfBoundary(); ++q)
				host_nVortices.push_back(W.getBoundary(s).virtualWake.devNWake);
					   			
			dev_ptr_nVortices = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_nVortices.data());

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

			dev_ptr_ptr_vtx = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_vtx.data());
			dev_ptr_ptr_rad = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_rad.data());
			dev_ptr_ptr_vel = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_vel.data());
			dev_ptr_ptr_i0 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_i0.data());
			dev_ptr_ptr_i1 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_i1.data());
			dev_ptr_ptr_i2 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_i2.data());
			dev_ptr_ptr_i3 = ReserveDevMemAndCopyFixedArray(W.getNumberOfBoundary(), host_ptr_i3.data());			
		}//if (n_CUDA_virtWake[s] < W.getBoundary(s).virtualWake.vtx.size())


		for (size_t s = 0; s < W.getNumberOfBoundary(); ++s)
		{
			//Обнуляем память, выделенную выше для хранения виртуального следа
			cuClearWakeMem(n_CUDA_virtWake[s], W.getBoundary(s).virtualWake.devVtxPtr);

			//Копирование виртуального следа на видеокарту
			cuCopyWakeToDev(W.getBoundary(s).virtualWake.vtx.size(), W.getBoundary(s).virtualWake.vtx.data(), W.getBoundary(s).virtualWake.devVtxPtr);			
		}
	}
}



#endif