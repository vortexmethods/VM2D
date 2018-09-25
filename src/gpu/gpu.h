/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.3    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/09/26     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: gpu.h                                                            |
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
\brief Заголовочный файл с описанием класса gpu
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/


#ifndef CUDA_CUH
#define CUDA_CUH

#include <memory>

#include "cuLib.cuh"
#include "gpudefs.h"

class World2D;

/*!
\brief Класс, обеспечивающий возможность выполнения вычислений на GPU по технологии Nvidia CUDA

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/
class gpu
{
private:
	/// Константная ссылка на решаемую задачу
	const World2D& W;

public:
	
#if defined(__CUDACC__) || defined(USE_CUDA)
	
	/// \brief Резервирование видеопамяти на графической карте
	///
	/// Резервирует на видеокарте память под хранение массива заданной размерности
	/// \n Производится округление длины массива вверх до длины, кратной размеру блока CUBLOCK
	///
	/// \tparam T тип данных
	/// \tparam dim внутренняя размерность массива
	/// \param[in] n длина массива (внешняя размерность), для которого требуется зарезервировать память
	/// \param[out] new_n ссылка на длину (внешнюю размерность) массива, под который выделена память
	/// \return адрес на графической карте - указатель на начало массива
	template<typename T, size_t dim>
	T* ReserveDevMem(size_t n, size_t& new_n) 
	{
		size_t nBlocks = n / CUBLOCK;
		if (n % CUBLOCK)
			nBlocks++;

		new_n = nBlocks*CUBLOCK;

		void* ptr;

		cuReserveDevMem(ptr, CUBLOCK * nBlocks * dim * sizeof(T));

		return (T*)ptr;
	} //ReserveDevMem(...)

	
	/// \brief Резервирование видеопамяти на графической карте и копирование туда данных
	///
	/// Резервирует на видеокарте память под хранение массива фиксированной длины и копирует туда этот массив
	/// \n Не производится округление длины массива
	///
	/// \tparam T тип данных
	/// \param[in] n длина массива, для которого требуется зарезервировать память
	/// \param[in] host_src указатель на начало массива в памяти хоста
	/// \return адрес на графической карте - указатель на начало массива, скопированного с хоста
	template<typename T>
	T* ReserveDevMemAndCopyFixedArray(size_t n, T* host_src)
	{
		void* dev_ptr;

		cuReserveDevMem(dev_ptr, n * sizeof(T));
		cuCopyFixedArray(dev_ptr, host_src, sizeof(size_t) * n);	

		return (T*)dev_ptr;
	}//ReserveDevMemAndCopyFixedArray(...)

	
	/// \brief Копирование данных из видеопамяти на графической карте на хост
	///
	/// \tparam T тип данных
	/// \tparam dim внутренняя размерность массива
	/// \param[in] n длина массива (внешняя размерность), который требуется скопировать
	/// \param[in] dev_ptr адрес на графической карте - указатель на начало массива
	/// \param[in] host_ptr адрес на хосте, куда требуется скопировать массив
	template<typename T, size_t dim>
	void CopyMemFromDev(size_t n, T* dev_ptr, T* host_ptr) const
	{
		cuCopyMemFromDev((void*)host_ptr, (void*)dev_ptr, sizeof(T) * n * dim);
	};//CopyMemFromDev(...)


	/// Синхронизация состояния следа с графической картой
	void RefreshWake();
	
	
	/// Синхронизация состояния профилей с графической картой
	void RefreshAfls();
	
	// Данные для вычисления скоростей	
	size_t n_CUDA_afls;

	size_t* dev_nPanels;
	size_t* dev_nVortices;
	

	double** dev_ptr_ptr_vtx;
	double** dev_ptr_ptr_i0;
	double** dev_ptr_ptr_i1;
	double** dev_ptr_ptr_i2;
	double** dev_ptr_ptr_i3;
	double** dev_ptr_ptr_rad;
	double** dev_ptr_ptr_vel;
	
	double** dev_ptr_ptr_r;
	double** dev_ptr_ptr_rhs;

	std::vector<size_t> n_CUDA_vtxs;
	std::vector<size_t> n_CUDA_i0s;
	std::vector<size_t> n_CUDA_i1s;
	std::vector<size_t> n_CUDA_i2s;
	std::vector<size_t> n_CUDA_i3s;
	std::vector<size_t> n_CUDA_rads;
	std::vector<size_t> n_CUDA_vels;

	std::vector<size_t> n_CUDA_rs;
	std::vector<size_t> n_CUDA_rhss;

	size_t n_CUDA_vel;
	size_t n_CUDA_rad;
	size_t n_CUDA_i0;
	size_t n_CUDA_i1;
	size_t n_CUDA_i2;
	size_t n_CUDA_i3;
		
	size_t n_CUDA_mesh;
	size_t n_CUDA_nei;	

#endif

	/// \brief Конструктор
	/// 
	/// \param[in] W_ константная ссылка на решаемую задачу	
	gpu(const World2D& W_);

	//// \brief Деструктор
	~gpu();
};

#endif