/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
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
| (at your option) any later version.	                                      |
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
\version 1.1
\date 2 апреля 2018 г.
*/


#ifndef CUDA_CUH
#define CUDA_CUH

#include <memory>
#include <vector>


#include "defs.h"
#include "Vortex2D.h"
#include "Parallel.h"
#include "Boundary.h"
#include "cuLib.cuh"

/*!
\brief Класс, обеспечиваюий возможность выполнения вычислений на GPU по технологии Nvidia CUDA

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/
class gpu
{
private:
	/// Константная ссылка на список умных казателей на граничные условия
	const std::vector<std::unique_ptr<Boundary>>& boundary;

	/// Константная ссылка на вихревой след
	const Wake& wake;

	/// Константная ссылка на параметры исполнения в параллельном режиме
	const Parallel& parallel;
	
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

	
	/// \brief Резервирование видеопамяти на графической карте
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
	void CopyMemFromDev(size_t n, T* dev_ptr, T* host_ptr)
	{
		cuCopyMemFromDev((void*)host_ptr, (void*)dev_ptr, sizeof(T) * n * dim);
	};//CopyMemFromDev(...)


	/// Синхронизация состояния следа с гафической картой
	void RefreshWake();
	
	
	/// Синхронизация состояния профилей с гафической картой
	void RefreshAfls();
	
	/// Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
	void ExpCalcDiffVeloI1I2ToSetOfPoints(
		const size_t npt, double* dev_ptr_pt,
		double* dev_ptr_dr,
		const size_t nvt, double* dev_ptr_vt,
		std::vector<double>& I1, std::vector<Point2D>& I2,
		std::vector<double>& loci1, std::vector<Point2D>& loci2,
		double* dev_ptr_i1, double* dev_ptr_i2,
		bool useMesh);

	/// Вычисление влияния вихревого следа в заданном наборе точек - расчет правых частей системы
	void ExpGetWakeInfluence(
		const size_t npt, double* dev_ptr_pt, 
		const size_t nvt, double* dev_ptr_vt, 
		std::vector<double>& rhs, std::vector<double>& locrhs, double* dev_ptr_rhs);

	/// Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
	void ExpGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(
		const size_t npt, double* dev_ptr_pt, double* dev_ptr_rad,
		const size_t nr, double* dev_ptr_r,
		std::vector<double>& I0, std::vector<Point2D>& I3,
		std::vector<double>& loci0, std::vector<Point2D>& loci3,
		double* dev_ptr_i0, double* dev_ptr_i3);

	/// Вычисление конвективных скоростей и радиусов вихревых доменов в заданном наборе точек
	void ExpCalcConvVeloToSetOfPoints(
		const size_t npt, double* dev_ptr_pt,
		const size_t nvt, double* dev_ptr_vt,
		const size_t nbou, size_t* dev_nPanels, double** dev_ptr_ptr_pnl,
		std::vector<Point2D>& Vel, std::vector<double>& Rad,
		std::vector<Point2D>& locvel, std::vector<double>& locrad,
		double* dev_ptr_vel, double* dev_ptr_rad,
		double minRad, double eps2);

	/// Вычисление конвективных скоростей вихревых доменов в заданном наборе точек		
	void ExpGetConvVelocityToSetOfPointsFromVirtualVortexes(
		const size_t npt, double* dev_ptr_pt,
		const size_t nvt, double* dev_ptr_vt,
		std::vector<Point2D>& Vel, std::vector<Point2D>& locvel, double* dev_ptr_vel,
		double eps2);
	
	/// Поиск соседних вихревых элементов
	void ExpGetPairs(int type, std::vector<int>& NEIB);
	
	// Данные для вычисления скоростей
	std::vector<double*> host_ptr_ptr_pnl;
	std::vector<double*> host_ptr_ptr_r;
	std::vector<double*> host_ptr_ptr_i0;
	std::vector<double*> host_ptr_ptr_i1;
	std::vector<double*> host_ptr_ptr_i2;
	std::vector<double*> host_ptr_ptr_i3;
	std::vector<double*> host_ptr_ptr_rad;
	std::vector<double*> host_ptr_ptr_vel;
	std::vector<double*> host_ptr_ptr_rhs;

	size_t* dev_nPanels;
	double** dev_ptr_ptr_pnl;
	double** dev_ptr_ptr_r;
	double** dev_ptr_ptr_i0;
	double** dev_ptr_ptr_i1;
	double** dev_ptr_ptr_i2;
	double** dev_ptr_ptr_i3;
	double** dev_ptr_ptr_rad;
	double** dev_ptr_ptr_vel;
	double** dev_ptr_ptr_rhs;

	//double* dev_ptr_wake;
	double* dev_ptr_vel;
	double* dev_ptr_rad;
	double* dev_ptr_i0;
	double* dev_ptr_i1;
	double* dev_ptr_i2;
	double* dev_ptr_i3;

	size_t n_CUDA_afls;
	std::vector<size_t> n_CUDA_pnls;
	std::vector<size_t> n_CUDA_rs;
	std::vector<size_t> n_CUDA_i0s;
	std::vector<size_t> n_CUDA_i1s;
	std::vector<size_t> n_CUDA_i2s;
	std::vector<size_t> n_CUDA_i3s;
	std::vector<size_t> n_CUDA_rads;
	std::vector<size_t> n_CUDA_vels;
	std::vector<size_t> n_CUDA_rhss;

	//size_t n_CUDA_wake;
	size_t n_CUDA_vel;
	size_t n_CUDA_rad;
	size_t n_CUDA_i0;
	size_t n_CUDA_i1;
	size_t n_CUDA_i2;
	size_t n_CUDA_i3;

	std::vector<Point2D> vels;
	std::vector<double> rads;
	std::vector<double> i0;
	std::vector<double> i1;
	std::vector<Point2D> i2;
	std::vector<Point2D> i3;
	std::vector<int> nei;

	std::vector<std::vector<Point2D>> virtvels;
	std::vector<std::vector<double>> virtrads;
	std::vector<std::vector<double>> virti0;
	std::vector<std::vector<double>> virti1;
	std::vector<std::vector<Point2D>> virti2;
	std::vector<std::vector<Point2D>> virti3;
	std::vector<std::vector<double>> virtrhs;

	//Данные для коллапса
	int* dev_ptr_mesh;
	int* dev_ptr_nei;

	size_t n_CUDA_mesh;
	size_t n_CUDA_nei;

	//std::vector<numvector<int, 2>> mesh;

#endif

	/// \brief Конструктор
	/// 
	/// \param[in] boundary_ константная ссылка на вектор из указателей на все граничные условия
	/// \param[in] wake_ константная ссылка на вихревой след
	/// \param[in] parallel_ константная ссылка на параметры параллельного исполнения
	gpu(const std::vector<std::unique_ptr<Boundary>>& boundary_, const Wake& wake_, const Parallel& parallel_);

	//// \brief Деструктор
	~gpu();
};

#endif