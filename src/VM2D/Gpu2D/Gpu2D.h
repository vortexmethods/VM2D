/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.6    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/10/28     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Gpu2D.h                                                          |
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
\brief Заголовочный файл с описанием класса Gpu
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.6   
\date 28 октября 2019 г.
*/


#ifndef GPU_H
#define GPU_H

#include <memory>

#include "cuLib2D.cuh"
#include "Gpudefs.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс, обеспечивающий возможность выполнения вычислений на GPU по технологии Nvidia CUDA

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.6
	\date 28 октября 2019 г.
	*/
	class Gpu
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
		/// \param[out] new_n ссылка на длину (внешнюю размерность) массива, под который выделена память (округляется вверх до длины, кратной размеру блока CUBLOCK)
		/// \return адрес на графической карте - указатель на начало массива
		template<typename T, size_t dim>
		T* ReserveDevMem(size_t n, size_t& new_n)
		{
			size_t nBlocks = n / CUBLOCK;
			if (n % CUBLOCK)
				nBlocks++;

			new_n = nBlocks * CUBLOCK;

			void* ptr;

			cuReserveDevMem(ptr, new_n * dim * sizeof(T));

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
		void RefreshVirtualWakes();

		// Ниже - данные для вычисления скоростей	

		/// Число профилей в задаче
		size_t n_CUDA_afls;

		/// Указатели на массивы, которые хранятся на видеокарте и содержат число панелей и виртуальных вихрей на всех профилях
		size_t* dev_ptr_nPanels;
		size_t* dev_ptr_nVortices;


		//Переменная, которая лежит на хосте и хранит адрес на видеокарте массива, в котором хранятся указатели на соответствующие массивы
		double** dev_ptr_ptr_vtx;
		double** dev_ptr_ptr_i0;
		double** dev_ptr_ptr_i1;
		double** dev_ptr_ptr_i2;
		double** dev_ptr_ptr_i3;
		double** dev_ptr_ptr_rad;
		double** dev_ptr_ptr_vel;

		double** dev_ptr_ptr_r;
		double** dev_ptr_ptr_rhs;

		double** dev_ptr_ptr_freeVortexSheet;
		double** dev_ptr_ptr_attachedVortexSheet;
		double** dev_ptr_ptr_attachedSourceSheet;


		/// Массив на хосте, содержащий число виртуальных вихрей на профилях; его длина равна числу профилей
		std::vector<size_t> n_CUDA_virtWake;

		/// Массив на хосте, содержащий число вершин на профилях (на единицу больше, чем число панелей); его длина равна числу профилей
		std::vector<size_t> n_CUDA_r;

		/// Массив на хосте, содержащий число панелей на профилях; его длина равна числу профилей
		std::vector<size_t> n_CUDA_panel;

		/// Длина массивов на видеокарте, зарезервированных на хранение сведений о следе (wake)
		size_t n_CUDA_wake;

		/// Длина массивов на видеокарте, зарезервированных на хранение сведений о следе (source)
		size_t n_CUDA_source;


#endif

		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу	
		Gpu(const World2D& W_);

		//// \brief Деструктор
		~Gpu();


		/// \brief Установка коэффициента разгона потока
		///
		/// \param[in] cft_ множитель, соответствующий степени разгона потока
		void setAccelCoeff(double cft_)
		{
#if defined(__CUDACC__) || defined(USE_CUDA)		
			cuSetAccelCoeff(cft_);
#endif
		}


		/// \brief Установка правой границы самого правого профиля (для организации увеличения радиуса коллапса)
		///
		/// \param[in] pos_ абсцисса правой границы самого правого профиля
		/// \param[in] refLength_ характерная длина, на которой происходит увеличение радиуса коллапса
		void setCollapseCoeff(double pos_, double refLength_)
		{
#if defined(__CUDACC__) || defined(USE_CUDA)		
			cuSetCollapseCoeff(pos_, refLength_);
#endif
		}


		/// \brief Установка максимально допустимой циркуляции вихря
		///
		/// \param[in] gam_ максимально допустимая циркуляция вихря
		void setMaxGamma(double gam_)
		{
#if defined(__CUDACC__) || defined(USE_CUDA)		
			cuSetMaxGamma(gam_);
#endif
		}
	};

}//namespace VM2D

#endif