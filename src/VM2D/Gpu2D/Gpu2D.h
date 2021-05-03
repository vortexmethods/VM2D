/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
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
\brief Заголовочный файл с описанием класса Gpu
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
*/


#ifndef GPU_H
#define GPU_H

#include <limits>
#include <memory>

#include "cuLib2D.cuh"
#include "Gpudefs.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс, обеспечивающий возможность выполнения вычислений на GPU по технологии Nvidia CUDA

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.10
	\date 17 мая 2021 г.
	*/
	class Gpu
	{
	private:
		/// Константная ссылка на решаемую задачу
		const World2D& W;

	public:

		//static int nReserve; //Для контроля паритета выделения и освобождения памяти

#if defined(__CUDACC__) || defined(USE_CUDA)

		/// \brief Освобождение видеопамяти на графической карте
		/// \todo Откомментировать
		template<typename T>
		void ReleaseDevMem(T* ptr, int code)
		{
			cuDeleteFromDev(ptr, code);
			//--nReserve;
		}


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
			//++nReserve;

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
			//++nReserve;

			cuCopyFixedArray(dev_ptr, host_src, sizeof(T) * n);

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
		void CopyMemFromDev(size_t n, T* dev_ptr, T* host_ptr, int code = 0) const
		{
			cuCopyMemFromDev((void*)host_ptr, (void*)dev_ptr, sizeof(T) * n * dim, code);
		};//CopyMemFromDev(...)

		/// \brief Копирование данных в видеопамять на графической карте с хоста
		///
		/// \tparam T тип данных
		/// \tparam dim внутренняя размерность массива
		/// \param[in] n длина массива (внешняя размерность), который требуется скопировать
		/// \param[in] dev_ptr адрес на графической карте - указатель на начало массива
		/// \param[in] host_ptr адрес на хосте, куда требуется скопировать массив
		template<typename T, size_t dim>
		void CopyMemToDev(size_t n, T* host_ptr, T* dev_ptr) const
		{
			cuCopyFixedArray((void*)dev_ptr, (void*)host_ptr, sizeof(T) * n * dim);
		};//CopyMemToDev(...)


		/// Синхронизация состояния следа с графической картой
		void RefreshWake(int code = 0);
		
		/// Синхронизация состояния профилей с графической картой
		void RefreshAfls(int code = 0);
		void RefreshVirtualWakes(int code = 0);

		/// Обновление состояния сетки для вычисления VP
		void RefreshVP(int code = 0);


		// Ниже - данные для вычисления скоростей	

		/// Число профилей в задаче
		size_t n_CUDA_afls;

		/// Указатели на массивы, которые хранятся на видеокарте и содержат число панелей и число виртуальных вихрей на всех профилях
		size_t* dev_ptr_nPanels;
		size_t* dev_ptr_nVortices;


		//Переменная, которая лежит на хосте и хранит адрес на видеокарте массива, в котором хранятся указатели на соответствующие массивы
		double** dev_ptr_ptr_vtx;
		double** dev_ptr_ptr_vel;
		double** dev_ptr_ptr_rad;
		double** dev_ptr_ptr_i0;
		double** dev_ptr_ptr_i1;
		double** dev_ptr_ptr_i2;
		double** dev_ptr_ptr_i3;

		double** dev_ptr_ptr_r;
		double** dev_ptr_ptr_rhs;

		double** dev_ptr_ptr_freeVortexSheet;
		double** dev_ptr_ptr_attachedVortexSheet;
		double** dev_ptr_ptr_attachedSourceSheet;

		double** dev_ptr_ptr_meanEpsOverPanel;

		double** dev_ptr_ptr_viscousStresses;

		/// Массив на хосте, содержащий число виртуальных вихрей на профилях; его длина равна числу профилей
		std::vector<size_t> n_CUDA_virtWake;
		size_t n_CUDA_totalVirtWake;

		/// Массив на хосте, содержащий число панелей на профилях; его длина равна числу профилей
		std::vector<size_t> n_CUDA_panel;

		/// Длина массивов на видеокарте, зарезервированных на хранение сведений о следе (wake)
		size_t n_CUDA_wake;

		/// Длина массивов на видеокарте, зарезервированных на хранение сведений о следе (source)
		size_t n_CUDA_source;

		/// Длина массивов на видеокарте, зарезервированных на хранение сведений о точках вычисления VP
		size_t n_CUDA_velVP;


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


		/// \brief Установка переключателя расчетных схем 
		///
		/// \param[in] schemeSwitcher_ тип схемы
		/// - schemeSwitcher = 1 -- кусочно-постоянная схема
		/// - schemeSwitcher = 2 -- кусочно-линейная схема
		/// - schemeSwitcher = 11 -- схема типа МДВ
		void setSchemeSwitcher(int schemeSwitcher_)
		{
#if defined(__CUDACC__) || defined(USE_CUDA)		
			cuSetSchemeSwitcher(schemeSwitcher_, 1);
#endif
		}
	};

}//namespace VM2D




namespace VM2D 
{
	template <class T>
	class MyAlloc {
	public:
		// type definitions
		typedef T              value_type;
		typedef T*             pointer;
		typedef const T*       const_pointer;
		typedef T&             reference;
		typedef const T&       const_reference;
		typedef std::size_t    size_type;
		typedef std::ptrdiff_t difference_type;

		// rebind allocator to type U
		template <class U>
		struct rebind {
			typedef MyAlloc<U> other;
		};

		// return address of values
		pointer address(reference value) const {
			return &value;
		}
		const_pointer address(const_reference value) const {
			return &value;
		}

		/* constructors and destructor
		 * - nothing to do because the allocator has no state
		 */
		MyAlloc() throw() {
		}
		MyAlloc(const MyAlloc&) throw() {
		}
		template <class U>
		MyAlloc(const MyAlloc<U>&) throw() {
		}
		~MyAlloc() throw() {
		}

		// return maximum number of elements that can be allocated
		size_type max_size() const throw() {
			return std::numeric_limits<std::size_t>::max() / sizeof(T);
		}

		// allocate but don't initialize num elements of type T
		pointer allocate(size_type num, const void* = 0) {
			// print message and allocate memory with global new
			//std::cerr << "allocate " << num << " element(s)" 	<< " of size " << sizeof(T) << std::endl;
			
			pointer ret = (pointer)(::operator new(num * sizeof(T)));

			//pointer ret;
			//cudaHostAlloc((void**)&ret, num * sizeof(T), cudaHostAllocDefault);
			cuAlloc((void**)&ret, num * sizeof(T));

			//std::cerr << " allocated at: " << (void*)ret << std::endl;
			return ret;
		}

		// initialize elements of allocated storage p with value value
		void construct(pointer p, const T& value) {
			// initialize memory with placement new
			new((void*)p)T(value);
			//std::cerr << " construct " << std::endl;
		}

		// destroy elements of initialized storage p
		void destroy(pointer p) {
			// destroy objects by calling their destructor
			p->~T();
			//std::cerr << " destroy " << std::endl;
		}

		// deallocate storage p of deleted elements
		void deallocate(pointer p, size_type num) {
			// print message and deallocate memory with global delete
			
			//std::cerr << "deallocate " << num << " element(s)" << " of size " << sizeof(T) << " at: " << (void*)p << std::endl;
			
			//::operator delete((void*)p);
			cuDalloc((void*)p);
		}
	};

	// return that all specializations of this allocator are interchangeable
	template <class T1, class T2>
	bool operator== (const MyAlloc<T1>&,
		const MyAlloc<T2>&) throw() {
		return true;
	}
	template <class T1, class T2>
	bool operator!= (const MyAlloc<T1>&,
		const MyAlloc<T2>&) throw() {
		return false;
	}


	/*{
		// create a vector, using MyAlloc<> as allocator
		std::vector<int, VM2D::MyAlloc<int> > v(5, 37);

		v.resize(0);

		// insert elements
		// - causes reallocations
		v.push_back(42);
		v.push_back(56);

		v.reserve(10);

		v.push_back(11);
		v.push_back(22);
		v.push_back(33);
		v.push_back(44);
	}
	*/

}
















#endif