/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.5    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2024/06/19     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2024 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: defs.h                                                           |
| Info: Source code of BH                                                     |
|                                                                             |
| This file is part of BH.                                                    |
| BH is free software: you can redistribute it and/or modify it               |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| BHcu is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with BH.  If not, see <http://www.gnu.org/licenses/>.                 |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Вспомогательные функции
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 19 июня 2024 г.
*/

#pragma once

#include <iostream>
#include "PointsCopy.h"

#ifdef calcOp
	#define ADDOP(n) BH::op += ((int)n)
#else
	#define ADDOP(n) 
#endif

namespace BH
{
	extern long long op;

	//Мировые константы
	static const double PI = 3.1415926535897932384626;
	static const double DPI = 2.0 * 3.1415926535897932384626;
	static const double IPI = 1.0 / PI;
	static const double IDPI = 0.5 / PI;

	/// Длина мортоновского кода для каждой координаты (не более половины длины int в битах)
	static const int codeLength = 14;

	/// 2 в степени длины мортоновского кода (на каждую координату)
	static const int twoPowCodeLength = (1 << codeLength);

	
	/// Вспомогательная функция корректировки capacity вектора (при необходимости - удваивает)
	inline void SizeCheck(std::vector<Point2D>& i00)
	{
		if (i00.capacity() == i00.size())
			i00.reserve(i00.size() * 2);
	}


	/// Умножение комплексных чисел
	inline Point2D multz(const Point2D& a, const Point2D& b)
	{
		ADDOP(4);
		return Point2D({ a[0] * b[0] - a[1] * b[1], a[1] * b[0] + a[0] * b[1] });
	}

	/// Возведение в степень комплексных чисел
	inline Point2D powz(const Point2D& z, double n)
	{
		double phi, R;
		ADDOP(10);
		phi = n * atan2(z[1], z[0]);
		R = pow(z.length2(), 0.5*n);
		return Point2D({ R * cos(phi), R * sin(phi) });
	}

	/// Умножение a на комплексно сопряженноe к b
	inline Point2D multzA(const Point2D& a, const Point2D& b)
	{
		ADDOP(4);
		return Point2D({ a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1] });
	}

	/// Шаблонная функция возведения в квадрат
	template <typename T>
	inline T sqr(T x)
	{
		ADDOP(1);
		return x * x;
	}

	/// \brief Шаблонная функция знака числа
	/// Написана оптимизированная версия, которая работает самым быстрым возможным образом, 
	/// т.к. не содержит ни одного оператора условного перехода
	template <typename T>
	inline int sign(T val)
	{
		return (T(0) < val) - (val < T(0));
	}


	//Округление "в потолок" результата деления x на y
	//(годится для любых натуральных, но работает довольно медленно, здесь не нужна)
	//unsigned int ceil(unsigned int x, unsigned int y)
	//{
	//    return x / y + (x % y != 0);
	//}

	/// \brief Округление "в потолок" результата деления на степень двойки, эквивалент ceil(x / (2^p))
	///
	/// \param[in] x делимое
	/// \param[in] p показатель степени двойки в делителе
	inline int ceilpow2(unsigned int x, unsigned int p) 
	{
		return (x >> p) + !!(x & ((1 << p) - 1));
	}

	/// \brief Округление "в потолок" результата деления пополам, эквивалент ceil(x / 2)
	inline int ceilhalf(unsigned int x) 
	{
		return (x >> 1) + (x & 1);
	}

	/// \brief Шаблонная функция вычисления евклидовой нормы вектора или списка
	///
	/// \param[in] b константная ссылка на вектор или список
	template<typename T>
	inline double norm(const T& b)
	{
		double norm = 0;
#ifndef OLD_OMP
#pragma omp simd reduction(+:norm)
#endif
			for (size_t i = 0; i < b.size(); i++)
				norm += (b[i] * b[i]);
			ADDOP(2*b.size() + 1);
			return sqrt(norm);
	}

	/// Шаблонная функция сложения двух векторов
	template<typename T>
	inline std::vector<T> operator+(const std::vector<T>& x, const std::vector<T>& y)
	{
		std::vector<T> c(x);
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (size_t i = 0; i < x.size(); ++i)
			c[i] += y[i];
		return c;
	}


	/// Шаблонная функция прибавления к одному вектору другого
	template<typename T>
	inline std::vector<T>& operator+=(std::vector<T>& x, const std::vector<T>& y)
	{
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (size_t i = 0; i < x.size(); ++i)
			x[i] += y[i];
		return x;
	}

	/// Шаблонная функция вычитания векторов
	template<typename T>
	inline std::vector<T> operator-(const std::vector<T>& x, const std::vector<T>& y)
	{
		std::vector<T> c(x);
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (size_t i = 0; i < x.size(); ++i)
			c[i] -= y[i];
		return c;
	}

	/// Шаблонная функция вычитания из одного вектора другого
	template<typename T>
	inline std::vector<T>& operator-=(std::vector<T>& x, const std::vector<T>& y)
	{
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (size_t i = 0; i < x.size(); ++i)
			x[i] -= y[i];
		return x;
	}

	/// Шаблонная функция умножения числа на вектор
	template<typename T>
	inline std::vector<T> operator*(const T lambda, const std::vector<T>& x)
	{
		std::vector<T> c(x);
		c.resize(x.size());
	
#ifndef OLD_OMP
#pragma omp simd
#endif
		for (size_t i = 0; i < x.size(); ++i)
			c[i] *= lambda;
		ADDOP(x.size());
		return c;
	}


	/// Шаблонная функция вычисления скалярного произведения двух векторов
	template<typename T>
	inline T operator&(const std::vector<T>& x, const std::vector<T>& y)
	{
		T c = 0;
#ifndef OLD_OMP
#pragma omp simd reduction(+:c)
#endif
		for (size_t i = 0; i < x.size(); ++i)
			c += x[i] * y[i];
		ADDOP(x.size());
		return c;
	}

}//namespace BH
