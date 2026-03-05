/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: operations.cuh                                                   |
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
\brief Вспомогательные операции
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/


#ifndef OPERATIONS_H
#define OPERATIONS_H

#include <cuda.h>

#include "Gpudefs.h"


__device__ const double ifac[] = { 0.0, 1.0, 1.0, 1.0 / 2.0, 1.0 / 6.0, 1.0 / 24.0, 1.0 / 120.0, 1.0 / 720.0, 1.0 / 5040.0, 1.0 / 40320.0, 1.0 / 362880.0, 1.0 / 3628800.0, 1.0 / 39916800.0, 1.0 / 479001600.0, 1.0 / 6227020800.0, 1.0 / 87178291200.0 };



namespace BHcu
{

#ifndef WIN32
#define __forceinline inline
#endif

	__device__ __forceinline double2 multz(const double2& a, const double2& b) {
		return make_double2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
	}

	__device__ __forceinline double2 multz(double ax, double ay, const double2& b) {
		return make_double2(ax * b.x - ay * b.y, ax * b.y + ay * b.x);
	}

	__device__ __forceinline double2 multzA(const double2& a, const double2& b)	{		
		return make_double2(a.x * b.x + a.y * b.y, -a.x * b.y + a.y * b.x);
	}	
	
	__device__ __forceinline double2 multzA(const double2& a, double bx, double by) {
		return make_double2(a.x * bx + a.y * by, a.y * bx - a.x * by);
	}

	__device__ __forceinline double2 operator*(double a, const double2& b) {
		return make_double2(a * b.x, a * b.y);
	}

	__device__ __forceinline float2 operator*(float a, const float2& b) {
		return make_float2(a * b.x, a * b.y);
	}

	__device__ __forceinline double2 operator*(int a, const double2& b)	{
		return make_double2(a * b.x, a * b.y);
	}

	__device__ __forceinline float2 operator*(int a, const float2& b) {
		return make_float2(a * b.x, a * b.y);
	}

	__device__ __forceinline double2 operator*(const double2& b, double a) {
		return make_double2(a * b.x, a * b.y);
	}

	__device__ __forceinline float2 operator*(const float2& b, float a) {
		return make_float2(a * b.x, a * b.y);
	}

	__device__ __forceinline double2 operator*(const double2& b, int a) {
		return make_double2(a * b.x, a * b.y);
	}

	__device__ __forceinline float2 operator*(const float2& b, int a) {
		return make_float2(a * b.x, a * b.y);
	}

	__device__ __forceinline double2 operator/(const double2& b, double a) {
		return make_double2(b.x / a, b.y / a);		
	}

	__device__ __forceinline float2 operator/(const float2& b, float a) {
		return make_float2(b.x / a, b.y / a);
	}

	__device__ __forceinline double2& operator+=(double2& a, const double2& b) {
		a.x += b.x;
		a.y += b.y;
		return a;
	}

	__device__ __forceinline float2& operator+=(float2& a, const float2& b) {
		a.x += b.x;
		a.y += b.y;
		return a;
	}

	__device__ __forceinline double2& operator-=(double2& a, const double2& b) {
		a.x -= b.x;
		a.y -= b.y;
		return a;
	}

	__device__ __forceinline float2& operator-=(float2& a, const float2& b) {
		a.x -= b.x;
		a.y -= b.y;
		return a;
	}

	__device__ __forceinline float2& operator*=(float2& b, float a) {
		b.x *= a;
		b.y *= a;
		return b;
	}

	__device__ __forceinline double2& operator*=(double2& b, double a) {
		b.x *= a;
		b.y *= a;
		return b;
	}


	__device__ __forceinline double2 operator-(const double2& a, const double2& b) {
		return make_double2(a.x - b.x, a.y - b.y);
	}

	__device__ __forceinline float2 operator-(float2 a, float2 b) {
		return make_float2(a.x - b.x, a.y - b.y);
	}

	__device__ __forceinline double2 operator+(const double2& a, const double2& b) {
		return make_double2(a.x + b.x, a.y + b.y);
	}

	__device__ __forceinline float2 operator+(float2 a, float2 b) {
		return make_float2(a.x + b.x, a.y + b.y);
	}

	__device__ __forceinline double operator&(const double2& a, const double2& b) {
		return a.x * b.x + a.y * b.y;
	}

	__device__ __forceinline float operator&(float2 a, float2 b) {
		return a.x * b.x + a.y * b.y;
	}

	__device__ __forceinline double length(const double2& a) {
		return sqrt(a & a);
	}

	__device__ __forceinline float length(float2 a)	{
		return sqrtf(a & a);
	}
	
	__device__ __forceinline double2 normalize(const double2& a) {
		double inrm = rsqrt(a & a);
		return make_double2(a.x * inrm, a.y * inrm);		
	}

	__device__ __forceinline float2 normalize(float2 a)	{
		float inrm = rsqrtf(a & a);
		return make_float2(a.x * inrm, a.y * inrm);
	}

	__device__ __forceinline unsigned int MExpandBits(unsigned int v)
	{
		// вставит 1 нуль
		v = (v | (v << 8)) & 0x00FF00FF;      //  00000000`00000000`abcdefgh`ijklmnop 
		//                                      | 00000000`abcdefgh`ijklmnop`00000000
		//                                      = 00000000`abcdefgh`XXXXXXXX`ijklmnop
		//                                      & 00000000`11111111`00000000`11111111
		//                                      = 00000000`abcdefgh`00000000`ijklmnop

		v = (v | (v << 4)) & 0x0F0F0F0F;      //  00000000`abcdefgh`00000000`ijklmnop 
		//                                      | 0000abcd`efgh0000`0000ijkl`mnop0000
		//                                      = 0000abcd`XXXXefgh`0000ijkl`XXXXmnop
		//                                      & 00001111`00001111`00001111`00001111
		//                                      = 0000abcd`0000efgh`0000ijkl`0000mnop

		v = (v | (v << 2)) & 0x33333333;      //  0000abcd`0000efgh`0000ijkl`0000mnop 
		//                                      | 00abcd00`00efgh00`00ijkl00`00mnop00
		//                                      = 00abXXcd`00efXXgh`00ijXXkl`00mnXXop
		//                                      & 00110011`00110011`00110011`00110011
		//                                      = 00ab00cd`00ef00gh`00ij00kl`00mn00op

		v = (v | (v << 1)) & 0x55555555;      //  00ab00cd`00ef00gh`00ij00kl`00mn00op 
		//                                      | 0ab00cd0`0ef00gh0`0ij00kl0`0mn00op0
		//                                      = 0aXb0cXd`0eXf0gXh`0iXj0kXl`0mXn0oXp
		//                                      & 01010101`01010101`01010101`01010101
		//                                      = 0a0b0c0d`0e0f0g0h`0i0j0k0l`0m0n0o0p
		return v;
	}


	__device__ __forceinline unsigned int MShrinkBits(unsigned int x)
	{
		x = x & 0x55555555;
		x = (x | (x >> 1)) & 0x33333333;
		x = (x | (x >> 2)) & 0x0F0F0F0F;
		x = (x | (x >> 4)) & 0x00FF00FF;
		x = (x | (x >> 8)) & 0x0000FFFF;
		return x;
	}



	//Знак числа, написана оптимизированная версия, 
	//которая работает самым быстрым возможным образом, т.к. не содержит ни одного оператора условного перехода
	__device__ __forceinline
	int sign(int val) {
		return (0 < val) - (val < 0);
	}

	__device__ __forceinline
	int Delta(const int i, const int j, const int nbodies, const int* __restrict MmortonCodesKeyd) {
		if ((j < 0) || (j > nbodies - 1))
			return -1;

		int count = __clz(MmortonCodesKeyd[i] ^ MmortonCodesKeyd[j]);

		if ((count == 32) && (i != j))
		{
			int add = __clz((i + 1) ^ (j + 1));
			return (2 * codeLength) + (add + (2 * codeLength) - 32);
		}
		return count + (2 * codeLength) - 32;
	}//Delta(...)

	__device__ __forceinline
		int Delta(const int codei, const int i, const int j, const int nbodies, const int* __restrict MmortonCodesKeyd)
	{
		if ((j < 0) || (j > nbodies - 1))
			return -1;

		int count = __clz(codei ^ MmortonCodesKeyd[j]);

		if ((count == 32) && (i != j))
		{
			int add = __clz((i + 1) ^ (j + 1));
			return (2 * codeLength) + (add + (2 * codeLength) - 32);
		}
		return count + (2 * codeLength) - 32;
	}//Delta(...)

	__device__ __forceinline
		int Delta(const int codei, const int codej, const int i, const int j, const int nbodies)
	{
		if ((j < 0) || (j > nbodies - 1))
			return -1;

		int count = __clz(codei ^ codej);

		if ((count == 32) && (i != j))
		{
			int add = __clz((i + 1) ^ (j + 1));
			return (2 * codeLength) + (add + (2 * codeLength) - 32);
		}
		return count + (2 * codeLength) - 32;
	}//Delta(...)


	//Округление "в потолок" результата деления на степень двойки
	__device__ __forceinline
	int ceilpow2(int x, int p) //  =ceil(x / 2^p)
	{
		return (x >> p) + !!(x & ((1 << p) - 1));
	}

	//Округление "в потолок" результата деления пополам
	__device__ __forceinline
	int ceilhalf(int x) //  =ceil(x / 2), т.е. предыдущая функция при p=1
	{
		return (x >> 1) + (x & 1);
	}
	__device__ __forceinline double myMax(double x, double y)
	{
		return (x > y) ? x : y;
	}

	__device__ __forceinline double myMin(double x, double y)
	{
		return (x > y) ? y : x;
	}

	__device__ __forceinline float myMax(float x, float y)
	{
		return (x > y) ? x : y;
	}

	__device__ __forceinline float myMin(float x, float y)
	{
		return (x > y) ? y : x;
	}

	__device__ __forceinline int sqr(int x)
	{
		return x * x;
	}

	__device__ __forceinline double sqr(double x)
	{
		return x * x;
	}

	__device__ __forceinline float sqrf(float x)
	{
		return x * x;
	}

	//Функция вычисления длины общей части префиксов всех частиц в конкретной ячейке
	//__device__ __forceinline
	//int PrefixLength(const int2& range)
	//{		
	//	return min(Delta(range.x, range.y), 2 * 14);
	//}//PrefixLength(...)


	//Функция вычисления общего префикса двух частиц
	//std::pair<unsigned int, int> MortonTree::Prefix(int cell) const
	//{
	//	int length = PrefixLength(cell);
	//	unsigned int el = mortonCodes[mortonTree[cell].range[0]].key;
	//	return { el >> (2 * codeLength - length), length };
	//}//Prefix(...)



}//namespace BHcu

#endif