/*--------------------------------*- BHgpu -*----------------*---------------*\
| #####   ##  ##                |                            | Version 1.5    |
| ##  ##  ##  ##   ####  ##  ## |  BHgpu: Barnes-Hut method  | 2023/08/29     |
| #####   ######  ##     ##  ## |  for 2D vortex particles   *----------------*
| ##  ##  ##  ##  ##     ##  ## |  Open Source Code                           |
| #####   ##  ##   ####   ####  |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2023 I. Marchevsky, E. Ryatina, A. Kolganova             |
| Copyright (C) 2013, Texas State University-San Marcos. All rights reserved. |
*-----------------------------------------------------------------------------*
| File name: operations.cuh                                                   |
| Info: Source code of BHgpu                                                  |
|                                                                             |
| This file is part of BHgpu.                                                 |
| BHcu is free software: you can redistribute it and/or modify it             |
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
| along with BHgpu.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Вспомогательные операции
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 29 августа 2023 г.
*/

#ifndef OPERATIONS_H_
#define OPERATIONS_H_

#include <cuda.h>


namespace BHcu
{

#ifndef WIN32
#define __forceinline inline
#endif

	__device__ __forceinline real2 multz(real2 a, real2 b)
	{
		real2 res;
		res.x = a.x * b.x - a.y * b.y;
		res.y = a.x * b.y + a.y * b.x;
		return res;
	}

	__device__ __forceinline real2 multz(real ax, real ay, real2 b)
	{
		real2 res;
		res.x = ax * b.x - ay * b.y;
		res.y = ax * b.y + ay * b.x;
		return res;
	}
	__device__ __forceinline real2 multzA(real2 a, real2 b)
	{
		real2 res;
		res.x =  a.x * b.x + a.y * b.y;
		res.y = -a.x * b.y + a.y * b.x;
		return res;
	}	
	
	__device__ __forceinline real2 multzA(real2 a, real bx, real by)
	{
		real2 res;
		res.x = a.x * bx + a.y * by;
		res.y = a.y * bx - a.x * by;
		return res;
	}

	__device__ __forceinline real2 operator*(real a, real2 b)
	{
		real2 res;
		res.x = a * b.x;
		res.y = a * b.y;
		return res;
	}

	__device__ __forceinline float2 operator*(float a, float2 b)
	{
		float2 res;
		res.x = a * b.x;
		res.y = a * b.y;
		return res;
	}


	__device__ __forceinline real2 operator*(int a, real2 b)
	{
		real2 res;
		res.x = a * b.x;
		res.y = a * b.y;
		return res;
	}

	__device__ __forceinline float2 operator*(int a, float2 b)
	{
		float2 res;
		res.x = a * b.x;
		res.y = a * b.y;
		return res;
	}

	__device__ __forceinline real2 operator*(real2 b, real a)
	{
		real2 res;
		res.x = a * b.x;
		res.y = a * b.y;
		return res;
	}



	__device__ __forceinline float2 operator*(float2 b, float a)
	{
		float2 res;
		res.x = a * b.x;
		res.y = a * b.y;
		return res;
	}

	__device__ __forceinline real2 operator*(real2 b, int a)
	{
		real2 res;
		res.x = a * b.x;
		res.y = a * b.y;
		return res;
	}


	__device__ __forceinline float2 operator*(float2 b, int a)
	{
		float2 res;
		res.x = a * b.x;
		res.y = a * b.y;
		return res;
	}

	__device__ __forceinline real2 operator/(real2 b, real a)
	{
		real2 res;
		res.x = b.x / a;
		res.y = b.y / a;
		return res;
	}

	__device__ __forceinline float2 operator/(float2 b, float a)
	{
		float2 res;
		res.x = b.x / a;
		res.y = b.y / a;
		return res;
	}

	__device__ __forceinline real2& operator+=(real2& a, real2 b)
	{
		a.x += b.x;
		a.y += b.y;
		return a;
	}

	__device__ __forceinline float2& operator+=(float2& a, float2 b)
	{
		a.x += b.x;
		a.y += b.y;
		return a;
	}

	__device__ __forceinline real2& operator-=(real2& a, real2 b)
	{
		a.x += b.x;
		a.y += b.y;
		return a;
	}

	__device__ __forceinline float2& operator-=(float2& a, float2 b)
	{
		a.x += b.x;
		a.y += b.y;
		return a;
	}


	__device__ __forceinline real2 operator-(real2 a, real2 b)
	{
		real2 res;
		res.x = a.x - b.x;
		res.y = a.y - b.y;
		return res;
	}

	__device__ __forceinline float2 operator-(float2 a, float2 b)
	{
		float2 res;
		res.x = a.x - b.x;
		res.y = a.y - b.y;
		return res;
	}

	__device__ __forceinline real2 operator+(real2 a, real2 b)
	{
		real2 res;
		res.x = a.x + b.x;
		res.y = a.y + b.y;
		return res;
	}

	__device__ __forceinline float2 operator+(float2 a, float2 b)
	{
		float2 res;
		res.x = a.x + b.x;
		res.y = a.y + b.y;
		return res;
	}

	__device__ __forceinline real operator&(real2 a, real2 b)
	{
		return a.x * b.x + a.y * b.y;
	}

	__device__ __forceinline float operator&(float2 a, float2 b)
	{
		return a.x * b.x + a.y * b.y;
	}

	__device__ __forceinline real length(real2 a)
	{
		return sqrt(a & a);
	}

	__device__ __forceinline float length(float2 a)
	{
		return sqrtf(a & a);
	}
	
	__device__ __forceinline real2 normalize(real2 a)
	{
		real nrm = sqrt(a & a);
		real2 res{ a.x / nrm, a.y / nrm };
		return res;
	}


	__device__ __forceinline float2 normalize(float2 a)
	{
		float nrm = sqrtf(a & a);
		float2 res{ a.x / nrm, a.y / nrm };
		return res;
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
	int sign(int val)
	{
		return (0 < val) - (val < 0);
	}

	__device__ __forceinline
	int Delta(const int i, const int j, const int nbodies, const int* __restrict MmortonCodesKeyd)
	{
		if ((j < 0) || (j > nbodies - 1))
			return -1;

		int count = __clz(MmortonCodesKeyd[i] ^ MmortonCodesKeyd[j]);

		if ((count == 32) && (i != j))
		{
			int add = __clz((i + 1) ^ (j + 1));
			return 28 + (add + 28 - 32);
		}
		return count + 28 - 32;
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