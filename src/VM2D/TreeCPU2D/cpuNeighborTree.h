/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: cpuNeighborTree.h                                                |
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
\brief Простое дерево для точного поиска соседей в заданном радиусе на CPU

Точки упорядочиваются по мортоновским кодам и разбиваются на листья по leafSize точек;
над листьями строится полное бинарное дерево габаритных прямоугольников.
Используется для расчета сумм с компактным носителем (диффузионные скорости I1, I2),
результат совпадает с прямым перебором с точностью до порядка суммирования.

\version 1.14
\date 6 марта 2026 г.
*/

#ifndef CPUNEIGHBORTREE_H
#define CPUNEIGHBORTREE_H

#include <algorithm>
#include <cstdint>
#include <limits>
#include <vector>

#include "Point2D.h"

namespace VM2D
{
	/// Габаритный прямоугольник
	struct NbBox
	{
		double minx, miny, maxx, maxy;

		void SetEmpty()
		{
			minx = miny = std::numeric_limits<double>::infinity();
			maxx = maxy = -std::numeric_limits<double>::infinity();
		}

		void Add(const Point2D& p)
		{
			minx = std::min(minx, p[0]); maxx = std::max(maxx, p[0]);
			miny = std::min(miny, p[1]); maxy = std::max(maxy, p[1]);
		}

		void Add(const NbBox& b)
		{
			minx = std::min(minx, b.minx); maxx = std::max(maxx, b.maxx);
			miny = std::min(miny, b.miny); maxy = std::max(maxy, b.maxy);
		}

		bool IsEmpty() const { return minx > maxx; }

		/// Квадрат расстояния до другого прямоугольника (для пустого --- бесконечность)
		double Dist2(const NbBox& b) const
		{
			if (IsEmpty() || b.IsEmpty())
				return std::numeric_limits<double>::infinity();
			const double dx = std::max(0.0, std::max(minx - b.maxx, b.minx - maxx));
			const double dy = std::max(0.0, std::max(miny - b.maxy, b.miny - maxy));
			return dx * dx + dy * dy;
		}

		/// Квадрат расстояния до точки
		double Dist2(const Point2D& p) const
		{
			if (IsEmpty())
				return std::numeric_limits<double>::infinity();
			const double dx = std::max(0.0, std::max(minx - p[0], p[0] - maxx));
			const double dy = std::max(0.0, std::max(miny - p[1], p[1] - maxy));
			return dx * dx + dy * dy;
		}
	};


	class NeighborTree
	{
	public:
		static const int leafSize = 32;

		int n = 0;
		int nLeaves = 0;
		int firstLeaf = 1;            ///< номер первого листа в неявном бинарном дереве (степень двойки)

		std::vector<int> idx;         ///< исходные номера точек в порядке сортировки
		std::vector<Point2D> pos;     ///< координаты точек в порядке сортировки
		std::vector<NbBox> box;       ///< прямоугольники узлов: корень 1, потомки узла k --- 2k и 2k+1

		int LeafBegin(int leaf) const { return leaf * leafSize; }
		int LeafEnd(int leaf) const { return std::min(n, (leaf + 1) * leafSize); }
		const NbBox& LeafBox(int leaf) const { return box[firstLeaf + leaf]; }

		/// Построение дерева; getPos(i) возвращает координаты i-й точки
		template<typename F>
		void Build(int n_, F getPos)
		{
			n = n_;
			idx.resize(n);
			pos.resize(n);

			NbBox all;
			all.SetEmpty();
			for (int i = 0; i < n; ++i)
				all.Add(getPos(i));

			std::vector<uint64_t> keys(n), tmp(n);
			const double scale = std::max(all.maxx - all.minx, all.maxy - all.miny);
			const double iscale = (scale > 0.0) ? (65535.0 / scale) : 0.0;

#pragma omp parallel for
			for (int i = 0; i < n; ++i)
			{
				const Point2D& p = getPos(i);
				const uint32_t ix = (uint32_t)((p[0] - all.minx) * iscale);
				const uint32_t iy = (uint32_t)((p[1] - all.miny) * iscale);
				keys[i] = ((uint64_t)((Spread(ix) << 1) | Spread(iy)) << 32) | (uint32_t)i;
			}

			//поразрядная сортировка по старшим 32 битам (мортоновский код)
			for (int shift = 32; shift < 64; shift += 8)
			{
				size_t cnt[257] = { 0 };
				for (int i = 0; i < n; ++i)
					++cnt[((keys[i] >> shift) & 0xFF) + 1];
				for (int b = 0; b < 256; ++b)
					cnt[b + 1] += cnt[b];
				for (int i = 0; i < n; ++i)
					tmp[cnt[(keys[i] >> shift) & 0xFF]++] = keys[i];
				keys.swap(tmp);
			}

			for (int i = 0; i < n; ++i)
			{
				idx[i] = (int)(keys[i] & 0xFFFFFFFFu);
				pos[i] = getPos(idx[i]);
			}

			nLeaves = (n + leafSize - 1) / leafSize;
			firstLeaf = 1;
			while (firstLeaf < std::max(nLeaves, 1))
				firstLeaf <<= 1;

			box.resize(2 * firstLeaf);
#pragma omp parallel for
			for (int l = 0; l < firstLeaf; ++l)
			{
				NbBox& b = box[firstLeaf + l];
				b.SetEmpty();
				if (l < nLeaves)
					for (int k = LeafBegin(l); k < LeafEnd(l); ++k)
						b.Add(pos[k]);
			}
			for (int k = firstLeaf - 1; k >= 1; --k)
			{
				box[k] = box[2 * k];
				box[k].Add(box[2 * k + 1]);
			}
		}

		/// Список листьев, прямоугольники которых находятся на расстоянии меньше R от прямоугольника q
		void CollectLeaves(const NbBox& q, double R, std::vector<int>& leaves) const
		{
			leaves.clear();
			if (nLeaves == 0)
				return;
			const double R2 = R * R;
			int stack[64];
			int top = 0;
			stack[top++] = 1;
			while (top > 0)
			{
				const int k = stack[--top];
				if (box[k].Dist2(q) >= R2)
					continue;
				if (k >= firstLeaf)
					leaves.push_back(k - firstLeaf);
				else
				{
					stack[top++] = 2 * k + 1;
					stack[top++] = 2 * k;
				}
			}
		}

	private:
		static uint32_t Spread(uint32_t v)
		{
			v &= 0x0000FFFF;
			v = (v | (v << 8)) & 0x00FF00FF;
			v = (v | (v << 4)) & 0x0F0F0F0F;
			v = (v | (v << 2)) & 0x33333333;
			v = (v | (v << 1)) & 0x55555555;
			return v;
		}
	};

}//namespace VM2D

#endif
