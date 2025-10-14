/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.5    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2024/06/19     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2024 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: Tree.cpp                                                         |
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
\brief Основные операции работы с деревом
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 19 июня 2024 г.
*/

#include <algorithm>
#include <cstring>

#include "Tree.h"

namespace BH
{
	extern long long op;

	//Конструктор
	MortonTree::MortonTree(const params& prm_, int maxTreeLevel_, std::vector<PointsCopy>& points)
		: prm(prm_), pointsCopy(points), maxTreeLevel(maxTreeLevel_)
	{
		mortonTree.resize(2 * points.size());
		for (auto& c : mortonTree)
		{
			c.mom.resize(prm.order + 1);
			c.E.resize(prm.order);
		}

		mortonCodes.resize(points.size());
		mortonLowCells.reserve(points.size());
		

		//Вычисляем факториалы и биномиальные коэффиценты
		iFact.resize(prm.order + 1);
		binomCft.resize((prm.order + 1) * (prm.order + 1));
		
		iFact[0] = 1.0;
		for (int i = 1; i <= prm.order; ++i)
			iFact[i] = iFact[i - 1] / i;
		getBinom(0, 0) = 1.0;
		for (int i = 1; i <= prm.order; ++i) 
		{
			getBinom(i, 0) = 1.0;
			getBinom(i, i) = 1.0;
			for (int j = 1; j < i; ++j)
				getBinom(i, j) = getBinom(i - 1, j - 1) + getBinom(i - 1, j);
		}//for i
	}//MortonTree(...)


	//Поиск габаритного прямоугольника системы точек
	std::pair<Point2D, Point2D> MortonTree::FindEnclosingRectangleOld(const std::vector<PointsCopy>& points)
	{
		double shared_minx = 1e+10, shared_maxx = -1e+10, shared_miny = 1e+10, shared_maxy = -1e+10;
#pragma omp parallel 
		{
			double minx = 1e+10, maxx = -1e+10, miny = 1e+10, maxy = -1e+10;
#pragma omp for nowait
			for (int ii = 0; ii < (int)points.size(); ++ii)
			{
				const Point2D& r = points[ii].r();
				minx = std::min(r[0], minx);
				maxx = std::max(r[0], maxx);
				miny = std::min(r[1], miny);				
				maxy = std::max(r[1], maxy);
			}//for ii
#pragma omp critical 
			{
				shared_minx = std::min(shared_minx, minx);
				shared_maxx = std::max(shared_maxx, maxx);
				shared_miny = std::min(shared_miny, miny);
				shared_maxy = std::max(shared_maxy, maxy);
			}//critical
		}//parallel
		
		return { {shared_minx, shared_miny}, {shared_maxx, shared_maxy} };
	}//FindEnclosingRectangleOld(...)


/*
//Поиск габаритного прямоугольника системы точек
#pragma omp declare reduction(min : struct PointsCopy : \
        omp_out.r()[0] = omp_in.r()[0] > omp_out.r()[0]  ? omp_out.r()[0] : omp_in.r()[0], \
        omp_out.r()[1] = omp_in.r()[1] > omp_out.r()[1]  ? omp_out.r()[1] : omp_in.r()[1] ) \
        initializer( omp_priv = Vortex2D{ { 1e+10, 1e+10 } , 0.0} )

#pragma omp declare reduction(max : struct PointsCopy : \
        omp_out.r()[0] = omp_in.r()[0] < omp_out.r()[0]  ? omp_out.r()[0] : omp_in.r()[0],  \
        omp_out.r()[1] = omp_in.r()[1] < omp_out.r()[1]  ? omp_out.r()[1] : omp_in.r()[1] ) \
        initializer( omp_priv = Vortex2D{ { -1e+10, -1e+10 } , 0.0} )

	
	std::pair<Point2D, Point2D> MortonTree::FindEnclosingRectangle(const std::vector<PointsCopy>& points)
	{
		PointsCopy minp(Vortex2D{ { 1e+10, 1e+10 } , 0.0}), maxp(Vortex2D{ { -1e+10, -1e+10 } , 0.0 });
				
#pragma omp parallel for reduction(min:minp) reduction(max:maxp)
		for (int i = 0; i < (int)points.size(); ++i) {
			if (points[i].r()[0] < minp.r()[0]) minp.r()[0] = points[i].r()[0];
			if (points[i].r()[1] < minp.r()[1]) minp.r()[1] = points[i].r()[1];
			if (points[i].r()[0] > maxp.r()[0]) maxp.r()[0] = points[i].r()[0];
			if (points[i].r()[1] > maxp.r()[1]) maxp.r()[1] = points[i].r()[1];
		}//for i
		return { minp.r(), maxp.r() };		
	}//FindEnclosingRectangle(...)
*/


	//Построение корня дерева и задание его общих параметров
	void MortonTree::MakeRootMortonTree()
	{
		auto gabs = FindEnclosingRectangleOld(pointsCopy);

		double iLeft = gabs.first[0];
		double iBottom = gabs.first[1];
		
		double iRight = gabs.second[0];
		double iTop = gabs.second[1];

		Point2D posCentre = { 0.5 * (iLeft + iRight), 0.5 * (iBottom + iTop) };

		double quadSide = std::max(iRight - iLeft, iTop - iBottom);
		iQuadSideVar = 1.0 / (quadSide * (1.0 + 1.0 / ((1 << codeLength) - 1)));
		lowLeft = posCentre - 0.5*quadSide*Point2D{ 1.0, 1.0 };
				
		//iQuadSideVar = 1.0;
		//lowLeft = Point2D{ 0.0, 0.0 };
				
#pragma omp parallel for schedule(dynamic, 1000)
		for (int i = 0; i < (int)pointsCopy.size(); ++i)
		{
			Point2D locR = (pointsCopy[i].r() - lowLeft) * iQuadSideVar;
			unsigned int mortonCode = Morton2D(locR);
				
			auto& mc = mortonCodes[i];
			mc.key = mortonCode;
			mc.originNumber = i;
			mc.r = pointsCopy[i].r();
			mc.g = pointsCopy[i].g();
		}//for i

		//std::cout << "omp_threads = " << omp_get_max_threads() << std::endl;

		//RSort_Node3(mortonCodes.data(), mortonCodes_temp.data(), (int)pointsCopy.size());				
		
		//Временные массивы для сортировки		
		std::vector<unsigned int> s(256 * omp_get_max_threads());
		std::vector<TParticleCode> mortonCodes_temp(pointsCopy.size());
		RSort_Parallel(mortonCodes.data(), mortonCodes_temp.data(), (int)pointsCopy.size(), s.data());

		mortonTree[0].range = { 0, (int)pointsCopy.size() - 1 };
		mortonTree[0].particle = false;
	}//MakeRootMortonTree()


	//Функция вычисления длины общей части префиксов двух(а значит - и диапазона) частиц
	int MortonTree::Delta(int i, int j) const
	{
		if ((j < 0) || (j > (int)pointsCopy.size() - 1))
			return -1;

		if (i > j)
			std::swap(i, j);

		//if ((i < 0) || (j > n-1))
		//    exit(111);

		const unsigned int& ki = mortonCodes[i].key;
		const unsigned int& kj = mortonCodes[j].key;
		
		//Поиск номера самого старшего ненулевого бита в числе c 
		int count = 0;
		for (unsigned int c = (ki ^ kj); c; c >>= 1, ++count);

		if ((!count) && (i != j))
		{
			int addCount = 0;
			//единички к номерам i и j добавлены для совместимости с Wolfram Mathematica, 
			//для кода они не важны, но дерево без них почти наверняка построится по-другому        
			for (unsigned int add = ((i + 1) ^ (j + 1)); add; add >>= 1, ++addCount);
			return 2 * codeLength + (2 * codeLength - addCount);
		}//if ((!count) && (i != j))

		return (2 * codeLength - count);
	}//Delta(...)


	//Функция вычисления длины общей части префиксов всех частиц в конкретной ячейке
	int MortonTree::PrefixLength(int cell) const
	{
		const auto& range = mortonTree[cell].range;
		return std::min(Delta(range[0], range[1]), 2 * codeLength);
	}//PrefixLength(...)

	
	//Функция вычисления общего префикса двух частиц
	std::pair<unsigned int, int> MortonTree::Prefix(int cell) const
	{
		int length = PrefixLength(cell);
		unsigned int el = mortonCodes[mortonTree[cell].range[0]].key;
		return { el >> (2 * codeLength - length), length };
	}//Prefix(...)


	//Функция вычисления геометрических параметров внутренней ячейки дерева
	void MortonTree::SetCellGeometry(int cell)
	{
		int prLength;
		unsigned int pr;
		std::tie<unsigned int, int>(pr, prLength) = Prefix(cell);
		prLength -= PrefixLength(0);

		Point2D sz = { 1.0 / (double)(1 << ceilhalf(prLength)), 1.0 / (1 << (prLength / 2)) };

		Point2D pos = 0.5 * sz;

		for (int i = 0; i < prLength; i += 2)
			if (pr & (1 << (prLength - i - 1)))
				pos[0] += 1.0 / (1 << (1 + i / 2));

		for (int i = 1; i < prLength; i += 2)
			if (pr & (1 << (prLength - i - 1)))
				pos[1] += 1.0 / (1 << (1 + i / 2));

		mortonTree[cell].center = pos * (1.0 / iQuadSideVar) + lowLeft;
		mortonTree[cell].size = sz * (1.0 / iQuadSideVar);
		mortonTree[cell].level = prLength;

		mortonTree[cell].leaf = (prLength >= maxTreeLevel);
	}//SetCellGeometry(...)


	//Функция построения i-й внутренней ячейки дерева
	void MortonTree::BuildInternalTreeCell(int i)
	{
		const int n = (int)pointsCopy.size();
		int d = sign(Delta(i, i + 1) - Delta(i, i - 1));
		int delta_min = Delta(i, i - d);
		
		int Lmax = 2;
		while (Delta(i, i + Lmax * d) > delta_min)
			Lmax *= 2;

		int L = 0;
		for (int t = (Lmax >> 1); t >= 1; t >>= 1)		
			if (Delta(i, i + (L + t) * d) > delta_min)
				L += t;
		
		int j = i + L * d;
		int delta_node = Delta(i, j);
		int s = 0;
		for (int p = 1, t = ceilhalf(L); L > (1 << (p - 1)); ++p, t = ceilpow2(L, p))
		{
			int dl = Delta(i, i + (s + t) * d);
			if (dl > delta_node)
				s += t;
		}//for p
		int gamma = i + s * d + std::min(d, 0);

		auto min_max = std::minmax(i, j);

		const int& left = gamma;
		const int& right = gamma + 1;

		treeCellT& treeCell = mortonTree[i];

		// Левый потомок - лист или внутренний узел
		bool ifLeftParticle = (min_max.first == gamma);
		treeCell.child[0] = ifLeftParticle ? n + left : left;
		mortonTree[treeCell.child[0]].range = { min_max.first, gamma };
		mortonTree[treeCell.child[0]].particle = ifLeftParticle;
		mortonTree[treeCell.child[0]].parent = i;

		// Правый потомок - лист или внутренний узел
		bool ifRightParticle = (min_max.second == gamma + 1);
		treeCell.child[1] = ifRightParticle ? n + right : right;
		mortonTree[treeCell.child[1]].range = { gamma + 1, min_max.second };
		mortonTree[treeCell.child[1]].particle = ifRightParticle;
		mortonTree[treeCell.child[1]].parent = i;
	}//BuildInternalTreeCell(...)


	//Построение внутренних ячеек дерева
	void MortonTree::BuildMortonInternalTree()
	{
#pragma omp parallel for schedule(dynamic, 1000)
		for (int q = 0; q < (int)pointsCopy.size() - 1; ++q)
			BuildInternalTreeCell(q);
		
		mortonTree[pointsCopy.size() - 1].leaf = false;
		mortonTree[pointsCopy.size() - 1].particle = false;

#pragma omp parallel for schedule(dynamic, 1000)
		for (int q = 0; q < (int)pointsCopy.size() - 1; ++q)
			SetCellGeometry(q);
	}//BuildMortonInternalTree()


	//Построение верхушек дерева --- отдельных частиц
	void MortonTree::BuildMortonParticlesTree()
	{
#pragma omp parallel for		
		for (int q = 0; q < (int)pointsCopy.size(); ++q)
		{
			auto& pt = mortonCodes[q];

			auto& cell = mortonTree[pointsCopy.size() + q];
			cell.center = pt.r;

			cell.leaf = true;
			cell.size = { 0.0, 0.0 };
		}//for q
	}//BuildMortonParticlesTree(...)


	//Заполнение списка нижних вершин: 
	//рекурсивный алгоритм, быстро работает в последовательном варианте		
	void MortonTree::FillMortonLowCells(int cell)
	{
		if (mortonTree[cell].leaf)
			mortonLowCells.push_back(cell);
		else
		{
			FillMortonLowCells(mortonTree[cell].child[0]);
			FillMortonLowCells(mortonTree[cell].child[1]);
		}//else
	}//FillMortonLowCells(...)


	//Заполнение списка нижних вершин: 
	//нерекурсивный алгоритм, хорошо распараллеливается, но неэффективен в последовательном варианте
	void MortonTree::FillMortonLowCellsA()
	{
#pragma omp parallel
		{
			std::vector<int> good_matches_private;
#pragma omp for nowait
			for (int cell = 0; cell < (int)mortonTree.size(); ++cell)			
				if (!mortonTree[mortonTree[cell].parent].leaf && mortonTree[cell].leaf)				
					good_matches_private.push_back(cell);				
#pragma omp critical
			mortonLowCells.insert(mortonLowCells.end(), good_matches_private.begin(), good_matches_private.end());
		}//parallel
	}//FillMortonLowCellsA()


	/// Вычисление параметров дерева (циркуляций и высших моментов)
	void MortonTree::CalculateMortonTreeParams(int cell, int omplevel)
	{
		auto& cl = mortonTree[cell];
		if (!cl.particle)
		{
#pragma omp parallel /*default(none)*/ shared(cl, omplevel) num_threads(2) if (omplevel < prm.maxLevelOmp + 1)
			{				
				std::vector<Point2D> h(prm.order);
#pragma omp for 
				for (int s = 0; s < 2; ++s)
					CalculateMortonTreeParams(cl.child[s], omplevel + 1);

#pragma omp single
				{					
					for (auto& m : cl.mom)
						m.toZero();
										
					for (int s = 0; s < 2; ++s)//цикл по двум потомкам
					{
						auto& chld = mortonTree[cl.child[s]];
						const Point2D& g = chld.mom[0];
						cl.mom[0] += g;

						
						h[0] = (chld.center) - (cl.center);
						for (int i = 1; i < prm.order; ++i)
							h[i] = multz(h[i - 1], h[0]);
						
						for (int p = 1; p <= prm.order; ++p) 
						{
							Point2D shiftMom = chld.mom[p];
							for (int k = 1; k <= p; ++k)
							{
								shiftMom += getBinom(p, k) *  multz(chld.mom[p - k], h[k - 1]);
								//ADDOP(2);
							}//for k
							cl.mom[p] += shiftMom;
						}//for m
						
					}//for s
				}//single				
			}//parallel
		}//if(!particle)
		else
		{
			//расчет моментов ячейки нижнего уровня для вихрей или для панелей
				//C учетом того, что дерево строится до частиц --- только монопольный момент отличен от нуля
				cl.mom[0][0] = mortonCodes[cl.range[0]].g;
				cl.mom[0][1] = 0.0;

				for (int p = 1; p <= prm.order; ++p)
					cl.mom[p].toZero();
		}//else
	}//CalculateMortonTreeParams(...)		

	//Расчет коэффициентов разложения в ряд Тейлора внутри ячейки нижнего уровня
	void MortonTree::CalcLocalCoeffToLowLevel(int lowCell, std::unique_ptr<MortonTree>& treeInf, int fromWho, bool calcCloseTrees)
	{		
		auto& lt = mortonTree[lowCell];
		//if (calcCloseTrees)
		//	lt.closeCells.resize(0);

		if (treeInf.get() != this || lowCell != fromWho)
		{
			auto& wh = treeInf->mortonTree[fromWho];

			double h, h0;
			h = wh.size[0] + wh.size[1];
			h0 = lt.size[0] + lt.size[1];

			Point2D& r0 = lt.center;
			Point2D& r1 = wh.center;

			Point2D rr = r0 - r1;

			double crit2 = rr.length2();
			//ADDOP(3);
			//ADDOP(3);
			// если выполнен критерий дальности => считаем коэффициенты
			if ((crit2 >= sqr((h0 + h + 2.0 * prm.eps) / prm.theta)) && (wh.level <= prm.NumOfLevelsVortex))
			{
				if (wh.range[1] == wh.range[0])
					lt.closeCells.push_back(fromWho);
				else
				{
					if ((lt.range[0] <= prm.outVortexCoded) && (prm.outVortexCoded <= lt.range[1]))
					{
						//std::cout << fromWho << " " << wh.level << std::endl;
						prm.farCells.push_back(fromWho);
					}

					Point2D  rr = r0 - r1;
					rr *= 1.0 / rr.length2();
					//ADDOP(3);

					Point2D varTheta = rr;
					for (int diag = 0; diag < prm.order; ++diag)
					{
						for (int q = 0; q <= diag; ++q)
						{
							lt.E[q] += ((q % 2 ? -1.0 : 1.0) * iFact[diag - q]) * multzA(varTheta, wh.mom[diag - q]);
							//ADDOP(3);
						}
						varTheta = (diag + 1) * multz(varTheta, rr);
						//ADDOP(2);
					}
					lt.E[0] += iFact[prm.order] * multzA(varTheta, wh.mom[prm.order]);
					//ADDOP(2);
				}
			}//if crit2
			else // если не выполнен критерий, то рекурсия 
			{
				if ((!wh.leaf) && (wh.level < prm.NumOfLevelsVortex))
				{
					CalcLocalCoeffToLowLevel(lowCell, treeInf, wh.child[0], calcCloseTrees);
					CalcLocalCoeffToLowLevel(lowCell, treeInf, wh.child[1], calcCloseTrees);
				}
				else if (calcCloseTrees)
				{	
					lt.closeCells.push_back(fromWho);					
				}
			}//else if
		}//if (lowTree != this)
		else if (calcCloseTrees)
		{
			lt.closeCells.push_back(lowCell); //себя тоже добавляем в ближнюю зону 			
		}
	}//CalcLocalCoeffToLowLevel(...)



	//Расчет влияния от ближних ячеек по формуле Био - Савара
	void MortonTree::CalcVeloBiotSavart(int lowCell, std::unique_ptr<MortonTree>& treeInf, bool calcRadius)
	{
		double d_1, d_2, d_3, dst23, dst12;

		auto& lc = mortonTree[lowCell];
		for (int i = lc.range[0]; i <= lc.range[1]; ++i)		
		{
			//Локальные переменные для цикла
			Point2D velI, velLinI;
			double dst2eps;
			double r2;
						
			d_1 = 1e+5;
			d_2 = 1e+5;
			d_3 = 1e+5;			

			for (size_t k = 0; k < lc.closeCells.size(); ++k)
			{
				velI.toZero();
				velLinI.toZero();			

				const Point2D& posI = mortonCodes[i].r; 

				auto& rg = treeInf->mortonTree[lc.closeCells[k]].range;
				for (int j = rg[0]; j <= rg[1]; ++j)
				{
					auto& pt = treeInf->mortonCodes[j];

					const Point2D& posJ = pt.r;
					const double& gamJ = pt.g;
					r2 = (posI - posJ).length2();
					dst2eps = std::max(r2, prm.eps2);
					velI += (gamJ / dst2eps) * (posI - posJ).kcross();

					if (calcRadius)
					{
						if ((r2 < d_3) && (r2 > 0))
						{
							//dst23 = realmin(r2, d_2);
							dst23 = (r2 < d_2) ? r2 : d_2;

							//d_3 = realmax(r2, d_2);
							d_3 = (r2 > d_2) ? r2 : d_2;

							//dst12 = realmin(dst23, d_1);
							dst12 = (dst23 < d_1) ? dst23 : d_1;

							//d_2 = realmax(dst23, d_1);
							d_2 = (dst23 > d_1) ? dst23 : d_1;

							d_1 = dst12;

							//if (mortonCodes[i].originNumber == 0)
							//	std::cout << "d: " << d_1 << " " << d_2 << " " << d_3 << std::endl;
						}
					}
					//ADDOP(5);					
				}//for j

				pointsCopy[mortonCodes[i].originNumber].veloCopy += velI;		
			}//for k 

			if (calcRadius)
			{
				//if (mortonCodes[i].originNumber == 0)
				//	std::cout << "fin: " << "lowCell = " << lowCell << " " << d_1 << " " << d_2 << " " << d_3 << " " << sqrt((d_1 + d_2 + d_3) / 3) << std::endl;
				pointsCopy[mortonCodes[i].originNumber].epsast = sqrt((d_1 + d_2 + d_3) / 3);
			}
		}//for i
	}//CalcVeloBiotSavart(...)


	// Расчет влияния от дальней зоны при помощи Тейлоровского разложения
	void MortonTree::CalcVeloTaylorExpansion(int lowCell)
	{
		auto& lc = mortonTree[lowCell];
		for (int i = lc.range[0]; i <= lc.range[1]; ++i)
		{
			Point2D deltaPos = mortonCodes[i].r - lc.center;
			Point2D v = lc.E[0];


			Point2D hh = deltaPos;
			for (int k = 1; k < prm.order - 1; ++k) 
			{
				v += iFact[k] * multzA(lc.E[k], hh);
				hh = multz(hh, deltaPos);
				//ADDOP(2);
			}//for k
			v += iFact[prm.order-1] * multzA(lc.E[prm.order-1], hh);
			//ADDOP(2);

			pointsCopy[mortonCodes[i].originNumber].veloCopy += Point2D({ -v[1], v[0] });	

		}//for i
	}//CalcVeloTaylorExpansion(...)
	

	//"Разрежение" двоичного представления беззнакового целого, вставляя по одному нулику между всеми битами
	unsigned int MortonTree::ExpandBits(unsigned int v)
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


	//Мортоновский код для пары из чисел типа double
	//Исходное число - строго в диапазоне [0, 1 / (1.0 + 1/(2^codeLength - 1) ))
	unsigned int MortonTree::Morton2D(const Point2D& r)
	{
		/*
		if (x < 0) std::cout << "x*... < 0" << std::endl;
		if (y < 0) std::cout << "y*... < 0" << std::endl;

		if (x * twoPowCodeLength > twoPowCodeLengthMinus)
			std::cout << "x*... > ..." << std::endl;

		if (y * twoPowCodeLength > twoPowCodeLengthMinus)
			std::cout << "y*... > ..." << std::endl;

		x = std::min(std::max(x * twoPowCodeLength, 0.0), twoPowCodeLengthMinus);
		y = std::min(std::max(y * twoPowCodeLength, 0.0), twoPowCodeLengthMinus);

		unsigned int xx = expandBits((unsigned int)x);
		unsigned int yy = expandBits((unsigned int)y);
		*/
		const Point2D& rscale = twoPowCodeLength * r;
		const unsigned int& xx = ExpandBits((unsigned int)(rscale[0]));
		const unsigned int& yy = ExpandBits((unsigned int)(rscale[1]));
		return yy | (xx << 1);
	}



	//========================================================
	void MortonTree::RSort_step3(TParticleCode* source, TParticleCode* dest, unsigned int n, unsigned int* offset, unsigned char sortable_bit)
	{
		//unsigned char* b = (unsigned char*)&source[n].key + sortable_bit;

		for (unsigned int i = 0; i < n; ++i)
		{
			TParticleCode* src = &source[i];
			unsigned int off = (src->key >> (sortable_bit * 8)) & 0xFF;
			dest[offset[off]++] = *src;
		}
	}
	//========================================================
	void MortonTree::RSort_Node3(TParticleCode* m, TParticleCode* m_temp, unsigned int n)
	{
		// Заводим массив корзин
		unsigned int s[sizeof(m->key) * 256] = { 0 };
		// Заполняем массив корзин для всех разрядов
		for (unsigned int i = 0; i < n; ++i)
		{
			unsigned int key = m[i].key;
			for (unsigned int digit = 0; digit < sizeof(m->key); digit++)
				++s[((key >> (digit << 3)) & 0xFF) + (digit << 8)];
		}

		// Пересчитываем смещения для корзин
		for (unsigned int digit = 0; digit < sizeof(m->key); digit++)
		{
			unsigned int off = 0;
			for (unsigned int i = 0; i < 256; i++)
			{
				auto& rs = s[i + (digit << 8)];
				unsigned int value = rs;
				rs = off;
				off += value;
			}
		}

		// Вызов сортировки по битам от младших к старшим (LSD)
		for (unsigned int digit = 0; digit < sizeof(m->key); digit++)
		{
			RSort_step3(m, m_temp, n, &s[digit << 8], digit);
			std::swap(m, m_temp);
			//TParticleCode* temp = m;
			//m = m_temp;
			//m_temp = temp;
		}

		// Если ключ структуры однобайтовый, копируем отсортированное в исходный массив
		if (sizeof(m->key) == 1)
		{
			std::swap(m, m_temp);
			//TParticleCode* temp = m;
			//m = m_temp;
			//m_temp = temp;
			memcpy(m, m_temp, n * sizeof(TParticleCode));
		}
	}


	void MortonTree::RSort_Parallel(TParticleCode* m, TParticleCode* m_temp, unsigned int n, unsigned int* s)
	{
		if (n == 0)
			return;
		// Количество задействованных потоков
		unsigned char threads = omp_get_max_threads();
		//std::cout << "threads = " << (int)threads << std::endl;

#pragma omp parallel num_threads(threads)
		{
			TParticleCode* source = m;
			TParticleCode* dest = m_temp;
			unsigned int l = omp_get_thread_num();
			unsigned int div = n / omp_get_num_threads();
			unsigned int mod = n % omp_get_num_threads();
			unsigned int left_index = l < mod ? (div + (mod == 0 ? 0 : 1)) * l : n - (omp_get_num_threads() - l) * div;
			unsigned int right_index = left_index + div - (mod > l ? 0 : 1);

			for (unsigned int digit = 0; digit < sizeof(m->key); ++digit)
			{
				unsigned int s_sum[256] = { 0 };
				unsigned int s0[256] = { 0 };
				unsigned char* b1 = (unsigned char*)&source[right_index].key;
				unsigned char* b2 = (unsigned char*)&source[left_index].key;
				while (b1 >= b2)
				{
					++s0[*(b1 + digit)];
					b1 -= sizeof(TParticleCode);
				}
				for (unsigned int i = 0; i < 256; i++)
				{
					s[i + 256 * l] = s0[i];
				}

#pragma omp barrier
				for (unsigned int j = 0; j < threads; j++)
				{
					for (unsigned int i = 0; i < 256; i++)
					{
						s_sum[i] += s[i + 256 * j];
						if (j < l)
						{
							s0[i] += s[i + 256 * j];
						}
					}
				}

				for (unsigned int i = 1; i < 256; ++i)
				{
					s_sum[i] += s_sum[i - 1];
					s0[i] += s_sum[i - 1];
				}
				unsigned char* b = (unsigned char*)&source[right_index].key + digit;
				TParticleCode* v1 = &source[right_index];
				TParticleCode* v2 = &source[left_index];
				while (v1 >= v2)
				{
					dest[--s0[*b]] = *v1--;
					b -= sizeof(TParticleCode);
				}
#pragma omp barrier
				std::swap(source, dest);
			}
		}

		// Если ключ структуры однобайтовый, просто копируем в исходный массив
		if (sizeof(m->key) == 1)
		{
			memcpy(m, m_temp, n * sizeof(TParticleCode));
		}
	}

	void MortonTree::getStatistics(int& numParticles, int& treshold, std::pair<int, int>& partLevel, std::pair<int, int>& leafsLevel, int& lowLevelCells) const
	{
		numParticles = (int)pointsCopy.size();
		treshold = maxTreeLevel;

		int leafsLevelMin = (int)1e+9;
		int leafsLevelMax = 0;
		for (int q = 0; q < numParticles - 1; ++q)
		{
			if (!mortonTree[mortonTree[q].parent].leaf)
			{
				if (mortonTree[q].leaf)
				{
					int level = mortonTree[q].level;

					if (level > leafsLevelMax)
						leafsLevelMax = level;

					if (level < leafsLevelMin)
						leafsLevelMin = level;

					continue;
				}

				if ((mortonTree[mortonTree[q].child[0]].particle) || (mortonTree[mortonTree[q].child[1]].particle))
				{
					int level = mortonTree[q].level + 1;
					if (level > leafsLevelMax)
						leafsLevelMax = level;

					if (level < leafsLevelMin)
						leafsLevelMin = level;

					continue;
				}
			}
		}
		leafsLevel = { leafsLevelMin, leafsLevelMax };


		int partLevelMin = (int)1e+9;
		int partLevelMax = 0;
		for (int q = 0; q < numParticles - 1; ++q)
		{
			if ((mortonTree[mortonTree[q].child[0]].particle) || (mortonTree[mortonTree[q].child[1]].particle))
			{
				int level = mortonTree[q].level;
				if (level > partLevelMax)
					partLevelMax = level;

				if (level < partLevelMin)
					partLevelMin = level;
			}
		}
		partLevel = { partLevelMin + 1, partLevelMax + 1 };
		lowLevelCells = (int)mortonLowCells.size();
	}



	void MortonTree::MinMaxCellsSearch(int numCell)
	{
		auto& cell = mortonTree[numCell];
		if (cell.leaf)
		{
			MinMaxCellsCalc(cell);
			cell.hasGabs = true;
		}
		else
		{
			std::cout << "!!!!!!!!!!!!" << std::endl;

			int ch0 = cell.child[0];
			int ch1 = cell.child[1];
			MinMaxCellsSearch(ch0);
			MinMaxCellsSearch(ch1);
			cell.hasGabs = false;
		}
	}

	void MortonTree::MinMaxCellsCalc(BH::treeCellT& cell)
	{
		int spoint = cell.range[0];
		int epoint = cell.range[1];
		
		double minx = 1e+10, maxx = -1e+10, miny = 1e+10, maxy = -1e+10;

		for (int ii = spoint; ii <= epoint; ++ii)
		{
			const Point2D& r = mortonCodes[ii].r;
			minx = std::min(r[0], minx);
			maxx = std::max(r[0], maxx);
			miny = std::min(r[1], miny);
			maxy = std::max(r[1], maxy);
		}
		cell.minx_miny[0] = minx;
		cell.minx_miny[1] = miny;
		cell.maxx_maxy[0] = maxx;
		cell.maxx_maxy[1] = maxy;

		cell.center = (cell.minx_miny + cell.maxx_maxy) * 0.5;
		cell.size = cell.maxx_maxy - cell.minx_miny;

	}

	void MortonTree::ShareGabs(int numCell/*, Point2D& minx_miny, Point2D& maxx_maxy*/)
	{
		auto& lcell = mortonTree[mortonTree[numCell].child[0]]; //левый ребенок 	
		auto& rcell = mortonTree[mortonTree[numCell].child[1]]; //правый ребенок 
		auto& pcell = mortonTree[numCell];
		{
			pcell.minx_miny[0] = std::min(lcell.minx_miny[0], rcell.minx_miny[0]);
			pcell.minx_miny[1] = std::min(lcell.minx_miny[1], rcell.minx_miny[1]);
			pcell.maxx_maxy[0] = std::max(lcell.maxx_maxy[0], rcell.maxx_maxy[0]);
			pcell.maxx_maxy[1] = std::max(lcell.maxx_maxy[1], rcell.maxx_maxy[1]);
		}

		pcell.center = (pcell.minx_miny + pcell.maxx_maxy) * 0.5;
		pcell.size = pcell.maxx_maxy - pcell.minx_miny;
	}

	void MortonTree::GabsPyr(int numCell)
	{
		auto& cell = mortonTree[numCell];
		if (cell.hasGabs)
			return;

		int ch0 = cell.child[0];
		int ch1 = cell.child[1];
		if (!mortonTree[ch0].hasGabs)
			GabsPyr(ch0);
		if (!mortonTree[ch1].hasGabs)
			GabsPyr(ch1);
		
		ShareGabs(numCell/*, cell.minx_miny, cell.maxx_maxy*/);
	}


}//namespace BH

