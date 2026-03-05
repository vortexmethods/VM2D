/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: BarnesHut.cpp                                                    |
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
\brief Основные операции метода Барнса-Хата
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#include <fstream>

#include "BarnesHut.h"
#include "knnCPU.h"

namespace BH
{	
	//Конструктор для вычисления скоростей частиц pointsVrt, влияющих самих на себя
	BarnesHut::BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt)
		: prm(prm_) 
	{
		pointsCopyVrt.insert(pointsCopyVrt.end(), pointsVrt.begin(), pointsVrt.end());
	}//BarnesHut(...)

	

	//Построение одного дерева
	void BarnesHut::BuildOneTree(std::unique_ptr<MortonTree>& tree, int maxTreeLevel, std::vector<PointsCopy>& pointsCopy, double& time)
	{
		tree = std::make_unique<MortonTree>(prm, maxTreeLevel, pointsCopy);

		double t1 = omp_get_wtime();
		tree->MakeRootMortonTree();
		tree->BuildMortonInternalTree();
		tree->BuildMortonParticlesTree();

		// Один из двух вариантов заполнения нижних вершин
		//tree->FillMortonLowCellsA(); // Рекурсивный
		tree->FillMortonLowCells();  // Нерекурсивный

		double t2 = omp_get_wtime();

		time += t2-t1;
	}


	//Построение всех нужных деревьев на основе заданных точек pointsCopy  
	void BarnesHut::BuildNecessaryTrees(double& time)
	{
#ifdef needTreeVrt
		BuildOneTree(treeVrt, prm.NumOfLevelsVortex, pointsCopyVrt, time); // дерево вихрей
#endif		
	}

	//габаритные прямоугольники для листовых ячеек и объединение их при подъеме
	void BarnesHut::BuildEnclosingRectangle(double& time)
	{
		double t1 = omp_get_wtime();

#ifdef needTreeVrt
		//treeVrt->MinMaxCellsSearch(0); //габаритные прямоугольники для листовых ячеек

		for (auto& cell : treeVrt->mortonLowCells)
			treeVrt->MinMaxCellsSearch(cell);

		treeVrt->GabsPyr(0); //расчет габаритных прямоугольников при подъеме от листов к корню (но на самом деле наоборот)
#endif

		double t2 = omp_get_wtime();
		time += t2 - t1;

		
		
	}

	// Расчет влияния 
#ifndef CALCSHEET
	void BarnesHut::InfluenceComputation(std::vector<Point2D>& result, std::vector<double>& epsast, double& timeParams, double& timeInfl, bool calcRadius)
	{
		double tTreeParamsStart = omp_get_wtime();

#ifdef OLD_OMP
		omp_set_nested(1);
#else
		omp_set_max_active_levels(3);
		std::cout << "OMP_LEVEL = 2" << std::endl;
		//omp_set_max_active_levels(prm.maxLevelOmp + 1);
#endif


		auto& treeContr = treeVrt;
		auto& pointsCopy = pointsCopyVrt;


		treeVrt->CalculateMortonTreeParams(0, 0);

		double tTreeParamsFinish = omp_get_wtime();
		timeParams += tTreeParamsFinish - tTreeParamsStart;

		double tInflStart = omp_get_wtime();

#pragma omp parallel for schedule(dynamic, 10) //reduction(+:t1, t2, t3)
		for (int i = 0; i < (int)treeContr->mortonLowCells.size(); ++i)
		{
			auto& lci = treeContr->mortonLowCells[i];
			auto& lowCell = treeContr->mortonTree[lci];

			for (auto& e : lowCell.E)
				e.toZero();

			treeContr->CalcLocalCoeffToLowLevel(lci, treeVrt, 0, true);
			treeContr->CalcVeloBiotSavart(lci, treeVrt);
			treeContr->CalcVeloTaylorExpansion(lci);
		}


		//CPU - neib
		const size_t knb = 3;
		std::vector<std::vector<std::pair<double, size_t>>> initdist(pointsCopy.size());
		for (auto& d : initdist)
			d.resize(2 * knb, { -1.0, -1 });
		double timeKnn = -omp_get_wtime();
		WakekNNnewForEpsast(pointsCopy, knb, initdist);//CPU

#pragma omp parallel for
		for (int j = 0; j < pointsCopy.size(); ++j)
		{
			double sd2 = (initdist[j][0].first + initdist[j][1].first + initdist[j][2].first) / 3.0;
			pointsCopy[j].epsast = (sd2 > 0) ? sqrt(sd2) : 1000.0;
		}
		timeKnn += omp_get_wtime();

		double tInflStop = omp_get_wtime();
		timeInfl += tInflStop - tInflStart;

		int n = (int)pointsCopy.size();

#pragma omp parallel for 			
		for (int i = 0; i < n; ++i)
			result[i] = IDPI * pointsCopy[i].veloCopy;

		if (calcRadius)
		{
#pragma omp parallel for 			
			for (int i = 0; i < n; ++i)
				epsast[i] = pointsCopy[i].epsast;
		}
	}//InfluenceComputation(...)
#endif

}//namespace BH