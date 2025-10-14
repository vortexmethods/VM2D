/*---------------------------------*- BH -*------------------*---------------*\
|        #####   ##  ##         |                            | Version 1.5    |
|        ##  ##  ##  ##         |  BH: Barnes-Hut method     | 2024/06/19     |
|        #####   ######         |  for 2D vortex particles   *----------------*
|        ##  ##  ##  ##         |  Open Source Code                           |
|        #####   ##  ##         |  https://www.github.com/vortexmethods/fastm |
|                                                                             |
| Copyright (C) 2020-2024 I. Marchevsky, E. Ryatina, A. Kolganova             |
*-----------------------------------------------------------------------------*
| File name: BarnesHut.cpp                                                    |
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
\brief Основные операции метода Барнса-Хата
\author Марчевский Илья Константинович
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\version 1.5
\date 19 июня 2024 г.
*/

#include <fstream>

#include "BarnesHut.h"

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
	omp_set_max_active_levels(prm.maxLevelOmp + 1);
#endif


	auto& treeContr = treeVrt;
	auto& pointsCopy = pointsCopyVrt;


		treeVrt->CalculateMortonTreeParams(0, 0);

		double tTreeParamsFinish = omp_get_wtime();
		timeParams += tTreeParamsFinish - tTreeParamsStart;

		double tInflStart = omp_get_wtime();
		
		
		prm.outVortexCoded = -1;
		for (size_t s = 0; s < treeContr->mortonCodes.size(); ++s)
			if (treeContr->mortonCodes[s].originNumber == prm.outVortex)
				prm.outVortexCoded = (int)s;

//		double t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0;


#pragma omp parallel for schedule(dynamic, 10) //reduction(+:t1, t2, t3)
		for (int i = 0; i < (int)treeContr->mortonLowCells.size(); ++i)
		{
//			double tau1, tau2, tau3, tau4;
			auto& lci = treeContr->mortonLowCells[i];
			auto& lowCell = treeContr->mortonTree[lci];

			for (auto& e : lowCell.E)
				e.toZero();					

			if ((lowCell.range[0] <= prm.outVortexCoded) && (lowCell.range[1] >= prm.outVortexCoded))
				prm.outCell = lci;
//			tau1 = omp_get_wtime();
			treeContr->CalcLocalCoeffToLowLevel(lci, treeVrt, 0, true);
//			tau2 = omp_get_wtime();
			treeContr->CalcVeloBiotSavart(lci, treeVrt, calcRadius);
//			tau3 = omp_get_wtime();
			treeContr->CalcVeloTaylorExpansion(lci);
//			tau4 = omp_get_wtime();

//			t1 += (tau2 - tau1);
//			t2 += (tau3 - tau2);
//			t3 += (tau4 - tau3);
		}

//		t4 = -omp_get_wtime();

		/*
		std::ofstream WolfFile("visio.txt");
		const auto& lowCell = treeContr->mortonTree[prm.outCell];
		//std::cout << "lci = " << lci << std::endl;
		WolfFile << "ll:\n" << (lowCell.center - 0.5 * lowCell.size)[0] << " " << (lowCell.center - 0.5 * lowCell.size)[1] << std::endl;
		WolfFile << "ru:\n" << (lowCell.center + 0.5 * lowCell.size)[0] << " " << (lowCell.center + 0.5 * lowCell.size)[1] << std::endl;
		WolfFile << lowCell.range[1] - lowCell.range[0] + 1 << std::endl;
		for (size_t p = lowCell.range[0]; p <= lowCell.range[1]; ++p)
			WolfFile << treeVrt->mortonCodes[p].r[0] << " " << treeVrt->mortonCodes[p].r[1] << std::endl;
		
		
		WolfFile << prm.farCells.size() << std::endl;
		for (const auto c : prm.farCells)
		{
			const auto& cell = treeVrt->mortonTree[c];
			WolfFile << (cell.center - 0.5 * cell.size)[0] << " " << (cell.center - 0.5 * cell.size)[1] << " " <<
				(cell.center + 0.5 * cell.size)[0] << " " << (cell.center + 0.5 * cell.size)[1] << " " << cell.level << std::endl;
		}
		
		


		WolfFile.close();
		//*/

		double tInflStop = omp_get_wtime();
		timeInfl += tInflStop - tInflStart;

		int n = (int)pointsCopy.size();

#pragma omp parallel for 			
		for (int i = 0; i < n; ++i)
		{
			result[i] = IDPI * pointsCopy[i].veloCopy + prm.velInf;
//			ADDOP(2);
		}//for i

		if (calcRadius)
		{
#pragma omp parallel for 			
			for (int i = 0; i < n; ++i)
				epsast[i] = pointsCopy[i].epsast;
		}

//		t4 += omp_get_wtime();

//		std::cout << "t: " << t1 << " " << t2 << " " << t3 << " " << t4 << std::endl;

		}//InfluenceComputation(...)
#endif

}//namespace BH