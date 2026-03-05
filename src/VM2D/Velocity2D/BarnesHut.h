/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: BarnesHut.h                                                      |
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
\brief Заголовок основного класса BarnesHut
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/


#ifndef BARNESHUT_H
#define BARNESHUT_H

#include "TreeBH.h"

namespace BH
{
	/*!
	\brief Класс, определяющий основной алгоритм модификации метода Барнса - Хата

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна

	\Version 1.14
	\date 6 марта 2026 г.
	*/
	class BarnesHut
	{
	public:
		///Ссылка на параметры, считываемые из файла
		const params& prm;
		
		///Список оберток положений вихрей
		std::vector<PointsCopy> pointsCopyVrt;
				
		///Умный yказатель на дерево вихрей
		mutable std::unique_ptr<MortonTree> treeVrt;
		
		/// \brief Конструктор для решения задачи NBODY о вычислении скоростей вихревых частиц
		///
		/// \param[in] prm константная ссылка на параметры, считываемые из файла
		/// \param[in] pointsVrt константная ссылка на список частиц
		BarnesHut(const params& prm_, const std::vector<Vortex2D>& pointsVrt);
		
		/// Деструктор
		~BarnesHut() {};

		/// \brief Построение одного дерева tree на основе заданных точек pointsCopy  
		///
		/// \param[in,out] tree умный указатель на дерево
		/// \param[in] maxTreeLevel максимальная глубина при обходе дерева
		/// \param[in] pointsCopy неконстантная ссылка на список данных (оберток) точек, по которым строится дерево
		/// \param[in,out] time время, затрачиваемое на построение дерева (накопительный итог)
		void BuildOneTree(std::unique_ptr<MortonTree>& tree, int maxTreeLevel, std::vector<PointsCopy>& pointsCopy, double& time);

		/// \brief Построение всех нужных деревьев на основе заданных точек pointsCopy  
		///
		/// \param[in,out] time время, затрачиваемое на построение дерева (накопительный итог)
		void BuildNecessaryTrees(double& time);

		//габаритные прямоугольники для листовых ячеек
		void BuildEnclosingRectangle(double& time);

		/// \brief Расчет влияния в точках дерева, характерных для решаемой задачи (определяется внутри функции) 
		/// 
		/// \param[out] result ссылка на вектор, в который сохраняются вычисленные скорости
		/// \param[in,out] timeParams время расчета параметров деревьев (накопительный итог)
		///	\param[in,out] timeInfl время расчета влияния (накопительный итог)
		void InfluenceComputation(std::vector<Point2D>& result, std::vector<double>& epsast, double& timeParams, double& timeInfl, bool calcRadius);
	};

}//namespace BH

#endif
