/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Tree2D.h                                                         |
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
\brief Заголовочный файл с описанием класса Tree
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7
\date 22 ноября 2019 г.
*/

#ifndef TREE_H
#define TREE_H

#include "Cell2D.h"

namespace VM2D
{

	class World2D;
	struct PointsCopy;


	/*!
	\brief Класс, описывающий дерево для вычисления скоростей методом Барнса --- Хата

	// Корень дерева

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.7
	\date 22 ноября 2019 г.
	*/
	class Tree
	{
	public:

		/// Константная ссылка на решаемую задачу
		const World2D& W;

		/// Параметр точности 
		const double theta = 0.3;

		/// Максимальное количество уровней дерева
		const double numOfLevels = 10;

		/// Минимальный размер ячейки (сумма длины и ширины)
		double minDist;

		/// Минимальное количество вихрей в ячейке
		const double minNumOfVort = 1;

		/// Вектор умных указателей на ячейки нижнего уровня 
		std::vector<std::shared_ptr<Cell>> lowCells; 

		/// Умный указатель на корень
		std::shared_ptr<Cell> rootCell;	/// ????	

		/// \brief Конструктор
		///
		/// \param[in] W_ константная ссылка на решаемую задачу
		Tree(const World2D& W_);

		/// Деструктор
		virtual ~Tree();

		///  \brief Построение корня дерева на основе заданных вихрей
		///
		/// \param[in] pointsCopy вектор из копий вихрей, на основе которых строится дерево
		void BuildTree(std::vector<PointsCopy>& pointsCopy);
	};
}//namespace VM2D

#endif