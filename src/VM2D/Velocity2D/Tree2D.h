/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.8    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/03/09     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
\version 1.8
\date 09 марта 2020 г.
*/

#ifndef TREE_H
#define TREE_H

#include "Cell2D.h"

namespace VM2D
{

	class World2D;


	/*!
	\brief Класс, описывающий дерево для вычисления скоростей методом Барнса --- Хата

	// Корень дерева

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.8
	\date 09 марта 2020 г.
	*/
	class Tree
	{
	public:

		/// Константная ссылка на решаемую задачу
		const World2D& W;

		/// Параметр точности 
		const double theta = 0.1;

		/// Максимальное количество уровней дерева
		const double numOfLevels = 10;

		/// Минимальный размер ячейки (сумма длины и ширины)
		double minDist;

		/// Минимальное количество вихрей в ячейке
		const double minNumOfVort = 1;

		/// Вектор указателей на ячейки нижнего уровня 
		std::vector<Cell*> lowCells; 
		
		/// Умный указатель на корень
		Cell rootCell;		

		/// Ссылка на объект PointsCopy, содержащий все точки, их скорости, eps*, type
		PointsCopy& allPnt;

		/// \brief Конструктор
		///
		/// Включает в себя построение корня дерева на основе заданных вихрей
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in, out] pointsCopy_ ссылка на вектор из элементов PointsCopy
		Tree(const World2D& W_, PointsCopy& pointsCopy_);

		/// Деструктор
		virtual ~Tree();

	};
}//namespace VM2D

#endif