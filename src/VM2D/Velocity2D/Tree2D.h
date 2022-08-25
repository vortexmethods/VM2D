/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.11   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2022/08/07     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
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
\brief Заголовочный файл с описанием класса Tree
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.11
\date 07 августа 2022 г.
*/

#ifndef TREE_H
#define TREE_H

#include "Cell2D.h"
#include "WakeDataBase2D.h"

namespace VM2D
{

	class World2D;


	/*!
	\brief Класс, описывающий дерево для вычисления скоростей методом Барнса --- Хата

	// Корень дерева

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.11
	\date 07 августа 2022 г.
	*/
	class Tree
	{
	public:

		/// Константная ссылка на решаемую задачу
		const World2D& W;

		/// Параметр точности 
		const double theta = 0.5;

		/// Максимальное количество уровней дерева
		double numOfLevels; 

		/// Минимальный размер ячейки (сумма длины и ширины)
		double minDist;

		/// Минимальное количество вихрей в ячейке
		const double minNumOfVort = 1;

		/// Вектор указателей на ячейки нижнего уровня 
		std::vector<Cell*> lowCells; 
		
		/// Корень дерева
		Cell rootCell;		

		/// Объект PointsCopy, содержащий все точки, их скорости, eps*, type
		PointsCopy allPnt;

		/// \brief Конструктор
		///
		/// Включает в себя построение корня дерева на основе заданных вихрей
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] points константная ссылка на набор точек, по которым строится дерево (treeWake, treeWakeVP, treeSourcesWake)
		/// \param[in] type тип точек, по которым строится дерево
		/// \param[in] numOfLevels_ максимальное количество уровней дерева
		Tree(const World2D& W_, const WakeDataBase& points, const PointType type, int numOfLevels_);
		
		/// \brief Конструктор
		///
		/// Включает в себя построение корня дерева на основе заданных вихрей
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] type тип точек, по которым строится дерево (treeSheetsGam, treeSheetsSource)
		/// \param[in] numOfLevels_ максимальное количество уровней дерева
		Tree(const World2D& W_, const PointType type, int numOfLevels_);

		/// Деструктор
		virtual ~Tree();

		/// \brief Создание копии данных для точек (для сортировки) 
		///
		/// Используется для заполнения allPnt в конструкторе для деревьев treeWake, treeWakeVP, treeSourceWake
		///
		/// \param[in] points массив точек, указатели на которые необходимо запомнить
		/// \param[in] type тип массива точек (wake, sourceWake или wakeVP)
		void CreatePointsCopy(const WakeDataBase& points, const PointType type);

		/// \brief Создание массива указателей на массив точек (для сортировки)
		///
		/// Используется для заполнения allPnt в конструкторе для деревьев treeSheetsGam, treeSource
		///
		/// \param[in] type тип массива точек (sheetGam, source)
		void CreateSheetsCopy(const PointType type);
	};
}//namespace VM2D

#endif