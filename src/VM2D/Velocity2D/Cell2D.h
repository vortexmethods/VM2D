/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Cell2D.h                                                         |
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
\brief Заголовочный файл с описанием класса Cell
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7
\date 22 ноября 2019 г.
*/

#ifndef CELL_H
#define CELL_H

#include <memory>

#include "Point2D.h"
#include "Vortex2D.h"

namespace VM2D
{

	class Tree;
	struct PointsCopy;
	/*!
	\brief Класс, описывающий ячейку дерева для вычисления скоростей методом Барнса --- Хата

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.7
	\date 22 ноября 2019 г.
	*/
	class Cell
	{
	private:
		/// Габариты прямоугольника
		Point2D leftDown, rightUp;

		/// Признак нижнего уровня (листа дерева)
		bool lowLevel = false;

		/// Массив из двух умных указателей на дочерние ячейки
		std::shared_ptr<Cell> mChildren[2];

		
		/// Суммарные циркуляции вихрей, находящихся внутри данного прямоугольника
		//  sumGam[0] --- сумма положительных циркуляций
		//  sumGam[1] --- сумма отрицательных циркуляций
		double sumGam[2];

		/// Центры завихренностей
		//  centerGam[0] --- центр положительной завихренности 
		//  centerGam[1] --- центр отрицательной завихренности 
		Point2D centerGam[2];

		/// Коэффициенты линейного разложения для вычисления скоростей
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		double d = 0.0;

	public:

		/// Номер уровня данной ячейки
		int level;

		/// Вектор умных указателей на ячейки в ближней зоне (там, где надо считать влияние "напрямую") 
		//имеет смысл только для ячеек нижнего уровня
		std::vector<std::shared_ptr<Cell>> closeCells;

		/// Итераторы начала и конца вектора PointsCopy
		std::vector<PointsCopy>::iterator itBegin;
		std::vector<PointsCopy>::iterator itEnd;


		/// Ссылка на дерево целиком
		Tree& tree;

		/// Конструктор инициализации
		Cell(Tree& Tree_);
		Cell(const Point2D leftDown_, const Point2D rightUp_, Tree& Tree_);
		
		/// Деструктор
		virtual ~Cell();

		///  \brief Построение корня дерева на основе заданных вихрей
		///
		/// \param[in] pointsCopy вектор копий вихрей, на основе которых строится прямоугольник
		void MakeRootCell(std::vector<PointsCopy>& pointsCopy);

		///  \brief Построение иерархической структуры ячеек (дерева)
		///
		/// Рекурсивная процедура
		/// \param[in] itBegin ссылка на итератор начала списка вихрей в данной ячейке
		/// \param[in] itEnd ссылка на итератор конца списка вихрей в данной ячейке
		void CreateAllCells(std::vector<PointsCopy>::iterator &itBegin, std::vector<PointsCopy>::iterator &itEnd);

		/// \brief Вычисление параметров всех ячеек дерева (циркуляций и центров завихренности)			
		///
		/// Рекурсивная процедура
		// для нижних уровней считается по набору вихрей, для остальных --- по потомкам
		void CalculateCellsParams();

		/// \brief Расчет коэффициентов a,b,c,d для ячейки нижнего уровня
		///
		/// Рекурсивная процедура
		/// Вызывается внутри цикла по ячейкам нижнего уровня
		/// \param[in] infCell указатель на влияющую ячейку (при первом витке рекурсии передается rootCell, затем - ее потомки)
		void CalcCoeffToLowLevel(std::shared_ptr<Cell> infCell);

		/// Расчет конвективных скоростей вихрей внутри одной ячейки от вихрей внутри всех ближних ячеек
		void CalcConvVeloByBiotSavart();

		/// \brief Расчет приближенного влияния от вихрей на элемент it внутри ячейки нижнего уровня
		/// Вызывается внутри цикла по вихрям ячеек нижнего уровня
		///
		/// \param[in] it итератор на конкретный элемент, на который рассчитывается влияние
		void CalcInfluenceFromVortexFarCells(std::vector<PointsCopy>::iterator it);

		/// \brief Расчет приближенного влияния от источников на элемент it внутри ячейки нижнего уровня
		/// Вызывается внутри цикла по вихрям ячеек нижнего уровня
		///
		/// \param[in] it итератор на конкретный элемент, на который рассчитывается влияние
		void CalcInfluenceFromSourceFarCells(std::vector<PointsCopy>::iterator it);

		void PrintTree();
		void PrintMyVortexes();
	};
}//namespace VM2D

#endif