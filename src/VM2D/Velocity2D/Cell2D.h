/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.9    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/07/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
\brief Заголовочный файл с описанием класса Cell
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.9
\date 22 июля 2020 г.
*/

#ifndef CELL_H
#define CELL_H

#include <memory>

#include "Point2D.h"
#include "Vortex2D.h"
#include "PointsCopy2D.h"

namespace VM2D
{

	class Tree;
	/*!
	\brief Класс, описывающий ячейку дерева для вычисления скоростей методом Барнса --- Хата

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.9
	\date 22 июля 2020 г.
	*/
	class Cell
	{
	private:
		/// Габариты прямоугольника
		Point2D leftDown, rightUp;

		/// Признак нижнего уровня (листа дерева)
		bool lowLevel = false;

		/// Массив из двух умных указателей на дочерние ячейки
		std::unique_ptr<Cell> mChildren[2];
		
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

		/// Расстояния до трех ближайших ячеек (необходимо для корректного вычисления eps*)
		numvector<double, 3> closestCellDist = { 10000.0, 10000.0, 10000.0 };

	public:

		/// Номер уровня данной ячейки
		int level;

		/// Вектор указателей на ячейки в ближней зоне (там, где надо считать влияние "напрямую") 
		//имеет смысл только для ячеек нижнего уровня
		std::vector<Cell*> closeCells;

		/// Итераторы начала и конца объекта PointsCopy
		PointsCopy::iterator itBegin;
		PointsCopy::iterator itEnd;

		/// Ссылка на дерево целиком
		Tree& tree;

		/// \brief Конструктор инициализации
		///
		/// Включает в себя построение корня дерева на основе заданных вихрей
		/// \param[in] Tree_ ссылка на дерево целиком
		Cell(Tree& Tree_);
		Cell(Tree& Tree_, const Point2D leftDown_, const Point2D rightUp_);

		/// Деструктор
		virtual ~Cell();

		void CreateRootCell(PointsCopy& pointsCopy);

		///  \brief Построение иерархической структуры ячеек (дерева)
		///
		/// Рекурсивная процедура
		/// \param[in] itBegin ссылка на итератор начала списка вихрей в данной ячейке
		/// \param[in] itEnd ссылка на итератор конца списка вихрей в данной ячейке
		void CreateAllCells(PointsCopy::iterator& itBegin, PointsCopy::iterator& itEnd);

		/// \brief Вычисление параметров всех ячеек дерева (циркуляций и центров завихренности)			
		///
		/// Рекурсивная процедура
		// для нижних уровней считается по набору вихрей, для остальных --- по потомкам
		void CalculateCellsParams();

		/// \brief Обнуление коэффициентов a,b,c,d для ячейки нижнего уровня 
		///
		/// Рекурсивная процедура
		/// Вызывается перед накоплением коэффициентов процедурой CalcABCDandCloseCellsToLowLevel
		void ClearABCD();


		/// \brief Расчет коэффициентов a,b,c,d для ячейки нижнего уровня и определение ячеек ближней зоны 
		///
		/// Рекурсивная процедура
		/// Вызывается внутри цикла по ячейкам нижнего уровня
		/// \param[in] infCell указатель на влияющую ячейку (при первом витке рекурсии передается rootCell, затем - ее потомки)
		/// \param[in] calcCoef = true, если нужно считать коэффициенты a,b,c,d (не нужно только для виртуальных вихрей)
		void CalcABCDandCloseCellsToLowLevel(Cell& infCell, bool calcCoef);
		
		void CalcCloseZone(Cell& infCell, double distance);

		/// Расчет конвективных скоростей вихрей внутри одной ячейки от вихрей внутри всех ближних ячеек
		///
		/// Вызывается дважды: для wake и VP
		/// \param[in] calcRadius = true, если нужно вычислять eps*
		/// \param[out] savedEe2 - вектор, сохраняющий расстояние до трех ближайших вихрей (имеет смысл при calcRadius = true)
		void CalcConvVeloByBiotSavartFromVortices(bool calcRadius, std::vector<numvector<double, 3>>& savedEe2);
		
		/// Расчет конвективных скоростей вихрей внутри одной ячейки от вихревых слоев (присоединенный + свободный) внутри всех ближних ячеек
		///
		/// Вызывается дважды: для wake и VP
		/// \param[in] calcRadius = true, если нужно вычислять eps* от виртуальных вихрей
		/// \param[in, out] savedEe2 - вектор, сохраняющий расстояние до трех ближайших вихрей (имеет смысл при calcRadius = true)
		void CalcConvVeloByBiotSavartFromSheets(bool calcRadius, std::vector<numvector<double, 3>>& savedEe2);

		/// Расчет eps* по сохраненным расстояниям до трех ближайших элементов
		///
		/// \param[in] savedEe2 - вектор, содержащий расстояние до трех ближайших вихрей
		void SetDomainRadius(const std::vector<numvector<double, 3>>& savedEe2);

		/// Расчет eps* для виртуальных вихрей от всех вихрей (в пелене и виртуальных) внутри всех ближних ячеек
		/// 
		/// Проход по ближней зоне от treeWake и сохранение расстояний до трех ближайших вихрей
		/// Проход по дереву treeSheet, определение ближних ячеек, расчет расстояний до трех ближайших вихрей
		/// Определение eps* 
		void CalcDomainRadiusForVirtualWake();

		/// Расчет конвективных скоростей вихрей внутри одной ячейки от источников внутри всех ближних ячеек
		void CalcConvVeloByBiotSavartFromSources();

		/// \brief Расчет приближенного влияния от вихрей на элемент it внутри ячейки нижнего уровня
		/// Вызывается внутри цикла по вихрям ячеек нижнего уровня
		///
		/// \param[in, out] it итератор на конкретный элемент, на который рассчитывается влияние
		void CalcInfluenceFromVortexFarCells(PointsCopy::iterator it);

		/// \brief Расчет приближенного влияния от источников на элемент it внутри ячейки нижнего уровня
		/// Вызывается внутри цикла по вихрям ячеек нижнего уровня
		///
		/// \param[in] it итератор на конкретный элемент, на который рассчитывается влияние
		void CalcInfluenceFromSourceFarCells(PointsCopy::iterator it);

		Point2D GetR() const
		{
			return 0.5 * (rightUp + leftDown);
		};

		Point2D GetDiag() const
		{
			return rightUp - leftDown;
		};

		void PrintTree();
	};
}//namespace VM2D

#endif