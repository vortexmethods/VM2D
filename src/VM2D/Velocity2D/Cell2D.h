/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
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
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
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
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.10
	\date 17 мая 2021 г.
	*/
	class Cell
	{
	private:
		/// Габариты прямоугольника
		Point2D leftDown, rightUp;

		/// Центр прямоугольника
		Point2D posCentre;

		/// Номер моего уровня (глубина дерева)
		int level;
		/// Номер потомка внутри уровня
		int p;

		/// Признак нижнего уровня (листа дерева)
		bool lowLevel = false;

		/// Массив из двух умных указателей на дочерние ячейки
		std::unique_ptr<Cell> mChildren[2];
		
		/// Монопольный момент ячейки 
		double mon;
		/// Дипольный, квадрупольный, октупольный и гексадекапольный моменты ячейки 
		Point2D dip, qua, oct, hex;

		/// \brief Коэффициенты для вычисления скоростей
		numvector<Point2D, 4> E;

		/// Расстояния до трех ближайших ячеек (необходимо для корректного вычисления eps*)
		numvector<double, 3> closestCellDist = { 10000.0, 10000.0, 10000.0 };

	public:
		///Ссылка на вектор всех исходных точек 
		std::vector<PointsCopy>& points;

		/// Вектор указателей на ячейки в ближней зоне (там, где надо считать влияние "напрямую") 
		//имеет смысл только для ячеек нижнего уровня
		std::vector<Cell*> closeCells;

		/// Индексы элементов в исходном массиве точек
		std::vector<int> index;

		/// Ссылка на дерево целиком
		Tree& tree;

		/// \brief Конструктор инициализации
		///
		/// Включает в себя построение корня дерева на основе заданных вихрей
		/// \param[in] Tree_ ссылка на дерево целиком
		Cell(Tree& Tree_);
		
		/// Деструктор
		virtual ~Cell();

		void CreateRootCell(PointsCopy& pointsCopy);

		///  \brief Построение иерархической структуры ячеек (дерева)
		///
		/// Рекурсивная процедура
		void CreateAllCells();

		/// \brief Вычисление параметров всех ячеек дерева (циркуляций и центров завихренности)			
		///
		/// Рекурсивная процедура
		void CalculateCellsParams();

		/// Вычисление параметров нижней ячейки по набору вихрей (для расчета скоростей)
		void CalcPointsParams();
		/// Вычисление параметров нижней ячейки по набору панелей (для решения СЛАУ)
		void CalcPanelsParams();

		/// \brief Обнуление коэффициентов локальных разложений		
		///
		/// Рекурсивная процедура
		/// Вызывается перед накоплением коэффициентов процедурой CalcABCDandCloseCellsToLowLevel
		void ClearCoeff();

		/// \brief Расчет коэффициентов локальных разложений для ячейки нижнего уровня и определение ячеек ближней зоны 
		///
		/// Рекурсивная процедура
		/// Вызывается внутри цикла по ячейкам нижнего уровня
		/// \param[in] infCell указатель на влияющую ячейку (при первом витке рекурсии передается rootCell, затем - ее потомки)
		/// \param[in] calcCoef = true, если нужно считать коэффициенты локальных разложений (не нужно только для виртуальных вихрей)
		/// \param[in] calcCloseTrees = true (по умолчанию), если нужно заполнить вектор ближних ячеек
		void CalcLocalCoeffandCloseCellsToLowLevel(Cell& infCell, bool calcCoef= true, bool calcCloseTrees = true);
		
		/// \brief Нахождение ближних ячеек для расчета диффузионных скоростей - другой критерий
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

		/// \brief Расчет приближенного влияния от вихрей на один элемент внутри ячейки нижнего уровня
		/// Вызывается внутри цикла по вихрям ячеек нижнего уровня
		///
		void CalcInfluenceFromVortexFarCells();

		/// \brief Расчет приближенного влияния от источников на один элемент внутри ячейки нижнего уровня
		/// Вызывается внутри цикла по вихрям ячеек нижнего уровня
		void CalcInfluenceFromSourceFarCells();

		/// Расчет влияния от слоев внутри одного уровня от панелей всех ближних уровней (для решения СЛАУ)
		void CalcInfluenceFromPanels();
		void CalcInfluenceFromPanelsToPoints();

		/// Обновление влияния от слоев внутри одного уровня от панелей всех ближних уровней (для решения СЛАУ)
		/// 
		/// Вызывается внутри функции IterativeInfluenceComputation для всех итераций, кроме первой 
		/// Коэффициенты СЛАУ не рассчитываются напрямую, а берутся из сохраненного массива velSave
		void UpdateInfluence();

		void PushbackLowTrees();

		Point2D GetR() const
		{
			return 0.5 * (rightUp + leftDown);
		};

		Point2D GetDiag() const
		{
			return rightUp - leftDown;
		};
	};
}//namespace VM2D

//умножение комплексных чисел
inline Point2D multz(const Point2D& a, const Point2D& b)
{
#ifdef calcOp
	op += 4;
#endif
	return Point2D({ a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0] });
}

// умножение a на комплексно сопряженноe к b
inline Point2D multzA(const Point2D& a, const Point2D& b)
{
#ifdef calcOp
	op += 4;
#endif
	return Point2D({ a[0] * b[0] + a[1] * b[1], a[1] * b[0] - a[0] * b[1] });
}
#endif