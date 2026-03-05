/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: TreeBH.h                                                         |
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
\brief Заголовок класса Tree для метода Барнса - Хата для CPU
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.14
\date 6 марта 2026 г.
*/

#ifndef TREEBH_H
#define TREEBH_H


#include <memory>

#include "defsBH.h"


namespace BH
{
    /// Глобальная переменная - счетчик количества операций
	extern long long op;

	/*!
	\brief Структура, соответствующая частице и ее мортоновскому коду

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна

	\Version 1.14
	\date 6 марта 2026 г.
	*/
	struct TParticleCode
	{
		/// Мортоновский код частицы
		unsigned int key;

		/// Индекс частицы в глобальном массиве частиц
		int originNumber;

		/// Положение частицы (заносится в эту структуру для облегчения доступа)
		Point2D r;

		/// Циркуляция частицы (заносится в эту структуру для облегчения доступа)
		double g;
	};


	/*!
	\brief Структура, соответствующая ячейке мортоновского дерева

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна

	\Version 1.14
	\date 6 марта 2026 г.
	*/
	struct treeCellT
	{
		/// Индекс родительской ячейки (0 --- для корня дерева)
		int parent;

		/// Индексы потомков (если это не particle): 0 - левый, 1 - правый
		numvector<int, 2> child;

		/// Диапазон индексов частиц в ячейке: 0 - первая; 1 - последняя
		numvector<int, 2> range; 

		/// Признак того, что данная вершина - частица, а не кластер 
		/// (нужно не для обхода дерева, а для "адресации" частиц и кластеров)
		bool particle; 

		/// Признак того, что при обходе дерева на более глубокие уровни спускаться не нужно
		bool leaf;     

		/// Размер ячейки по x и по y (для листа (0;0)
		Point2D size;   

		/// Положение центра ячейки
		Point2D center;

		/// Уровень ячейки в дереве (не заполняется для частиц)
		int level;   

		/// Мультипольные моменты ячейки 
		std::vector<Point2D> mom;

		/// Коэффициенты локального разложения в ячейке для вычисления скоростей
		std::vector<Point2D> E;

		/// Вектор индексов ячеек в ближней зоне 
		/// (там, где надо считать влияние "напрямую"), 
		/// имеет смысл только для нижних уровней
		std::vector<int> closeCells;

		//габариты ячейки
		Point2D minx_miny, maxx_maxy;

		//посчитан или нет габаритный прямоугольник
		bool hasGabs;
	};


	/*!
	\brief Структура, соответствующая мортоновскому дереву
	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна
	\Version 1.14
	\date 6 марта 2026 г.
	*/
	class MortonTree
	{
	private:
		/// Ссылка на параметры, загружаемые из файла
		const params& prm;

		/// Ссылка на список обобщенных вихрей, по которым строится дерево
		std::vector<PointsCopy>& pointsCopy;

	public:
		/// Список мортоновских кодов, вычисляемых по положениям частиц
		std::vector<TParticleCode> mortonCodes;

	private:
		/// Масштабный фактор дерева (на сколько нужно умножить сторону квадрата 
		/// содержащего все частицы, чтобы новый квадрат имел размер
		/// [0, U]x[0, U], где U выбирается так, чтобы все мортоновские коды 
		/// заданной длины были в диапазоне [0, 1), т.е. положение частицы должно
		/// быть в диапазоне [0, 1 / (1.0 + 1/(2^codeLength - 1) ))
		double iQuadSideVar;

		/// Положение нижнего левого угла дерева
		Point2D lowLeft;

		/// Функция вычисления длины общей части префиксов двух 
		/// (а значит - и диапазона) частиц
		/// \param[in] i индекс одной из частиц в отсортированном по мортоновским кодам списке 
		/// \param[in] j индекс второй из частиц в отсортированном по мортоновским кодам списке 
		/// \return длину общей части мортоновского кода
		int Delta(int i, int j) const;

		/// Функция вычисления длины общей части префиксов всех частиц в конкретной ячейке
		/// \param[in] cell индекс ячейки дерева в линейном массиве
		/// \return длину общей части мортоновских кодов всех частиц в данной ячейке
		int PrefixLength(int cell) const;

		/// Функция вычисления общего префикса и его длины для всех частиц в конкретной ячейке
		/// \param[in] cell индекс ячейки дерева в линейном массиве
		/// \return пару из общей части мортоновских кодов всех частиц в данной ячейке и ее длины
		std::pair<unsigned int, int> Prefix(int cell) const;
		
		/// Функция вычисления геометрических параметров внутренней ячейки дерева
		void SetCellGeometry(int cell);

		/// Функция построения i-й внутренней ячейки дерева
		void BuildInternalTreeCell(int i);
		
		/// "Разрежение" двоичного представления беззнакового целого, вставляя по одному нулику между всеми битами
		static unsigned int ExpandBits(unsigned int v);

		/// Мортоновский код для пары из чисел типа double
		/// Исходные числа (r[0], r[1]) - строго в диапазоне [0, 1 / (1.0 + 1/(2^codeLength - 1) ))
		static unsigned int Morton2D(const Point2D& r);

		/// Поиск габаритного прямоугольника системы точек (для OpenMP 2.0)
		static std::pair<Point2D, Point2D> FindEnclosingRectangleOld(const std::vector<PointsCopy>& points);
		
		//Поиск габаритного прямоугольника системы точек
		//работает под OpenMP 4.0
		//static std::pair<Point2D, Point2D> FindEnclosingRectangle(const std::vector<PointsCopy>& points);

		/// Поразрядная сортировка мортоновских кодов
		void RSort_Node3(TParticleCode* m, TParticleCode* m_temp, unsigned int n);

		/// Вспомогательная функция к RSort_Node3
		inline void RSort_step3(TParticleCode* source, TParticleCode* dest, unsigned int n, unsigned int* offset, unsigned char sortable_bit);
		
		/// Параллельная версия поразрядной сортировки мортоновских кодов
		inline void RSort_Parallel(TParticleCode* m, TParticleCode* m_temp, unsigned int n, unsigned int* s);

	public:

		/// Спикок обратных факториалов целых чисел
		std::vector<double> iFact;

		/// список биномиальных коэффициентов 
		/// (имеет смысл только нижняя треугольная часть матрицы, остальные = 0)
		std::vector<double> binomCft;

		/// Максимальная глубина дерева при его обходе
		const int maxTreeLevel;

		/// Сбор статистики по дереву
		/// \param[out] numParticles число частиц в дереве
		/// \param[out] treshold уровень, до которого производится обход
		/// \param[out] partLevel наименьший и наибольший уровни, на которых находятся частицы
		/// \param[out] leafsLevel наименьший и наибольший уровни, на которых ячейки, до которых производится обход
		/// \param[out] lowLevelCells число ячеек "нижнего" уровня (до которых производится обход)
		void getStatistics(int& numParticles, int& treshold, std::pair<int, int>& partLevel, std::pair<int, int>& leafsLevel, int& lowLevelCells) const;

		/// Извлечение биномиального коэффициента из массива binomCft
		/// \return ссылку на C_n^k
		double& getBinom(int n, int k)
		{
			return binomCft[n * (prm.order + 1) + k];
		}
				
		//Конструктор
		MortonTree(const params& prm_, int maxTreeLevel_, std::vector<PointsCopy>& points);
		
		/// Вектор индесков ячеек нижнего уровня в дереве 
		std::vector<int> mortonLowCells;
		//std::vector<std::vector<int>> mortonLowCells;

		/// \brief Мортоновское дерево
		/// Размер дерева равен сумме числа внутренних вершин 
		/// (число частиц без единицы) и числа частиц
		/// резервируется на 1 штуку больше для удобства работы
		/// ячейки [0 ... pos.size()-2] отвечают внутренним узлам,
		/// ячейка с идексом pos.size()-1 не используется 
		/// ячейки [pos.size() ... 2*pos.size()-1] отвечают частицам
		std::vector<treeCellT> mortonTree;

		/// Построение корня дерева и задание его общих параметров
		void MakeRootMortonTree();

		/// Построение внутренних ячеек дерева
		void BuildMortonInternalTree();

		/// Построение верхушек дерева --- отдельных частиц
		void BuildMortonParticlesTree();

		/// Функция заполнения списка нижних ячеек дерева (до которых производится обход)
		// Два варианта заполнения списка нижних вершин: 
		// первый рекурсивный, быстро работает в последовательном варианте
		// второй нерекурсивный, хорошо распараллеливается, но неэффективен в последовательном варианте
		void FillMortonLowCells(int cell = 0);
		void FillMortonLowCellsA();

		/// Расчет мультипольных моментов всех ячеек дерева
		/// \param[in] cell ячейка для вычисления мультипольных моментов
		/// \param[in] omplevel текущий уровень вложенности OpenMP  
		void CalculateMortonTreeParams(int cell, int omplevel);
		
		/// Расчет коэффицентов локального разложения для нижних ячеек 
		/// \param[in] lowCell индекс ячейки нижнего уровня, в которой будет рассчитываться влияние
		/// \param[in] treeInf ссылка на влияющее дерево
		/// \param[in] fromWho индекс влияющей ячейки во влияющем дереве
		/// \param[in] calcCloseTrees признак заполнения списка ближних ячеек (нижнего уровня) к ячейке lowCell
		void CalcLocalCoeffToLowLevel(int lowCell, std::unique_ptr<MortonTree>& treeInf, int fromWho, bool calcCloseTrees);
		
		/// Расчет влияния от ближних ячеек по формуле Био - Савара
		/// \param[in] lowCell индекс ячейки нижнего уровня, в которой рассчитывается влияние
		/// \param[in] treeInf влияющее дерево (чтобы брать ближние ячейки по индексу из нужного массива)
		void CalcVeloBiotSavart(int lowCell, std::unique_ptr<MortonTree>& treeInf);//, bool calcRadius);
		
				
		/// Расчет влияния от дальней зоны при помощи Тейлоровского разложения
		/// \param[in] lowCell --- индекс ячейки нижнего уровня, в которой рассчитывается влияние
		void CalcVeloTaylorExpansion(int lowCell);

		/// Добавление влияния дальней зоны для правой части
		/// \param[in] lowCell --- индекс ячейки нижнего уровня, в которой рассчитывается влияние
		void AddVelo(int lowCell);


		//спускается до листьев рекурсивно и запускает расчет габаритного прямоугольника
		void MinMaxCellsSearch(int numCell);

		//расчет габаритного прямоугольника
		void MinMaxCellsCalc(BH::treeCellT& cell);

		//расчет габаритных прямоугольников при подъеме от листов к корню
		void GabsPyr(int numCell);

		//
		void ShareGabs(int numCell/*, Point2D& minx_miny, Point2D& maxx_maxy*/);

	};
	   
}//namespace BH


#endif
