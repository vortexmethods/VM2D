/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Velocity2DBarnesHut.h                                            |
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
\brief Заголовочный файл с описанием класса VelocityBiotSavart
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7
\date 22 ноября 2019 г.
*/

#ifndef VELOCITY2DBARNESHUT_H
#define VELOCITY2DBARNESHUT_H

#include "Velocity2D.h"
#include "Tree2D.h"

namespace VM2D
{

	class World2D;

	/// Класс перечисления для определения типа набора точек (пелена/виртуальные вихри/точки для вычисления скорости и давления)
	enum class PointType{ wake,  wakeVP, sheetGam, sourceWake, source };
	   	 

	/// Структура, содержащая копию вихря, его тип, номер (откуда пришел), скорость и domainRadius 
	/// Не содержит константных членов, так как вектор из PointsCopy сортируется
	struct PointsCopy : public Vortex2D
	{		
		int aflN; //если wake или wakeVP, то -1
		int num; //если wake, sourceWake или wakeVP, то num - номер в wake или wakeVP; если sheet, то num - номер в соответствующем профиле aflN
		PointType type;
		Point2D veloCopy;
		double domainRadiusCopy;

		PointsCopy(const Vortex2D& vtx_, int afl_, int num_, PointType type_) :
			Vortex2D(vtx_), aflN(afl_), num(num_), type(type_), veloCopy({ 0.0, 0.0 }), domainRadiusCopy(0.0)
		{};
	};

	/*!
	\brief Класс, определяющий способ вычисления скоростей

	Способ вычисления скоростей
	- быстрый метод типа Барнса - Хата

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.7
	\date 22 ноября 2019 г.
	*/
	class VelocityBarnesHut : public Velocity
	{
	private:

		/// \brief Векторы из копий точек:
		/// pointsCopyWake - вихрей в пелене и источников,
		/// pointsCopyVP - в которых считается скорость и давление,	
		/// sheetsGamCopy - маркеров для вихревых слоев,
		/// sourcesCopy -  источников в области течения + маркеры для присоединенного слоя источников
		std::vector<PointsCopy> pointsCopyWake, pointsCopyVP, sheetsGamCopy, sourcesCopy;


	public:		

		/// \brief Умные указатели на деревья: 
		/// treeWake содержит вихри из вихревого следа и источники,  
		/// treeVP - точки, в которых считается скорость и давление, 
		/// treeSheetsGam - маркеры для вихревых слоев,
		/// treeSources - источники в области течения + маркеры для присоединенного слоя источников 
		std::unique_ptr<Tree> treeWake, treeVP, treeSheetsGam, treeSources;

		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу 
		VelocityBarnesHut(const World2D& W_);

		/// Деструктор
		virtual ~VelocityBarnesHut();

		/// \brief Создание копии данных для точек (для сортировки) 
		///
		/// Используется для массивов pointsCopyVP и pointsCopyWake
		///
		/// \param[out] pointsCopy вектор из указателей на массив точек (заполняется)
		/// \param[in] points массив точек, указатели на которые необходимо запомнить
		/// \param[in] type тип массива точек (wake, sourceWake или wakeVP)
		void CreatePointsCopy(std::vector<PointsCopy>& pointsCopy, const WakeDataBase& points, const PointType type);

		/// \brief Создание массива указателей на массив точек (для сортировки)
		///
		/// Используется для массива sheetsCopy
		///
		void CreateSheetsCopy(std::vector<PointsCopy>& pointsCopy, const PointType type);

		/// Построение дерева treeWake и расчет его параметров
		void BuildTrees(PointType type);

		//реализация виртуальных функций
		//virtual void CalcDiffVeloI1I2ToSetOfPoints(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2) override {};
		//virtual void CalcDiffVeloI1I2ToSetOfPointsFromPanels(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2) override {};
		virtual void CalcConvVelo(timePeriod& convWakeTime, timePeriod& convVPTime) override;

		virtual void FillRhs(Eigen::VectorXd& rhs) const override;

		
#if defined (USE_CUDA)
		void GPUCalcConvVeloToSetOfPoints(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius) {};
#endif
		/// \brief Генерация вектора влияния вихревого следа на профиль
		///
		/// Генерирует вектор влияния вихревого следа на профиль, используемый затем для расчета вектора правой части.
		/// 
		/// \param[in] afl константная ссылка на профиль, правая часть для которого вычисляется
		/// \param[out] wakeRhs ссылка на вектор влияния вихревого следа на ВСЕ профили
		/// \warning Использует OMP, MPI
		/// \ingroup Parallel
		void GetWakeInfluenceToRhs(std::vector<double>& wakeRhs) const;
#if defined(USE_CUDA)
		void GPUGetWakeInfluenceToRhs(std::vector<double>& wakeRhs) const {};
#endif

	};



}//namespace VM2D

#endif