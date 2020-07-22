/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.9    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/07/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
\brief Заголовочный файл с описанием класса VelocityBiotSavart
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.9
\date 22 июля 2020 г.
*/

#ifndef VELOCITY2DBARNESHUT_H
#define VELOCITY2DBARNESHUT_H

#include "Velocity2D.h"
#include "Tree2D.h"
#include "Vortex2D.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс, определяющий способ вычисления скоростей

	Способ вычисления скоростей
	- быстрый метод типа Барнса - Хата

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.9
	\date 22 июля 2020 г.
	*/
	class VelocityBarnesHut : public Velocity
	{
	private:

		/// \brief Объекты из точек:
		/// pointsCopyWake - вихрей в пелене и источников,
		/// pointsCopyVP - в которых считается скорость и давление,	
		/// sheetsGamCopy - маркеров для вихревых слоев,
		/// sourcesCopy -  источников в области течения + маркеры для присоединенного слоя источников
		//PointsCopy pointsCopyWake, pointsCopyVP, sheetsGamCopy, sourcesCopy;

	public:		

		/// \brief Умные указатели на деревья: 
		/// treeWake содержит вихри из вихревого следа и источники,  
		/// treeVP - точки, в которых считается скорость и давление, 
		/// treeSheetsGam - маркеры для вихревых слоев,
		/// treeSources - источники в области течения + маркеры для присоединенного слоя источников 
		//std::unique_ptr<Tree> treeWake, treeVP, treeSheetsGam, treeSources;

		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу 
		VelocityBarnesHut(const World2D& W_);

		/// Деструктор
		virtual ~VelocityBarnesHut();

		//реализация виртуальных функций
		virtual void CalcConvVelo() override;
		virtual void FillRhs(Eigen::VectorXd& rhs) const override;

		
#if defined (USE_CUDA)
		void GPUCalcConvVeloToSetOfPoints(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius, bool calcVelo, bool calcRadius) {};
#endif
		/// \brief Генерация вектора влияния вихревого следа на профиль
		///
		/// Генерирует вектор влияния вихревого следа на профиль, используемый затем для расчета вектора правой части.
		/// 
		/// \param[out] wakeRhs ссылка на вектор влияния вихревого следа на ВСЕ профили
		/// \warning Использует OMP, MPI
		/// \ingroup Parallel
		void GetWakeInfluenceToRhs(std::vector<double>& wakeRhs) const;

		void GetWakeInfluenceToRhsBS(const Airfoil& afl, std::vector<double>& wakeRhs) const;
#if defined(USE_CUDA)
		void GPUGetWakeInfluenceToRhs(const Airfoil& afl, std::vector<double>& wakeRhs) const;
#endif
		void CalcDiffVeloI1I2ToWakeFromSheets(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2) override;
		void CalcDiffVeloI1I2ToWakeFromWake(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2) override;
	};



}//namespace VM2D

#endif