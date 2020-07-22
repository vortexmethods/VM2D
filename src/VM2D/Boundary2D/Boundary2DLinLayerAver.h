/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.9    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2020/07/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2020 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Boundary2DLinLayerAver.h                                         |
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
\brief Заголовочный файл с описанием класса BoundaryLinLayerAver
\warning Пока реализованы не все методы
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.9
\date 22 июля 2020 г.
*/
#ifndef BOUNDARYLINLAYERAVER_H
#define BOUNDARYLINLAYERAVER_H

#include "Boundary2D.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс, определяющий способ удовлетворения граничного условия на обтекаемом профиле

	Способ удовлетворения граничного условия:
	- вихревой слой представляется кусочно-линейном виде с разрывами в каждом узле профиля;
	- условие ортогональности невязки граничного условия константным и линейным проекционным функциям.

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.9
	\date 22 июля 2020 г.
	*/
	class BoundaryLinLayerAver : public Boundary
	{

	public:

		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] numberInPassport_ номер профиля в паспорте задачи
		BoundaryLinLayerAver(const World2D& W_, size_t numberInPassport_) :
			Boundary(W_, numberInPassport_, 2)
		{ };

		/// Деструктор
		virtual ~BoundaryLinLayerAver() {};

		//далее -- реализации виртуальных функций
		virtual void FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol) override;
		virtual void FillIQSelf(std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ) override;
		virtual void FillMatrixFromOther(const Boundary& otherBoundary, Eigen::MatrixXd& matr) override;
		virtual void FillIQFromOther(const Boundary& otherBoundary, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ) override;
		virtual void SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol) override;

		virtual void CalcConvVelocityToSetOfPointsFromSheets(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const override;
#if defined(USE_CUDA)
		virtual void GPUCalcConvVelocityToSetOfPointsFromSheets(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const override;
#endif

		virtual void CalcConvVelocityAtVirtualVortexes(std::vector<Point2D>& velo) const override;

		virtual void ComputeAttachedSheetsIntensity() override;		
	
		virtual void GetInfluenceFromVorticesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const override;
		virtual void GetInfluenceFromVorticesToCurvPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const override {};

		virtual void GetInfluenceFromSourcesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const override;
		virtual void GetInfluenceFromSourcesToCurvPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const override {};

		virtual void GetInfluenceFromVortexSheetAtRectPanelToVortex(size_t panel, const Vortex2D& vtx, Point2D& vel) const override;
		virtual void GetInfluenceFromSourceSheetAtRectPanelToVortex(size_t panel, const Vortex2D& vtx, Point2D& vel) const override;
	
		virtual void GetInfluenceFromVInfToRectPanel(std::vector<double>& vInfRhs) const override;
		virtual void GetInfluenceFromVInfToCurvPanel(std::vector<double>& vInfRhs) const override;

        };

}//namespace VM2D

#endif


