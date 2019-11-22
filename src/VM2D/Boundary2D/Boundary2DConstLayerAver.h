/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Boundary2DConstLayerAver.h                                       |
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
\brief Заголовочный файл с описанием класса BoundaryConstLayerAver
\warning Пока реализованы не все методы
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7   
\date 22 ноября 2019 г.
*/

#ifndef BOUNDARYCONSTLAYERAVER_H
#define BOUNDARYCONSTLAYERAVER_H

#include "Boundary2D.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс, определяющий способ удовлетворения граничного условия на обтекаемом профиле

	Способ удовлетворения граничного условия:
	- генерация вихревого слоя постоянной интенсивности на панелях;
	- условие ортогональности невязки граничного условия константе (выполнение граничного условия в среднем по панели).

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.7
	\date 22 ноября 2019 г.
	*/
	class BoundaryConstLayerAver : public Boundary
	{

	public:

		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] numberInPassport_ номер профиля в паспорте задачи
		BoundaryConstLayerAver(const World2D& W_, size_t numberInPassport_) :
			Boundary(W_, numberInPassport_, 1) {};

		/// Деструктор
		virtual ~BoundaryConstLayerAver() {};

		//далее -- реализации виртуальных функций
		virtual void FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol) override;
		virtual void FillIQSelf(std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ) override;
		virtual void FillMatrixFromOther(const Boundary& otherBoundary, Eigen::MatrixXd& matr) override;
		virtual void FillIQFromOther(const Boundary& otherBoundary, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& IQ) override;
		virtual void SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol) override;

		virtual void GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const override;
		virtual void GetConvVelocityToSetOfPointsFromVirtualVortexes(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const override;
#if defined(USE_CUDA)
		virtual void GPUGetConvVelocityToSetOfPointsFromVirtualVortexes(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const override;
#endif

		virtual void GetConvVelocityAtVirtualVortexes(std::vector<Point2D>& velo) const override;

		virtual void ComputeAttachedSheetsIntensity() override;
		virtual void FillRhsFromOther(const Airfoil& otherAirfoil, Eigen::VectorXd& rhs) override;

		virtual void GetInfluenceFromVorticesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const override;
		virtual void GetInfluenceFromVorticesToCurvPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const override {};
		
		virtual void GetInfluenceFromSourcesToRectPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const override;
		virtual void GetInfluenceFromSourcesToCurvPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& wakeRhs) const override {};


	};

	   	 	   
}//namespace VM2D

#endif
