/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.5    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/02/20     |
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
\version 1.5   
\date 20 февраля 2019 г.
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

	\version 1.5
	\date 20 февраля 2019 г.
	*/
	class BoundaryConstLayerAver : public Boundary
	{
	private:
		/// \brief Вспомогательная функция вычисления угла между векторами
		///
		/// \param[in] p константная ссылка на первый вектор
		/// \param[in] s константная ссылка на второй вектор
		/// \return угол между векторами в диапазоне \f$ (-\pi; \pi] \f$
		double Alpha(const Point2D& p, const Point2D& s) const
		{
			return atan2(cross3(p, s), p*s);
		}

		/// \brief Вспомогательная функция вычисления логарифма отношения норм векторов
		///
		/// \param[in] p константная ссылка на первый вектор
		/// \param[in] s константная ссылка на второй вектор
		/// \return логарифм отношения норм векторов
		double Lambda(const Point2D& p, const Point2D& s) const
		{
			return 0.5*log((s*s) / (p*p));
		}

		/// \brief Вспомогательная функция вычисления величины \f$ (\vec a \cdot \vec b) \cdot \vec c + (\vec a \times \vec b) \times \vec c \f$
		///
		/// Для оптимизации все векторы считаются двумерными
		///
		/// \param[in] a константная ссылка на первый вектор
		/// \param[in] b константная ссылка на второй вектор
		/// \param[in] c константная ссылка на третий вектор
		/// \return логарифм отношения норм векторов
		Point2D Omega(const Point2D& a, const Point2D& b, const Point2D& c) const
		{
			return (a * b) * c + (Point2D({ -c[1], c[0] }))*cross3(a, b);
		}

	public:

		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] numberInPassport_ номер профиля в паспорте задачи
		BoundaryConstLayerAver(const World2D& W_, size_t numberInPassport_) :
			Boundary(W_, numberInPassport_, 1)
		{ };

		/// Деструктор
		virtual ~BoundaryConstLayerAver() {};

		//далее -- реализации виртуальных функций
		virtual void FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol);
		virtual void FillMatrixFromOther(const Boundary& otherBoundary, Eigen::MatrixXd& matr);
		virtual void FillRhs(const Point2D& V0, Eigen::VectorXd& rhs, double* lastRhs, bool move, bool deform);
		virtual size_t GetUnknownsSize() const;
		virtual void SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol);

		virtual void GetWakeInfluence(std::vector<double>& wakeVelo) const;
#if defined(USE_CUDA)
		virtual void GPUGetWakeInfluence(std::vector<double>& wakeVelo) const;
#endif

		virtual void GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const;
		virtual void GetConvVelocityToSetOfPointsFromVirtualVortexes(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const;
#if defined(USE_CUDA)
		virtual void GPUGetConvVelocityToSetOfPointsFromVirtualVortexes(const WakeDataBase& pointsDb, std::vector<Point2D>& velo) const;
#endif

		virtual void GetConvVelocityAtVirtualVortexes(std::vector<Point2D>& velo) const;

		virtual void ComputeAttachedSheetsIntensity();
		virtual void FillRhsFromOther(const Airfoil& otherAirfoil, Eigen::VectorXd& rhs);
	};

}//namespace VM2D

#endif
