/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.10   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2021/05/17     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2021 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Airfoil2DRect.h                                                  |
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
\brief Заголовочный файл с описанием класса AirfoilRect
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.10
\date 17 мая 2021 г.
*/

#ifndef AIRFOILRECT_H
#define AIRFOILRECT_H

#include "Airfoil2D.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс, определяющий тип обтекаемого профиля

	Тип профиля:
	- профиль с прямолинейными панелями.

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.10
	\date 17 мая 2021 г.
	*/
	class AirfoilRect
		: public Airfoil
	{
	public:

		/// Конструктор
		AirfoilRect(const World2D& W_, const size_t numberInPassport_)
			:Airfoil(W_, numberInPassport_)
		{ };

		AirfoilRect(const Airfoil& afl) : Airfoil(afl) {};

		/// Деструктор
		virtual ~AirfoilRect() { };
		
		/// Вычисление нормалей, касательных и длин панелей по текущему положению вершин
		void CalcNrmTauLen();

		///Вычисляет габаритный прямоугольник профиля
		virtual void GetGabarits(double gap = 0.02) override;

		//далее -- реализация виртуальных функций
		virtual void ReadFromFile(const std::string& dir) override;
		virtual void Rotate(double alpha) override;
		virtual void Scale(double factor) override;
		virtual void Move(const Point2D& dr) override;
		
		virtual std::vector<double> getA(size_t p, size_t i, const Airfoil& airfoil, size_t j) const override;
		virtual void calcIQ(size_t p, const Airfoil& otherAirfoil, std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& matrPair) const override;

		virtual bool IsPointInAirfoil(const Point2D& point) const override;

		virtual void GetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3) override;
#if defined(USE_CUDA)
		virtual void GPUGetDiffVelocityI0I3ToSetOfPointsAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3) override;
#endif
		virtual void GetDiffVelocityI0I3ToWakeAndViscousStresses(const WakeDataBase& pointsDb, std::vector<double>& domainRadius, std::vector<double>& I0, std::vector<Point2D>& I3) override;

		virtual void GetInfluenceFromVorticesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const override;
		virtual void GetInfluenceFromSourcesToPanel(size_t panel, const Vortex2D* ptr, ptrdiff_t count, std::vector<double>& panelRhs) const override;

		virtual void GetInfluenceFromSourceSheetToVortex(size_t panel, const Vortex2D& vtx, Point2D& vel) const override;
		virtual void GetInfluenceFromVortexSheetToVortex(size_t panel, const Vortex2D& vtx, Point2D& vel) const override;

		virtual void GetInfluenceFromVInfToPanel(std::vector<double>& vInfRhs) const override;
	};

}//namespace VM2D

#endif