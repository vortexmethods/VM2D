/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DDeformable.h                                          |
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
\brief Заголовочный файл с описанием класса MechanicsDeformable
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#ifndef MECHANICSDEFORMABLE_H
#define MECHANICSDEFORMABLE_H

#include "Mechanics2D.h"

namespace VM2D
{

	class World2D;
	class Beam
	{
		double x0;
		double L;

	public:
		Beam(double x0_, double L_):		
			x0(x0_),
			L(L_)
		{};
		double getDisp(double x, double t) 
		{ 
			return 0.2 * 2.84 * (x - x0) * (x - x0) * sin(DPI * t); 
		}
	};


	struct ChordPanel
	{
		Point2D beg, end;
		std::pair<size_t, size_t> infPanels;
		double rightWidth;
	};



	/*!
	\brief Класс, определяющий вид механической системы

	Деформируемое

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
\author Колганова Александра Олеговна

	\Version 1.12
	\date 14 января 2024 г.
	*/

	class MechanicsDeformable :
		public Mechanics
	{
	private:

		

	public:
		/// текущая скорость профиля
		Point2D& getVcm() { return Vcm; };		

		/// текущее отклонение профиля
		Point2D& getRcm() { return Rcm; };

		/// текущая угловая скорость профиля
		double& getWcm() { return Wcm; };

		/// текущий угол поворота профиля
		double& getPhicm() { return Phi; };

		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] numberInPassport_ номер профиля в паспорте задачи	
		MechanicsDeformable(const World2D& W_, size_t numberInPassport_);

		/// Деструктор
		~MechanicsDeformable() {};

		//далее -- реализации виртуальных функций
		virtual void GetHydroDynamForce() override;
		virtual Point2D VeloOfAirfoilRcm(double currTime) override;
		virtual Point2D PositionOfAirfoilRcm(double currTime) override;
		virtual double AngularVelocityOfAirfoil(double currTime) override;
		virtual double AngleOfAirfoil(double currTime) override;
		virtual void VeloOfAirfoilPanels(double currTime) override;
		virtual void ReadSpecificParametersFromDictionary() override;
		virtual void Move() override;

		/// хорда
		size_t indexOfUpperRightAngle;
		size_t indexOfUpperLeftAngle;
		size_t indexOfLowerRightAngle;
		size_t indexOfLowerLeftAngle;

		std::vector<ChordPanel> chord;
		std::vector<double> upperShifts;
		std::vector<double> lowerShifts;

		std::unique_ptr<Beam> beam;

	};

}//namespace VM2D

#endif