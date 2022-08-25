/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.11   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2022/08/07     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2022 Ilia Marchevsky, Kseniia Sokol, Evgeniya Ryatina    |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DRigidGivenLaw.h                                       |
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
\brief Заголовочный файл с описанием класса MechanicsRigidGivenLaw
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.11
\date 07 августа 2022 г.
*/

#ifndef MECHANICSRIGIDGIVENLAW_H
#define MECHANICSRIGIDGIVENLAW_H

#include <functional>

#include "Mechanics2D.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс, определяющий вид механической системы

	Жесткое тело, движущееся по заданному закону

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.11
	\date 07 августа 2022 г.
	*/


	class MechanicsRigidGivenLaw :
		public Mechanics
	{
	private:

		const double timeOfMovement = 5.0;
		const double refVel = 2.0;
		//const Point2D maxDisplacementOfCM = { 2.0, 3.0 };
		//const double maxRotationAngle = PI / 3.0;

		std::function<Point2D(double)> PositionOfCenterOfMass = [=](double t) -> Point2D
		{
			if (numberInPassport < 4)
			{

				if (t < timeOfMovement)
					return { -25.0 + 0.5 * (t * t / timeOfMovement) * refVel - 1.04 * numberInPassport, 0.0 };
				else
					return { -25.0 + (t - 0.5 * timeOfMovement) * refVel - 1.04 * numberInPassport, 0.0 };
			}
			else
			{

				if (t < timeOfMovement)
					return { 25.0 - 0.5 * (t * t / timeOfMovement) * refVel + 1.04 * (numberInPassport - 4), -0.164 };
				else
					return { 25.0 - (t - 0.5 * timeOfMovement) * refVel + 1.04 * (numberInPassport - 4), -0.164 };
			}

			

		};

		std::function<Point2D(double)> VelocityOfCenterOfMass = [=](double t) -> Point2D
		{
			if (numberInPassport < 4)
			{
				if (t < timeOfMovement)
					return { (t / timeOfMovement) * refVel, 0.0 };
				else
					return { refVel, 0.0 };
			}
			else
			{
				if (t < timeOfMovement)
					return { -(t / timeOfMovement) * refVel, 0.0 };
				else
					return { -refVel, 0.0 };
			}

		};

		std::function<double(double)> RotationAngle = [=](double t) -> double
		{
			return 0.0;
		};

		std::function<double(double)> AngularVelocity = [=](double t) -> double
		{
			return 0.0;
		};

		
	public:
		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] numberInPassport_ номер профиля в паспорте задачи	
		MechanicsRigidGivenLaw(const World2D& W_, size_t numberInPassport_);


		/// Деструктор
		~MechanicsRigidGivenLaw() { };
					
		

		//далее -- реализации виртуальных функций
		virtual void GetHydroDynamForce() override;
		virtual Point2D VeloOfAirfoilRcm(double currTime) override;
		virtual Point2D PositionOfAirfoilRcm(double currTime) override;
		virtual double AngularVelocityOfAirfoil(double currTime) override;
		virtual double AngleOfAirfoil(double currTime) override;
		virtual void VeloOfAirfoilPanels(double currTime) override;
		virtual void ReadSpecificParametersFromDictionary() override;

		virtual void Move() override;
	};

}//namespace VM2D

#endif