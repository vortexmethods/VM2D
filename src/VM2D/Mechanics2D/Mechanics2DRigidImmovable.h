/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DRigidImmovable.h                                      |
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
\brief Заголовочный файл с описанием класса MechanicsRigidImmovable
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#ifndef MECHANICSRIGIDIMMOVABLE_H
#define MECHANICSRIGIDIMMOVABLE_H

#include "Mechanics2D.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс, определяющий вид механической системы

	Жесткое неподвижное тело

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
\author Колганова Александра Олеговна

	\Version 1.12
	\date 14 января 2024 г.
	*/

	class MechanicsRigidImmovable :
		public Mechanics
	{
	public:
		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] numberInPassport_ номер профиля в паспорте задачи	
		MechanicsRigidImmovable(const World2D& W_, size_t numberInPassport_);

		/// Деструктор
		~MechanicsRigidImmovable() {};


		//далее -- реализации виртуальных функций
		virtual void GetHydroDynamForce() override;
		virtual Point2D VeloOfAirfoilRcm(double currTime) override;
		virtual Point2D PositionOfAirfoilRcm(double currTime) override;
		virtual double AngularVelocityOfAirfoil(double currTime) override;
		virtual double AngleOfAirfoil(double currTime) override;

		virtual void VeloOfAirfoilPanels(double currTime) override;
		virtual void ReadSpecificParametersFromDictionary() override {};

		virtual void Move() override {};
	};

}//namespace VM2D

#endif