/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Mechanics2DRigidOscillPart.h                                     |
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
\brief Заголовочный файл с описанием класса MechanicsRigidOscillPart
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.7   
\date 22 ноября 2019 г.
*/

#ifndef MECHANICSRIGIDOSCILLPART_H
#define MECHANICSRIGIDOSCILLPART_H

#include "Mechanics2D.h"

namespace VM2D
{

	class World2D;

	/*!
	\brief Класс, определяющий вид механической системы

	Упруго закрепленное тело, метод расщепления

	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна

	\version 1.7
	\date 22 ноября 2019 г.
	*/

	class MechanicsRigidOscillPart :
		public Mechanics
	{
	private:

		/// начальная скорость профиля
		const double u0;

		/// начальное отклонение профиля
		const double y0;

		/// текущая скорость профиля
		double u;

		/// текущее отклонение профиля
		double y;

		/// скорость профиля с предыдущего шага
		double uOld;

		/// отклонение профиля с предыдущего шага
		double yOld;

		/// масса профиля
		double m;

		/// параметр демпфирования механической системы
		const double b;

		/// параметр жесткости механической системы
		double k;

	public:

		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] numberInPassport_ номер профиля в паспорте задачи	
		MechanicsRigidOscillPart(const World2D& W_, size_t numberInPassport_)
			: Mechanics(W_, numberInPassport_, 0, true, false, false, { 0.0, 0.0 }, { 0.0, 0.0 }, 0.0, 0.0), u0(0.0), y0(0.0), b(0.0*0.731)
		{
			u = u0;
			y = y0;
			uOld = u0;
			yOld = y0;
			ReadSpecificParametersFromDictionary();
		};

		/// Деструктор
		~MechanicsRigidOscillPart() {};


		//далее -- реализации виртуальных функций
		virtual void GetHydroDynamForce(timePeriod& time) override;
		virtual Point2D VeloOfAirfoilRcm(double currTime) override;
		virtual Point2D PositionOfAirfoilRcm(double currTime) override;
		virtual void VeloOfAirfoilPanels(double currTime) override;
		virtual void ReadSpecificParametersFromDictionary() override;
		virtual void Move() override;
	};

}//namespace VM2D

#endif