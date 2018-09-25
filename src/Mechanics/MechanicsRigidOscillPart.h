/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MechanicsRigidOscillPart.h                                       |
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
\brief Заголовочный файл с описанием класса MechanicsRigidOscillPart
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/

#ifndef MECHANICSRIGIDOSCILLPART_H
#define MECHANICSRIGIDOSCILLPART_H

#include "Mechanics.h"

/*!
\brief Класс, определяющий вид механической системы

Упруго закрепленное тело, метод расщепления

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.1
\date 2 апреля 2018 г.
*/

class MechanicsRigidOscillPart :
	public Mechanics
{
private:
	//начальная скорость и отклонение
	const double u0;
	const double y0;

	//плотность и площадь круга (радиус 0.5) для вычисления массы
//	const double rhoAfl = 30.0;
//	const double SAfl = 0.7853981633974483;
	const double m;
	const double b;
	const double k;

	//текущие скорость и отклонение
	double u, y;
public:

	/// \brief Конструктор
	/// 
	/// \param[in] passport_ константная ссылка на паспорт
	/// \param[in] afl_ ссылка на профиль;
	/// \param[in] boundary_ константная ссылка на граничное условие;
	/// \param[in] virtVortParams_ константная ссылка на параметры виртуального вихревого следа для профиля;
	/// \param[in] parallel_ константная ссылка на параметры параллельного исполнения.
	MechanicsRigidOscillPart(const Passport& passport_, Airfoil& afl_, const Boundary& boundary_, const VortexesParams& virtVortParams_, const Parallel& parallel_) :
		Mechanics(passport_, afl_, boundary_, virtVortParams_, parallel_, 0, true, false), u0(0.0), y0(0.0), m(39.15), b(0.731), k(39.15 * 4.0 * PI * PI * passport_.param * passport_.param * 3.0 * 3.0)
	{
		u = u0;
		y = y0;
	};

	/// Деструктор
	~MechanicsRigidOscillPart() {};

	
	//далее -- реализации виртуальных функций
	virtual void GetHydroDynamForce(timePeriod& time);
	virtual Point2D VeloOfAirfoilRcm(double currTime);
	virtual Point2D PositionOfAirfoilRcm(double currTime);
	virtual void VeloOfAirfoilPanels(double currTime);

	//TODO реализовать
	virtual void FillMechanicsRowsAndCols(Eigen::MatrixXd& row, Eigen::MatrixXd& col) {};
	virtual void FillMechanicsRhs(std::vector<double>& rhs) {};
	virtual void Move();
	
};

#endif