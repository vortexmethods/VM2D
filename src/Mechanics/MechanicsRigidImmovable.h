/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MechanicsRigidImmovable.h                                        |
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
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/

#ifndef MECHANICSRIGIDIMMOVABLE_H
#define MECHANICSRIGIDIMMOVABLE_H

#include "Mechanics.h"

/*!
\brief Класс, определяющий вид механической системы

Жесткое неподвижное тело

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.1
\date 2 апреля 2018 г.
*/

class MechanicsRigidImmovable :
	public Mechanics
{
public:
	/// \brief Конструктор
	/// 
	/// \param[in] passport_ константная ссылка на паспорт
	/// \param[in] afl_ ссылка на профиль;
	/// \param[in] boundary_ константная ссылка на граничное условие;
	/// \param[in] virtVortParams_ константная ссылка на параметры виртуального вихревого следа для профиля;
	/// \param[in] parallel_ константная ссылка на параметры параллельного исполнения.
	MechanicsRigidImmovable(const Passport& passport_, Airfoil& afl_, const Boundary& boundary_, const VortexesParams& virtVortParams_, const Parallel& parallel_) :
		Mechanics(passport_, afl_, boundary_, virtVortParams_, parallel_, 0, false, false)
	{};

	/// Деструктор
	~MechanicsRigidImmovable() {};


	//далее -- реализации виртуальных функций
	virtual void GetHydroDynamForce(timePeriod& time);
	virtual Point2D VeloOfAirfoilRcm(double currTime);
	virtual Point2D PositionOfAirfoilRcm(double currTime);
	virtual void VeloOfAirfoilPanels(double currTime);

	//TODO реализовать
	virtual void FillMechanicsRowsAndCols(Eigen::MatrixXd& row, Eigen::MatrixXd& col) {};
	virtual void FillMechanicsRhs(std::vector<double>& rhs) {};
	virtual void ComputeDisplacementRcm() {};
	virtual void Move() {};
};

#endif