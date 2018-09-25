/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.3    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/09/26     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MechanicsRigidOscillStronglyCoupled.h                            |
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
\brief Заголовочный файл с описанием класса MechanicsRigidOscillStronglyCoupled
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/

#ifndef MECHANICSRIGIDOSCILLSTRONGLYCOUPLED_H
#define MECHANICSRIGIDOSCILLSTRONGLYCOUPLED_H

#include "Mechanics.h"

class World2D;

/*!
\brief Класс, определяющий вид механической системы

Жесткое тело, движущееся по заданному закону

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.3
\date 26 сентября 2018 г.
*/

class MechanicsRigidOscillStronglyCoupled :
	public Mechanics
{
private:
	//плотность и площадь круга (радиус 0.5) для вычисления массы
	//	const double rhoAfl = 30.0;
	//	const double SAfl = 0.7853981633974483;
	double m;
	const double b;

	double k;
	
	
public:

	/// \brief Конструктор
	/// 
	/// \param[in] W_ константная ссылка на решаемую задачу
	/// \param[in] numberInPassport_ номер профиля в паспорте задачи	
	MechanicsRigidOscillStronglyCoupled(const World2D& W_, size_t numberInPassport_)
		: Mechanics(W_, numberInPassport_, 1, true, false, false, { 0.0, 0.0 }, { 0.0, 0.0 }, 0.0, 0.0), m(39.15), b(0.0 * 0.731)
	{
		
		ReadSpecificParametersFromDictionary();
	};

	/// Деструктор
	~MechanicsRigidOscillStronglyCoupled() {};


	//далее -- реализации виртуальных функций
	virtual void GetHydroDynamForce(timePeriod& time);
	virtual Point2D VeloOfAirfoilRcm(double currTime);
	virtual Point2D PositionOfAirfoilRcm(double currTime);
	virtual void VeloOfAirfoilPanels(double currTime);
	virtual void ReadSpecificParametersFromDictionary();

	/// \todo реализовать
	virtual void FillMechanicsRowsAndCross(Eigen::MatrixXd& row, Eigen::MatrixXd& cross);
	virtual void FillMechanicsRhs(std::vector<double>& rhs);
	virtual void FillAtt(Eigen::MatrixXd& row, Eigen::MatrixXd& rhs);
	virtual void SolutionToMechanicalSystem(Eigen::VectorXd& col);
	virtual void Move();

};

#endif