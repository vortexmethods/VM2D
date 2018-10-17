/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MechanicsRigidRotateMon.h                                        |
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
\brief Заголовочный файл с описанием класса MechanicsRigidRotateMon
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#ifndef MECHANICSRIGIDROTATEMON_H
#define MECHANICSRIGIDROTATEMON_H

#include "Mechanics.h"

class World2D;

/*!
\brief Класс, определяющий вид механической системы

Жесткое тело, движущееся по заданному закону

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.4
\date 16 октября 2018 г.
*/

class MechanicsRigidRotateMon :
	public Mechanics
{
private:
	//момент инерции тела и жесткость пружины
	double J;
	double b;
	double k;

	///время, в течение которого ротор крутится без нагрузки
	double tRotateAccel;

	///время, в течение которого постепенно увеличивается момент, которым нагружается ротор
	double tMomentAccel;

	///Момент, который снимается с ротора
	double Mz;

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

	double GetMz(double t)
	{
		if (t > (tMomentAccel + tRotateAccel))
			return Mz;
		else if (t < tRotateAccel)
			return 0.0;
		else return ((Mz / tMomentAccel) * t - (Mz * tRotateAccel / tMomentAccel));
	}


public:

	/// \brief Конструктор
	/// 
	/// \param[in] W_ константная ссылка на решаемую задачу
	/// \param[in] numberInPassport_ номер профиля в паспорте задачи	
	MechanicsRigidRotateMon(const World2D& W_, size_t numberInPassport_)
		: Mechanics(W_, numberInPassport_, 1, true, false, true, { 0.0, 0.0 }, {0.0, 0.0}, 0.0, 0.0), b(0.0 * 0.731)
	{
		ReadSpecificParametersFromDictionary();
	};

	/// Деструктор
	~MechanicsRigidRotateMon() {};


	//далее -- реализации виртуальных функций
	virtual void GetHydroDynamForce(timePeriod& time);
	virtual Point2D VeloOfAirfoilRcm(double currTime);
	virtual Point2D PositionOfAirfoilRcm(double currTime);
	virtual void VeloOfAirfoilPanels(double currTime);
	virtual void ReadSpecificParametersFromDictionary();

	virtual void FillMechanicsRowsAndCross(Eigen::MatrixXd& row, Eigen::MatrixXd& cross);
	virtual void FillMechanicsRhs(std::vector<double>& rhs);
	virtual void FillAtt(Eigen::MatrixXd& row, Eigen::MatrixXd& rhs);
	virtual void SolutionToMechanicalSystem(Eigen::VectorXd& col);
	virtual void Move();

};

#endif