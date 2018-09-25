/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.2    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/06/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: BoundaryVortColl.h                                               |
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
\brief Заголовочный файл с описанием класса BoundaryVortColl
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.2
\date 14 июня 2018 г.
*/

#ifndef BOUNDARYVORTCOLL_H
#define BOUNDARYVORTCOLL_H

#include "Boundary.h"

class World2D;

/*!
\brief Класс, определяющий способ удовлетворения граничного условия на обтекаемом профиле

Способ удовлетворения граничного условия:
- рождение вихрей в центрах панелей;
- условие коллокации в центрах панелей для касательной составляющей скорости.

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.2
\date 14 июня 2018 г.
*/
class BoundaryVortColl : public Boundary
{
private:
	/// Контрольные точки (устанавливаются в центры панелей)
	std::vector<Point2D> KK;

	/// \brief Вычисляет скос на точку R от точки X
	///
	/// Вектор вихревого влияния на точку R от вихря единичной интенсивности, находящегося в точке X
	/// \n (сглаживание поля скоростей не производится)
	/// 
	/// \param[in] R точка наблюдения
	/// \param[in] X точка, где находится вихрь
	/// \return вектор скорости
	Point2D Skos(const Point2D& R, const Point2D& X);

public:
	/// \brief Конструктор
	/// 
	/// \param[in] W_ константная ссылка на решаемую задачу
	/// \param[in] numberInPassport_ номер профиля в паспорте задачи
	BoundaryVortColl(const World2D& W_, size_t numberInPassport_);

	/// Деструктор
	virtual ~BoundaryVortColl() { };

	//далее -- реализации виртуальных функций
	virtual void FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol);
	
	virtual void GetWakeInfluence(std::vector<double>& wakeVelo) const;
#if defined(USE_CUDA)
	virtual void GPUGetWakeInfluence(std::vector<double>& wakeVelo) const;
#endif

	virtual void GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const;
	virtual void FillRhs(const Point2D& V0, Eigen::VectorXd& rhs, double* lastRhs, bool move, bool deform);
	virtual size_t GetUnknownsSize() const;
	virtual void SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol);
	virtual void ComputeAttachedSheetsIntensity();
	virtual void FillRhsFromOther(const Airfoil& otherAirfoil, Eigen::VectorXd& rhs);
};

#endif
