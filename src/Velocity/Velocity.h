/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Velocity.h                                                       |
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
\brief Заголовочный файл с описанием класса Velocity
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#ifndef VELOCITY_H
#define VELOCITY_H

#include <vector>

#include "Point2D.h"
#include "WakeDataBase.h"
#include "Boundary.h"

//#include "Gpu.h"

class World2D;

/*!
\brief Структура, определяющая параметры виртуальных вихрей для отдельного профиля
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
struct VortexesParams
{
	/// Вектор конвективных скоростей вихрей
	std::vector<Point2D> convVelo;

	/// Вектор диффузионных скоростей вихрей
	std::vector<Point2D> diffVelo;

	/// Вектор знаменателей (I0) диффузионных скоростей вихрей (обусловленных профилем)
	std::vector<double> I0;

	/// Вектор знаменателей (I1) диффузионных скоростей вихрей (обусловленных завихренностью)
	std::vector<double> I1;

	/// Вектор числителей (I2) диффузионных скоростей вихрей (обусловленных завихренностью)
	std::vector<Point2D> I2;

	/// Вектор числителей (I3) диффузионных скоростей вихрей (обусловленных профилем)
	std::vector<Point2D> I3;

	/// Вектор характерных радиусов вихревых доменов (eps*)
	std::vector<double> epsastWake;
};

/*!
\brief Абстрактный класс, определяющий способ вычисления скоростей
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
class Velocity
{
protected:
	/// Константная ссылка на решаемую задачу
	const World2D& W;	
	
public:
	/// Струтура, определяющая параметры вихрей в следе
	VortexesParams wakeVortexesParams;

	/// Вектор струтур, определяющий параметры виртуальных вихрей для профилей
	std::vector<VortexesParams> virtualVortexesParams;
	

	/// \brief Конструктор
	/// 
	/// \param[in] W_ константная ссылка на решаемую задачу 	
	Velocity(const World2D& W_) :
		W(W_)
	{			
		virtualVortexesParams.resize(0);
	};
	
	/// \brief Вычисление конвективных скоростей и радиусов вихревых доменов в заданном наборе точек от следа
	///
	/// \param[in] pointsDb константная ссылка на базу данных пелены из вихрей, в которых надо сосчитать конвективные скорости
	/// \param[out] velo ссылка на вектор скоростей в требуемых точках
	/// \param[out] domainRadius ссылка на вектор радиусов вихревых доменов
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	virtual void CalcConvVeloToSetOfPoints(const WakeDataBase& pointsDb, std::vector<Point2D>& velo,  std::vector<double>& domainRadius) = 0;
#if defined(USE_CUDA)	
	virtual void GPUCalcConvVeloToSetOfPoints(const WakeDataBase& pointsDb, std::vector<Point2D>& velo, std::vector<double>& domainRadius) = 0;
#endif

	/// \brief Вычисление конвективных скоростей вихрей и виртуальных вихрей в вихревом следе
	///
	/// Вызывает 2 раза функцию CalcConvVeloToSetOfPoints
	///
	/// \warning скорости приплюсовываются к уже имеющимся
	void CalcConvVelo();

	/// \brief Вычисление числителей и знаменателей диффузионных скоростей в заданном наборе точек
	///
	/// \param[in] pointsDb константная ссылка на базу данных пелены из вихрей, в которых надо сосчитать диффузионные скорости
	/// \param[in] domainRadius константная ссылка на вектор радиусов вихревых доменов
	/// \param[in] vorticesDb константная ссылка на на базу данных пелены из вихрей,от которых надо сосчитать влияния на points
	/// \param[out] I1 ссылка на вектор величин I1 (знаменателей в диффузионных скоростях) в требуемых точках
	/// \param[out] I2 ссылка на вектор величин I2 (числителей в диффузионных скоростях) в требуемых точках
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	virtual void CalcDiffVeloI1I2ToSetOfPoints(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2) = 0;
	virtual void CalcDiffVeloI1I2ToSetOfPointsFromPanels(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2) = 0;
	
#if defined(USE_CUDA)
	virtual void GPUCalcDiffVeloI1I2ToSetOfPoints(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const WakeDataBase& vorticesDb, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh = false) = 0;
	virtual void GPUCalcDiffVeloI1I2ToSetOfPointsFromPanels(const WakeDataBase& pointsDb, const std::vector<double>& domainRadius, const Boundary& bnd, std::vector<double>& I1, std::vector<Point2D>& I2, bool useMesh = false) = 0;
#endif

	/// \brief Вычисление диффузионных скоростей вихрей и виртуальных вихрей в вихревом следе
	///
	/// Вызывает 4 раза функцию CalcDiffVeloToSetOfPoints
	///
	/// \warning скорости приплюсовываются к уже имеющимся
	void CalcDiffVelo();

	/// Деструктор
	virtual ~Velocity() { };
};

#endif
