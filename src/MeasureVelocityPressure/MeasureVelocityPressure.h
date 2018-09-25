/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.3    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/09/26     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: MeasureVelocityPressure.h                                        |
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
\brief Заголовочный файл с описанием класса MeasureVelocityPressure
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.3
\date 26 сентября 2018 г.
*/

#ifndef MEASUREVELOCITYPRESSURE_H
#define MEASUREVELOCITYPRESSURE_H

#include <memory>
#include <vector>
#include "Eigen/Dense"

#include "defs.h"
//#include "Airfoil.h"
//#include "Boundary.h"
#include "Point2D.h"
#include "StreamParser.h"
#include "Times.h"
#include "WakeDataBase.h"
//#include "Velocity.h"

/*!
\brief Класс, отвечающий за вычисление поля скорости и давления в заданых точках для вывода

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.3
\date 26 сентября 2018 г.
*/

class MeasureVelocityPressure
{
private:
	///количество точек
	size_t nPts;

	/// Точки, которые считываются из файла
	std::vector<Point2D> initialPoints;

	/// Точки, в которых нужно вычислять в данный момент времени (исключены точки внутри профиля)
	//записываем их в мнимую систему вихрей, чтобы воспользоваться методом CalcConvVeloToSetOfPoints(...)
	WakeDataBase additionalWake;

	/// Скорости в нужных точках
	std::vector<Point2D> velocity;

	/// Радиусы вихрей в нужных точках
	std::vector<double> domainRadius;

	/// Давление в нужных точках
	std::vector<double> pressure;

	///???
	const double pInf = 1.0;

protected:
	/// Константная ссылка на решаемую задачу
	const World2D& W;

public:
	/// \brief Конструктор
	/// 
	/// \param[in] W_ константная ссылка на решаемую задачу
	MeasureVelocityPressure(const World2D& W_) : W(W_) {};

	/// Деструктор
	~MeasureVelocityPressure() { };

	/// \brief Чтение точек, в которых нужно посчитать давление и скорость
	///
	/// \param[in] dir константная ссылка на строку --- имя каталога, где лежит cчитываемый файл
	void ReadPointsFromFile(const std::string& dir);

	/// \brief Расчет и сохранение в файл поля скоростей и давления
	///
	/// \param[in] dir константная ссылка на строку, задающую каталог, куда сохранять файл с вихревым следом
	/// \param[in] step номер кадра для сохранения
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void CalcSaveVelocityPressure(const std::string& dir, size_t step, timePeriod& time);

};

#endif