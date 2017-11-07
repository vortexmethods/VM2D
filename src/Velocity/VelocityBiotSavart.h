/*!
\file
\brief Заголовочный файл с описанием класса VelocityBiotSavart
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef VELOCITYBIOTSAVART_H
#define VELOCITYBIOTSAVART_H

#include "Velocity.h"

/*!
\brief Класс, определяющий способ вычисления скоростей

Способ вычисления скоростей
- напрямую по закону Био --- Савара

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.0
\date 1 декабря 2017 г.
*/
class VelocityBiotSavart : public Velocity
{
public:
	/// \brief Конструктор
	/// 
	/// \param[in] parallel_ константная ссылка на параметры исполнения задачи в параллельном MPI-режиме
	/// \param[in] wake_ константная ссылка на вихревой след	
	/// \param[in] boundary_ константная ссылка на вектор указателей на граничные условия 	
	VelocityBiotSavart(const Parallel& parallel_, const Wake& wake_, const std::vector<std::unique_ptr<Boundary>>& boundary_) :
		Velocity(parallel_, wake_, boundary_)
	{ };

	/// Деструктор
	virtual ~VelocityBiotSavart(){};

	//реализация виртуальной функции
	virtual void CalcConvVeloToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo, std::vector<double>& domainRadius);
	virtual void CalcDiffVeloToSetOfPoints(const std::vector<Vortex2D>& points, const std::vector<double>& domainRadius, const std::vector<Vortex2D>& vortices, std::vector<Point2D>& velo);
};

#endif