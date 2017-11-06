/*!
\file
\brief Заголовочный файл с описанием класса VelocityDirect
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef VELOCITYDIRECT_H
#define VELOCITYDIRECT_H

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
class VelocityDirect : public Velocity
{
public:
	/// \brief Конструктор
	/// 
	/// \param[in] parallel_ константная ссылка на параметры исполнения задачи в параллельном MPI-режиме
	/// \param[in] wake_ константная ссылка на вихревой след	
	VelocityDirect(const Parallel& parallel_, const Wake& wake_) :
		Velocity(parallel_, wake_)
	{ };

	/// Деструктор
	virtual ~VelocityDirect(){};

	//реализация виртуальной функции
	virtual void CalcConvVelo(double dt = 1.0);	
};

#endif