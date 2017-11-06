/*!
\file
\brief Заголовочный файл с описанием класса Velocity
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef VELOCITY_H
#define VELOCITY_H

#include "Wake.h"

/*!
\brief Абстрактный класс, определяющий способ вычисления скоростей
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class Velocity
{
protected:
	/// Константная ссылка на вихревой след
	const Wake& wake;

	/// Константная ссылка на параметры исполнения задачи в параллельном MPI-режиме
	const Parallel& parallel;	

public:
	/// Вектор конвективных скоростей вихрей в следе
	std::vector<Point2D> convVelo;
	
	/// \brief Конструктор
	/// 
	/// \param[in] parallel_ константная ссылка на параметры исполнения задачи в параллельном MPI-режиме
	/// \param[in] wake_ константная ссылка на вихревой след	
	Velocity(const Parallel& parallel_, const Wake& wake_) :
		wake(wake_), parallel(parallel_)
	{		
		convVelo.resize(0);
	};
	
	/// \brief Вычисление конвективных скоростей вихрей в вихревом следе
	///
	/// \param[in] dt величина шага по времени, на который умножаются скорости, чтобы сразу получить перемемещения вихрей (по умолчанию 1.0)
	/// \warning Использует OMP, MPI
	virtual void CalcConvVelo(double dt = 1.0) = 0;
	
	/// Деструктор
	virtual ~Velocity() { };
};

#endif
