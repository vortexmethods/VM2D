/*!
\file
\brief Заголовочный файл с описанием класса World2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef WORLD2D_H
#define WORLD2D_H

#include "BoundaryVortColl.h"
#include "BoundaryConstLayerAver.h"

#include "VelocityDirect.h"
#include "VelocityFourier.h"

#include "AirfoilRect.h"

class Queue; //дальнее описание

/*!
\brief Класс, опеделяющий текущую решаемую задачу
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class World2D
{
private:
	/// Константная ссылка на очередь задачи
	const Queue& queue;

	/// Константная ссылка на параметры исполнения задачи в параллельном MPI-режиме
	const Parallel& parallel;

	// Ссылка на поток для вывода временной статистики
	std::ostream& teleFile;

	/// Список умных казателей на обтекаемые профили
	std::vector<std::unique_ptr<Airfoil>> airfoil;

	/// Список умных указателей на формирователи граничных условий на профилях
	std::vector<std::unique_ptr<Boundary>> boundary;

	/// Вихревой след
	Wake wake;

	/// Матрица системы
	Eigen::MatrixXd matr;

	/// Правая часть системы
	Eigen::VectorXd rhs; 
	
	/// Решение системы
	Eigen::VectorXd sol;
	
	/// Умный укзатель на объект, определяющий методику вычисления скоростей
	std::unique_ptr<Velocity> velocity;

	/// \brief Функция, возвраящающая константную ссылку на паспорт конкретного расчета
	///
	/// \return константная ссылка на паспорт задачи
	const Passport& Passport() const;

public:
	/// \brief Конструктор
	///
	/// \param[in] _queue константная ссылка на очередь задач
	/// \param[in,out] _telefile ссылка на поток для вывода временной статистики 
	World2D(const Queue& _queue, std::ostream& _telefile);
	
	/// Деструктор
	~World2D() { };

	/// Текущее время в решаемой задаче
	double currentTime; 

	/// Текущий номер шага в решаемой задаче
	int currentStep; 

	/// \brief Функция, возвращающая номер текущей решаемой задачи (для справки)
	///
	/// \return номер решаемой задачи
	int problemNumber() const;    

	/// \brief Функция, возвращающая признак завершения счета в решаемой задаче
	///
	/// true, если задача решена и выполнен признак останова; false если требуется еще выполнять шаги по времени
	bool isFinished() const;      

	/// Основная функция выполнения одного шага по времени
	void Step(); 
};

#endif

