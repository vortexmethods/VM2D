/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.0    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2017/12/01     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina       |
*-----------------------------------------------------------------------------*
| File name: World2D.h                                                        |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.	                                      |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.	                                                          |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/


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

#include "BoundaryMDV.h"
#include "BoundaryVortColl.h"
#include "BoundaryConstLayerAver.h"


#include "MechanicsRigidImmovable.h"

#include "VelocityBiotSavart.h"
#include "VelocityFourier.h"

#include "AirfoilRect.h"

#include "Times.h"

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

	/// Список умных указателей на типы механической системы для каждого профиля
	std::vector<std::unique_ptr<Mechanics>> mechanics;

	/// Вихревой след
	Wake wake;

	/// Матрица системы
	Eigen::MatrixXd matr;

	/// Обратная матрица
	Eigen::MatrixXd invMatr;

	/// Правая часть системы
	Eigen::VectorXd rhs; 
	
	/// Решение системы
	Eigen::VectorXd sol;
	
	/// Умный укзатель на объект, определяющий методику вычисления скоростей
	std::unique_ptr<Velocity> velocity;

	/// \brief Решение системы линейных алгебраических уравнений
	///
	/// Вызывается в Step()
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void SolveLinearSystem(timePeriod& time);

	/// \brief Заполнение матрицы системы для всех профилей
	///
	/// Вызывается в Step()
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void FillMatrixAndRhs(timePeriod& time);

	/// \brief Вычисляем размер матрицы и резервируем память под нее и под правую часть
	///
	/// Вызывается в Step()
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void ReserveMemoryForMatrixAndRhs(timePeriod& time);
	
	/// \brief Вычисляем скорости вихрей (в пелене и виртуальных)
	///
	/// Вызывается в Step()
	/// \param[in] dt шаг по времени
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void CalcVortexVelo(double dt, timePeriod& time);

	/// \brief Вычисляем новые положения вихрей (в пелене и виртуальных)
	///
	/// Вызывается в Step()
	/// \param[in] dt шаг по времени
	/// \param[out] newPos новые позиции вихрей
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void MoveVortexes(double dt, std::vector<Point2D>& newPos, timePeriod& time);

	/// \brief Проверка проникновения вихрей внутрь профиля
	///
	/// Вызывается в Step()
	/// \param[in] newPos новые позиции вихрей
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void CheckInside(std::vector<Point2D>& newPos, timePeriod& time);

public:
	/// \brief Конструктор
	///
	/// \param[in] _queue константная ссылка на очередь задач
	/// \param[in,out] _telefile ссылка на поток для вывода временной статистики 
	World2D(const Queue& _queue, std::ostream& _telefile);
	
	/// Деструктор
	~World2D() { };

	/// \brief Функция, возвраящающая константную ссылку на паспорт конкретного расчета
	///
	/// \return константная ссылка на паспорт задачи
	const Passport& GetPassport() const;

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

	/// Сведения о временах выполнения основных операций
	Times timestat;

	/// Основная функция выполнения одного шага по времени
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void Step(timePeriod& time);

	/// Метод-обертка для вызова метода генерации заголовка файла нагрузок
	/// \param mechanicsNumber номер профиля, для которого генерируется заголовок файла
	void GenerateMechanicsHeader(int mechanicsNumber)
	{
		mechanics[mechanicsNumber]->GenerateForcesHeader(mechanicsNumber);
	}
};

#endif

