/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: World2D.h                                                        |
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
\brief Заголовочный файл с описанием класса World2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#ifndef WORLD2D_H
#define WORLD2D_H

#include <iostream>
#include <memory>
#include <vector>
#include "Eigen/Dense"

#include "defs.h"
#include "Gpu.h"

#include "Airfoil.h"
#include "AirfoilRect.h"

#include "Boundary.h"
#include "BoundaryMDV.h"
#include "BoundaryVortColl.h"
#include "BoundaryConstLayerAver.h"

#include "Mechanics.h"
#include "MechanicsRigidImmovable.h"
#include "MechanicsRigidGivenLaw.h"
#include "MechanicsRigidOscillPart.h"
#include "MechanicsRigidOscillMon.h"
#include "MechanicsRigidRotateMon.h"

#include "Velocity.h"
#include "VelocityBiotSavart.h"
#include "VelocityFourier.h"

#include "MeasureVP.h"

#include "LogStream.h"
#include "Parallel.h"
#include "Passport.h"
#include "Point2D.h"
#include "Times.h"
#include "Wake.h"

/*!
\brief Класс, опеделяющий текущую решаемую задачу
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
class World2D
{
private:
	/// Поток для вывода логов и сообщений об ошибках
	mutable LogStream info;
	
	/// Список умных казателей на обтекаемые профили
	std::vector<std::unique_ptr<Airfoil>> airfoil;

	/// Список умных казателей на обтекаемые профили для сохранения старого положения
	std::vector<std::unique_ptr<Airfoil>> oldAirfoil;

	/// Список умных указателей на формирователи граничных условий на профилях
	std::vector<std::unique_ptr<Boundary>> boundary;

	/// Список умных указателей на типы механической системы для каждого профиля
	std::vector<std::unique_ptr<Mechanics>> mechanics;

	/// Умный укзатель на объект, определяющий методику вычисления скоростей
	std::unique_ptr<Velocity> velocity;

	/// Вихревой след
	Wake wake;

	/// Поле скоростей и давления (для сохранения в файл)
	std::unique_ptr<MeasureVP> measureVP;

	/// Матрица системы
	Eigen::MatrixXd matr;

	/// Обратная матрица
	Eigen::MatrixXd invMatr;

	/// Правая часть системы
	Eigen::VectorXd rhs;

	/// Решение системы
	Eigen::VectorXd sol;

	/// Константная ссылка на паспорт конкретного расчета
	const Passport& passport;

	/// Константная ссылка на параметры исполнения задачи в параллельном MPI-режиме
	const Parallel& parallel;

	/// Объект, управляющий графическим ускорителем
	Gpu cuda;

	/// Текущий номер шага в решаемой задаче
	size_t currentStep;

	/// Сведения о временах выполнения основных операций
	Times timestat;
	
public:	
	/// \todo откомментировать
	LogStream& getInfo() const { return info; };

	/// \todo откомментировать
	std::ostream& getInfo(char x) const { return info(x); };

	/// \todo откомментировать
	const Airfoil& getAirfoil(size_t i) const { return *airfoil[i]; };
	
	/// \todo откомментировать
	Airfoil& getNonConstAirfoil(size_t i) const { return *airfoil[i]; };
	
	/// \todo откомментировать
	size_t getNumberOfAirfoil() const { return airfoil.size(); };
	

	/// \todo откомментировать
	const Boundary& getBoundary(size_t i) const { return *boundary[i]; };
	
	/// \todo откомментировать
	Boundary& getNonConstBoundary(size_t i) const { return *boundary[i]; };
	
	/// \todo откомментировать
	size_t getNumberOfBoundary() const { return boundary.size(); };


	/// \todo откомментировать
	const Wake& getWake() const { return wake; };

	/// \todo откомментировать
	const Velocity& getVelocity() const { return *velocity; };

	/// \todo откомментировать
	Velocity& getNonConstVelocity() const { return *velocity; };

	/// \todo откомментировать
	const Passport& getPassport() const { return passport; };
	
	/// \todo откомментировать
	const Parallel& getParallel() const { return parallel; };

	/// \todo откомментировать
	const Gpu& getCuda() const { return cuda; };	

	/// \todo откомментировать
	size_t getCurrentStep() const { return currentStep; };

	/// \todo откомментировать
	const Times& getTimestat() const { return timestat; };

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
	/// \param[out] convTime ссылка на промежуток времени для вычисления конвективных скоростей --- пару чисел (время начала и время конца операции)
	/// \param[out] diffTime ссылка на промежуток времени для вычисления диффузионных скоростей --- пару чисел (время начала и время конца операции)
	void CalcVortexVelo(double dt, timePeriod& convTime, timePeriod& diffTime);

	/// \brief Вычисляем новые положения вихрей (в пелене и виртуальных)
	///
	/// Вызывается в Step()
	/// \param[in] dt шаг по времени
	/// \param[out] newPos новые позиции вихрей
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void MoveVortexes(double dt, std::vector<Point2D>& newPos, timePeriod& time);

	/// \brief Проверка проникновения вихрей внутрь  профиля
	///
	/// Вызывается в Step()
	/// \param[in] newPos новые позиции вихрей
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	void CheckInside(std::vector<Point2D>& newPos, timePeriod& time);
	
	/// \brief Конструктор
	///
	/// \param[in] passport_ константная ссылка на паспорт расчета
	/// \param[in] parallel_ коенстантная ссылка на параметры исполнения задачи в параллельном MPI-режиме
	World2D(const Passport& passport_, const Parallel& parallel_);
	
	/// Деструктор
	~World2D() { };

	/// \brief Функция, возвращающая признак завершения счета в решаемой задаче
	///
	/// true, если задача решена и выполнен признак останова; false если требуется еще выполнять шаги по времени
	bool isFinished() const;      

	/// Основная функция выполнения одного шага по времени	
	void Step();

	/// Метод-обертка для вызова метода генерации заголовка файла нагрузок и заголовка файла положения (последнее --- если профиль движется) 
	/// \param[in] mechanicsNumber номер профиля, для которого генерируется заголовок файла
	void GenerateMechanicsHeader(size_t mechanicsNumber)
	{
		mechanics[mechanicsNumber]->GenerateForcesHeader();
		mechanics[mechanicsNumber]->GeneratePositionHeader();
	}
};

#endif

