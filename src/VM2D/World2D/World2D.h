/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.7    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2019/11/22     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2019 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
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
| VM is distributed in the hope that it will be useful, but WITHOUT           |
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
\version 1.7   
\date 22 ноября 2019 г.
*/

#ifndef WORLD2D_H
#define WORLD2D_H

#include "Gpu2D.h"
#include "Times2D.h"
#include "WorldGen.h"

namespace VMlib
{
	class Parallel;
}

namespace VM2D
{
	class Airfoil;
	class Boundary;
	class MeasureVP;
	class Mechanics;
	class Passport;
	class Times;
	class Velocity;	
	class Wake;
	class WakeDataBase;

	/*!
	\brief Класс, опеделяющий текущую решаемую задачу
	\author Марчевский Илья Константинович
	\author Кузьмина Ксения Сергеевна
	\author Рятина Евгения Павловна
	\version 1.7
	\date 22 ноября 2019 г.
	*/
	class World2D : public VMlib::WorldGen
	{
	private:
		/// Список умных казателей на обтекаемые профили
		std::vector<std::unique_ptr<Airfoil>> airfoil;

		/// Список умных казателей на обтекаемые профили для сохранения старого положения
		std::vector<std::unique_ptr<Airfoil>> oldAirfoil;

		/// Список умных указателей на формирователи граничных условий на профилях
		std::vector<std::unique_ptr<Boundary>> boundary;

		/// Список номеров, с которых начинаются элементы правой части (или матрицы) системы для профилей
		std::vector<size_t> dispBoundaryInSystem;

		/// Список умных указателей на типы механической системы для каждого профиля
		std::vector<std::unique_ptr<Mechanics>> mechanics;

		/// Умный укзатель на объект, определяющий методику вычисления скоростей
		std::unique_ptr<Velocity> velocity;

		/// Умный указатель на вихревой след
		std::unique_ptr<Wake> wake;

		/// Умный указатель на источники
		std::unique_ptr<WakeDataBase> source;

		/// Умный указатель на алгоритм вычисления полей скоростей и давления (для сохранения в файл)
		std::unique_ptr<MeasureVP> measureVP;

		/// Матрица системы
		Eigen::MatrixXd matr;

		/// Матрица, состоящая из пар матриц, в которых хранятся касательные и нормальные компоненты интегралов от ядра
		std::vector<std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>>> IQ;

		/// Обратная матрица
		Eigen::MatrixXd invMatr;

		/// Правая часть системы
		Eigen::VectorXd rhs;

		/// Решение системы
		Eigen::VectorXd sol;

		/// Константная ссылка на паспорт конкретного расчета
		const Passport& passport;

		/// Объект, управляющий графическим ускорителем
		mutable Gpu cuda;

	public:
		/// \brief Возврат константной ссылки на объект профиля
		///
		/// \param[in] i номер профиля, константная ссылка на который возвращается
		/// \return константную ссылку на i-й профиль
		const Airfoil& getAirfoil(size_t i) const { return *airfoil[i]; };

		/// \brief Возврат неконстантной ссылки на объект профиля
		///
		/// \param[in] i номер профиля, неконстантная ссылка на который возвращается
		/// \return неконстантную ссылку на i-й профиль
		Airfoil& getNonConstAirfoil(size_t i) const { return *airfoil[i]; };

		/// \brief Возврат количества профилей в задаче
		///
		/// \return количество профилей в задаче
		size_t getNumberOfAirfoil() const { return airfoil.size(); };

		/// \brief Возврат константной ссылки на объект граничного условия
		///
		/// \param[in] i номер граничного условия, константная ссылка на которое возвращается
		/// \return константную ссылку на i-е граничное условие
		const Boundary& getBoundary(size_t i) const { return *boundary[i]; };

		/// \brief Возврат неконстантной ссылки на объект граничного условия
		///
		/// \param[in] i номер граничного условия, неконстантная ссылка на которое возвращается
		/// \return неконстантную ссылку на i-е граничное условие
		Boundary& getNonConstBoundary(size_t i) const { return *boundary[i]; };

		/// \brief Возврат количества граничных условий в задаче
		///
		/// \return количество граничных условий в задаче
		size_t getNumberOfBoundary() const { return boundary.size(); };

		/// \brief Возврат смещения в системе dispBoundaryInSystem
		///
		/// \param[in] i номер граничного условия, константная ссылка на которое возвращается
		/// \return константную ссылку на i-е граничное условие
		size_t getDispBoundaryInSystem(size_t i) const { return dispBoundaryInSystem[i]; };

		/// \brief Возврат константной ссылки на measureVP
		///
		/// \return константную ссылку на вихревой след
		const MeasureVP& getMeasureVP() const { return *measureVP; };

		/// \brief Возврат неконстантной ссылки на measureVP
		///
		/// \return константную ссылку на вихревой след
		MeasureVP& getNonConstMeasureVP() const { return *measureVP; };

		/// \brief Возврат константной ссылки на объект механики
		///
		/// \param[in] i номер механики, константная ссылка на который возвращается
		/// \return константную ссылку на i-ю механику
		const Mechanics& getMechanics(size_t i) const { return *mechanics[i]; };

		/// \brief Возврат константной ссылки на вихревой след
		///
		/// \return константную ссылку на вихревой след
		const Wake& getWake() const { return *wake; };

		/// \brief Возврат константной ссылки на источники в области течения
		///
		/// \return константную ссылку на источники в области течения
		const WakeDataBase& getSource() const { return *source; };

		/// \brief Возврат константной ссылки на объект для вычисления скоростей
		///
		/// \return константную ссылку на объект для вычисления скоростей
		const Velocity& getVelocity() const { return *velocity; };

		/// \brief Возврат неконстантной ссылки на объект для вычисления скоростей
		///
		/// \return неконстантную ссылку на объект для вычисления скоростей
		Velocity& getNonConstVelocity() const { return *velocity; };

		/// \brief Возврат константной ссылки на паспорт
		///
		/// \return константную ссылку на паспорт
		const Passport& getPassport() const { return passport; };

		/// \brief Возврат константной ссылки на объект, связанный с видеокартой (GPU)
		///
		/// \return константную ссылку на объект, связанный с видеокартой (GPU)
		const Gpu& getCuda() const { return cuda; };

		/// \brief Возврат неконстантной ссылки на объект, связанный с видеокартой (GPU)
		///
		/// \return неконстантную ссылку на объект, связанный с видеокартой (GPU)
		Gpu& getNonConstCuda() const { return cuda; };

		/// \brief Возврат константной ссылки на объект, связанный с матрицей интегралов от (r-xi)/|r-xi|^2
		///
		/// \return константную ссылку на объект, связанный с матрицей интегралов от (r-xi)/|r-xi|^2
		const std::pair<Eigen::MatrixXd, Eigen::MatrixXd>& getIQ(size_t i, size_t j) const { return IQ[i][j]; };

		/// \brief Возврат ссылки на временную статистику выполнения шага расчета по времени
		///
		/// \return ссылку на временную статистику выполнения шага расчета по времени
		Times& getTimestat() { return dynamic_cast<Times&>(*timestat); };

		/// \brief Решение системы линейных алгебраических уравнений
		///
		/// Вызывается в Step()
		/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
		void SolveLinearSystem(timePeriod& time);

		/// \brief Заполнение матрицы, состоящей из интегралов от (r-xi) / |r-xi|^2
		///
		/// Вызывается в Step()
		void FillIQ();

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
		void CalcVortexVelo(double dt, timePeriod& convTime, timePeriod& diffTime, timePeriod& vpTime);

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
		World2D(const VMlib::PassportGen& passport_, const VMlib::Parallel& parallel_);

		/// Деструктор
		~World2D() { };
		
		/// Метод-обертка для вызова метода генерации заголовка файла нагрузок и заголовка файла положения (последнее --- если профиль движется) 
		/// \param[in] mechanicsNumber номер профиля, для которого генерируется заголовок файла
		void GenerateMechanicsHeader(size_t mechanicsNumber);



		// Реализация виртуальных функций
		virtual void Step() override;
	};

}//namespace VM2D

#endif

