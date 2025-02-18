/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Mechanics2D.h                                                    |
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
\brief Заголовочный файл с описанием класса Mechanics
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/

#ifndef MECHANICS_H
#define MECHANICS_H

#include <memory>

#include "defs.h"

namespace VMlib
{
	class StreamParser;
}

namespace VM2D
{
	class Airfoil;
	class Boundary;
	class World2D;
	struct VortexesParams;

	/*!
	\brief Абстрактный класс, определяющий вид механической системы

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
\author Колганова Александра Олеговна

	\Version 1.12
	\date 14 января 2024 г.
	*/

	class Mechanics
	{
	private:
		/// Парсинг списка параметров механической системы
		void ReadParametersFromDictionary();

	protected:
		/// Константная ссылка на решаемую задачу
		const World2D& W;

		/// Номер профиля в паспорте
		const size_t numberInPassport;

		/// Константная ссылка на профиль
		/// \n инициализируется автоматом в конструкторе, при помощи const_cast
		/// \warning использует const_cast для получения неконстантной ссылки
		Airfoil& afl;

		/// Константная ссылка на граничное условие
		/// \n инициализируется автоматом в конструкторе
		const Boundary& boundary;

		/// Константная ссылка на структуру с параметрами виртуального вихревого слоя для профиля
		/// \n инициализируется автоматом в конструкторе
		const VortexesParams& virtVortParams;

		/// Умный указатель на парсер параметров механической системы
		std::unique_ptr<VMlib::StreamParser> mechParamsParser;

	public:

		/// Начальная скорость центра и угловая скорость
		Point2D Vcm0; double Wcm0;

		/// Начальное положение профиля
		Point2D Rcm0; double Phi0;

		/// Текущие скорость центра и угловая скорость
		Point2D Vcm; double Wcm;	

		/// Текущие положение профиля
		Point2D Rcm; double Phi;

		/// Скорость и отклонение с предыдущего шага
		Point2D VcmOld;	double WcmOld;

		/// Текущие положение профиля
		Point2D RcmOld;	double PhiOld;

		/// Текущая циркуляция скорости по границе профиля
		double circulation;

		/// Циркуляция скорости по границе профиля с предыдущего шага
		double circulationOld;

	public:
		/// Переменная, отвечающая за то, двигается профиль или нет
		const bool isMoves;

		/// Переменная, отвечающая за то, деформируется профиль или нет
		const bool isDeform;

		/// Вектор гидродинамической силы и момент, действующие на профиль
		Point2D hydroDynamForce;
		double hydroDynamMoment;

		/// Вектор силы и момент вязкого трения, действующие на профиль
		Point2D viscousForce;
		double viscousMoment;


		/// \brief Чтение параметров конкретной механической системы
		///
		virtual void ReadSpecificParametersFromDictionary() = 0;

		/// \brief Конструктор
		/// 
		/// \param[in] W_ константная ссылка на решаемую задачу
		/// \param[in] numberInPassport_ номер профиля в паспорте задачи		
		/// \param[in] isMoves_ является ли профиль подвижным (1 - является, 0 - не является)
		/// \param[in] isDeform_ является ли профиль деформируемым (1 - является, 0 - не является)
		Mechanics(const World2D& W_, size_t numberInPassport_, bool isMoves_, bool isDeform_);

		/// Деструктор
		virtual ~Mechanics() { };


		/// \brief Задание начального положения и начальной скорости
		/// 
		/// \param[in] Vcm0_ - скорость центра масс
		/// \param[in] Rcm0_ - положение центра масс
		/// \param[in] Wcm0_ - угловая скорость центра масс 
		/// \param[in] Phi0_ - угол поворота центра масс
		void Initialize(Point2D Vcm0_, Point2D Rcm0_, double Wcm0_, double Phi0_);


		/// Вычисление гидродинамической силы, действующей на профиль
		virtual void GetHydroDynamForce() = 0;

		/// Генерация заголовка файла нагрузок
		void GenerateForcesHeader();

		/// Генерация заголовка файла положения профиля
		void GeneratePositionHeader();

		/// Сохранение строки со статистикой в файл нагрузок
		void GenerateForcesString();

		/// Сохранение строки со статистикой в файл нагрузок
		void GeneratePositionString();


		/// \brief Вычисление скорости центра масс профиля
		///
		/// \param[in] currTime текущее время
		virtual Point2D VeloOfAirfoilRcm(double currTime) = 0;

		/// \brief Вычисление положения центра масс профиля
		///
		/// \param[in] currTime текущее время
		virtual Point2D PositionOfAirfoilRcm(double currTime) = 0;

		/// \brief Вычисление угловой скорости профиля
		///
		/// \param[in] currTime текущее время
		virtual double AngularVelocityOfAirfoil(double currTime) = 0;

		/// \brief Вычисление угла поворота профиля
		///
		/// \param[in] currTime текущее время
		virtual double AngleOfAirfoil(double currTime) = 0;

		/// \brief Вычисление скоростей начал панелей
		/// \param[in] currTime текущее время
		virtual void VeloOfAirfoilPanels(double currTime) = 0;

		/// \brief Перемещение профиля в соответствии с законом
		///
		virtual void Move() = 0;
	};

}//namespace VM2D

#endif