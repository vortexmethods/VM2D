/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Mechanics.h                                                      |
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
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/

#ifndef MECHANICS_H
#define MECHANICS_H

#include <memory>
#include <vector>
#include "Eigen/Dense"

#include "defs.h"
#include "Airfoil.h"
#include "Boundary.h"
#include "Point2D.h"
#include "StreamParser.h"
#include "Times.h"
#include "Velocity.h"

class World2D;

/*!
\brief Абстрактный класс, определяющий вид механической системы

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.4
\date 16 октября 2018 г.
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
	std::unique_ptr<StreamParser> mechParamsParser;

	/// Начальная скорость центра и угловая скорость
	const Point2D Vcm0;	const double Wcm0;

	/// Начальное положение профиля
	const Point2D Rcm0; double Phi0;
	
	/// Текущие скорость центра и угловая скорость
	Point2D Vcm; double Wcm;

	/// Текущие положение профиля
	Point2D Rcm; double Phi;

	/// Скорость и отклонение с предыдущего шага
	Point2D VcmOld;	double WcmOld;

	/// Текущие положение профиля
	Point2D RcmOld;	double PhiOld;

public:
	/// Переменная, отвечающая за то, двигается профиль или нет
	const bool isMoves;

	/// Переменная, отвечающая за то, деформируется профиль или нет
	const bool isDeform;

	/// Переменная, отвечающая за то, может профиль вращаться или нет
	const bool isRotate;

	/// Количество степеней свободы
	const size_t degOfFreedom;

	/// Вектор гидродинамической силы и момент, действующие на профиль
	Point2D hydroDynamForce;
	double hydroDynamMoment;


	/// \brief Чтение параметров конкретной механической системы
	///
	virtual void ReadSpecificParametersFromDictionary() = 0;

	/// \todo Проверить описание последних параметров конструктора Mechanics

	/// \brief Конструктор
	/// 
	/// \param[in] W_ константная ссылка на решаемую задачу
	/// \param[in] numberInPassport_ номер профиля в паспорте задачи
	/// \param[in] degOfFreedom_ количество степеней свободы
	/// \param[in] isMoves_ является ли профиль подвижным (1 - является, 0 - не является)
	/// \param[in] isDeform_ является ли профиль деформируемым (1 - является, 0 - не является)
	/// \param[in] isRotate_ является ли профиль вращающимся (1 - является, 0 - не является)
	/// \param[in] Vcm0_ - скорость центра масс
	/// \param[in] Rcm0_ - положение центра масс
	/// \param[in] Wcm0_ - угловая скорость центра масс 
	/// \param[in] Phi0_ - угол поворота центра масс
	Mechanics(const World2D& W_, size_t numberInPassport_, int degOfFreedom_, bool isMoves_, bool isDeform_, bool isRotate_, Point2D Vcm0_, Point2D Rcm0_, double Wcm0_, double Phi0_);
		
	
	/// Деструктор
	virtual ~Mechanics() { };

	/// \brief Вычисление гидродинамической силы, действующей на профиль
	///
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	virtual void GetHydroDynamForce(timePeriod& time) = 0;


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

	/// \brief Вычисление скоростей начал панелей
	/// \param[in] currTime текущее время
	virtual void VeloOfAirfoilPanels(double currTime) = 0;

	/// \brief Заполнение строк в матрице, соответствующих механической системе. А также пересечений строк и столбцов, соответствующих мех. системе
	///
	/// \param[in] row строки (количество строк = количество степеней свободы, количество элементов в строках = количество неизвестных в локальной матрице)
	/// \param[in] cross пересечения (количество столбцов = количество степеней свободы, количество элементов в столбцах = количество степеней свободы)
	virtual void FillMechanicsRowsAndCross(Eigen::MatrixXd& row, Eigen::MatrixXd& cross) = 0;

	/// \brief Заполнение правых частей, соответствующих механической системе
	///
	/// \param[in] rhs вектор правой части для механики одного профиля (количество элементов = количеству степеней свободы)
	virtual void FillMechanicsRhs(std::vector<double>& rhs) = 0;

	/// \brief Заполнение правых частей (расщепленная схема) либо коэффициентов в матрице при скорости (монолитный подход),
	///			соответствующих присоединенным слоям
	///
	/// \param[in] row строки (количество строк = количество степеней свободы, количество элементов в строках = количество неизвестных в локальной матрице)
	/// \param[in] rhs вектор правой части для механики одного профиля (количество элементов = количеству степеней свободы)
	virtual void FillAtt(Eigen::MatrixXd& row, Eigen::MatrixXd& rhs) = 0;

	/// \brief Извлечение из решения параметров механической системы
	///
	/// \param[in] col
	virtual void SolutionToMechanicalSystem(Eigen::VectorXd& col) = 0;

	/// \brief Перемещение профиля в соответствии с законом
	///
	virtual void Move() = 0;
};

#endif