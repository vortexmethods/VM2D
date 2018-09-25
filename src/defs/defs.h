/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.1    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/04/02     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: defs.h                                                           |
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
\brief Описание базовых вспомогательных функций
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.1
\date 2 апреля 2018 г.
*/


#ifndef DEFS_H
#define DEFS_H

#include <iostream>
#include <fstream>
#include <string>
#include <ctime>

#include "Eigen/Dense"

#include "gpudefs.h"
#include "Point2D.h"
   
/// Тип для хранения начала и конца промежутка времени
typedef std::pair<double, double> timePeriod;

/// Число \f$ \pi \f$
const double   PI = 3.1415926535897932384626433832795;

/// Число \f$ \frac{1}{2\pi} \f$
const double IDPI = 0.1591549430918953357688837633725;

/// Число \f$ 2\pi \f$
const double  DPI = 6.2831853071795864769252867665590;

/// \brief Глобальные параметры по умолчанию
///
/// Они могут не указываться ни в параметрах задач, ни в файле defaults
namespace defaults
{
	/// Начало расчета
	const double defaultTimeStart = 0.0;

	/// Время разгона
	const double defaultTimeAccel = 0.0;
	
	/// Шаг сохранения в текстовый файл
	const int defaultDeltacntText = 1;

	/// Шаг сохранения в бинарный файл
	const int defaultDeltacntBinary = 0;

	/// Радиус убивания дальнего следа
	const double defaultDistKill = 10.0;

	/// Расстояние, на которое рождаемый вихрь отодвигается от профиля
	const double defaultDelta = 1.e-5;

	/// Каталог с файлами профилей
	const std::string defaultAirfoilsDir = "./airfoils/";

	/// Каталог с файлами вихревых следов
	const std::string defaultWakesDir = "./wakes/";

	/// Список профилей
	const std::vector<std::string> defaultFileAirfoil({});

	/// Файл со следом
	const std::string defaultFileWake("");

	/// Список профилей
	const std::vector<std::string> defaultAirfoil({});

	/// Базовое смещение профиля
	const Point2D defaultBasePoint = { 0.0, 0.0 };
	
	/// Коэффициент масштабирования профиля
	const double defaultScale = 1.0;

	/// Угол атаки
	const double defaultAngle = 0.0;
	
	/// Тип панелей
	const int defaultPanelsType = 0;
	
	/// Способ удовлетворения граничного условия
	const int defaultBoundaryCondition = 0;

	/// Тип механической системы
	const int defaultMechanicalSystem = 0;

	/// Имя файла с паспортом задачи
	const std::string defaultPspFile = "passport";
	
	/// Необходимое число процессоров для решения задачи
	const int defaultNp = 1;

	/// Имя файла с паспортом задачи для копирования в новые каталоги
	const std::string defaultCopyPspFile = "";

	/// Поток вывода логов
	static std::ostream* defaultPinfo = &std::cout;
	
	/// Поток вывода ошибок
	static std::ostream* defaultPerr = &std::cout;

	/// Поток вывода телеметрии
	static std::ostream* defaultPtele = &std::cout;
	
} //namespace defaults


/// \brief Формирование заголовка файла программы VM2D
///
/// Печатает в шапку файла заголовок программы VM2D
///
/// \param[out] str ссылка на файл для вывода, должен быть открыт
/// \param[in] fileName константная ссылка на строку с именем файла
/// \param[in] descr константная ссылка на строку с описанием файла
void PrintLogoToTextFile(std::ofstream& str, const std::string& fileName, const std::string& descr);


/// \brief Формирование подзаголовка в текстовом файле вывода программы VM2D
///
/// Печатает в файл вывода программы VM2D подзаголовок и подчеркивает его
///
/// \param[out] str ссылка на файл для вывода, должен быть открыт
/// \param[in] header константная ссылка на строку с заголовком
void PrintHeaderToTextFile(std::ofstream& str, const std::string& header);


/// \brief Переопределение оператора "<<" для вывода в поток вектора std::vector
///
/// Компоненты вектора выводятся в фигурных скобках с разделителем ";"
///
/// \tparam T тип элементов в векторе
/// \param[in,out] _stream ссылка на поток для вывода
/// \param[in] _vec константная ссылка на выводимый вектор
/// \return ссылка на поток
template <typename T>
std::ostream& operator<< (std::ostream& _stream, const std::vector<T>& _vec)
{
	size_t n = _vec.size();
	_stream << "{ ";
	if (n > 0)
	{		
		for (size_t j = 0; j < n - 1; ++j)
			_stream << _vec[j] << ", ";
		_stream << _vec[n - 1];
	}
	_stream << " }";
	return _stream;
}//operator<<(...)

/// \brief Проверка существования файла
///
/// \param[in] fileName константная ссылка на имя проверяемого файла
/// \param[in, out] Pinfo указатель на поток вывода логов
/// \param[in, out] Perr указатель на поток вывода ошибок
/// \param[in] _str константная ссылка на текстовую метку, которая выводится в потоки
/// \return true, если файл существует; false если файл отсутствует
inline bool fileExistTest(const std::string& fileName, std::ostream *Pinfo, std::ostream *Perr, const std::string& _str)
{
	std::ifstream f(fileName.c_str());
	if (f.good())
	{
		f.close();
		f.clear();

		if (Pinfo != nullptr)
			*Pinfo << _str << " info: FILE " << fileName << " IS FOUND" << std::endl;
		return true;
	}
	else
	{
		if (Perr != nullptr)
			*Perr << _str << " error: FILE " << fileName << " IS NOT FOUND" << std::endl;
		exit(-1);
		return false;
	}
}


/// \brief Возведение числа в квадрат
/// 
/// \tparam T тип данных
/// \param[in] x аргумент
/// \return квадрат аргумента \f$ x^2 \f$
template<typename T>
inline T sqr(T x) { return x*x; }

/// \brief Сохранение матрицы в поток
/// 
/// \param[in] matr константная ссылка на сохраняемую матрицу
/// \param[in,out] str ссылка на поток для сохранения
void SaveToStream(const Eigen::MatrixXd& matr, std::ostream& str);

/// \brief Сохранение комплекснозначной матрицы в поток
/// 
/// \param[in] matr константная ссылка на сохраняемую комплекснозначную матрицу
/// \param[in,out] str ссылка на поток для сохранения
void SaveToStream(const Eigen::MatrixXcd& matr, std::ostream& str);

/// \brief Сохранение вектора в поток
/// 
/// \param[in] vec константная ссылка на сохраняемый вектор
/// \param[in,out] str ссылка на поток для сохранения
void SaveToStream(const Eigen::VectorXd& vec, std::ostream& str);

/// \brief Сохранение списка из двумерных векторов (точек) в поток
/// 
/// \param[in] vec константная ссылка на сохраняемый список из двумерных векторов (точек)
/// \param[in,out] str ссылка на поток для сохранения
void SaveToStream(const std::vector<Point2D>& vec, std::ostream& str);

/// \brief Ядро сглаживания (Монагана)
///
/// \param[in] t аргумент ядра
/// \return Значение ядра сглаживания в конкретной точке
double M4(double t);


#endif