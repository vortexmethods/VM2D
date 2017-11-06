/*!
\file
\brief Описание базовых вспомогательных функций
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/


#ifndef DEFS_H
#define DEFS_H

#include <iostream>
#include <fstream>
#include <string>

#include "Ext_Lib/Eigen/Dense"

#include "Point2D.h"


/// Число \f$ \pi \f$
const double   PI = 3.1415926535897932384626433832795;

/// Число \f$ \frac{1}{2\pi} \f$
const double IDPI = 0.1591549430918953357688837633725;

/// Число \f$ 2\pi \f$
const double  DPI = 6.2831853071795864769252867665590;

/// \brief Глобальные параметры по умолчанию
///
/// Они могут не указываться ни в параметрах задач, ни в файле defaults.txt
namespace defaults
{
	/// Начало расчета
	const double defaultTimeStart = 0.0;
	
	/// Шаг сохранения в текстовый файл
	const int defaultDeltacntText = 1;

	/// Шаг сохранения в бинарный файл
	const int defaultDeltacntBinary = 0;

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

	/// Имя файла с паспортом задачи
	const std::string defaultPspFile = "passport.txt";
	
	/// Необходимое число процессоров для решения задачи
	const int defaultNp = 1;

	/// Поток вывода логов
	static std::ostream* defaultPinfo = &std::cout;
	
	/// Поток вывода ошибок
	static std::ostream* defaultPerr = &std::cout;

	/// ПОток вывода телеметрии
	static std::ostream* defaultPtele = &std::cout;
	
} //namespace defaults


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
	int n = _vec.size();
	_stream << "{ ";
	if (n > 0)
	{
		
		for (int j = 0; j < n - 1; ++j)
			_stream << _vec[j] << ", ";
		_stream << _vec[n - 1];
	}
	_stream << " }";
	return _stream;
}//operator<<(...)

/// \brief Проверка существования файла
///
/// \param[in] fileName константная ссылка на имя проверяемого файла
/// \param[in] _info ссылка на поток вывода логов
/// \param[in] _err ссылка на поток вывода ошибок
/// \param[in] _str константная ссылка на текстовую метку, которая выводится в потоки
/// \return true, если файл существует; false если файл отсутствует
inline bool fileExistTest(const std::string& fileName, std::ostream& _info, std::ostream& _err, const std::string& _str)
{
	std::ifstream f(fileName.c_str());
	if (f.good())
	{
		f.close();
		f.clear();
		_info << _str << " info: FILE " << fileName << " IS FOUND" << std::endl;
		return true;
	}
	else
	{
		_err << _str << " error: FILE " << fileName << " IS NOT FOUND" << std::endl;
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