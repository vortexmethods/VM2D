/*--------------------------------*- VMlib -*----------------*---------------*\
| ##  ## ##   ## ##   ##  ##    |                            | Version 1.12   |
| ##  ## ### ### ##       ##    |  VMlib: VM2D/VM3D Library  | 2024/01/14     |
| ##  ## ## # ## ##   ##  ####  |  Open Source Code          *----------------*
|  ####  ##   ## ##   ##  ## ## |  https://www.github.com/vortexmethods/VM2D  |
|   ##   ##   ## #### ### ####  |  https://www.github.com/vortexmethods/VM3D  |
|                                                                             |
| Copyright (C) 2017-2024 Ilia Marchevsky                                     |
*-----------------------------------------------------------------------------*
| File name: defs.h                                                           |
| Info: Source code of VMlib                                                  |
|                                                                             |
| This file is part of VMlib.                                                 |
| VMLib is free software: you can redistribute it and/or modify it            |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VMlib is distributed in the hope that it will be useful, but WITHOUT        |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VMlib.  If not, see <http://www.gnu.org/licenses/>.              |
\*---------------------------------------------------------------------------*/


/*!
\file
\brief Описание базовых вспомогательных функций
\author Марчевский Илья Константинович
\Version 1.12
\date 14 января 2024 г.
*/


#ifndef DEFS_H
#define DEFS_H

#if defined(_WIN32)
	#include <direct.h>
#endif

#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>

#include <ctime>
#include <iostream>
#include <list>

#include "Eigen/Dense"

#include "Gpudefs.h"
#include "LogStream.h"

#include "Vortex2D.h"
#include "v3D.h"
#include "PairInt.h"

const int multipoleOrder = 12;
const double multipoleTheta = 1.2;


/// параметр расписания для распараллеливания циклов по OpenMP
#define DYN_SCHEDULE 20
   
/// Класс перечисления для определения типа набора точек (пелена/виртуальные вихри/точки для вычисления скорости и давления)
enum class PointType { wake, wakeVP, sheetGam, sourceWake, source };

/// Тип для хранения начала и конца промежутка времени
typedef std::pair<double, double> timePeriod;

/// Число \f$ \pi \f$
const double PI = 3.1415926535897932384626433832795;

/// Число \f$ \frac{1}{2\pi} \f$
const double IDPI = 0.1591549430918953357688837633725;

/// Число \f$ 2\pi \f$
const double DPI = 6.2831853071795864769252867665590;

/// Число \f$ \frac{1}{4\pi} \f$
const double IQPI = 0.07957747154594766788444188168626;

/// Число \f$ 4\pi \f$
const double QPI = 12.566370614359172953850573533118;

/// \brief Глобальные параметры по умолчанию
///
/// Они могут не указываться ни в параметрах задач, ни в файле defaults
namespace defaults
{
	/// Начало расчета
	const double defaultTimeStart = 0.0;

	/// Время разгона
	const std::pair<std::pair<std::string, int>, std::string> defaultVelAccel = { {"RampLin", 1}, "" };
	const double defaultTimeAccel = 0.0;
	
	/// Шаг подсчета поля скорости и давления
	const std::pair<std::pair<std::string, int>, std::string> defaultSaveVP = { {"text", 0}, "" };
	const int defaultSaveVPstep = 0;

	/// Шаг подсчета поля скорости и давления
	const std::pair<std::pair<std::string, int>, std::string> defaultSaveVtx = { {"text", 0}, "" };
	const int defaultSaveVtxStep = 100;

	/// Шаг подсчета поля скорости и давления
	const int defaultSaveVisStress = 0;

	/// Число разрядов в имени файла
	const int defaultNameLength = 5;

	/// Радиус убивания дальнего следа
	const double defaultDistFar = 10.0;

	/// Расстояние, на которое рождаемый вихрь отодвигается от профиля
	const double defaultDelta = 1.e-5;

	/// Число вихрей, рождаемых на одной панели
	const int defaultVortexPerPanel = 1;

	/// Число вихрей, рождаемых на одной панели
	const double defaultMaxGamma = 0.0;
		
	/// Референсная скорость, равная нулю, что означает ее равенство скорости набегающего потока
	const double defaultVRef = 0.0;

	/// Желаемое число панелей для разбиения геометрии
	const size_t defaultRequiredNPanels = 0;

	/// Способ удовлетворения граничного условия
	const std::pair<std::string, int> defaultBoundaryCondition = { "boundaryConstantLayerAver", 1 };

	/// Способ решения линейной системы
	const std::pair<std::string, int> defaultLinearSystemSolver = { "linearSystemGauss", 0 };


	/// Способ вычисления скоростей вихрей
	const std::pair<std::string, int> defaultVelocityComputation{ "velocityBiotSavart", 0 };

	/// Признак работы в "географической" системе координат
	const bool defaultGeographicalAngles = false;

	/// Признак поворота сил в профильную систему координат
	const bool defaultRotateForces = false;

	/// Признак расчета безразмерных коэффициентов вместо сил
	const bool defaultCalcCoefficients = false;

	/// Угол поворота точек VP
	const double rotateAngleVpPoints = 0;

	/// Каталог с файлами профилей
	const std::string defaultAirfoilsDir = "../settings/airfoils/";
    const std::string defaultBodiesDir = "../settings/bodies/";

	/// Каталог с файлами вихревых следов
	const std::string defaultWakesDir = "../settings/wakes/";

	/// Список профилей
	const std::vector<std::string> defaultFileAirfoil({});
    const std::vector<std::string> defaultFileBody({});

	/// Файл со следом
	const std::string defaultFileWake("");

	/// Файл с источниками
	const std::string defaultFileSource("");

	/// Список профилей
	const std::vector<std::string> defaultAirfoil({});
    const std::vector<std::string> defaultBody({});

	/// Базовое смещение профиля
	const VMlib::Point2D defaultBasePoint = { 0.0, 0.0 };
    const VMlib::v3D defaultBasePoint3D = { 0.0, 0.0, 0.0 };
	
	/// Коэффициент масштабирования профиля
	const Point2D defaultScale = { 1.0, 1.0 };

	/// Угол атаки
	const double defaultAngle = 0.0;

	/// Хорда
	const double defaultChord = 1.0;

	/// Присоединенная масса
	const Point2D defaultAddedMass = { 0.0, 0.0 };
	
	/// Признак разворота нормалей (для расчета внутреннего течения)
	const bool defaultInverse = false;
	
	/// Тип механической системы
	const int defaultMechanicalSystemType = 0;
	const std::string defaultMechanicalSystem = "";

	/// Имя файла с паспортом задачи
	const std::string defaultPspFile = "passport";

	/// Необходимое число процессоров для решения задачи
	const int defaultNp = 1;

	/// Путь к каталогу с задачей для копирования в новые каталоги
	const std::string defaultCopyPath = "";

	/// Поток вывода логов и ошибок очереди
	static std::ostream* defaultQueueLogStream = &std::cout;

	/// Поток вывода логов и ошибок задачи
	static std::ostream* defaultWorld2DLogStream = &std::cout;
	
	/// Расчет присоединенной массы
	static bool defaultAddMass = false;
	static v3D defaultAddMassVcm = { 0.0, 0.0, 0.0 };
	static v3D defaultAddMassWcm = { 0.0, 0.0, 0.0 };

	/// Для профиля на упругих связях - начальные отклонения и скорости	
	static Point2D defaultInitDisplacement = { 0.0, 0.0 };
	static double defaultInitAngularDisplacement = 0.0;
	static Point2D defaultInitVelocity = {0.0, 0.0};
	static double defaultInitAngularVelocity = 0.0;

} //namespace defaults


namespace VMlib
{

	/// \brief Формирование строки с текущем временем и датой
	///
	/// \return Возвращает строку с текущей датой и временем в формате XX Month 20XX at XX:XX:XX
	std::string CurrentDataTime();


	/// \brief Передача в поток вывода шапки программы VM2D/VM3D
	/// 
	/// Печатает в поток вывода заголовок программы VM2D/VM3D
	///
	/// \param[out] str ссылка на поток вывода
	void PrintLogoToStream(std::ostream& str);

	/// \brief Передача в поток вывода универсальной шапки программы VM2D/VM3D
	/// 
	/// Печатает в поток вывода заголовок программы VM2D/VM3D
	///
	/// \param[out] str ссылка на поток вывода
	void PrintUniversalLogoToStream(std::ostream& str);


	/// \brief Формирование заголовка файла программы VM2D/VM3D
	///
	/// Печатает в шапку файла заголовок программы VM2D/VM3D
	///
	/// \param[out] str ссылка на файл для вывода, должен быть открыт
	/// \param[in] fileName константная ссылка на строку с именем файла
	/// \param[in] descr константная ссылка на строку с описанием файла
	void PrintLogoToTextFile(std::ofstream& str, const std::string& fileName, const std::string& descr);


	/// \brief Формирование подзаголовка в текстовом файле вывода программы VM2D/VM3D
	///
	/// Печатает в файл вывода программы VM2D/VM3D подзаголовок и подчеркивает его
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


	/// \brief Переопределение оператора "<<" для вывода в поток пары ((строка, целое число), строка)
	///
	/// Компоненты пары выводятся в скобках с разделителем ","
	///
	/// \param[in,out] _stream ссылка на поток для вывода
	/// \param[in] _pair константная ссылка на выводимую пару
	/// \return ссылка на поток
	inline std::ostream& operator<< (std::ostream& _stream, const std::pair<std::pair<std::string, int>, std::string>& _pair)
	{
		_stream << "( ( " << _pair.first.first << ", " << _pair.first.second << " ), " << _pair.second << " )";	
		return _stream;
	}//operator<<(...)


	/// \brief Переопределение оператора "<<" для вывода в поток пары (строка, целое число)
	///
	/// Компоненты пары выводятся в скобках с разделителем ","
	///
	/// \param[in,out] _stream ссылка на поток для вывода
	/// \param[in] _pair константная ссылка на выводимую пару
	/// \return ссылка на поток
	inline std::ostream& operator<< (std::ostream& _stream, const std::pair<std::string, int>& _pair)
	{
		_stream << "( " << _pair.first << ", " << _pair.second << " )";
		return _stream;
	}//operator<<(...)



	/// \brief Проверка существования файла
	///
	/// \param[in, out] fileName ссылка на имя проверяемого файла, заменяется на имя + то из расширений, с которым файл присутствует
	/// \param[in, out] info ссылка на поток вывода логов/ошибок
	/// \param[in] extList константная ссылка на список проверяемых расширений
	/// \return true, если файл существует; false если файл отсутствует
	inline bool fileExistTest(std::string& fileName, LogStream& info, const std::list <std::string>& extList = {})
	{
		std::list<std::string> fullExtList(extList);
		fullExtList.insert(fullExtList.begin(), "");

		for (const auto& ext : fullExtList)
		{
			std::string newFileName = ((ext == "") ? fileName : (fileName + "." + ext));
			std::ifstream f(newFileName.c_str());
			if (f.good())
			{
				f.close();
				f.clear();

				info('i') << "file " << newFileName << " is found" << std::endl;
				fileName = newFileName;
				return true;
			}
		}

		info('e') << "file " << fileName << " is not found" << std::endl;
		exit(-1);
		return false;
	}


	/// \brief Формирование имени файла
	///
	/// \param[in] name константная ссылка на префикс имени файла
	/// \param[in] length количество знаков под номер
	/// \param[in] number номер для имени файла
	/// \param[in] ext константная ссылка на раширение (дописывается, если непустое)
	/// \return строку --- имя текстового файла
	inline std::string fileNameStep(const std::string& name, int length, size_t number, const std::string& ext)
	{
		std::string fname(name);

		size_t dec = 1;

		for (int i = 1; i < length; ++i)
		{
			dec *= 10;
			if (number < dec)
				fname += "0";
		}

		std::ostringstream ss;
		ss << number;
		fname += ss.str();
		
		if (ext.size() > 0)
		{
			fname += ".";
			fname += ext;
		}

		return fname;
	}


	/// \brief Копирование файла
	///
	/// \param[in] fileNameFrom константная ссылка на имя исходного файла
	/// \param[in] fileNameTo константная ссылка на имя нового файла	
	inline void copyFile(const std::string& fileNameFrom, const std::string& fileNameTo) 
	{
		std::string buf;
		buf.resize(BUFSIZ);

		FILE *in, *out;
		size_t n;

		
#pragma warning (push)
#pragma warning (disable: 4996)
		in = fopen(fileNameFrom.c_str(), "rb");
		out = fopen(fileNameTo.c_str(), "wb");
#pragma warning (pop)

		while ((n = fread((void*)buf.data(), 1, BUFSIZ, in)) != 0) 
		{
			fwrite((void*)buf.data(), 1, n, out);
		}
	}//copyFile


	/// \brief Создание каталога
	///
	/// \param[in] dir константная ссылка на имя текущего каталога (со слешем на конце)
	/// \param[in] name константная ссылка на имя создаваемого подкаталога
	inline void CreateDirectory(const std::string& dir, const std::string& name)
	{
#if defined(_WIN32)
		_mkdir((dir + name).c_str());
#else
		mkdir((dir + name).c_str(), S_IRWXU | S_IRGRP | S_IROTH);
#endif
	}
	

	/// \brief Возведение числа в квадрат
	/// 
	/// \tparam T тип данных
	/// \param[in] x аргумент
	/// \return квадрат аргумента \f$ x^2 \f$
	template<typename T>
	inline T sqr(T x) { return x * x; }


	/// \brief Усовершенствованный аркосинус
	/// 
	/// \tparam T тип данных
	/// \param[in] x аргумент
	/// \return аркосинус аргумента, если |x|<1, или 0, или Pi,
	template <typename T>
	inline double macos(const T x)
	{
		double res = abs(x) > 1.0 ? 0.5*PI*(1 - sign(x)) : acos(x);
		return res;
	}


	/// \brief Усовершенствованный сигнум
	/// 
	/// \tparam T тип данных
	/// \param[in] x аргумент
	/// \return знак числа в обычном математическом смысле,
	template <typename T>
	int sign(T x)
	{
		double dx = (double)x;
		if (dx > 0) return 1;
		if ( (dx < 0.0) || (dx < -0.0) ) return -1;
		return 0;
	}



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

	/// \brief Модифицирует массив квадратов расстояний до ближайших вихрей из wake
	///
	/// \param[out] ee2 массив расстояний до трех ближайших вихрей
	/// \param[in] dst2 новое расстояние (заменит одно из расстояний в ee2, если меньше)
	void ModifyE2(double* ee2, double dst2);


	/// \brief Вспомогательная функция вычисления угла между векторами
	///
	/// \param[in] p константная ссылка на первый вектор
	/// \param[in] s константная ссылка на второй вектор
	/// \return угол между векторами в диапазоне \f$ (-\pi; \pi] \f$
	double Alpha(const Point2D& p, const Point2D& s);

	/// \brief Вспомогательная функция вычисления логарифма отношения норм векторов
	///
	/// \param[in] p константная ссылка на первый вектор
	/// \param[in] s константная ссылка на второй вектор
	/// \return логарифм отношения норм векторов
	double Lambda(const Point2D& p, const Point2D& s);

	/// \brief Вспомогательная функция вычисления величины \f$ (\vec a \cdot \vec b) \cdot \vec c + (\vec a \times \vec b) \times \vec c \f$
	///
	/// Для оптимизации все векторы считаются двумерными
	///
	/// \param[in] a константная ссылка на первый вектор
	/// \param[in] b константная ссылка на второй вектор
	/// \param[in] c константная ссылка на третий вектор
	/// \return логарифм отношения норм векторов
	Point2D Omega(const Point2D& a, const Point2D& b, const Point2D& c);

	/// \brief Вспомогательная функция перестановки байт местами (нужно для сохранения бинарных VTK)
	///
	///
	/// \param[in] var ссылка на переменную, у которой байты нужно поменять местами
	template<typename T>
	void SwapEnd(T& var)
	{
		char* varArray = reinterpret_cast<char*>(&var);
		for (long i = 0; i < static_cast<long>(sizeof(var) / 2); ++i)
			std::swap(varArray[sizeof(var) - 1 - i], varArray[i]);
	}

	/// \brief Способ сглаживания скорости вихря (вихрь Рэнкина или вихрь Ламба)
	inline double boundDenom(double r2, double eps2)
	{
#ifndef LAMBVORTEX
		return std::max(r2, eps2);
#else
		if (r2 > eps2)
			return std::max(r2, eps2);
		else
			return (r2 < 1e-10) ? 1e-10 : r2 / (1.0 - exp(-6.0*r2 / eps2));
#endif
	}

}//namespace VMlib

using VMlib::sqr;

#endif