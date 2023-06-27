/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.12   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2024/01/14     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2024 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Passport2D.h                                                     |
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
\brief Заголовочный файл с описанием класса Passport (двумерный) и cоответствующими структурами

Помимо класса Passport определены следующие структуры:
- PhysicalProperties --- физические свойства задачи
- TimeDiscretizationProperties --- параметры процесса интегрирования по времени
- WakeDiscretizationProperties --- параметры дискретизации вихревого следа
- NumericalSchemes --- используемые численные схемы
- AirfoilParams --- параметры профиля

\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\Version 1.12
\date 14 января 2024 г.
*/


#ifndef PASSPORT_H
#define PASSPORT_H

#include "PassportGen.h"
#include "Point2D.h"

namespace VM2D
{

	/*!
	\brief Структура, задающая физические свойства задачи

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
	\Version 1.12
	\date 14 января 2024 г.
	*/
	struct PhysicalProperties
	{
	private:
		const VMlib::TimeDiscretizationProperties& timeProp;

	public:
		/// Плотность потока
		double rho;

		/// Скоростью набегающего потока	
		Point2D vInf;

		/// Референсная скорость
		double vRef;

		/// Способ разгона потока
		std::pair<std::string, int> typeAccel;

		/// Время разгона потока
		double timeAccel;
		
		/// Функция-множитель, позволяющая моделировать разгон
		double accelCft() const;
		
		/// Функция скорости набегающего потока с учетом разгона
		Point2D V0() const
		{
			return static_cast<Point2D>(accelCft()*vInf);
		};

		/// Коэффициент кинематической вязкости среды
		double nu;

		/// Возвращает текуще время
		double getCurrTime() const
		{
			return timeProp.currTime;
		}

		/// Установка текущего времени
		void setCurrTime(double t_) const
		{
			timeProp.currTime = t_;
		}

		/// Добавление шага к текущему времени
		void addCurrTime(double dt_) const
		{
			timeProp.currTime += dt_;
		}

		PhysicalProperties(const VMlib::TimeDiscretizationProperties& timeProp_)
			:timeProp(timeProp_)
		{};

	};//PhysicalProperties


	/*!
	\brief Структура, задающая параметры параметры дискретизации вихревого следа

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
	\Version 1.12
	\date 14 января 2024 г.
	*/
	struct WakeDiscretizationProperties
	{
		/// Радиус вихря
		double eps;

		/// Квадрат радиуса вихря
		double eps2;

		/// Радиус коллапса
		double epscol;
		
		/// Расстояние от центра самого подветренного (правого) профиля, на котором вихри уничтожаются
		double distFar;

		/// Расстояние, на которое рождаемый вихрь отодвигается от профиля
		double delta;

		/// Минимальное число вихрей, рождаемых на каждой панели профииля
		int minVortexPerPanel;

		/// Максимально допустимая циркуляция вихря
		double maxGamma;

		/// Имя файла с начальным состоянием вихревого следа (без полного пути)
		std::string fileWake;

		/// Имя файла с положениями источников (без полного пути)
		std::string fileSource;

		/// Функция минимально возможного значения для epsAst
		double getMinEpsAst() const
		{
			return 2.0 * epscol;
		};

	};//WakeDiscretizationProperties


	/*!
	\brief Структура, задающая используемые численные схемы

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
    \author Колганова Александра Олеговна
	\Version 1.12
	\date 14 января 2024 г.
	*/
	struct NumericalSchemes
	{
		//Решатель СЛАУ
		std::pair<std::string, int> linearSystemSolver;

		//Метод вычисления скоростей вихрей и в отдельных точках
		std::pair<std::string, int> velocityComputation;

		//Метод решения ОДУ движения вихрей
		//std::pair<std::string, int> wakeMotionIntegrator;

		/// Метод аппроксимации граничных условий
		std::pair<std::string, int> boundaryCondition;
	};//NumericalSchemes


	/*!
	\brief Структура, задающая параметры профиля

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
	\Version 1.12
	\date 14 января 2024 г.
	*/
	struct AirfoilParams
	{
		/// Имя файла с начальным состоянием профилей (без полного пути)
		std::string fileAirfoil;

		/// Желаемое число панелей для разбиения геометрии
		size_t requiredNPanels;

		/// Смещение центра масс (перенос профиля)
		Point2D basePoint;

		/// Коэффициент масштабирования
		Point2D scale;

		/// Хорда
		double chord;

		/// Угол поворота (угол атаки)
		double angle;

		/// Присоединенная масса
		Point2D addedMass;

		/// Признак разворота нормалей (для расчета внутреннего течения)
		bool inverse;

		/// Тип механической системы
		int mechanicalSystemType;
		std::string mechanicalSystem;
		std::string mechanicalSystemParameters;

	};//AirfoilParams



	/*!
	\brief Класс, опеделяющий паспорт двумерной задачи

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
	\Version 1.12
	\date 14 января 2024 г.
	*/
	class Passport : public VMlib::PassportGen
	{
	private:
		//далее -- реализации виртуальных функций
		virtual void GetAllParamsFromParser
		(
			std::istream& mainStream,
			std::istream& mechanicsStream,
			std::istream& defaultStream,
			std::istream& switcherStream,
			std::istream& varsStream
		) override;
		virtual void PrintAllParams() override;

	public:
		/// Каталог с файлами профилей
		std::string airfoilsDir;

		/// Каталог с файлами вихревых следов
		std::string wakesDir;

		/// Список структур с параметрами профилей 
		std::vector<AirfoilParams> airfoilParams;

		/// Признак работы в "географической" системе координат
		bool geographicalAngles;

		/// Признак поворота вычисляемых сил в профильную систему координат
		bool rotateForces;

		/// Признак вычисления коэффициентов вместо сил
		bool calcCoefficients;

		/// Угол поворота точек VP
		double rotateAngleVpPoints;


		/// Структура с физическими свойствами задачи
		PhysicalProperties physicalProperties;

		/// Структура с параметрами дискретизации вихревого следа
		WakeDiscretizationProperties wakeDiscretizationProperties;

		/// Структура с используемыми численными схемами
		NumericalSchemes numericalSchemes;

		/// \brief Конструктор
		///
		/// Осуществляет чтение всех данных из соответствующих потоков, полностью инициализирует паспорт
		///
		/// \param[in, out] infoStream базовый поток для вывода логов
		/// \param[in] _problemName константная ссылка наназвание задачи
		/// \param[in] _problemNumber номер (по счету) решаемой задачи
		/// \param[in] _filePassport константная ссылка на файл (без пути) с паспортом задачи
		/// \param[in] _mechanics константная ссылка на файл (c путем) со словарем механических систем
		/// \param[in] _defaults константная ссылка на имя файла (с путем) с параметрами по умолчанию
		/// \param[in] _switchers константная ссылка на имя файла (с путем) со значениями параметров-переключателей
		/// \param[in] vars константная ссылка на список переменных, заданных в виде строк
		Passport
		(
			VMlib::LogStream& infoStream,
			const std::string& _problemName,
			const size_t _problemNumber,
			const std::string& _filePassport,
			const std::string& _mechanics,
			const std::string& _defaults,
			const std::string& _switchers,
			const std::vector<std::string>& vars
		);

		/// Деструктор
		virtual ~Passport() { };
	};

}//namespace VM2D

#endif

