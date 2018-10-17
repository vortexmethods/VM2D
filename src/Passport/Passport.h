/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.4    |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2018/10/16     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2018 Ilia Marchevsky, Kseniia Kuzmina, Evgeniya Ryatina  |
*-----------------------------------------------------------------------------*
| File name: Passport.h                                                       |
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
\brief Заголовочный файл с описанием класса Passport и cоответствующими структурами

Помимо класса Passport определены следующие структуры:
- PhysicalProperties --- физические свойства задачи
- TimeDiscretizationProperties --- параметры процесса интегрирования по времени
- WakeDiscretizationProperties --- параметры дискретизации вихревого следа
- NumericalSchemes --- используемые численные схемы
- AirfoilParams --- параметры профиля

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/


#ifndef PASSPORT_H
#define PASSPORT_H

#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Point2D.h"
#include "LogStream.h"


/*!
\brief Структура, задающая физические свойства задачи

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
struct PhysicalProperties
{
// TODO: После реализации синхронизации паспортов сделать его private
public: 
	/// Текущее время
	mutable double currTime;

public:
	/// Плотность потока
	double rho;     

	/// Скоростью набегающего потока	
	Point2D vInf;

	/// Время разгона
	double timeAccel;

	/// Функция скорости набегающего потока с учетом разгона
	Point2D V0() const
	{
		return (currTime < timeAccel) ? static_cast<Point2D>(vInf * (currTime / timeAccel)) : vInf;
	};

	/// Коэффициент кинематической вязкости среды
	double nu;

	/// Возвращает текуще время
	double getCurrTime() const
	{
		return currTime;
	}

	/// Установка текущего времени
	void setCurrTime(double t_) const
	{
		currTime = t_;
	}

	/// Добавление шага к текущему времени
	void addCurrTime(double dt_) const
	{
		currTime += dt_;
	}



};//PhysicalProperties


/*!
\brief Структура, задающая параметры процесса интегрирования по времени 

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
struct TimeDiscretizationProperties
{
	/// Начальное время 
	double timeStart;    
	
	/// Конечное время
	double timeStop;
	
	/// Шаг по времени
	double dt;    
	
	/// Шаг сохранения кадров в текстовые файлы
	int saveTXT;
	
	/// Шаг сохранения кадров в бинарные файлы	
	int saveVTK; 

	/// Шаг вычисления и сохранения скорости и давления
	int saveVP;
};//TimeDiscretizationProperties


/*!
\brief Структура, задающая параметры параметры дискретизации вихревого следа

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
struct WakeDiscretizationProperties
{	
	/// Имя файла с начальным состоянием вихревого следа (без полного пути)
	std::string fileWake;
	
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

	/// Число вихрей, рождаемых на каждой панели профииля
	int vortexPerPanel;

};//WakeDiscretizationProperties


/*!
\brief Структура, задающая используемые численные схемы

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
struct NumericalSchemes
{
	//Решатель СЛАУ
	int linearSystemSolver;     

	//Метод вычисления скоростей вихрей и в отдельных точках
	int velocityComputation;		
	
	//Метод решения ОДУ движения вихрей
	int wakeMotionIntegrator;            
};//NumericalSchemes


/*!
\brief Структура, задающая параметры профиля

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
struct AirfoilParams
{
	/// Имя файла с начальным состоянием профилей (без полного пути)
	std::string fileAirfoil;
	
	/// Смещение центра масс (перенос профиля)
	Point2D basePoint;  
	
	/// Коэффициент масштабирования
	double scale;        
	
	/// Угол поворота (угол атаки)
	double angle;        

	/// Тип панелей (0 --- прямые)
	int panelsType;	
		
	/// Метод аппроксимации граничных условий
	int boundaryCondition;    

	/// Тип механической системы
	int mechanicalSystemType;
	std::string mechanicalSystem;
	std::string mechanicalSystemParameters;

};//AirfoilParams



/*!
\brief Класс, опеделяющий паспорт задачи

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.4
\date 16 октября 2018 г.
*/
class Passport
{
private:	
	/// \brief Считывание всех параметров расчета из соответствующих потоков
	/// 
	/// \param[in] mainStream ссылка на основной поток
	/// \param[in] mechanicsStream ссылка на поток со словарем механических систем
	/// \param[in] defaultStream ссылка на поток с параметрами по умолчанию
	/// \param[in] switcherStream ссылка на поток со значениями параметров-переключателей
	/// \param[in] varsStream ссылка на поток с параметрами конкретной задачи и переменными
	void GetAllParamsFromParser
		(
		std::istream& mainStream,
		std::istream& mechanicsStream,
		std::istream& defaultStream,
		std::istream& switcherStream,
		std::istream& varsStream
	);

	/// Поток для вывода логов и сообщений об ошибках
	mutable LogStream info;
	
	/// Печать всех параметров расчета в поток логов
	void PrintAllParams();

public:
	/// Рабочий каталог задачи
	std::string dir;

	/// Название задачи
	std::string problemName;

	/// Номер задачи
	size_t problemNumber;

	/// Каталог с файлами профилей
	std::string airfoilsDir;

	/// Каталог с файлами вихревых следов
	std::string wakesDir;

	/// Список структур с параметрами профилей 
	std::vector<AirfoilParams> airfoilParams;

	/// Структура с физическими свойствами задачи
	PhysicalProperties physicalProperties;

	/// Структура с параметрами процесса интегрирования по времени
	TimeDiscretizationProperties timeDiscretizationProperties;

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
		LogStream& infoStream,
		const std::string& _problemName,
		const size_t _problemNumber,
		const std::string& _filePassport,
		const std::string& _mechanics,
		const std::string& _defaults, 
		const std::string& _switchers,
		const std::vector<std::string>& vars
	);

	/// Деструктор
	virtual ~Passport(){ };			
};


#endif

