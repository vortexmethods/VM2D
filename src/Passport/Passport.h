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
\version 1.0
\date 1 декабря 2017 г.
*/


#ifndef PASSPORT_H
#define PASSPORT_H

#include <memory>

#include "Preprocessor.h"
#include "StreamParser.h"


/*!
\brief Структура, задающая физические свойства задачи

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
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
		return (currTime < timeAccel) ? (vInf * (currTime / timeAccel)) : vInf;
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
\version 1.0
\date 1 декабря 2017 г.
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
	int deltacntText;
	
	/// Шаг сохранения кадров в бинарные файлы	
	int deltacntBinary;  
};//TimeDiscretizationProperties


/*!
\brief Структура, задающая параметры параметры дискретизации вихревого следа

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
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
	double distKill; 

	/// Расстояние, на которое рождаемый вихрь отодвигается от профиля
	double delta;
};//WakeDiscretizationProperties


/*!
\brief Структура, задающая используемые численные схемы

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
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
\version 1.0
\date 1 декабря 2017 г.
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
	int mechanicalSystem;
};//AirfoilParams



/*!
\brief Класс, опеделяющий паспорт задачи

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class Passport
{
private:	
	/// \brief Считывание всех параметров расчета из соответствующих потоков
	/// 
	/// \param[in] mainStream ссылка на основной поток
	/// \param[in] defaultStream ссылка на поток с параметрами по умолчанию
	/// \param[in] switcherStream ссылка на поток со значениями параметров-переключателей
	/// \param[in] varsStream ссылка на поток с параметрами конкретной задачи и переменными
	void GetAllParamsFromParser
	(
		std::istream& mainStream, 
		std::istream& defaultStream, 
		std::istream& switcherStream,
		std::istream& varsStream
	);
	
	/// Печать всех параметров расчета в поток логов
	void PrintAllParams();

public:
	/// Рабочий каталог задачи
	std::string dir;

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
	/// \param[in] _dir константная ссылка на рабочий каталог задачи
	/// \param[in] _filePassport константная ссылка на файл (без пути) с паспортом задачи
	/// \param[in] _defaults константная ссылка на имя файла (с путем) с параметрами по умолчанию
	/// \param[in] _switchers константная ссылка на имя файла (с путем) со значениями параметров-переключателей
	/// \param[in] vars константная ссылка на список переменных, заданных в виде строк
	Passport
	(
		const std::string& _dir, 
		const std::string& _filePassport, 
		const std::string& _defaults, 
		const std::string& _switchers,
		const std::vector<std::string>& vars);

	/// Деструктор
	virtual ~Passport(){ };			
};


#endif

