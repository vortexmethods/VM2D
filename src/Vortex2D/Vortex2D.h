/*!
\file
\brief Заголовочный файл с описанием класса Vortex2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef VORTEX2D_H_
#define VORTEX2D_H_

#include <iostream>

#include "Point2D.h"

/*!
\brief Класс, опеделяющий двумерный вихревой элемент
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class Vortex2D
{
private:
	/// Радиус-вектор вихря
	Point2D pos;	
	
	/// Циркуляция вихря
	double gam;	

public:
	/// MPI-описатель типа
	static MPI_Datatype mpiVortex2D;

	/// Пустой конструктор
	Vortex2D() {};

	/// \brief Конструктор инициализации
	///
	/// \param[in] _r константная ссылка на радиус-вектор положения вихря
	/// \param[in] _g циркуляция (интенсивность) вихря
	Vortex2D(const Point2D& _r, const double _g)
		: pos(_r), gam(_g)	{ }
	
	/// Деструктор
	~Vortex2D() {};

	/// \brief Функция для доступа к радиус-вектору вихря
	/// \return ссылка на радиус-вектор вихря
	Point2D& r() { return pos; }

	/// \brief Функция для доступа для чтения к радиус-вектору вихря
	/// \return константная ссылка на радиус-вектор вихря
	const Point2D& r() const { return pos; }

	/// \brief Функция для доступа к циркуляции вихря
	/// \return ссылка на циркуляцию вихря
	double& g() { return gam; }

	/// \brief Функция для доступа для чтения к циркуляции вихря
	/// \return константная ссылка на циркуляцию вихря
	const double& g() const { return gam; }

	/// Cоздание MPI-описателя типа
	static void CreateMpiType();	
};

#endif
 