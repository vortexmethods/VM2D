/*!
\file
\brief Заголовочный файл с описанием класса Point2D
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef POINT2D_H_
#define POINT2D_H_

#include <cmath>

#include "mpi.h"

#include "numvector.h"

/*!
\brief Класс, опеделяющий двумерный вектор

Наследуется от numvector<double, 2>, имеет дополнительные возможности:
- поворота на заданный угол против часовой стрелки;
- генерируется MPI-описатель для возможности его пересылки как единичного объекта.

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class Point2D 
	: public numvector<double, 2>
{
public:	
	/// MPI-описатель типа
	static MPI_Datatype mpiPoint2D;
	
	/// Пустой конструктор
	Point2D() { };

	/// \brief Конструктор и приведение типа из numvector<double, 2>
	///
	/// \param[in] _r константная ссылка на копируемый объект типа numvector<double, 2>
	Point2D(const numvector<double, 2>& _r);
	
	/// \brief Конструктор копирования
	///
	/// \param[in] _r константная ссылка на копируемый вектор
	Point2D(const Point2D& _r);
	
	/// \brief Конструктор инициализации списком
	///
	/// \param[in] z константная ссылка на список инициализации из чисел типа double
	/// \warning Длина списка инициализации не проверяется, от него берутся только 2 первых элемента
	Point2D(const std::initializer_list<double>& z);
	
	/// Деструктор
	~Point2D() { };

	/// \brief Поворот вектора на произвольный угол против часовой стрелки (по умолчанию 90 градусов)
	///
	/// \param[in] angle угол поворота в радианах (по умолчанию \f$ \frac{\pi}{2} \f$)
	/// \return новый вектор, полученный поворотом старого 
	Point2D rotated(const double angle = 1.5707963267948966192313216916398) const;
	
	/// Cоздание MPI-описателя типа
	static void CreateMpiType();	
};


#endif
 