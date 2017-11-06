/*!
\file
\brief Заголовочный файл с описанием класса AirfoilRect
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef AIRFOILRECT_H
#define AIRFOILRECT_H

#include "Airfoil.h"

/*!
\brief Класс, определяющий тип обтекаемого профиля

Тип профиля:
- профиль с прямолинейными панелями.

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.0
\date 1 декабря 2017 г.
*/
class AirfoilRect 
	: public Airfoil
{
public:

	/// Конструктор
	AirfoilRect() { };

	/// Деструктор
	virtual ~AirfoilRect() { };

	//далее -- реализация виртуальной функции
	virtual void ReadFromFile(const std::string& dir, const AirfoilParams& param);
};

#endif