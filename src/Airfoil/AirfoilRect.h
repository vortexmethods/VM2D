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
	AirfoilRect(const Passport& passport_, const int numberInPassport_, const Parallel& parallel_)
		:Airfoil(passport_, numberInPassport_, parallel_)
	{ };

	/// Деструктор
	virtual ~AirfoilRect() { };

	//далее -- реализация виртуальной функции
	virtual void ReadFromFile(const std::string& dir);

	virtual void GetDiffVelocityToSetOfPointsAndViscousStresses(const std::vector<Vortex2D>& points, std::vector<double>& domainRadius, std::vector<Point2D>& velo, double epscol);
};

#endif