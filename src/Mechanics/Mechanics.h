/*!
\file
\brief Заголовочный файл с описанием класса Mechanics
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef MECHANICS_H
#define MECHANICS_H

#include "Velocity.h"

/*!
\brief Абстрактный класс, определяющий вид механической системы

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.0
\date 1 декабря 2017 г.
*/

class Mechanics
{
protected:
	/// Константная ссылка на Passport
	const Passport& passport;

	/// Ссылка на профиль
	Airfoil& afl;

	/// Константная ссылка на граничное условие
	const Boundary& boundary;

	/// Константная ссылка на параметры исполнения в параллельном режиме
	const Parallel& parallel;

	/// Константная ссылка на структуру с параметрами виртуального вихревого слоя для профиля
	const VortexesParams& virtVortParams;

public:
	/// Вектор гидродинамической силы, действующей на профиль
	Point2D hydroDynamForce;

	/// \brief Конструктор
	/// 
	/// \param[in] passport_ константная ссылка на паспорт
	/// \param[in] afl_ ссылка на профиль;
	/// \param[in] boundary_ константная ссылка на граничное условие;
	/// \param[in] virtVortParams_ константная ссылка на параметры виртуального вихревого следа для профиля;
	/// \param[in] parallel_ константная ссылка на параметры параллельного исполнения.
	Mechanics(const Passport& passport_, Airfoil& afl_, const Boundary& boundary_, const VortexesParams& virtVortParams_, const Parallel& parallel_)
		: passport(passport_), afl(afl_), boundary(boundary_), parallel(parallel_), virtVortParams(virtVortParams_) {};

	/// Деструктор
	virtual ~Mechanics() { };

	/// \brief Вычисление гидродинамической силы, действующей на профиль
	///
	/// \param[in] dt шаг по времени
	virtual void GetHydroDynamForce() = 0;

};

#endif