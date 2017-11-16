/*!
\file
\brief Заголовочный файл с описанием класса MechanicsRigidImmovable
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef MECHANICSRIGIDIMMOVABLE_H
#define MECHANICSRIGIDIMMOVABLE_H

#include "Mechanics.h"

/*!
\brief Класс, определяющий вид механической системы

Жесткое неподвижное тело

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.0
\date 1 декабря 2017 г.
*/

class MechanicsRigidImmovable :
	public Mechanics
{
public:
	/// \brief Конструктор
	/// 
	/// \param[in] passport_ константная ссылка на паспорт
	/// \param[in] afl_ ссылка на профиль;
	/// \param[in] boundary_ константная ссылка на граничное условие;
	/// \param[in] virtVortParams_ константная ссылка на параметры виртуального вихревого следа для профиля;
	/// \param[in] parallel_ константная ссылка на параметры параллельного исполнения.
	MechanicsRigidImmovable(const Passport& passport_, Airfoil& afl_, const Boundary& boundary_, const VortexesParams& virtVortParams_, const Parallel& parallel_) :
		Mechanics(passport_, afl_, boundary_, virtVortParams_, parallel_)
	{};

	/// Деструктор
	~MechanicsRigidImmovable() {};

	/// Вычисление гидродинамической силы, действующей на профиль
	virtual void GetHydroDynamForce(timePeriod& time);

	///
};

#endif