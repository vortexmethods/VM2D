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
	/// \param[out] time ссылка на промежуток времени --- пару чисел (время начала и время конца операции)
	virtual void GetHydroDynamForce(timePeriod& time) = 0;


	/// Генерация заголовка файла нагрузок
	/// \param[in] airfoilNumber номер профиля, для которого сохраняются силы
	void GenerateForcesHeader(int airfoilNumber)
	{
		std::stringstream forceFileName;
		forceFileName << passport.dir << "forces-airfoil-" << airfoilNumber << ".txt";

		std::ofstream newForcesFile(forceFileName.str());
		PrintLogoToTextFile(newForcesFile, forceFileName.str(), "Hydrodynamic loads for the airfoil " + passport.airfoilParams[airfoilNumber].fileAirfoil);

		PrintHeaderToTextFile(newForcesFile, "currentStep     currentTime     Fx     Fy");

		newForcesFile.close();
		newForcesFile.clear();
	}//GenerateForcesHeader(...)


	/// Сохранение строки со статистикой в файл нагрузок
	/// \param[in] currentStep номер текущего шага
	/// \param[in] airfoilNumber номер профиля, для которого сохраняются силы
	void GenerateForcesString(int currentStep, int airfoilNumber)
	{
		std::stringstream forceFileName;
		forceFileName << passport.dir << "forces-airfoil-" << airfoilNumber << ".txt";

		std::ofstream forcesFile(forceFileName.str(), std::ios::app);
		forcesFile << std::endl << currentStep << "	" << passport.physicalProperties.getCurrTime() << "	" << hydroDynamForce[0] << "	" << hydroDynamForce[1];
		forcesFile.close();
	}//GenerateForcesString(...)

};

#endif