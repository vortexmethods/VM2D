/*!
\file
\brief Заголовочный файл с описанием класса Velocity
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef VELOCITY_H
#define VELOCITY_H

#include "Boundary.h"


/*!
\brief Структура, определяющая параметры виртуальных вихрей для отдельного профиля
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
struct VortexesParams
{
	/// Вектор конвективных скоростей виртуальных вихрей в следе
	std::vector<Point2D> convVelo;

	/// Вектор диффузионных скоростей виртуальных вихрей в следе
	std::vector<Point2D> diffVelo;

	/// Вектор характерных радиусов вихревых доменов (eps*) виртуальных вихрей
	std::vector<double> epsastWake;
};

/*!
\brief Абстрактный класс, определяющий способ вычисления скоростей
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/
class Velocity
{
protected:
	/// Константная ссылка на вихревой след
	const Wake& wake;

	/// Константная ссылка на параметры исполнения задачи в параллельном MPI-режиме
	const Parallel& parallel;	

	/// Константная ссылка на вектор указателей на граничные условия
	const std::vector<std::unique_ptr<Boundary>>& boundary;

public:
	/// Струтура, определяющая параметры вихрей в следе
	VortexesParams wakeVortexesParams;

	/// Вектор струтур, определяющий параметры виртуальных вихрей для профилей
	std::vector<VortexesParams> virtualVortexesParams;
	

	/// \brief Конструктор
	/// 
	/// \param[in] parallel_ константная ссылка на параметры исполнения задачи в параллельном MPI-режиме
	/// \param[in] wake_ константная ссылка на вихревой след	
	/// \param[in] boundary_ константная ссылка на вектор указателей на граничные условия 	
	Velocity(const Parallel& parallel_, const Wake& wake_, const std::vector<std::unique_ptr<Boundary>>& boundary_) :
		wake(wake_), parallel(parallel_), boundary(boundary_)
	{			
		virtualVortexesParams.resize(0);
	};
	
	/// \brief Вычисление конвективных скоростей и радиусов вихревых доменов в заданном наборе точек
	///
	/// \param[in] points константная ссылка на вектор из вихрей, в которых надо сосчитать конвективные скорости
	/// \param[out] velo ссылка на вектор скоростей в требуемых точках
	/// \param[out] domainRadius ссылка на вектор радиусов вихревых доменов
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	virtual void CalcConvVeloToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo, std::vector<double>& domainRadius) = 0;

	/// \brief Вычисление конвективных скоростей вихрей и виртуальных вихрей в вихревом следе
	///
	/// Вызывает 2 раза функцию CalcConvVeloToSetOfPoints
	///
	/// \warning скорости приплюсовываются к уже имеющимся
	void CalcConvVelo();

	/// \brief Вычисление диффузионных скоростей в заданном наборе точек
	///
	/// \param[in] points константная ссылка на вектор из вихрей, в которых надо сосчитать диффузионные скорости
	/// \param[in] domainRadius константная ссылка на вектор радиусов вихревых доменов
	/// \param[in] vortices константная ссылка на вектор из вихрей, от которых надо сосчитать влияния на points
	/// \param[out] velo ссылка на вектор скоростей в требуемых точках
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	virtual void CalcDiffVeloToSetOfPoints(const std::vector<Vortex2D>& points, const std::vector<double>& domainRadius, const std::vector<Vortex2D>& vortices, std::vector<Point2D>& velo) = 0;

	/// \brief Вычисление диффузионных скоростей вихрей и виртуальных вихрей в вихревом следе
	///
	/// Вызывает 4 раза функцию CalcDiffVeloToSetOfPoints
	///
	/// \param[in] nu кинематическая вязкость среды
	///
	/// \warning скорости приплюсовываются к уже имеющимся
	void CalcDiffVelo(double nu);

	/// Деструктор
	virtual ~Velocity() { };
};

#endif
