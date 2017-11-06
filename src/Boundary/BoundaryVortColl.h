/*!
\file
\brief Заголовочный файл с описанием класса BoundaryVortColl
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef BOUNDARYVORTCOLL_H
#define BOUNDARYVORTCOLL_H

#include "Boundary.h"

/*!
\brief Класс, определяющий способ удовлетворения граничного условия на обтекаемом профиле

Способ удовлетворения граничного условия:
- рождение вихрей в центрах панелей;
- условие коллокации в центрах панелей для касательной составляющей скорости.

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.0
\date 1 декабря 2017 г.
*/
class BoundaryVortColl : public Boundary
{
private:
	/// Контрольные точки (устанавливаются в центры панелей)
	std::vector<Point2D> KK;

	/// \brief Вычисляет скос на точку R от точки X
	///
	/// Вектор вихревого влияния на точку R от вихря единичной интенсивности, находящегося в точке X
	/// \n (сглаживание поля скоростей не производится)
	/// 
	/// \param[in] R точка наблюдения
	/// \param[in] X точка, где находится вихрь
	/// \return вектор скорости
	Point2D Skos(const Point2D& R, const Point2D& X);

public:
	/// \brief Конструктор
	/// 
	/// \param[in] afl_ константная ссылка на указатель на профиль
	/// \param[in] wake_ константная ссылка на вихревой след
	/// \param[in] parallel_ константная ссылка на параметры параллельного исполнения
	BoundaryVortColl(const std::unique_ptr<Airfoil>& afl_, const Wake& wake_, const Parallel& parallel_);
	
	/// Деструктор
	virtual ~BoundaryVortColl() { };

	//далее -- реализации виртуальных функций
	virtual void FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol);
	virtual void GetWakeInfluence(std::vector<double>& wakeVelo) const;
	virtual void GetWakeVelocity(std::vector<Point2D>& wakeVelo, double dt = 1.0) const;
	virtual void FillRhs(const Point2D& V0, Eigen::VectorXd& rhs);
	virtual int GetUnknownsSize() const;
	virtual void SolutionToFreeVortexSheet(const Eigen::VectorXd& sol);
};

#endif
