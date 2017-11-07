/*!
\file
\brief Заголовочный файл с описанием класса BoundaryVortColl
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef BOUNDARYMDV_H
#define BOUNDARYMDV_H

#include "Boundary.h"

/*!
\brief Класс, определяющий способ удовлетворения граничного условия на обтекаемом профиле

Способ удовлетворения граничного условия:
- рождение вихрей на концах панелей;
- условие коллокации в центрах панелей для нормальной составляющей скорости.

\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна

\version 1.0
\date 1 декабря 2017 г.
*/

class BoundaryMDV :
	public Boundary
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
	/// \param[in] afl_ константная ссылка на профиль;
	/// \param[in] allBoundry_ константная ссылка на вектор из указателей на все граничные условия
	/// \param[in] wake_ константная ссылка на вихревой след;
	/// \param[in] parallel_ константная ссылка на параметры параллельного исполнения.
	BoundaryMDV(const Passport& passport_, const Airfoil& afl_, const std::vector<std::unique_ptr<Boundary>>& allBoundary_, const Wake& wake_, const Parallel& parallel_);

	/// Деструктор
	virtual ~BoundaryMDV() { };

	//далее -- реализации виртуальных функций
	virtual void FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol);
	virtual void GetWakeInfluence(std::vector<double>& wakeVelo) const;
	virtual void GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const;
	virtual void FillRhs(const Point2D& V0, Eigen::VectorXd& rhs, double* lastRhs);
	virtual int GetUnknownsSize() const;
	virtual void SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol);
};

#endif