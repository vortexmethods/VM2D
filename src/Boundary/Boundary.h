/*!
\file
\brief Заголовочный файл с описанием класса Boundary
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Sheet.h"
#include "Wake.h"

/*!
\brief Абстрактный класс, определяющий способ удовлетворения граничного условия на обтекаемом профиле
\author Марчевский Илья Константинович
\author Кузьмина Ксения Сергеевна
\author Рятина Евгения Павловна
\version 1.0
\date 1 декабря 2017 г.
*/

class Boundary
{
protected:
	/// Константная ссылка на паспорт
	const Passport& passport;

	/// Константная ссылка на вихревой след
	const Wake& wake;

	/// Константная ссылка на список указателей на все граничные условия
	const std::vector<std::unique_ptr<Boundary>>& allBoundary;

	/// Константная ссылка на параметры исполнения в параллельном режиме
	const Parallel& parallel;

public:	
	/// Константная ссылка на профиль
	const Airfoil& afl;

	/// Константная ссылка на вектор из начал (и концов) панелей
	const std::vector<Point2D>& CC;

	/// Виртуальный вихревой след конкретного профиля
	std::vector<Vortex2D> virtualWake;

	/// \brief Размерность параметров каждого из слоев на каждой из панелей
	///
	/// Указывает, сколькими числами задается интенсивность каждого из слоев на каждой панели:
	/// - 1 --- одно число --- задается только среднее значение;
	/// - 2 --- два числа --- задается среднее значение и "наклон".
	int sheetDim;

	/// Слои на профиле
	Sheet sheets;


	/// \brief Конструктор
	/// 
	/// \param[in] passport_ константная ссылка на паспорт расчета
	/// \param[in] afl_ константная ссылка на профиль
	/// \param[in] allBoundary_ константная ссылка на вектор из указателей на все граничные условия
	/// \param[in] sheetDim_ размерность параметра слоя на каждой панели профиля (сколько чисел используется, чтобы описать слой на каждой панели);
	/// \param[in] wake_ константная ссылка на вихревой след
	/// \param[in] parallel_ константная ссылка на параметры параллельного исполнения
	Boundary(const Passport& passport_, const Airfoil& afl_, const std::vector<std::unique_ptr<Boundary>>& allBoundary_, int sheetDim_, const Wake& wake_, const Parallel& parallel_);
	
	/// Деструктор
	virtual ~Boundary() { };
	
	/// \brief Генерация блока матрицы
	///
	/// Генерирует следующие компоненты матрицы:
	/// - диагональный блок матрицы --- влияние данного профиля на самого себя;
	/// - нижнюю строку для матрицы для данного профиля;
	/// - правый столбец матрицы для данного профиля.
	/// 
	/// \param[out] matr ссылка на генерируемую матрицу
	/// \param[out] lastLine ссылка на нижнюю строку
	/// \param[out] lactCol ссылка на правый столбец
	virtual void FillMatrixSelf(Eigen::MatrixXd& matr, Eigen::VectorXd& lastLine, Eigen::VectorXd& lactCol) = 0;



	/// \brief Генерация блока матрицы влияния от другого профиля того же типа
	///
	/// \todo ДОДЕЛАТЬ ОПИСАНИЕ
	virtual void FillMatrixFromOther(const Boundary& otherBoundary, Eigen::MatrixXd& matr) = 0;




	/// \brief Генерация вектора влияния вихревого следа на профиль
	///
	/// Генерирует вектор влияния вихревого следа на профиль, используемый затем для расчета вектора правой части.
	/// 
	/// \param[out] wakeVelo ссылка на вектор влияния вихревого следа на профиль
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	virtual void GetWakeInfluence(std::vector<double>& wakeVelo) const = 0;
	
	/// \brief Вычисление конвективных скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле как от слоев
	///
	/// Вычисляет конвективные скорости в наборе точек, которые вызваны влиянием завихренности и источников на профиле как от слоев
	/// 
	/// \param[in] points константная ссылка на набор точек, в которых вычисляются скорости
	/// \param[out] velo ссылка на вектор скоростей, которые приобретают точки из-за влияния завихренности и источников на профиле
	/// 
	/// \warning velo --- накапливается!
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	virtual void GetConvVelocityToSetOfPoints(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const = 0;

	/// \brief Вычисление конвективных скоростей в наборе точек, вызываемых наличием завихренности и источников на профиле как виртуальных вихрей
	///
	/// Вычисляет конвективные скорости в наборе точек, которые вызваны влиянием завихренности и источников на профиле как виртуальных вихрей
	/// 
	/// \param[in] points константная ссылка на набор точек, в которых вычисляются скорости
	/// \param[out] velo ссылка на вектор скоростей, которые приобретают точки из-за влияния завихренности и источников на профиле
	/// 
	/// \warning velo --- накапливается!
	/// \warning Использует OMP, MPI
	/// \ingroup Parallel
	virtual void GetConvVelocityToSetOfPointsFromVirtualVortexes(const std::vector<Vortex2D>& points, std::vector<Point2D>& velo) const = 0;

	
	/// \brief Заполнение правой части
	///
	/// Заполняет блок правой части матрицы СЛАУ, обеспечивающей удовлетворение граничного условия, соответствующий данному профилю.
	/// 	
	/// \param[in] V0 константная ссылка на вектор набегающего потока
	/// \param[out] rhs ссылка на блок вектора правой части
	/// \param[out] lastRhs указатель на последний элемент вектора правой части
	virtual void FillRhs(const Point2D& V0, Eigen::VectorXd& rhs, double* lastRhs) = 0;


	/// \brief Возврат размерности вектора решения 
	///
	/// (без учета регуляризирующей переменной)
	///
	/// \return размерность вектора решения
	virtual int GetUnknownsSize() const = 0;

	/// \brief Пересчет решения на интенсивность вихревого слоя и на рождаемые вихри на конкретном профиле
	///
	/// 1) Приводит решение к интенсивности вихревого слоя и записывает его в sheets.freeVortexSheet:
	///
	/// - если неизвестное --- интенсивность вихря, то он "размазывается" по панели;
	/// - если неизвестное --- интенсивность слоя, то она передается непостредственно.
	///
	/// 2) Приводит интенсивность вихревого слоя к рождаемым вихрям, а также вычисляет их положения
	/// \param[in] sol вектор решения СЛАУ
	virtual void SolutionToFreeVortexSheetAndVirtualVortex(const Eigen::VectorXd& sol) = 0;
};

#endif
