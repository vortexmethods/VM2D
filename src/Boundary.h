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
	
	/// Константная ссылка на указатель на профиль
	const std::unique_ptr<Airfoil>& afl;

	/// Константная ссылка на вихревой след
	const Wake& wake;

	/// \brief Размерность параметров каждого из слоев на каждой из панелей
	///
	/// Указывает, сколькими числами задается интенсивность каждого из слоев на каждой панели:
	/// - 1 --- одно число --- задается только среднее значение;
	/// - 2 --- два числа --- задается среднее значение и "наклон".
	int sheetDim;

	/// Слои на профиле
	Sheet sheets;

	/// Константная ссылка на параметры исполнения в параллельном режиме
	const Parallel& parallel;

	/// Константная ссылка на вектор из начал (и концов) панелей
	const std::vector<Point2D>& CC;

public:
	/// \brief Конструктор
	/// 
	/// \param[in] afl_ константная ссылка на указатель на профиль;
	/// \param[in] sheetDim_ размерность параметра слоя на каждой панели профиля (сколько чисел используется, чтобы описать слой на каждой панели);
	/// \param[in] wake_ константная ссылка на вихревой след;
	/// \param[in] parallel_ константная ссылка на параметры параллельного исполнения.
	Boundary(const std::unique_ptr<Airfoil>& afl_, int sheetDim_, const Wake& wake_, const Parallel& parallel_);
	
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

	/// \brief Генерация вектора влияния вихревого следа на профиль
	///
	/// Генерирует вектор влияния вихревого следа на профиль, используемый затем для расчета вектора правой части.
	/// 
	/// \param[out] wakeVelo ссылка на вектор влияния вихревого следа на профиль
	/// \warning Использует OMP, MPI
	virtual void GetWakeInfluence(std::vector<double>& wakeVelo) const = 0;
	
	/// \brief Вычисление скоростей вихрей в вихревом следе, вызываемых наличием завихренности и источников на профиле
	///
	/// Вычисляет скорости всех вихрей, уже имеющихся в вихревом следе, которые вызваны влиянием завихренности и источников на профиле 
	/// 
	/// \param[out] wakeVelo ссылка на вектор скоростей, которые приобретают вихри в вихревом следе из-за влияния завихренности и источников на профиле
	/// \param[in] dt шаг по времени (по умолчанию 1.0), на который нужно домножить скорость, чтобы вычислить не скорости, а перемещения вихрей
	/// \warning Использует OMP, MPI
	virtual void GetWakeVelocity(std::vector<Point2D>& wakeVelo, double dt = 1.0) const = 0;

	/// \brief Заполнение правой части
	///
	/// Заполняет блок правой части матрицы СЛАУ, обеспечивающей удовлетворение граничного условия, соответствующий данному профилю.
	/// 	
	/// \param[in] V0 константная ссылка на вектор набегающего потока
	/// \param[out] rhs ссылка на блок вектора правой части
	virtual void FillRhs(const Point2D& V0, Eigen::VectorXd& rhs) = 0;

	/// \brief Возврат размерности вектора решения 
	///
	/// (без учета регуляризирующей переменной)
	///
	/// \return размерность вектора решения
	virtual int GetUnknownsSize() const = 0;

	/// \brief Пересчет решения на интенсивность вихревого слоя
	///
	/// Приводит решение к интенсивности вихревого слоя и записывает его в sheets.freeVortexSheet:
	///
	/// - если неизвестное --- интенсивность вихря, то он "размазывается" по панели;
	/// - если неизвестное --- интенсивность слоя, то она передается непостредственно.
	///
	/// \param[in] sol вектор решения СЛАУ
	virtual void SolutionToFreeVortexSheet(const Eigen::VectorXd& sol) = 0;
};

#endif
