#ifndef GMRES2D_H_
#define GMRES2D_H_

#include <vector>
#include <Point2D.h>
#include <iostream>
#include "omp.h"

#include "Airfoil2D.h"
#include "World2D.h"

#include <string>


/// \brief Шаблонная функция вычисления евклидовой нормы вектора или списка
///
/// \param[in] b константная ссылка на вектор
template<typename T>
inline double norm(const std::vector<T>& b)
{
	double norm = 0;
//#ifndef OLD_OMP
//#pragma omp simd reduction(+:norm)
//#endif
	for (size_t i = 0; i < b.size(); i++)
		norm += (b[i] * b[i]);
	return sqrt(norm);
}


/// Шаблонная функция умножения числа на вектор
template<typename T>
inline std::vector<T> operator*(const T lambda, const std::vector<T>& x)
{
	std::vector<T> c(x);
	c.resize(x.size());

//#ifndef OLD_OMP
//#pragma omp simd
//#endif
	for (size_t i = 0; i < x.size(); ++i)
		c[i] *= lambda;
	return c;
}


/// \brief Контроль невязки после выполнения очередной итерации GMRES
///
/// \param[in,out] H ссылка на матрицу вращений Гивенса
/// \param[in] rhs константная ссылка на вектор правой части решаемой СЛАУ
/// \param[in,out] gs ссылка на правую часть решаемой СЛАУ
/// 
bool IterRot(std::vector<std::vector<double>>& H, const std::vector<double>& rhs, double& gs, std::vector<double>& c, std::vector<double>& s, int m, int n, double epsGMRES, int iter, bool residualShow);

void SolCircleRundirect(const std::vector<double>& A, const std::vector<double>& rhs, size_t startRow, size_t startRowReg, size_t np, std::vector<double>& res);

void SolMdirect(const std::vector<double>& A, const std::vector<double>& rhs, size_t startRow, size_t startRowReg, size_t np, std::vector<double>& res, bool lin);

void GMRES_Direct(
	const VM2D::World2D& W,
	int nAllVars,
	int nafl,
	const std::vector<double>& mtrDir,
	const std::vector<double>& rhsDir,
	const std::vector<int>& pos,
	const std::vector<int>& vsize,
	std::vector<std::vector<double>>& gam,
	std::vector<double>& R);


///////////////////// FAST //////////////////////

void SolCircleRun(std::vector<double>& AX, const std::vector<double>& rhs, const VM2D::Airfoil& afl, const std::vector<std::vector<Point2D>>& prea1, const std::vector<std::vector<Point2D>>& prec1, int p, int n);

void SolM(const VM2D::World2D& W, std::vector<double>& AX, const std::vector<double>& rhs,
	const std::vector<std::vector<Point2D>>& prea, const std::vector<std::vector<Point2D>>& prec, const std::vector<std::vector<Point2D>>& prea1, const std::vector<std::vector<Point2D>>& prec1, int p, int n, bool linScheme);

#ifdef USE_CUDA
void GMRES(const VM2D::World2D& W,
	std::vector<std::vector<double>>& X,
	std::vector<double>& R,
	const std::vector<std::vector<double>>& rhs,
	const std::vector<double> rhsReg,
	int& niter,
	bool linScheme);
#endif


#endif