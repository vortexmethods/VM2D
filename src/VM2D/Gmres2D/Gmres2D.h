/*--------------------------------*- VM2D -*-----------------*---------------*\
| ##  ## ##   ##  ####  #####   |                            | Version 1.14   |
| ##  ## ### ### ##  ## ##  ##  |  VM2D: Vortex Method       | 2026/03/06     |
| ##  ## ## # ##    ##  ##  ##  |  for 2D Flow Simulation    *----------------*
|  ####  ##   ##   ##   ##  ##  |  Open Source Code                           |
|   ##   ##   ## ###### #####   |  https://www.github.com/vortexmethods/VM2D  |
|                                                                             |
| Copyright (C) 2017-2026 I. Marchevsky, K. Sokol, E. Ryatina, A. Kolganova   |
*-----------------------------------------------------------------------------*
| File name: Gmres2D.h                                                        |
| Info: Source code of VM2D                                                   |
|                                                                             |
| This file is part of VM2D.                                                  |
| VM2D is free software: you can redistribute it and/or modify it             |
| under the terms of the GNU General Public License as published by           |
| the Free Software Foundation, either version 3 of the License, or           |
| (at your option) any later version.                                         |
|                                                                             |
| VM2D is distributed in the hope that it will be useful, but WITHOUT         |
| ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       |
| FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License       |
| for more details.                                                           |
|                                                                             |
| You should have received a copy of the GNU General Public License           |
| along with VM2D.  If not, see <http://www.gnu.org/licenses/>.               |
\*---------------------------------------------------------------------------*/

/*!
\file
\brief Заголовочный файл с функциями для метода GMRES
\author Марчевский Илья Константинович
\author Сокол Ксения Сергеевна
\author Рятина Евгения Павловна
\author Колганова Александра Олеговна
\author Кобзарь Дарья Юрьевна
\Version 1.14
\date 6 марта 2026 г.
*/


#ifndef GMRES2D_H
#define GMRES2D_H

#include <Point2D.h>
#include <iostream>
#include "omp.h"

#include <string>

namespace VM2D
{
	/*!
	\brief Структура, определяющий необходимые массивы для рализации метода прогонки

	\author Марчевский Илья Константинович
	\author Сокол Ксения Сергеевна
	\author Рятина Евгения Павловна
	\author Колганова Александра Олеговна
	\author Кобзарь Дарья Юрьевна

	\Version 1.14
	\date 6 марта 2026 г.
	*/
	struct sweepVectors
	{
		std::vector<double> alpha;
		std::vector<double> beta;
		std::vector<double> gamma;
		std::vector<double> delta;
		std::vector<double> phi;
		std::vector<double> xi;
		std::vector<double> psi;

		void resize(size_t n)
		{
			alpha.resize(n);
			beta.resize(n);
			gamma.resize(n);
			delta.resize(n);
			phi.resize(n);
			xi.resize(n);
			psi.resize(n);
		}
	};

	class Airfoil;
	class World2D;

	/// \brief Шаблонная функция вычисления евклидовой нормы вектора или списка
	///
	/// \param[in] b константная ссылка на вектор
	template<typename T>
	inline double norm(const std::vector<T>& b)
	{
		double norm = 0;
		#ifndef OLD_OMP
		#pragma omp simd reduction(+:norm)
		#endif
		for (size_t i = 0; i < b.size(); i++)
			norm += (b[i] * b[i]);
		return sqrt(norm);
	}

	/// \brief Шаблонная функция вычисления евклидовой нормы вектора или списка
	///
	/// \param[in] b константная ссылка на вектор
	template<typename T>
	inline double norm2(const std::vector<T>& b)
	{
		double norm2 = 0;
		//#ifndef OLD_OMP
		//#pragma omp simd reduction(+:norm)
		//#endif
		for (size_t i = 0; i < b.size(); i++)
			norm2 += (b[i] * b[i]);
		return norm2;
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
	bool IterRot(std::vector<std::vector<double>>& H, const double nrmRhs, double& gs, std::vector<double>& c, std::vector<double>& s, int m, int n, double epsGMRES, int iter, bool residualShow);

	void SolCircleRundirect(const std::vector<double>& A, const std::vector<double>& rhs, size_t startRow, size_t startRowReg, size_t np, std::vector<double>& res);

	void SolMdirect(const std::vector<double>& A, const std::vector<double>& rhs, size_t startRow, size_t startRowReg, size_t np, std::vector<double>& res, bool lin);

	void GMRES_Direct(
		const World2D& W,
		int nAllVars,
		int nafl,
		const std::vector<double>& mtrDir,
		const std::vector<double>& rhsDir,
		const std::vector<int>& pos,
		const std::vector<int>& vsize,
		std::vector<std::vector<double>>& gam,
		std::vector<double>& R);


	///////////////////// FAST //////////////////////

	void SolCircleRun(std::vector<double>& AX, const std::vector<double>& rhs, const Airfoil& afl, const std::vector<std::vector<Point2D>>& prea1, const std::vector<std::vector<Point2D>>& prec1, int p, int n, sweepVectors& sw);

	void SolM(std::vector<double>& AX, const std::vector<double>& rhs, const Airfoil& afl,
		const std::vector<std::vector<Point2D>>& prea, const std::vector<std::vector<Point2D>>& prec, int p, int n, bool linScheme, sweepVectors& sw);

#ifdef USE_CUDA
	void GMRES(const World2D& W,
		std::vector<std::vector<double>>& X,
		std::vector<double>& R,
		const std::vector<std::vector<double>>& rhs,
		const std::vector<double> rhsReg,
		int& niter,
		bool linScheme);
#endif


#endif

}