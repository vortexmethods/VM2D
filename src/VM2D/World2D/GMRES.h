#pragma once

#include <vector>
#include <Point2D.h>
#include <iostream>
#include "omp.h"

#include <string>

#define ADDOP(n) 


/// \brief Шаблонная функция вычисления евклидовой нормы вектора или списка
///
/// \param[in] b константная ссылка на вектор или список
template<typename T>
inline double norm(const T& b)
{
	double norm = 0;
//#ifndef OLD_OMP
//#pragma omp simd reduction(+:norm)
//#endif
	for (size_t i = 0; i < b.size(); i++)
		norm += (b[i] * b[i]);
	ADDOP(2 * b.size() + 1);
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
	ADDOP(x.size());
	return c;
}




bool IterRot(std::vector<std::vector<double>>& H, const std::vector<double>& rhs, double& gs, std::vector<double>& c, std::vector<double>& s, int m, int n, double epsGMRES, int iter, bool residualShow)
{
	bool fl;
	double buf;

	for (int i = 0; i < m - 1; ++i)
	{
		buf = H[i][m - 1];

		H[i][m - 1] = c[i] * buf + s[i] * H[i + 1][m - 1];
		H[i + 1][m - 1] = -s[i] * buf + c[i] * H[i + 1][m - 1];
		ADDOP(4);
	}

	double zn = sqrt(H[m - 1][m - 1] * H[m - 1][m - 1] + H[m][m - 1] * H[m][m - 1]);

	c.push_back(H[m - 1][m - 1] / zn);
	s.push_back(H[m][m - 1] / zn);
	gs *= -s[m - 1];
	H[m - 1][m - 1] = c[m - 1] * H[m - 1][m - 1] + s[m - 1] * H[m][m - 1];
	H[m][m - 1] = 0.0;
	ADDOP(8);

	if (residualShow)
		std::cout << "Iteration: " << iter << ", residual = " << (fabs(gs) / norm(rhs)) << std::endl;

	//double npprmRhs = norm(rhs);

	fl = ((fabs(gs) / norm(rhs)) < epsGMRES);

	return fl;
}



void SolCircleRundirect(const std::vector<double>& A, const std::vector<double>& rhs, size_t startRow, size_t startRowReg, size_t np, std::vector<double>& res)
{
	size_t nAllVars = rhs.size();
	std::vector<double> alpha(np + 2), beta(np + 2), gamma(np + 2);
	std::vector<double> u(np), v(np);

	//for (size_t i = 0; i < np; ++i)
	//	rhs[startRow + np + i] = i + 1;

	double zn = A[(startRow + np + 1) * nAllVars + (startRow + np + 1)];
	alpha[2] = -A[(startRow + np + 1) * nAllVars + (startRow + np + 2)] / zn;
	beta[2] = rhs[startRow + np + 1] / zn;
	gamma[2] = -A[(startRow + np + 1) * nAllVars + (startRow + np + 0)] / zn;
	ADDOP(3);

	for (size_t i = 2; i <= np; ++i)
	{
		zn = A[(startRow + np + (i % np)) * nAllVars + (startRow + np + (i % np))] + alpha[i] * A[(startRow + np + (i % np)) * nAllVars + (startRow + np + ((i - 1) % np))];
		alpha[i + 1] = -A[(startRow + np + (i % np)) * nAllVars + (startRow + np + ((i + 1) % np))] / zn;
		beta[i + 1] = (rhs[startRow + np + (i % np)] - beta[i] * A[(startRow + np + (i % np)) * nAllVars + (startRow + np + ((i - 1) % np))]) / zn;
		gamma[i + 1] = -(gamma[i] * A[(startRow + np + (i % np)) * nAllVars + (startRow + np + ((i - 1) % np))]) / zn;
		ADDOP(6);
	}

	u[np - 1] = beta[np];
	v[np - 1] = alpha[np] + gamma[np];
	for (size_t i = np - 2; i >= 1; --i)
	{
		u[i] = alpha[i + 1] * u[i + 1] + beta[i + 1];
		v[i] = alpha[i + 1] * v[i + 1] + gamma[i + 1];
		ADDOP(2);
	}


	res[startRow + np] = (beta[np + 1] + alpha[np + 1] * u[1]) / (1.0 - gamma[np + 1] - alpha[np + 1] * v[1]);
	for (size_t i = 1; i < np; ++i)
		res[startRow + np + i] = u[i] + res[startRow + np] * v[i];
	ADDOP(np + 2);


	//std::cout << "Hello!" << std::endl;
	return;

}




void SolMdirect(const std::vector<double>& A, const std::vector<double>& rhs, size_t startRow, size_t startRowReg, size_t np, std::vector<double>& res, bool lin)
{
	std::vector<double> alpha(np), beta(np), gamma(np), delta(np), phi(np), xi(np), psi(np);
	double yn, ynn;
	double lam1, lam2, mu1, mu2, xi1, xi2;
	double a;

	//double len2 = pnt[0].len * 2.0;

	size_t nAllVars = rhs.size();

	double izn = 1.0 / A[(startRow + 0) * nAllVars + (startRow + 0)];
	alpha[1] = -A[(startRow + 0) * nAllVars + (startRow + 1)] * izn;
	beta[1] = -A[(startRow + 0) * nAllVars + (startRow + np - 1)] * izn;
	gamma[1] = -A[(startRow + 0) * nAllVars + (startRowReg)] * izn;
	delta[1] = rhs[startRow] * izn;
	ADDOP(5);
	double zn;

	//std::cout << "coeff[1]" << std::endl;

	for (size_t i = 1; i < np - 1; ++i)
	{
		//std::cout << "i = " << i << std::endl;
		//std::cout << "A.size = " << A.size() << std::endl;
		//std::cout << "index = " << (startRow + i) * nAllVars + (startRow + i-1) << std::endl;
		a = A[(startRow + i) * nAllVars + (startRow + i - 1)];
		//std::cout << "a " << std::endl;

		zn = a * alpha[i] + A[(startRow + i) * nAllVars + (startRow + i)];
		//std::cout << "zn " << std::endl;

		alpha[i + 1] = -A[(startRow + i) * nAllVars + (startRow + i + 1)] / zn;
		//std::cout << "alpha " << std::endl;

		beta[i + 1] = -a * beta[i] / zn;
		//std::cout << "beta " << std::endl;

		gamma[i + 1] = -(A[(startRow + i) * nAllVars + (startRowReg)] + a * gamma[i]) / zn;
		//std::cout << "gamma " << std::endl;

		delta[i + 1] = (rhs[startRow + i] - a * delta[i]) / zn;
		//std::cout << "delta " << std::endl;
		ADDOP(8);
	}

	//std::cout << "coeff[n]" << std::endl;

	a = A[(startRow + (np - 2)) * nAllVars + (startRow + (np - 3))];
	zn = alpha[np - 2] * a + A[(startRow + (np - 2)) * nAllVars + (startRow + (np - 2))];

	phi[np - 1] = -(A[(startRow + np - 2) * nAllVars + (startRow + np - 1)] + a * beta[np - 2]) / zn;
	psi[np - 1] = -(A[(startRow + np - 2) * nAllVars + (startRowReg)] + gamma[np - 2] * a) / zn;
	xi[np - 1] = (rhs[startRow + np - 2] - delta[np - 2] * a) / zn;
	ADDOP(7);

	//std::cout << "phi[n]" << std::endl;

	for (size_t i = np - 2; i > 0; --i)
	{
		phi[i] = alpha[i] * phi[i + 1] + beta[i];
		psi[i] = alpha[i] * psi[i + 1] + gamma[i];
		xi[i] = alpha[i] * xi[i + 1] + delta[i];
		ADDOP(3);
	}


	//std::cout << "phi[1]" << std::endl;

	double e = A[(startRow + np - 1) * nAllVars + (startRow + 0)];
	a = A[(startRow + np - 1) * nAllVars + (startRow + np - 2)];
	lam1 = e * phi[1] + a * phi[np - 1] + A[(startRow + np - 1) * nAllVars + (startRow + np - 1)];
	mu1 = e * psi[1] + a * psi[np - 1] + A[(startRow + np - 1) * nAllVars + (startRowReg)];
	xi1 = rhs[startRow + np - 1] - e * xi[1] - a * xi[np - 1];
	ADDOP(6);

	//std::cout << "lam1" << std::endl;


	lam2 = mu2 = xi2 = 0.0;
	for (size_t j = 0; j < np - 1; ++j)
	{
		lam2 += phi[j + 1];
		mu2 += psi[j + 1];
		xi2 -= xi[j + 1];
	}
	lam2 += 1.0;

	//std::cout << "lam2" << std::endl;



//#ifndef linScheme
//	xi2 += rhs[n];
//#else
//	xi2 += rhs[2 * n];
//#endif

	xi2 += rhs[startRowReg];

	zn = lam1 * mu2 - lam2 * mu1;
	yn = (xi1 * mu2 - xi2 * mu1) / zn;
	ynn = -(xi1 * lam2 - xi2 * lam1) / zn;
	ADDOP(8);

	res[startRowReg] = ynn;

	//#ifndef linScheme
	//	AX[n] = ynn;
	//#else
	//	AX[2 * n] = ynn;
	//#endif

		//std::cout << "Reg" << std::endl;


		//std::cout << "res.size = " << res.size() << std::endl;
		//std::cout << "index = " << startRow + np - 1 << std::endl;
	res[startRow + np - 1] = yn;
	for (size_t i = np - 2; i + 1 > 0; --i)
	{
		//std::cout << "i = " << i << std::endl;
		//std::cout << "sizes = " << res.size() << " " << phi.size() << " " << psi.size() << " " << xi.size() << std::endl;
		//std::cout << "index = " << startRow + i << " " << i+1 << " " << i+1 << " " << i+1 << std::endl;
		res[startRow + i] = phi[i + 1] * yn + psi[i + 1] * ynn + xi[i + 1];
		ADDOP(2);
	}


	//std::cout << "Sol" << std::endl;

//Циклическая прогонка для линейной схемы
//#ifdef linScheme 
	if (lin)
		SolCircleRundirect(A, rhs, startRow, startRowReg, np, res);
//#endif 

	//std::cout << "Sweep finished!" << std::endl;
}

void GMRES_Direct(
	int nAllVars,
	int nafl,
	const std::vector<double>& mtrDir,
	const std::vector<double>& rhsDir,
	const std::vector<int>& n,
	const std::vector<int>& pos,
	const std::vector<int>& vsize,
	std::vector<std::vector<double>>& gam, 
	std::vector<double>& R,	
	int currentStep
	//const std::vector<double>& Matr,
	//const std::vector<std::vector<double>>& Len,
	//std::vector<std::vector<Point2D>>& Tau,
	//std::vector<std::vector<double>>& X, std::vector<double>& R, std::vector<std::vector<double>>& rhs, const std::vector<int>& n/*, const params& prm*/, double& timing2, double& timing3, int& niter
)
{
	double tt3;
	
	
	
	const size_t maxIter = nAllVars+1;
	std::vector<double> x(nAllVars, 0.0);
	std::vector<double> y(nAllVars, 0.0);
	std::vector<double> residual(nAllVars, 0.0);
	std::vector<double> r(nAllVars, 0.0);
	std::vector<double> w(nAllVars, 0.0);
	std::vector<double> wOld(nAllVars, 0.0);
	std::vector<double> g(nAllVars, 0.0);
	std::vector<double> c(nAllVars, 0.0);
	std::vector<double> s(nAllVars, 0.0);
	std::vector<double> hisGs(nAllVars + 1, 0.0);
	std::vector<double> v(nAllVars * (maxIter + 1), 0.0);

	std::vector<double> vecFromV(nAllVars, 0.0);

	std::vector<double> h((nAllVars + 1) * maxIter, 0.0);
	double beta = 0.0, gs, normb = 0.0;
	size_t m;

	std::cout << "memory allocated" << std::endl;

	tt3 = -omp_get_wtime();
	for (size_t i = 0; i < nAllVars; ++i)
		normb += rhsDir[i] * rhsDir[i];
	normb = sqrt(normb);

	ADDOP(nAllVars + 1);

	for (size_t i = 0; i < nAllVars; ++i)
	{
		residual[i] = rhsDir[i];
		//for (size_t j = 0; j < nAllVars; ++j)
/*!!!*/	//	residual[i] -= mtrDir[i * nAllVars + j] * x[j];

	}

	//cblas_dgemv(CblasRowMajor, CblasNoTrans, nAllVars, nAllVars, 1.0, mtrDir.data(), nAllVars, x.data(), 1, 0.0, wOld.data(), -1);
	for (int i = 0; i < nAllVars; ++i)
	{
		wOld[i] = 0.0;
		for (int j = 0; j < nAllVars; ++j)
		{
			wOld[i] += mtrDir[i * nAllVars + j] * x[j];
		}
	}



	//std::cout << "1-blas" << std::endl;

	ADDOP(nAllVars * nAllVars);

	//posI = 0;
#pragma omp parallel for 
	for (int aflI = 0; aflI < nafl; ++aflI)
	{
//		SolMdirect(mtrDir, residual, pos[aflI], nAllVars - nafl + aflI, n[aflI], r, (vsize[aflI] == 2*n[aflI]));
		r = residual;
	}

	//std::cout << "1-sweep" << std::endl;

	for (size_t i = 0; i < nAllVars; ++i)
		beta += r[i] * r[i];
	beta = sqrt(beta);
	gs = beta;
	g[0] = beta;
	ADDOP(nAllVars + 1);

	for (size_t i = 0; i < nAllVars; ++i)
		v[i * (maxIter + 1) + 0] = r[i] / beta;
	ADDOP(nAllVars);

	for (size_t j = 0; j < nAllVars; ++j)
	{
		//std::cout << "j = " << j << std::endl;
		double timer1 = omp_get_wtime();

		//for (size_t p = 0; p < nAllVars; ++p)
		//{
		//	wOld[p] = 0.0;
		//	for (size_t q = 0; q < nAllVars; ++q)
		//		wOld[p] += mtrDir[p * nAllVars + q] * v[q * (maxIter + 1) + j];
		//}//for p	

		for (size_t p = 0; p < nAllVars; ++p)
			vecFromV[p] = v[p * (maxIter + 1) + j];


		//cblas_dgemv(CblasRowMajor, CblasNoTrans, nAllVars, nAllVars, 1.0, mtrDir.data(), nAllVars, &v[0 * (maxIter + 1) + j], (nAllVars + 1), 0.0, wOld.data(), 1);
		//cblas_dgemv(CblasRowMajor, CblasNoTrans, nAllVars, nAllVars, 1.0, mtrDir.data(), nAllVars, vecFromV.data(), 1, 0.0, wOld.data(), 1);
		for (int i = 0; i < nAllVars; ++i)
		{
			wOld[i] = 0.0;
			for (int j = 0; j < nAllVars; ++j)
			{
				wOld[i] += mtrDir[i * nAllVars + j] * vecFromV[j];
			}
		}

		//{
		//	std::string sA = "d:\\Marchevsky\\VM2D\\build\\circle200\\dbg\\mul";
		//	std::string sB = "_old.txt";
		//	std::ofstream bufnewfile(sA + std::to_string(currentStep) + "_" + std::to_string(j) + sB);
		//	for (auto& q : vecFromV)
		//		bufnewfile << q << '\n';
		//	bufnewfile.close();
		//}
		
		//{			
		//	std::string sA = "d:\\Marchevsky\\VM2D\\build\\circle200\\dbg\\s";
		//	std::string sB = "_old.txt";
		//	std::ofstream bufnewfile(sA + std::to_string(currentStep) + "_" + std::to_string(j) + sB);
		//	for (auto& q : wOld)
		//		bufnewfile << q << '\n';
		//	bufnewfile.close();
		//}
		ADDOP(nAllVars * nAllVars);

		//std::cout << "2-blas" << std::endl;

		double timer2 = omp_get_wtime();

		//posI = 0;
#pragma omp parallel for
		for (int aflI = 0; aflI < nafl; ++aflI)
		{
//			SolMdirect(mtrDir, wOld, pos[aflI], nAllVars - nafl + aflI, n[aflI], w, (vsize[aflI] == 2 * n[aflI]));
			w = wOld;
		}

		//std::cout << "2-sweep" << std::endl;


		double timer3 = omp_get_wtime();

		for (size_t i = 0; i <= j; ++i)
		{
			h[i * maxIter + j] = 0.0;
			for (size_t q = 0; q < nAllVars; ++q)
				h[i * maxIter + j] += w[q] * v[q * (maxIter + 1) + i];
			ADDOP(nAllVars);

			for (size_t q = 0; q < nAllVars; ++q)
				w[q] -= h[i * maxIter + j] * v[q * (maxIter + 1) + i];
			ADDOP(nAllVars);
		}//for i

		double timer4 = omp_get_wtime();

		h[(j + 1) * maxIter + j] = 0.0;
		for (size_t p = 0; p < nAllVars; ++p)
			h[(j + 1) * maxIter + j] += w[p] * w[p];
		h[(j + 1) * maxIter + j] = sqrt(h[(j + 1) * maxIter + j]);
		ADDOP(nAllVars + 1);


		double timer5 = omp_get_wtime();

		//std::cout << "|w| = " << h[(j + 1) * maxIter + j] << std::endl;
		//std::cout << "j = " << j << ", gs = " << fabs(gs)/normb << std::endl;

		for (size_t p = 0; p < nAllVars; ++p)
			v[p * (maxIter + 1) + (j + 1)] = w[p] / h[(j + 1) * maxIter + j];
		ADDOP(nAllVars);

		double timer6 = omp_get_wtime();

		for (size_t i = 0; i < j; ++i)
		{
			double buf = h[i * maxIter + j];
			h[i * maxIter + j] = c[i] * buf + s[i] * h[(i + 1) * maxIter + j];
			h[(i + 1) * maxIter + j] = -s[i] * buf + c[i] * h[(i + 1) * maxIter + j];
			ADDOP(4);
		}

		double timer7 = omp_get_wtime();

		double zn = sqrt(h[j * maxIter + j] * h[j * maxIter + j] + h[(j + 1) * maxIter + j] * h[(j + 1) * maxIter + j]);
		c[j] = fabs(h[j * maxIter + j] / zn);
		s[j] = c[j] * h[(j + 1) * maxIter + j] / h[j * maxIter + j];

		gs *= -s[j];

		h[j * maxIter + j] = c[j] * h[j * maxIter + j] + s[j] * h[(j + 1) * maxIter + j];
		h[(j + 1) * maxIter + j] = 0.0;
		hisGs[j] = fabs(gs) / normb;
		ADDOP(10);

		printf("res = %e\n", fabs(gs));

		if (fabs(gs) / normb < 1e-14)
		{
			m = j + 1;
			break;
		}

		double timer8 = omp_get_wtime();

		//std::cout << "timer: " <<
		//	timer2 - timer1 << " " <<
		//	timer3 - timer2 << " " <<
		//	timer4 - timer3 << " " <<
		//	timer5 - timer4 << " " <<
		//	timer6 - timer5 << " " <<
		//	timer7 - timer6 << " " <<
		//	timer8 - timer7 << " " << std::endl;
	}//for j

	//Восстанавливаем решение
	for (size_t i = 0; i < m; ++i)
	{
		double buf = g[i];
		g[i] *= c[i];
		if (i < (m - 1))
			g[i + 1] = -s[i] * buf;

	}
	y[m - 1] = g[m - 1] / h[(m - 1) * maxIter + (m - 1)];
	ADDOP(2 * m + 1);

	for (size_t i = m - 2; i + 1 > 0; --i)
	{
		double sum = 0;
		for (size_t j = i + 1; j < m; ++j)
			sum += h[i * maxIter + j] * y[j];

		y[i] = (g[i] - sum) / h[i * maxIter + i];
		ADDOP(m + 1);
	}

	for (size_t p = 0; p < nAllVars; ++p)
	{
		for (size_t q = 0; q < m; ++q)
			x[p] += v[p * (maxIter + 1) + q] * y[q];
		ADDOP(m);
	}
	tt3 += omp_get_wtime();

	std::cout << "#iterations:     " << m << std::endl;
	//std::cout << "estimated error: " << fabs(gs) / normb << std::endl;

	//posI = 0;
	for (size_t aflI = 0; aflI < nafl; ++aflI)
	{
		for (size_t i = 0; i < vsize[aflI]; ++i)
			gam[aflI][i] = x[pos[aflI] + i];// *Len[aflI][i % n[aflI]];
		//posI += vsize[aflI];
	}

	for (size_t aflI = 0; aflI < nafl; ++aflI)
		R[aflI] = x[x.size() - nafl + aflI];
}



