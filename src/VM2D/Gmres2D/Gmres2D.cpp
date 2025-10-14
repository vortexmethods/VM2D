#include "Gmres2D.h"

#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "Passport2D.h"
#include "Preprocessor.h"
#include "StreamParser.h"
#include "Velocity2D.h"
#include "WakeDataBase2D.h"
#include "World2D.h"
#include "Wake2D.h"


using namespace VM2D;

bool IterRot(std::vector<std::vector<double>>& H, const std::vector<double>& rhs, double& gs, std::vector<double>& c, std::vector<double>& s, int m, int n, double epsGMRES, int iter, bool residualShow)
{
	bool fl;
	double buf;

	for (int i = 0; i < m - 1; ++i)
	{
		buf = H[i][m - 1];

		H[i][m - 1] = c[i] * buf + s[i] * H[i + 1][m - 1];
		H[i + 1][m - 1] = -s[i] * buf + c[i] * H[i + 1][m - 1];
	}

	double zn = sqrt(H[m - 1][m - 1] * H[m - 1][m - 1] + H[m][m - 1] * H[m][m - 1]);

	c.push_back(H[m - 1][m - 1] / zn);
	s.push_back(H[m][m - 1] / zn);
	gs *= -s[m - 1];
	H[m - 1][m - 1] = c[m - 1] * H[m - 1][m - 1] + s[m - 1] * H[m][m - 1];
	H[m][m - 1] = 0.0;

	if (residualShow)
		std::cout << "Iteration: " << iter << ", residual = " << (fabs(gs) / norm(rhs)) << std::endl;

	fl = ((fabs(gs) / norm(rhs)) < epsGMRES);
	return fl;
}


void SolCircleRundirect(const std::vector<double>& A, const std::vector<double>& rhs, size_t startRow, size_t startRowReg, size_t np, std::vector<double>& res)
{
	size_t nAllVars = rhs.size();
	std::vector<double> alpha(np + 2), beta(np + 2), gamma(np + 2);
	std::vector<double> u(np), v(np);

	double zn = A[(startRow + np + 1) * nAllVars + (startRow + np + 1)];
	alpha[2] = -A[(startRow + np + 1) * nAllVars + (startRow + np + 2)] / zn;
	beta[2] = rhs[startRow + np + 1] / zn;
	gamma[2] = -A[(startRow + np + 1) * nAllVars + (startRow + np + 0)] / zn;

	for (size_t i = 2; i <= np; ++i)
	{
		zn = A[(startRow + np + (i % np)) * nAllVars + (startRow + np + (i % np))] + alpha[i] * A[(startRow + np + (i % np)) * nAllVars + (startRow + np + ((i - 1) % np))];
		alpha[i + 1] = -A[(startRow + np + (i % np)) * nAllVars + (startRow + np + ((i + 1) % np))] / zn;
		beta[i + 1] = (rhs[startRow + np + (i % np)] - beta[i] * A[(startRow + np + (i % np)) * nAllVars + (startRow + np + ((i - 1) % np))]) / zn;
		gamma[i + 1] = -(gamma[i] * A[(startRow + np + (i % np)) * nAllVars + (startRow + np + ((i - 1) % np))]) / zn;
	}

	u[np - 1] = beta[np];
	v[np - 1] = alpha[np] + gamma[np];
	for (size_t i = np - 2; i >= 1; --i)
	{
		u[i] = alpha[i + 1] * u[i + 1] + beta[i + 1];
		v[i] = alpha[i + 1] * v[i + 1] + gamma[i + 1];
	}

	res[startRow + np] = (beta[np + 1] + alpha[np + 1] * u[1]) / (1.0 - gamma[np + 1] - alpha[np + 1] * v[1]);
	for (size_t i = 1; i < np; ++i)
		res[startRow + np + i] = u[i] + res[startRow + np] * v[i];

	return;
}




void SolMdirect(const std::vector<double>& A, const std::vector<double>& rhs, size_t startRow, size_t startRowReg, size_t np, std::vector<double>& res, bool lin)
{
	std::vector<double> alpha(np), beta(np), gamma(np), delta(np), phi(np), xi(np), psi(np);
	double yn, ynn;
	double lam1, lam2, mu1, mu2, xi1, xi2;
	double a;

	size_t nAllVars = rhs.size();

	double izn = 1.0 / A[(startRow + 0) * nAllVars + (startRow + 0)];
	alpha[1] = -A[(startRow + 0) * nAllVars + (startRow + 1)] * izn;
	beta[1] = -A[(startRow + 0) * nAllVars + (startRow + np - 1)] * izn;
	gamma[1] = -A[(startRow + 0) * nAllVars + (startRowReg)] * izn;
	delta[1] = rhs[startRow] * izn;
	double zn;


	for (size_t i = 1; i < np - 1; ++i)
	{
		a = A[(startRow + i) * nAllVars + (startRow + i - 1)];
		zn = a * alpha[i] + A[(startRow + i) * nAllVars + (startRow + i)];
		alpha[i + 1] = -A[(startRow + i) * nAllVars + (startRow + i + 1)] / zn;
		beta[i + 1] = -a * beta[i] / zn;
		gamma[i + 1] = -(A[(startRow + i) * nAllVars + (startRowReg)] + a * gamma[i]) / zn;
		delta[i + 1] = (rhs[startRow + i] - a * delta[i]) / zn;
	}

	a = A[(startRow + (np - 2)) * nAllVars + (startRow + (np - 3))];
	zn = alpha[np - 2] * a + A[(startRow + (np - 2)) * nAllVars + (startRow + (np - 2))];

	phi[np - 1] = -(A[(startRow + np - 2) * nAllVars + (startRow + np - 1)] + a * beta[np - 2]) / zn;
	psi[np - 1] = -(A[(startRow + np - 2) * nAllVars + (startRowReg)] + gamma[np - 2] * a) / zn;
	xi[np - 1] = (rhs[startRow + np - 2] - delta[np - 2] * a) / zn;

	for (size_t i = np - 2; i > 0; --i)
	{
		phi[i] = alpha[i] * phi[i + 1] + beta[i];
		psi[i] = alpha[i] * psi[i + 1] + gamma[i];
		xi[i] = alpha[i] * xi[i + 1] + delta[i];
	}


	double e = A[(startRow + np - 1) * nAllVars + (startRow + 0)];
	a = A[(startRow + np - 1) * nAllVars + (startRow + np - 2)];
	lam1 = e * phi[1] + a * phi[np - 1] + A[(startRow + np - 1) * nAllVars + (startRow + np - 1)];
	mu1 = e * psi[1] + a * psi[np - 1] + A[(startRow + np - 1) * nAllVars + (startRowReg)];
	xi1 = rhs[startRow + np - 1] - e * xi[1] - a * xi[np - 1];

	lam2 = mu2 = xi2 = 0.0;
	for (size_t j = 0; j < np - 1; ++j)
	{
		lam2 += phi[j + 1];
		mu2 += psi[j + 1];
		xi2 -= xi[j + 1];
	}
	lam2 += 1.0;

	xi2 += rhs[startRowReg];

	zn = lam1 * mu2 - lam2 * mu1;
	yn = (xi1 * mu2 - xi2 * mu1) / zn;
	ynn = -(xi1 * lam2 - xi2 * lam1) / zn;

	res[startRowReg] = ynn;

	res[startRow + np - 1] = yn;
	for (size_t i = np - 2; i + 1 > 0; --i)
		res[startRow + i] = phi[i + 1] * yn + psi[i + 1] * ynn + xi[i + 1];

//Циклическая прогонка для линейной схемы
	if (lin)
		SolCircleRundirect(A, rhs, startRow, startRowReg, np, res);
}



void GMRES_Direct(
	const VM2D::World2D& W,
	int nAllVars,
	int nafl,
	const std::vector<double>& mtrDir,
	const std::vector<double>& rhsDir,
	const std::vector<int>& pos,
	const std::vector<int>& vsize,
	std::vector<std::vector<double>>& gam,
	std::vector<double>& R)
{
	const size_t maxIter = nAllVars + 1;
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

	for (size_t i = 0; i < nAllVars; ++i)
		normb += rhsDir[i] * rhsDir[i];
	normb = sqrt(normb);

	for (size_t i = 0; i < nAllVars; ++i)	
		residual[i] = rhsDir[i];

#pragma omp parallel for 
	for (int aflI = 0; aflI < nafl; ++aflI)
	{
		// SolMdirect(mtrDir, residual, pos[aflI], nAllVars - nafl + aflI, W.getAirfoil(aflI).getNumberOfPanels(), r, (vsize[aflI] == 2*n[aflI]));
		r = residual;
	}

	for (size_t i = 0; i < nAllVars; ++i)
		beta += r[i] * r[i];
	beta = sqrt(beta);
	gs = beta;
	g[0] = beta;

	for (size_t i = 0; i < nAllVars; ++i)
		v[i * (maxIter + 1) + 0] = r[i] / beta;

	for (size_t j = 0; j < nAllVars; ++j)
	{
		double timer1 = omp_get_wtime();

		for (size_t p = 0; p < nAllVars; ++p)
			vecFromV[p] = v[p * (maxIter + 1) + j];

		//cblas_dgemv(CblasRowMajor, CblasNoTrans, nAllVars, nAllVars, 1.0, mtrDir.data(), nAllVars, &v[0 * (maxIter + 1) + j], (nAllVars + 1), 0.0, wOld.data(), 1);
		//cblas_dgemv(CblasRowMajor, CblasNoTrans, nAllVars, nAllVars, 1.0, mtrDir.data(), nAllVars, vecFromV.data(), 1, 0.0, wOld.data(), 1);
		for (int i = 0; i < nAllVars; ++i)
		{
			wOld[i] = 0.0;
			for (int jj = 0; jj < nAllVars; ++jj)
			{
				wOld[i] += mtrDir[i * nAllVars + jj] * vecFromV[jj];
			}
		}

		double timer2 = omp_get_wtime();

#pragma omp parallel for
		for (int aflI = 0; aflI < nafl; ++aflI)
		{
			// SolMdirect(mtrDir, wOld, pos[aflI], nAllVars - nafl + aflI, W.getAirfoil(aflI).getNumberOfPanels(), w, (vsize[aflI] == 2 * n[aflI]));
			w = wOld;
		}

		for (size_t i = 0; i <= j; ++i)
		{
			h[i * maxIter + j] = 0.0;
			for (size_t q = 0; q < nAllVars; ++q)
				h[i * maxIter + j] += w[q] * v[q * (maxIter + 1) + i];

			for (size_t q = 0; q < nAllVars; ++q)
				w[q] -= h[i * maxIter + j] * v[q * (maxIter + 1) + i];
		}//for i

		h[(j + 1) * maxIter + j] = 0.0;
		for (size_t p = 0; p < nAllVars; ++p)
			h[(j + 1) * maxIter + j] += w[p] * w[p];
		h[(j + 1) * maxIter + j] = sqrt(h[(j + 1) * maxIter + j]);

		for (size_t p = 0; p < nAllVars; ++p)
			v[p * (maxIter + 1) + (j + 1)] = w[p] / h[(j + 1) * maxIter + j];

		for (size_t i = 0; i < j; ++i)
		{
			double buf = h[i * maxIter + j];
			h[i * maxIter + j] = c[i] * buf + s[i] * h[(i + 1) * maxIter + j];
			h[(i + 1) * maxIter + j] = -s[i] * buf + c[i] * h[(i + 1) * maxIter + j];
		}


		double zn = sqrt(h[j * maxIter + j] * h[j * maxIter + j] + h[(j + 1) * maxIter + j] * h[(j + 1) * maxIter + j]);
		c[j] = fabs(h[j * maxIter + j] / zn);
		s[j] = c[j] * h[(j + 1) * maxIter + j] / h[j * maxIter + j];

		gs *= -s[j];

		h[j * maxIter + j] = c[j] * h[j * maxIter + j] + s[j] * h[(j + 1) * maxIter + j];
		h[(j + 1) * maxIter + j] = 0.0;
		hisGs[j] = fabs(gs) / normb;

		if ((fabs(gs) / normb) < W.getPassport().numericalSchemes.gmresEps)
		{
			m = j + 1;
			break;
		}
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

	for (size_t i = m - 2; i + 1 > 0; --i)
	{
		double sum = 0;
		for (size_t j = i + 1; j < m; ++j)
			sum += h[i * maxIter + j] * y[j];

		y[i] = (g[i] - sum) / h[i * maxIter + i];
	}

	for (size_t p = 0; p < nAllVars; ++p)
	{
		for (size_t q = 0; q < m; ++q)
			x[p] += v[p * (maxIter + 1) + q] * y[q];
	}

	std::cout << "#iterations:     " << m << std::endl;

	for (size_t aflI = 0; aflI < nafl; ++aflI)
	{
		for (size_t i = 0; i < vsize[aflI]; ++i)
			gam[aflI][i] = x[pos[aflI] + i];
	}

	for (size_t aflI = 0; aflI < nafl; ++aflI)
		R[aflI] = x[x.size() - nafl + aflI];
}









/////////////////////// FAST ///////////////////////////



void SolCircleRun(std::vector<double>& AX, const std::vector<double>& rhs, const Airfoil& afl, const std::vector<std::vector<Point2D>>& prea1, const std::vector<std::vector<Point2D>>& prec1, int p, int n)
{
	std::vector<double> alpha(n), beta(n), delta(n), phi(n), xi(n);
	double yn, a;

	double len24 = afl.len[0] * 24.0;

	alpha[1] = (afl.tau[0] & prec1[p][0]) * len24;
	beta[1] = (afl.tau[0] & prea1[p][0]) * len24;
	delta[1] = -len24 * rhs[n + 0];
	double zn;

	for (int i = 1; i < n - 1; ++i)
	{
		a = (afl.tau[i] & prea1[p][i]);
		zn = a * alpha[i] - 1.0 / (24.0 * afl.len[i]);
		alpha[i + 1] = -(afl.tau[i] & prec1[p][i]) / zn;
		beta[i + 1] = -a * beta[i] / zn;
		delta[i + 1] = (rhs[n + i] - a * delta[i]) / zn;
	}

	a = afl.tau[n - 2] & prea1[p][n - 2];
	zn = alpha[n - 2] * a - 1.0 / (24.0 * afl.len[n - 2]);
	phi[n - 1] = -((afl.tau[n - 2] & prec1[p][n - 2]) + beta[n - 2] * a) / zn;
	xi[n - 1] = (rhs[n + n - 2] - delta[n - 2] * a) / zn;

	for (int i = n - 2; i > 0; --i)
	{
		phi[i] = alpha[i] * phi[i + 1] + beta[i];
		xi[i] = alpha[i] * xi[i + 1] + delta[i];
	}


	double e = (afl.tau[n - 1] & prec1[p][n - 1]);
	a = (afl.tau[n - 1] & prea1[p][n - 1]);

	yn = (rhs[n + n - 1] - e * xi[1] - a * xi[n - 1]) / (e * phi[1] + a * phi[n - 1] - 1.0 / (24.0 * afl.len[n - 1]));

	AX[n + n - 1] = yn;
	for (int i = n - 2; i >= 0; --i)
		AX[n + i] = phi[i + 1] * yn + xi[i + 1];
}



void SolM(const World2D& W, std::vector<double>& AX, const std::vector<double>& rhs,
	const std::vector<std::vector<Point2D>>& prea, const std::vector<std::vector<Point2D>>& prec, const std::vector<std::vector<Point2D>>& prea1, const std::vector<std::vector<Point2D>>& prec1, int p, int n, bool linScheme)
{
	const Airfoil& afl = W.getAirfoil(p);

	std::vector<double> alpha(n), beta(n), gamma(n), delta(n), phi(n), xi(n), psi(n);
	double yn, ynn;
	double lam1, lam2, mu1, mu2, xi1, xi2;
	double a;

	double len2 = afl.len[0] * 2.0;

	alpha[1] = (afl.tau[0] & prec[p][0]) * len2;
	beta[1] = (afl.tau[0] & prea[p][0]) * len2;
	gamma[1] = len2;
	delta[1] = -len2 * rhs[0];

	double zn;

	for (int i = 1; i < n - 1; ++i)
	{
		a = (afl.tau[i] & prea[p][i]);
		zn = a * alpha[i] - 0.5 / afl.len[i];
		alpha[i + 1] = -(afl.tau[i] & prec[p][i]) / zn;
		beta[i + 1] = -a * beta[i] / zn;
		gamma[i + 1] = -(1.0 + a * gamma[i]) / zn;
		delta[i + 1] = (rhs[i] - a * delta[i]) / zn;
	}

	a = afl.tau[n - 2] & prea[p][n - 2];
	zn = alpha[n - 2] * a - 0.5 / afl.len[n - 2];
	phi[n - 1] = -((afl.tau[n - 2] & prec[p][n - 2]) + beta[n - 2] * a) / zn;
	psi[n - 1] = -(1.0 + gamma[n - 2] * a) / zn;
	xi[n - 1] = (rhs[n - 2] - delta[n - 2] * a) / zn;

	for (int i = n - 2; i > 0; --i)
	{
		phi[i] = alpha[i] * phi[i + 1] + beta[i];
		psi[i] = alpha[i] * psi[i + 1] + gamma[i];
		xi[i] = alpha[i] * xi[i + 1] + delta[i];
	}

	double e = (afl.tau[n - 1] & prec[p][n - 1]);
	a = (afl.tau[n - 1] & prea[p][n - 1]);
	lam1 = e * phi[1] + a * phi[n - 1] - 0.5 / afl.len[n - 1];
	mu1 = e * psi[1] + a * psi[n - 1] + 1.0;
	xi1 = rhs[n - 1] - e * xi[1] - a * xi[n - 1];
	lam2 = mu2 = xi2 = 0.0;
	for (int j = 0; j < n - 1; ++j)
	{
		lam2 += phi[j + 1];
		mu2 += psi[j + 1];
		xi2 -= xi[j + 1];
	}
	lam2 += 1.0;

	if (!linScheme)
		xi2 += rhs[n];
	else if (linScheme)
		xi2 += rhs[2 * n];

	zn = lam1 * mu2 - lam2 * mu1;
	yn = (xi1 * mu2 - xi2 * mu1) / zn;
	ynn = -(xi1 * lam2 - xi2 * lam1) / zn;

	if (!linScheme)
		AX[n] = ynn;
	else if (linScheme)
		AX[2 * n] = ynn;

	AX[n - 1] = yn;
	for (int i = n - 2; i >= 0; --i)
		AX[i] = phi[i + 1] * yn + psi[i + 1] * ynn + xi[i + 1];

	//Циклическая прогонка для линейной схемы
	if (linScheme)
		SolCircleRun(AX, rhs, afl, prea1, prec1, p, n);

}

// вычисление коэффициентов для применения предобуславливателя до начала итераций
void PreCalculateCoef(
const VM2D::Airfoil& aflP,
size_t nPanelsP,
std::vector<Point2D> & prea,
std::vector<Point2D>& prec,
std::vector<Point2D>& prea1, // для линейной схемы
std::vector<Point2D>& prec1, // для линейной схемы
bool linScheme) 
{
	Point2D p1, s1, p2, s2, di, dj;
	numvector<double, 3> alpha, lambda;
	numvector<Point2D, 3> v00, v11;
	Point2D i00, i11;

	for (size_t i = 0; i < nPanelsP; ++i)
	{
		std::array<int, 2> vecV = { ((int)i - 1 >= 0) ? int(i - 1) : int(nPanelsP - 1), (i + 1 < nPanelsP) ? int(i + 1) : 0 };
		for (int j : vecV)
		{
			p1 = aflP.getR(i + 1) - aflP.getR(j + 1);
			s1 = aflP.getR(i + 1) - aflP.getR(j);
			p2 = aflP.getR(i) - aflP.getR(j + 1);
			s2 = aflP.getR(i) - aflP.getR(j);
			di = aflP.getR(i + 1) - aflP.getR(i);
			dj = aflP.getR(j + 1) - aflP.getR(j);

			alpha = { \
				(aflP.isAfter(j, i)) ? 0.0 : VMlib::Alpha(s2, s1), \
				VMlib::Alpha(s2, p1), \
				(aflP.isAfter(i, j)) ? 0.0 : VMlib::Alpha(p1, p2) \
			};

			lambda = { \
				(aflP.isAfter(j, i)) ? 0.0 : VMlib::Lambda(s2, s1), \
				VMlib::Lambda(s2, p1), \
				(aflP.isAfter(i, j)) ? 0.0 : VMlib::Lambda(p1, p2) \
			};

			v00 = {
				VMlib::Omega(s1, di.unit(), dj.unit()),
				-VMlib::Omega(di, di.unit(), dj.unit()),
				VMlib::Omega(p2, di.unit(), dj.unit())
			};

			i00 = IDPI / di.length() * (-(alpha[0] * v00[0] + alpha[1] * v00[1] + alpha[2] * v00[2]).kcross() \
				+ (lambda[0] * v00[0] + lambda[1] * v00[1] + lambda[2] * v00[2]));

			if (aflP.isAfter(j, i))
				prec[i] = (i00.kcross()) * (1.0 / dj.length());
			if (aflP.isAfter(i, j))
				prea[i] = (i00.kcross()) * (1.0 / dj.length());

			if (linScheme) {

				v11 = { 1.0 / (12.0 * di.length() * dj.length()) * \
				(2.0 * (s1 & Omega(s1 - 3.0 * p2, di.unit(), dj.unit())) * VMlib::Omega(s1, di.unit(), dj.unit()) - \
					s1.length2() * (s1 - 3.0 * p2)) - 0.25 * VMlib::Omega(s1, di.unit(), dj.unit()),
				-di.length() / (12.0 * dj.length()) * Omega(di, dj.unit(), dj.unit()),
				-dj.length() / (12.0 * di.length()) * Omega(dj, di.unit(), di.unit()) };

				i11 = (IDPI / di.length()) * ((alpha[0] + alpha[2]) * v11[0] + (alpha[1] + alpha[2]) * v11[1] + alpha[2] * v11[2]\
					+ ((lambda[0] + lambda[2]) * v11[0] + (lambda[1] + lambda[2]) * v11[1] + lambda[2] * v11[2] \
						+ 1.0 / 12.0 * (dj.length() * di.unit() + di.length() * dj.unit() - 2.0 * VMlib::Omega(s1, di.unit(), dj.unit()))).kcross());
				if (aflP.isAfter(j, i))
					prec1[i] = (i11) * (1.0 / dj.length());
				if (aflP.isAfter(i, j))
					prea1[i] = (i11) * (1.0 / dj.length());
			}
		}
	}
}


#ifdef USE_CUDA
void GMRES(const VM2D::World2D& W,
	std::vector<std::vector<double>>& X,
	std::vector<double>& R,
	const std::vector<std::vector<double>>& rhs,
	const std::vector<double> rhsReg,
	int& niter,
	bool linScheme)
{
	const size_t nAfl = W.getNumberOfAirfoil();

	size_t totalVsize = 0;
	for (size_t i = 0; i < nAfl; ++i)
		totalVsize += W.getAirfoil(i).getNumberOfPanels();

	size_t nTotPan = 0;
	for (size_t s = 0; s < nAfl; ++s)
		nTotPan += W.getAirfoil(s).getNumberOfPanels();

	if (linScheme)
		totalVsize *= 2;

	size_t m;
	const size_t iterSize = 50;
	double beta;
	std::vector<double> w(totalVsize + nAfl), c, s;
	c.reserve(iterSize);
	s.reserve(iterSize);

	std::vector<std::vector<double>> V, H, diag(nAfl);
	V.reserve(iterSize);
	V.resize(1);
	H.resize(iterSize + 1);
	for (int i = 0; i < iterSize + 1; ++i)
		H[i].resize(iterSize);

	for (int p = 0; p < (int)nAfl; ++p)
	{
		size_t nPanelsP = W.getAirfoil(p).getNumberOfPanels();
		if (!linScheme)
			diag[p].resize(nPanelsP);
		else
			diag[p].resize(2 * nPanelsP);

		for (int i = 0; i < nPanelsP; ++i)
		{
			double lenI = W.getAirfoil(p).len[i];
			diag[p][i] = 0.5 / lenI;

			if (linScheme)
				diag[p][i + nPanelsP] = (1.0 / 24.0) / lenI;
		}
	}


	std::vector<std::vector<double>> residual(nAfl);
	for (size_t p = 0; p < nAfl; ++p)
		residual[p] = rhs[p];

	for (size_t i = 0; i < nAfl; ++i)
		V[0].insert(V[0].end(), residual[i].begin(), residual[i].end());

	for (size_t i = 0; i < nAfl; ++i)
		V[0].push_back(rhsReg[i]); // суммарная гамма

	// PRECONDITIONER
	std::vector<std::vector<Point2D>> prea(nAfl), prec(nAfl), prea1(nAfl), prec1(nAfl);

	for (size_t p = 0; p < nAfl; ++p)
	{
		const auto& aflP = W.getAirfoil(p);
		size_t nPanelsP = aflP.getNumberOfPanels();

		prea[p].resize(nPanelsP);
		prec[p].resize(nPanelsP);
		prea1[p].resize(nPanelsP);
		prec1[p].resize(nPanelsP);

		PreCalculateCoef(aflP, nPanelsP, prea[p], prec[p], prea1[p], prec1[p], linScheme);
	}

	std::vector<double> vbuf1;
	size_t np = 0;

	for (size_t p = 0; p < nAfl; ++p)
	{
		size_t nPanelsP = W.getAirfoil(p).getNumberOfPanels();

		if (!linScheme)
			vbuf1.resize(nPanelsP + 1);
		else
			vbuf1.resize(2 * nPanelsP + 1);

		if (!linScheme)
			np += ((p == 0) ? 0 : W.getAirfoil(p - 1).getNumberOfPanels());
		else
			np += ((p == 0) ? 0 : 2 * W.getAirfoil(p - 1).getNumberOfPanels());

		for (size_t i = 0; i < nPanelsP; ++i)
		{
			vbuf1[i] = V[0][np + i];
			if (linScheme)
				vbuf1[i + nPanelsP] = V[0][np + nPanelsP + i];
		}

		if (!linScheme)
			vbuf1[nPanelsP] = V[0][V[0].size() - (nAfl - p)];
		else
			vbuf1[2 * nPanelsP] = V[0][V[0].size() - (nAfl - p)];

		SolM(W, vbuf1, vbuf1, prea, prec, prea1, prec1, (int)p, (int)nPanelsP, linScheme);

		for (size_t j = np; j < np + nPanelsP; ++j)
		{
			V[0][j] = vbuf1[j - np];
			if (linScheme)
				V[0][j + nPanelsP] = vbuf1[j - np + nPanelsP];
		}

		V[0][V[0].size() - (nAfl - p)] = vbuf1[vbuf1.size() - 1];
	}

	beta = norm(V[0]);
	if (beta > 0)
		V[0] = (1.0 / beta) * V[0];
	else
	{
		for (auto& xp : X)
			for (auto& x : xp)
				x = 0.0;
		return;
	}

	double gs = beta;

	std::vector<double> bufnewSol((linScheme ? 2 : 1) * nTotPan), bufcurrentSol(totalVsize, 0.0);

	W.getNonConstCuda().AllocateSolution(W.getNonConstCuda().dev_sol, nTotPan);
	if (linScheme)
		W.getNonConstCuda().AllocateSolution(W.getNonConstCuda().dev_solLin, nTotPan);
	else
		W.getNonConstCuda().dev_solLin = nullptr;

	double*& dev_ptr_pt = W.getAirfoil(0).devRPtr;

	double*& dev_ptr_rhs = W.getAirfoil(0).devRhsPtr;
	double*& dev_ptr_rhsLin = W.getAirfoil(0).devRhsLinPtr;
	double* linPtr = (double*)(!linScheme ? nullptr : dev_ptr_rhsLin);

	double timingsToRHS[7];

	BHcu::wrapperMatrixToVector wrapper(
		(double*)dev_ptr_pt,    //начала и концы панелей
		(double*)dev_ptr_rhs,   //куда сохранить результат 
		linPtr,
		W.getNonConstCuda().CUDAptrsAirfoilVrt[0],  //указатели на дерево вихрей
		true,                   //признак перестроения дерева вихрей				
		(int)nTotPan,           //общее число панелей на всех профилях
		timingsToRHS,           //засечки времени				
		W.getPassport().numericalSchemes.fastGmresTheta,
		W.getPassport().numericalSchemes.multipoleOrderGmres,
		W.getPassport().numericalSchemes.boundaryCondition.second
	);

	//Iterations
	for (int j = 0; j < totalVsize - 1; ++j) //+ n.size() 
	{
		size_t npred = 0;
		if (!linScheme)
			for (size_t i = 0; i < W.getNumberOfAirfoil(); ++i)
			{
				for (size_t p = 0; p < W.getAirfoil(i).getNumberOfPanels(); ++p)
					bufcurrentSol[npred + p] = V[j][npred + p] / W.getAirfoil(i).len[p];
				npred += W.getAirfoil(i).getNumberOfPanels();
			}
		else
		{
			for (size_t i = 0; i < W.getNumberOfAirfoil(); ++i)
			{
				for (size_t p = 0; p < W.getAirfoil(i).getNumberOfPanels(); ++p)
				{
					bufcurrentSol[npred + p] = V[j][2 * npred + p] / W.getAirfoil(i).len[p];
					bufcurrentSol[nTotPan + npred + p] = V[j][2 * npred + W.getAirfoil(i).getNumberOfPanels() + p] / W.getAirfoil(i).len[p];
				}
				npred += W.getAirfoil(i).getNumberOfPanels();
			}
			W.getNonConstCuda().SetSolution(bufcurrentSol.data() + nTotPan, W.getNonConstCuda().dev_solLin, nTotPan);
		}
		W.getNonConstCuda().SetSolution(bufcurrentSol.data(), W.getNonConstCuda().dev_sol, nTotPan);


		wrapper.calculate(
			(double*)W.getNonConstCuda().dev_sol,
			(double*)W.getNonConstCuda().dev_solLin, j);

		if (!linScheme)
			W.getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_rhs, w.data(), 22);

		if (linScheme)
		{
			W.getCuda().CopyMemFromDev<double, 1>(nTotPan, dev_ptr_rhs, bufnewSol.data(), 22);
			W.getCuda().CopyMemFromDev<double, 1>(nTotPan, linPtr, bufnewSol.data() + nTotPan, 22);
		}

		npred = 0;
		if (linScheme)
			for (size_t i = 0; i < W.getNumberOfAirfoil(); ++i)
			{
				for (size_t j = 0; j < W.getAirfoil(i).getNumberOfPanels(); ++j)
				{
					w[2 * npred + j] = bufnewSol[npred + j];
					w[2 * npred + W.getAirfoil(i).getNumberOfPanels() + j] = bufnewSol[nTotPan + npred + j];
				}
				npred += W.getAirfoil(i).getNumberOfPanels();
			}

		size_t cntr = 0;
		for (size_t pi = 0; pi < nAfl; ++pi)
		{
			size_t nPanelsPi = W.getAirfoil(pi).getNumberOfPanels();
#pragma omp parallel for
			for (int i = 0; i < (int)nPanelsPi; ++i)
			{
				w[cntr + i] += V[j][totalVsize + pi] - V[j][cntr + i] * diag[pi][i];
				if (linScheme)
					w[cntr + i + nPanelsPi] += -V[j][cntr + i + nPanelsPi] * diag[pi][i + nPanelsPi];
			}

			cntr += nPanelsPi;
			if (linScheme)
				cntr += nPanelsPi;
		}

		for (size_t i = 0; i < nAfl; ++i)
			w[totalVsize + i] = 0.0;

		cntr = 0;
		for (size_t i = 0; i < nAfl; ++i)
		{
			size_t nPanelsI = W.getAirfoil(i).getNumberOfPanels();
			for (size_t k = 0; k < nPanelsI; ++k)
				w[totalVsize + i] += V[j][cntr + k];
			cntr += nPanelsI;
			if (linScheme)
				cntr += nPanelsI;
		}


		// PRECONDITIONER
		std::vector<double> vbuf2;
		np = 0;

		for (size_t p = 0; p < nAfl; ++p)
		{
			size_t nPanelsP = W.getAirfoil(p).getNumberOfPanels();

			if (!linScheme)
				vbuf2.resize(nPanelsP + 1);
			else
				vbuf2.resize(2 * nPanelsP + 1);

			if (!linScheme)
				np += ((p == 0) ? 0 : W.getAirfoil(p - 1).getNumberOfPanels());
			else
				np += ((p == 0) ? 0 : 2 * W.getAirfoil(p - 1).getNumberOfPanels());

			for (size_t i = 0; i < nPanelsP; ++i)
			{
				vbuf2[i] = w[np + i];
				if (linScheme)
					vbuf2[i + nPanelsP] = w[np + nPanelsP + i];
			}

			if (!linScheme)
				vbuf2[nPanelsP] = w[w.size() - (nAfl - p)];
			else
				vbuf2[2 * nPanelsP] = w[w.size() - (nAfl - p)];

			SolM(W, vbuf2, vbuf2, prea, prec, prea1, prec1, (int)p, (int)nPanelsP, linScheme);

			for (size_t j = np; j < np + nPanelsP; ++j)
			{
				w[j] = vbuf2[j - np];
				if (linScheme)
					w[j + nPanelsP] = vbuf2[j - np + nPanelsP];
			}

			w[w.size() - (nAfl - p)] = vbuf2[vbuf2.size() - 1];
		}

		if (j == H.size() - 1)
		{
			H.resize(j * 2 + 1);
			for (int i = 0; i < j * 2 + 1; ++i)
				H[i].resize(j * 2);
			V.reserve(j * 2);
			c.reserve(j * 2);
			s.reserve(j * 2);
		}

		for (int i = 0; i <= j; ++i)
		{
			double scal = 0.0;
#pragma omp parallel 
			{
#pragma omp for reduction(+:scal)
				for (int q = 0; q < (int)w.size(); ++q)
					scal += w[q] * V[i][q];

				H[i][j] = scal;
#pragma omp for				
				for (int q = 0; q < (int)w.size(); ++q)
					w[q] -= scal * V[i][q];
			}
		}
		m = j + 1;
		H[m][j] = norm(w);
		V.push_back((1 / H[m][j]) * w);

		if (IterRot(H, rhs[0], gs, c, s, (int)m, (int)totalVsize, W.getPassport().numericalSchemes.gmresEps, (int)m, false))
		{
			std::cout << "iterations in GMRES = " << j + 1 << std::endl;
			break;
		}
	} // end of iterations

	W.getNonConstCuda().ReleaseSolution(W.getNonConstCuda().dev_sol);
	if (linScheme)
		W.getNonConstCuda().ReleaseSolution(W.getNonConstCuda().dev_solLin);

	niter = (int)m;

	std::vector<double> g(m + 1);
	g[0] = beta;

	//GivensRotations
	double oldValue;
	for (int i = 0; i < m; i++)
	{
		oldValue = g[i];
		g[i] = c[i] * oldValue;
		g[i + 1] = -s[i] * oldValue;
	}
	//end of GivensRotations

	// Solve HY=g
	std::vector<double> Y(m);
	Y[m - 1] = g[m - 1] / H[m - 1][m - 1];

	double sum;
	for (int k = (int)m - 2; k >= 0; --k)
	{
		sum = 0.0;
		for (int s = k + 1; s < m; ++s)
			sum += H[k][s] * Y[s];
		Y[k] = (g[k] - sum) / H[k][k];
	}
	// end of Solve HY=g

	size_t cntr = 0;
	for (size_t p = 0; p < nAfl; p++) {

		size_t nPanelsP = W.getAirfoil(p).getNumberOfPanels();

		if (!linScheme)
			for (size_t i = 0; i < nPanelsP; i++)
			{
				sum = 0.0;
				for (size_t j = 0; j < m; j++)
					sum += V[j][i + cntr] * Y[j];
				X[p][i] += sum;
			}
		else {
			for (size_t i = 0; i < 2 * nPanelsP; i++)
			{
				sum = 0.0;
				for (size_t j = 0; j < m; j++)
					sum += V[j][i + cntr] * Y[j];
				X[p][i] += sum;
			}
		}
		cntr += nPanelsP;

		if (linScheme)
			cntr += nPanelsP;

	}
	sum = 0.0;
	for (size_t p = 0; p < nAfl; p++)
	{
		for (size_t j = 0; j < m; j++)
			sum += V[j][totalVsize + p] * Y[j];
		R[p] += sum;
	}
}

#endif