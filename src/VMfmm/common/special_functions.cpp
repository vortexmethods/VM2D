#include "special_functions.h"
#include "simple_math.h"
#include <cmath>
#include <utility>
#include <oneapi/tbb/parallel_for.h>

namespace fmm {

treevector<double> ComputeKnm(int N)
{
    treevector<double> Knm(N);
	Knm(0, 0) = 1.0;
	for (int n = 1; n < N; ++n)
	{
		for (int m = 0; m <= n - 1; ++m)
			Knm(n, m) = Knm(n - 1, m) * (n - m) / (n + m);
		Knm(n, n) = Knm(n, n - 1) / (2 * n);
	}
	for (int n = 1; n < N; ++n)
	{
		for (int m = 0; m <= n; ++m)
			Knm(n, m) = sqrt(Knm(n, m));
	}
    return Knm;
}

treevector<double> ComputeAnm(int N)
{
    treevector<double> Anm(2 * N);
    Anm(0, 0) = 1.0;
    for (int n = 1; n < 2 * N; ++n)
    {
        for (int m = 0; m <= n - 1; ++m)
            Anm(n, m) = Anm(n - 1, m) * (n - m) * (n + m);
        Anm(n, n) = Anm(n, n - 1) * 2 * n;
    }
    for (int n = 1; n < 2 * N; ++n)
    {
        for (int m = 0; m <= n; ++m)
            Anm(n, m) = ni(n) / sqrt(Anm(n, m));
    }
    return Anm;
}

//std::pair<std::vector<std::vector<double>>,std::vector<std::vector<double>>> compute_wignerd_coefs(int N)
//{
//    std::vector<std::vector<double>> res1(N);
//    std::vector<std::vector<double>> res2(N);
//    int idx1, idx2;
//    for (int n = 0; n < N; ++n)
//	{
//        res1[n].resize((1 + 2 * n) * (1 + n));
//        res2[n].resize((1 + 2 * n) * (1 + n));
//        idx1 = 2 * n + 1;
//        auto& r1 = res1[n];
//        auto& r2 = res2[n];
//        for (int m = 0; m <= n - 1; ++m)
//		{
//			idx2 = idx1 * (m + 1) + n;
//            r2[idx2 - n] = - 1.0 * (m - n) / sqrt(n * (n + 1) - m * (m + 1));
//            for (int k = -n + 1; k <= n; ++k)
//            {
//                r1[idx2 + k] = 1.0 * sqrt((n * (n + 1.0) - k * (k - 1)) / (n * (n + 1) - m * (m + 1)));
//                r2[idx2 + k] = - 1.0 * (m + k) / sqrt(n * (n + 1) - m * (m + 1));
//            } 
//        }
//    }
//    return {res1, res2};
//}

//const auto wignerdcoefs = compute_wignerd_coefs(100);

//std::vector<double> ComputeWignerD(int N, double theta)
//{
//    assert(std::abs(theta) <= std::numbers::pi);
//    // length can be computed as 
//    // \sum_{n=0}^{N-1} \sum_{m=0}^n \sum_{k=-n}^n (1) = N * (1 + N) * (4 * N - 1) / 6
//    std::vector<double> D(N * (1 + N) * (4 * N - 1) / 6);
//
//    // theta must be in [0;pi/2], other cases are considered using properties of d-matrix
//    // see properties here: https://en.wikipedia.org/wiki/Wigner_D-matrix
//    bool isnegative = theta < 0;
//    theta = std::abs(theta);
//
//    bool closetopi = false;
//    if (theta > std::numbers::pi / 2) 
//    {
//        theta = std::numbers::pi - theta;
//        closetopi = true;
//    }
//
//    if (fabs(theta) < 1e-3)
//    {    
//        int idx0, idx1, idx2;
//        for (int n = 0; n < N; ++n)
//        {
//            idx0 = n * (1 + n) * (4 * n - 1) / 6;
//            idx1 = 2 * n + 1;
//            for (int m = 0; m <= n; ++m)
//            {
//                idx2 = idx0 + idx1 * m + n;
//                D[idx2 + m] = 1.0;
//            }
//        }
//    }
//    else
//    {
//        const auto& [wignerdcoef_1,wignerdcoef_2] = wignerdcoefs;
//        auto Pnm = ComputePnm(N, theta);
//        
//        int idx0, idx1, idx2, idx_pred;
//        for (int n = 0; n < N; ++n)
//        {
//            idx1 = n + n * (1 + n) * (4 * n - 1) / 6;
//            for (int k = -n; k < 0; ++k)
//            {
//                D[idx1 + k] = ni(k) * Pnm(n, std::abs(k)) * detail::Knm(n,std::abs(k));
//            }
//            for (int k = 0; k <= n; ++k)
//            {
//                D[idx1 + k] = Pnm(n, std::abs(k)) * detail::Knm(n,std::abs(k));
//            }
//        }
//
//        double x = sin(theta) / (1.0 + cos(theta));
//        for (int n = 0; n < N; ++n)
//        {
//            idx0 = n * (1 + n) * (4 * n - 1) / 6;
//            idx1 = 2 * n + 1;
//            const auto& c1 = wignerdcoef_1[n];
//            const auto& c2 = wignerdcoef_2[n];
//            for (int m = 0; m <= n - 1; ++m)
//            {
//                int idxc = idx1 * (m + 1) + n;
//                idx2 = idx0 + idxc;
//                idx_pred = idx0 + idx1 * m + n;
//
//                D[idx2 - n] = x * c2[idxc - n] * D[idx_pred - n];
//                for (int k = -n + 1; k <= n; ++k)
//                {
//                    D[idx2 + k] = c1[idxc + k] * D[idx_pred + k - 1] + c2[idxc + k] * x * D[idx_pred + k];
//                }
//            
//            }
//        }
//    }
//
//    int idx0, idx1, idx2;
//    auto Dcopy = D;
//    if (closetopi) // D^n_{m,k}(pi - theta) = (-1)^(n+m) * D^n_{m,-k}(theta)
//    {
//        for (int n = 1; n < N; ++n)
//        {
//            idx0 = n * (1 + n) * (4 * n - 1) / 6;
//            idx1 = 2 * n + 1;
//            for (int m = 0; m <= n; ++m)
//            {
//                idx2 = idx0 + idx1 * m + n;
//                for (int k = -n; k <= n; ++k)
//                {
//                    D[idx2 + k] = ni(n + m) * Dcopy[idx2 - k];
//                }
//            }
//        }
//    }
//
//    for (int n = 1; n < N; ++n)
//    {
//        idx0 = n * (1 + n) * (4 * n - 1) / 6;
//        idx1 = 2 * n + 1;
//        for (int m = 0; m <= n; ++m)
//        {
//            idx2 = idx0 + idx1 * m + n;
//            for (int k = -n; k <= n; ++k)
//            {
//                // define eps(m) = (-1)^(m-|m|)/2, 1/eps(m) = eps(m),
//                // because Ynm in Laplace FMM is different from Ynm in physics (denote it Ynm_p)
//                // and they are connected as Ynm_p = eps(m) * sqrt((2n+1)/4pi) * Ynm,
//                // Wigner D-matrix are introduced for Ynm_p such that Ynm_p = sum_k Ynk_p*Dnkm,
//                // or equivalently  eps(m)*Ynm = sum_k eps(k)*Ynk*Dnkm <=> Ynm = sum_k Ynk*[eps(m)*eps(k)*Dnkm]
//                // so Wigner D-matrix in fmm should have multiplier eps(m)*eps(k)
//                // since we consider only m>=0, we have only eps(k) left (eps(m) = 1 for m >= 0)
//                D[idx2 + k] *= ni((k - std::abs(k)) / 2); 
//
//                if (isnegative) // D^n_{m,k}(-theta) = (-1)^(m - k) * D^n_{m,k}(theta)
//                    D[idx2 + k] *= ni(m - k);
//            }
//        }
//    }
//
//    return D;
//}

std::vector<double> compute_wignerd_coef(int N)
{
    std::vector<double> g(N * (N + 1) / 2);
    g[0] = 1;
    for (int n = 1; n < N; ++n)
    {
        g[n * (n + 1) / 2] = sqrt((2 * n - 1.) / (2 * n)) * g[n * (n - 1) / 2];
        for (int m = 1; m <= n; ++m)
        {
            g[n * (n + 1) / 2 + m] = sqrt((n - m + 1.) / (n + m)) * g[n * (n + 1) / 2 + m - 1];
        }
    }
    return g;
}

const auto wignerd_auxillary_coef = compute_wignerd_coef(100);

std::vector<double> ComputeWignerD(int N, double theta)
{
    assert(std::abs(theta) <= std::numbers::pi);
    // size can be computed as 
    // \sum_{n=0}^{N-1} \sum_{m=0}^n \sum_{k=-n}^n (1) = N * (1 + N) * (4 * N - 1) / 6
    // linear index (n,m,k) can be found as 
    // \sum_{n'=0}^{n-1} \sum_{m=0}^{n'} \sum_{k=-n'}^{n'} 1 +
    // + \sum_{m'=0}^{m-1} \sum_{k=-n}^n 1 +
    // + \sum_{k'=-n}^{k-1} 1 =
    // = 1/6 n (1 + n) (-1 + 4 n) + m (1 + 2 n) + (k + n)
    std::vector<double> D(N * (1 + N) * (4 * N - 1) / 6);

    // theta must be in [0;pi/2], other cases are considered using properties of d-matrix
    // see properties here: https://en.wikipedia.org/wiki/Wigner_D-matrix
    bool isnegative = theta < 0;
    theta = std::abs(theta);

    bool closetopi = false;
    if (theta > std::numbers::pi / 2)
    {
        theta = std::numbers::pi - theta;
        closetopi = true;
    }

    if (fabs(theta) < 1e-3) // theta = 0 => d^n_mm = 1, other zero
    {
        int idx0, idx1, idx2;
        for (int n = 0; n < N; ++n)
        {
            idx0 = n * (1 + n) * (4 * n - 1) / 6 + n;
            idx1 = 2 * n + 1;
            for (int m = 0; m <= n; ++m)
            {
                idx2 = idx0 + idx1 * m;
                D[idx2 + m] = 1.0;
            }
        }
    }
    else
    {
        const auto& g = wignerd_auxillary_coef;
        int idx0, idx1, idx2, idx3;
        double s = sin(theta);
        double c = 1.0 + cos(theta);
        double f = s / c;
        std::vector<double> sn(2 * N), cn(2 * N);
        sn[N] = cn[N] = 1.0;
        for (int i = 1; i < N; ++i)
        {
            sn[N + i] = sn[N + i - 1] * s;
            sn[N - i] = sn[N - i + 1] / s;
            cn[N + i] = cn[N + i - 1] * c;
            cn[N - i] = cn[N - i + 1] / c;
        }
        

        // 1) compute d^n_mn
        for (int n = 0; n < N; ++n)
        {
            idx0 = n * (1 + n) * (4 * n - 1) / 6 + 2 * n;
            idx1 = 1 + 2 * n;
            for (int m = 0; m <= n; ++m)
            {
                D[idx0 + m * idx1] = ni(n+m) * g[n * (n + 1) / 2 + m] * cn[N + m] * sn[N + n - m];
            }
        }

        // 2) compute d^n_n,m-1
        for (int n = 0; n < N; ++n)
        {
            idx0 = n * (1 + n) * (4 * n - 1) / 6 + n * (2 * n + 1) + n;
            for (int m = n; m > -n; --m)
            {
                D[idx0 + m - 1] = (n + m) / sqrt(n * (n + 1.) - m * (m - 1.)) * f * D[idx0 + m];
            }
        }

        // 3) compute d^n_m,k-1
        for (int n = 0; n < N; ++n)
        {
            idx0 = n * (1 + n) * (4 * n - 1) / 6 + n;
            idx1 = 2 * n + 1;
            for (int m = n - 1; m >= 0; --m)
            {
                idx2 = m * idx1;
                idx3 = (m + 1) * idx1;
                for (int k = n; k > -n; --k)
                {
                    D[idx0 + idx2 + k - 1] = sqrt((n * (n + 1.) - m * (m + 1.)) / (n * (n + 1.) - k * (k - 1.))) * D[idx0 + idx3 + k] +
                        (m + k) / sqrt(n * (n + 1.) - k * (k - 1.)) * f * D[idx0 + idx2 + k];
                }
            }
        }

    }

    int idx0, idx1, idx2;
    auto Dcopy = D;
    if (closetopi) // D^n_{m,k}(pi - theta) = (-1)^(n+m) * D^n_{m,-k}(theta)
    {
        for (int n = 1; n < N; ++n)
        {
            idx0 = n * (1 + n) * (4 * n - 1) / 6;
            idx1 = 2 * n + 1;
            for (int m = 0; m <= n; ++m)
            {
                idx2 = idx0 + idx1 * m + n;
                for (int k = -n; k <= n; ++k)
                {
                    D[idx2 + k] = ni(n + m) * Dcopy[idx2 - k];
                }
            }
        }
    }

    for (int n = 1; n < N; ++n)
    {
        idx0 = n * (1 + n) * (4 * n - 1) / 6;
        idx1 = 2 * n + 1;
        for (int m = 0; m <= n; ++m)
        {
            idx2 = idx0 + idx1 * m + n;
            for (int k = -n; k <= n; ++k)
            {
                // define eps(m) = (-1)^(m-|m|)/2, 1/eps(m) = eps(m),
                // because Ynm in Laplace FMM is different from Ynm in physics (denote it Ynm_p)
                // and they are connected as Ynm_p = eps(m) * sqrt((2n+1)/4pi) * Ynm,
                // Wigner D-matrix are introduced for Ynm_p such that Ynm_p = sum_k Ynk_p*Dnkm,
                // or equivalently  eps(m)*Ynm = sum_k eps(k)*Ynk*Dnkm <=> Ynm = sum_k Ynk*[eps(m)*eps(k)*Dnkm]
                // so Wigner D-matrix in fmm should have multiplier eps(m)*eps(k)
                // since we consider only m>=0, we have only eps(k) left (eps(m) = 1 for m >= 0)
                D[idx2 + k] *= ni((k - std::abs(k)) / 2);

                if (isnegative) // D^n_{m,k}(-theta) = (-1)^(m - k) * D^n_{m,k}(theta)
                    D[idx2 + k] *= ni(m - k);
            }
        }
    }

    return D;
}

std::vector<std::complex<double>> ComputeExp(int N, double phi)
{
    std::vector<std::complex<double>> res(N);
    for (int i = 0; i < N; ++i)
        res[i] = std::complex<double>(cos(i * phi), sin(i * phi));
    return res;
}

namespace detail {

void InitMathConstants(int N)
{
    for (int i = -3; i < 4; ++i)
        for (int j = -3; j < 4; ++j)
            for (int k = -3; k < 4; ++k) {
                if (std::abs(i) <= 1 && std::abs(j) <= 1 && std::abs(k) <= 1)
                    continue;
                auto rho = DecToSph(Vector3d(i, j, k));
                auto theta = rho[1];
                auto phi = rho[2];
                auto key = doublehash(theta);
                dmatrix[key] = ComputeWignerD(N, theta);
                dmatrix[-key] = ComputeWignerD(N, -theta);
                key = doublehash(phi);
                rotation_exponents[key] = ComputeExp(N, phi);
                rotation_exponents[-key] = ComputeExp(N, -phi);
            }
    //std::cout << dmatrix.size() << std::endl;
}

} // detail

} // fmm