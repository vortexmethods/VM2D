#pragma once
#include <vector>
#include <complex>
#include "utils.h"
#include "simple_math.h"
#include <unordered_map>
//#include "cuda_utils.h"

namespace fmm {

// associated Legendre(n,m), n=0..N-1, m=0..n
__HOST__ __DEVICE__ inline void ComputePnm(int N, double theta, double* Pnm)
{
#ifdef __NVCC__
	double x, y;
	sincos(theta, &y, &x);
#else
	auto x = cos(theta);
	auto y = sin(theta);
#endif
	Pnm[0] = 1;
	for (int n = 1; n < N; ++n)
	{
		for (int m = 0; m <= n - 2; ++m)
		{
			Pnm[n * (n + 1) / 2 + m] = ((2 * n - 1) * x * Pnm[(n - 1) * n / 2 + m] - (n + m - 1) * Pnm[(n - 2) * (n - 1) / 2 + m]) / (n - m);
		}
		Pnm[n * (n + 1) / 2 + n - 1] = Pnm[(n - 1) * n / 2 + n - 1] * x * (2 * n - 1);
		Pnm[n * (n + 1) / 2 + n] = Pnm[(n - 1) * n / 2 + n - 1] * y * (1 - 2 * n);
	}
}
// Knm = sqrt{(n - |m|)! / (n + |m|)!} - normalization coefficient for Spherical harmonic Ynm
// computes only m >= 0 
[[deprecated("Use only when constexpr Knm unavailable")]]
treevector<double> ComputeKnm(int N);

// Anm = (-1)^n / sqrt{(n-m)!(n+m)!}
[[deprecated("Use only when constexpr Anm unavailable")]]
treevector<double> ComputeAnm(int N);

// WignerD(n,m,k), n=0..N-1, m=0..n, k=-n..n
std::vector<double> ComputeWignerD(int N, double theta);

namespace detail {

	inline double constexpr sqrtNewton(double x, double curr, double prev)
	{
		return (prev - curr) / curr < 1.e-16 && (prev - curr) / curr > -1.e-16 ? curr : sqrtNewton(x, 0.5 * (curr + x / curr), curr);
	}

	inline double constexpr constexpr_sqrt(double x)
	{
		return sqrtNewton(x, x, 0);
	}

	template <int N>
	struct knm_wrapper
	{
		std::array<double, N * (N + 1) / 2> Knm;
		FMM_CONSTEXPR knm_wrapper() : Knm()
		{
			Knm[0] = 1.0;
			for (int n = 1; n < N; ++n)
			{
				for (int m = 0; m <= n - 1; ++m)
					Knm[n * (n + 1) / 2 + m] = Knm[n * (n - 1) / 2 + m] * (n - m) / (n + m);
				Knm[n * (n + 1) / 2 + n] = Knm[n * (n + 1) / 2 + n - 1] / (2 * n);
			}
			for (int n = 1; n < N; ++n)
			{
				for (int m = 0; m <= n; ++m)
#ifdef FMM_CONSTEXPR_MATH
					Knm[n * (n + 1) / 2 + m] = constexpr_sqrt(Knm[n * (n + 1) / 2 + m]);
#else
					Knm[n * (n + 1) / 2 + m] = sqrt(Knm[n * (n + 1) / 2 + m]);
#endif
			}
		}
		FMM_CONSTEXPR double operator()(int n, int m) const
		{
			return Knm[n * (n + 1) / 2 + m];
		}
	};

	inline FMM_CONSTEXPR knm_wrapper<_3d_MAX_MULTIPOLE_NUM> Knm;
	
	template <int N>
	struct anm_wrapper
	{
		std::array<double, 2 * N * (2 * N + 1) / 2> Anm;
		FMM_CONSTEXPR anm_wrapper() : Anm()
		{
			Anm[0] = 1.0;
			for (int n = 1; n < 2 * N; ++n)
			{
				for (int m = 0; m <= n - 1; ++m)
					Anm[n * (n + 1) / 2 + m] = Anm[(n - 1) * n / 2 + m] * (n - m) * (n + m);
				Anm[n * (n + 1) / 2 + n] = Anm[n * (n + 1) / 2 + n - 1] * 2 * n;
			}
			for (int n = 1; n < 2 * N; ++n)
			{
				for (int m = 0; m <= n; ++m)
#ifdef FMM_CONSTEXPR_MATH
					Anm[n * (n + 1) / 2 + m] = ni(n) / constexpr_sqrt(Anm[n * (n + 1) / 2 + m]);
#else
					Anm[n * (n + 1) / 2 + m] = ni(n) / sqrt(Anm[n * (n + 1) / 2 + m]);
#endif
			}
		}
		FMM_CONSTEXPR double operator()(int n, int m) const
		{
			return Anm[n * (n + 1) / 2 + m];
		}
	};

	inline FMM_CONSTEXPR anm_wrapper<_3d_MAX_MULTIPOLE_NUM> Anm;

	template <int N>
	struct m2lcoef_wrapper
	{
		std::array<double, N * (N + 1) * (N + 2) / 2 + 1> m2lcoef;
		FMM_CONSTEXPR m2lcoef_wrapper() : m2lcoef()
		{
			for (int j = 0; j < N; ++j)
			{
				for (int k = 0; k <= j; ++k)
				{
					for (int n = k; n < N; ++n)
					{
						int idx = N * (j + 2) * (j + 1) / 2 + N * (k + 1) + n + 1;
						m2lcoef[idx] = Anm(n, k) * Anm(j, k) * ni(n + k) / Anm(j + n, 0);
					}
				}
			}
		}
	};
	
	inline FMM_CONSTEXPR m2lcoef_wrapper<_3d_MAX_MULTIPOLE_NUM> m2lcoef;

#ifndef FMM_CONSTEXPR_MATH
	inline double* dev_Knm; // deprecated, substituted with __constant__ version
	inline double* dev_Anm; // deprecated, substituted with __constant__ version
	inline double* dev_m2lcoef; // deprecated, substituted with __constant__ version
#endif
	inline std::unordered_map<int, std::vector<double>> dmatrix;
	inline std::unordered_map<int, std::vector<std::complex<double>>> rotation_exponents;
	inline double* dev_dmatrix;
	// if dev_dm_map[i] = key, use dev_dmatrix[i] for given key
	// dmatrices on device are sorted by key,
	// so it is possible to use binary search for given key to find i such that dev_dm_map[i] = key
	inline int* dev_dm_map; 


	void InitMathConstants(int N);
	void cudaCopyMathConstants();
	void cudaClearMathConstants();
}

} // fmm