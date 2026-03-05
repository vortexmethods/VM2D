#pragma once
#ifndef __NVCC__
#include <oneapi/tbb/parallel_for.h>
#endif
#include <vector>
#include <complex>
#include <numbers>
#include <functional>
#include <fstream>
//#include "cuda_utils.h"
#include "mpi_utils.h"
#include <iostream>



namespace fmm {

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

namespace detail {

template<int n>
struct BinomNewton_wrapper
{
	std::array<double, n * n> cft;
	int dim;

	constexpr BinomNewton_wrapper() : cft(), dim(n)
	{
		cft[0 * dim + 0] = 1.0;

		for (int i = 1; i < n; ++i)
		{
			cft[i * dim + 0] = 1.0;
			cft[i * dim + i] = 1.0;
			for (int j = 1; j < i; ++j)
				cft[i * dim + j] = cft[(i - 1) * dim + j] + cft[(i - 1) * dim + (j - 1)];
		}
	}
	constexpr double operator()(int p, int q) const
	{
		return cft[p * dim + q];
	}
};

constexpr BinomNewton_wrapper<2 * _2d_MAX_MULTIPOLE_NUM> binom;

}

class BinomNewton
{
public:
	std::vector<double> cft;
	int dim;

public:
	BinomNewton(int n)
		: dim(n + 1)
	{
		cft.resize(dim * dim, 0.0);
		cft[0 * dim + 0] = 1.0;

		for (int i = 1; i <= n; ++i)
		{
			cft[i * dim + 0] = 1.0;
			cft[i * dim + i] = 1.0;
			for (int j = 1; j < i; ++j)
				cft[i * dim + j] = cft[(i - 1) * dim + j] + cft[(i - 1) * dim + (j - 1)];
		}
	}

	double operator()(int p, int q) const
	{
		return cft[p * dim + q];
	}
};

template <typename T>
__DEVICE__ __HOST__ constexpr T MyPow(T base, unsigned int exp)
{
	T res = 1;
	while (exp) {
		if (exp & 1)
		{
			res *= base;
		}
		exp >>= 1;
		base *= base;
	}
	return res;
}

__DEVICE__ __HOST__ constexpr inline double ni(int n)
{
	return ((n & 1) == 1) ? -1.0 : 1.0;
}

__DEVICE__ __HOST__ inline Vector3d DecToSph(const Vector3d& vec)
{
	const double eps = 1e-12;
	Vector3d c;
	c[0] = abs(vec) + eps;
	c[1] = acos(vec[2] / c[0]);
	if (fabs(vec[0]) + fabs(vec[1]) < eps) {
		c[2] = 0;
	}
	else if (fabs(vec[0]) < eps) {
		c[2] = vec[1] / fabs(vec[1]) * /*std::numbers::pi*/ 3.14159265358979323846 * 0.5;
	}
	else {
		c[2] = atan2(vec[1], vec[0]);
	}
	return c;
}

inline double Potential2d(const particle2d& p1, const particle2d& p2)
{
	return p2.q * 0.5 * log(std::max(norm(p1.center - p2.center), FORCE_EPS2));
}

inline std::complex<double> Force2d(const particle2d& p1, const particle2d& p2)
{
	auto dz = p1.center - p2.center;
	return p2.q * std::conj(dz) / std::max(norm(dz), FORCE_EPS2);
}

inline double Potential3d(const particle3d& p1, const particle3d& p2)
{
	return p2.q / std::max(abs(p1.center - p2.center), FORCE_EPS);
}

inline Vector3d Force3d(const particle3d& p1, const particle3d& p2)
{
	auto dr = p1.center - p2.center;
	return p2.q * dr / MyPow(std::max(abs(dr), FORCE_EPS), 3);
}

inline Vector3d Force3d(const particle3d3& p1, const particle3d3& p2)
{
	auto dr = p1.center - p2.center;
	return cross(p2.q, dr) / MyPow(std::max(abs(dr), FORCE_EPS), 3);
}

// ── Kernel interpolation tables ──────────────────────────────────────────────
// lambdaEta / lambdaEtaA are called for every near-field pair (argLE < 3.0).
// Replacing std::erf + std::exp with a linear-interpolation table gives a
// ~10–20× speedup for the leaf direct computation.
//
// Table range [0, XI_MAX_INTERP] with TABLE_SIZE points; step = XI_MAX_INTERP/(N-1).
// For argLE outside the table the original analytical formula is used as fallback.

namespace detail_interp {

static constexpr int    TABLE_SIZE    = 2049;
static constexpr double XI_MAX_INTERP = 5.0;
static constexpr double INV_STEP      = (TABLE_SIZE - 1) / XI_MAX_INTERP;

struct KernelTables {
    double lambda[TABLE_SIZE];
    double eta[TABLE_SIZE];
    double lambdaA[TABLE_SIZE];
    double etaA[TABLE_SIZE];

    KernelTables() {
        constexpr double sqrtPi = 1.7724538509055159; // sqrt(π)
        for (int i = 0; i < TABLE_SIZE; ++i) {
            double xi = i / INV_STEP;
            if (xi < 1e-10) {
                lambda[i]  = 8.0 / (3.0 * sqrtPi);
                eta[i]     = 8.0 / (5.0 * sqrtPi);
                lambdaA[i] = 0.0;
                etaA[i]    = 0.0;
            } else {
                double xi2 = xi * xi, xi3 = xi2 * xi, xi4 = xi2 * xi2;
                double xi5 = xi3 * xi2, xi6 = xi4 * xi2;
                double e = std::erf(xi), ex = std::exp(-xi2);
                lambda[i]  = -e/xi3 + (2/sqrtPi)*(2 + 1/xi2)*ex;
                eta[i]     =  3*e/xi5 - (2/sqrtPi/xi2)*(2 + 3/xi2)*ex;
                lambdaA[i] =  3*e/xi4 - (8/sqrtPi)*(3/(4*xi3) + 1/(2*xi) + xi)*ex;
                etaA[i]    = -15*e/xi6 + (1/sqrtPi)*(8/xi + 20/xi3 + 30/xi5)*ex;
            }
        }
    }
};

inline const KernelTables& tables() {
    static const KernelTables t;
    return t;
}

inline std::pair<double,double> lookup2(const double* tA, const double* tB, double xi) {
    double idx_f = xi * INV_STEP;
    int    idx   = static_cast<int>(idx_f);
    if (idx >= TABLE_SIZE - 1) idx = TABLE_SIZE - 2;
    double frac  = idx_f - idx;
    double a = tA[idx] * (1.0 - frac) + tA[idx + 1] * frac;
    double b = tB[idx] * (1.0 - frac) + tB[idx + 1] * frac;
    return {a, b};
}

} // namespace detail_interp

// Generic template — used for non-double types (e.g. future SIMD)
template <typename T>
std::pair<T, T> lambdaEta(T xi)
{
	T sqrtPi = sqrt((T)3.141592653589793);
	if (fabs(xi) < 1e-10)
		return { 8.0 / (3.0 * sqrtPi), 8.0 / (5.0 * sqrtPi) };
	T xi2 = xi * xi, xi3 = xi2 * xi, xi5 = xi3 * xi2;
	T e = std::erf(xi), ex = std::exp(-xi2);
	return {
		   -e/xi3 + (2/sqrtPi)      *(2 + 1/xi2)*ex,
		3 * e/xi5 - (2/sqrtPi/xi2)*(2 + 3/xi2)*ex
	};
};

// double specialisation — uses interpolation tables (no erf/exp)
template <>
inline std::pair<double, double> lambdaEta<double>(double xi)
{
    if (xi >= detail_interp::XI_MAX_INTERP) {
        // far-field fallback (power-law, no erf/exp needed)
        double xi3 = xi*xi*xi, xi5 = xi3*xi*xi;
        return { -1.0/xi3, 3.0/xi5 };
    }
    const auto& t = detail_interp::tables();
    return detail_interp::lookup2(t.lambda, t.eta, xi);
}

template <typename T>
std::pair<T, T> lambdaEtaA(T xi)
{
	T sqrtPi = sqrt((T)3.141592653589793);
	if (fabs(xi) < 1e-10) return { 0.0, 0.0 };
	T xi2 = xi*xi, xi3 = xi2*xi, xi4 = xi2*xi2, xi5 = xi3*xi2, xi6 = xi4*xi2;
	T e = std::erf(xi), ex = std::exp(-xi2);
	return {
		 3*e/xi4 - (8/sqrtPi)*(3/(4*xi3) + 1/(2*xi) + xi)*ex,
		-15*e/xi6 + (1/sqrtPi)*(8/xi + 20/xi3 + 30/xi5)*ex
	};
};

// double specialisation — uses interpolation tables (no erf/exp)
template <>
inline std::pair<double, double> lambdaEtaA<double>(double xi)
{
    if (xi >= detail_interp::XI_MAX_INTERP) {
        double xi4 = xi*xi*xi*xi, xi6 = xi4*xi*xi;
        return { 3.0/xi4, -15.0/xi6 };
    }
    const auto& t = detail_interp::tables();
    return detail_interp::lookup2(t.lambdaA, t.etaA, xi);
}


inline Vector3d velDipole3d(const particle3d3& p1, const particle3d3& p2)
{
	auto dr = (p1.center - p2.center);
	auto Ldr = abs(dr);

	auto argLE = (Ldr / FORCE_EPS);
	if (argLE < 3.0)
	{
		auto [lambda, eta] = lambdaEta(argLE);
		auto result = p2.q * lambda / FORCE_EPS3 + dot(p2.q, dr) * dr * eta / FORCE_EPS5;
		//if (std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2]))
		//	std::cout << "AAA" << std::endl;
		return result;
	}
	else
	{
		auto Ldr2 = Ldr * Ldr;
		auto Ldr3 = Ldr2 * Ldr;
		auto Ldr5 = Ldr3 * Ldr2;
		auto result = p2.q / (-Ldr3) + 3 * dot(p2.q, dr) * dr / Ldr5;
		//if (std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2]))
		//	std::cout << "BBB" << std::endl;
		return result;
	}
}

inline Vector3d momDipole3d(const particle3d3& p1, const particle3d3& p2)
{
	//p1 - �����������, p2 - ��������
	auto dr = (p1.center - p2.center);
	auto Ldr = abs(dr);

	auto argLE = (Ldr / FORCE_EPS);
	if (argLE < 3.0)
	{
		auto [lambda, eta] = lambdaEta(argLE);
		auto [lambdaA, etaA] = lambdaEtaA(argLE);

		auto result = -1 * ((dot(p1.q, p2.q) * lambdaA / FORCE_EPS4 + dot(p2.q, dr) * dot(p1.q, dr) * etaA / FORCE_EPS6) / std::max(Ldr, 1e-10) * dr + \
			(p2.q * dot(p1.q, dr) + p1.q * dot(p2.q, dr)) * eta / FORCE_EPS5);

		//if (std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2]))
		//	std::cout << "AAA" << std::endl;
		return result;
	}
	else
	{
		auto Ldr2 = Ldr * Ldr;
		auto Ldr3 = Ldr2 * Ldr;
		auto Ldr4 = Ldr2 * Ldr2;
		auto Ldr5 = Ldr3 * Ldr2;
		auto Ldr6 = Ldr3 * Ldr3;
		auto result = -1 * ((dot(p1.q, p2.q) * (3 / Ldr4) + dot(p2.q, dr) * dot(p1.q, dr) * (-15 / Ldr6)) / std::max(Ldr, 1e-10) * dr + \
			(p2.q * dot(p1.q, dr) + p1.q * dot(p2.q, dr)) * 3 / Ldr5);

		//if (std::isnan(result[0]) || std::isnan(result[1]) || std::isnan(result[2]))
		//	std::cout << "BBB" << std::endl;
		return result;
	}
}


inline void Potential3dMutual(const particle3d& p1, const particle3d& p2, double& potential1, double& potential2)
{
	auto invdr = 1.0 / std::max(abs(p1.center - p2.center), FORCE_EPS);
	potential1 += p2.q * invdr;
	potential2 += p1.q * invdr;
}

inline void Force3dMutual(const particle3d& p1, const particle3d& p2, Vector3d& force1, Vector3d& force2)
{
	auto dr = p1.center - p2.center;
	auto invdr = dr / MyPow(std::max(abs(dr), FORCE_EPS), 3);
	force1 += p2.q * invdr;
	force2 -= p1.q * invdr;
}

inline void Force3dMutual(const particle3d3& p1, const particle3d3& p2, Vector3d& force1, Vector3d& force2)
{
	auto dr = p1.center - p2.center;
	auto invdr = dr / MyPow(std::max(abs(dr), FORCE_EPS), 3);
	force1 += cross(p2.q, invdr);
	force2 -= cross(p1.q, invdr);
}

inline void Potential2dMutual(const particle2d& p1, const particle2d& p2, double& potential1, double& potential2)
{
	double dz = 0.5 * log(std::max(norm(p1.center - p2.center), FORCE_EPS2));
	potential1 += p2.q * dz;
	potential2 += p1.q * dz;
}

inline void Force2dMutual(const particle2d& p1, const particle2d& p2, std::complex<double>& force1, std::complex<double>& force2)
{
	auto dz = p1.center - p2.center;
	auto invdz = std::conj(dz) / std::max(norm(dz), FORCE_EPS2);
	force1 += p2.q * invdz;
	force2 -= p1.q * invdz;
}

#ifdef __NVCC__
__device__ inline double Potential2d(const gpu::particle2d& p1, const gpu::particle2d& p2)
{
	return p2.q * 0.5 * log(max(cuda::std::norm(p1.center - p2.center), CUDA_FORCE_EPS2));
}

__device__ inline gpu::cuda_complex Force2d(const gpu::particle2d& p1, const gpu::particle2d& p2)
{
	auto dz = p1.center - p2.center;
	return p2.q * cuda::std::conj(dz) / max(cuda::std::norm(dz), CUDA_FORCE_EPS2);
}

__device__ inline double Potential3d(const gpu::particle3d& p1, const gpu::particle3d& p2)
{
	double dr2 = max(norm(p1.center - p2.center), CUDA_FORCE_EPS2);
	return p2.q * rsqrt(dr2);
}

__device__ inline Vector3d Force3d(const gpu::particle3d& p1, const gpu::particle3d& p2)
{
	auto dr = p1.center - p2.center;
	auto dr2 = max(norm(dr), CUDA_FORCE_EPS2);
	return p2.q * dr * MyPow(rsqrt(dr2), 3);
}

__device__ inline Vector3d Force3d(const gpu::particle3d3& p1, const gpu::particle3d3& p2)
{
	auto dr = p1.center - p2.center;
	auto dr2 = max(norm(dr), CUDA_FORCE_EPS2);
	return cross(p2.q, dr) * MyPow(rsqrt(dr2), 3);
}
#endif

#ifndef __NVCC__
template <typename point_type, typename interaction_type, typename value_type>
void ComputeExact(const std::vector<particle<point_type, value_type>>& source_particles,
	const std::vector<particle<point_type, value_type>>& target_particles,
	std::function<interaction_type(const particle<point_type, value_type>&, const particle<point_type, value_type>&)> func, std::string filename)
{
	size_t num_particles = target_particles.size();
#ifdef FMM_MPI
	auto [shift, end_part] = LocalPart(0, num_particles);
	size_t local_size = end_part - shift;
#else
	size_t shift = 0, local_size = num_particles;
#endif
	std::vector<interaction_type> exact(local_size);

	auto particles_ptr = target_particles.data() + shift;
	tbb::parallel_for(size_t(0), local_size, [&](size_t i)
	{
		const auto& p1 = particles_ptr[i];
		for (const auto& p2 : source_particles)
		{
			exact[i] += func(p1, p2);
		}
	});

#ifdef FMM_MPI
	std::vector<int> sizes(NProc());
	std::vector<int> displs(NProc());
	std::vector<interaction_type> buf(num_particles);
	for (int j = 0; j < NProc(); ++j)
	{
		sizes[j] = LocalPart(0, num_particles, j);
	}
	for (int j = 1; j < NProc(); ++j)
	{
		displs[j] = displs[j - 1] + sizes[j - 1];
	}
	if constexpr (std::is_same_v<interaction_type,double>)
		MPI_Allgatherv(exact.data(), local_size, MPI_DOUBLE,
			buf.data(), sizes.data(), displs.data(), MPI_DOUBLE, MPI_COMM_WORLD);
	if constexpr (std::is_same_v<interaction_type,std::complex<double>>)
		MPI_Allgatherv(exact.data(), local_size, MPI_COMPLEX16,
			buf.data(), sizes.data(), displs.data(), MPI_COMPLEX16, MPI_COMM_WORLD);
	if constexpr (std::is_same_v<interaction_type, Vector3d>)
	{
		MPI_Datatype MPI_VECTOR3;
		MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_VECTOR3);
		MPI_Type_commit(&MPI_VECTOR3);
		MPI_Allgatherv(exact.data(), local_size, MPI_VECTOR3,
			buf.data(), sizes.data(), displs.data(), MPI_VECTOR3, MPI_COMM_WORLD);
	}
#endif

	std::ofstream fout(filename);
	fout.precision(12);
#ifdef FMM_MPI
	for (const auto& x : buf)
		fout << x << "\n";
#else
	for (const auto& x : exact)
		fout << x << "\n";
#endif
}
#endif

/*
namespace gpu {

namespace detail {

template <InteractionType it, typename interaction_type, typename point_type, typename value_type, typename gpu_interaction_type, typename gpu_point_type>
struct ComputeExact {
	static void Compute(const std::vector<fmm::particle<point_type, value_type>>& source_particles, const std::vector<fmm::particle<point_type, value_type>>& target_particles);
}; }

template <InteractionType it, typename point_type, typename value_type>
void ComputeExact(const std::vector<fmm::particle<point_type, value_type>>& source_particles, const std::vector<fmm::particle<point_type, value_type>>& target_particles)
{
	if constexpr (it == fmm::InteractionType::Potential2d)
		fmm::gpu::detail::ComputeExact<fmm::InteractionType::Potential2d, double, fmm::point2d, value_type, double, fmm::gpu::point2d>::Compute(source_particles, target_particles);
	if constexpr (it == fmm::InteractionType::Force2d)
		fmm::gpu::detail::ComputeExact<fmm::InteractionType::Force2d, fmm::point2d, fmm::point2d, value_type, fmm::gpu::point2d, fmm::gpu::point2d>::Compute(source_particles, target_particles);
	if constexpr (it == fmm::InteractionType::Potential3d)
		fmm::gpu::detail::ComputeExact<fmm::InteractionType::Potential3d, double, fmm::point3d, value_type, double, fmm::gpu::point3d>::Compute(source_particles, target_particles);
	if constexpr (it == fmm::InteractionType::Force3d)
		fmm::gpu::detail::ComputeExact<fmm::InteractionType::Force3d, fmm::point3d, fmm::point3d, value_type, fmm::gpu::point3d, fmm::gpu::point3d>::Compute(source_particles, target_particles);
}

template <InteractionType it, typename point_type, typename value_type>
void ComputeExact(const std::vector<fmm::particle<point_type, value_type>>& particles)
{
	fmm::gpu::ComputeExact<it, point_type, value_type>(particles, particles);
}

}
*/

template <InteractionType it, typename point_type, typename value_type>
void ComputeExact(const std::vector<particle<point_type, value_type>>& source_particles, const std::vector<particle<point_type, value_type>>& target_particles)
{
	if constexpr (it == fmm::InteractionType::Potential2d) {
		auto foo = [](const particle2d& p1, const particle2d& p2) {return Potential2d(p1, p2); };
		ComputeExact<point2d, double, value_type>(source_particles, target_particles, foo, "potential2d_exact.txt");
	}
	if constexpr (it == fmm::InteractionType::Force2d) {
		auto foo = [](const particle2d& p1, const particle2d& p2) {return Force2d(p1, p2); };
		ComputeExact<point2d, point2d, value_type>(source_particles, target_particles, foo, "force2d_exact.txt");
	}
	if constexpr (it == fmm::InteractionType::Potential3d)
	{
		auto foo = [](const particle3d& p1, const particle3d& p2) {return Potential3d(p1, p2); };
		ComputeExact<point3d, double, value_type>(source_particles, target_particles, foo, "potential3d_exact.txt");
	}
	if constexpr (it == fmm::InteractionType::Force3d)
	{
		auto foo = [](const particle<point_type, value_type>& p1, const particle<point_type, value_type>& p2) {return Force3d(p1, p2); };
		ComputeExact<point3d, Vector3d, value_type>(source_particles, target_particles, foo, "force3d_exact.txt");
	}
	if constexpr (it == fmm::InteractionType::VelDipole3d)
	{
		auto foo = [](const particle<point_type, value_type>& p1, const particle<point_type, value_type>& p2) {return velDipole3d(p1, p2); };
		ComputeExact<point3d, Vector3d, value_type>(source_particles, target_particles, foo, "velDipole3d_exact.txt");
	}
	if constexpr (it == fmm::InteractionType::MomDipole3d)
	{
		auto foo = [](const particle<point_type, value_type>& p1, const particle<point_type, value_type>& p2) {return momDipole3d(p1, p2); };
		ComputeExact<point3d, Vector3d, value_type>(source_particles, target_particles, foo, "momDipole3d_exact.txt");
	}
}

template <InteractionType it, typename point_type, typename value_type>
void ComputeExact(const std::vector<particle<point_type, value_type>>& particles)
{
	ComputeExact<it, point_type, value_type>(particles, particles);
}

template <InteractionType it, typename T>
void ReadError(const std::vector<T>& res)
{
	size_t num_particles = res.size();
	std::vector<T> exact(num_particles);

	std::string filename;
	switch (it)
	{
	case fmm::InteractionType::Potential2d:
		filename = "potential2d_exact.txt";
		break;
	case fmm::InteractionType::Force2d:
		filename = "force2d_exact.txt";
		break;
	case fmm::InteractionType::Potential3d:
		filename = "potential3d_exact.txt";
		break;
	case fmm::InteractionType::Force3d:
		filename = "force3d_exact.txt";
		break;
	case fmm::InteractionType::VelDipole3d:
		filename = "velDipole3d_exact.txt";
		break;
	default:
		break;
	}

	std::ifstream fin(filename);
	for (auto& x : exact)
		fin >> x;

	double max_error = 0.0;
	double err1 = 0.0, err2 = 0.0;
	double l2error = 0.0;
//#pragma omp parallel for reduction(+: err1, err2, l2error) reduction(max: max_error)
	for (int i = 0; i < num_particles; ++i)
	{
		double err = 0;
		if constexpr ((it != InteractionType::Force3d) && (it != InteractionType::VelDipole3d))
			err = std::abs(exact[i] - res[i]);
		else
			err = abs(exact[i] - res[i]);
		max_error = std::max(max_error, err);
		
		err1 += err;
		
		if constexpr ((it != InteractionType::Force3d) && (it != InteractionType::VelDipole3d))
		{
			err2 += std::abs(exact[i]);
			l2error += (err / std::abs(exact[i])) * (err / std::abs(exact[i]));
		}
		else
		{
			err2 += abs(exact[i]);
			l2error += (err / abs(exact[i])) * (err / abs(exact[i]));
		}

	}
	l2error /= num_particles;
	//std::cout << "--------------------------------" << std::endl;
	//std::cout << "max error = " << max_error << std::endl;
	//std::cout << "l2 error = " << sqrt(l2error) << std::endl;
	std::cout << "FMM relative error = " << err1 / err2 << std::endl;
	//std::cout << "--------------------------------" << std::endl;
}


inline void ReadErrorHyb(const std::vector<std::array<double, 2>>& res)
{
	size_t num_particles = res.size();
	std::vector<std::complex<double>> exact(num_particles);

	std::string filename = "force2d_exact.txt";

	std::ifstream fin(filename);
	for (auto& x : exact)
		fin >> x;

	double max_error = 0.0;
	double err1 = 0.0, err2 = 0.0;
	double l2error = 0.0;
	//#pragma omp parallel for reduction(+: err1, err2, l2error) reduction(max: max_error)

	for (int i = 0; i < num_particles; ++i)
	{
		double err = 0;

		std::array<double, 2> hyb = { res[i][1] * 6.283185307179586476925286766559, res[i][0] * 6.283185307179586476925286766559 };

		err = sqrt((exact[i].real() - hyb[0]) * (exact[i].real() - hyb[0]) +
			(exact[i].imag() - hyb[1]) * (exact[i].imag() - hyb[1]));
		
		max_error = std::max(max_error, err);

		err1 += err;


			err2 += std::abs(exact[i]);
			l2error += (err / std::abs(exact[i])) * (err / std::abs(exact[i]));
	

	}
	l2error /= num_particles;
	//std::cout << "--------------------------------" << std::endl;
	//std::cout << "max error = " << max_error << std::endl;
	//std::cout << "l2 error = " << sqrt(l2error) << std::endl;
	std::cout << "Hybrid relative error = " << err1 / err2 << std::endl;
	std::cout << "--------------------------------" << std::endl;
}

} // fmm
