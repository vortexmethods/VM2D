#pragma once

#include <limits>

// FMM_MPI is controlled by CMake (passed as -DFMM_MPI via target_compile_definitions)
// Uncomment below to force-enable MPI regardless of CMake:
//#define FMM_MPI
#define FMM_CONSTEXPR_MATH // fast CUDA solver, but limits max multipole num (min error 10^-7)

#ifdef FMM_CONSTEXPR_MATH
#define FMM_CONSTEXPR constexpr
#else
#define FMM_CONSTEXPR
#endif

#ifdef __NVCC__
#define __HOST__ __host__
#define __DEVICE__ __device__
#else
#define __HOST__
#define __DEVICE__
#endif

namespace fmm {

	const int FMM_AUTO = std::numeric_limits<int>::max();

	constexpr double FORCE_EPS = 0.01;
	constexpr double FORCE_EPS2 = FORCE_EPS * FORCE_EPS;
	constexpr double FORCE_EPS3 = FORCE_EPS2 * FORCE_EPS;
	constexpr double FORCE_EPS4 = FORCE_EPS2 * FORCE_EPS2;
	constexpr double FORCE_EPS5 = FORCE_EPS3 * FORCE_EPS2;
	constexpr double FORCE_EPS6 = FORCE_EPS3 * FORCE_EPS3;

	//__DEVICE__ constexpr double CUDA_FORCE_EPS = FORCE_EPS;
	//__DEVICE__ constexpr double CUDA_FORCE_EPS2 = CUDA_FORCE_EPS * CUDA_FORCE_EPS;
	//__DEVICE__ constexpr double CUDA_FORCE_EPS3 = CUDA_FORCE_EPS2 * CUDA_FORCE_EPS;
	//__DEVICE__ constexpr double CUDA_FORCE_EPS5 = CUDA_FORCE_EPS3 * CUDA_FORCE_EPS2;

namespace detail {

	const int _3d_MAX_MULTIPOLE_NUM = 18; // set max 18 if constexpr math
	const int _2d_MAX_MULTIPOLE_NUM = 9; // set max 30 if constexpr math
	


} // detail

} // fmm
