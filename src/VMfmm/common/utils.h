#pragma once
#include "defs.h"
#include <complex>
#include <array>
#include <vector>
#include <cassert>
//#include <cuda/std/complex>
#include "omp.h"

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

namespace fmm {

enum class InteractionType {
	Potential2d,
	Force2d,
	Potential3d,
	Force3d,
	VelDipole3d,
	MomDipole3d
};

template <typename point, typename value>
struct particle
{
	point center;
	value q;
	size_t morton_code;
};

template <typename T>
class Vector3 : public std::array<T, 3>
{
public:

	constexpr Vector3() = default;
	constexpr Vector3(const Vector3<T>& other) = default;
	constexpr Vector3(Vector3<T>&& other) = default;
	constexpr Vector3<T>& operator=(const Vector3<T>& other) = default;
	constexpr Vector3<T>& operator=(Vector3<T>&& other) = default;
	constexpr __HOST__ __DEVICE__ Vector3(T x, T y, T z) : std::array<T, 3> {x, y, z} {};

	constexpr __HOST__ __DEVICE__ Vector3<T> operator+(const Vector3<T>& rhs) const {return Vector3<T>(std::array<T, 3>::operator[](0)+rhs[0], std::array<T, 3>::operator[](1)+rhs[1], std::array<T, 3>::operator[](2)+rhs[2]);}
	constexpr __HOST__ __DEVICE__ Vector3<T> operator-(const Vector3<T>& rhs) const {return Vector3<T>(std::array<T, 3>::operator[](0)-rhs[0], std::array<T, 3>::operator[](1)-rhs[1], std::array<T, 3>::operator[](2)-rhs[2]);}
	constexpr __HOST__ __DEVICE__ Vector3<T> operator*(const double& x) const {return Vector3<T>(std::array<T, 3>::operator[](0)*x, std::array<T, 3>::operator[](1)*x, std::array<T, 3>::operator[](2)*x);}
	constexpr __HOST__ __DEVICE__ Vector3<T> operator/(const double& x) const {return Vector3<T>(std::array<T, 3>::operator[](0)/x, std::array<T, 3>::operator[](1)/x, std::array<T, 3>::operator[](2)/x);}

	constexpr __HOST__ __DEVICE__ void operator+=(const Vector3<T>& rhs) { std::array<T, 3>::operator[](0)+=rhs[0]; std::array<T, 3>::operator[](1)+=rhs[1]; std::array<T, 3>::operator[](2)+=rhs[2];}
	constexpr __HOST__ __DEVICE__ void operator-=(const Vector3<T>& rhs) { std::array<T, 3>::operator[](0)-=rhs[0]; std::array<T, 3>::operator[](1)-=rhs[1]; std::array<T, 3>::operator[](2)-=rhs[2];}
	constexpr __HOST__ __DEVICE__ void operator*=(const double& x) { std::array<T, 3>::operator[](0)*=x; std::array<T, 3>::operator[](1)*=x; std::array<T, 3>::operator[](2)*=x;}
	constexpr __HOST__ __DEVICE__ void operator/=(const double& x) { std::array<T, 3>::operator[](0)/=x; std::array<T, 3>::operator[](1)/=x; std::array<T, 3>::operator[](2)/=x;}

	constexpr __HOST__ __DEVICE__ bool operator==(const Vector3<T>& rhs) const { return std::array<T, 3>::operator[](0) == rhs[0] && std::array<T, 3>::operator[](1) == rhs[1] && std::array<T, 3>::operator[](2) == rhs[2]; }
	constexpr __HOST__ __DEVICE__ bool operator!=(const Vector3<T>& rhs) const { return !(*this == rhs); }
};

using Vector3d = Vector3<double>;
using Matrix3d = std::array<Vector3d, 3>;
using Vector3cd = Vector3<std::complex<double>>;
//using cuda_Vector3cd = Vector3<cuda::std::complex<double>>;

template <typename T>
inline __HOST__ __DEVICE__ constexpr Vector3<T> operator*(const double x, const Vector3<T>& rhs) {return Vector3<T>(rhs[0]*x, rhs[1]*x, rhs[2]*x);}

inline __HOST__ __DEVICE__ constexpr Vector3cd operator*(const std::complex<double>& x, const Vector3d& rhs) { return Vector3cd(rhs[0] * x, rhs[1] * x, rhs[2] * x); }
inline __HOST__ __DEVICE__ constexpr Vector3cd operator*(const Vector3d& rhs, const std::complex<double>& x) { return Vector3cd(rhs[0] * x, rhs[1] * x, rhs[2] * x); }
inline __HOST__ __DEVICE__ constexpr Vector3cd operator*(const std::complex<double>& x, const Vector3cd& rhs) { return Vector3cd(rhs[0] * x, rhs[1] * x, rhs[2] * x); }
inline __HOST__ __DEVICE__ constexpr Vector3cd operator*(const Vector3cd& rhs, const std::complex<double>& x) { return Vector3cd(rhs[0] * x, rhs[1] * x, rhs[2] * x); }

inline __HOST__ __DEVICE__ /*constexpr*/ void operator*=(Vector3cd& lhs, const std::complex<double>& x) { lhs[0] *= x; lhs[1] *= x; lhs[2] *= x; }

//inline __HOST__ __DEVICE__ cuda_Vector3cd operator*(const cuda::std::complex<double>& x, const Vector3d& rhs) { return cuda_Vector3cd(rhs[0] * x, rhs[1] * x, rhs[2] * x); }
//inline __HOST__ __DEVICE__ cuda_Vector3cd operator*(const Vector3d& rhs, const cuda::std::complex<double>& x) { return cuda_Vector3cd(rhs[0] * x, rhs[1] * x, rhs[2] * x); }
//inline __HOST__ __DEVICE__ cuda_Vector3cd operator*(const cuda::std::complex<double>& x, const cuda_Vector3cd& rhs) { return cuda_Vector3cd(rhs[0] * x, rhs[1] * x, rhs[2] * x); }
//inline __HOST__ __DEVICE__ cuda_Vector3cd operator*(const cuda_Vector3cd& rhs, const cuda::std::complex<double>& x) { return cuda_Vector3cd(rhs[0] * x, rhs[1] * x, rhs[2] * x); }

//inline __HOST__ __DEVICE__ void operator*=(cuda_Vector3cd& lhs, const cuda::std::complex<double>& x) { lhs[0] *= x; lhs[1] *= x; lhs[2] *= x; }

inline __HOST__ __DEVICE__ /*constexpr*/ void operator+=(Vector3cd& lhs, const Vector3d& rhs) {
	lhs[0].real(lhs[0].real() + rhs[0]);
	lhs[1].real(lhs[1].real() + rhs[1]);
	lhs[2].real(lhs[2].real() + rhs[2]);
}

//inline __HOST__ __DEVICE__ void operator+=(cuda_Vector3cd& lhs, const Vector3d& rhs) {
//	lhs[0].real(lhs[0].real() + rhs[0]);
//	lhs[1].real(lhs[1].real() + rhs[1]);
//	lhs[2].real(lhs[2].real() + rhs[2]);
//}

template <typename complex>
inline __HOST__ __DEVICE__ constexpr Vector3d real(const Vector3cd& rhs) { return Vector3d(rhs[0].real(), rhs[1].real(), rhs[2].real()); }
//template <typename complex>
//inline __HOST__ __DEVICE__ constexpr Vector3d real(const cuda_Vector3cd& rhs) { return Vector3d(rhs[0].real(), rhs[1].real(), rhs[2].real()); }

template <typename complex>
inline __HOST__ __DEVICE__ constexpr Vector3d imag(const Vector3cd& rhs) { return Vector3d(rhs[0].imag(), rhs[1].imag(), rhs[2].imag()); }
//template <typename complex>
//inline __HOST__ __DEVICE__ constexpr Vector3d imag(const cuda_Vector3cd& rhs) { return Vector3d(rhs[0].imag(), rhs[1].imag(), rhs[2].imag()); }

template <typename T>
inline __HOST__ __DEVICE__ constexpr Vector3<T> cross(const Vector3<T>& lhs, const Vector3<T>& rhs) { return Vector3<T>(lhs[1] * rhs[2] - lhs[2] * rhs[1], lhs[2] * rhs[0] - lhs[0] * rhs[2], lhs[0] * rhs[1] - lhs[1] * rhs[0]); }

inline __HOST__ __DEVICE__ constexpr double dot(const Vector3d& lhs, const Vector3d& rhs) { return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2]; }

template <typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector3<T>& rhs) { os << rhs[0] << " " << rhs[1] << " " << rhs[2]; return os; }

template <typename T>
inline std::istream& operator>>(std::istream& is, Vector3<T>& rhs) {is >> rhs[0] >> rhs[1] >> rhs[2]; return is;}

template <typename T>
inline __HOST__ __DEVICE__ double norm(const Vector3<T>& vec) {return vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];}

template <typename T>
inline __HOST__ __DEVICE__ double abs(const Vector3<T>& vec) {return sqrt(norm(vec));}

inline __HOST__ __DEVICE__ int doublehash(double val)
{
	return int(val * 100);
}

inline __HOST__ __DEVICE__ int doublehash(double val1, double val2)
{
	return 100 * int(val1 * 100) + int(val2 * 100);
}

template<typename T>
struct treevector : public std::vector<T>
{
	treevector() = default;
	treevector(const treevector&) = default;
	treevector(treevector&&) = default;
	treevector& operator=(const treevector&) = default;
	treevector& operator=(treevector&&) = default;
	constexpr treevector(int n) : std::vector<T>(n * (n + 1) / 2) {}
	T& operator()(int n, int m)
	{
		return this->operator[](n * (n + 1) / 2 + m);
	}
	const T& operator()(int n, int m) const
	{
		return this->operator[](n * (n + 1) / 2 + m);
	}
};

using point2d = std::complex<double>;
using point3d = Vector3d;
using particle2d = particle<point2d, double>;
using particle3d = particle<point3d, double>;
using particle3d3 = particle<point3d, point3d>;

} // fmm
