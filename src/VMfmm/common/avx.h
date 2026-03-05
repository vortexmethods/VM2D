#pragma once
#include <immintrin.h>
#include <type_traits>

namespace fmm {

	#define AVX256

	namespace detail {
#if defined(AVX128)
		inline const int avx_vec_length = 2;
#elif defined(AVX256)
		inline const int avx_vec_length = 4;
#elif defined(AVX512)
		inline const int avx_vec_length = 8;
#endif
	}

	template <typename vec_type>
	vec_type avx_zero_vec()
	{
		if constexpr (std::is_same_v<vec_type, __m128d>)
			return _mm_setzero_pd();
		if constexpr (std::is_same_v<vec_type, __m256d>)
			return _mm256_setzero_pd();
		if constexpr (std::is_same_v<vec_type, __m512d>)
			return _mm512_setzero_pd();
	}

	template <typename vec_type>
	inline void avx_store(double* data, vec_type& v)
	{
		if constexpr (std::is_same_v<vec_type, __m128d>)
			_mm_store_pd(data, v);
		if constexpr (std::is_same_v<vec_type, __m256d>)
			_mm256_store_pd(data, v);
		if constexpr (std::is_same_v<vec_type, __m512d>)
			_mm512_store_pd(data, v);
	}

	template <typename vec_type>
	inline vec_type avx_add(const vec_type& v1, const vec_type& v2)
	{
		if constexpr (std::is_same_v<vec_type, __m128d>)
			return _mm_add_pd(v1, v2);
		if constexpr (std::is_same_v<vec_type, __m256d>)
			return _mm256_add_pd(v1, v2);
		if constexpr (std::is_same_v<vec_type, __m512d>)
			return _mm512_add_pd(v1, v2);
	}

	template <typename vec_type>
	inline vec_type avx_sub(const vec_type& v1, const vec_type& v2)
	{
		if constexpr (std::is_same_v<vec_type, __m128d>)
			return _mm_sub_pd(v1, v2);
		if constexpr (std::is_same_v<vec_type, __m256d>)
			return _mm256_sub_pd(v1, v2);
		if constexpr (std::is_same_v<vec_type, __m512d>)
			return _mm512_sub_pd(v1, v2);
	}

	template <typename vec_type>
	inline vec_type avx_mul(const vec_type& v1, const vec_type& v2)
	{
		if constexpr (std::is_same_v<vec_type, __m128d>)
			return _mm_mul_pd(v1, v2);
		if constexpr (std::is_same_v<vec_type, __m256d>)
			return _mm256_mul_pd(v1, v2);
		if constexpr (std::is_same_v<vec_type, __m512d>)
			return _mm512_mul_pd(v1, v2);
	}

	inline double avx_hsum(__m128d v) {
		return  _mm_cvtsd_f64(_mm_add_sd(v, _mm_unpackhi_pd(v, v)));
	}

	inline double avx_hsum(__m256d v) {
		__m128d vlow = _mm256_castpd256_pd128(v);
		vlow = _mm_add_pd(vlow, _mm256_extractf128_pd(v, 1));
		return  _mm_cvtsd_f64(_mm_add_sd(vlow, _mm_unpackhi_pd(vlow, vlow)));
	}

	inline double avx_hsum(__m512d v) {
		__m256d vlow = _mm512_castpd512_pd256(v);
		return avx_hsum(_mm256_add_pd(vlow, _mm512_extractf64x4_pd(v, 1)));
	}

	template <typename vec_type>
	inline void avx_set_range(vec_type& vec, double* val)
	{
		if constexpr (std::is_same_v<vec_type, __m128d>)
			vec = _mm_set_pd(val[1], val[0]);
		if constexpr (std::is_same_v<vec_type, __m256d>)
			vec = _mm256_set_pd(val[3], val[2], val[1], val[0]);
		if constexpr (std::is_same_v<vec_type, __m512d>)
			vec = _mm512_set_pd(val[7], val[6], val[5], val[4], val[3], val[2], val[1], val[0]);
	}

	template <typename vec_type>
	inline void avx_set_constant(vec_type& vec, double val)
	{
		if constexpr (std::is_same_v<vec_type, __m128d>)
			vec = _mm_set_pd(val, val);
		if constexpr (std::is_same_v<vec_type, __m256d>)
			vec = _mm256_set_pd(val, val, val, val);
		if constexpr (std::is_same_v<vec_type, __m512d>)
			vec = _mm512_set_pd(val, val, val, val, val, val, val, val);
	}

	inline void avx_inv_dr(const __m512d& oneVec, const __m512d& epsVec,
		const __m512d& x_target, const  __m512d& y_target, const __m512d& z_target,
		__m512d& temp1, __m512d& temp2, __m512d& temp3,
		__m512d& dx, __m512d& dy, __m512d& dz, __m512d& invdr, __m512d& invdr2)
	{
		dx = _mm512_sub_pd(x_target, temp1);
		dy = _mm512_sub_pd(y_target, temp2);
		dz = _mm512_sub_pd(z_target, temp3);
		temp1 = _mm512_mul_pd(dx, dx);
		temp2 = _mm512_mul_pd(dy, dy);
		temp3 = _mm512_mul_pd(dz, dz);
		invdr = _mm512_add_pd(temp1, temp2);
		invdr = _mm512_add_pd(invdr, temp3);
		invdr = _mm512_max_pd(invdr, epsVec);
		invdr2 = _mm512_div_pd(oneVec, invdr);
		invdr = _mm512_sqrt_pd(invdr2);
	}

	inline void avx_inv_dr(const __m256d& oneVec, const __m256d& epsVec,
		const __m256d& x_target, const __m256d& y_target, const __m256d& z_target,
		__m256d& temp1, __m256d& temp2, __m256d& temp3,
		__m256d& dx, __m256d& dy, __m256d& dz, __m256d& invdr, __m256d& invdr2)
	{
		dx = _mm256_sub_pd(x_target, temp1);
		dy = _mm256_sub_pd(y_target, temp2);
		dz = _mm256_sub_pd(z_target, temp3);
		temp1 = _mm256_mul_pd(dx, dx);
		temp2 = _mm256_mul_pd(dy, dy);
		temp3 = _mm256_mul_pd(dz, dz);

		invdr = _mm256_add_pd(temp1, temp2);
		invdr = _mm256_add_pd(invdr, temp3);

		invdr = _mm256_max_pd(invdr, epsVec);
		invdr2 = _mm256_div_pd(oneVec, invdr);
		invdr = _mm256_sqrt_pd(invdr2);
	}

	inline void avx_inv_dr(const __m128d& oneVec, const __m128d& epsVec,
		const __m128d& x_target, const __m128d& y_target, const __m128d& z_target,
		__m128d& temp1, __m128d& temp2, __m128d& temp3,
		__m128d& dx, __m128d& dy, __m128d& dz, __m128d& invdr, __m128d& invdr2)
	{
		dx = _mm_sub_pd(x_target, temp1);
		dy = _mm_sub_pd(y_target, temp2);
		dz = _mm_sub_pd(z_target, temp3);
		temp1 = _mm_mul_pd(dx, dx);
		temp2 = _mm_mul_pd(dy, dy);
		temp3 = _mm_mul_pd(dz, dz);

		invdr = _mm_add_pd(temp1, temp2);
		invdr = _mm_add_pd(invdr, temp3);
		invdr = _mm_max_pd(invdr, epsVec);
		invdr2 = _mm_div_pd(oneVec, invdr);
		invdr = _mm_sqrt_pd(invdr2);
	}
}
