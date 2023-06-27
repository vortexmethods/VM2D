#ifndef LBVH_UTILITY_CUH
#define LBVH_UTILITY_CUH
#include <vector_types.h>
#include <math_constants.h>

namespace lbvh
{

template<typename T> struct vector_of;
template<> struct vector_of<float>  {using type = float4;};
template<> struct vector_of<double> {using type = double4;};

template<typename T>
using vector_of_t = typename vector_of<T>::type;

template<typename T> struct vector_of2;
template<> struct vector_of2<float> { using type = float2; };
template<> struct vector_of2<double> { using type = double2; };

template<typename T>
using vector_of2_t = typename vector_of2<T>::type;


//template<typename T> struct vector_ofS;
//template<> struct vector_ofS<float> { using type = float2; };

//template<typename T>
//using vector_ofS_t = typename vector_ofS<T>::type;



template<typename T>
__device__
inline T infinity() noexcept;

template<>
__device__
inline float  infinity<float >() noexcept {return CUDART_INF_F;}
template<>
__device__
inline double infinity<double>() noexcept {return CUDART_INF;}


template<typename T>
struct segment
{
    vector_of2_t<T> p1, p2;
};

template<typename T>
inline __host__ __device__ 
segment<T> make_segment(vector_of2_t<T> p1_, vector_of2_t<T> p2_)
{
    return segment<T>{ p1_, p2_ };
}


} // lbvh
#endif// LBVH_UTILITY_CUH
