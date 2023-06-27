#ifndef LBVH_MORTON_CODE_CUH
#define LBVH_MORTON_CODE_CUH
#include <vector_types.h>
#include <cuda_runtime.h>
#include <cstdint>

namespace lbvh
{

__device__ __host__
inline std::uint32_t expand_bits(std::uint32_t v) noexcept
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

//Аналогично для 2d случая
//https://stackoverflow.com/questions/30539347/2d-morton-code-encode-decode-64bits
//Expands a 10 - bit integer into 20 bits
// by inserting 1 zero after each bit.
__device__ __host__
inline std::uint32_t expand_bits2(std::uint32_t v) noexcept
{
    //v = (v | (v << 16)) & 0x0000FFFF;
    v = (v | (v << 8)) & 0x00FF00FF;
    v = (v | (v << 4)) & 0x0F0F0F0F;
    v = (v | (v << 2)) & 0x33333333;
    v = (v | (v << 1)) & 0x55555555;
    return v;
}
// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
__device__ __host__
inline std::uint32_t morton_code(float4 xyz, float resolution = 1024.0f) noexcept
{
    xyz.x = ::fminf(::fmaxf(xyz.x * resolution, 0.0f), resolution - 1.0f);
    xyz.y = ::fminf(::fmaxf(xyz.y * resolution, 0.0f), resolution - 1.0f);
    xyz.z = ::fminf(::fmaxf(xyz.z * resolution, 0.0f), resolution - 1.0f);
    const std::uint32_t xx = expand_bits(static_cast<std::uint32_t>(xyz.x));
    const std::uint32_t yy = expand_bits(static_cast<std::uint32_t>(xyz.y));
    const std::uint32_t zz = expand_bits(static_cast<std::uint32_t>(xyz.z));
    return xx * 4 + yy * 2 + zz;
}

__device__ __host__
inline std::uint32_t morton_code2(float2 xy, float resolution = 1024.0f) noexcept
{
    xy.x = ::fminf(::fmaxf(xy.x * resolution, 0.0f), resolution - 1.0f);
    xy.y = ::fminf(::fmaxf(xy.y * resolution, 0.0f), resolution - 1.0f);
    
    const std::uint32_t xx = expand_bits2(static_cast<std::uint32_t>(xy.x));
    const std::uint32_t yy = expand_bits2(static_cast<std::uint32_t>(xy.y));
    
    return xx * 2 + yy;  //xx * 4 + yy * 2 + zz;
}

__device__ __host__
inline std::uint32_t morton_code(double4 xyz, double resolution = 1024.0) noexcept
{
    xyz.x = ::fmin(::fmax(xyz.x * resolution, 0.0), resolution - 1.0);
    xyz.y = ::fmin(::fmax(xyz.y * resolution, 0.0), resolution - 1.0);
    xyz.z = ::fmin(::fmax(xyz.z * resolution, 0.0), resolution - 1.0);
    const std::uint32_t xx = expand_bits(static_cast<std::uint32_t>(xyz.x));
    const std::uint32_t yy = expand_bits(static_cast<std::uint32_t>(xyz.y));
    const std::uint32_t zz = expand_bits(static_cast<std::uint32_t>(xyz.z));
    return xx * 4 + yy * 2 + zz;
}


__device__ __host__
inline std::uint32_t morton_code2(double2 xy, double resolution = 1024.0) noexcept
{
    xy.x = ::fmin(::fmax(xy.x * resolution, 0.0), resolution - 1.0);
    xy.y = ::fmin(::fmax(xy.y * resolution, 0.0), resolution - 1.0);
    
    const std::uint32_t xx = expand_bits2(static_cast<std::uint32_t>(xy.x));
    const std::uint32_t yy = expand_bits2(static_cast<std::uint32_t>(xy.y));
    
    return xx * 2 + yy;  //xx * 4 + yy * 2 + zz;
}


__device__
inline int common_upper_bits(const unsigned int lhs, const unsigned int rhs) noexcept
{
    return ::__clz(lhs ^ rhs);
}
__device__
inline int common_upper_bits(const unsigned long long int lhs, const unsigned long long int rhs) noexcept
{
    return ::__clzll(lhs ^ rhs);
}

} // lbvh
#endif// LBVH_MORTON_CODE_CUH
