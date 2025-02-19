#ifndef LBVH_PREDICATOR_CUH
#define LBVH_PREDICATOR_CUH
#include "aabb.cuh"

namespace lbvh
{

template<typename Real>
struct query_overlap
{
    __device__ __host__
    query_overlap(const aabb<Real>& tgt): target(tgt) {}

    query_overlap()  = default;
    ~query_overlap() = default;
    query_overlap(const query_overlap&) = default;
    query_overlap(query_overlap&&)      = default;
    query_overlap& operator=(const query_overlap&) = default;
    query_overlap& operator=(query_overlap&&)      = default;

    __device__ __host__
    inline bool operator()(const aabb<Real>& box) noexcept
    {
        return intersects(box, target);
    }

    aabb<Real> target;
};

template<typename Real>
__device__ __host__
query_overlap<Real> overlaps(const aabb<Real>& region) noexcept
{
    return query_overlap<Real>(region);
}

template<typename Real>
struct query_nearest
{
    // float4/double4
    using vector_type = typename vector_of<Real>::type;

    __device__ __host__
    query_nearest(const vector_type& tgt): target(tgt) {}

    query_nearest()  = default;
    ~query_nearest() = default;
    query_nearest(const query_nearest&) = default;
    query_nearest(query_nearest&&)      = default;
    query_nearest& operator=(const query_nearest&) = default;
    query_nearest& operator=(query_nearest&&)      = default;

    vector_type target;
};


template<typename Real>
struct query_overlap2
{
    __device__ __host__
    query_overlap2(const aabb2<Real>& tgt): target(tgt) {}

    query_overlap2()  = default;
    ~query_overlap2() = default;
    query_overlap2(const query_overlap2&) = default;
    query_overlap2(query_overlap2&&)      = default;
    query_overlap2& operator=(const query_overlap2&) = default;
    query_overlap2& operator=(query_overlap2&&)      = default;

    __device__ __host__
    inline bool operator()(const aabb2<Real>& box) noexcept
    {
        return intersects2(box, target);
    }

    aabb2<Real> target;
};

template<typename Real>
__device__ __host__
query_overlap2<Real> overlaps2(const aabb2<Real>& region) noexcept
{
    return query_overlap2<Real>(region);
}

template<typename Real>
struct query_nearest2
{
    // float4/double4
    using vector_type2 = typename vector_of2<Real>::type;

    __device__ __host__
    query_nearest2(const vector_type2& tgt): target(tgt) {}

    query_nearest2()  = default;
    ~query_nearest2() = default;
    query_nearest2(const query_nearest2&) = default;
    query_nearest2(query_nearest2&&)      = default;
    query_nearest2& operator=(const query_nearest2&) = default;
    query_nearest2& operator=(query_nearest2&&)      = default;

    vector_type2 target;
};

template<typename Real>
struct query_nearest2ray
{
    // float4/double4
    using vector_type2 = typename vector_of2<Real>::type;

    __device__ __host__
        query_nearest2ray(const vector_type2& start_, const vector_type2& finish_) : start(start_), finish(finish_) {}

    query_nearest2ray() = default;
    ~query_nearest2ray() = default;
    query_nearest2ray(const query_nearest2ray&) = default;
    query_nearest2ray(query_nearest2ray&&) = default;
    query_nearest2ray& operator=(const query_nearest2ray&) = default;
    query_nearest2ray& operator=(query_nearest2ray&&) = default;

    vector_type2 start;
    vector_type2 finish;
};


__device__ __host__
inline query_nearest<float> nearest(const float4& point) noexcept
{
    return query_nearest<float>(point);
}

__device__ __host__
inline query_nearest<float> nearest(const float3& point) noexcept
{
    return query_nearest<float>(make_float4(point.x, point.y, point.z, 0.0f));
}

__device__ __host__
inline query_nearest2<float> nearest2(const float2& point) noexcept
{
    return query_nearest2<float>(make_float2(point.x, point.y));
}

__device__ __host__
inline query_nearest<double> nearest(const double4& point) noexcept
{
    return query_nearest<double>(point);
}

__device__ __host__
inline query_nearest<double> nearest(const double3& point) noexcept
{
    return query_nearest<double>(make_double4(point.x, point.y, point.z, 0.0));
}

__device__ __host__
inline query_nearest2<double> nearest2(const double2& point) noexcept
{
    return query_nearest2<double>(make_double2(point.x, point.y));
}

} // lbvh
#endif// LBVH_PREDICATOR_CUH
