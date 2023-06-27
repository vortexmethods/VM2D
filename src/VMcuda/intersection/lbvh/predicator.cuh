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

__device__ __host__
inline query_nearest<float> nearest(const float2& point) noexcept
{
    return query_nearest<float>(point);
}

__device__ __host__
inline query_nearest<double> nearest(const double2& point) noexcept
{
    return query_nearest<double>(point);
}


} // lbvh
#endif// LBVH_PREDICATOR_CUH
