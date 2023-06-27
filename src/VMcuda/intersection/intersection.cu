#include "lbvh.cuh"
#include <random>
#include <vector>
#include <thrust/random.h>
#include <cmath>

#include <omp.h>

#include "Point2D.h"


struct aabbS_getter
{
    __device__
        lbvh::aabb2<float> operator()(const lbvh::segment<float> S) const noexcept
    {
        lbvh::aabb2<float> retval;
        retval.upper.x = ::fmaxf(S.p1.x, S.p2.x);
        retval.upper.y = ::fmaxf(S.p1.y, S.p2.y);
        retval.lower.x = ::fminf(S.p1.x, S.p2.x);
        retval.lower.y = ::fminf(S.p1.y, S.p2.y);
        return retval;
    }
    __device__
        lbvh::aabb2<double> operator()(const lbvh::segment<double> S) const noexcept
    {
        lbvh::aabb2<double> retval;
        retval.upper.x = ::fmax(S.p1.x, S.p2.x);
        retval.upper.y = ::fmax(S.p1.y, S.p2.y);
        retval.lower.x = ::fmin(S.p1.x, S.p2.x);
        retval.lower.y = ::fmin(S.p1.y, S.p2.y);
        return retval;
    }
};


struct distance_calculatorS
{
    //определяет квадрат расстояния от точки до отрезка
    __device__
        thrust::pair<float, int> operator()(const float2 point, const lbvh::segment<float> object) const noexcept
    {
        int location;

        float a = object.p2.y - object.p1.y;
        float b = object.p2.x - object.p1.x;

        float distanceSegment;
        float r_numerator = (point.x - object.p1.x) * b + (point.y - object.p1.y) * a;
        float r_denomenator = b * b + a * a;
        float r = r_numerator / r_denomenator;

        float s = ((object.p1.y - point.y) * b - (object.p1.x - point.x) * a);

        if ((r >= 0) && (r <= 1))
        {
            distanceSegment = s * s / r_denomenator;
            location = 0;
        }
        else
        {
            float dist1 = (point.x - object.p1.x) * (point.x - object.p1.x) +
                (point.y - object.p1.y) * (point.y - object.p1.y);
            float dist2 = (point.x - object.p2.x) * (point.x - object.p2.x) +
                (point.y - object.p2.y) * (point.y - object.p2.y);
            if (dist1 < dist2)
            {
                distanceSegment = dist1;
                location = -1;
            }
            else
            {
                distanceSegment = dist2;
                location = 1;
            }
        }

        return thrust::make_pair(distanceSegment, location);
    }

    //определяет квадрат расстояния от точки до отрезка
    __device__
        thrust::pair<double, int> operator()(const double2 point, const lbvh::segment<double> object) const noexcept
    {
        int location;

        double a = object.p2.y - object.p1.y;
        double b = object.p2.x - object.p1.x;

        double distanceSegment;
        double r_numerator = (point.x - object.p1.x) * b + (point.y - object.p1.y) * a;
        double r_denomenator = b * b + a * a;
        double r = r_numerator / r_denomenator;

        double s = ((object.p1.y - point.y) * b - (object.p1.x - point.x) * a);

        if ((r >= 0) && (r <= 1))
        {
            distanceSegment = s * s / r_denomenator;
            location = 0;
        }
        else
        {
            double dist1 = (point.x - object.p1.x) * (point.x - object.p1.x) +
                (point.y - object.p1.y) * (point.y - object.p1.y);
            double dist2 = (point.x - object.p2.x) * (point.x - object.p2.x) +
                (point.y - object.p2.y) * (point.y - object.p2.y);
            if (dist1 < dist2)
            {
                distanceSegment = dist1;
                location = -1;
            }
            else
            {
                distanceSegment = dist2;
                location = 1;
            }
        }

        return thrust::make_pair(distanceSegment, location);
    }
};





//нормаль к отрезку
inline __device__ float2 normal2d(const lbvh::segment<float> object)
{
    float dx = object.p2.x - object.p1.x;
    float dy = object.p2.y - object.p1.y;
    float len = sqrtf(dx * dx + dy * dy);
    return  { dy / len, -dx / len };
}


inline __device__ double2 normal2d(const lbvh::segment<double> object)
{
    double dx = object.p2.x - object.p1.x;
    double dy = object.p2.y - object.p1.y;
    double len = sqrt(dx * dx + dy * dy);
    return  { dy / len, -dx / len };
}

template <typename T> 
__host__ __device__
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

std::vector<int> lbvh_check_inside(
    int timeStep,
    //const std::vector<Point2D>& vecPositions,
    int nPoints,
    double* devNewpos,
    int nPanels,
    double* devRPtr
)
{     
    double t1 = omp_get_wtime();
    
    thrust::device_vector<lbvh::segment<float>> ps(nPanels);
    thrust::transform(
       thrust::make_counting_iterator<unsigned int>(0),
       thrust::make_counting_iterator<unsigned int>(nPanels),
       ps.begin(),
       [devRPtr] __device__(const unsigned int idx) {
        const float x1 = devRPtr[4 * idx + 0];
        const float y1 = devRPtr[4 * idx + 1];
        const float x2 = devRPtr[4 * idx + 2];
        const float y2 = devRPtr[4 * idx + 3];
        return lbvh::segment<float>{ {x1, y1}, {x2, y2} };
    }
    );
    lbvh::segment<float>* psptr = thrust::raw_pointer_cast(&ps[0]);

    double t2 = omp_get_wtime();

    thrust::device_vector<float2> norms;
    norms.resize(nPanels);
    thrust::transform(
        ps.begin(),
        ps.end(),
        norms.begin(),
        []__device__(const lbvh::segment<float>& s){   
        return normal2d(s);
    }
    );
    float2* normptr = thrust::raw_pointer_cast(&norms[0]);
   
    double t3 = omp_get_wtime();

    lbvh::bvhS<float, lbvh::segment<float>, aabbS_getter> bvhS(ps.begin(), ps.end(), true);
    const auto bvh_dev = bvhS.get_device_reprS();

    double t4 = omp_get_wtime();

    thrust::device_vector<float2> random_points(nPoints);
    thrust::transform(
        thrust::make_counting_iterator<unsigned int>(0),
        thrust::make_counting_iterator<unsigned int>(nPoints),
        random_points.begin(),
        [devNewpos] __device__(const unsigned int idx) {
        const float x = devNewpos[2 * idx + 0];
        const float y = devNewpos[2 * idx + 1];
        return make_float2(x, y);
    }
    );
    float2* ptr = thrust::raw_pointer_cast(&random_points[0]);

    double t5 = omp_get_wtime();

    /////////////////////////////////////

    thrust::device_vector<thrust::pair<int, int>> num;
    num.resize(nPoints);
    auto numptr = thrust::raw_pointer_cast(&num[0]);

    thrust::for_each(thrust::device,
        thrust::make_counting_iterator<unsigned int>(0),
        thrust::make_counting_iterator<unsigned int>(nPoints),
        [bvh_dev, ptr, numptr] __device__(const unsigned int idx) {
        const auto self = ptr[idx];
        const auto nest = lbvh::query_device2(bvh_dev, lbvh::nearest2(self),
            distance_calculatorS());

        if (nest.first == 0xFFFFFFFF)
        {
            printf("???????????????????????????????????????????????????????????????\n");
            printf("x,y = {%.16f, %.16f}\n", self.x, self.y);
        }
        numptr[idx] = nest;
        return;
    });
    
    double t6 = omp_get_wtime();

    thrust::device_vector<int> flag(random_points.size());
    auto flagptr = thrust::raw_pointer_cast(&flag[0]);
    
    thrust::transform(thrust::device,
        thrust::make_counting_iterator<unsigned int>(0),
        thrust::make_counting_iterator<unsigned int>(nPoints),
        flag.begin(),
        [ptr, normptr, numptr, psptr, nPanels]__device__(const unsigned int idx) {
        int pnl = numptr[idx].first;
        int loc = numptr[idx].second;
        float2 norm, vp;

        if (loc > 0)  //1
        {
            norm.x = normptr[pnl].x + normptr[pnl == nPanels - 1 ? 0 : pnl + 1].x;
            norm.y = normptr[pnl].y + normptr[pnl == nPanels - 1 ? 0 : pnl + 1].y;
            vp.x = ptr[idx].x - psptr[pnl].p2.x;
            vp.y = ptr[idx].y - psptr[pnl].p2.y;
        }
        else if (loc == 0)
        {
            norm = normptr[pnl];
            vp.x = ptr[idx].x - (psptr[pnl].p2.x + psptr[pnl].p1.x) * 0.5;
            vp.y = ptr[idx].y - (psptr[pnl].p2.y + psptr[pnl].p1.y) * 0.5;
        }
        else //-1
        {
            norm.x = normptr[pnl].x + normptr[pnl == 0 ? nPanels - 1 : pnl - 1].x;
            norm.y = normptr[pnl].y + normptr[pnl == 0 ? nPanels - 1 : pnl - 1].y;
            vp.x = ptr[idx].x - psptr[pnl].p1.x;
            vp.y = ptr[idx].y - psptr[pnl].p1.y;
        }

        return sgn(norm.x * vp.x + norm.y * vp.y);
    }
    );

    double t7 = omp_get_wtime();

    thrust::device_vector<int> fitpanel(random_points.size());
    thrust::transform(thrust::device,
        thrust::make_counting_iterator<unsigned int>(0),
        thrust::make_counting_iterator<unsigned int>(nPoints),
        fitpanel.begin(),
        [flagptr, numptr]__device__(const unsigned int idx) {
            if (flagptr[idx] > 0)
                return -1;
            else
                return (numptr[idx].first);            
        }
    );

    double t8 = omp_get_wtime();

    ///////////////////////////////////////

    thrust::host_vector<int> fitpanelHost = fitpanel;
    
    double t9 = omp_get_wtime();

    //std::cout << 1000*(t2 - t1) << " " << 1000 * (t3 - t2) << " " << 1000 * (t4 - t3) << " " << 1000 * (t5 - t4) << " " << 1000 * (t6 - t5) << " " << 1000 * (t7 - t6) << " " << 1000 * (t8 - t7) << " " << 1000 * (t9 - t8) << std::endl;

    return std::vector<int>(fitpanelHost.begin(), fitpanelHost.end());
}
