#include "lbvh.cuh"
#include <random>
#include <vector>
#include <thrust/random.h>
#include <cmath>

#include <omp.h>

#include "Point2D.h"


#include <memory>



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



std::shared_ptr<lbvh::bvhS<float, lbvh::segment<float>, aabbS_getter>> ptrTree(nullptr);
std::shared_ptr<thrust::device_vector<float2>> ptrNorm(nullptr);
std::shared_ptr < thrust::device_vector<lbvh::segment<float>>> ptrPs(nullptr);


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




template <typename T>
__device__ int get_line_intersection(T p0_x, T p0_y, T p1_x, T p1_y,
    T p2_x, T p2_y, T p3_x, T p3_y, T* i_x, T* i_y)

{
    T s1_x, s1_y, s2_x, s2_y;
    s1_x = p1_x - p0_x;     s1_y = p1_y - p0_y;
    s2_x = p3_x - p2_x;     s2_y = p3_y - p2_y;

    T s, t;
    T den = (-s2_x * s1_y + s1_x * s2_y);

    T dx = (p0_x - p2_x) / den;
    T dy = (p0_y - p2_y) / den;

    s = (-s1_y * dx + s1_x * dy);
    t = (s2_x * dy - s2_y * dx);

    //printf("s = %f, t = %f\n", s, t);

    if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
    {
        // Collision detected
        if (i_x != NULL)
            *i_x = p0_x + (t * s1_x);
        if (i_y != NULL)
            *i_y = p0_y + (t * s1_y);
        return 1;
    }

    return 0; // No collision
}


struct check_hit_obj_Calculator
{
    //определяет расстояние по лучу до точки пересечения этим лучом отрезка object (или infinity)
    __device__
        float operator()(const float2 start, const float2 finish, const lbvh::segment<float> object) const noexcept
    {
        //printf("D3\n");
        float xint, yint;

        //printf("st.x = %f, st.y = %f, f.x = %f, f.y = %f, op1.x = %f, op1.y = %f, op2.x = %f, op2.y = %f\n", 
        //    start.x, start.y, finish.x, finish.y, object.p1.x, object.p1.y, object.p2.x, object.p2.y);

        int intersect = get_line_intersection(start.x, start.y, finish.x, finish.y, object.p1.x, object.p1.y, object.p2.x, object.p2.y, &xint, &yint);
        if (intersect)
            return (xint - start.x) * (xint - start.x) + (yint - start.y) * (yint - start.y);
        else
            return lbvh::infinity<float>();
    }

    __device__
        float operator()(const double2 start, const double2 finish, const lbvh::segment<double> object) const noexcept
    {
        //printf("D3\n");
        double xint, yint;

        //printf("st.x = %f, st.y = %f, f.x = %f, f.y = %f, op1.x = %f, op1.y = %f, op2.x = %f, op2.y = %f\n", 
        //    start.x, start.y, finish.x, finish.y, object.p1.x, object.p1.y, object.p2.x, object.p2.y);

        int intersect = get_line_intersection(start.x, start.y, finish.x, finish.y, object.p1.x, object.p1.y, object.p2.x, object.p2.y, &xint, &yint);
        if (intersect)
            return (xint - start.x) * (xint - start.x) + (yint - start.y) * (yint - start.y);
        else
            return lbvh::infinity<double>();
    }
};


struct check_hit_node_Calculator
{
    //определяет расстояние по лучу до точки пересечения этим лучом прямоугольника (или infinity)
    __device__
        float operator()(const float2 start, const float2 finish, const lbvh::aabb2<float> object) const noexcept
    {
        if ( (object.lower.x <= start.x) && (start.x <= object.upper.x) && (object.lower.y <= start.y) && (start.y <= object.upper.y) &&
            (object.lower.x <= finish.x) && (finish.x <= object.upper.x) && (object.lower.y <= finish.y) && (finish.y <= object.upper.y) )
            return 0;
        
        float xint, yint, d1, d2, d3, d4;
        int intersect1 = get_line_intersection(start.x, start.y, finish.x, finish.y, object.upper.x,
            object.upper.y, object.upper.x, object.lower.y, &xint, &yint);
        d1 = (intersect1) ? (xint - start.x) * (xint - start.x) + (yint - start.y) * (yint - start.y) : lbvh::infinity<float>();

        int intersect2 = get_line_intersection(start.x, start.y, finish.x, finish.y, object.upper.x,
            object.lower.y, object.lower.x, object.lower.y, &xint, &yint);
        d2 = (intersect2) ? (xint - start.x) * (xint - start.x) + (yint - start.y) * (yint - start.y) : lbvh::infinity<float>();

        int intersect3 = get_line_intersection(start.x, start.y, finish.x, finish.y, object.lower.x, object.lower.y,
            object.lower.x, object.upper.y, &xint, &yint);
        d3 = (intersect3) ? (xint - start.x) * (xint - start.x) + (yint - start.y) * (yint - start.y) : lbvh::infinity<float>();

        int intersect4 = get_line_intersection(start.x, start.y, finish.x, finish.y, object.lower.x, object.upper.y,
            object.upper.x, object.upper.y, &xint, &yint);
        d4 = (intersect4) ? (xint - start.x) * (xint - start.x) + (yint - start.y) * (yint - start.y) : lbvh::infinity<float>();

        //printf("d: %f %f %f %f\n",d1, d2, d3, d4);

        d1 = (d1 < d2) ? d1 : d2;
        d1 = (d1 < d3) ? d1 : d3;
        d1 = (d1 < d4) ? d1 : d4;

        return d1;
    }
};


__host__ __device__
bool aabbContainsSegment(const float2 start, const float2 finish, const lbvh::aabb2<float> object) noexcept
{
    // Completely outside.
    if ((start.x <= object.lower.x && finish.x <= object.lower.x) || (start.y <= object.lower.y && finish.y <= object.lower.y) ||
        (start.x >= object.upper.x && finish.x >= object.upper.x) || (start.y >= object.upper.y && finish.y >= object.upper.y))
        return false;

    float m = (finish.y - start.y) / (finish.x - start.x);

    float y = m * (object.lower.x - start.x) + start.y;
    if (y > object.lower.y && y < object.upper.y) return true;

    y = m * (object.upper.x - start.x) + start.y;
    if (y > object.lower.y && y < object.upper.y) return true;

    float x = (object.lower.y - start.y) / m + start.x;
    if (x > object.lower.x && x < object.upper.x) return true;

    x = (object.upper.y - start.y) / m + start.x;
    if (x > object.lower.x && x < object.upper.x) return true;

    // Start or end inside.
    if ((start.x > object.lower.x && start.x < object.upper.x && start.y > object.lower.y && start.y < object.lower.y)
        || (finish.x > object.lower.x && finish.x < object.upper.x && finish.y > object.lower.y && finish.y < object.upper.y)) return true;

    return false;
}


__host__ __device__
bool aabbContainsSegment(const double2 start, const double2 finish, const lbvh::aabb2<double> object) noexcept
{
    // Completely outside.
    if ((start.x <= object.lower.x && finish.x <= object.lower.x) || (start.y <= object.lower.y && finish.y <= object.lower.y) ||
        (start.x >= object.upper.x && finish.x >= object.upper.x) || (start.y >= object.upper.y && finish.y >= object.upper.y))
        return false;

    double m = (finish.y - start.y) / (finish.x - start.x);

    double y = m * (object.lower.x - start.x) + start.y;
    if (y > object.lower.y && y < object.upper.y) return true;

    y = m * (object.upper.x - start.x) + start.y;
    if (y > object.lower.y && y < object.upper.y) return true;

    double x = (object.lower.y - start.y) / m + start.x;
    if (x > object.lower.x && x < object.upper.x) return true;

    x = (object.upper.y - start.y) / m + start.x;
    if (x > object.lower.x && x < object.upper.x) return true;

    // Start or end inside.
    if ((start.x > object.lower.x && start.x < object.upper.x && start.y > object.lower.y && start.y < object.lower.y)
        || (finish.x > object.lower.x && finish.x < object.upper.x && finish.y > object.lower.y && finish.y < object.upper.y)) return true;

    return false;
}




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
    double* devRPtr,
    bool rebuildTree
)
{     
    double t1 = omp_get_wtime();
    //double t2 = t1, t3 = t1;

    bool flagRebuild = (ptrTree == nullptr || rebuildTree);

    if (flagRebuild)
    {
        ptrTree.reset();
        ptrNorm.reset();
        ptrPs.reset();

        std::shared_ptr < thrust::device_vector<lbvh::segment<float>>> ps = std::make_shared<thrust::device_vector<lbvh::segment<float>>>(nPanels);
        thrust::transform(
            thrust::make_counting_iterator<unsigned int>(0),
            thrust::make_counting_iterator<unsigned int>(nPanels),
            ps->begin(),
            [devRPtr] __device__(const unsigned int idx) {
            const float x1 = devRPtr[4 * idx + 0];
            const float y1 = devRPtr[4 * idx + 1];
            const float x2 = devRPtr[4 * idx + 2];
            const float y2 = devRPtr[4 * idx + 3];
            return lbvh::segment<float>{ {x1, y1}, { x2, y2 } };
        }
        );
        ptrPs = ps;

        //t2 = omp_get_wtime();

        std::shared_ptr<thrust::device_vector<float2>> norms = std::make_shared<thrust::device_vector<float2>>(nPanels);
        thrust::transform(
            ps->begin(),
            ps->end(),
            norms->begin(),
            []__device__(const lbvh::segment<float>&s) {
            return normal2d(s);
        }
        );

        ptrNorm = norms;


        //t3 = omp_get_wtime();              

        std::shared_ptr<lbvh::bvhS<float, lbvh::segment<float>, aabbS_getter>> bvhS = std::make_shared<lbvh::bvhS<float, lbvh::segment<float>, aabbS_getter>>(ps->begin(), ps->end(), true);
        ptrTree = bvhS;
    }

    //const auto bvh_dev = ptrTree->get_device_reprS();
    const auto bvh_dev = ptrTree->get_device_reprS();
    float2* normptr = thrust::raw_pointer_cast(&((*ptrNorm)[0]));
    lbvh::segment<float>* psptr = thrust::raw_pointer_cast(&((*ptrPs)[0]));
    

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





std::vector<int> lbvh_check_inside_ray(
    int timeStep,
    //const std::vector<Point2D>& vecPositions,
    int nPoints,
    double* devOldpos,
    double* devNewpos,
    int nPanels,
    double* devRPtr,
    bool rebuildTree
)
{
	double t1 = omp_get_wtime();
    double t2 = t1;
    
    bool flag = (ptrTree == nullptr || rebuildTree);

    if (flag)
    {
        ptrTree.reset();

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
            return lbvh::segment<float>{ {x1, y1}, { x2, y2 } };
        }
        );
        lbvh::segment<float>* psptr = thrust::raw_pointer_cast(&ps[0]);

        t2 = omp_get_wtime();

        std::shared_ptr<lbvh::bvhS<float, lbvh::segment<float>, aabbS_getter>> bvhS = std::make_shared<lbvh::bvhS<float, lbvh::segment<float>, aabbS_getter>>(ps.begin(), ps.end(), true);
        ptrTree = bvhS;
    }
    
    //const auto bvh_dev = ptrTree->get_device_reprS();
    const auto bvh_dev = ptrTree->get_device_reprS();

    double t3 = omp_get_wtime();

    thrust::device_vector<lbvh::segment<float>> random_rays(nPoints);
    thrust::transform(
        thrust::make_counting_iterator<unsigned int>(0),
        thrust::make_counting_iterator<unsigned int>(nPoints),
        random_rays.begin(),
        [devOldpos, devNewpos] __device__(const unsigned int idx) {
        const float x1 = devOldpos[3 * idx + 0];
        const float y1 = devOldpos[3 * idx + 1];
        const float x2 = devNewpos[2 * idx + 0];
        const float y2 = devNewpos[2 * idx + 1];
        return lbvh::segment<float>{ {x1, y1}, { x2, y2 } };
    }
    );

    lbvh::segment<float>* ptr = thrust::raw_pointer_cast(&random_rays[0]);

    double t4 = omp_get_wtime();

    /////////////////////////////////////

    thrust::device_vector<int> fitpanel(random_rays.size());
    
    auto numptr = thrust::raw_pointer_cast(&fitpanel[0]);

    thrust::for_each(thrust::device,
        thrust::make_counting_iterator<unsigned int>(0),
        thrust::make_counting_iterator<unsigned int>(nPoints),
        [bvh_dev, ptr, numptr] __device__(const unsigned int idx) {
        const auto self = lbvh::query_nearest2ray<float>{ ptr[idx].p1, ptr[idx].p2 };
        const auto nest = lbvh::query_device2ray(bvh_dev, self, check_hit_obj_Calculator(), check_hit_node_Calculator());
        numptr[idx] = nest;
        return;
    });
    
    thrust::host_vector<int> fitpanelHost = fitpanel;
    
    double t5 = omp_get_wtime();

    std::cout << 1000*(t2 - t1) << " " << 1000 * (t3 - t2) << " " << 1000 * (t4 - t3) << " " << 1000 * (t5 - t4)  << std::endl;

    return std::vector<int>(fitpanelHost.begin(), fitpanelHost.end());
}



