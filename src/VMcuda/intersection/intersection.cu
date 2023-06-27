#include "lbvh.cuh"
#include <random>
#include <vector>
#include <thrust/random.h>
#include <cmath>

#include <omp.h>

#include "Point2D.h"

struct segment
{
    double2 beg;
    double2 end;
};

struct aabb_getter
{
    __device__
        lbvh::aabb<double> operator()(const segment f) const noexcept
    {
        lbvh::aabb<double> retval;
        retval.upper.x = ::fmax(f.beg.x, f.end.x);
        retval.upper.y = ::fmax(f.beg.y, f.end.y);
        retval.lower.x = ::fmin(f.beg.x, f.end.x);
        retval.lower.y = ::fmin(f.beg.y, f.end.y);
        return retval;
    }
};

struct distance_calculator
{
    __device__
        thrust::pair<double,double> operator()(const double2 point, const segment object) const noexcept
    {
        const double L2 = (object.end.x - object.beg.x) * (object.end.x - object.beg.x) \
            + (object.end.y - object.beg.y) * (object.end.y - object.beg.y);
        if (L2 == 0)
            return thrust::make_pair( (object.end.x - point.x) * (object.end.x - point.x) \
            + (object.end.y - point.y) * (object.end.y - point.y), 0.0 );

        const double t = ::fmax(0.0, ::fmin(1.0, ((point.x - object.beg.x) * (object.end.x - object.beg.x) + (point.y - object.beg.y) * (object.end.y - object.beg.y)) / L2));

        const double2 projection = { object.beg.x + t * (object.end.x - object.beg.x), \
            object.beg.y + t * (object.end.y - object.beg.y) };
        
        return thrust::make_pair( (point.x - projection.x) * (point.x - projection.x) +
            (point.y - projection.y) * (point.y - projection.y), t );
    }
};

std::vector<unsigned int> lbvh_check_inside(
    int timeStep,
    //const std::vector<Point2D>& vecPositions,
    int nPoints,
    double* devNewpos,
    int nPanels,
    double* devRPtr
)
{    
    //auto tm1 = omp_get_wtime();

    thrust::device_vector<segment> ps(nPanels);
    thrust::transform(
       thrust::make_counting_iterator<unsigned int>(0),
       thrust::make_counting_iterator<unsigned int>(nPanels),
       ps.begin(),
       [devRPtr] __device__(const unsigned int idx) {
        const double x1 = devRPtr[4 * idx + 0];
        const double y1 = devRPtr[4 * idx + 1];
        const double x2 = devRPtr[4 * idx + 2];
        const double y2 = devRPtr[4 * idx + 3];
        return segment{ {x1, y1}, {x2, y2} };
    }
    );
   
    //auto tm2 = omp_get_wtime();

    //lbvh::bvh<double, segment, aabb_getter> bvh(psHost.begin(), psHost.end(), true);
    lbvh::bvh<double, segment, aabb_getter> bvh(ps);

    const auto bvh_dev = bvh.get_device_repr();
    
    //auto tm3 = omp_get_wtime();

    thrust::device_vector<double2> random_points(nPoints);
    thrust::transform(
        thrust::make_counting_iterator<unsigned int>(0),
        thrust::make_counting_iterator<unsigned int>(nPoints),
        random_points.begin(),
        [devNewpos] __device__(const unsigned int idx) {
        const double x = devNewpos[2 * idx + 0];
        const double y = devNewpos[2 * idx + 1];
        return make_double2(x, y);
    }
    );

    //auto tm4 = omp_get_wtime();

    thrust::device_vector<unsigned int> fitpanel(nPoints);
    double2* ptr_random_points = thrust::raw_pointer_cast(&random_points[0]);
    thrust::transform(
        thrust::make_counting_iterator<unsigned int>(0),
        thrust::make_counting_iterator<unsigned int>(nPoints),
        fitpanel.begin(),        
        [bvh_dev, nPanels, devRPtr, ptr_random_points, timeStep] __device__(const unsigned int idx) {
        const double2 pos = *(ptr_random_points + idx);
        const auto calc = distance_calculator();

        //if ((timeStep >= 1401) && (idx == 29441))
        //    printf("pos = {%f, %f}\n,", pos.x, pos.y);

        const auto nest = lbvh::query_device((int)idx, (int)timeStep, bvh_dev, lbvh::nearest(pos), calc);

        const unsigned int pnl = nest.first.first;
        const double frac = nest.first.second;

        if (pnl >= nPanels)
            return (unsigned int)(-1);

        //if (pnl < 0 || pnl >= nPanels)
        //    printf("QA! idx = %d, nPanels = %d, pnl = %d, frac = %f\n", (int)idx, (int)nPanels, (int)pnl, (double)frac);
        //
        //if (frac > 1.0 || frac < 0.0)
        //    printf("QB! idx = %d, nPanels = %d, pnl = %d, frac = %f\n", (int)idx, (int)nPanels, (int)pnl, (double)frac);

        const double2 ptBeg = { devRPtr[4 * pnl + 0], devRPtr[4 * pnl + 1] };
        const double2 ptEnd = { devRPtr[4 * pnl + 2], devRPtr[4 * pnl + 3] };
        const double2 ptMid = { ptBeg.x + frac * (ptEnd.x - ptBeg.x), ptBeg.y + frac * (ptEnd.y - ptBeg.y) };

        const double2 vec = { ptMid.x - pos.x, ptMid.y - pos.y };

        double2 nrmBase = { ptEnd.y - ptBeg.y, -(ptEnd.x - ptBeg.x) };

        if (frac == 0.0)
        {
            unsigned int pnlPrev = (pnl > 0) ? pnl - 1 : nPanels - 1;
            const double2 ptPrev = { devRPtr[4 * pnlPrev + 0], devRPtr[4 * pnlPrev + 1] };

            const double2 nrmAdd = { ptBeg.y - ptPrev.y, -(ptBeg.x - ptPrev.x)};
            nrmBase.x += nrmAdd.x;
            nrmBase.y += nrmAdd.y;
        }

        if (frac == 1.0)
        {
            unsigned int pnlNext = (pnl + 1 < nPanels) ? pnl + 1 : 0;
            const double2 ptNext = { devRPtr[4 * pnlNext + 2], devRPtr[4 * pnlNext + 3] };

            const double2 nrmAdd = { ptNext.y - ptEnd.y, -(ptNext.x - ptEnd.x) };
            nrmBase.x += nrmAdd.x;
            nrmBase.y += nrmAdd.y;
        }

        if (vec.x * nrmBase.x + vec.y * nrmBase.y > 0)
            return (unsigned int)pnl;
        else 
            return (unsigned int)(-1);
    }
        );

    thrust::host_vector<unsigned int> fitpanelHost = fitpanel;

    //auto tm5 = omp_get_wtime();
        
    //printf("TimeInt = %f, %f, %f, %f\n", tm2 - tm1, tm3 - tm2, tm4 - tm3, tm5 - tm4);
     
    return std::vector<unsigned int>(fitpanelHost.begin(), fitpanelHost.end());
}
