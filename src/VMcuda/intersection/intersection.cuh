#ifndef INTERSECTION_CUH
#define INTERSECTION_CUH

#include <vector>
#include "Point2D.h"
#include <memory>



struct aabbS_getter;

class lbvh_find_nearest_panel
{
public:
    lbvh_find_nearest_panel();
    ~lbvh_find_nearest_panel();

    std::shared_ptr<void> baseTreeVoid;

    //std::shared_ptr<lbvh::bvhS<double, lbvh::segment<double>, aabbS_getter>> baseTree;
    //std::shared_ptr<thrust::device_vector<lbvh::segment<double>>> ptrPs;

    
    std::vector<std::pair<int, double>> get_nearest_panels(
        int timeStep,
        int nPoints,
        double* devNewpos, //«десь это вихри, т.е. (x,y,G)
        int nPanels,
        double* devRPtr,
        bool rebuildTree
    );
    
};


class lbvh_penetration_control
{
public:
    //std::shared_ptr<lbvh::bvhS<double, lbvh::segment<double>, aabbS_getter>> ptrTree(nullptr);
    //std::shared_ptr<thrust::device_vector<double2>> ptrNorm(nullptr);
    //std::shared_ptr<thrust::device_vector<lbvh::segment<double>>> ptrPs(nullptr);
    std::shared_ptr<void> ptrTreeVoid;
    std::shared_ptr<void> ptrNormVoid;
    std::shared_ptr<void> ptrPsVoid;

    lbvh_penetration_control();
    ~lbvh_penetration_control();

    std::vector<int> lbvh_check_inside(
        int timeStep,
        int nPoints,
        double* devNewpos,
        int nPanels,
        double* devRPtr,
        bool rebuildTree
    );


    std::vector<int> lbvh_check_inside_ray(
        int timeStep,
        int nPoints,
        double* devOldpos,
        double* devNewpos,
        int nPanels,
        double* devRPtr,
        bool rebuildTree
    );
};




#endif //INTERSECTION_CUH