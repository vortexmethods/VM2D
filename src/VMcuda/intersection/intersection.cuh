#ifndef INTERSECTION_CUH
#define INTERSECTION_CUH

#include <vector>
#include "Point2D.h"

std::vector<int> lbvh_check_inside(
    int timeStep,
    //const std::vector<Point2D>& vecPositions,
    int nPoints,
    double* devNewpos,
    int nPanels,
    double* devRPtr
);

#endif //INTERSECTION_CUH