#ifndef OPTIMIZED_VELOCITY_2D_H
#define OPTIMIZED_VELOCITY_2D_H

#include <memory>

#include "defs.h"
#include "cudaTreeInfo.h"


class OptimizedVelocity
{
private:
    std::vector<std::vector<Point2D>> allStepsVelocities;
    std::vector<std::vector<double>> allStepsPressures;
    bool isAveraged = true;

public:
    void AddVelocities(const std::vector<Point2D>& velocities);
    std::vector<Point2D> GetAverageVelocity() const;
    Point2D GetGlobalAverageVelocity(const std::vector<Point2D>& averageVelocities) const;
    void PerformAveragingVelo(const std::string& path);

    /*void AddPressures(const std::vector<double>& pressures);
    std::vector<double> GetAveragePressure() const;
    double GetGlobalAveragePressure(const std::vector<double>& averagePressures) const;\
    void PerformAveragingPress();*/
};
#endif