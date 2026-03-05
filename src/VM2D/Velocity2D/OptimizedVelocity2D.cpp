#include "OptimizedVelocity2D.h"

#include "Airfoil2D.h"
#include "Boundary2D.h"
#include "MeasureVP2D.h"
#include "Mechanics2D.h"
#include "StreamParser.h"
#include "Wake2D.h"
#include "World2D.h"

#include "BarnesHut.h"

using namespace VM2D;

std::vector<Point2D> OptimizedVelocity::GetAverageVelocity() const
{
	size_t numSteps = allStepsVelocities.size();
	size_t numPoints = allStepsVelocities[0].size();

	std::vector<Point2D> averageVelocities(numPoints, { 0.0, 0.0 });

	for (size_t step = 0; step < numSteps; ++step) {
		for (size_t i = 0; i < numPoints; ++i) {
			averageVelocities[i][0] += allStepsVelocities[step][i][0];
			averageVelocities[i][1] += allStepsVelocities[step][i][1];
		}
	}

	for (size_t i = 0; i < numPoints; ++i) {
		averageVelocities[i][0] /= numSteps;
		averageVelocities[i][1] /= numSteps;
	}

	return averageVelocities;
}

Point2D OptimizedVelocity::GetGlobalAverageVelocity(const std::vector<Point2D>& averageVelocities) const
{
	double sumX = 0.0;
	double sumY = 0.0;
	size_t numPoints = averageVelocities.size();

	for (size_t i = 0; i < numPoints; ++i) {
		sumX += averageVelocities[i][0];
		sumY += averageVelocities[i][1];
	}

	Point2D globalAverage;
	globalAverage[0] = sumX / numPoints;
	globalAverage[1] = sumY / numPoints;

	return globalAverage;
}

void OptimizedVelocity::AddVelocities(const std::vector<Point2D>& velocities)
{
	allStepsVelocities.push_back(velocities);
}

void OptimizedVelocity::PerformAveragingVelo(const std::string& path)
{
	if (isAveraged) {
		std::vector<Point2D> averageVelocities = GetAverageVelocity();
		Point2D globalAverage = GetGlobalAverageVelocity(averageVelocities);
		std::cout << "Velocity value: (" << globalAverage[0] << ", " << globalAverage[1] << ")\n";
		std::ofstream psiFile(path + "psi.txt");
		psiFile << "Velocity value: (" << globalAverage[0] << ", " << globalAverage[1] << ")\n";
		psiFile.close();

		isAveraged = false;
	}
}

//void OptimizedVelocity::AddPressures(const std::vector<double>& pressures)
//{
//	allStepsPressures.push_back(pressures);
//}
//
//std::vector<double> OptimizedVelocity::GetAveragePressure() const
//{
//	size_t numSteps = allStepsPressures.size();
//	size_t numPoints = allStepsPressures[0].size();
//
//	std::vector<double> averagePressures(numPoints, 0.0);
//
//	for (size_t step = 0; step < numSteps; ++step) {
//		for (size_t i = 0; i < numPoints; ++i) {
//			averagePressures[i] += allStepsPressures[step][i];
//		}
//	}
//
//	for (size_t i = 0; i < numPoints; ++i) {
//		averagePressures[i] /= numSteps;
//	}
//
//	return averagePressures;
//}
//
//double OptimizedVelocity::GetGlobalAveragePressure(const std::vector<double>& averagePressures) const
//{
//	double sumPressure = 0.0;
//	size_t numPoints = averagePressures.size();
//
//	for (size_t i = 0; i < numPoints; ++i) {
//		sumPressure += averagePressures[i];
//	}
//
//	return sumPressure / numPoints;
//}
//
//void OptimizedVelocity::PerformAveragingPress()
//{
//	if (isAveraged) {
//		std::vector<double> averagePressure = GetAveragePressure();
//		double globalAverage = GetGlobalAveragePressure(averagePressure);
//		std::cout << "Значение давления: (" << globalAverage << ")\n";
//		isAveraged = false;
//	}
//}