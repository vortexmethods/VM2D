#ifndef KNNCPU_H
#define KNNCPU_H

#include "Vortex2D.h"

void WakekNN(const std::vector<Vortex2D>& vtx, const size_t k, std::vector<std::vector<std::pair<double, size_t>>>& initdist,
	double cSP, double cRBP, double maxG, double epsCol, int type);


#endif