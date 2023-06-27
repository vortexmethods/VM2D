#ifndef KNN_CUH
#define KNN_CUH

#include <vector>
#include "Vortex2D.h"


class VectorsForKnn;

void kNNcuda(const std::vector<Vortex2D>& vtx,
	const size_t k,
	std::vector<std::pair<double, size_t>>& initdist,
	VectorsForKnn*& vecForKnn,
	double cSP, double cRBP, double maxG, double epsCol, int typ);


#endif //KNN_CUH