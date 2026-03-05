#pragma once

#include "morton_tree.h"
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <ranges>
#include "../common/simple_math.h"

namespace fmm {

class FastMultipole
{	
public:
	FastMultipole(const std::vector<particle2d>& particles, double eps = 1.e-8, int N = FMM_AUTO, int tree_depth = FMM_AUTO);
	FastMultipole(std::vector<particle2d>&& particles, double eps = 1.e-8, int N = FMM_AUTO, int tree_depth = FMM_AUTO);
	FastMultipole(const std::vector<particle2d>& source_particles, const std::vector<particle2d>& target_particles, double eps = 1.e-8, int N = FMM_AUTO, int tree_depth = FMM_AUTO);

	std::vector<std::complex<double>> forces;
	std::vector<double> potentials;
private:
	void Solve(std::vector<particle2d>&& particles, double eps = 1.e-8, int N = FMM_AUTO, int tree_depth = FMM_AUTO);
	void Multipole(std::complex<double>* a, const std::pair<size_t, size_t>& particle_range, const std::complex<double>& z0);
	void M2M(const std::complex<double>* a, const std::complex<double>& z0, std::complex<double>* b, std::vector<std::complex<double>>& work);
	void M2L(const std::complex<double>* a, const std::complex<double>& z0, std::complex<double>* b, std::vector<std::complex<double>>& work);
	void L2L(const std::complex<double>* a, const std::complex<double>& z0, std::complex<double>* b);

	void Upward();
	void Downward();

	void ComputePotentials();
	void ComputeForces();

	std::shared_ptr<MortonTree2d> tree;
	std::vector<std::vector<std::complex<double>>> outer_expansions; // for all tree levels
	std::vector<std::vector<std::complex<double>>> inner_expansions; // for all tree levels
		
	size_t tree_depth;
	int N; // multipole_num
	size_t num_particles;
	size_t targets_num;

	tbb::enumerable_thread_specific<std::vector<std::complex<double>>> local_work;	
};

} // fmm