#pragma once

#include "morton_tree.h"
#include <oneapi/tbb/enumerable_thread_specific.h>
#include <oneapi/tbb/concurrent_unordered_map.h>
#include "../common/utils.h"

namespace fmm {

template <typename value>
class FastMultipole3d
{
public:
	using complex_value = std::conditional_t<std::is_same_v<value, double>, std::complex<double>, Vector3cd>;
	using particle_t = particle<Vector3d, value>;
	using MortonTree_t = MortonTree<3, point3d, value>;
	using TreeCell_t = typename MortonTree_t::TreeCell_t;

	FastMultipole3d(const std::vector<particle_t>& particles, double eps = 1.e-8, int N = FMM_AUTO, int tree_depth = FMM_AUTO);
	FastMultipole3d(std::vector<particle_t>&& particles, double eps = 1.e-8, int N = FMM_AUTO, int tree_depth = FMM_AUTO);
	FastMultipole3d(const std::vector<particle_t>& source_particles, const std::vector<particle_t>& target_particles, double eps = 1.e-8, int N = FMM_AUTO, int tree_depth = FMM_AUTO);

	std::vector<Vector3d> forces;
	std::vector<Vector3d> moments;
	std::vector<double> potentials;
	std::vector<Vector3d> dipolevelocities;
private:
	void Solve(std::vector<particle_t>&& particles, double eps = 1.e-8, int N = FMM_AUTO, int tree_depth = FMM_AUTO);
	void Multipole(complex_value* a, const std::pair<size_t, size_t>& particle_range, const Vector3d& z0, double* Pnm, std::complex<double>* eim);
	void M2M(const complex_value* a, const Vector3d& p0, complex_value* b, std::vector<complex_value>& work1, std::vector<complex_value>& work2);
	void M2L(const complex_value* a, const Vector3d& p0, complex_value* b, std::vector<complex_value>& work1, std::vector<complex_value>& work2);
	void L2L(const complex_value* a, const Vector3d& p0, complex_value* b, std::vector<complex_value>& work1, std::vector<complex_value>& work2);

	void RotateY(const complex_value* a, const std::vector<double>& dmatrix, complex_value* res);
	void RotateZ(const complex_value* a, const std::vector<std::complex<double>>& rotation_exp, complex_value* b);
	void RotateZ(complex_value* a, const std::vector<std::complex<double>>& rotation_exp);

	void Upward();
	void Downward();
	void ComputeLeaves();


	void Prepare();

	std::shared_ptr<MortonTree_t> tree;
	std::vector<std::vector<complex_value>> outer_expansions; // for all tree levels
	std::vector<std::vector<complex_value>> inner_expansions; // for all tree levels

	size_t tree_depth;
	int N; // multipole_num
	int Nx2;
	size_t num_particles;
	size_t targets_num;
	//std::vector<double> m2lcoef;

	tbb::enumerable_thread_specific<std::vector<complex_value>> local_work1;
	tbb::enumerable_thread_specific<std::vector<complex_value>> local_work2;

};

}