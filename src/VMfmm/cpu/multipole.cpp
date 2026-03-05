#include "multipole.h"
#include "../common/simple_math.h"

#include <algorithm>
#include <numeric>
#include <iostream>
#include "omp.h"
#include <oneapi/tbb/parallel_for_each.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/blocked_range.h>
#include <cmath>
#include <fstream>
#include <ranges>
#include "../common/mpi_utils.h"


namespace fmm {



using std::complex;
using std::vector;
using namespace std::complex_literals;

void FastMultipole::Multipole(std::complex<double>* a, const std::pair<size_t, size_t>& particle_range, const std::complex<double>& z0)
{
	auto p_begin = tree->particles.begin() + particle_range.first;
	auto p_end = tree->particles.begin() + particle_range.second;
	for (int i = 1; i < N; ++i)
	{
		auto temp = std::accumulate(p_begin, p_end, 0.0i, [&z0, i](const complex<double>& x, const particle2d& y)
			{
				return x + y.q * MyPow(y.center - z0, i);
			});
		a[i] = -temp / double(i);
	}

	a[0] = std::accumulate(p_begin, p_end, 0.0i, [](const complex<double>& x, const particle2d& y) { return x + y.q; });
}

void FastMultipole::M2M(const std::complex<double>* a, const complex<double>& z0, std::complex<double>* b, vector<complex<double>>& work)
{
	work[0] = 1.0;
	for (int i = 1; i < N; ++i)
		work[i] = work[i - 1] * z0;

	for (int j = 1; j < N; ++j)
	{
		for (int i = 1; i <= j; ++i)
		{
			b[j] += a[i] * work[j - i] * detail::binom(j - 1, i - 1);
		}
		b[j] -= a[0] * work[j] / double(j);
	}
	b[0] += a[0];
}

inline std::complex<double> divComp(const std::complex<double>& a, const std::complex<double>& b)
{
	const double& ax = a.real();
	const double& bx = b.real();
	const double& ay = a.imag();
	const double& by = b.imag();

	double zn = bx * bx + by * by;
	return { (ax * bx + ay * by) / zn, (ay * bx - ax * by) / zn };
}

void FastMultipole::M2L(const std::complex<double>* a, const complex<double>& z0, std::complex<double>* b, vector<complex<double>>& work)
{
	work[0] = a[0];
	b[0] += a[0] * std::conj(log(-z0));
	for (int i = 1; i < N; ++i)
	{
		work[i] = ni(i) * a[i] / MyPow(z0, i);
		b[0] += work[i];
	}
	for (int i = 1; i < N; ++i)
	{
		complex<double> tmp{ 0.0 };
		for (int j = 1; j < N; ++j)
		{
			tmp += work[j] * detail::binom(i + j - 1, j - 1);
		}
		b[i] += divComp(-a[0] / double(i) + tmp, MyPow(z0, i));
	}
}

void FastMultipole::L2L(const std::complex<double>* a, const complex<double>& z0, std::complex<double>* b)
{
	for (int i = 0; i < N; ++i)
		b[i] = a[i];
	for (int j = 0; j < N; ++j)
	{
		for (int k = N - j - 1; k < N - 1; ++k)
		{
			b[k] -= z0 * b[k + 1];
		}
	}
}

void FastMultipole::Upward()
{
	double t1 = omp_get_wtime();

	auto& leaves = tree->levels.back();
	auto& leaves_outer = outer_expansions.back();
#ifdef FMM_MPI
	auto [Begin, End] = LocalPart(0, tree->level_sizes[tree->tree_depth - 1]);
#else
	size_t Begin = 0, End = tree->level_sizes[tree->tree_depth - 1];
#endif
	tbb::parallel_for(Begin, End, [&](size_t i) {
		TreeCell2d& cell = leaves[i];
		Multipole(leaves_outer.data() + N * i, cell.source_range, cell.center);
	});
#ifdef FMM_MPI
	std::vector<std::complex<double>> buf(leaves_outer.size());
	std::vector<int> sizes(NProc());
	std::vector<int> displs(NProc());
	for (int i = 0; i < NProc(); ++i)
	{
		sizes[i] = N * LocalPart(0, tree->level_sizes[tree->tree_depth - 1], i);
	}multipole time:
	for (int i = 1; i < NProc(); ++i)
	{
		displs[i] = displs[i - 1] + sizes[i - 1];
	}
	MPI_Allgatherv(leaves_outer.data() + displs[MyID()], sizes[MyID()], MPI_COMPLEX16,
		buf.data(), sizes.data(), displs.data(), MPI_COMPLEX16, MPI_COMM_WORLD);
	std::swap(buf, leaves_outer);
#endif
//	if (IAmRoot()) std::cout << "multipole time: " << omp_get_wtime() - t1 << std::endl;

	t1 = omp_get_wtime();

	for (int i = (int)(tree->tree_depth - 2); i >= 0; --i)
	{
		auto& top_level = tree->levels[i];
		auto& top_outer = outer_expansions[i];
		const auto& bottom_level = tree->levels[i + 1];
		const auto& bottom_outer = outer_expansions[i + 1];
#ifdef FMM_MPI
		auto [Begin, End] = LocalPart(0, tree->level_sizes[i]);
#else
		size_t Begin = 0, End = tree->level_sizes[i];
#endif
		tbb::parallel_for(Begin, End, [&](size_t j)
			{
				TreeCell2d& cell = top_level[j];
				auto& work = local_work.local();
				const auto& [begin, end] = cell.source_range;

				for (auto k = begin; k < end; ++k)
				{
					const auto& children_cell = bottom_level[k];
					M2M(bottom_outer.data() + k * N, children_cell.center - cell.center, top_outer.data() + j * N, work);
				}
			});

#ifdef FMM_MPI
		{
			std::vector<std::complex<double>> buf(top_outer.size());
			for (int j = 0; j < NProc(); ++j)
			{
				sizes[j] = N * LocalPart(0, tree->level_sizes[i], j);
			}
			for (int j = 1; j < NProc(); ++j)
			{
				displs[j] = displs[j - 1] + sizes[j - 1];
			}
			MPI_Allgatherv(top_outer.data() + displs[MyID()], sizes[MyID()], MPI_COMPLEX16,
				buf.data(), sizes.data(), displs.data(), MPI_COMPLEX16, MPI_COMM_WORLD);
			std::swap(buf, top_outer);
		}
#endif
	}

	//if (IAmRoot()) std::cout << "m2m time: " << omp_get_wtime() - t1 << std::endl;
}

void FastMultipole::Downward()
{
#ifdef FMM_MPI
	std::vector<int> sizes(NProc());
	std::vector<int> displs(NProc());
#endif

	for (int i = 2; i < tree->tree_depth; ++i)
	{
		const auto& top_level = tree->levels[i - 1];
		const auto& top_inner = inner_expansions[i - 1];
		auto& bottom_level = tree->levels[i];
		auto& bottom_inner = inner_expansions[i];
		const auto& bottom_outer = outer_expansions[i];

#ifdef FMM_MPI
		auto [Begin, End] = LocalPart(0, tree->level_sizes[i]);
#else
		size_t Begin = 0, End = tree->level_sizes[i];
#endif

		tbb::parallel_for(Begin, End, [&](size_t j)
			{
				TreeCell2d& cell = bottom_level[j];
				auto& work = local_work.local();
				auto& parent_cell = top_level[cell.parent];
				L2L(top_inner.data() + N * cell.parent, parent_cell.center - cell.center, bottom_inner.data() + N * j);

				for (const auto& idx : cell.farneighbours)
				{
					const auto& far_neighbour = bottom_level[idx];
					M2L(bottom_outer.data() + idx * N, far_neighbour.center - cell.center, bottom_inner.data() + N * j, work);
				}
			});

#ifdef FMM_MPI
		{
			std::vector<std::complex<double>> buf(bottom_inner.size());
			for (int j = 0; j < NProc(); ++j)
			{
				sizes[j] = N * LocalPart(0, tree->level_sizes[i], j);
			}
			for (int j = 1; j < NProc(); ++j)
			{
				displs[j] = displs[j - 1] + sizes[j - 1];
			}
			MPI_Allgatherv(bottom_inner.data() + displs[MyID()], sizes[MyID()], MPI_COMPLEX16,
				buf.data(), sizes.data(), displs.data(), MPI_COMPLEX16, MPI_COMM_WORLD);
			std::swap(buf, bottom_inner);
		}
#endif
	}
}

void FastMultipole::ComputePotentials()
{
	tbb::enumerable_thread_specific<std::complex<double>> dz_local;
	tbb::enumerable_thread_specific<double> phi_local;
	tbb::enumerable_thread_specific<std::vector<double>> buf_potentials_local((std::vector<double>(num_particles)));

	const auto& leaves = tree->levels.back();
	const auto& leaves_inner = inner_expansions.back();
	auto& particles = tree->particles;

#ifdef FMM_MPI
	auto [Begin, End] = LocalPart(0, tree->level_sizes.back());
#else
	size_t Begin = 0, End = tree->level_sizes.back();
#endif

	decltype(&TreeCell2d::source_range) range_ptr;
	if (tree->targets_num == 0) // targets = sources
		range_ptr = &TreeCell2d::source_range;
	else
		range_ptr = &TreeCell2d::target_range;
	const auto& index_mapping = tree->positions_map;

	tbb::parallel_for(Begin, End, [&](size_t i) {
		auto& dz = dz_local.local();
		auto& phi = phi_local.local();
		const auto& cell = leaves[i];
		const auto& inner = leaves_inner.data() + i * N;
		auto& buf_potentials = buf_potentials_local.local();

		const auto& [p1_begin, p1_end] = cell.*range_ptr;

		for (auto j = p1_begin; j < p1_end; ++j)
		{
			auto& particle = particles[j];

			if (range_ptr == &TreeCell2d::target_range && particle.q == 0)
			{
				phi = 0.0;

				const auto& [p2_begin, p2_end] = cell.source_range;
				for (auto k = p2_begin; k < p2_end; ++k)
					phi += Potential2d(particle, particles[k]);

				for (const auto& idx : cell.closeneighbours)
				{
					const auto& neighbour_cell = leaves[idx];
					const auto& [p2_begin, p2_end] = neighbour_cell.source_range;

					for (auto k = p2_begin; k < p2_end; ++k)
						phi += Potential2d(particle, particles[k]);
				}
			}
			else
			{
				phi = Potential2d(particle, particle);

				for (auto k = j + 1; k < p1_end; ++k)
				{
					Potential2dMutual(particle, particles[k], phi, buf_potentials[k]);
				}

				for (const auto& idx : cell.closeneighbours)
				{
					const auto& neighbour_cell = leaves[idx];
					const auto& [p2_begin, p2_end] = neighbour_cell.source_range;

					for (auto k = p2_begin; k < p2_end; ++k)
						Potential2dMutual(particle, particles[k], phi, buf_potentials[k]);
				}
			}
			buf_potentials[j] += phi;
		}
		
		for (auto j = p1_begin; j < p1_end; ++j)
		{
			auto& particle = particles[j];
			phi = 0.0;
			dz = particle.center - cell.center;
			complex<double> dzpow{ 1.0 };
			for (int k = 0; k < N; ++k)
			{
				phi += inner[k].real() * dzpow.real() - inner[k].imag() * dzpow.imag();
				dzpow *= dz;
			}
			buf_potentials[j] += phi;
		}
	});

	auto& buf_potentials = *buf_potentials_local.begin();
	for (auto it = buf_potentials_local.begin() + 1; it < buf_potentials_local.end(); it++)
	{
		const auto& x = *it;
		tbb::parallel_for(tbb::blocked_range<size_t>(0, num_particles),
			[&](tbb::blocked_range<size_t> r) {
				for (auto i = r.begin(); i < r.end(); ++i)
					buf_potentials[i] += x[i];
			});
	}

#ifdef FMM_MPI
	AllReduce(buf_potentials.data(), buf_potentials.size());
#endif

	tbb::parallel_for(size_t(0), targets_num, [&](size_t i){
		potentials[i] = buf_potentials[index_mapping[i]];
	});
}

void FastMultipole::ComputeForces()
{
	tbb::enumerable_thread_specific<std::complex<double>> dz_local;
	tbb::enumerable_thread_specific<std::complex<double>> force_local;
	tbb::enumerable_thread_specific<std::vector<std::complex<double>>> buf_forces_local((std::vector<std::complex<double>>(num_particles)));

	const auto& leaves = tree->levels.back();
	const auto& leaves_inner = inner_expansions.back();
	auto& particles = tree->particles;

#ifdef FMM_MPI
	auto [Begin, End] = LocalPart(0, tree->level_sizes.back());
#else
	size_t Begin = 0, End = tree->level_sizes.back();
#endif

	decltype(&TreeCell2d::source_range) range_ptr;
	if (tree->targets_num == 0) // targets = sources
		range_ptr = &TreeCell2d::source_range;
	else
		range_ptr = &TreeCell2d::target_range;
	const auto& index_mapping = tree->positions_map;

	tbb::parallel_for(Begin, End, [&](size_t i) {
		auto& dz = dz_local.local();
		auto& force = force_local.local();
		const auto& cell = leaves[i];
		const auto& inner = leaves_inner.data() + i * N;
		auto& buf_forces = buf_forces_local.local();

		const auto& [p1_begin, p1_end] = cell.*range_ptr;

		for (auto j = p1_begin; j < p1_end; ++j)
		{
			auto& particle = particles[j];

			if (range_ptr == &TreeCell2d::target_range && particle.q == 0)
			{
				force = 0.0;

				const auto& [p2_begin, p2_end] = cell.source_range;
				for (auto k = p2_begin; k < p2_end; ++k)
					force += Force2d(particle, particles[k]);

				for (const auto& idx : cell.closeneighbours)
				{
					const auto& neighbour_cell = leaves[idx];
					const auto& [p2_begin, p2_end] = neighbour_cell.source_range;

					for (auto k = p2_begin; k < p2_end; ++k)
						force += Force2d(particle, particles[k]);
				}
			}
			else
			{
				force = Force2d(particle, particle);

				for (auto k = j + 1; k < p1_end; ++k)
				{
					Force2dMutual(particle, particles[k], force, buf_forces[k]);
				}

				for (const auto& idx : cell.closeneighbours)
				{
					const auto& neighbour_cell = leaves[idx];
					const auto& [p2_begin, p2_end] = neighbour_cell.source_range;

					for (auto k = p2_begin; k < p2_end; ++k)
						Force2dMutual(particle, particles[k], force, buf_forces[k]);
				}
			}
			buf_forces[j] += force;
		}

		for (auto j = p1_begin; j < p1_end; ++j)
		{
			auto& particle = particles[j];
			force = 0.0;
			dz = particle.center - cell.center;
			complex<double> dzpow{ 1.0 };
			for (int k = 1; k < N; ++k)
			{
				force += double(k) * inner[k] * dzpow;
				dzpow *= dz;
			}
			buf_forces[j] += force;
		}
	});

	auto& buf_forces = *buf_forces_local.begin();
	for (auto it = buf_forces_local.begin() + 1; it < buf_forces_local.end(); it++)
	{
		const auto& x = *it;
		tbb::parallel_for(tbb::blocked_range<size_t>(0, num_particles),
			[&](tbb::blocked_range<size_t> r) {
				for (auto i = r.begin(); i < r.end(); ++i)
					buf_forces[i] += x[i];
			});
	}

#ifdef FMM_MPI
	AllReduce(buf_forces.data(), buf_forces.size());
#endif

	tbb::parallel_for(size_t(0), targets_num, [&](size_t i){
		forces[i] = buf_forces[index_mapping[i]];
	});
}

FastMultipole::FastMultipole(const std::vector<particle2d>& particles, double eps, int N_, int tree_depth_)
{
	auto particles_copy = particles;
	targets_num = particles.size();
	Solve(std::move(particles_copy), eps, N_, tree_depth_);
}

FastMultipole::FastMultipole(std::vector<particle2d>&& particles, double eps, int N_, int tree_depth_)
{
	targets_num = particles.size();
	Solve(std::move(particles), eps, N_, tree_depth_);
}

FastMultipole::FastMultipole(const std::vector<particle2d>& source_particles, const std::vector<particle2d>& target_particles, double eps, int N_, int tree_depth_)
{
	std::vector<particle2d> particles(source_particles.size() + target_particles.size());
	memcpy(particles.data(), target_particles.data(), target_particles.size() * sizeof(particle2d));
	memcpy(particles.data() + target_particles.size(), source_particles.data(), source_particles.size() * sizeof(particle2d));
	targets_num = target_particles.size();
	Solve(std::move(particles), eps, N_, tree_depth_);
}

void FastMultipole::Solve(std::vector<particle2d>&& particles, double eps, int N_, int tree_depth_)
{
	double T = omp_get_wtime();
	//if (IAmRoot()) std::cout << "\n***************** Start FMM **************" << std::endl;
	num_particles = particles.size();
	if (N_ == FMM_AUTO)
		N = std::min(fmm::detail::_2d_MAX_MULTIPOLE_NUM, int(-1.18 * log(1.65 * eps)));
	else
		N = std::min(fmm::detail::_2d_MAX_MULTIPOLE_NUM, N_);
	//if (IAmRoot()) std::cout << "multipole num = " << N << std::endl;

	size_t initial_depth;
	if (tree_depth_ == FMM_AUTO)
		initial_depth = std::max((size_t)4, size_t(log(double(num_particles)) / log(4.0) + 0.5)) - 1;
		///size_t(std::max(3.0, ))//7;// (size_t)std::max(3.0, 2 + log(double(num_particles) / 100.0) / log(4.0));
	else
		initial_depth = tree_depth_;
	//if (IAmRoot()) std::cout << "tree depth = " << initial_depth << std::endl;
	
	tree = std::make_shared<MortonTree2d>(std::move(particles), initial_depth);
	tree_depth = tree->tree_depth;

	//if (tree_depth < initial_depth)
	//	std::cout << "Warning: tree_depth was reduced from " << initial_depth << " to " << tree_depth << " to ensure H_CELL < EPS_FORCE" << std::endl;


	double t = omp_get_wtime();
	outer_expansions.resize(tree_depth);
	inner_expansions.resize(tree_depth);
	for (int i = 0; i < tree->tree_depth; ++i)
	{
		auto& outer = outer_expansions[i];
		auto& inner = inner_expansions[i];
		outer.resize(tree->level_sizes[i] * N);
		inner.resize(tree->level_sizes[i] * N);
	}
	forces.resize(targets_num);
	potentials.resize(targets_num);
	//if (IAmRoot()) std::cout << "allocation time: " << omp_get_wtime() - t << std::endl;
	local_work = tbb::enumerable_thread_specific<vector<complex<double>>>{ vector<complex<double>>(N) };

	Upward();
	
	t = omp_get_wtime();
	Downward();
	//if (IAmRoot()) std::cout << "m2l+l2l time: " << omp_get_wtime() - t << std::endl;

	t = omp_get_wtime();
	ComputePotentials();
	//if (IAmRoot()) std::cout << "leaf potential time: " << omp_get_wtime() - t << std::endl;

	t = omp_get_wtime();
	ComputeForces();
	//if (IAmRoot()) std::cout << "leaf force time: " << omp_get_wtime() - t << std::endl;

	//if (IAmRoot()) std::cout << "total fmm time: " << omp_get_wtime() - T << std::endl;
	if (IAmRoot()) std::cout << "FMM time ConvVelo Wake: " << (omp_get_wtime() - T) * 1000.0 << std::endl;
	//if (IAmRoot()) std::cout << "***************** End FMM **************\n" << std::endl;
}


} // fmm