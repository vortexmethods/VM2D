#include "morton_tree.h"
#include "../common/simple_math.h"

#include <numeric>
#include <algorithm>
#include <execution>
#include <cmath>
#include "omp.h"
#include <iostream>
#include <oneapi/tbb/parallel_sort.h>
#include <oneapi/tbb/parallel_for_each.h>
#include <oneapi/tbb/scalable_allocator.h>
#include <oneapi/tbb/concurrent_vector.h>
#include <oneapi/tbb/parallel_reduce.h>
#include <atomic>
#include <cassert>
//#include <boost/sort/sort.hpp>
#include <iterator>
#include <oneapi/tbb/enumerable_thread_specific.h>
#include "../common/mpi_utils.h"

namespace fmm {

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

#define exec_policy std::execution::par

template <typename point, typename value>
std::pair<point, double> box_extent(const std::vector<particle<point, value>>& particles);

template <int dim>
inline size_t ExpandBits(double x, int tree_depth);

template <int dim, typename point>
inline size_t MortonCode(point val, int tree_depth, point shift, double scale);

template <typename point>
point PointFromMortonCode(size_t val, int depth, point shift, double scale);

//#define deterministic_behaviour

template <>
inline size_t ExpandBits<2>(double x, int tree_depth)
{
	assert(x > 0.0 && x < 1.0);
#ifdef deterministic_behaviour
	size_t bit_shift = (1ull << 30);
#else
	size_t bit_shift = (1ull << tree_depth);
#endif
	size_t val = (size_t)(bit_shift * x);
	val = ((val << 16) | val) & 0xffff0000ffff;
	val = ((val << 8) | val) & 0xff00ff00ff00ff;
	val = ((val << 4) | val) & 0xf0f0f0f0f0f0f0f;
	val = ((val << 2) | val) & 0x3333333333333333;
	val = ((val << 1) | val) & 0x5555555555555555;
	return val;
}

template <>
inline size_t ExpandBits<3>(double x, int tree_depth)
{
	assert(x > 0.0 && x < 1.0);
#ifdef deterministic_behaviour
	size_t bit_shift = (1ull << 20);
#else
	size_t bit_shift = (1ull << tree_depth);
#endif
	size_t val = (size_t)(bit_shift * x);
	val = ((val << 32) | val) & 0xffff00000000ffff;
    val = ((val << 16) | val) & 0xff0000ff0000ff;
    val = ((val << 8) | val) & 0xf00f00f00f00f00f;
    val = ((val << 4) | val) & 0x30c30c30c30c30c3;
    val = ((val << 2) | val) & 0x9249249249249249;
	return val;
}

template <>
inline size_t MortonCode<2, point2d>(point2d val, int tree_depth, point2d shift, double scale)
{
	double x = (val.real() - shift.real()) / scale;
	double y = (val.imag() - shift.imag()) / scale;
	return (ExpandBits<2>(x, tree_depth) << 1) | ExpandBits<2>(y, tree_depth);
}

template <>
inline size_t MortonCode<3, point3d>(point3d val, int tree_depth, point3d shift, double scale)
{
	double x = (val[0] - shift[0]) / scale;
	double y = (val[1] - shift[1])  / scale;
	double z = (val[2] - shift[2])  / scale;
	return (ExpandBits<3>(x, tree_depth) << 2) | (ExpandBits<3>(y, tree_depth) << 1) | ExpandBits<3>(z, tree_depth);
}

template <>
point2d PointFromMortonCode(size_t val, int depth, point2d shift, double scale)
{
	double step = 1.0 / (1 << (depth + 1));
	double x = step, y = step;
	step *= 2;
	for (int i = 0; i < depth; ++i)
	{
		y += step * (val & 1);
		val >>= 1;
		x += step * (val & 1);
		val >>= 1;
		step *= 2;
	}
	x *= scale;
	x += shift.real();
	y *= scale;
	y += shift.imag();
	return { x, y };
}

template <>
point3d PointFromMortonCode(size_t val, int depth, point3d shift, double scale)
{
	double step = 1.0 / (1 << (depth + 1));
	double x = step, y = step, z = step;
	step *= 2;
	for (int i = 0; i < depth; ++i)
	{
		z += step * (val & 1);
		val >>= 1;
		y += step * (val & 1);
		val >>= 1;
		x += step * (val & 1);
		val >>= 1;
		step *= 2;
	}
	x *= scale;
	x += shift[0];
	y *= scale;
	y += shift[1];
	z *= scale;
	z += shift[2];
	return { x, y, z };
}

template <int dim, typename point, typename value>
void MortonTree<dim, point, value>::FillNeighbours()
{
	if (tree_depth > 1)
	for (int i = 0; i < level_sizes[1]; ++i)
	{
		auto& cell = levels[1][i];
		auto& [s, e] = levels[0][0].source_range;
		std::vector<size_t> indexes(e - s);
		std::iota(indexes.begin(), indexes.end(), s);
		indexes.erase(indexes.begin() + i);
		cell.closeneighbours = indexes;
	}
	for (int i = 2; i < tree_depth; ++i)
	{
		double dw = sqrt(dim + 1.e-5) * scale / (1 << i);
		dw *= dw;
		const auto& parent_level = levels[i - 1];
		auto& current_level = levels[i];

		tbb::parallel_for_each(current_level.begin(), current_level.begin() + level_sizes[i], [&](auto& cell)
			{
				if constexpr (dim == 2)
				{
					cell.closeneighbours.reserve(9);
					cell.farneighbours.reserve(27);
				}
				else
				{
					cell.closeneighbours.reserve(27);
					cell.farneighbours.reserve(189);
				}

				const auto& parent_cell = parent_level[cell.parent];
				for (const auto& parent_neighbours : parent_cell.closeneighbours)
				{
					const auto& [s, e] = parent_level[parent_neighbours].source_range;
					for (auto k = s; k < e; ++k)
					{
						if (norm(current_level[k].center - cell.center) < dw)
							cell.closeneighbours.push_back(k);
						else
							cell.farneighbours.push_back(k);
					}
				}

				const auto& [s, e] = parent_level[cell.parent].source_range;
				for (auto k = s; k < e; ++k)
				{
					if (&cell != &(current_level[k]))
					{
						if (norm(current_level[k].center - cell.center) < dw)
							cell.closeneighbours.push_back(k);
						else
							cell.farneighbours.push_back(k);
					}
				}
				cell.farneighbours.shrink_to_fit();
				cell.closeneighbours.shrink_to_fit();
			});
	}
}

template <int dim, typename value>
struct small_particle
{
	double coords[dim];
	value q;
	size_t morton_code;
};

template<typename value>
std::pair<point2d, double> box_extent(const std::vector<particle<point2d, value>>& particles)
{
 	tbb::enumerable_thread_specific<double> xmin(particles[0].center.real()), xmax(particles[0].center.real()),
											ymin(particles[0].center.imag()), ymax(particles[0].center.imag());
	tbb::parallel_for(tbb::blocked_range<size_t>(0, particles.size()), [&](const tbb::blocked_range<size_t>& r)
		{
			auto& Xmin = xmin.local();
			auto& Xmax = xmax.local();
			auto& Ymin = ymin.local();
			auto& Ymax = ymax.local();
			for (size_t i = r.begin(); i < r.end(); ++i)
			{
				double x = particles[i].center.real();
				double y = particles[i].center.imag();
				Xmin = std::min(Xmin, x);
				Xmax = std::max(Xmax, x);
				Ymin = std::min(Ymin, y);
				Ymax = std::max(Ymax, y);
			}
		});
	double Xmin = *std::min_element(xmin.begin(), xmin.end());
	double Ymin = *std::min_element(ymin.begin(), ymin.end());
	double Xmax = *std::max_element(xmax.begin(), xmax.end());
	double Ymax = *std::max_element(ymax.begin(), ymax.end());

	double dw = std::max(Xmax - Xmin, Ymax - Ymin);
	return { {Xmin,Ymin}, dw };
}

template<typename value>
std::pair<point3d, double> box_extent(const std::vector<particle<point3d, value>>& particles)
{
 	tbb::enumerable_thread_specific<double> xmin(particles[0].center[0]), xmax(particles[0].center[0]),
											ymin(particles[0].center[1]), ymax(particles[0].center[1]),
											zmin(particles[0].center[2]), zmax(particles[0].center[2]);
	tbb::parallel_for(tbb::blocked_range<size_t>(0, particles.size()), [&](const tbb::blocked_range<size_t>& r)
		{
			auto& Xmin = xmin.local();
			auto& Xmax = xmax.local();
			auto& Ymin = ymin.local();
			auto& Ymax = ymax.local();
			auto& Zmin = zmin.local();
			auto& Zmax = zmax.local();
			for (size_t i = r.begin(); i < r.end(); ++i)
			{
				double x = particles[i].center[0];
				double y = particles[i].center[1];
				double z = particles[i].center[2];
				Xmin = std::min(Xmin, x);
				Xmax = std::max(Xmax, x);
				Ymin = std::min(Ymin, y);
				Ymax = std::max(Ymax, y);
				Zmin = std::min(Zmin, z);
				Zmax = std::max(Zmax, z);
			}
		});
	double Xmin = *std::min_element(xmin.begin(), xmin.end());
	double Ymin = *std::min_element(ymin.begin(), ymin.end());
	double Zmin = *std::min_element(zmin.begin(), zmin.end());
	double Xmax = *std::max_element(xmax.begin(), xmax.end());
	double Ymax = *std::max_element(ymax.begin(), ymax.end());
	double Zmax = *std::max_element(zmax.begin(), zmax.end());

	double dw = std::max({Xmax - Xmin, Ymax - Ymin, Zmax - Zmin});
	return { {Xmin, Ymin, Zmin}, dw };
}

template <int dim, typename point, typename value>
MortonTree<dim, point, value>::MortonTree(std::vector<particle_t>&& particles_, size_t tree_depth_) : tree_depth(tree_depth_), particles(particles_)
{
	double T = omp_get_wtime();
	//if (IAmRoot()) std::cout << "\n************ Start Building Tree ************" << std::endl;
	const int dim2 = MyPow(2, dim);
	double t = omp_get_wtime();
	{
		auto be = box_extent(particles);
		shift = be.first; scale = be.second;

		while (scale  < FORCE_EPS * (1 << tree_depth))
		{
			--tree_depth;
		}

		if constexpr (dim == 2)
			shift -= std::complex<double> {1.e-7, 1.e-7}; // min particles shouldnt have coordinates equal zero
		else
			shift -= Vector3d {1.e-7, 1.e-7, 1.e-7}; // min particles shouldnt have coordinates equal zero
		scale += 1.e-6; // max particles shouldnt have coordinate equal 1
	}
	//if (IAmRoot()) std::cout << "compute extent time: " << omp_get_wtime() - t << std::endl;
	//if (IAmRoot()) std::cout << "scale factors: shift = " << shift << ", extent = " << scale << std::endl;

	size_t num_particles = particles.size();
	// if (IAmRoot()) std::cout << "particles num = " << num_particles << std::endl;

	levels.resize(tree_depth);
	tbb::parallel_for(size_t(0), tree_depth, [&](size_t i) {
		levels[i].resize(std::min(num_particles, (size_t)MyPow(size_t(dim2), (unsigned int)i)));
		});
	level_sizes.resize(tree_depth);
	positions_map.resize(num_particles);

	t = omp_get_wtime();
	tbb::parallel_for_each(particles.begin(), particles.end(), [&](particle_t& x) {
		x.morton_code = MortonCode<dim, point>(x.center, (int)(tree_depth - 1), shift, scale);
		assert(x.morton_code < MyPow(dim2, tree_depth - 1));
		});
	//if (IAmRoot()) std::cout << "compute morton codes time: " << omp_get_wtime() - t << std::endl;

	assert(MyPow(dim2, tree_depth) > std::max_element(particles.begin(), particles.end(), [&](const particle_t& x, const particle_t& y) {return x.morton_code < y.morton_code; })->morton_code);
	t = omp_get_wtime();
	{
#ifdef deterministic_behaviour
		tbb::parallel_sort(particles.begin(), particles.end(), [](const auto& x, const auto& y) {return x.morton_code < y.morton_code; });
		//boost::sort::block_indirect_sort(particles.begin(), particles.end(), [](const auto& x, const auto& y) {return x.morton_code < y.morton_code; }, 1);
#else
		if (IAmRoot())
		{
			// counting sort, using known number of unique numbers
			tbb::concurrent_vector<std::atomic_size_t> counts(MyPow(dim2, (unsigned int)(tree_depth - 1)) + 1);
			tbb::parallel_for_each(particles.begin(), particles.end(), [&](const auto& x) {counts[1 + x.morton_code]++; });
			for (int i = 1; i < counts.size(); ++i)
				counts[i] += counts[i - 1];
			assert(counts.back() == num_particles);

			tbb::scalable_allocator<small_particle<dim, value>> temp;
			constexpr size_t memsize = sizeof(small_particle<dim, value>);
			auto ptr = temp.allocate(num_particles);

			//t = omp_get_wtime();
			tbb::parallel_for(size_t(0), num_particles, [&](size_t i) {
				const auto& x = particles[i];
				size_t pos = counts[x.morton_code]++;
				assert(pos < num_particles);
				positions_map[i] = pos;
				memcpy(&ptr[pos], &x, memsize);
				});

			tbb::parallel_for(size_t(0), num_particles, [&](size_t i) {
				memcpy(&particles[i], &ptr[i], memsize);
				});

			temp.deallocate(ptr, num_particles);

			assert(std::is_sorted(particles.begin(), particles.end(), [](const auto& x, const auto& y) {return x.morton_code < y.morton_code; }));
		}
#ifdef FMM_MPI
		MPI_Bcast(particles.data(), num_particles * sizeof(particle_t), MPI_BYTE, RootID(), MPI_COMM_WORLD);
		MPI_Bcast(positions_map.data(), num_particles * sizeof(size_t), MPI_BYTE, RootID(), MPI_COMM_WORLD);
#endif

#endif
	}
	//if (IAmRoot()) std::cout << "sort time: " << omp_get_wtime() - t << std::endl;
	t = omp_get_wtime();
#ifdef deterministic_behaviour
	tbb::parallel_for_each(particles.begin(), particles.end(), [&](auto& x) {x.morton_code >>= dim * (30 - tree_depth); });
#endif
	{
		auto& top_level = levels.back();
		auto ptr = std::unique_copy(exec_policy, particles.begin(), particles.end(), top_level.begin(), [&](const auto& x, const auto& y) {return (x.morton_code) == (y.morton_code); });
		level_sizes.back() = ptr - top_level.begin();

		int idx_bottom = 0, idx_top = 0;
		for (int j = 0; j < num_particles; ++j)
		{
			if ((particles[j].morton_code) != (top_level[idx_top].morton_code))
			{
				auto& cell = top_level[idx_top];
				cell.source_range = { idx_bottom, j };
				cell.center = PointFromMortonCode(cell.morton_code, (int)(tree_depth - 1), shift, scale);
				++idx_top;
				idx_bottom = j;
			}
		}
		auto& cell = top_level[idx_top];
		cell.source_range = { idx_bottom, num_particles };
		cell.center = PointFromMortonCode(cell.morton_code, (int)(tree_depth - 1), shift, scale);

		//if (IAmRoot()) std::cout << "level " << tree_depth - 1 << " size: " <<  level_sizes.back() << std::endl;
	}

	for (int i = (int)(tree_depth - 2); i >= 0; --i)
	{
		auto& top_level = levels[i];
		auto& bottom_level = levels[i + 1];

		auto ptr = std::unique_copy(exec_policy, bottom_level.begin(), bottom_level.begin() + level_sizes[i + 1], top_level.begin(), [&](const auto& x, const auto& y) {return (x.morton_code >> dim) == (y.morton_code >> dim); });
		level_sizes[i] = ptr - top_level.begin();

		//if (IAmRoot()) std::cout << "level " << i << " size: " << level_sizes[i] << std::endl;

		int idx_bottom = 0, idx_top = 0;
		for (int j = 0; j < level_sizes[i + 1]; ++j)
		{
			if ((bottom_level[j].morton_code >> dim) != (top_level[idx_top].morton_code >> dim))
			{
				auto& cell = top_level[idx_top];
				cell.source_range = { idx_bottom, j };
				cell.morton_code >>= dim;
				cell.center = PointFromMortonCode(cell.morton_code, i, shift, scale);
				++idx_top;
				idx_bottom = j;
			}
			bottom_level[j].parent = idx_top;
		}
		auto& cell = top_level[idx_top];
		cell.source_range = { idx_bottom, level_sizes[i + 1] };
		cell.morton_code >>= dim;
		cell.center = PointFromMortonCode(cell.morton_code, i, shift, scale);
	}
	//if (IAmRoot()) std::cout << "create levels time: " << omp_get_wtime() - t << std::endl;

	t = omp_get_wtime();
	{
		std::vector<size_t> inverse_positions_map(num_particles);
		tbb::parallel_for(size_t(0), num_particles, [&](size_t i){
			inverse_positions_map[positions_map[i]] = i;
		});

		value Zero;
		if constexpr (std::is_same_v<value, double>)
			Zero = 0;
		if constexpr (std::is_same_v<value, Vector3d>)
			Zero = { 0,0,0 };
		auto& leaves = levels.back();
		tbb::enumerable_thread_specific<size_t> local_num(0);
		tbb::parallel_for(size_t(0), level_sizes.back(), [&](size_t k){
			auto& leaf = leaves[k];
			auto& [s_begin, s_end] = leaf.source_range;
			auto& [t_begin, t_end] = leaf.target_range;
			t_end = s_end;
			t_begin = s_begin;
			for (size_t i = s_begin; i < s_end; ++i)
			{
				if (particles[i].q != Zero)
				{
					std::swap(particles[t_begin], particles[i]);
					std::swap(inverse_positions_map[t_begin++], inverse_positions_map[i]);
				}
			}
			s_end = t_begin;
			local_num.local() += t_end - t_begin;

		});

		tbb::parallel_for(size_t(0), num_particles, [&](size_t i){
				positions_map[inverse_positions_map[i]] = i;
			});

		for (const auto& x : local_num)
			targets_num += x;
		//if (IAmRoot()) std::cout << "targets num: " << targets_num << std::endl;

	}
	//if (IAmRoot()) std::cout << "source/target partition time: " << omp_get_wtime() - t << std::endl;

	t = omp_get_wtime();
	FillNeighbours();
	//if (IAmRoot()) std::cout << "fill neighbours time: " << omp_get_wtime() - t << std::endl;

	//������ ������ ������� ������������
	/*
	t = omp_get_wtime();
	if (targets_num == 0)
	{
		auto& leaves = levels.back();
		for (int i = 0; i < level_sizes.back(); ++i)
		{
			const auto& cell = leaves[i];
			for (const auto& idx : cell.closeneighbours)
			{
				auto& neighbour_cell = leaves[idx];
				auto ptr = std::find(neighbour_cell.closeneighbours.begin(), neighbour_cell.closeneighbours.end(), i);
				neighbour_cell.closeneighbours.erase(ptr);
			}
		}
	}
	if (IAmRoot()) std::cout << "remove duplicate nbr indices time:" << omp_get_wtime() - t << std::endl;
	*/

	//if (IAmRoot()) std::cout << "total build time: " << omp_get_wtime() - T << std::endl;
	//if (IAmRoot()) std::cout << "************ End Building Tree ************\n\n";
}

template class MortonTree<2, point2d, double>;
template class MortonTree<3, point3d, double>;
template class MortonTree<3, point3d, point3d>;

} // fmm
