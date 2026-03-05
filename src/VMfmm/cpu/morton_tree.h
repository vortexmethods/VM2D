#pragma once
#include <utility>
#include <vector>
#include <complex>
#include "../common/utils.h"

namespace fmm {

template <typename point, typename value>
struct TreeCell
{
	using particle_t = particle<point, value>;

	TreeCell& operator=(const particle_t& particle) { morton_code = particle.morton_code; return *this; }
	point center;

	std::pair<size_t, size_t> source_range;
	std::pair<size_t, size_t> target_range;
	size_t parent;
	size_t morton_code;
	std::vector<size_t> closeneighbours, farneighbours;
};

template <int dim, typename point, typename value>
class MortonTree
{
public:
	using particle_t = particle<point, value>;
	using TreeCell_t = TreeCell<point, value>;

	MortonTree(std::vector<particle_t>&& particles, size_t tree_depth);

	std::vector<std::vector<TreeCell_t>> levels;
	std::vector<size_t> level_sizes;
	size_t tree_depth;
	size_t targets_num = 0; // number of target points in tree
	std::vector<particle_t> particles;
	std::vector<size_t> positions_map; // maps sorted to initial

private:
	void FillNeighbours();
	point shift;
	double scale;

};

using TreeCell2d = TreeCell<point2d, double>;
using TreeCell3d = TreeCell<point3d, double>;
using TreeCell3d3 = TreeCell<point3d, point3d>;
using MortonTree2d = MortonTree<2, point2d, double>;
using MortonTree3d = MortonTree<3, point3d, double>;
using MortonTree3d3 = MortonTree<3, point3d, point3d>;

} // fmm