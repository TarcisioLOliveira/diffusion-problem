/*
 *   Copyright (C) 2023 Tarcísio Ladeia de Oliveira.
 *
 *   This file is part of diffusion-problem
 *
 *   diffusion-problem is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   diffusion-problem is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with diffusion-problem.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifndef DPLIB_MESH_HPP
#define DPLIB_MESH_HPP

#include <cstddef>
#include <functional>
#include <vector>
#include "lib/eigen.hpp"
#include "lib/sparse_matrix.hpp"

namespace dplib{

struct Point{
    double x = 0, y = 0, z = 0;

    inline double distance(const Point& p) const{
        const double dx = x - p.x;
        const double dy = y - p.y;
        const double dz = z - p.z;
        return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
};

void reverse_cuthill_mckee(std::vector<size_t>& element_nodes, std::vector<size_t>& old_position_mapping, const size_t nodes_per_element, const size_t number_of_nodes);

class RectangularMesh{
    public:
    RectangularMesh(size_t W, size_t H, double t, double elem_size);

    // Node range
    void apply_Dirichlet(double d, Point begin, Point end);
    // Node range
    void apply_Neumann(double d, Point begin, Point end);
    void generate_K(const double K_MIN);
    void solve();

    std::vector<double> get_result();

    dplib::SparseMatrix K;
    inline size_t matrix_size(){
        return this->load.size();
    }
    private:
    struct NeumannBoundary{
        double d;
        Point begin, end;
    };
    const size_t W, H;
    const size_t dof_per_node = 1;
    const size_t nodes_per_element = 4;
    const double element_size;
    const double t;
    std::vector<size_t> element_nodes;
    std::vector<long> node_vector_mapping;
    std::vector<double> load;
    std::vector<double> dirichlet;
    std::vector<NeumannBoundary> neumann;
    std::vector<double> psi;
    std::vector<size_t> old_position_mapping;
    dplib::EigenCholesky solver;

    double ring(const Point& p, double min);
};

//Mesh cantilever(size_t W, size_t H, double element_size, double fx, double fy, double f_len);
//void cantilever_convolution_map(size_t W, size_t H, double r, std::vector<std::vector<size_t>>& neighbors, std::vector<std::vector<double>>& distances, std::vector<double>& dist_total);

}

#endif
