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

#include <cblas.h>
#include <cmath>
#include <set>
#include <queue>
#include <algorithm>
#include "lib/mesh.hpp"
#include "lib/print.hpp"
#include "lib/Q4.hpp"

namespace dplib{

RectangularMesh::RectangularMesh(size_t W, size_t H, double t, double elem_size):
    W(W), H(H), element_size(elem_size), t(t), element_nodes(W*H*nodes_per_element, 0),
    node_vector_mapping((W+1)*(H+1)*dof_per_node, 0){
    dplib::print_line("Mesh: generating mesh...");
    for(size_t y = 0; y < H; ++y){
        for(size_t x = 0; x < W; ++x){
            const size_t e = (y*W + x)*this->nodes_per_element;
            const size_t n = x + y*(W+1);
            // Ordering matters! (because of Q4 assumptions when generating matrices)
            element_nodes[e+0] = n+W+1;
            element_nodes[e+1] = n+W+2;
            element_nodes[e+2] = n+1;
            element_nodes[e+3] = n;
        }
    }

    dplib::print_line("Mesh: running RCM...");
    this->old_position_mapping.resize((W+1)*(H+1), 0);
    reverse_cuthill_mckee(element_nodes, old_position_mapping, nodes_per_element, (W+1)*(H+1));
}

void RectangularMesh::apply_Neumann(double d, Point begin, Point end){
    if(begin.x == end.x && begin.x == W+1){
        begin.x -= 1;
        end.x -= 1;
    }
    if(begin.y == end.y && begin.y == H+1){
        begin.y -= 1;
        end.y -= 1;
    }
    this->neumann.push_back({d, begin, end});
}

void RectangularMesh::apply_Dirichlet(double d, Point begin, Point end){
    dplib::print_line("Mesh: generating mesh...");
    long id = this->dirichlet.size() + 1;
    size_t extension = std::max(std::abs(end.x - begin.x), std::abs(end.y - begin.y));
    this->dirichlet.reserve(this->dirichlet.size() + extension);
    if(begin.x == end.x && begin.x == W+1){
        begin.x -= 1;
        end.x -= 1;
    }
    if(begin.y == end.y && begin.y == H+1){
        begin.y -= 1;
        end.y -= 1;
    }
    for(size_t x = begin.x; x < end.x || (x == begin.x && x == end.x); ++x){
        for(size_t y = begin.y; y < end.y || (y == begin.y && y == end.y); ++y){
            const size_t n = this->old_position_mapping[y*(W+1) + x];
            for(size_t i = 0; i < this->dof_per_node; ++i){
                node_vector_mapping[n*dof_per_node+i] = -id;
            }
            this->dirichlet.push_back(d);
            ++id;
        }
    }
}

void RectangularMesh::generate_K(const double K_MIN){
    long id = 0;
    for(auto& n:node_vector_mapping){
        if(n > -1){
            n = id;
            ++id;
        }
    }

    dplib::print_line("Mesh: generating Neumann vector...");
    this->load.resize(id, 0);
    this->psi.resize(id, 0);
    for(const auto& n:this->neumann){
        const Point& begin = n.begin;
        const Point& end = n.end;
        const double de = n.d/(begin.distance(end)*this->element_size);
        bool first = true;
        bool last = false;
        for(size_t x = begin.x; x < end.x || (x == begin.x && x == end.x); ++x){
            if(begin.x != end.x && x + 1 == end.x){
                last = true;
            }
            for(size_t y = begin.y; y < end.y || (y == begin.y && y == end.y); ++y){
                if(begin.y != end.y && y + 1 == end.y){
                    last = true;
                }
                const size_t n = this->old_position_mapping[y*(W+1) + x];
                for(size_t i = 0; i < this->dof_per_node; ++i){
                    const long u1_id = node_vector_mapping[n*dof_per_node+i];
                    if(u1_id > -1){
                        if(first || last){
                            this->load[u1_id] += de*this->element_size/2;
                        } else {
                            this->load[u1_id] += de*this->element_size;
                        }
                    }
                }
            }
        }
    }

    dplib::print_line("Mesh: generating global matrix and Dirichlet vector...");
    std::vector<double> A{1.0, 0.0,
                          0.0, 1.0};
    const auto k = dplib::Q4::get_diffusion_2D(this->t, this->element_size/2, this->element_size/2, A);
    std::vector<double> rho_k(k);
    std::vector<long> u_pos(this->nodes_per_element*this->dof_per_node, 0);
    for(size_t y = 0; y < H; ++y){
        for(size_t x = 0; x < W; ++x){
            const size_t e = (y*W + x);
            const Point p{static_cast<double>(x), static_cast<double>(y), 0.0};
            for(size_t n = 0; n < this->nodes_per_element; ++n){
                const size_t node_id = this->element_nodes[e*this->nodes_per_element + n];
                for(size_t i = 0; i < this->dof_per_node; ++i){
                    const size_t dof_id = node_id*this->dof_per_node + i;
                    u_pos[n*this->dof_per_node + i] = this->node_vector_mapping[dof_id];
                }
            }
            std::copy(k.begin(), k.end(), rho_k.begin());
            const double rho = this->ring(p, K_MIN);
            cblas_dscal(rho_k.size(), rho, rho_k.data(), 1);
            this->K.insert_matrix_symmetric_mumps(rho_k, u_pos);
            // Add Dirichlet boundary conditions
            if(x == 0 || x == W-1 || y == 0 || y == H-1){
                for(size_t i = 0; i < u_pos.size(); ++i){
                    if(u_pos[i] < 0){
                        continue;
                    }
                    for(size_t j = 0; j < u_pos.size(); ++j){
                        if(u_pos[j] < 0){
                            long dirich_id = -(u_pos[j]+1);
                            this->load[u_pos[i]] -= this->dirichlet[dirich_id]*rho_k[i*u_pos.size() + j];
                        }
                    }
                }
            }
        }
    }
}

void RectangularMesh::solve(){
    solver.set_K(this->K, this->load.size());

    solver.compute();
    solver.solve(this->psi, this->load);
}
    
std::vector<double> RectangularMesh::get_result(){
    std::vector<double> result(W*H, 0);

    for(size_t y = 0; y < H; ++y){
        for(size_t x = 0; x < W; ++x){
            const size_t e = (y*W + x);
            for(size_t n = 0; n < this->nodes_per_element; ++n){
                const size_t node_id = this->element_nodes[e*this->nodes_per_element + n];
                for(size_t i = 0; i < this->dof_per_node; ++i){
                    const size_t dof_id = node_id*this->dof_per_node + i;
                    const long pos = this->node_vector_mapping[dof_id];
                    if(pos > -1){
                        result[e] += this->psi[pos]/4;
                    } else {
                        const long d_pos = -(pos+1);
                        result[e] += this->dirichlet[d_pos]/4;
                    }
                }
            }
        }
    }

    return result;
}


double RectangularMesh::ring(const Point& p, double min){
    const Point center{W/2.0, H/2.0, 0};

    const double dist = center.distance(p);
    const double ri = std::min(W, H)/6.0;
    const double ro = 2*ri;
    if(dist >= ri && dist <= ro){
        return min;
    } else {
        return 1;
    }
}

void reverse_cuthill_mckee(std::vector<size_t>& element_nodes, std::vector<size_t>& old_position_mapping, const size_t nodes_per_element, const size_t number_of_nodes){
    dplib::print_line("Mesh: RCM: creating adjacency matrix...");
    // Create adjacency "matrix"
    std::vector<std::set<size_t>> adjacents(number_of_nodes);
    const size_t number_of_elements = element_nodes.size()/nodes_per_element;
    for(size_t e = 0; e < number_of_elements; ++e){
        for(size_t i = 0; i < nodes_per_element; ++i){
            for(size_t j = i+1; j < nodes_per_element; ++j){
                const size_t ni = element_nodes[e*nodes_per_element+i];
                const size_t nj = element_nodes[e*nodes_per_element+j];
                adjacents[ni].insert(nj);
                adjacents[nj].insert(ni);
            }
        }
    }
    std::vector<bool> added(number_of_nodes, false);

    // Get node with least degree
    dplib::print_line("Mesh: RCM: getting node with least degree...");
    size_t min_degree = adjacents[0].size();
    size_t min_node = 0;
    for(size_t i = 1; i < adjacents.size(); ++i){
        if(adjacents[i].size() < min_degree){
            min_node = i;
            min_degree = adjacents[i].size();
        }
    }

    // Generate Cuthill-McKee
    dplib::print_line("Mesh: RCM: reorganizing nodes...");
    auto comp = [&](size_t n1, size_t n2){
        return adjacents[n1].size() < adjacents[n2].size();
    };
    std::queue<size_t> queue;
    queue.push(min_node);
    added[min_node] = true;
    std::vector<size_t> result;
    result.reserve(number_of_nodes);
    while(!queue.empty()){
        size_t node = queue.front();
        auto& adj = adjacents[node];

        std::vector<size_t> new_nodes;
        new_nodes.reserve(adj.size());
        for(auto& n:adj){
            if(!added[n]){
                added[n] = true;
                auto upper = std::upper_bound(new_nodes.begin(), new_nodes.end(), n, comp);
                new_nodes.insert(upper, n);
            }
        }
        for(auto& n:new_nodes){
            queue.push(n);
        }

        result.push_back(node);
        queue.pop();
    }

    // Reorder node list
    dplib::print_line("Mesh: RCM: finishing up...");
    std::vector<size_t> new_node_mapping(number_of_nodes, 0);
    size_t pos = 0;
    for(auto i = result.rbegin(); i < result.rend(); ++i){
        new_node_mapping[*i] = pos;
        ++pos;
    }
    std::copy(new_node_mapping.begin(), new_node_mapping.end(), old_position_mapping.begin());
    for(auto& n:element_nodes){
        n = new_node_mapping[n];
    }
}

}
