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

#include <algorithm>
#include <numeric>
#include <iomanip>
#include <cblas.h>
#include "lib/print.hpp"
#include "lib/Q4.hpp"
#include "lib/window.hpp"
#include "lib/mesh.hpp"
#include "lib/sparse_matrix.hpp"
#include "lib/eigen.hpp"

int main(){
    Eigen::initParallel();

    dplib::print_line("Launching window...");
    const size_t window_width = 600;
    const size_t window_height = 500;

    const size_t W = 400;
    const size_t H = 400;

    const double E_SIZE = 1;

    const double K_MIN = 0;
    const double I_MIN = 1e-9;

    dplib::Window window(window_width, window_height, W, H, "test3 - psi");

    dplib::print_line("Creating mesh...");
    dplib::RectangularMesh mesh(W, H, 1.0, E_SIZE);

    mesh.apply_Dirichlet(1, {0,0,0}, {W+1,0,0});
    mesh.apply_Dirichlet(1, {W+1,0,0}, {W+1,H+1,0});
    mesh.apply_Dirichlet(1, {0,0,0}, {0,H+1,0});
    mesh.apply_Dirichlet(1, {0,H+1,0}, {W+1,H+1,0});

    dplib::print_line("Generating global matrix...");
    mesh.generate_K(K_MIN);

    for(size_t i = 0; i < mesh.matrix_size(); ++i){
        mesh.K.add(i, i, I_MIN);
    }

    dplib::print_line("Solving linear equation...");
    mesh.solve();

    dplib::print_line("Displaying results...");
    auto result = mesh.get_result();
    double maxx = *std::max_element(result.begin(), result.end());
    double minx = *std::min_element(result.begin(), result.end());

    std::cout << minx << " " << maxx << std::endl;
    window.update(result, 0, 1);
    do{
        window.update();
    } while(window.is_open());

    return 0;
}
