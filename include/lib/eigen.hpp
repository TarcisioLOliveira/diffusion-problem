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

#ifndef DPLIB_EIGEN_HPP
#define DPLIB_EIGEN_HPP

#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/src/OrderingMethods/Ordering.h>
#include <Eigen/src/SparseCholesky/SimplicialCholesky.h>
#include <cstddef>
#include "lib/sparse_matrix.hpp"

namespace dplib{

class EigenPCG{
    public:
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, std::ptrdiff_t> Mat;

    void set_K(SparseMatrix& M, size_t L);
    void compute();
    void solve(std::vector<double>& x, std::vector<double>& b);

    inline void reset(){
        this->first_time = true;
    }

    private:
    bool first_time = true;
    Mat K;
    Eigen::ConjugateGradient<Mat, Eigen::Lower|Eigen::Upper> cg;
};

class EigenCholesky{
    public:
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor, std::ptrdiff_t> Mat;

    void set_K(SparseMatrix& M, size_t L);
    void compute();
    void solve(std::vector<double>& x, std::vector<double>& b);

    inline void reset(){
        this->first_time = true;
    }

    private:
    bool first_time = true;
    Mat K;
    Eigen::SimplicialCholesky<Mat, Eigen::Lower, Eigen::AMDOrdering<std::ptrdiff_t>> solver;
};


}

#endif
