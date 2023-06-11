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

#include "lib/eigen.hpp"

namespace dplib{

void EigenPCG::set_K(SparseMatrix& M, size_t L){
    std::fill(this->K.valuePtr(), this->K.valuePtr() + this->K.nonZeros(), 0);
    if(this->first_time){
        this->K = Mat(L, L);
    }
    M.to_eigen_sparse(this->K);
}

void EigenPCG::compute(){
    this->cg.compute(this->K);
}

void EigenPCG::solve(std::vector<double>& x, std::vector<double>& b){
    Eigen::VectorXd f = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), b.size());
    Eigen::VectorXd u = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(x.data(), x.size());

    u = this->cg.solveWithGuess(f, u);

    std::copy(u.cbegin(), u.cend(), x.begin());
}

// CHOLESKY

void EigenCholesky::set_K(SparseMatrix& M, size_t L){
    std::fill(this->K.valuePtr(), this->K.valuePtr() + this->K.nonZeros(), 0);
    if(this->first_time){
        this->K = Mat(L, L);
    }
    M.to_eigen_sparse(this->K);
}

void EigenCholesky::compute(){
    this->solver.compute(this->K);
}

void EigenCholesky::solve(std::vector<double>& x, std::vector<double>& b){
    Eigen::VectorXd f = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), b.size());

    Eigen::VectorXd u = this->solver.solve(f);

    std::copy(u.cbegin(), u.cend(), x.begin());
}

}
