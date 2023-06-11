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

#include "lib/Q4.hpp"

namespace dplib::Q4{

std::vector<double> get_diffusion_2D(double t, double a, double b, const std::vector<double>& A){
    std::vector<double> k{
    A[1]*t/4 + A[2]*t/4 + (A[0]*b*b*t + A[3]*a*a*t)/(3*a*b)
    ,
    A[1]*t/4 - A[2]*t/4 + (-2*A[0]*b*b*t + A[3]*a*a*t)/(6*a*b)
    ,
    -A[1]*t/4 - A[2]*t/4 + (-A[0]*b*b*t - A[3]*a*a*t)/(6*a*b)
    ,
    -A[1]*t/4 + A[2]*t/4 + (A[0]*b*b*t - 2*A[3]*a*a*t)/(6*a*b)
    ,
    -A[1]*t/4 + A[2]*t/4 + (-2*A[0]*b*b*t + A[3]*a*a*t)/(6*a*b)
    ,
    -A[1]*t/4 - A[2]*t/4 + (A[0]*b*b*t + A[3]*a*a*t)/(3*a*b)
    ,
    A[1]*t/4 - A[2]*t/4 + (A[0]*b*b*t - 2*A[3]*a*a*t)/(6*a*b)
    ,
    A[1]*t/4 + A[2]*t/4 + (-A[0]*b*b*t - A[3]*a*a*t)/(6*a*b)
    ,
    -A[1]*t/4 - A[2]*t/4 + (-A[0]*b*b*t - A[3]*a*a*t)/(6*a*b)
    ,
    -A[1]*t/4 + A[2]*t/4 + (A[0]*b*b*t - 2*A[3]*a*a*t)/(6*a*b)
    ,
    A[1]*t/4 + A[2]*t/4 + (A[0]*b*b*t + A[3]*a*a*t)/(3*a*b)
    ,
    A[1]*t/4 - A[2]*t/4 + (-2*A[0]*b*b*t + A[3]*a*a*t)/(6*a*b)
    ,
    A[1]*t/4 - A[2]*t/4 + (A[0]*b*b*t - 2*A[3]*a*a*t)/(6*a*b)
    ,
    A[1]*t/4 + A[2]*t/4 + (-A[0]*b*b*t - A[3]*a*a*t)/(6*a*b)
    ,
    -A[1]*t/4 + A[2]*t/4 + (-2*A[0]*b*b*t + A[3]*a*a*t)/(6*a*b)
    ,
    -A[1]*t/4 - A[2]*t/4 + (A[0]*b*b*t + A[3]*a*a*t)/(3*a*b)
    };

    return k;
}

}
