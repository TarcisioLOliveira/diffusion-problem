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

#ifndef DPLIB_PRINT_HPP
#define DPLIB_PRINT_HPP

#include <iostream>

namespace dplib{

inline void print_line(std::string s){
    std::cout << s << std::endl;
}

}

#endif
