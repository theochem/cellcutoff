// CellList is a 3D domain decomposition library.
// Copyright (C) 2011-2015 The CellList Development Team
//
// This file is part of CellList.
//
// CellList is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// CellList is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// aint with this program; if not, see <http://www.gnu.org/licenses/>
//
//--


#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include "common.h"
#include "celllists/vec3.h"


//! Fills an array of doubles with random numbers in range ]-0.5*scale, 0.5*scale]
int fill_random_double(unsigned int seed, double* array, size_t size,
    double low, double high) {

    if (size <= 0)
        throw std::domain_error("Array size must be strictly positive."); // IGNORE_COVERAGE

    srand(seed);
    for (size_t i=0; i < size; i++)
        array[i] = (rand() + 1.0)/(RAND_MAX + 1.0)*(high - low) + low;
    return rand();
}

//! Fills an array of int with random numbers in range [-range, range]
int fill_random_int(unsigned int seed, int* array, size_t size,
    int begin, int end) {

    if (size <= 0)
        throw std::domain_error("Array size must be strictly positive."); // IGNORE_COVERAGE
    if (begin > end)
        throw std::domain_error("Begin cannot be larger than end."); // IGNORE_COVERAGE
    srand(seed);
    for (size_t i=0; i < size; i++)
        array[i] = (rand() % (end - begin)) + begin;
    return rand();
}

int myrandom(int i) { return rand()%i; }

//! Fills and array of int with a random permutation
int fill_random_permutation(unsigned int seed, int* array, size_t size) {

    if (size <= 0)
        throw std::domain_error("Array size must be strictly positive."); // IGNORE_COVERAGE
    for (size_t i=0; i < size; i++)
        array[i] = i;
    srand(seed);
    std::random_shuffle(array, array+size, myrandom);
    return rand();
}

//! Random cell with a volume larger than 0.01
Cell* create_random_cell_nvec(unsigned int seed, int nvec, double scale, bool cuboid) {
    if ((nvec <= 0) || (nvec > 3)) {
        throw std::domain_error("A random cell must be 1D, 2D or 2D periodic."); // IGNORE_COVERAGE
    }
    double rvecs[nvec*3];
    while (true) {
        seed = fill_random_double(seed, rvecs, nvec*3, -scale, +scale);
        if (cuboid) {
            rvecs[1] = 0.0;
            rvecs[2] = 0.0;
            if (nvec > 1) {
                rvecs[3] = 0.0;
                rvecs[5] = 0.0;
            }
            if (nvec > 2) {
                rvecs[6] = 0.0;
                rvecs[7] = 0.0;
            }
        }
        try {
            Cell* cell = new Cell(rvecs, nvec);
            if (cell->get_volume() > pow(0.1*scale, nvec))
                return cell;
            delete cell;
        } catch (singular_cell_vectors) { }
    }
}

//! Compute a random point in a cubic box centered around center. Also computes distance.
void random_point(unsigned int seed, double* point, double rcut, const double* center,
    double &norm) {
    fill_random_double(seed, point, 3, -rcut, rcut);
    norm = vec3::norm(point);
    vec3::iadd(point, center);
}
