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
#include <cstdlib>
#include "common.h"


//! Fills an array of doubles with random numbers in range ]-0.5*scale, 0.5*scale]
int fill_random_double(unsigned int seed, double* array, size_t size, double scale) {
    srand(seed);
    for (int i=0; i<size; i++)
        array[i] = (rand() + 1.0)/(RAND_MAX + 1.0) - 0.5;
    return rand();
}

//! Fills an array of int with random numbers in range [-range, range]
int fill_random_int(unsigned int seed, int* array, size_t size, int range) {
    srand(seed);
    for (int i=0; i<size; i++)
        array[i] = (rand() % (2*range+1))-range;
    return rand();
}

//! Random cell with a volume larger than 0.01
Cell* create_random_cell_nvec(unsigned int seed, int nvec, double scale, bool cuboid) {
    if (nvec == 0) {
        throw std::domain_error("A random cell must be at least 1D periodic.");
    }
    double rvecs[nvec*3];
    while (true) {
        seed = fill_random_double(seed, rvecs, nvec*3, scale);
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
            if (cell->get_volume() > 0.01)
                return cell;
            delete cell;
        } catch (singular_cell_vectors) { }
    }
}
