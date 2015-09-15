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


#ifndef CELLLISTS_TESTS_COMMON_H_
#define CELLLISTS_TESTS_COMMON_H_

#include <memory>

#include "celllists/cell.h"

namespace cl = celllists;

/* Some notes:

    - Usually, the output argument is last. The exception in this module is the size
      of the output argument, which comes after the actual output argument.

*/

#define NREP 100
#define NPOINT 1000

unsigned int fill_random_double(unsigned int seed, double* array, int size,
    double low = -0.5, double high = 0.5);
unsigned int fill_random_int(unsigned int seed, int* array, int size,
    int begin, int end);
unsigned int fill_random_permutation(unsigned int seed, int* array, int size);
std::unique_ptr<cl::Cell> create_random_cell_nvec(unsigned int seed, int nvec,
    double scale = 1.0, bool cuboid = false);
unsigned int random_point(unsigned int seed,  const double* center, double rcut,
    double* point, double* norm);

#endif  // CELLLISTS_TESTS_COMMON_H_
