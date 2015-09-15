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


#ifndef CELLLISTTEST_COMMON_H_
#define CELLLISTTEST_COMMON_H_

#include "celllists/cell.h"

#define NREP 100
#define NPOINT 1000

unsigned int fill_random_double(unsigned int seed, double* array, int size,
    double low=-0.5, double high=0.5);
unsigned int fill_random_int(unsigned int seed, int* array, int size,
    int begin, int end);
unsigned int fill_random_permutation(unsigned int seed, int* array, int size);
Cell* create_random_cell_nvec(unsigned int seed, int nvec, double scale=1, bool cuboid=false);
unsigned int  random_point(unsigned int seed, double* point, double rcut, const double* center,
    double &norm);

#endif
