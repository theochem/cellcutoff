// CellCutoff is a library for periodic boundary conditions and real-space cutoff calculations.
// Copyright (C) 2017 The CellCutoff Development Team
//
// This file is part of CellCutoff.
//
// CellCutoff is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// CellCutoff is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
// --


/** @file

    Low-level helper functions to facilitate python wrapping of cellcutoff.

  */

#ifndef CELLCUTOFF_WRAPPING_H_
#define CELLCUTOFF_WRAPPING_H_


#include <vector>

#include "cellcutoff/iterators.h"


namespace cl = cellcutoff;

/** @brief
        Put all results from a cutoff calculation into C++ vectors.

    @param bsp
        A cl::BoxSortedPoints instance.

    @param center
        The center of the cutoff sphere.

    @param radius
        The radius of the cutoff sphere.

    @param dds
        The vector to store relative vectors and distances.

    @param ipoints
        The vector to store the piont indexes in the original points array.

    @return npoint
        The number of points within the cutoff sphere.
 */
int box_cutoff_points(const cl::BoxSortedPoints* bsp, const double* center,
                      const double radius, std::vector<double> *dds,
                      std::vector<size_t> *ipoints);


#endif  // CELLCUTOFF_WRAPPING_H_

// vim: textwidth=90 et ts=2 sw=2
