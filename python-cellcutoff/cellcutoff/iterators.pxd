# CellCutoff is a library for periodic boundary conditions and real-space cutoff calculations.
# Copyright (C) 2017 The CellCutoff Development Team
#
# This file is part of CellCutoff.
#
# CellCutoff is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# CellCutoff is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
# --
# cython: linetrace=True, embedsignature=True, language_level=3


cimport cellcutoff.cell as cell

cdef extern from "cellcutoff/iterators.h" namespace "cellcutoff":
    size_t cutoff_ranges(const cell.Cell* cell, const double* center,
        double cutoff, int* ranges_begin, int* ranges_end);

    cdef cppclass BoxSortedPoints:
        BoxSortedPoints(double* points, int npoint, cell.Cell* cell, double threshold) except +

        const double* points()
        size_t npoint()
        const cell.Cell* subcell()
        const int* shape()
        const size_t* ipoints()
