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


from libcpp cimport bool

cdef extern from "cellcutoff/cell.h" namespace "cellcutoff":
    cdef cppclass Cell:
        Cell(double* vecs, int nvec) except +
        Cell()

        Cell* create_subcell(const double threshold, int* shape) except +
        Cell* create_reciprocal()

        int nvec()
        const double* vecs()
        const double* gvecs()
        double volume()
        double gvolume()
        const double* lengths()
        const double* glengths()
        const double* spacings()
        const double* gspacings()
        const bool cubic()
        const bool cuboid()

        void to_frac(double* cart, double* frac)
        void to_cart(double* frac, double* cart)

        void iwrap_mic(double* delta);
        void iwrap_box(double* delta);

        size_t ranges_cutoff(const double* center, double cutoff, int* ranges_begin,
            int* ranges_end) const;

    Cell* create_random_cell(const unsigned int seed, const int nvec,
        const double scale, const double ratio, const bool cuboid)
