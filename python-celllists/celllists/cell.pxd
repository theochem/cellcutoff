#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CellList is a 3D domain decomposition library.
# Copyright (C) 2011-2015 The CellList Development Team
#
# This file is part of CellList.
#
# CellList is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# CellList is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--


from libcpp cimport bool

cdef extern from "celllists/cell.h" namespace "celllists":
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
