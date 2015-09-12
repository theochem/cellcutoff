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
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--

/** @file */

#ifndef CELLLIST_SPHERE_SLICE_H_
#define CELLLIST_SPHERE_SLICE_H_

class SphereSlice {
    private:
        const double* center;
        const double* normals;
        double radius;
        //double norms[2];
        double cut_begin[2];
        double cut_end[2];

    public:
        SphereSlice(const double* center, const double* normals, double radius);
        void solve_range(int ncut, double &begin, double &end) const;
        void set_cut_begin_end(int icut, double new_begin, double new_end);

        void solve_sphere(const double* axis, double &begin,
            double &end, double* point_begin, double* point_end) const;
        bool solve_circle(const double* axis, const double* cut_normal,
            double cut, double &begin, double &end, double* point_begin,
            double* point_end) const;
        void solve_range_0(double &begin, double &end) const;
        void solve_range_1(double &begin, double &end) const;
        void solve_range_2(double &begin, double &end) const;
};

#endif
