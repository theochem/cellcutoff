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


#ifndef __CELLLISTS_SPHERE_SLICE_H__
#define __CELLLISTS_SPHERE_SLICE_H__


namespace celllists {


/** @brief.
        An exception for when solve_range(_*) can not find a solution.
 */
class no_solution_found : public std::domain_error {
    public:
        explicit no_solution_found (const std::string& what_arg) : std::domain_error(what_arg) {};
};

class SphereSlice {
    private:
        // Constant independent data members
        const double* center;
        const double* normals;
        const double radius;

        // Configurable data members
        double cut_begin[2];
        double cut_end[2];

        // Derived from constant data members upon construction
        double radius_sq;
        double norms_sq[3];
        double norms[3];
        double frac_radii[3];
        double frac_center[3];
        double sphere_frac_begin[3];
        double sphere_frac_end[3];
        double radius_normals[9];
        double sphere_point_begin[9];
        double sphere_point_end[9];
        double dots[9];
        double denoms[9];
        double cut_ortho[27];
    public:
        SphereSlice(const double* center, const double* normals, double radius);

        // Main API
        void solve_range(int ncut, double &begin, double &end) const;
        void set_cut_begin_end(int icut, double new_begin, double new_end);

        // Auxiliary API, could also be useful and there is no need to really
        // make this private. Having it public also facilitates testing.
        bool inside_cuts(int id_cut, double* point) const;

        void solve_full_low(int id_axis, double &begin, double &end,
            double* point_begin, double* point_end) const;
        void solve_full(int id_axis, double &begin, double &end,
            int id_cut0, int id_cut1) const;

        void solve_plane_low(int id_axis, int id_cut, double frac_cut,
            double &begin, double &end, double* point_begin, double* point_end)
            const;
        void solve_plane(int id_axis, int id_cut0, double frac_cut,
            double &begin, double &end, int id_cut1) const;

        double compute_plane_intersection(int id_cut0, int id_cut1,
            double cut0, double cut1, double* other_center) const;
        void solve_line_low(int id_axis, int id_cut0, int id_cut1,
            double frac_cut0, double frac_cut1, double &begin, double &end,
            double* point_begin, double* point_end) const;
        void solve_line(int id_axis, int id_cut0, int id_cut1,
            double frac_cut0, double frac_cut1, double &begin, double &end) const;

        void solve_range_0(double &begin, double &end) const;
        void solve_range_1(double &begin, double &end) const;
        void solve_range_2(double &begin, double &end) const;
};

void compute_begin_end(const double* other_center, const double* ortho,
    const double* axis, double &begin, double &end,
    double* point_begin, double* point_end);

void update_begin_end(double work_begin, double work_end, double &begin, double &end);


}

#endif
