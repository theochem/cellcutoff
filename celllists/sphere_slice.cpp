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


#include <cmath>
#include <stdexcept>
#include "celllists/sphere_slice.h"


SphereSlice::SphereSlice(const double* center, const double* gvecs, double rcut) :
    center(center), normals(normals), rcut(rcut) {

    //norms[0] = sqrt(normals[0]*normals[0] + normals[1]*normals[1] + normals[2]*normals[2]);
    //norms[1] = sqrt(normals[3]*normals[3] + normals[4]*normals[4] + normals[5]*normals[5]);
    begin[0] = 0.0;
    begin[1] = 0.0;
    end[0] = 0.0;
    end[1] = 0.0;

}


void SphereSlice::solve_range_0(const double* normal, double &begin, double &end) const {
    double norm = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
    // Convert cutoff to fractional coordinates
    double frac_rcut = rcut*norm;
    // Convert the center to fractional coordinates
    double frac_center = center[0]*normal[0] + center[1]*normal[1] + center[2]*normal[2];
    // Find the ranges in fractional coordinates that encloses the cutoff
    begin = frac_center - frac_rcut;
    end = frac_center + frac_rcut;
}


void SphereSlice::solve_range_1(const double* normal, double &begin, double &end) const {
    throw std::logic_error("TODO");
}


void SphereSlice::solve_range_2(const double* normal, double &begin, double &end) const{
    throw std::logic_error("TODO");
}


void SphereSlice::solve_range(int ncut, const double* normal, double &begin, double &end) const {
    switch (ncut) {
        case 0:
            solve_range_0(normal, begin, end);
            break;
        case 1:
            solve_range_1(normal, begin, end);
            break;
        case 2:
            solve_range_2(normal, begin, end);
            break;
        default:
            throw std::domain_error("ncut must be 0, 1, or 2.");
            break;
    }
}


void SphereSlice::set_begin_end(int icut, double new_begin, double new_end) {
    if ((icut < 0) || (icut >= 2))
        throw std::domain_error("icut must be 0 or 1.");

    begin[icut] = new_begin;
    end[icut] = new_end;
}
