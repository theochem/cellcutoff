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
#include <vector>
#include "celllists/cell.h"
#include "celllists/vec3.h"


Cell::Cell(const double* _rvecs, int _nvec): nvec(_nvec) {
    // check if nvec is sensible
    if ((_nvec < 0) || (_nvec > 3)) {
        throw std::domain_error("The number of cell vectors must be 0, 1, 2 or 3.");
    }

    // copy the given _rvecs and _nvec:
    for (int ivec=nvec*3-1; ivec>=0; ivec--) {
        rvecs[ivec] = _rvecs[ivec];
    }

    // compute the volume
    switch(nvec) {
        case 0:
            volume = 0.0;
            break;
        case 1:
            volume = vec3::norm(rvecs);
            break;
        case 2:
            double tmp;
            tmp = vec3::dot(rvecs, rvecs+3);
            tmp = vec3::normsq(rvecs)*vec3::normsq(rvecs+3) - tmp*tmp;
            if (tmp > 0) {
                volume = sqrt(tmp);
            } else {
                volume = 0.0;
            }
            break;
        case 3:
            volume = fabs(vec3::triple(rvecs, rvecs+3, rvecs+6));
            break;
    }

    // If the volume is zero and nvec > 0, raise an error. In this case, the
    // reciprocal cell vectors can not be computed.
    if ((volume == 0.0) && (nvec > 0)) {
        throw singular_cell_vectors("The cell vectors are degenerate");
    }

    // complete the list of rvecs in case nvec < 3
    switch(nvec) {
        case 0:
            // Just put in the identity matrix.
            rvecs[0] = 1.0;
            rvecs[1] = 0.0;
            rvecs[2] = 0.0;
            rvecs[3] = 0.0;
            rvecs[4] = 1.0;
            rvecs[5] = 0.0;
            rvecs[6] = 0.0;
            rvecs[7] = 0.0;
            rvecs[8] = 1.0;
            break;
        case 1: {
            // Add two rvecs that are orthogonal to the given rvec, orthogonal
            // to each other and normalized. The three vectors will be
            // right-handed.
            // 1) find the component of the given vector with the smallest
            // absolute value
            int ismall = 0;
            if (fabs(rvecs[1]) < fabs(rvecs[0])) {
                ismall = 1;
                if (fabs(rvecs[2]) <= fabs(rvecs[1])) {
                    ismall = 2;
                }
            } else if (fabs(rvecs[2]) < fabs(rvecs[0])) {
                ismall = 2;
            }
            // 2) store a temporary vector in position 3
            rvecs[6] = 0.0;
            rvecs[7] = 0.0;
            rvecs[8] = 0.0;
            rvecs[ismall+6] = 1.0;
            // 3) compute the cross product of vector 3 and 1
            vec3::cross(rvecs+6, rvecs, rvecs+3);
            // 4) normalize
            double norm = vec3::norm(rvecs+3);
            vec3::iscale(rvecs+3, 1.0/norm);
            // the rest is done in case 2, so no break here!
        }
        case 2:
            // Add one rvec that is normalized and orthogonal to the two given
            // rvecs. The three vectors will be right-handed.
            // 1) compute the cross product of vector 1 and 2
            vec3::cross(rvecs, rvecs+3, rvecs+6);
            // 2) normalize
            double norm = vec3::norm(rvecs+6);
            vec3::iscale(rvecs+6, 1.0/norm);
    }

    // Now we assume that rvecs contains a set of three well-behaved
    // non-degenerate vectors. Cramer's rule is used to compute the reciprocal
    // space vectors. This is fairly ugly in terms of numerical stability but
    // it keeps things simple.
    vec3::cross(rvecs+3, rvecs+6, gvecs);
    vec3::cross(rvecs+6, rvecs, gvecs+3);
    vec3::cross(rvecs, rvecs+3, gvecs+6);
    // inverse of determinant
    double denom = 1.0/vec3::dot(gvecs, rvecs);
    // inverse
    vec3::iscale(gvecs, denom);
    vec3::iscale(gvecs+3, denom);
    vec3::iscale(gvecs+6, denom);

    // compute the spacings and the lengths of the cell vectors
    for (int ivec=2; ivec>=0; ivec--) {
        rlengths[ivec] = vec3::norm(rvecs+3*ivec);
        glengths[ivec] = vec3::norm(gvecs+3*ivec);
        rspacings[ivec] = 1.0/glengths[ivec];
        gspacings[ivec] = 1.0/rlengths[ivec];
    }
}


void Cell::wrap(double* delta) const {
    // Wrap the relative vector back into the cell in the range ]-0.5, 0.5].
    double x;
    if (nvec == 0) return;
    // Compute the first fractional coordinates, subtract one half and ceil. The round
    // founction is intentionally not used here! The half-ways case is always up instead
    // of away from zero.
    x = ceil(vec3::dot(gvecs, delta) - 0.5);
    vec3::iadd(delta, rvecs, -x);
    if (nvec == 1) return;
    // Compute the second fractional coordinates, subtract one half and ceil.
    x = ceil(vec3::dot(gvecs+3, delta) - 0.5);
    vec3::iadd(delta, rvecs+3, -x);
    if (nvec == 2) return;
    // Compute the third fractional coordinates, subtract one half and ceil.
    x = ceil(vec3::dot(gvecs+6, delta) - 0.5);
    vec3::iadd(delta, rvecs+6, -x);
}


void Cell::to_rfrac(const double* rcart, double* rfrac) const {
    // Transfrom to real-space fractional coordinates
    vec3::matvec(gvecs, rcart, rfrac);
}


void Cell::to_rcart(const double* rfrac, double* rcart) const {
    // Transfrom to real-space Cartesian coordinates
    vec3::tmatvec(rvecs, rfrac, rcart);
}


void Cell::to_gfrac(const double* gcart, double* gfrac) const {
    // Transform to reciprocal space Cartesian coordinates
    vec3::matvec(rvecs, gcart, gfrac);
}


void Cell::to_gcart(const double* gfrac, double* gcart) const {
    // Transform to reciprocal-space Cartesian coordinates
    vec3::tmatvec(gvecs, gfrac, gcart);
}


void Cell::add_rvec(double* delta, const int* coeffs) const {
    // Simply adds an linear combination of real cell vectors to delta.
    if (nvec == 0) return;
    vec3::iadd(delta, rvecs, coeffs[0]);
    if (nvec == 1) return;
    vec3::iadd(delta, rvecs+3, coeffs[1]);
    if (nvec == 2) return;
    vec3::iadd(delta, rvecs+6, coeffs[2]);
}


double Cell::get_rvec(int ivec, int icomp) const {
    if ((ivec < 0) || (ivec >= 3)) {
        throw std::domain_error("ivec must be 0, 1 or 2.");
    }
    if ((icomp < 0) || (icomp >= 3)) {
        throw std::domain_error("icomp must be 0, 1 or 2.");
    }
    return rvecs[3*ivec + icomp];
}


double Cell::get_gvec(int ivec, int icomp) const {
    if ((ivec < 0) || (ivec >= 3)) {
        throw std::domain_error("ivec must be 0, 1 or 2.");
    }
    if ((icomp < 0) || (icomp >= 3)) {
        throw std::domain_error("icomp must be 0, 1 or 2.");
    }
    return gvecs[3*ivec + icomp];
}


double Cell::get_rlength(int ivec) const {
    if ((ivec < 0) || (ivec >= 3)) {
        throw std::domain_error("ivec must be 0, 1 or 2.");
    }
    return rlengths[ivec];
}


double Cell::get_glength(int ivec) const {
    if ((ivec < 0) || (ivec >= 3)) {
        throw std::domain_error("ivec must be 0, 1 or 2.");
    }
    return glengths[ivec];
}


double Cell::get_rspacing(int ivec) const {
    if ((ivec < 0) || (ivec >= 3)) {
        throw std::domain_error("ivec must be 0, 1 or 2.");
    }
    return rspacings[ivec];
}


double Cell::get_gspacing(int ivec) const {
    if ((ivec < 0) || (ivec >= 3)) {
        throw std::domain_error("ivec must be 0, 1 or 2.");
    }
    return gspacings[ivec];
}


bool Cell::is_cubic() const {
    if (!is_cuboid()) return false;
    if (nvec < 2) return true;
    if (rvecs[0] != rvecs[4]) return false;
    if (nvec < 3) return true;
    if (rvecs[0] != rvecs[8]) return false;
    return true;
}


bool Cell::is_cuboid() const {
    if (nvec < 1) return true;
    if (rvecs[1] != 0.0) return false;
    if (rvecs[2] != 0.0) return false;
    if (nvec < 2) return true;
    if (rvecs[3] != 0.0) return false;
    if (rvecs[5] != 0.0) return false;
    if (nvec < 3) return true;
    if (rvecs[6] != 0.0) return false;
    if (rvecs[7] != 0.0) return false;
    return true;
}


int Cell::set_ranges_rcut(const double* center, double rcut, int* ranges_begin,
    int* ranges_end) const {
    if (rcut <= 0) {
        throw std::domain_error("rcut must be strictly positive.");
    }
    double frac[3];
    int ncell = 1;
    to_rfrac(center, frac);
    for (int ivec=nvec-1; ivec>=0; ivec--) {
        // Use spacings between planes to find first plane before cutoff sphere and last
        // plane after cutoff sphere. To this end, we must divide rcut by the spacing
        // between planes.
        double frac_rcut = rcut/rspacings[ivec];
        ranges_begin[ivec] = static_cast<int>(floor(frac[ivec]-frac_rcut));
        ranges_end[ivec] = static_cast<int>(ceil(frac[ivec]+frac_rcut));
        ncell *= (ranges_end[ivec] - ranges_begin[ivec]);
    }
    return ncell;
}


void Cell::select_inside_low(SphereSlice* slice, const int* shape,
    const bool* pbc, std::vector<int> &bars, int* prefix, int ivec) const {

    double begin_exact = 0.0;
    double end_exact = 0.0;
    // Solve the hard problem elsewhere.
    slice->solve_range(ivec, begin_exact, end_exact);
    int begin = static_cast<int>(floor(begin_exact));
    int end = static_cast<int>(ceil(end_exact));
    // Truncate this range if there are non-periodic bounds
    if (!pbc[ivec]) {
        if (begin < 0) begin = 0;
        if (end > shape[ivec]) end = shape[ivec];
    }

    if (ivec == nvec - 1) {
        // If we are dealing with the last recursion, just store the bar.
        for (int i=0; i<ivec; i++)
            bars.push_back(prefix[i]);
        bars.push_back(begin);
        bars.push_back(end);
    } else {
        // If this is not yet the last recursion, iterate over the range of integer
        // fractional coordinates, and go one recursion deeper in each iteration.
        for (int i = begin; i < end; i++) {
            // Make sure the following recursion knows the indices of the current bar.
            prefix[ivec] = i;
            // Make a new cut in the spere slice.
            slice->set_cut_begin_end(ivec, i, i+1);
            // Make recursion
            select_inside_low(slice, shape, pbc, bars, prefix, ivec+1);
        }
    }
}


size_t Cell::select_inside_rcut(const double* center, double rcut,
    const int* shape, const bool* pbc, std::vector<int> &bars) const {
    if (nvec == 0) {
        throw std::domain_error("The cell must be at least 1D periodic for select_inside_rcut.");
    }
    if (rcut <= 0) {
        throw std::domain_error("rcut must be strictly positive.");
    }

    // For all the heavy work, precomputes a lot.
    SphereSlice sphere_slice(center, gvecs, rcut);

    // Prefix is used to keep track of current bar indices while going into recursion.
    int prefix[nvec-1];
    // Compute bars and return the number of bars
    select_inside_low(&sphere_slice, shape, pbc, bars, prefix, 0);
    return bars.size()/(nvec+1);
}


int smart_wrap(int i, int shape, bool pbc) {
    if (pbc) {
        i %= shape;
        if (i < 0) i += shape;
        return i;
    } else if (i < 0) {
        return -1;
    } else if (i >= shape) {
        return -1;
    } else {
        return i;
    }
}
