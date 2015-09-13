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
#include "celllists/vec3.h"


#define CHECK_ID(arg) if ((arg < 0) || (arg >= 3)) throw std::domain_error("arg must be 0, 1 or 2.")
#define UPDATE_BEGIN(found, work, dest) if (found) {if (work < dest) dest = work;} else {dest = work; found = true;}
#define UPDATE_END(found, work, dest)   if (found) {if (work > dest) dest = work;} else {dest = work; found = true;}


SphereSlice::SphereSlice(const double* center, const double* normals, double radius) :
    center(center), normals(normals), radius(radius) {

    if (radius <= 0) {
        throw std::domain_error("radius must be strictly positive.");
    }
    // Initialize variable data members
    cut_begin[0] = 0.0;
    cut_begin[1] = 0.0;
    cut_end[0] = 0.0;
    cut_end[1] = 0.0;
    // Compute from derived data members
    radius_sq = radius*radius;
    for (int i=0; i < 3; i++) {
        norms_sq[i] = vec3::normsq(normals + 3*i);
        norms[i] = sqrt(norms_sq[i]);
        frac_radii[i] = radius*norms[i];
        frac_center[i] = vec3::dot(center, normals + 3*i);
        vec3::copy(normals + 3*i, radius_normals + 3*i);
        vec3::iscale(radius_normals + 3*i, radius/norms[i]);
    }
}


void SphereSlice::solve_sphere(int id_axis, double &begin,
    double &end, double* point_begin, double* point_end) const {

    // TODO: everything in this function can be precomputed
    // Get the axis
    CHECK_ID(id_axis);
    // Find the ranges in fractional coordinates that encloses the cutoff
    begin = frac_center[id_axis] - frac_radii[id_axis];
    end = frac_center[id_axis] + frac_radii[id_axis];

    const double* radius_normal = radius_normals + 3*id_axis;
    if (point_begin != NULL) {
        point_begin[0] = center[0] - radius_normal[0];
        point_begin[1] = center[1] - radius_normal[1];
        point_begin[2] = center[2] - radius_normal[2];
    }
    if (point_end != NULL) {
        point_end[0] = center[0] + radius_normal[0];
        point_end[1] = center[1] + radius_normal[1];
        point_end[2] = center[2] + radius_normal[2];
    }
}

bool SphereSlice::solve_circle(int id_axis, int id_cut, double frac_cut,
    double &begin, double &end, double* point_begin, double* point_end) const {

    // Get the axis
    CHECK_ID(id_axis);
    const double* axis = normals + 3*id_axis;
    // Get the cut_normal
    CHECK_ID(id_cut);
    const double* cut_normal = normals + 3*id_cut;

    /* Define the parameters of a circle that is the intersection of
         - the sphere
         - the plane defined by cut_normal and cut: dot(r, cut_normal) = cut
     */
    // The difference in reduced coordinate between the center of the sphere
    // and the center of the circle.
    double delta_cut = frac_cut - frac_center[id_cut];
    // The amount lost from the total radius squared.
    double lost_radius_sq = delta_cut*delta_cut/norms_sq[id_cut];
    // The rest of the radius squared is for the size of the circle.
    double circle_radius_sq = radius_sq - lost_radius_sq;
    // Check if an intersecting circle exists, if not return false;
    if (circle_radius_sq < 0) return false;
    // Compute the circle radius
    double circle_radius = sqrt(circle_radius_sq);

    // Compute the center of the circle
    double circle_center[3];
    vec3::copy(center, circle_center);
    vec3::iadd(circle_center, cut_normal, delta_cut/norms_sq[id_cut]);

    /* Define a vector orthogonal to cut_normal, in the plane of axis and
       cut_normal. The length of the vector is such that, when added to the
       center of the circle, it just ends on the circle edge. The direction
       is chosen to either minimize or maximise the projection on axis. */
    double ortho[3]; // TODO: maybe precompute six of such ortho vectors
    // Copy of axis -> in plane of axis
    vec3::copy(axis, ortho);
    // Subtract projection on cut_normal
    //  -> in plane of axis and cut_normal
    //  -> orthogonal to cut_normal
    vec3::iadd(ortho, cut_normal, -vec3::dot(axis, cut_normal)/norms_sq[id_cut]);
    // Normalize to circle_radius
    vec3::iscale(ortho, circle_radius/vec3::norm(ortho));

    // Compute projection on axis, optionally compute points;
    if (point_begin==NULL) {
        begin = (circle_center[0] - ortho[0])*axis[0] +
                (circle_center[1] - ortho[1])*axis[1] +
                (circle_center[2] - ortho[2])*axis[2];
    } else {
        vec3::copy(circle_center, point_begin);
        vec3::iadd(point_begin, ortho, -1);
        begin = vec3::dot(point_begin, axis);
    }
    if (point_end==NULL) {
        end = (circle_center[0] + ortho[0])*axis[0] +
              (circle_center[1] + ortho[1])*axis[1] +
              (circle_center[2] + ortho[2])*axis[2];
    } else {
        vec3::copy(circle_center, point_end);
        vec3::iadd(point_end, ortho, +1);
        end = vec3::dot(point_end, axis);
    }

    return true;
}

bool SphereSlice::solve_line(int id_axis, int id_cut0, int id_cut1,
    double frac_cut0, double frac_cut1, double &begin, double &end,
    double* point_begin, double* point_end) const {
    throw std::logic_error("TODO");
}

void SphereSlice::solve_range_0(double &begin, double &end) const {
    // The first normal serves as axis on which begin and end is defined.
    solve_sphere(0, begin, end, NULL, NULL);
}


void SphereSlice::solve_range_1(double &begin, double &end) const {
    // work
    double work_begin;
    double work_end;
    double frac_tmp;
    double point_begin[3];
    double point_end[3];
    bool found_begin = false;
    bool found_end = false;

    // cut_normal
    const double* cut_normal = normals;

    // Whole-sphere solution
    solve_sphere(1, work_begin, work_end, point_begin, point_end);
    frac_tmp = vec3::dot(point_begin, cut_normal);
    if ((frac_tmp > cut_begin[0]) && (frac_tmp < cut_end[0])) {
        UPDATE_BEGIN(found_begin, work_begin, begin);
    }
    frac_tmp = vec3::dot(point_end, cut_normal);
    if ((frac_tmp > cut_begin[0]) && (frac_tmp < cut_end[0])) {
        UPDATE_END(found_end, work_end, end);
    }

    // Cut circle begin
    if (solve_circle(1, 0, cut_begin[0], work_begin, work_end, NULL, NULL)) {
        UPDATE_BEGIN(found_begin, work_begin, begin);
        UPDATE_END(found_end, work_end, end);
    }

    // Cut circle end
    if (solve_circle(1, 0, cut_end[0], work_begin, work_end, NULL, NULL)) {
        UPDATE_BEGIN(found_begin, work_begin, begin);
        UPDATE_END(found_end, work_end, end);
    }

    if (!(found_end && found_begin))
        throw std::logic_error("No solution found");
}


void SphereSlice::solve_range_2(double &begin, double &end) const{
    throw std::logic_error("TODO");
}


void SphereSlice::solve_range(int ncut, double &begin, double &end) const {
    switch (ncut) {
        case 0:
            solve_range_0(begin, end);
            break;
        case 1:
            solve_range_1(begin, end);
            break;
        case 2:
            solve_range_2(begin, end);
            break;
        default:
            throw std::domain_error("ncut must be 0, 1, or 2.");
            break;
    }
}


void SphereSlice::set_cut_begin_end(int icut, double new_begin, double new_end) {
    if ((icut < 0) || (icut >= 2))
        throw std::domain_error("icut must be 0 or 1.");
    if (new_begin >= new_end)
        throw std::domain_error("begin must be strictly smaller than end.");

    cut_begin[icut] = new_begin;
    cut_end[icut] = new_end;
}
