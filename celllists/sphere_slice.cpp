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


SphereSlice::SphereSlice(const double* center, const double* normals, double radius) :
    center(center), normals(normals), radius(radius) {

    if (radius <= 0) {
        throw std::domain_error("radius must be strictly positive.");
    }

    cut_begin[0] = 0.0;
    cut_begin[1] = 0.0;
    cut_end[0] = 0.0;
    cut_end[1] = 0.0;

}


void SphereSlice::solve_sphere(const double* axis, double &begin,
    double &end, double* point_begin, double* point_end) const {

    // Convert the radius to reduced coordinates
    double axis_norm = vec3::norm(axis);
    double proj_radius = radius*axis_norm;
    // Convert the center to reduced coordinates
    double proj_center = vec3::dot(center, axis);
    // Find the ranges in fractional coordinates that encloses the cutoff
    begin = proj_center - proj_radius;
    end = proj_center + proj_radius;

    if (point_begin != NULL) {
        point_begin[0] = center[0] - radius*axis[0]/axis_norm;
        point_begin[1] = center[1] - radius*axis[1]/axis_norm;
        point_begin[2] = center[2] - radius*axis[2]/axis_norm;
    }
    if (point_end != NULL) {
        point_end[0] = center[0] + radius*axis[0]/axis_norm;
        point_end[1] = center[1] + radius*axis[1]/axis_norm;
        point_end[2] = center[2] + radius*axis[2]/axis_norm;
    }
}

bool SphereSlice::solve_circle(const double* axis, const double* cut_normal,
    double cut, double &begin, double &end, double* point_begin,
    double* point_end) const {

    // Compute the circle radius
    double delta_cut = cut - vec3::dot(center, cut_normal);
    double cut_normal_sq = vec3::normsq(cut_normal);
    double lost_radius_sq = delta_cut*delta_cut/cut_normal_sq;
    double circle_radius_sq = radius*radius - lost_radius_sq;
    // Check if an intersection circle exists, if not return false;
    if (circle_radius_sq < 0) return false;
    double circle_radius = sqrt(circle_radius_sq);

    // Compute the center of the circle
    double circle_center[3];
    vec3::copy(center, circle_center);
    vec3::iadd(circle_center, cut_normal, delta_cut/cut_normal_sq);

    // Then add a vector orthogonal to normal, in the plane of axis and normal
    // that brings us to the surface of the sphere
    double ortho[3];
    vec3::copy(axis, ortho);
    vec3::iadd(ortho, cut_normal, -vec3::dot(axis, cut_normal)/cut_normal_sq);
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

void SphereSlice::solve_range_0(double &begin, double &end) const {
    // The first normal serves as axis on which begin and end is defined.
    solve_sphere(normals, begin, end, NULL, NULL);
}


#define UPDATE_BEGIN(found, work, dest) if (found) {if (work < dest) dest = work;} else {dest = work; found = true;}
#define UPDATE_END(found, work, dest)   if (found) {if (work > dest) dest = work;} else {dest = work; found = true;}


void SphereSlice::solve_range_1(double &begin, double &end) const {
    // The first normal is used to define the slice.
    const double* cut_normal = normals;
    // The second normal serves as axis on which begin and end is defined.
    const double* axis = normals + 3;
    // work
    double work_begin;
    double work_end;
    double cut_tmp;
    double point_begin[3];
    double point_end[3];
    bool found_begin = false;
    bool found_end = false;

    // Whole-sphere solution
    solve_sphere(axis, work_begin, work_end, point_begin, point_end);
    cut_tmp = vec3::dot(point_begin, cut_normal);
    if ((cut_tmp > cut_begin[0]) && (cut_tmp < cut_end[0])) {
        UPDATE_BEGIN(found_begin, work_begin, begin);
    }
    cut_tmp = vec3::dot(point_end, cut_normal);
    if ((cut_tmp > cut_begin[0]) && (cut_tmp < cut_end[0])) {
        UPDATE_END(found_end, work_end, end);
    }

    // Cut begin
    if (solve_circle(axis, cut_normal, cut_begin[0], work_begin, work_end, NULL, NULL)) {
        UPDATE_BEGIN(found_begin, work_begin, begin);
        UPDATE_END(found_end, work_end, end);
    }

    // Cut end
    if (solve_circle(axis, cut_normal, cut_end[0], work_begin, work_end, NULL, NULL)) {
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
