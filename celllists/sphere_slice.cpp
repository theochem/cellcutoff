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


#include "celllists/sphere_slice.h"

#include <cmath>
#include <memory>

#include "celllists/vec3.h"


namespace celllists {


#define CHECK_ID(ARG) \
    if ((ARG < 0) || (ARG >= 3)) throw std::domain_error(#ARG " must be 0, 1 or 2.")


SphereSlice::SphereSlice(const double* center, const double* normals, double radius)
    : center_(center), normals_(normals), radius_(radius) {
  // Check sanity of arguments
  if (radius_ <= 0)
    throw std::domain_error("radius must be strictly positive.");
  if (vec3::triple(normals_, normals_ + 3, normals_ + 6) == 0.0)
    throw std::domain_error("The three normals must be linearly independent.");
  // Initialize variable data members
  cut_begin[0] = 0.0;
  cut_begin[1] = 0.0;
  cut_end[0] = 0.0;
  cut_end[1] = 0.0;
  // Compute from derived data members
  radius_sq_ = radius_*radius_;
  for (int id_axis=0; id_axis < 3; ++id_axis) {
    const double* axis = normals_ + 3*id_axis;
    for (int id_cut=0; id_cut < 3; ++id_cut) {
      const double* cut_normal = normals_ + 3*id_cut;
      dots_[id_axis + 3*id_cut] = vec3::dot(axis, cut_normal);
    }
    norms_sq_[id_axis] = vec3::normsq(axis);
    norms_[id_axis] = sqrt(norms_sq_[id_axis]);
    frac_radii_[id_axis] = radius_*norms_[id_axis];
    frac_center_[id_axis] = vec3::dot(center_, axis);
    sphere_frac_begin_[id_axis] = frac_center_[id_axis] - frac_radii_[id_axis];
    sphere_frac_end_[id_axis] = frac_center_[id_axis] + frac_radii_[id_axis];
    vec3::copy(axis, radius_normals_ + 3*id_axis);
    vec3::iscale(radius_normals_ + 3*id_axis, radius_/norms_[id_axis]);
    vec3::copy(center_, sphere_point_begin_ + 3*id_axis);
    vec3::iadd(sphere_point_begin_ + 3*id_axis, radius_normals_ + 3*id_axis, -1);
    vec3::copy(center_, sphere_point_end_ + 3*id_axis);
    vec3::iadd(sphere_point_end_ + 3*id_axis, radius_normals_ + 3*id_axis);
  }
  for (int id_axis=0; id_axis < 3; ++id_axis) {
    for (int id_cut=0; id_cut < 3; ++id_cut) {
      denoms[id_axis + 3*id_cut] = (
          dots_[id_axis + 3*id_cut]*dots_[id_axis + 3*id_cut] -
          dots_[id_axis + 3*id_axis]*dots_[id_cut + 3*id_cut]);
    }
  }
  for (int id_axis=0; id_axis < 3; ++id_axis) {
    const double* axis = normals_ + 3*id_axis;
    for (int id_cut=0; id_cut < 3; ++id_cut) {
      /* Define a vector orthogonal to cut_normal, in the plane of axis
         and cut_normal. The length of the vector is such that, when added
         to the center of the circle, it just ends on the circle edge. The
         direction is chosen to either minimize or maximise the projection
         on axis. */
      const double* cut_normal = normals_ + 3*id_cut;
      double ortho[3] = {0.0, 0.0, 0.0};
      if (id_cut != id_axis) {
        // Copy of axis -> in plane of axis
        vec3::copy(axis, ortho);
        // Subtract projection on cut_normal
        //  -> in plane of axis and cut_normal
        //  -> orthogonal to cut_normal
        vec3::iadd(ortho, cut_normal,
                   -vec3::dot(axis, cut_normal)/norms_sq_[id_cut]);
        // Normalize
        vec3::iscale(ortho, 1.0/vec3::norm(ortho));
      }
      // Store solution
      vec3::copy(ortho, cut_ortho + 3*id_axis + 9*id_cut);
    }
  }
}


void SphereSlice::solve_range(const int ncut, double* begin, double* end) const {
  switch (ncut) {
    case 0: {
      solve_range_0(begin, end);
      break;
    }
    case 1: {
      solve_range_1(begin, end);
      break;
    }
    case 2: {
      solve_range_2(begin, end);
      break;
    }
    default: {
      throw std::domain_error("ncut must be 0, 1, or 2.");
    }
  }
}


void SphereSlice::set_cut_begin_end(const int icut, double new_begin, double new_end) {
  if ((icut < 0) || (icut >= 2))
    throw std::domain_error("icut must be 0 or 1.");
  if (new_begin >= new_end)
    throw std::domain_error("begin must be strictly smaller than end.");

  cut_begin[icut] = new_begin;
  cut_end[icut] = new_end;
}


void SphereSlice::solve_range_0(double* begin, double* end) const {
  // Start out with NaNs in begin and end, as to indicate that they are not
  // found yet.
  *begin = NAN;
  *end = NAN;
  solve_full(0, begin, end);
}


void SphereSlice::solve_range_1(double* begin, double* end) const {
  // Start out with NaNs in begin and end, as to identify the unsolvable case.
  // (See end of function.)
  *begin = NAN;
  *end = NAN;

  // Find solutions for all sensible combinations of constraints
  // * Case A: cut_begin[0]
  solve_plane(1, 0, cut_begin[0], begin, end);
  // * Case B: cut_end[0]
  solve_plane(1, 0, cut_end[0], begin, end);
  // * Case C: Whole-sphere solution, only accepted when within cut0_normal bounds
  solve_full(1, begin, end, 0, -1);

  if (std::isnan(*begin) || std::isnan(*end))
    throw no_solution_found("No solution found in solve_range_1.");
}


void SphereSlice::solve_range_2(double* begin, double* end) const {
  // Start out with NaNs in begin and end, as to identify the unsolvable case.
  // (See end of function.)
  *begin = NAN;
  *end = NAN;

  // Find solutions for all sensible combinations of constraints
  // * Case A: cut_begin[0], cut_begin[1]
  solve_line(2, 0, 1, cut_begin[0], cut_begin[1], begin, end);
  // * Case B: cut_begin[0], cut_end[1]
  solve_line(2, 0, 1, cut_begin[0], cut_end[1], begin, end);
  // * Case C: cut_end[0], cut_begin[1]
  solve_line(2, 0, 1, cut_end[0], cut_begin[1], begin, end);
  // * Case D: cut_end[0], cut_end[1]
  solve_line(2, 0, 1, cut_end[0], cut_end[1], begin, end);
  // * Case E: cut_begin[0]
  solve_plane(2, 0, cut_begin[0], begin, end, 1);
  // * Case F: cut_end[0]
  solve_plane(2, 0, cut_end[0], begin, end, 1);
  // * Case G: cut_begin[1]
  solve_plane(2, 1, cut_begin[1], begin, end, 0);
  // * Case H: cut_end[1]
  solve_plane(2, 1, cut_end[1], begin, end, 0);
  // * Case I: Whole-sphere solution, only if within cut0_normal and cut1_normal bounds
  solve_full(2, begin, end, 0, 1);

  if (std::isnan(*begin) || std::isnan(*end))
    throw no_solution_found("No solution found in solve_range_2.");
}


void SphereSlice::solve_full(const int id_axis, double* begin, double* end,
    const int id_cut0, const int id_cut1) const {
  double work_begin, work_end;
  if ((id_cut0 == -1) && (id_cut1 == -1)) {
    solve_full_low(id_axis, &work_begin, &work_end);
  } else {
    double point_begin[3];
    double point_end[3];
    solve_full_low(id_axis, &work_begin, &work_end, point_begin, point_end);
    // Reject solution if not between cut0 and/or cut1 planes
    if (!inside_cuts(id_cut0, point_begin))
      work_begin = NAN;
    if (!inside_cuts(id_cut1, point_begin))
      work_begin = NAN;
    if (!inside_cuts(id_cut0, point_end))
      work_end = NAN;
    if (!inside_cuts(id_cut1, point_end))
      work_end = NAN;
  }
  update_begin_end(work_begin, work_end, begin, end);
}


void SphereSlice::solve_full_low(const int id_axis, double* begin,
    double* end, double* point_begin, double* point_end) const {
  // Check the axis
  CHECK_ID(id_axis);
  // Everything is precomputed...
  *begin = sphere_frac_begin_[id_axis];
  *end = sphere_frac_end_[id_axis];
  if (point_begin != nullptr)
    vec3::copy(sphere_point_begin_ + 3*id_axis, point_begin);
  if (point_end != nullptr)
    vec3::copy(sphere_point_end_ + 3*id_axis, point_end);
}


void SphereSlice::solve_plane(const int id_axis, int const id_cut0,
    double const frac_cut0, double* begin, double* end, const int id_cut1) const {
  double work_begin, work_end;
  if (id_cut1 == -1) {
    solve_plane_low(id_axis, id_cut0, frac_cut0, &work_begin, &work_end);
  } else {
    double point_begin[3];
    double point_end[3];
    solve_plane_low(id_axis, id_cut0, frac_cut0, &work_begin, &work_end,
            point_begin, point_end);
    // Reject solution if not between cut1 planes
    if (std::isfinite(work_begin)) {
      if (!inside_cuts(id_cut1, point_begin))
        work_begin = NAN;
    }
    if (std::isfinite(work_end)) {
      if (!inside_cuts(id_cut1, point_end))
        work_end = NAN;
    }
  }
  update_begin_end(work_begin, work_end, begin, end);
}


void SphereSlice::solve_plane_low(const int id_axis, const int id_cut,
    const double frac_cut, double* begin, double* end,
    double* point_begin, double* point_end) const {
  // Get the axis
  CHECK_ID(id_axis);
  const double* axis = normals_ + 3*id_axis;
  // Get the cut_normal
  CHECK_ID(id_cut);
  const double* cut_normal = normals_ + 3*id_cut;

  /* Define the parameters of a circle that is the intersection of
       - the sphere
       - the plane defined by cut_normal and cut: dot(r, cut_normal) = cut
   */
  // The difference in reduced coordinate between the center of the sphere
  // and the center of the circle.
  double delta_cut = frac_cut - frac_center_[id_cut];
  // The amount lost from the total radius squared.
  double lost_radius_sq = delta_cut*delta_cut/norms_sq_[id_cut];
  // The rest of the radius squared is for the size of the circle.
  double circle_radius_sq = radius_sq_ - lost_radius_sq;
  // Check if an intersecting circle exists, if not return;
  if (circle_radius_sq < 0) {
    *begin = NAN;
    *end = NAN;
    return;
  }
  // Compute the circle radius
  double circle_radius = sqrt(circle_radius_sq);

  // Compute the center of the circle
  double circle_center[3];
  vec3::copy(center_, circle_center);
  vec3::iadd(circle_center, cut_normal, delta_cut/norms_sq_[id_cut]);

  // Get a vector orthogonal to cut_normal, in the plane of axis;
  double ortho[3];
  vec3::copy(cut_ortho + 3*id_axis + 9*id_cut, ortho);
  // Normalize to circle_radius
  vec3::iscale(ortho, circle_radius);
  // Compute projection on axis of two solutions, optionally compute points;
  compute_begin_end(circle_center, ortho, axis, begin, end, point_begin, point_end);
}


void SphereSlice::solve_line(const int id_axis, const int id_cut0, const int id_cut1,
    const double frac_cut0, const double frac_cut1, double* begin, double* end) const {

  double work_begin, work_end;
  solve_line_low(id_axis, id_cut0, id_cut1, frac_cut0, frac_cut1,
                 &work_begin, &work_end);
  update_begin_end(work_begin, work_end, begin, end);
}


void SphereSlice::solve_line_low(const int id_axis, const int id_cut0, const int id_cut1,
    const double frac_cut0, const double frac_cut1, double* begin, double* end,
    double* point_begin, double* point_end) const {

  // Run some checks on the ID arguments.
  CHECK_ID(id_axis);
  CHECK_ID(id_cut0);
  CHECK_ID(id_cut1);

  // Select the vectors
  const double* axis = normals_ + 3*id_axis;
  const double* cut0_normal = normals_ + 3*id_cut0;
  const double* cut1_normal = normals_ + 3*id_cut1;

  // Cuts relative to the center
  double delta_cut0 = frac_cut0 - frac_center_[id_cut0];
  double delta_cut1 = frac_cut1 - frac_center_[id_cut1];

  double line_center[3];
  double lost_radius_sq = compute_plane_intersection(id_cut0, id_cut1,
      delta_cut0, delta_cut1, line_center);
  vec3::iadd(line_center, center_);

  // Compute the remaining line radius
  double line_radius_sq = radius_sq_ - lost_radius_sq;
  if (line_radius_sq < 0) {
    *begin = NAN;
    *end = NAN;
    return;
  }
  double line_radius = sqrt(line_radius_sq);

  // Compute the basis vector (easy).
  double basis[3];
  vec3::cross(cut0_normal, cut1_normal, basis);
  double scale = line_radius/vec3::norm(basis);
  if (vec3::dot(axis, basis) < 0) scale *= -1;
  vec3::iscale(basis, scale);

  // Compute projection on axis, optionally compute points;
  compute_begin_end(line_center, basis, axis, begin, end, point_begin, point_end);
}


double SphereSlice::compute_plane_intersection(const int id_cut0, const int id_cut1,
    const double cut0, const double cut1, double* other_center) const {

  CHECK_ID(id_cut0);
  CHECK_ID(id_cut1);

  // Select the vectors
  const double* cut0_normal = normals_ + 3*id_cut0;
  const double* cut1_normal = normals_ + 3*id_cut1;

  // Find the nearest point where the two planes cross
  double dot00 = norms_sq_[id_cut0];
  double dot01 = dots_[id_cut0 + 3*id_cut1];
  double dot11 = norms_sq_[id_cut1];
  double denom = denoms[id_cut0 + 3*id_cut1];
  double ratio0 = (cut1*dot01 - cut0*dot11)/denom;
  double ratio1 = (cut0*dot01 - cut1*dot00)/denom;
  if (other_center != nullptr) {
    other_center[0] = cut0_normal[0]*ratio0 + cut1_normal[0]*ratio1;
    other_center[1] = cut0_normal[1]*ratio0 + cut1_normal[1]*ratio1;
    other_center[2] = cut0_normal[2]*ratio0 + cut1_normal[2]*ratio1;
  }

  // Compute the distance squared from the origin to the nearest point on the
  // intersection.
  return ratio0*ratio0*dot00 + 2*ratio0*ratio1*dot01 + ratio1*ratio1*dot11;
}


bool SphereSlice::inside_cuts(const int id_cut, const double* point) const {
  // if id_cut == -1, the test always passes, i.e. bounds are not imposed.
  if (id_cut == -1) return true;
  CHECK_ID(id_cut);
  const double* cut_normal = normals_ + 3*id_cut;
  double frac_cut = vec3::dot(point, cut_normal);
  return (frac_cut > cut_begin[id_cut]) && (frac_cut < cut_end[id_cut]);
}


void compute_begin_end(const double* other_center, const double* ortho,
    const double* axis, double* begin, double* end,
    double* point_begin, double* point_end) {
  // Compute projection on axis, optionally compute points;
  if (point_begin == nullptr) {
    *begin = (other_center[0] - ortho[0])*axis[0] +
             (other_center[1] - ortho[1])*axis[1] +
             (other_center[2] - ortho[2])*axis[2];
  } else {
    vec3::copy(other_center, point_begin);
    vec3::iadd(point_begin, ortho, -1);
    *begin = vec3::dot(point_begin, axis);
  }
  if (point_end == nullptr) {
    *end = (other_center[0] + ortho[0])*axis[0] +
           (other_center[1] + ortho[1])*axis[1] +
           (other_center[2] + ortho[2])*axis[2];
  } else {
    vec3::copy(other_center, point_end);
    vec3::iadd(point_end, ortho, 1);
    *end = vec3::dot(point_end, axis);
  }
}


void update_begin_end(const double work_begin, const double work_end,
    double* begin, double* end) {
  if (!std::isnan(work_begin)) {
    if (std::isnan(*begin)) {
      *begin = work_begin;
    } else if (work_begin < *begin) {
      *begin = work_begin;
    }
  }
  if (!std::isnan(work_end)) {
    if (std::isnan(*end)) {
      *end = work_end;
    } else if (work_end > *end) {
      *end = work_end;
    }
  }
}


}  // namespace celllists

// vim: textwidth=90 et ts=2 sw=2
