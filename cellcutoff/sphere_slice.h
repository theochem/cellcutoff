// CellCutoff is a library for periodic boundary conditions and real-space cutoff calculations.
// Copyright (C) 2017 The CellCutoff Development Team
//
// This file is part of CellCutoff.
//
// CellCutoff is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// CellCutoff is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
// --

/** @file

    This module implements the algorithm used by bars_cutoff to find all subcells
    that overlap with a given cutoff sphere. With hindsight, this might not be the
    simplest implementation, but it does the job relatively efficiently. Better algorithms
    should be possible in which even more gets precomputed in the Constructor of
    SphereSlice.

    The SphereSlice object can determine the begin and end of the cutoff sphere along one
    of the normals of the cell planes. It can also solve this problem for slices of the
    cutoff sphere, e.g. a slice contained between two parallel cutting planes, or a strip
    contained between two pairs of parallel cutting planes.

  */


#ifndef CELLCUTOFF_SPHERE_SLICE_H_
#define CELLCUTOFF_SPHERE_SLICE_H_

#include <stdexcept>
#include <string>


namespace cellcutoff {


/** @brief
        An exception for when solve_range(_*) can not find a solution.
 */
class no_solution_found : public std::domain_error {
 public:
  //! Create exception, default constructor
  explicit no_solution_found(const std::string& what_arg)
      : std::domain_error(what_arg) {}
};


/**
    @brief
        Describes a 0D, 1D or 2D slice of a cutoff sphere.

    In the 0D case, the slice is just the whole sphere.

    In the 1D case, a slice is just the part of a cutoff sphere between a pair of parallel
    cutting planes.

    In the the 2D, this object describes a strip from a cutoff sphere defined by two pairs
    of parallel cutting planes.

  */
class SphereSlice {
 public:
  /** @brief
          Create a SphereSlice

      This construcutor precomputes a lot of stuff that gets reused by the solve_range
      method, which is expected to be called many times.

      @params center
          The center of the cutoff sphere (len=3)

      @normals
          The normal vectors for the three sets of planes (len=9). These are typically the
          normals of the crystal cell planes.

      @radius
          The radius of the cutoff sphere
    */
  SphereSlice(const double* center, const double* normals, double radius);

  /* Move-constructor and assignment make no sense as the SphereSlice is
     constant after construction! Just pass a reference or a pointer instead. If needed,
     a copy can be made to guarantee that the data remains available. */
  //! Disable copy constructor
  SphereSlice(SphereSlice&&) = delete;
  //! Disable move constructor
  SphereSlice& operator=(const SphereSlice&) = delete;

  // Main API
  /** @brief Find the beginning and end of a (slice of) the cutoff sphere along a normal

      @param ncut
          The number of slices to consider (0, 1 or 2). 0 means no slices, and the
          beginning and end of the sphere along the normal 0 are found. 1 means the slice
          is defined by the pair of parallel cutting planes orthogonal to normal 0 and the
          begin and end along the second normal are found. 2 means the strip is defined by
          two pairs of parallel (within each set) cutting planes (normals 0 and 1) and the
          begin and end along normal 2 are found.

      @param begin
          Output argument for the begin of the range containing the sphere.

      @param end
          Output argument for the end of the range containing the sphere.
    */
  void solve_range(const int ncut, double* begin, double* end) const;
  /** @brief Set the positions of the cutting planes.

      @param icut
          The index of the plane to set (0, 1). This also refers to the normal
          vectors provided to the constructor.

      @param new_begin
          The position of the first plane along the axis defined by the corresponding
          normal.

      @param new_end
          The position of the second plane along the axis defined by the corresponding
          normal.
   */
  void set_cut_begin_end(const int icut, double new_begin, double new_end);

  // Auxiliary API, could also be useful and there is no need to really
  // make this private. Having it public also facilitates testing.

  //! Same is solve_range with ncut=0
  void solve_range_0(double* begin, double* end) const;
  //! Same is solve_range with ncut=1
  void solve_range_1(double* begin, double* end) const;
  //! Same is solve_range with ncut=2
  void solve_range_2(double* begin, double* end) const;

  /** @brief
          Solve the problem along an axis without imposing constraints planes. These are
          simply the extrema of the sphere along a normal.

      The solution will be rejected if the points on the sphere at the begin and end
      fall outside of any defined cutting planes (by default no other cutting planes
      are active).

      @param id_axis
          The axis (or normal) along which the begin and end must be found (0, 1, or 2).

      @param begin
          Output argument for solution.

      @param end
          Output argument for solution.

      @param id_cut0
          The index of the first set of cutting planes. Active if not negative. Solutions
          outside the planes are rejected.

      @param id_cut1
          The index of the second set of cutting planes. Active if not negative. Solutions
          outside the planes are rejected.
    */
  void solve_full(const int id_axis, double* begin, double* end,
      const int id_cut0 = -1, const int id_cut1 = -1) const;
  /** @brief Low-level function for solve_full, includes computation of points at begin and end.

      The points are returned as two extra output arguments, if no null-pointers are
      given. These can be useful to check if these points are valid solutions, i.e. not
      outside one or more pairs of parallel planes.
   */
  void solve_full_low(const int id_axis, double* begin, double* end,
      double* point_begin = nullptr, double* point_end = nullptr) const;

  /** @brief
          Solve the problem along an axis and imposing the solutions sit on a cutting
          plane.

      The solution will be rejected if the points on the sphere at the begin and end
      fall outside of the second pair of parallel planes.

      @param id_axis
          The axis along which the begin and end must be found (0, 1, or 2).

      @param frac_cut0
          The position of the cutting plane defining the constraint.

      @param begin
          Output argument for solution.

      @param end
          Output argument for solution.

      @param id_cut1
          The index of the "second" set of cutting planes. Active if not negative.
          Solutions outside the planes are rejected.
    */
  void solve_plane(const int id_axis, const int id_cut0, const double frac_cut0,
      double* begin, double* end, const int id_cut1 = -1) const;
  //! Low-level function for solve_plane, includes computation of points at begin and end.
  void solve_plane_low(const int id_axis, const int id_cut, const double frac_cut,
      double* begin, double* end,
      double* point_begin = nullptr, double* point_end = nullptr) const;

  /** @brief
          Solve the problem along an axis and imposing the solutions sit on the
          intersections of a two cutting planes.

      @param id_axis
          The axis along which the begin and end must be found (0, 1, or 2).

      @param frac_cut0
          The position of the first cutting plane defining the constraint.

      @param frac_cut1
          The position of the second cutting plane defining the constraint.

      @param begin
          Output argument for solution.

      @param end
          Output argument for solution.
    */
  void solve_line(const int id_axis, const int id_cut0, const int id_cut1,
      const double frac_cut0, const double frac_cut1, double* begin, double* end) const;
  //! Low-level function for solve_line, includes computation of points at begin and end.
  void solve_line_low(const int id_axis, const int id_cut0, const int id_cut1,
      const double frac_cut0, const double frac_cut1, double* begin, double* end,
      double* point_begin = nullptr, double* point_end = nullptr) const;
  /** @brief
          Compute the line defined by the intersection of two cutting planes.

      @param id_cut0
          Index for the first cutting plane, refers to one of the normals.

      @param id_cut1
          Index for the second cutting plane, refers to one of the normals.

      @param cut0
          The fractional coordinate for the position of the first cutting plane, relative
          to the center.

      @param id_cut1
          The fractional coordinate for the position of the second cutting plane, relative
          to the center.
    */
  double compute_plane_intersection(const int id_cut0, int const id_cut1,
      const double cut0, const double cut1, double* other_center) const;

  //! Test if a point lies between a pair of parallel cutting planes.
  bool inside_cuts(const int id_cut, const double* point) const;

 private:
  // Constant independent data members
  const double* center_;    //!< Center of the cutoff sphere
  const double* normals_;   //!< Normals of the cutting (or cell) planes
  const double radius_;     //!< Radius of the cutoff sphere

  // Configurable data members
  double cut_begin[2];  //!< Fractional coordinates defining the lower-bound cutting planes
  double cut_end[2];    //!< Fractional coordinates defining the upper-bound cutting planes

  // Derived from constant data members upon construction
  double radius_sq_;    //!< Square radius of the cutoff sphere
  double norms_sq_[3];  //!< Square norms of the normal vectors
  double norms_[3];     //!< Norms of the normal vectors
  double frac_radii_[3];    //!< Radius of the sphere projected on the normals
  double frac_center_[3];   //!< Center of the cutoff sphere in fractional coordinates
  double sphere_frac_begin_[3];   //!< Begin of the sphere in fractional coordinates
  double sphere_frac_end_[3];     //!< End of the sphere in fractional coordinates
  double radius_normals_[9];      //!< Vectors parallel to normals with length equal to radius
  double sphere_point_begin_[9];  //!< Position of beginnings of the sphere along each normal
  double sphere_point_end_[9];    //!< Position of endings of the sphere along each normal
  double dots_[9];        //!< Dot products between normal vectors
  double denoms_[9];      //!< Auxiliary used in denominator of intersection formula
  double cut_ortho_[27];  //!< Complicated. See comment in constructor.
};


/** @brief
        Auxiliary function for computing begin and end from typical intermeiate results.

    This is used by solve_*_low functions. It is present in the header to facilitate unit
    testing but it is intended for internal usage only.

    @param other_center
        The center of the cutoff sphere, or the center projected on the active
        constraints.

    @orhto
        A vector from the center of the cutoff sphere orthogonal to the constraints.

    @axis
        The axis along which the begin and end are to be found.

    @param begin
        Output argument for solution, begin of the sphere in fractional coordinates along
        the given axis.

    @param end
        Output argument for solution, end of the sphere in fractional coordinates along
        the given axis.

    @param point_begin
        Output argument for solution, position of the begin of the sphere

    @param point_end
        Output argument for solution, position of the end of the sphere.
 */
void compute_begin_end(const double* other_center, const double* ortho,
    const double* axis, double* begin, double* end,
    double* point_begin, double* point_end);

/** @brief
        Updates begin and end solution provided new candidate solutions

    work_begin will replace begin if it is not Nan and if it is lower than begin or when
    begin was NaN. Similar for end (lower->higher).

    @param work_begin
        A new candidate solution for begin. Nan of no valid solution was found.

    @param work_end
        A new candidate solution for begin. Nan of no valid solution was found.

    @param begin
        In/Out argument. Gets replaced by work_begin if it is a better solution.

    @param work_end
        In/Out argument. Gets replaced by work_end if it is a better solution.
  */
void update_begin_end(const double work_begin, const double work_end,
    double* begin, double* end);


}  // namespace cellcutoff


#endif  // CELLCUTOFF_SPHERE_SLICE_H_

// vim: textwidth=90 et ts=2 sw=2
