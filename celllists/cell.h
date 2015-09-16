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


#ifndef CELLLISTS_CELL_H_
#define CELLLISTS_CELL_H_

#include <exception>
#include <string>
#include <vector>

#include "celllists/sphere_slice.h"


namespace celllists {


/** @brief
        An exception for singular cell vectors.
 */
class singular_cell_vectors : public std::domain_error {
 public:
  explicit singular_cell_vectors(const std::string& what_arg)
      : std::domain_error(what_arg) {}
};


/** @brief
        3D/2D/1D/0D cell and derived quantities in a 3D space.

    Upon construction, an object of this class acts as a read-only representation of the
    cell. Reciprocal cell vectors, vector lengths and spacings between planes are computed
    immediately. All sorts of manipulations on fractional or Cartesian coordinates are
    supported.

    Even though lower-dimensional periodic boundary conditions are supported, this class
    is specific for 3D systems. In case of 0D, 1D or 2D PBC, the cell vectors are
    internally extended with orthonormal basis vectors to guarantee an invertible
    transformation between Cartesian and fractional coordinates. In that case, the
    fractional coordinates are actually also Cartesian coordinates in directions
    orthogonal to the available cell vectors. The extra basis vectors are always such that
    the complete set of vectors is right-handed.
 */
class Cell {
 public:
  /** @brief
          Construct a Cell object.

      @param _rvecs
          A pointer to `3*nvec` doubles that represent the Cartesian real-space vectors in
          row-major ordering. The vectors may not have a linear dependency. Each vector is
          one row in a 3x3 matrix.

      @param _nvec
          The number of cell vectors. This corresponds to the dimensionality of the cell.
          `nvec` must be 0, 1, 2 or 3.
  */
  Cell(const double* rvecs, const int nvec);

  // Copy-constructor, move-constructor and assignment make no sense as the Cell is
  // constant after construction! Just pass references or pointers instead.
  Cell(const Cell& that) = delete;
  Cell(Cell&&) = delete;
  Cell& operator=(const Cell&) = delete;

  /** @brief
          Create a Cell object with the reciprocal cell. The caller owns the pointer. This
          practically equivalent to `Cell(cell.gvecs(), cell.nvec())`, but more efficient
          and without precision loss.
   */
  Cell* create_reciprocal() const;

  //! Returns the number of periodic dimensions.
  int nvec() const { return nvec_; }

  //! Returns all real-space vectors.
  const double* rvecs() const { return rvecs_; }

  //! Returns a real-space vector.
  const double* rvec(const int ivec) const;

  //! Returns all reciprocal-space vectors.
  const double* gvecs() const { return gvecs_; }

  //! Returns a reciprocal-space vector.
  const double* gvec(const int ivec) const;

  //! Returns the volume (or area or length) of the cell.
  double volume() const { return volume_; }

  //! Returns the volume (or area or length) of the reciprocal cell.
  double gvolume() const { return gvolume_; }

  //! Returns the lengths of the real-space vectors.
  const double* rlengths() const { return rlengths_; }

  //! Returns the lengths of the reciprocal-space vectors.
  const double* glengths() const { return glengths_; }

  //! Returns the spacings between the real-space crystal plane
  const double* rspacings() const { return rspacings_; }

  //! Returns the spacings between the reciprocal-space crystal plane
  const double* gspacings() const { return gspacings_; }


  /** @brief
          Test if cell is cubic

      The cell must also be aligned with Cartesian axes, i.e a to x, b to y and c to z. No
      small errors allowed.
    */
  bool cubic() const;


  /** @brief
          Test if cell is cuboid (orthorombic)

      The cell must also be aligned with Cartesian axes, i.e a to x, b to y and c to z. No
      small errors allowed.
    */
  bool cuboid() const;


  /** @brief
          Convert Cartesian real-space coordinates to fractional.

      @param cart
          A pointer to 3 doubles containing the input real-space vector.

      @param frac
          A pointer to 3 doubles in which the output is written.

      This effectively computes the dot products of the Cartesian vector with the
      reciprocal-space cell vectors.
   */
  void to_rfrac(const double* cart, double* frac) const;


  /** @brief
          Convert fractional real-space coordinates to Cartesian.

      @param frac
          A pointer to 3 doubles containing the input fractional coordinates.

      @param cart
          A pointer to 3 doubles to which the output is written

      This effectively computes a linear combination of real-space cell vectors.
   */
  void to_rcart(const double* frac, double* cart) const;


  /** @brief
          In-place wrap a (relative) vector back into the cell ]-0.5, 0.5].

      @param delta
          A pointer to 3 doubles with the (relative) vector. It will be modified in-place.

      After calling the wrap method, the fractional coordinates of delta will be in the
      range [-0.5, 0.5[.

      This is an approximate implementation of the minimum image convention that sometimes
      fails in very skewed cells, i.e. the wrapped vector is not always the shortest
      relative vector between a reference point and all of its periodic images. For more
      details see:
      http://scicomp.stackexchange.com/questions/3107/minimum-image-convention-for-triclinic-unit-cell
  */
  void iwrap(double* delta) const;


  /** @brief
          In-place addition of an integer linear combination of cell vectors to delta.

      @param delta
          A pointer to 3 doubles for the real-space vector to which the linear combination
          is added in-place.

      @param coeffs
          A pointer to 3 doubles with the coefficients of the linear combination.
   */
  void iadd_rvec(double* delta, const int* coeffs) const;


  /** @brief
          Get the ranges of cells within a cutoff radius.

      @param center
          A pointer to 3 doubles that specify the center of the cutoff sphere in
          real-space Cartesian coordinates.

      @param rcut
          The cutoff radius.

      @param ranges_begin
          A pointer to `nvec` ints, to which the begin of each range of cells along a cell
          vector is written. These integers are the highest indices of the crystal planes
          before/below the cutoff sphere.

      @param ranges_end
          A pointer to `nvec` ints to which the end of each range of cells along a cell
          vector is written. These integers are the lowest indices of the crystal planes
          after/above the cutoff sphere.

      This function effectively defines a supercell that is guaranteed to
      enclose the cutoff sphere.
   */
  int select_ranges_rcut(const double* center, const double rcut, int* ranges_begin,
      int* ranges_end) const;


  /** @brief
          Selects a cells inside or intersecting with a cutoff sphere. This function
          assumes space is divded in a regular grid of cells. The shape of one cell is
          by `rvecs` and `nvec`. This function finds all cells that contain a point within
          a cutoff sphere.

      @param center
          A pointer to 3 doubles that specify the center of the cutoff sphere in
          real-space Cartesian coordinates.

      @param rcut
          The cutoff radius.

      @param shape
          A pointer of 3 ints with the shape of the supercell.

      @param pbc
          A pointer to Boolean flags indicating the periodicity of the supercell along
          each cell vector.

      @param bars
          A std::vector<int> pointer in which the results, i.e. the cells overlapping with
          the cutoff sphere, are stored. In this output, a ranges of consecutive of cells
          along the last cell vector is called a bar and is represented by `(nvec + 1)`
          integers. Thus, the elements of the bars vectors should be used in groups of
          `(nvec + 1)`. For a single bar, the last two integers are the fractional
          coordinates of enclosing crystal planes along the last periodic vector. The
          preceding indexes in a bar are used to identify the position of the first cell
          in a bar along all but the last real-space vectors. The range in fractional
          coordinates along these all-but-last directions is `(i, i + 1)`. Finally, the
          union of all bars is a volume that completely contains the cutoff sphere but
          does not contain a single cell that does not overlap with the cutoff sphere.

      @return
          The number bars. The size of the bars vector is `nbar*(nvec+1)`.
    */
  size_t select_bars_rcut(const double* center, const double rcut, const int* shape,
      const bool* pbc, std::vector<int>* bars) const;

 protected:
  /** @brief
          Constructor that assumes the caller takes care of the consistency of all
          arguments. All arguments are copied.
   */
  Cell(const double* rvecs, const int nvec, const double* gvecs,
       const double volume, const double gvolume,
       const double* rlengths, const double* glenths,
       const double* rspacings, const double* gspacings);

  /** @brief
          TODO
   */
  void select_bars_low(SphereSlice* slice, const int* shape,
      const bool* pbc, std::vector<int>* prefix, std::vector<int>* bars)
      const;

 private:
  double rvecs_[9];       //!< real-space vectors,       one per row, row-major
  const int nvec_;        //!< number of defined cell vectors
  double gvecs_[9];       //!< reciprocal-space vectors, one per row, row-major
  double volume_;         //!< volume (or area or length) of the cell
  double gvolume_;        //!< volume of the reciprocal cell
  double rlengths_[3];    //!< real-space vector lengths
  double glengths_[3];    //!< reciprocal-space vector lengths
  double rspacings_[3];   //!< spacing between real-space crystal planes
  double gspacings_[3];   //!< spacing between reciprocal-space crystal planes
};

/**
    @brief
        A standardized modulo operation geared toward boundary conditions.

    @param i
        The numerator of the integer division.

    @param shape
        The denominator of the integer division.

    @param pbc
        Whether periodic boundary conditions apply.

    @return
        `i % shape` (guaranteed to be positive) if `pbc` is true. If `pbc` is false,
        `i` is returned if `i` lies in `[0,shape[`. If `pbc` is false and `i` falls out of
        that range, `-1` is returned.
 */
int smart_wrap(int i, const int shape, const bool pbc);


}  // namespace celllists


#endif  // CELLLISTS_CELL_H_


// vim: textwidth=90 et ts=2 sw=2
