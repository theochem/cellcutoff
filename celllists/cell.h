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
    cell. Reciprocal cell vectors, cell vector lengths and spacings between crystal planes
    are computed immediately. Several manipulations of fractional or Cartesian coordinates
    are implemented.

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

      @param vecs
          A pointer to `3*nvec` doubles that represent the Cartesian cell vectors in
          row-major ordering. The vectors must not have a linear dependency. Each vector
          is one row in a 3x3 matrix.

      @param nvec
          The number of cell vectors. This corresponds to the dimensionality of the cell.
          `nvec` must be 0, 1, 2 or 3.
  */
  Cell(const double* vecs, const int nvec);

  // Copy-constructor, move-constructor and assignment make no sense as the Cell is
  // constant after construction! Just pass a reference or a pointer instead.
  Cell(const Cell& that) = delete;
  Cell(Cell&&) = delete;
  Cell& operator=(const Cell&) = delete;

  /** @brief
          Create a Cell object with the reciprocal cell.

      The caller owns the memory of the returned pointer. This is practically equivalent
      to `Cell(cell.gvecs(), cell.nvec())`, but more efficient and without precision loss.

      @return
          A pointer to a newly allocated reciprocal cell.
   */
  Cell* create_reciprocal() const;

  //! Returns the number of periodic dimensions.
  int nvec() const { return nvec_; }

  //! Returns all cell vectors.
  const double* vecs() const { return vecs_; }

  //! Returns a cell vector.
  const double* vec(const int ivec) const;

  //! Returns all reciprocal cell vectors.
  const double* gvecs() const { return gvecs_; }

  //! Returns a reciprocal cell vector.
  const double* gvec(const int ivec) const;

  //! Returns the volume (or area or length) of the cell.
  double volume() const { return volume_; }

  //! Returns the volume (or area or length) of the reciprocal cell.
  double gvolume() const { return gvolume_; }

  //! Returns the lengths of the cell vectors.
  const double* lengths() const { return lengths_; }

  //! Returns the lengths of the reciprocal cell vectors.
  const double* glengths() const { return glengths_; }

  //! Returns the spacings between the crystal plane
  const double* spacings() const { return spacings_; }

  //! Returns the spacings between the reciprocal crystal plane
  const double* gspacings() const { return gspacings_; }


  /** @brief
          Returns `true` only when cell is cubic and aligned with Cartesian
          axes: a to x, b to y and c to z. No small errors allowed.
    */
  bool cubic() const;


  /** @brief
          Returns `true` only when cell is cuboid (orthorombic) and aligned with Cartesian
          axes: a to x, b to y and c to z. No small errors allowed.
    */
  bool cuboid() const;


  /** @brief
          Convert Cartesian coordinates to fractional.

      This effectively computes the dot products of the Cartesian vector with the
      reciprocal cell vectors.

      @param cart
          A pointer to 3 doubles containing the input cell vector.

      @param frac
          A pointer to 3 doubles in which the output is written.
   */
  void to_frac(const double* cart, double* frac) const;


  /** @brief
          Convert fractional coordinates to Cartesian.

      This effectively computes a linear combination of cell vectors.

      @param frac
          A pointer to 3 doubles containing the input fractional coordinates.

      @param cart
          A pointer to 3 doubles to which the output is written
   */
  void to_cart(const double* frac, double* cart) const;


  /** @brief
          In-place wrap a (relative) Cartesian vector back into the cell [-0.5, 0.5[.

      After calling the `iwrap_mic` method, the fractional coordinates of delta will be in
      the range [-0.5, 0.5[.

      This is an approximate implementation of the minimum image convention that sometimes
      fails in very skewed cells, i.e. the wrapped vector is not always the shortest
      relative vector between a reference point and all of its periodic images. For more
      details see:
      http://scicomp.stackexchange.com/questions/3107/minimum-image-convention-for-triclinic-unit-cell

      @param delta
          A pointer to 3 doubles with the (relative) Cartesian vector. It will be modified
          in-place.
  */
  void iwrap_mic(double* delta) const;


  /** @brief
          In-place wrap a (relative) Cartesian vector back into the cell [0.0, 1.0[.

      After calling the `iwrap_box` method, the fractional coordinates of delta will be in
      the range [0.0, 1.0[.

      @param delta
          A pointer to 3 doubles with the (relative) Cartesian vector. It will be modified
          in-place.
  */
  void iwrap_box(double* delta) const;


  /** @brief
          In-place addition of an integer linear combination of cell vectors to delta.

      @param delta
          A pointer to 3 doubles for the Cartesian vector to which the linear combination
          is added in-place.

      @param coeffs
          A pointer to 3 ints with the coefficients of the linear combination.
   */
  void iadd_vec(double* delta, const int* coeffs) const;


  /** @brief
          Get the ranges of cells within a cutoff radius.

      This function effectively defines a supercell that is guaranteed to enclose the
      cutoff sphere.

      @param center
          A pointer to 3 doubles that specify the center of the cutoff sphere in
          Cartesian coordinates.

      @param cutoff
          The cutoff radius.

      @param ranges_begin
          A pointer to `nvec` ints, to which the begin of each range of cells along a cell
          vector is written. These integers are the highest indices of the crystal planes
          below the cutoff sphere.

      @param ranges_end
          A pointer to `nvec` ints to which the end of each range of cells along a cell
          vector is written. These integers are the lowest indices of the crystal planes
          above the cutoff sphere.

      @return
          The number of cells contained in the supercell.
   */
  int ranges_cutoff(const double* center, const double cutoff, int* ranges_begin,
      int* ranges_end) const;


  /** @brief
          Selects a cells inside or intersecting with a cutoff sphere.

      This function assumes space is divded in a regular grid of cells. The shape of one
      cell is by `vecs` and `nvec`. This function finds all cells that contain a point
      within a cutoff sphere.

      @param center
          A pointer to 3 doubles that specify the center of the cutoff sphere in
          Cartesian coordinates.

      @param cutoff
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
          in a bar along all but the last cell vectors. The range in fractional
          coordinates along these all-but-last directions is `(i, i + 1)`. Finally, the
          union of all bars is a volume that completely contains the cutoff sphere but
          does not contain a single cell that does not overlap with the cutoff sphere.

      @return
          The number bars. The size of the bars vector is `nbar*(nvec+1)`.
    */
  size_t bars_cutoff(const double* center, const double cutoff, const int* shape,
      const bool* pbc, std::vector<int>* bars) const;

  /** @brief
          Helper to construct a subcell of a given cell.

          This partitions space into bins (with the size and shape of the subcell) that
          can be used to do a domain decomposition.

      @param shape
          A pointer to nvec integers, with the number subcells along the corresponding
          cell vector.

      @param spacings
          A pointer to (3-nvec) doubles, with the spacings between the
          inactive/non-periodic cell vectors.

      @param pbc
          A pointer to 3 bools. This is an auxiliary output argument, whose first nvec
          elements are set to true, while the remaining are set to false.

      @return
          A pointer to a `Cell` object with the subcell.
   */
  Cell* create_subcell(const int* shape, const double* spacings, bool* pbc);

 protected:
  /** @brief
          Constructor that assumes the caller takes care of the consistency of all
          arguments. All arguments are copied to data members.
   */
  Cell(const double* vecs, const int nvec, const double* gvecs,
       const double volume, const double gvolume,
       const double* lengths, const double* glenths,
       const double* spacings, const double* gspacings);

  /** @brief
          Low-level functions used by bars_cutoff.

      TODO. This method may change in future, so I'm not going to try explain it in
      detail. This can be fixed after the `sphere_slice` will be completely finalized.
      (The current implementation sphere_slice and bars_cutoff is good but not optimal.)

      This method goes recursively through all active cell vectors and divides space along
      this axis in cells that overlap with the cutoff sphere/circle/line, depending on
      the dimension at hand (i.e. the recursion depth). It makes use of the SphereSlice\
      object to find the begin-end range along each cell vector.
   */
  void bars_cutoff_low(SphereSlice* slice, const int* shape,
      const bool* pbc, std::vector<int>* prefix, std::vector<int>* bars)
      const;

 private:
  double vecs_[9];        //!< cell vectors, one per row, row-major
  const int nvec_;        //!< number of defined cell vectors
  double gvecs_[9];       //!< reciprocal cell vectors, one per row, row-major
  double volume_;         //!< volume (or area or length) of the cell
  double gvolume_;        //!< volume of the reciprocal cell
  double lengths_[3];     //!< cell vector lengths
  double glengths_[3];    //!< reciprocal cell vector lengths
  double spacings_[3];    //!< spacing between crystal planes
  double gspacings_[3];   //!< spacing between reciprocal crystal planes
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
