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

    The main module of the cellcutoff library: a convenient implementation of a Cell
    object describibg a triclinic box, e.g. to define periodic boundary conditions.

  */


#ifndef CELLCUTOFF_CELL_H_
#define CELLCUTOFF_CELL_H_

#include <exception>
#include <memory>
#include <string>
#include <vector>


namespace cellcutoff {


/** @brief
        An exception for singular cell vectors.
 */
class singular_cell_vectors : public std::domain_error {
 public:
  //! Create exception, default constructor
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
    is specific for 3D systems. In case of 0D, 1D or 2D PBC, the "active" cell vectors are
    internally extended with "inactive" cell vectors, i.e. an orthonormal set of basis
    vectors to guarantee an invertible transformation between Cartesian and fractional
    coordinates. In that case, the fractional coordinates are actually also Cartesian
    coordinates in directions orthogonal to the available cell vectors. The extra basis
    vectors are always such that the complete set of vectors is right-handed.
 */
class Cell {
 public:
  /** @brief
          Create a Cell object.

      @param vecs
          A pointer to `3*nvec` doubles that represent the Cartesian cell vectors in
          row-major ordering. The vectors must not have a linear dependency. Each vector
          is one row in a 3x3 matrix.

      @param nvec
          The number of cell vectors. This corresponds to the dimensionality of the cell.
          `nvec` must be 0, 1, 2 or 3.
  */
  Cell(const double* vecs, const int nvec);

  //! Construct a Cell object without cell vectors.
  Cell() : Cell(nullptr, 0) {}

  //! Copy constructor
  Cell(const Cell& cell) : Cell(cell.vecs(), cell.nvec()) {}

  /** @brief
          Create a Cell object with the reciprocal cell.

      The caller owns the memory of the returned pointer. This is practically equivalent
      to `Cell(cell.gvecs(), cell.nvec())`, but more efficient and without precision loss.

      @return
          A pointer to a newly allocated reciprocal cell.
   */
  Cell* create_reciprocal() const;


  /** @brief
          Helper to construct a subcell of a given cell.

      This partitions space into bins (with the size and shape of the subcell) that can be
      used to do a domain decomposition.

      @param threshold
          The spacing between the subcell crystal planes is guaranteed to be lower than
          or equal to this value. Within the constraints of an integer division of the
          active cell vectors and this threshold, the spacing between the subcell crystal
          planes as as large as possible. The spacing along inactive cell vectors is equal
          to the threshold.

      @param shape
          A point to three ints. For each active cell vector i, shape[i] will contain the
          integer number of subcells along cell vector i. For inactive cell vectors i,
          shape[i] is set to zero.

      @return
          A pointer to a `Cell` object with the subcell.
   */
  Cell* create_subcell(const double threshold, int* shape) const;


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

      This member function is deprecated since version 0.3 and will be removed in version
      1.0. Please use the free function ranges_cutoff from the iterators module instead.

   */
  size_t ranges_cutoff(const double* center, const double cutoff, int* ranges_begin,
      int* ranges_end) const;


  /** @brief
          Selects cells inside or at least partially overlapping with a cutoff sphere.

      This member function is deprecated since version 0.3 and will be removed in version
      1.0. Please use the free function ranges_cutoff from the iterators module instead.

    */
  void bars_cutoff(const double* center, const double cutoff, std::vector<int>* bars) const;

 protected:
  /** @brief
          Constructor that assumes the caller takes care of the consistency of all
          arguments. All arguments are copied to data members.
   */
  Cell(const double* vecs, const int nvec, const double* gvecs,
       const double volume, const double gvolume,
       const double* lengths, const double* glenths,
       const double* spacings, const double* gspacings);

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


// Free functions


/** @brief
        Return a randomized but sensible unit cell.

    @param seed
        The random seed used to initialize the cell vectors.

    @param nvec
        The dimensionality of the cell.

    @param scale
        Determines the overall size of the random cell. The components of the cell vectors
        are sampled from a uniform random distribution over the interval [-scale, scale].

    @param ratio
        Controls the minimal volume of the cell. Random cell vectors are tried until the
        volume is above (ratio*scale)**nvec. When a suitable volume is found, these
        cell vectors are used. The closer to 2, the more cubic the cell. High values
        of this parameter may make this function slow because many trials will be needed.

    @param cuboid
        When true, an orthorhombic cell is constructed.
  */
Cell* create_random_cell(const unsigned int seed, const int nvec, const double scale = 1.0,
  const double ratio = 0.1, const bool cuboid = false);


}  // namespace cellcutoff


#endif  // CELLCUTOFF_CELL_H_

// vim: textwidth=90 et ts=2 sw=2
