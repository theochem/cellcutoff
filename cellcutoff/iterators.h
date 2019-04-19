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

    Two iterators for efficient real-space cutoff implementations. Useful for when one
    is interested in looping over many points without a cutoff radius around any given
    center.

    Terminology used in the doc strings:

    - A periodic cell is a Cell object describing a triclinic unit cell.
    - A subcell is a Cell object describing a portion of the unit cell. The shape array
      shows how many times a subcell must be repeated in each direction to reconstruct
      the whole periodic cell.
  */


#ifndef CELLCUTOFF_ITERATORS_H_
#define CELLCUTOFF_ITERATORS_H_

#include <vector>
#include <array>

#include "cellcutoff/cell.h"
#include "cellcutoff/decomposition.h"


namespace cellcutoff {


/** @brief
        Loop over cells that have some overlap with a cutoff sphere.

    The cells to loop over must be provided as "bars", which can be generated with the
    function Cell::bars_cutoff. This class is mainly intended for usage by DeltaIterator,
    but in rare cases it might also be directly useful.
  */
class BarIterator {
 public:
  /** @brief
          Create a BarIterator.

      @param bars
          The output of the function Cell::bars_cutoff applied to a subcell.

      @param nvec
          The number of periodic dimensions. This is needed to interpret the bars argument
          correctly.

      @param shape
          Optional argument. When given and strictly positive, the iterator only loops
          over cells with indices 0 <= icell[i] < shape[i]. Cells outside this range
          are wrapped into this range, i.e. they are not ignored.
    */
  BarIterator(const std::vector<int>& bars, const int nvec, const int* shape);
  //! Initialize the BarIterator without shape argument.
  BarIterator(const std::vector<int>& bars, const int nvec)
      : BarIterator(bars, nvec, nullptr) {}
  //! Destruct the iterator.
  ~BarIterator();

  //! Return true if the iterator has not reached the last cell yet.
  bool busy() const { return busy_; }

  //! Iterate forward by one step.
  BarIterator& operator++();
  //! Disables post-increment.
  BarIterator operator++(int);

  //! The cell indices of the current iteration, wrapped.
  const int* icell() const { return icell_; }
  //! The cell indices of the periodic image in which the current (wrapped) cell sits.
  const int* coeffs() const { return coeffs_; }

 private:
  //! Consume one range of cells from the bars vector. increments ibar_ appropriately.
  void take_range(const int ivec);
  /** @brief
          Iterate to the next cell along a given dimension.

      When it reaches the end at a given dimension, it recurses and increments at a
      preceding dimensions.
    */
  void increment(const int ivec);

  const std::vector<int>& bars_;  //!< Vector produced by Cell::bars_cutoff
  const int nvec_;                //!< Number of dimensions encoded in the bars vector
  size_t ibar_;                   //!< Current position in the bars vector
  int* shape_;                    //!< Number of subcells in the periodic cell along each vector
  int* ranges_begin_;             //!< Begin of cell indec ranges to iterate over
  int* ranges_end_;               //!< Non-inclusive end of cell indec ranges to iterate over
  int* icell_unwrapped_;          //!< Cell indices at the current iterator before wrapping
  int* icell_;                    //!< Cell indices at the current iterator
  int* coeffs_;                   //!< Indices of the periodic image of the current cell
  bool busy_;                     //!< true as long as the iterator has not reached end
};


/** @brief
        Loop over points within a cutoff sphere around some center.
  */
class DeltaIterator {
 public:
  /** @brief
          Create a DeltaIterator.

      @param subcell
          A Cell object describing the volume of one grid cell in which points are sorted.

      @param shape
          The number of times the subcell has to be repeated in each direction to obtain
          the orthorhombic periodic unit cell.

      @param center
          The center of the cutoff sphere (size=3).

      @param cutoff
          The radius of the cutoff sphere.

      @param points
          An array of Point objects. These must be sorted by the caller with the function
          sort_by_icell.

      @param npoint
          The number of points.

      @param point_size
          The size in bytes of a Point object.

      @param cell_map
          An CellMap dictionary with icell vectors as keys and (start, size) pairs as
          values. This is the output of the create_cell_map function.
    */
  DeltaIterator(const Cell& subcell, const int* shape, const double* center,
      const double cutoff, const void* points, const size_t npoint,
      const size_t point_size, const CellMap& cell_map);
  //! See other DeltaIterator constructor. This one is useful for aperiodic systems.
  DeltaIterator(const Cell& subcell, const double* center,
      const double cutoff, const void* points, const size_t npoint,
      const size_t point_size, const CellMap& cell_map)
      : DeltaIterator(subcell, nullptr, center, cutoff, points, npoint, point_size,
        cell_map) {}
  //! Destruct a DeltaIterator
  ~DeltaIterator();

  //! Return true when not at the end yet.
  bool busy() const { return bar_iterator_->busy(); }
  //! Iterate by one point.
  DeltaIterator& operator++();
  //! Disables post-increment.
  DeltaIterator operator++(int);

  //! The relative vector from center to point, at the current iteration.
  const double* delta() const { return delta_; }
  //! The distance between center and point, at the current iteration.
  double distance() const { return distance_; }
  //! The index of the point in the points array.
  size_t ipoint() const { return ipoint_; }

 private:
  /** @brief
          Low-level increment implementation.

      The parameter initialization is only set to True for the first increment call. Ugly
      trick.
    */
  void increment(bool initialization);

  // Provided through constructor
  const Cell& subcell_;        //!< Description of volume in which points are grouped
  int* shape_;                 //!< Number of repetitions of subcell to obtain periodic cell
  const double center_[3];     //!< Center of the cutoff sphere
  const double cutoff_;        //!< Radius of the cutoff sphere
  const char* points_char_;    //!< Opaque pointer to array of sorted points
  const size_t npoint_;        //!< Number of points
  const size_t point_size_;    //!< Size of one point object
  const CellMap& cell_map_;    //!< Segments in sorted points for each cell index

  // Internal data
  std::vector<int> bars_;      //!< A bars vector created with subcell::bars_cutoff
  BarIterator* bar_iterator_;  //!< Bar iterator to loop over subcells
  const Point* point_;         //!< Current point object
  /** @brief
          Relative vector from the center to lower corner of the periodic cell
    */
  double cell_delta_[3];
  double delta_[3];            //!< Relative vector from the center to the current point
  double distance_;            //!< Distance between center and current point
  size_t ipoint_;              //!< Index of current point
  size_t ibegin_;              //!< First index of points in current subcell
  size_t iend_;                //!< Non-inclusive last index of points in current subcell
};


}  // namespace cellcutoff


#endif  // CELLCUTOFF_ITERATORS_H_

// vim: textwidth=90 et ts=2 sw=2
