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

#include <array>
#include <unordered_map>
#include <vector>

#include "cellcutoff/cell.h"
#include "cellcutoff/sphere_slice.h"
#include "cellcutoff/decomposition.h"


namespace cellcutoff {


/** @brief
        Get the ranges of cells within a cutoff radius.

    This function assumes the space is divided into boxes by crystal planes at integer
    indexes. For example, these planes for the first cell vector are parallel to the
    second and third cell vector and have one point in (first vector)*index where index
    is the integer index for these planes. Similar definitions are used for the other
    two cell vectors. The returned ranges are arrays referring to the integer indexes
    that demarcate the cutoff sphere. One could interpret the result as a supercell
    that contains the entire cutoff sphere.

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
size_t cutoff_ranges(const Cell* cell, const double* center, const double cutoff,
    int* ranges_begin, int* ranges_end);


/** @brief
        Selects cells inside or at least partially overlapping with a cutoff sphere.

    This function assumes space is divided in a regular grid of subcells. The shape of
    one subcell is defined by `vecs` and `nvec` (>= 1). This function then finds all
    subcells that overlap with a cutoff sphere.

    @param center
        A pointer to 3 doubles that specify the center of the cutoff sphere in
        Cartesian coordinates.

    @param cutoff
        The cutoff radius.

    @param bars
        A std::vector<int> pointer in which the results, i.e. the cells overlapping with
        the cutoff sphere, are stored.

        To keep the array compact, the following format is used to specify all cells
        that overlap with the cutoff sphere. The integers are always to be interpreted
        in (begin, end) pairs, corresponding to crystal planes just before and after the
        cutoff sphere. The first pair, (begin0, end0), corresponds to the planes along
        the [100] direction. If `nvec==2`, a list of pairs follows, corresponding to
        (begin1, end1) ranges along the [010] direction. One such pair is present for
        each slice of the cutoff sphere along the [100] direction, i.e. for (begin0,
        begin0 + 1), (begin0 + 1, begin0 + 2), etc. Similarly, if `nvec==3`, each pair
        for the [010] direction is followed with a set of pairs for the [001] direction.

        The above format assumes that one know `nvec` when parsing the list of integers.

    @return
        The number bars. The size of the bars vector is `nbar*(nvec+1)`.
 */
void cutoff_bars(const Cell* cell, const double* center, const double cutoff,
    std::vector<int>* bars);


/** @brief
        Low-level functions used by cutoff_bars.

    TODO. This method may change in future, so I'm not going to try explaining it in
    detail. This can be fixed after the `sphere_slice` will be completely finalized.
    (The current implementation sphere_slice and cutoff_bars is good but not optimal.)

    This method goes recursively through all active cell vectors and divides space along
    this axis in cells that overlap with the cutoff sphere/circle/line, depending on
    the dimension at hand (i.e. the recursion depth). It makes use of the SphereSlice\
    object to find the begin-end range along each cell vector.
 */
void cutoff_bars_low(const Cell* cell, SphereSlice* slice, int ivec, std::vector<int>* bars);



/** @brief
        Loop over cells that have some overlap with a cutoff sphere.

    The cells to loop over must be provided as "bars", which can be generated with the
    function cutoff_bars. This class is mainly intended for usage by DeltaIterator,
    but in rare cases it might also be directly useful.
 */
class BarIterator {
 public:
  /** @brief
          Create a BarIterator.

      @param bars
          The output of the function cutoff_bars applied to a subcell.

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

  const std::vector<int>& bars_;  //!< Vector produced by cutoff_bars
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

    The usage of this class is deprecated since version 0.3. It will be removed in version
    0.1. Consider using BoxSortedPoints and BoxCutoffIterator instead. These should be
    easier to work with.
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
  std::vector<int> bars_;      //!< A bars vector created with cutoff_bars for subcell
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


//! Map type to locate sorted points easily in BoxSortedPoints.
typedef std::unordered_map<size_t, std::array<size_t, 2>> RangesMap;


/** @brief
        Return a positive serial for given (sub)cell indices.

    @param icell
        Array with cell indexes.

    @return serial
        The positive integer identifying the cell.
 */
size_t serialize_icell(const int* icell);


/** @brief
        Return a positive serial for given (sub)cell indices.

    @param i0
        cell index 0.

    @param i1
        cell index 1.

    @param i2
        cell index 2.

    @return serial
        The positive integer identifying the cell.
 */
size_t serialize_icell(const int i0, const int i1, const int i2);


/** @brief
        Estimate a suitable threshold for binning the points into subcells.

    @param points
        C-contigiuous array of 3D Cartesian coordinates.

    @param npoint
        Number of points.

    @param cell
        Describing the periodic boundary conditions, if any.

    @return threshold
        The maximum distance between subsequent subcell planes.
 */
double sensible_threshold(const double* points, size_t npoint, const Cell* cell);


/** @brief
        Sort points and construct iterators over points within cutoff. Added in version 0.3.
  */
class BoxSortedPoints {
 public:
  /** @brief
          Create a BoxSortedPoints object.

      This would group all the points into subcells of the given cell. The size of the
      subcells is controlled by the threshold parameter.

      @param points
          C-contigiuous array of 3D Cartesian coordinates.

      @param npoint
          Number of points.

      @param cell
          Describing the periodic boundary conditions, if any.

      @param threshold
          The maximum spacing between opposite faces of the subcell. When not given or
          when not positive, a sensible default threshold is determined, which is the
          recommended usage.

   */
  BoxSortedPoints(const double* points, size_t npoint, const Cell* cell,
    double threshold = 0.0);
  //! Destruct a BoxSortedPoints
  ~BoxSortedPoints();

  //! Points after sorting and other sorts of manipulations.
  const double* points() const { return points_; }
  //! Number of points.
  size_t npoint() const { return npoint_; }

  //! Subcell used for binning.
  const Cell* subcell() const { return subcell_; }
  //! The number of repetitions of the subcell to obtain the periodic cell.
  const int* shape() const { return shape_; }
  //! The index of each point in the original array.
  const size_t* ipoints() const { return ipoints_; }
  //! The mapping from serials to ranges.
  const RangesMap* ranges() const { return &ranges_; }

 private:
  double* points_;    //!< Points after wrapping in periodic cell.
  size_t npoint_;     //!< Number of points.
  Cell* subcell_;     //!< Subcell used for binning the points.
  int shape_[3];      //!< Number of repetitions along each subcell vector.
  /** @brief The reordering of the points.

      First point in a cell is points_[ipoints_[ranges->at(serial)[0]]].
   */
  size_t* ipoints_;
  //! The begin and end indices (after sorting) of the points within each subcell.
  RangesMap ranges_;
};


/** @brief
        Iterator over points within a cutoff sphere. Added in version 0.3.
  */
class BoxCutoffIterator {
 public:
  /** @brief
          Create an iterator over points within a cutoff sphere.

      @param bsp
          A BoxSortedPoints instance.

      @param center
          The center of the cutoff sphere.

      @params radius
          Radius of the cutoff sphere.

      @return ci
          An iterator over the points within the cutoff.
   */
  BoxCutoffIterator(const BoxSortedPoints* bsp, const double* center, double radius);

  //! Destruct the iterator.
  ~BoxCutoffIterator();

  //! Distance between current point and center.
  double distance() const { return distance_; }
  //! Relative vector from center to current point.
  const double* delta() const { return delta_; }
  //! Index of the current point in the user-provided points array.
  size_t ipoint() const { return ipoint_; }

  //! Return true when not at the end yet.
  bool busy() const { return bar_iterator_->busy(); }
  //! Iterate by one point.
  BoxCutoffIterator& operator++();
  //! Disables post-increment.
  BoxCutoffIterator operator++(int);

 private:
  /** @brief
          Low-level increment implementation.

      The parameter initialization is only set to True for the first increment call in
      the constructor.
    */
  void increment(bool initialization);

  const BoxSortedPoints* bsp_;
  std::vector<int> bars_;      //!< A bars vector created with cutoff_bars for subcell
  BarIterator* bar_iterator_;  //!< Bar iterator to loop over subcells
  double center_[3];
  double radius_;
  size_t ibegin_;
  size_t iend_;
  size_t icurrent_;
  size_t ipoint_;              //!< Index of current point
  //! Relative vector from the center to lower corner of the periodic cell
  double cell_delta_[3];
  double delta_[3];            //!< Relative vector from the center to the current point
  double distance_;            //!< Distance between center and current point
};


}  // namespace cellcutoff


#endif  // CELLCUTOFF_ITERATORS_H_

// vim: textwidth=90 et ts=2 sw=2
