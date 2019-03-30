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
//
// --

/** @file

    Provides tools to sort grid points into cells.

  */


#ifndef CELLCUTOFF_DECOMPOSITION_H_
#define CELLCUTOFF_DECOMPOSITION_H_

#include <array>
#include <vector>
#include <unordered_map>
#include <string>

#include "cellcutoff/cell.h"


namespace cellcutoff {


/** @brief
        An exception for Point arrays that are not properly grouped by icell.
 */
class points_not_grouped : public std::domain_error {
 public:
  //! Create exception, default constructor
  explicit points_not_grouped(const std::string& what_arg)
      : std::domain_error(what_arg) {}
};


/** @brief
        A Point: its Cartesian coordinates and the indices of the cell in which it sits.

    A comparator is implemented to facilitate sorting of points.
 */
class Point {
 public:
  //! Initialize with just the Cartesian coordinates of the point alone.
  explicit Point(const double* cart);
  //! Initialize with Cartesian coordinates and cell indices.
  Point(const double* cart, const int* icell);
  //! Compare two points, only cell indices matter.
  bool operator<(const Point& other) const;

  double cart_[3];  //!< Cartesian coordinates of the point
  int icell_[3];    //!< Cell indices of the containing grid cell
};


/** @brief
        A hash for using an integer array of three elements as keys in a std::unordered_map.
 */
struct icell_hash {
  /** @brief
          Compute a unique (positive) integer from a vector of three integer.

      This is the usual folding trick, such that cells closer to the origin are assigned
      lower integer indices.
    */
  size_t operator()(const std::array<int, 3>& icell) const {
    const int small = 4*(icell[0] < 0) + 2*(icell[1] < 0) + (icell[2] < 0);
    const size_t x = abs(icell[0]) + (icell[0] >= 0);
    const size_t y = abs(icell[1]) + (icell[1] >= 0);
    const size_t z = abs(icell[2]) + (icell[2] >= 0);
    const size_t d0 = x + y + z;
    const size_t d1 = x + y;
    return (((d0-3)*(d0-2)*(d0-1))/6 + ((d1-2)*(d1-1))/2 + (x-1))*8 + small;
  }
};

//! Unordered map for sorting points into cells.
typedef std::unordered_map<std::array<int, 3>, std::array<size_t, 2>, icell_hash> CellMap;

/** @brief Assigns all cell indexes to an array of Points, not wrapping the cell indices.
 */
void assign_icell(const Cell &subcell, void* points, size_t npoint, size_t point_size);

//! Assigns all cell indexes to an array of points, wrapp
void assign_icell(const Cell &subcell, const int* shape, void* points, size_t npoint,
    size_t point_size);

//! Sort function for Point array
void sort_by_icell(void* points, size_t npoint, size_t point_size);

//! Create a mapping from cell indices to a list of points
CellMap* create_cell_map(const void* points, size_t npoint, size_t point_size);

//! Safe modulus operation with compatible division
inline int robust_wrap(int index, const int size, int* division) {
  if (size == 0) {
    *division = 0;
    return index;
  } else {
    *division = index/size;
    index %= size;
    if (index < 0) --*division;
    return (index + size) % size;
  }
}

//! Safe modulus operation
inline int robust_wrap(const int index, const int size) {
  if (size == 0) {
    return index;
  } else {
    return ((index % size) + size) % size;
  }
}


}  // namespace cellcutoff


#endif  // CELLCUTOFF_DECOMPOSITION_H_

// vim: textwidth=90 et ts=2 sw=2
