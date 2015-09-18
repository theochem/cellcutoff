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


#ifndef CELLLISTS_DECOMPOSITION_H_
#define CELLLISTS_DECOMPOSITION_H_

#include <vector>
#include <map>

#include "celllists/cell.h"


namespace celllists {


/** @brief
        An exception for Point arrays that are not properly grouped by icell.
 */
class points_not_grouped : public std::domain_error {
 public:
  explicit points_not_grouped(const std::string& what_arg)
      : std::domain_error(what_arg) {}
};


class Point {
 public:
  Point(const double* cart);
  Point(const double* cart, const int* icell);
  bool operator<(const Point& other) const;

  double cart_[3];
  int icell_[3];
};


// A typedef for cell_map objects
typedef std::map<std::array<int, 3>, std::array<size_t, 2>> CellMap;

//! Assigns all cell indexes
void assign_icell(const Cell &subcell, void* points, size_t npoint, size_t point_size);
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


}  // namespace celllists


#endif  // CELLLISTS_CELL_H_

// vim: textwidth=90 et ts=2 sw=2
