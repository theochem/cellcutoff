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


class Point {
 public:
  int index;
  std::array<double, 3> cart;
  std::array<int, 3> icell;

  Point(const int index, const double* _cart);
  Point(const int index, const double* _cart, const int* _icell);
  bool operator<(const Point& other) const;
};


// A typedef for cell_map objects
typedef std::map<std::array<int, 3>, std::array<int, 2>> CellMap;

//! Assigns all cell indexes
void assign_icell(const Cell &subcell, std::vector<Point>* points);
void assign_icell(const Cell &subcell, std::vector<Point>* points, const int* shape);

//! Create a mapping from cell indices to a list of points
CellMap* create_cell_map(const std::vector<Point> &points);


inline int robust_wrap(int num, const int denom, int* division) {
  if (denom == 0) {
    *division = 0;
    return num;
  } else {
    *division = num/denom;
    num %= denom;
    if (num < 0) --*division;
    return (num + denom) % denom;
  }
}


inline int robust_wrap(const int num, const int denom) {
  if (denom == 0) {
    return num;
  } else {
    return ((num % denom) + denom) % denom;
  }
}


}  // namespace celllists


#endif  // CELLLISTS_CELL_H_

// vim: textwidth=90 et ts=2 sw=2
