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


#include "cellcutoff/decomposition.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "cellcutoff/cell.h"
#include "cellcutoff/iterators.h"
#include "cellcutoff/vec3.h"


namespace cellcutoff {


Point::Point(const double* cart) {
  std::copy(cart, cart + 3, cart_);
  std::fill(icell_, icell_ + 3, 0);
}


Point::Point(const double* cart, const int* icell) {
  std::copy(cart, cart + 3, cart_);
  std::copy(icell, icell + 3, icell_);
}


bool Point::operator<(const Point& other) const {
  // Lexicographical less than (equivalent to std::lexicographical_compare but faster)
  if (icell_[0] < other.icell_[0]) return true;
  if (icell_[0] > other.icell_[0]) return false;
  if (icell_[1] < other.icell_[1]) return true;
  if (icell_[1] > other.icell_[1]) return false;
  if (icell_[2] < other.icell_[2]) return true;
  return false;
}


size_t icell_hash::operator()(const std::array<int, 3>& icell) const {
  return serialize_icell(icell.data());
}


void assign_icell(const Cell &subcell, void* points, size_t npoint, size_t point_size) {
  // Check args
  if (!(subcell.nvec() == 3))
    throw std::domain_error("Partitioning is only sensible for 3D subcells.");
  // Loop over all points, compute icell
  char* points_char = reinterpret_cast<char*>(points);  // Ugly sweet hack
  for (size_t ipoint = 0; ipoint < npoint; ++ipoint) {
    Point* point(reinterpret_cast<Point*>(points_char));  // Ugly sweet hack
    double frac[3];
    subcell.to_frac(point->cart_, frac);
    point->icell_[0] = static_cast<int>(floor(frac[0]));
    point->icell_[1] = static_cast<int>(floor(frac[1]));
    point->icell_[2] = static_cast<int>(floor(frac[2]));
    points_char += point_size;
  }
}


void assign_icell(const Cell &subcell, const int* shape, void* points, size_t npoint,
    size_t point_size) {
  // Check args
  if (!(subcell.nvec() == 3))
    throw std::domain_error("Partitioning is only sensible for 3D subcells.");
  // Loop over all points, compute icell and wrap if needed
  char* points_char = reinterpret_cast<char*>(points);  // Ugly sweet hack
  for (size_t ipoint = 0; ipoint < npoint; ++ipoint) {
    Point* point(reinterpret_cast<Point*>(points_char));  // Ugly sweet hack
    double frac[3];
    subcell.to_frac(point->cart_, frac);
    for (int ivec = 0; ivec < 3; ++ivec) {
      // Compute floored fractional coordinate, optionally wrapped.
      int i = static_cast<int>(floor(frac[ivec]));
      point->icell_[ivec] = robust_wrap(i, shape[ivec]);
      // Wrap point into box
      vec3::iadd(point->cart_, subcell.vec(ivec), point->icell_[ivec]-i);
    }
    points_char += point_size;
  }
}


static int cmp_points(const void *a, const void* b) {
  const Point* pa(reinterpret_cast<const Point*>(a));  // Ugly sweet hack
  const Point* pb(reinterpret_cast<const Point*>(b));  // Ugly sweet hack
  // Lexicographical comparison
  int result = pa->icell_[0] - pb->icell_[0];
  if (result != 0) return result;
  result = pa->icell_[1] - pb->icell_[1];
  if (result != 0) return result;
  return pa->icell_[2] - pb->icell_[2];
}


void sort_by_icell(void* points, size_t npoint, size_t point_size) {
  qsort(points, npoint, point_size, cmp_points);
}


static inline void _store_in_cell_map(const int* icell, size_t ibegin, size_t iend,
    CellMap* cell_map) {
  // Try to store the new range in the cell_map
  auto emplace_output = cell_map->emplace(
    std::array<int, 3>{{icell[0], icell[1], icell[2]}},
    std::array<size_t, 2>{{ibegin, iend}});
  // If the is already present, the input for create_cell_map was incorrect.
  if (!emplace_output.second) {
    delete cell_map;
    throw points_not_grouped("The given points are not grouped by icell.");
  }
}


CellMap* create_cell_map(const void* points, size_t npoint, size_t point_size) {
  auto cell_map(new CellMap);
  const char* points_char = reinterpret_cast<const char*>(points);  // Ugly sweet hack
  const Point* point(reinterpret_cast<const Point*>(points_char));  // Ugly sweet hack
  const int* icell_begin = point->icell_;
  size_t ibegin = 0;
  for (size_t ipoint = 0; ipoint < npoint; ++ipoint) {
    point = reinterpret_cast<const Point*>(points_char);  // Ugly sweet hack
    if ((icell_begin[0] != point->icell_[0]) ||
        (icell_begin[1] != point->icell_[1]) ||
        (icell_begin[2] != point->icell_[2])) {
      // Store
      _store_in_cell_map(icell_begin, ibegin, ipoint, cell_map);
      // New `begin` of a range of points
      icell_begin = point->icell_;
      ibegin = ipoint;
    }
    points_char += point_size;
  }
  // Storing the last range
  _store_in_cell_map(icell_begin, ibegin, npoint, cell_map);
  return cell_map;
}


}  // namespace cellcutoff

// vim: textwidth=90 et ts=2 sw=2
