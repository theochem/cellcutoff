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
// aint with this program; if not, see <http://www.gnu.org/licenses/>
//
//--


#include "celllists/decomposition.h"

#include <cmath>
#include <map>
#include <vector>

#include "celllists/cell.h"
#include "celllists/vec3.h"


namespace celllists {


Point::Point(const int index, const double* _cart) : index(index) {
  std::copy(_cart, _cart + 3, cart);
  std::fill(icell, icell + 3, 0.0);
}


Point::Point(const int index, const double* _cart, const int* _icell) : index(index) {
  std::copy(_cart, _cart + 3, cart);
  std::copy(_icell, _icell + 3, icell);
}


bool Point::operator<(const Point& other) const {
  if (icell[0] < other.icell[0]) return true;
  if (icell[0] > other.icell[0]) return false;
  if (icell[1] < other.icell[1]) return true;
  if (icell[1] > other.icell[1]) return false;
  if (icell[2] < other.icell[2]) return true;
  return false;
}


Cell* create_subcell(const Cell* cell, const int* shape, const double* spacings, bool* pbc) {\
  // Start by copying all cell vectors
  double new_vecs[9];
  std::copy(cell->vecs(), cell->vecs() + 9, new_vecs);
  // Divide the active dimensions according to shape
  for (int ivec = 0; ivec < cell->nvec(); ++ivec) {
    vec3::iscale(new_vecs + 3*ivec, 1.0/shape[ivec]);
    pbc[ivec] = true;
  }
  // Divide the remaining dimensions by given spacing
  for (int ivec = cell->nvec(); ivec < 3; ++ivec) {
    vec3::iscale(new_vecs + 3*ivec, spacings[ivec - cell->nvec()]);
    pbc[ivec] = false;
  }
  // Return the subcell
  return new Cell(new_vecs, 3);
}


void partition(std::vector<Point>* points, Cell* subcell) {
  if (!(subcell->nvec() == 3))
    throw std::domain_error("Partitioning is only sensible for 3D subcells.");
  for (auto& point : *points) {
    double frac[3];
    subcell->to_frac(point.cart, frac);
    point.icell[0] = static_cast<int>(floor(frac[0]));
    point.icell[1] = static_cast<int>(floor(frac[1]));
    point.icell[2] = static_cast<int>(floor(frac[2]));
  }
}


std::map<std::array<int, 3>, std::array<int, 2>>* create_map(const std::vector<Point>* points) {
  auto result(new std::map<std::array<int, 3>, std::array<int, 2>>);
  const int* icell_begin = points->at(0).icell;
  int ibegin = 0;
  int ipoint = 0;
  for (const auto& point : *points) {
    if ((icell_begin[0] != point.icell[0]) ||
        (icell_begin[1] != point.icell[1]) ||
        (icell_begin[2] != point.icell[2])) {
      // Store range
      result->emplace(
        std::array<int, 3>{icell_begin[0], icell_begin[1], icell_begin[2]},
        std::array<int, 2>{ibegin, ipoint}
      );
      // New `begin`
      icell_begin = point.icell;
      ibegin = ipoint;
    }
    ++ipoint;
  }
  // Storing the last range
  result->emplace(
    std::array<int, 3>{icell_begin[0], icell_begin[1], icell_begin[2]},
    std::array<int, 2>{ibegin, static_cast<int>(points->size())}
  );
  return result;
}


}  // namespace celllists

// vim: textwidth=90 et ts=2 sw=2
