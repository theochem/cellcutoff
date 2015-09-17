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


#include "celllists/decomposition.h"

#include <cmath>
#include <map>
#include <vector>

#include "celllists/cell.h"
#include "celllists/vec3.h"


namespace celllists {


Point::Point(const int index, const double* _cart) : index(index) {
  std::copy(_cart, _cart + 3, cart.data());
  std::fill(icell.begin(), icell.end(), 0);
}


Point::Point(const int index, const double* _cart, const int* _icell) : index(index) {
  std::copy(_cart, _cart + 3, cart.data());
  std::copy(_icell, _icell + 3, icell.data());
}


bool Point::operator<(const Point& other) const {
  return icell < other.icell;
}


void assign_icell(const Cell &subcell, std::vector<Point>* points) {
  if (!(subcell.nvec() == 3))
    throw std::domain_error("Partitioning is only sensible for 3D subcells.");
  for (auto& point : *points) {
    double frac[3];
    subcell.to_frac(&point.cart[0], frac);
    point.icell[0] = static_cast<int>(floor(frac[0]));
    point.icell[1] = static_cast<int>(floor(frac[1]));
    point.icell[2] = static_cast<int>(floor(frac[2]));
  }
}

void assign_icell(const Cell &subcell, std::vector<Point>* points, const int* shape,
    const bool* pbc){
  if (!(subcell.nvec() == 3))
    throw std::domain_error("Partitioning is only sensible for 3D subcells.");
  for (auto& point : *points) {
    double frac[3];
    subcell.to_frac(&point.cart[0], frac);
    for (int ivec = 0; ivec < 3; ++ivec) {
      int i = static_cast<int>(floor(frac[ivec]));
      if (pbc[ivec])
        point.icell[ivec] = robust_wrap(i, shape[ivec]);
      i -= point.icell[ivec];
      vec3::iadd(point.cart.data(), subcell.vec(ivec), -i);
    }
  }
}


CellMap* create_cell_map(const std::vector<Point> &points) {
  auto result(new std::map<std::array<int, 3>, std::array<int, 2>>);
  const int* icell_begin = &(points[0].icell[0]);
  int ibegin = 0;
  int ipoint = 0;
  for (const auto& point : points) {
    if ((icell_begin[0] != point.icell[0]) ||
        (icell_begin[1] != point.icell[1]) ||
        (icell_begin[2] != point.icell[2])) {
      // Store range
      result->emplace(
        std::array<int, 3>{icell_begin[0], icell_begin[1], icell_begin[2]},
        std::array<int, 2>{ibegin, ipoint}
      );
      // New `begin`
      icell_begin = &(point.icell[0]);
      ibegin = ipoint;
    }
    ++ipoint;
  }
  // Storing the last range
  result->emplace(
    std::array<int, 3>{icell_begin[0], icell_begin[1], icell_begin[2]},
    std::array<int, 2>{ibegin, static_cast<int>(points.size())}
  );
  return result;
}


}  // namespace celllists

// vim: textwidth=90 et ts=2 sw=2
