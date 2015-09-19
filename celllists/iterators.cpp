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


#include "celllists/iterators.h"

#include <vector>

#include "celllists/decomposition.h"


namespace celllists {


BarIterator3D::BarIterator3D(const std::vector<int>& bars)
    : bars_(bars), ibar_(0), nbar_(bars.size()/4), icell_({0, 0, 0}), shape_{0, 0, 0},
      coeffs_{0, 0, 0} {
  // First check if there is a consistent number of elements in the bars vector.
  if (bars.size() % 4 != 0)
    throw std::domain_error("The size of the bars vector is not consisent with nvec.");
  // If there is at least on bar, set icell using the first bar.
  if (nbar_ > 0) {
    icell_[0] = bars_[0];
    icell_[1] = bars_[1];
    icell_[2] = bars_[2];
  }
};


BarIterator3D::BarIterator3D(const std::vector<int>& bars, const int* shape)
    : bars_(bars), ibar_(0), nbar_(bars.size()/4), icell_({0, 0, 0}),
      shape_{shape[0], shape[1], shape[2]}, coeffs_{0, 0, 0} {
  // First check if there is a consistent number of elements in the bars vector.
  if (bars.size() % 4 != 0)
    throw std::domain_error("The size of the bars vector is not consisent with nvec.");
  // If there is at least on bar, set icell using the first bar.
  if (nbar_ > 0) {
    icell_[0] = robust_wrap(bars_[0], shape_[0], &coeffs_[0]);
    icell_[1] = robust_wrap(bars_[1], shape_[1], &coeffs_[1]);
    icell_[2] = robust_wrap(bars_[2], shape_[2], &coeffs_[2]);
  }
}


BarIterator3D& BarIterator3D::operator++() {
  // Check if continuing makes any sense.
  if (ibar_ == nbar_)
    throw std::range_error("Cannot iterate past the last bar.");
  ++icell_[2];
  if ((shape_[2] != 0) && (icell_[2] >= shape_[2])) {
    icell_[2] = 0;
    ++coeffs_[2];
  }
  if (icell_[2] + shape_[2]*coeffs_[2] == bars_[4*ibar_ + 3]) {
    ++ibar_;
    icell_[0] = robust_wrap(bars_[4*ibar_], shape_[0], &coeffs_[0]);
    icell_[1] = robust_wrap(bars_[4*ibar_ + 1], shape_[1], &coeffs_[1]);
    icell_[2] = robust_wrap(bars_[4*ibar_ + 2], shape_[2], &coeffs_[2]);
  }
  return *this;
}


BarIterator3D BarIterator3D::operator++(int) {
  throw std::logic_error("Don't use the post-increment operator if BarIterator3D.");
  return *this;
}


}  // namespace celllists

// vim: textwidth=90 et ts=2 sw=2
