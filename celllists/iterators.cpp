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

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "celllists/decomposition.h"


namespace celllists {


BarIterator::BarIterator(const std::vector<int>& bars, const int nvec, const int* shape)
    : bars_(bars), nvec_(nvec), ibar_(0), shape_(nullptr), icell_unwrapped_(nullptr),
      icell_(nullptr), coeffs_(nullptr), busy_(true) {
  // Argument checking
  if ((nvec < 1) || (nvec > 3))
    throw std::domain_error("BarIterator requires nvec to be 1, 2 or 3.");
  // Allocate arrays
  shape_ = new int[nvec];
  ranges_begin_ = new int[nvec];
  ranges_end_ = new int[nvec];
  icell_unwrapped_ = new int[nvec];
  icell_ = new int[nvec];
  coeffs_ = new int[nvec];
  // Initialize arrays
  if (shape == nullptr) {
    std::fill(shape_, shape_ + nvec, 0);
  } else {
    std::copy(shape, shape + nvec, shape_);
  }
  for (int ivec = 0; ivec < nvec_; ++ivec)
    take_range(ivec);
}


BarIterator::~BarIterator() {
  if (shape_ != nullptr) delete[] shape_;
  if (ranges_begin_ != nullptr) delete[] ranges_begin_;
  if (ranges_end_ != nullptr) delete[] ranges_end_;
  if (icell_unwrapped_ != nullptr) delete[] icell_unwrapped_;
  if (icell_ != nullptr) delete[] icell_;
  if (coeffs_ != nullptr) delete[] coeffs_;
}


BarIterator& BarIterator::operator++() {
  increment(nvec_-1);
  return *this;
}


BarIterator BarIterator::operator++(int) {
  throw std::logic_error("Don't use the post-increment operator of BarIterator.");
  return *this;
}


void BarIterator::take_range(const int ivec) {
  if (ibar_ > bars_.size() - 2)
    throw std::range_error("The bar vector is too short.");
  ranges_begin_[ivec] = bars_[ibar_];
  ranges_end_[ivec] = bars_[ibar_ + 1];
  icell_unwrapped_[ivec] = ranges_begin_[ivec];
  icell_[ivec] = robust_wrap(ranges_begin_[ivec], shape_[ivec], &coeffs_[ivec]);
  ibar_ += 2;
}


void BarIterator::increment(const int ivec) {
  ++icell_unwrapped_[ivec];
  ++icell_[ivec];
  if ((shape_[ivec] != 0) && (icell_[ivec] >= shape_[ivec])) {
    icell_[ivec] = 0;
    ++coeffs_[ivec];
  }
  if (icell_unwrapped_[ivec] >= ranges_end_[ivec]) {
    if (ivec > 0) {
      increment(ivec - 1);
      if (busy_)
        take_range(ivec);
    } else {
      busy_ = false;
      if (ibar_ != bars_.size())
        throw std::range_error("Cannot iterate past end of bar.");
    }
  }
}


}  // namespace celllists

// vim: textwidth=90 et ts=2 sw=2
