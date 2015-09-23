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

#include "celllists/cell.h"
#include "celllists/decomposition.h"
#include "celllists/vec3.h"


namespace celllists {


//
// BarIterator
//

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


// DeltaIterator

DeltaIterator::DeltaIterator(const Cell& subcell, const int* shape, const double* center,
    const double cutoff, const void* points, const size_t npoint, const size_t point_size,
    const CellMap& cell_map)
    : subcell_(subcell),
      shape_(nullptr),
      center_{center[0], center[1], center[2]},
      cutoff_(cutoff),
      points_char_(reinterpret_cast<const char*>(points)),
      npoint_(npoint),
      point_size_(point_size),
      cell_map_(cell_map),
      bar_iterator_(nullptr),
      point_(nullptr),
      cell_delta_{NAN, NAN, NAN},
      delta_{NAN, NAN, NAN},
      distance_(NAN),
      ipoint_(0),
      ibegin_(1),
      iend_(1) {

  int nvec = subcell_.nvec();
  shape_ = new int[nvec];
  // Initialize shape
  if (shape == nullptr) {
    std::fill(shape_, shape_ + nvec, 0);
  } else {
    std::copy(shape, shape + nvec, shape_);
  }
  // Set up the bar_iterator_
  subcell_.bars_cutoff(center_, cutoff_, &bars_);
  bar_iterator_ = new BarIterator(bars_, nvec, shape_);
  // Prepare first iteration
  increment(true);
}


DeltaIterator::~DeltaIterator() {
  if (shape_ != nullptr) delete[] shape_;
  if (bar_iterator_ != nullptr) delete bar_iterator_;
}


DeltaIterator& DeltaIterator::operator++() {
  increment(false);
  return *this;
}


DeltaIterator DeltaIterator::operator++(int) {
  throw std::logic_error("Don't use the post-increment operator of DeltaIterator.");
  return *this;
}


void DeltaIterator::increment(bool initialization) {
  do {
    // Just move one point further
    ++ipoint_;
    // If we moved past the last point in the cell, then do the outer loop
    if (ipoint_ == iend_) {
      do {
        // Move to the next cell, if any
        if (!initialization) {
          ++(*bar_iterator_);
        } else {
          initialization = false;
        }
        // Check if there is a next cell.
        if (!bar_iterator_->busy()) {
          delta_[0] = NAN;
          delta_[1] = NAN;
          delta_[2] = NAN;
          distance_ = NAN;
          return;
        }
        // Take all relevant data from bar_iterator_ and cell_map_.
        // - the cell index
        std::array<int, 3> key{
          bar_iterator_->icell()[0],
          bar_iterator_->icell()[1],
          bar_iterator_->icell()[2]};
        // - the next range of points, if any. If the next range of points is not in
        //   cell_map_, the while loop will try the next cell.
        auto it = cell_map_.find(key);
        if (it != cell_map_.end()) {
          // The range points + reset ipoint_ to the beginning if that range
          ibegin_ = it->second[0];
          iend_ = it->second[1];
          ipoint_ = ibegin_;
        }
      } while (ipoint_ == iend_);
      // When we get here, a new cell with some points is found.
      // Compute the relative vector of the cutoff center to the lower corner of the
      // periodic cell. (This is called the cell_delta_ vector, as it is, for a given
      // cell, a constant part of the relative vector from center to a point within a
      // given cell.)
      cell_delta_[0] = -center_[0];
      cell_delta_[1] = -center_[1];
      cell_delta_[2] = -center_[2];
      int translate_icell[3]{
        bar_iterator_->coeffs()[0]*shape_[0],
        bar_iterator_->coeffs()[1]*shape_[1],
        bar_iterator_->coeffs()[2]*shape_[2],
      };
      subcell_.iadd_vec(cell_delta_, translate_icell);
    }
    // When we reach this point, a new point is found, either in a new cell or not. Some
    // additional properties of that point are computed here. If the distance from the
    // center is beyond the cutoff, we just move to the next point.
    point_ = reinterpret_cast<const Point*>(points_char_ + ipoint_*point_size_);
    vec3::copy(point_->cart_, delta_);
    vec3::iadd(delta_, cell_delta_);
    distance_ = vec3::norm(delta_);
  } while (distance_ > cutoff_);
}


}  // namespace celllists

// vim: textwidth=90 et ts=2 sw=2
