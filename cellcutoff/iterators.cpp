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


#include "cellcutoff/iterators.h"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "cellcutoff/cell.h"
#include "cellcutoff/decomposition.h"
#include "cellcutoff/vec3.h"


namespace cellcutoff {


//
// Cutoff functions
//


size_t cutoff_ranges(const Cell* cell, const double* center, const double cutoff,
    int* ranges_begin, int* ranges_end) {
  if (cutoff <= 0) {
    throw std::domain_error("cutoff must be strictly positive.");
  }
  double frac[3];
  size_t ncell = 1;
  cell->to_frac(center, frac);
  for (int ivec = cell->nvec() - 1; ivec >= 0; --ivec) {
    // Use spacings between planes to find first plane before cutoff sphere and last
    // plane after cutoff sphere. To this end, we must divide cutoff by the spacing
    // between planes.
    double frac_cutoff = cutoff/cell->spacings()[ivec];
    ranges_begin[ivec] = static_cast<int>(floor(frac[ivec] - frac_cutoff));
    ranges_end[ivec] = static_cast<int>(ceil(frac[ivec] + frac_cutoff));

    ncell *= (ranges_end[ivec] - ranges_begin[ivec]);
  }
  return ncell;
}


void cutoff_bars(const Cell* cell, const double* center, const double cutoff,
    std::vector<int>* bars) {
  // Check arguments
  if (cell->nvec() == 0) {
    throw std::domain_error("The cell must be at least 1D periodic.");
  }
  if (cutoff <= 0) {
    throw std::domain_error("cutoff must be strictly positive.");
  }
  // For all the heavy work, a SphereSlice object is used, which precomputes a lot.
  SphereSlice sphere_slice(center, cell->gvecs(), cutoff);
  // Compute bars and return the number of bars
  cutoff_bars_low(cell, &sphere_slice, 0, bars);
}


void cutoff_bars_low(const Cell* cell, SphereSlice* slice, int ivec, std::vector<int>* bars) {
  // Use SphereSlice object to solve the tedious of problem of finding begin and end.
  double begin_exact = 0.0;
  double end_exact = 0.0;
  slice->solve_range(ivec, &begin_exact, &end_exact);
  int begin = static_cast<int>(floor(begin_exact));
  int end = static_cast<int>(ceil(end_exact));

  // Store the current range.
  bars->push_back(begin);
  bars->push_back(end);
  if (ivec < cell->nvec() - 1) {
    // If this is not yet the last recursion, iterate over the range of integer
    // fractional coordinates, and go one recursion deeper in each iteration.
    for (int i = begin; i < end; ++i) {
      // Define a slice (two cuts) in the sphere.
      slice->set_cut_begin_end(ivec, i, i + 1);
      // Recursive call in which the remaining details of the bar/bars is/are solved.
      cutoff_bars_low(cell, slice, ivec + 1, bars);
    }
  }
}


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


//
// DeltaIterator (DEPRECATED)
//

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
  cutoff_bars(&subcell_, center_, cutoff_, &bars_);
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


//
// BoxSortedPoints
//


size_t serialize_icell(const int* icell) {
  return serialize_icell(icell[0], icell[1], icell[2]);
}


size_t serialize_icell(const int i0, const int i1, const int i2) {
  const int small = 4*(i0 < 0) + 2*(i1 < 0) + (i2 < 0);
  const size_t x = abs(i0) + (i0 >= 0);
  const size_t y = abs(i1) + (i1 >= 0);
  const size_t z = abs(i2) + (i2 >= 0);
  const size_t d0 = x + y + z;
  const size_t d1 = x + y;
  return (((d0-3)*(d0-2)*(d0-1))/6 + ((d1-2)*(d1-1))/2 + (x-1))*8 + small;
}


BoxSortedPoints::BoxSortedPoints(const double* points, size_t npoint, const Cell* cell,
    const double threshold)
    : points_(nullptr), npoint_(npoint), subcell_(nullptr), shape_{-1, -1, -1},
      ipoints_(nullptr), ranges_() {
  if (npoint == 0) {
    throw std::domain_error("The number of points must be strictly positive.");
  }
  // Derive the subcell.
  subcell_ = cell->create_subcell(threshold, shape_);
  // Copy wrapped points and assign subcell serials for every point.
  points_ = new double[3*npoint];
  std::vector<size_t> serials(npoint);
  for (size_t ipoint = 0; ipoint < npoint; ++ipoint) {
    double cart[3];
    vec3::copy(points + 3*ipoint, cart);
    cell->iwrap_box(cart);
    vec3::copy(cart, points_ + 3*ipoint);
    double frac[3];
    subcell_->to_frac(cart, frac);
    int isubcell[3];
    std::transform(frac, frac + 3, isubcell, &floor);
    serials[ipoint] = serialize_icell(isubcell);
  }
  // Determine order for sorting.
  ipoints_ = new size_t[npoint];
  std::iota(ipoints_, ipoints_ + npoint, 0);
  std::sort(ipoints_, ipoints_ + npoint,
            [&serials](size_t i1, size_t i2) {return serials[i1] < serials[i2];});
  // Construct the ranges map, which will be used to look up points in a subcell.
  size_t begin = 0;
  size_t last_serial = serials[ipoints_[0]];
  for (size_t jpoint = 0; jpoint < npoint; ++jpoint) {
    size_t ipoint = ipoints_[jpoint];
    if (serials[ipoint] != last_serial) {
      ranges_.insert({last_serial, {begin, jpoint}});
      begin = jpoint;
      last_serial = serials[ipoint];
    }
  }
  ranges_.insert({last_serial, {begin, npoint}});
}


BoxSortedPoints::~BoxSortedPoints() {
  if (points_ != nullptr) delete[] points_;
  if (subcell_ != nullptr) delete subcell_;
  if (ipoints_ != nullptr) delete[] ipoints_;
}


//
// BoxCutoffIterator
//


BoxCutoffIterator::BoxCutoffIterator(const BoxSortedPoints* bsp, const double* center,
    double radius)
    : bsp_(bsp),
      bars_(),
      bar_iterator_(nullptr),
      center_{NAN, NAN, NAN},
      radius_(radius),
      ibegin_(1),
      iend_(1),
      icurrent_(0),
      ipoint_{0},
      cell_delta_{NAN, NAN, NAN},
      delta_{NAN, NAN, NAN},
      distance_(NAN) {
  // Create the bars vector, integers that encode the subcell indices to visit.
  cutoff_bars(bsp_->subcell(), center, radius, &bars_);
  bar_iterator_ = new BarIterator(bars_, bsp_->subcell()->nvec(), bsp_->shape());
  // Copy center
  vec3::copy(center, center_);
  // Prepare first iteration
  increment(true);
}


BoxCutoffIterator::~BoxCutoffIterator() {
  if (bar_iterator_ != nullptr) delete bar_iterator_;
}



BoxCutoffIterator& BoxCutoffIterator::operator++() {
  increment(false);
  return *this;
}


BoxCutoffIterator BoxCutoffIterator::operator++(int) {
  throw std::logic_error("Don't use the post-increment operator of BoxCutoffIterator.");
  return *this;
}


void BoxCutoffIterator::increment(bool initialization) {
  do {
    // Just move one point further
    ++icurrent_;
    // If we moved past the last point in the cell, then do the outer loop
    if (icurrent_ == iend_) {
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
        // Take all relevant data from bar_iterator_ and bsp_.
        // - the cell index, converting to serial
        size_t serial = serialize_icell(bar_iterator_->icell());
        // - the next range of points, if any. If the next range of points is not in
        //   bsp_->ranges(), the while loop will try the next cell.
        auto it = bsp_->ranges()->find(serial);
        if (it != bsp_->ranges()->end()) {
          // The range points + reset ipoint_ to the beginning if that range
          ibegin_ = it->second[0];
          iend_ = it->second[1];
          icurrent_ = ibegin_;
        }
      } while (icurrent_ == iend_);
      // When we get here, a new cell with some points is found.
      // Compute the relative vector of the cutoff center to the lower corner of the
      // periodic cell. (This is called the cell_delta_ vector, as it is, for a given
      // cell, a constant part of the relative vector from center to a point within a
      // given cell.)
      cell_delta_[0] = -center_[0];
      cell_delta_[1] = -center_[1];
      cell_delta_[2] = -center_[2];
      int translate_icell[3]{
        bar_iterator_->coeffs()[0]*bsp_->shape()[0],
        bar_iterator_->coeffs()[1]*bsp_->shape()[1],
        bar_iterator_->coeffs()[2]*bsp_->shape()[2],
      };
      bsp_->subcell()->iadd_vec(cell_delta_, translate_icell);
    }
    // When we reach this point, a new point is found, either in a new cell or not. Some
    // additional properties of that point are computed here. If the distance from the
    // center is beyond the cutoff, we just move to the next point.
    ipoint_ = bsp_->ipoints()[icurrent_];
    vec3::copy(bsp_->points() + 3*ipoint_, delta_);
    vec3::iadd(delta_, cell_delta_);
    distance_ = vec3::norm(delta_);
  } while (distance_ > radius_);
}


}  // namespace cellcutoff

// vim: textwidth=90 et ts=2 sw=2
