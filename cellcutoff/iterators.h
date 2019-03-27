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

/** @file */


#ifndef CELLCUTOFF_ITERATORS_H_
#define CELLCUTOFF_ITERATORS_H_

#include <vector>
#include <array>

#include "cellcutoff/cell.h"
#include "cellcutoff/decomposition.h"


namespace cellcutoff {


class BarIterator {
 public:
  BarIterator(const std::vector<int>& bars, const int nvec, const int* shape);
  BarIterator(const std::vector<int>& bars, const int nvec)
      : BarIterator(bars, nvec, nullptr) {}
  ~BarIterator();

  bool busy() const { return busy_; }
  BarIterator& operator++();
  BarIterator operator++(int);

  const int* icell() const { return icell_; }
  const int* coeffs() const { return coeffs_; }

 private:
  void take_range(const int ivec);
  void increment(const int ivec);

  const std::vector<int>& bars_;
  const int nvec_;
  size_t ibar_;
  int* shape_;
  int* ranges_begin_;
  int* ranges_end_;
  int* icell_unwrapped_;
  int* icell_;
  int* coeffs_;
  bool busy_;
};


class DeltaIterator {
 public:
  DeltaIterator(const Cell& subcell, const int* shape, const double* center,
      const double cutoff, const void* points, const size_t npoint,
      const size_t point_size, const CellMap& cell_map);
  DeltaIterator(const Cell& subcell, const double* center,
      const double cutoff, const void* points, const size_t npoint,
      const size_t point_size, const CellMap& cell_map)
      : DeltaIterator(subcell, nullptr, center, cutoff, points, npoint, point_size,
        cell_map) {}
  ~DeltaIterator();

  bool busy() const { return bar_iterator_->busy(); }
  DeltaIterator& operator++();
  DeltaIterator operator++(int);

  const double* delta() const { return delta_; }
  double distance() const { return distance_; }
  size_t ipoint() const { return ipoint_; }

 private:
  void increment(bool initialization);

  // Provided through constructor
  const Cell& subcell_;
  int* shape_;
  const double center_[3];
  const double cutoff_;
  const char* points_char_;
  const size_t npoint_;
  const size_t point_size_;
  const CellMap& cell_map_;

  // Internal data
  std::vector<int> bars_;
  BarIterator* bar_iterator_;
  const Point* point_;
  double cell_delta_[3];
  double delta_[3];
  double distance_;
  size_t ipoint_;
  size_t ibegin_;
  size_t iend_;
};


}  // namespace cellcutoff


#endif  // CELLCUTOFF_ITERATORS_H_

// vim: textwidth=90 et ts=2 sw=2
