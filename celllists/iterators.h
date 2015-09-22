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


#ifndef CELLLISTS_ITERATORS_H_
#define CELLLISTS_ITERATORS_H_

#include <vector>
#include <array>


namespace celllists {


class BarIterator {
 public:
  BarIterator(const std::vector<int>& bars, const int nvec, const int* shape);
  BarIterator(const std::vector<int>& bars, const int nvec)
      : BarIterator(bars, nvec, nullptr) {};
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


}  // namespace celllists


#endif  // CELLLISTS_ITERATORS_H_

// vim: textwidth=90 et ts=2 sw=2
