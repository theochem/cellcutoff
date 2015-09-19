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


class BarIterator3D {
 public:
  explicit BarIterator3D(const std::vector<int>& bars);
  BarIterator3D(const std::vector<int>& bars, const int* shape);

  bool busy() const { return ibar_ < nbar_; }
  BarIterator3D& operator++();
  BarIterator3D operator++(int);
  const std::array<int, 3> icell() const { return icell_; }
  const int* coeffs() const { return coeffs_; }
  size_t nbar() const { return nbar_; }

 private:
  const std::vector<int>& bars_;
  size_t ibar_;
  const size_t nbar_;
  std::array<int, 3> icell_;
  int shape_[3];
  int coeffs_[3];
};


}  // namespace celllists


#endif  // CELLLISTS_ITERATORS_H_

// vim: textwidth=90 et ts=2 sw=2
