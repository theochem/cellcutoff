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


#include <vector>

#include <gtest/gtest.h>

#include <celllists/cell.h>
#include <celllists/decomposition.h>
#include <celllists/iterators.h>

#include "common.h"


namespace cl = celllists;


class BarIteratorTestP : public ::testing::TestWithParam<int> {
 public:
  virtual void SetUp() {
    nvec = GetParam();
  }

  std::unique_ptr<cl::Cell> create_random_cell(const unsigned int seed,
      const double scale = 1.0, const double ratio = 0.1, const bool cuboid = false) {
    return create_random_cell_nvec(seed, nvec, scale, ratio, cuboid);
  }

  int nvec;
};


TEST(BarIteratorTest, exceptions) {
  std::vector<int> bars{1, 2, 3};
  int shape[3]{2, 3, 4};

  cl::BarIterator bit1(bars, 1);
  EXPECT_THROW(bit1++, std::logic_error);  // No post increment allowed
  EXPECT_THROW(++bit1, std::range_error);  // Cannot increment as ranges is too short

  cl::BarIterator bit1s(bars, 1, shape);
  EXPECT_THROW(bit1s++, std::logic_error);  // No post increment allowed
  EXPECT_THROW(++bit1s, std::range_error);  // Cannot increment as ranges is too short

  EXPECT_THROW(cl::BarIterator bit2(bars, 2), std::range_error);
  EXPECT_THROW(cl::BarIterator bit2s(bars, 2, shape), std::range_error);
  EXPECT_THROW(cl::BarIterator bit3(bars, 3), std::range_error);
  EXPECT_THROW(cl::BarIterator bit3s(bars, 3, shape), std::range_error);
}


TEST(BarIteratorTest, example_3) {
  const std::vector<int> bars{
    1, 2,
      -2, 0,
        5, 7,
        8, 10};
  cl::BarIterator it(bars, 3);
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], -2);
  EXPECT_EQ(it.icell()[2], 5);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], -2);
  EXPECT_EQ(it.icell()[2], 6);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], -1);
  EXPECT_EQ(it.icell()[2], 8);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], -1);
  EXPECT_EQ(it.icell()[2], 9);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_FALSE(it.busy());
}


TEST(BarIteratorTest, example_3_shape) {
  const std::vector<int> bars{
    1, 2,
      -2, 0,
        5, 7,
        8, 10};
  const int shape[3]{3, 4, 5};
  cl::BarIterator it(bars, 3, shape);
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 2);
  EXPECT_EQ(it.icell()[2], 0);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], -1);
  EXPECT_EQ(it.coeffs()[2], 1);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 2);
  EXPECT_EQ(it.icell()[2], 1);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], -1);
  EXPECT_EQ(it.coeffs()[2], 1);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 3);
  EXPECT_EQ(it.icell()[2], 3);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], -1);
  EXPECT_EQ(it.coeffs()[2], 1);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 3);
  EXPECT_EQ(it.icell()[2], 4);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], -1);
  EXPECT_EQ(it.coeffs()[2], 1);
  ++it;
  EXPECT_FALSE(it.busy());
}


TEST_P(BarIteratorTestP, example_3_random) {
  for (int irep = 0; irep < NREP; ++irep) {
    // Problem definition
    double cutoff = 1.0 + static_cast<double>(irep)/NREP;
    double center[3];
    fill_random_double(irep, center, 3, -cutoff, cutoff);

    std::unique_ptr<cl::Cell> cell(create_random_cell(irep+NREP, cutoff, 0.6, false));
    cell->iwrap_box(center);
    int shape[3] = {-1, -1, -1};
    std::unique_ptr<cl::Cell> subcell(cell->create_subcell(cutoff*0.2, shape));

    std::vector<int> bars;
    subcell->bars_cutoff(center, cutoff, &bars);
    cl::BarIterator bit(bars, 3, shape);

    int ibar = 2;
    int icell[3];
    int coeffs[3];
    for (int ifrac0 = bars[0]; ifrac0 < bars[1]; ++ifrac0) {
      icell[0] = cl::robust_wrap(ifrac0, shape[0], &coeffs[0]);
      int begin1 = bars[ibar];
      int end1 = bars[ibar+1];
      ibar += 2;
      for (int ifrac1 = begin1; ifrac1 < end1; ++ifrac1) {
        icell[1] = cl::robust_wrap(ifrac1, shape[1], &coeffs[1]);
        int begin2 = bars[ibar];
        int end2 = bars[ibar+1];
        ibar += 2;
        for (int ifrac2 = begin2; ifrac2 < end2; ++ifrac2) {
          icell[2] = cl::robust_wrap(ifrac2, shape[2], &coeffs[2]);
          EXPECT_EQ(icell[0], bit.icell()[0]);
          EXPECT_EQ(icell[1], bit.icell()[1]);
          EXPECT_EQ(icell[2], bit.icell()[2]);
          EXPECT_EQ(coeffs[0], bit.coeffs()[0]);
          EXPECT_EQ(coeffs[1], bit.coeffs()[1]);
          EXPECT_EQ(coeffs[2], bit.coeffs()[2]);
          EXPECT_TRUE(bit.busy());
          ++bit;
        }
      }
    }
    EXPECT_FALSE(bit.busy());
  }
}


// Instantiation of parameterized tests
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

INSTANTIATE_TEST_CASE_P(BarIteratorTest0123, BarIteratorTestP, ::testing::Range(0, 4));


// vim: textwidth=90 et ts=2 sw=2
