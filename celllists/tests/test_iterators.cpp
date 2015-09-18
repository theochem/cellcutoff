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


class BarIterator3DTestP : public ::testing::TestWithParam<int> {
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


TEST(BarIterator3DTest, exceptions) {
  std::vector<int> bars{1, 2, 3};
  int shape[3]{2, 3, 4};
  EXPECT_THROW(cl::BarIterator3D bit(bars), std::domain_error);
  EXPECT_THROW(cl::BarIterator3D bit(bars, shape), std::domain_error);
  bars.push_back(5);
  cl::BarIterator3D bit(bars, shape);
  EXPECT_THROW(bit++, std::logic_error);
  ++bit;  // First increment should be OK.
  ++bit;  // Second increment should be OK.
  EXPECT_THROW(++bit, std::range_error);  // Third increment goes too far.
}


TEST(BarIterator3DTest, example) {
  const std::vector<int> bars{1, 2, -2, 1, 5, -3, 7, 11};
  cl::BarIterator3D it(bars);
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 2);
  EXPECT_EQ(it.icell()[2], -2);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 2);
  EXPECT_EQ(it.icell()[2], -1);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 2);
  EXPECT_EQ(it.icell()[2], 0);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 5);
  EXPECT_EQ(it.icell()[1], -3);
  EXPECT_EQ(it.icell()[2], 7);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 5);
  EXPECT_EQ(it.icell()[1], -3);
  EXPECT_EQ(it.icell()[2], 8);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 5);
  EXPECT_EQ(it.icell()[1], -3);
  EXPECT_EQ(it.icell()[2], 9);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 5);
  EXPECT_EQ(it.icell()[1], -3);
  EXPECT_EQ(it.icell()[2], 10);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_FALSE(it.busy());
}


TEST(BarIterator3DTest, example_shape) {
  const std::vector<int> bars{1, 2, -2, 1, 5, -3, 7, 11};
  const int shape[3]{3, 4, 5};
  cl::BarIterator3D it(bars, shape);
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 2);
  EXPECT_EQ(it.icell()[2], 3);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], -1);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 2);
  EXPECT_EQ(it.icell()[2], 4);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], -1);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 2);
  EXPECT_EQ(it.icell()[2], 0);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  EXPECT_EQ(it.coeffs()[2], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 2);
  EXPECT_EQ(it.icell()[1], 1);
  EXPECT_EQ(it.icell()[2], 2);
  EXPECT_EQ(it.coeffs()[0], 1);
  EXPECT_EQ(it.coeffs()[1], -1);
  EXPECT_EQ(it.coeffs()[2], 1);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 2);
  EXPECT_EQ(it.icell()[1], 1);
  EXPECT_EQ(it.icell()[2], 3);
  EXPECT_EQ(it.coeffs()[0], 1);
  EXPECT_EQ(it.coeffs()[1], -1);
  EXPECT_EQ(it.coeffs()[2], 1);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 2);
  EXPECT_EQ(it.icell()[1], 1);
  EXPECT_EQ(it.icell()[2], 4);
  EXPECT_EQ(it.coeffs()[0], 1);
  EXPECT_EQ(it.coeffs()[1], -1);
  EXPECT_EQ(it.coeffs()[2], 1);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 2);
  EXPECT_EQ(it.icell()[1], 1);
  EXPECT_EQ(it.icell()[2], 0);
  EXPECT_EQ(it.coeffs()[0], 1);
  EXPECT_EQ(it.coeffs()[1], -1);
  EXPECT_EQ(it.coeffs()[2], 2);
  ++it;
  EXPECT_FALSE(it.busy());
}


TEST_P(BarIterator3DTestP, example_random) {
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
    size_t nbar = subcell->bars_cutoff(center, cutoff, &bars);
    cl::BarIterator3D bit(bars, shape);
    EXPECT_EQ(nbar, bit.nbar());
    for (size_t ibar = 0; ibar < nbar; ++ibar) {
      std::array<int, 3> icell;
      int coeffs[3];
      icell[0] = cl::robust_wrap(bars[4*ibar], shape[0], &coeffs[0]);
      icell[1] = cl::robust_wrap(bars[4*ibar + 1], shape[1], &coeffs[1]);
      int begin2 = bars[4*ibar + 2];
      int end2 = bars[4*ibar + 3];
      for (int icell2 = begin2; icell2 < end2; ++icell2) {
        icell[2] = cl::robust_wrap(icell2, shape[2], &coeffs[2]);
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
    EXPECT_FALSE(bit.busy());
  }
}


// Instantiation of parameterized tests
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

INSTANTIATE_TEST_CASE_P(BarIterator3DTest0123, BarIterator3DTestP, ::testing::Range(0, 4));


// vim: textwidth=90 et ts=2 sw=2
