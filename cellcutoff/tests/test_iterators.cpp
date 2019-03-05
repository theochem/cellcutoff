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


#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include <cellcutoff/cell.h>
#include <cellcutoff/decomposition.h>
#include <cellcutoff/iterators.h>

#include "common.h"


namespace cl = cellcutoff;


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

  // Tests on allowed nvec
  EXPECT_THROW(cl::BarIterator(bars, -1), std::domain_error);
  EXPECT_THROW(cl::BarIterator(bars, 0), std::domain_error);
  EXPECT_THROW(cl::BarIterator(bars, 4), std::domain_error);

  // Tests related to the size of the bars vector.
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


TEST(BarIteratorTest, example_1) {
  const std::vector<int> bars{1, 4};
  cl::BarIterator it(bars, 1);
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.coeffs()[0], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 2);
  EXPECT_EQ(it.coeffs()[0], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 3);
  EXPECT_EQ(it.coeffs()[0], 0);
  ++it;
  EXPECT_FALSE(it.busy());
}


TEST(BarIteratorTest, example_1_shape) {
  const std::vector<int> bars{1, 4};
  const int shape[1]{3};
  cl::BarIterator it(bars, 1, shape);
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.coeffs()[0], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 2);
  EXPECT_EQ(it.coeffs()[0], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 0);
  EXPECT_EQ(it.coeffs()[0], 1);
  ++it;
  EXPECT_FALSE(it.busy());
}


TEST(BarIteratorTest, example_2) {
  const std::vector<int> bars{
    1, 3,
      -2, 0,
       5, 7};
  cl::BarIterator it(bars, 2);
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], -2);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], -1);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 2);
  EXPECT_EQ(it.icell()[1], 5);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 2);
  EXPECT_EQ(it.icell()[1], 6);
  EXPECT_EQ(it.coeffs()[0], 0);
  EXPECT_EQ(it.coeffs()[1], 0);
  ++it;
  EXPECT_FALSE(it.busy());
}


TEST(BarIteratorTest, example_2_shape) {
  const std::vector<int> bars{
    4, 6,
      -2, 0,
       5, 7};
  const int shape[2]{3, 2};
  cl::BarIterator it(bars, 2, shape);
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 0);
  EXPECT_EQ(it.coeffs()[0], 1);
  EXPECT_EQ(it.coeffs()[1], -1);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 1);
  EXPECT_EQ(it.icell()[1], 1);
  EXPECT_EQ(it.coeffs()[0], 1);
  EXPECT_EQ(it.coeffs()[1], -1);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 2);
  EXPECT_EQ(it.icell()[1], 1);
  EXPECT_EQ(it.coeffs()[0], 1);
  EXPECT_EQ(it.coeffs()[1], 2);
  ++it;
  EXPECT_TRUE(it.busy());
  EXPECT_EQ(it.icell()[0], 2);
  EXPECT_EQ(it.icell()[1], 0);
  EXPECT_EQ(it.coeffs()[0], 1);
  EXPECT_EQ(it.coeffs()[1], 3);
  ++it;
  EXPECT_FALSE(it.busy());
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
    EXPECT_LT(-1, shape[0]);
    EXPECT_LT(-1, shape[1]);
    EXPECT_LT(-1, shape[2]);

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


// DeltaIterator
// ~~~~~~~~~~~~~

TEST(DeltaIteratorTest, exception) {
  double vecs[9]{2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell subcell(vecs, 3);
  const double center[3]{0.0, 0.0, 0.0};
  std::vector<cl::Point> points;
  cl::CellMap cell_map;
  cl::DeltaIterator dit1(subcell, nullptr, center, 1e-15, points.data(), points.size(),
      sizeof(cl::Point), cell_map);
  EXPECT_THROW(dit1++, std::logic_error);
}

TEST(DeltaIteratorTest, examples) {
  // Problem definition: points, cell and cell_map, center and cutoff
  // 1) points
  std::vector<cl::Point> points;
  double cart_0[3]{9.0, 9.0, 99.0};
  double cart_1[3]{22.5, 0.0, 0.0};
  double cart_2[3]{4.0, 5.1, -3.0};
  double cart_3[3]{4.0, 5.2, -3.0};
  points.emplace_back(cart_0);
  points.emplace_back(cart_1);
  points.emplace_back(cart_2);
  points.emplace_back(cart_3);
  // 2) cell and shape
  double vecs[9]{2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell subcell(vecs, 3);
  int shape[3]{10, 20, 20};
  // 3) cell_map (no need to sort, as they already are sorted)
  cl::assign_icell(subcell, shape, points.data(), points.size(), sizeof(cl::Point));
  std::unique_ptr<cl::CellMap> cell_map(
    cl::create_cell_map(points.data(), points.size(), sizeof(cl::Point)));
  // 4) center and cutoff
  double center[3]{1.0, 1.0, 1.0};
  double cutoff = 10.0;

  // Construct the iterator and go through every point, one by one
  cl::DeltaIterator dit1(subcell, shape, center, cutoff, points.data(), points.size(),
      sizeof(cl::Point), *cell_map);
  EXPECT_TRUE(dit1.busy());

  // First iteration
  EXPECT_NEAR(1.5, dit1.delta()[0], EPS);
  EXPECT_NEAR(-1.0, dit1.delta()[1], EPS);
  EXPECT_NEAR(-1.0, dit1.delta()[2], EPS);
  EXPECT_NEAR(sqrt(1.5*1.5 + 1.0 + 1.0), dit1.distance(), EPS);
  EXPECT_EQ(1, dit1.ipoint());
  ++dit1;
  EXPECT_TRUE(dit1.busy());
  // Second iteration
  EXPECT_NEAR(3.0, dit1.delta()[0], EPS);
  EXPECT_NEAR(4.1, dit1.delta()[1], EPS);
  EXPECT_NEAR(-4.0, dit1.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0*3.0 + 4.1*4.1 + 4.0*4.0), dit1.distance(), EPS);
  EXPECT_EQ(2, dit1.ipoint());
  ++dit1;
  EXPECT_TRUE(dit1.busy());
  // First iteration
  EXPECT_NEAR(3.0, dit1.delta()[0], EPS);
  EXPECT_NEAR(4.2, dit1.delta()[1], EPS);
  EXPECT_NEAR(-4.0, dit1.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0*3.0 + 4.2*4.2 + 4.0*4.0), dit1.distance(), EPS);
  EXPECT_EQ(3, dit1.ipoint());
  ++dit1;
  EXPECT_FALSE(dit1.busy());

  // Same thing, this time without shape argument
  points[0].cart_[0] = 0.0;
  points[0].cart_[1] = 0.0;
  points[0].cart_[2] = 0.0;
  points[1].cart_[0] = cart_1[0];
  points[2].cart_[2] = cart_2[2];
  points[3].cart_[2] = cart_3[2];
  cl::assign_icell(subcell, points.data(), points.size(), sizeof(cl::Point));
  EXPECT_EQ(11, points[1].icell_[0]);
  EXPECT_EQ(-3, points[2].icell_[2]);
  EXPECT_EQ(-3, points[3].icell_[2]);
  cell_map.reset(cl::create_cell_map(points.data(), points.size(), sizeof(cl::Point)));
  cl::DeltaIterator dit2(subcell, center, cutoff, points.data(), points.size(),
      sizeof(cl::Point), *cell_map);
  EXPECT_TRUE(dit2.busy());

  // First iteration
  EXPECT_NEAR(-1.0, dit2.delta()[0], EPS);
  EXPECT_NEAR(-1.0, dit2.delta()[1], EPS);
  EXPECT_NEAR(-1.0, dit2.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0), dit2.distance(), EPS);
  EXPECT_EQ(0, dit2.ipoint());
  ++dit2;
  EXPECT_TRUE(dit2.busy());
  // First iteration
  EXPECT_NEAR(3.0, dit2.delta()[0], EPS);
  EXPECT_NEAR(4.1, dit2.delta()[1], EPS);
  EXPECT_NEAR(-4.0, dit2.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0*3.0 + 4.1*4.1 + 4.0*4.0), dit2.distance(), EPS);
  EXPECT_EQ(2, dit2.ipoint());
  ++dit2;
  EXPECT_TRUE(dit2.busy());
  // Second iteration
  EXPECT_NEAR(3.0, dit2.delta()[0], EPS);
  EXPECT_NEAR(4.2, dit2.delta()[1], EPS);
  EXPECT_NEAR(-4.0, dit2.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0*3.0 + 4.2*4.2 + 4.0*4.0), dit2.distance(), EPS);
  EXPECT_EQ(3, dit2.ipoint());
  ++dit2;
  EXPECT_FALSE(dit2.busy());
}


// Instantiation of parameterized tests
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

INSTANTIATE_TEST_CASE_P(BarIteratorTest0123, BarIteratorTestP, ::testing::Range(0, 4));


// vim: textwidth=90 et ts=2 sw=2
