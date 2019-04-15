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


#include <cmath>

#include <vector>

#include <gtest/gtest.h>

#include <cellcutoff/cell.h>
#include <cellcutoff/decomposition.h>
#include <cellcutoff/iterators.h>
#include <cellcutoff/vec3.h>

#include "cellcutoff/tests/common.h"


namespace vec3 = cellcutoff::vec3;

// Cutoff functions
// ~~~~~~~~~~~~~---


class CutoffTestP : public ::testing::TestWithParam<int> {
 public:
  virtual void SetUp() {
    nvec = GetParam();
  }

  std::unique_ptr<cl::Cell> create_random_cell(const unsigned int seed,
      const double scale = 1.0, const double ratio = 0.1, const bool cuboid = false) {
    return std::unique_ptr<cl::Cell>(cl::create_random_cell(seed, nvec, scale, ratio, cuboid));
  }

  int nvec;
};


TEST_F(CellTest1, ranges_cutoff_example) {
  double center[3] = {6.3, 0.2, -0.8};
  int ranges_begin[1];
  int ranges_end[1];
  size_t ncell = 0;
  ncell = ranges_cutoff(mycell.get(), center, 1.0, ranges_begin, ranges_end);
  EXPECT_EQ(2, ncell);
  EXPECT_EQ(2, ranges_begin[0]);
  EXPECT_EQ(4, ranges_end[0]);
  ncell = ranges_cutoff(mycell.get(), center, 2.0, ranges_begin, ranges_end);
  EXPECT_EQ(3, ncell);
  EXPECT_EQ(2, ranges_begin[0]);
  EXPECT_EQ(5, ranges_end[0]);
  ncell = ranges_cutoff(mycell.get(), center, 3.0, ranges_begin, ranges_end);
  EXPECT_EQ(4, ncell);
  EXPECT_EQ(1, ranges_begin[0]);
  EXPECT_EQ(5, ranges_end[0]);
}


TEST_F(CellTest1, ranges_cutoff_edge) {
  double center[3] = {2.0, 0.2, -0.8};
  int ranges_begin[1];
  int ranges_end[1];
  size_t ncell = 0;
  ncell = ranges_cutoff(mycell.get(), center, 1.0, ranges_begin, ranges_end);
  EXPECT_EQ(2, ncell);
  EXPECT_EQ(0, ranges_begin[0]);
  EXPECT_EQ(2, ranges_end[0]);
  ncell = ranges_cutoff(mycell.get(), center, 2.0, ranges_begin, ranges_end);
  EXPECT_EQ(2, ncell);
  EXPECT_EQ(0, ranges_begin[0]);
  EXPECT_EQ(2, ranges_end[0]);
  ncell = ranges_cutoff(mycell.get(), center, 3.0, ranges_begin, ranges_end);
  EXPECT_EQ(4, ncell);
  EXPECT_EQ(-1, ranges_begin[0]);
  EXPECT_EQ(3, ranges_end[0]);
}


TEST_F(CellTest2, ranges_cutoff_example) {
  double center[3] = {6.3, 0.2, -5.0};
  int ranges_begin[2];
  int ranges_end[2];
  size_t ncell = 0;
  ncell = ranges_cutoff(mycell.get(), center, 1.1, ranges_begin, ranges_end);
  EXPECT_EQ(2*2, ncell);
  EXPECT_EQ(2, ranges_begin[0]);
  EXPECT_EQ(-2, ranges_begin[1]);
  EXPECT_EQ(4, ranges_end[0]);
  EXPECT_EQ(0, ranges_end[1]);
}


TEST_F(CellTest2, ranges_cutoff_edge) {
  double center[3] = {4.0, 0.2, -2.0};
  int ranges_begin[2];
  int ranges_end[2];
  size_t ncell = 0;
  ncell = ranges_cutoff(mycell.get(), center, 2.0, ranges_begin, ranges_end);
  EXPECT_EQ(2, ncell);
  EXPECT_EQ(1, ranges_begin[0]);
  EXPECT_EQ(-1, ranges_begin[1]);
  EXPECT_EQ(3, ranges_end[0]);
  EXPECT_EQ(0, ranges_end[1]);
}


TEST_F(CellTest3, ranges_cutoff_example) {
  double center[3] = {6.3, 2.2, -5.8};
  int ranges_begin[3];
  int ranges_end[3];
  size_t ncell = 0;
  ncell = ranges_cutoff(mycell.get(), center, 1.0, ranges_begin, ranges_end);
  EXPECT_EQ(2*3*1, ncell);
  EXPECT_EQ(2, ranges_begin[0]);
  EXPECT_EQ(1, ranges_begin[1]);
  EXPECT_EQ(-2, ranges_begin[2]);
  EXPECT_EQ(4, ranges_end[0]);
  EXPECT_EQ(4, ranges_end[1]);
  EXPECT_EQ(-1, ranges_end[2]);
}


TEST_F(CellTest3, ranges_cutoff_edge) {
  double center[3] = {10.0, -2.0, -6.0};
  int ranges_begin[3];
  int ranges_end[3];
  size_t ncell = 0;
  ncell = ranges_cutoff(mycell.get(), center, 2.0, ranges_begin, ranges_end);
  EXPECT_EQ(2*4*1, ncell);
  EXPECT_EQ(4, ranges_begin[0]);
  EXPECT_EQ(-4, ranges_begin[1]);
  EXPECT_EQ(-2, ranges_begin[2]);
  EXPECT_EQ(6, ranges_end[0]);
  EXPECT_EQ(0, ranges_end[1]);
  EXPECT_EQ(-1, ranges_end[2]);
}


TEST_P(CellTestP, ranges_cutoff_domain) {
  double center[3] = {6.3, 2.2, -5.8};
  int ranges_begin[3];
  int ranges_end[3];
  EXPECT_THROW(ranges_cutoff(mycell.get(), center, -1.0, ranges_begin, ranges_end),
               std::domain_error);
  EXPECT_THROW(ranges_cutoff(mycell.get(), center, 0.0, ranges_begin, ranges_end),
               std::domain_error);
}


TEST_P(CellTestP, ranges_cutoff_random) {
  int npoint_total = 0;
  for (int irep = 0; irep < NREP; ++irep) {
    std::unique_ptr<cl::Cell> cell(create_random_cell(irep));
    double center[3];
    int ranges_begin[3];
    int ranges_end[3];
    double cutoff = 0.3*(irep + 1);
    fill_random_double(irep + 2, center, 3, -5.0, 5.0);
    ranges_cutoff(cell.get(), center, cutoff, ranges_begin, ranges_end);
    for (int ipoint=0; ipoint < NPOINT; ++ipoint) {
      double point[3];
      double norm;
      random_point(ipoint + irep*NPOINT, center, cutoff, point, &norm);
      if (norm <= cutoff) {
        double frac[3];
        cell->to_frac(point, frac);
        for (int ivec=0; ivec < nvec; ++ivec) {
          EXPECT_LE(ranges_begin[ivec], frac[ivec]);
          EXPECT_GE(ranges_end[ivec], frac[ivec]);
        }
        ++npoint_total;
      }
    }
  }
  // Check sufficiency
  EXPECT_LT((NREP*NPOINT)/3, npoint_total);
}


// bars_cutoff
// ~~~~~~~~~~~

TEST_P(CellTestP, bars_cutoff_domain) {
  double center[3] = {2.5, 3.4, -0.6};
  std::vector<int> bars;
  EXPECT_THROW(bars_cutoff(mycell.get(), center, 0.0, &bars), std::domain_error);
  EXPECT_THROW(bars_cutoff(mycell.get(), center, -1.0, &bars), std::domain_error);
  cl::Cell zero_cell(nullptr, 0);
  EXPECT_THROW(bars_cutoff(&zero_cell, center, 1.0, &bars), std::domain_error);
}


TEST_F(CellTest1, bars_cutoff_example) {
  // All the parameters
  double cutoff = 5.0;
  double center[3] = {2.5, 3.4, -0.6};

  // Call
  std::vector<int> bars;
  bars_cutoff(mycell.get(), center, cutoff, &bars);
  EXPECT_EQ(2, bars.size());

  // Check results
  // lower end: -2.5 (-2 #> 8)
  // upper end:  7.5 (8 #> 4) non-inclusive
  EXPECT_EQ(-2, bars[0]);
  EXPECT_EQ(4, bars[1]);
}


TEST_F(CellTest2, bars_cutoff_example) {
  // All the parameters
  double cutoff = 5.0;
  double center[3] = {2.5, 3.4, -0.6};

  // Call
  std::vector<int> bars;
  bars_cutoff(mycell.get(), center, cutoff, &bars);
  EXPECT_EQ((1+6)*2, bars.size());

  // Test
  EXPECT_EQ(-2, bars[0]);
  EXPECT_EQ(4, bars[1]);
  EXPECT_EQ(-1, bars[2]);
  EXPECT_EQ(1, bars[3]);
  EXPECT_EQ(-2, bars[4]);
  EXPECT_EQ(1, bars[5]);
  EXPECT_EQ(-2, bars[6]);
  EXPECT_EQ(2, bars[7]);
  EXPECT_EQ(-2, bars[8]);
  EXPECT_EQ(2, bars[9]);
  EXPECT_EQ(-2, bars[10]);
  EXPECT_EQ(2, bars[11]);
  EXPECT_EQ(-2, bars[12]);
  EXPECT_EQ(1, bars[13]);
}


TEST_F(CellTest3, bars_cutoff_example) {
  // All the parameters
  double cutoff = 1.9;
  double center[3] = {2.0, 2.0, 2.0};

  // Call
  std::vector<int> bars;
  bars_cutoff(mycell.get(), center, cutoff, &bars);
  EXPECT_EQ((1+1+4+1+4)*2, bars.size());

  // Test
  EXPECT_EQ(0, bars[0]);
  EXPECT_EQ(2, bars[1]);
  EXPECT_EQ(0, bars[2]);
  EXPECT_EQ(4, bars[3]);
  EXPECT_EQ(0, bars[4]);
  EXPECT_EQ(1, bars[5]);
  EXPECT_EQ(0, bars[6]);
  EXPECT_EQ(1, bars[7]);
  EXPECT_EQ(0, bars[8]);
  EXPECT_EQ(1, bars[9]);
  EXPECT_EQ(0, bars[10]);
  EXPECT_EQ(1, bars[11]);
  EXPECT_EQ(0, bars[12]);
  EXPECT_EQ(4, bars[13]);
  EXPECT_EQ(0, bars[14]);
  EXPECT_EQ(1, bars[15]);
  EXPECT_EQ(0, bars[16]);
  EXPECT_EQ(1, bars[17]);
  EXPECT_EQ(0, bars[18]);
  EXPECT_EQ(1, bars[19]);
  EXPECT_EQ(0, bars[20]);
  EXPECT_EQ(1, bars[21]);
}


TEST_P(CellTestP, bars_cutoff_random) {
  size_t ncell_total = 0;
  for (int irep = 0; irep < NREP; ++irep) {
    // Test parameters:
    // - Random cell
    std::unique_ptr<cl::Cell> cell(create_random_cell(2*irep));
    // - Increasing cutoff
    double cutoff = (irep + 1)*0.1;
    // - Random center
    double center[3];
    fill_random_double(47332 + irep, center, 3, -1.0, 1.0);

    // Compute the bars.
    std::vector<int> bars;
    bars_cutoff(cell.get(), center, cutoff, &bars);

    // Construct a random vector in a cubic box around the cutoff sphere.
    double cart[3];
    fill_random_double(123 + irep, cart, 3, -cutoff*1.1, cutoff*1.1);
    double norm = vec3::norm(cart);
    // Center of the box must coincide with center of the sphere.
    cart[0] += center[0];
    cart[1] += center[1];
    cart[2] += center[2];
    // For the rest of the test, we need this random vector in fractional coordinates.
    double frac[3];
    cell->to_frac(cart, frac);

    // Does the fractional coordinate fit in one of the bars?
    int index[3] = {
      static_cast<int>(floor(frac[0])),
      static_cast<int>(floor(frac[1])),
      static_cast<int>(floor(frac[2]))
    };
    bool in_bar = false;
    for (cl::BarIterator bit(bars, nvec); bit.busy(); ++bit) {
      in_bar = true;
      for (int ivec = 0; ivec < nvec; ++ivec) {
        in_bar &= (bit.icell()[ivec] == index[ivec]);
      }
      if (in_bar) break;
      ++ncell_total;
    }

    // Does the relative vector sit in the cutoff sphere, taking into account
    // non-periodic boundaries that truncate the cutoff sphere.
    bool in_sphere = (norm < cutoff);

    if (in_sphere) {
      // First test: if the vector is in the cutoff and non-periodic boundaries, the
      //             point must be in a bar.
      EXPECT_TRUE(in_bar);
    } else if (!in_bar) {
      // Second test: if not in a bar, the norm must be larger than the cutoff,
      //              or the point is outside a non-periodic boundary.
      EXPECT_FALSE(in_sphere);
    }
    // Clean up
  }
  // Sufficiency check
  EXPECT_LE(NREP*((nvec - 1)*3 + 1), ncell_total);
}


TEST_F(CellTest1, bars_cutoff_corners) {
  for (int irep = 0; irep < NREP; ++irep) {
    // Test parameters:
    // - Random cell
    std::unique_ptr<cl::Cell> cell(create_random_cell(2*irep));
    // - Increasing cutoff
    double cutoff = (irep + 1)*0.1;
    // - Random center
    double center[3];
    fill_random_double(47332 + irep, center, 3, -2.0, 2.0);

    // Compute the bars.
    std::vector<int> bars;
    bars_cutoff(cell.get(), center, cutoff, &bars);
    EXPECT_EQ(2, bars.size());

    // Check the ranges
    double frac_corner[3] = {0, 0, 0};
    double cart_corner[3] = {0, 0, 0};
    cell->to_frac(center, frac_corner);
    frac_corner[0] = bars[0];
    cell->to_cart(frac_corner, cart_corner);
    EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
    frac_corner[0] = bars[1];
    cell->to_cart(frac_corner, cart_corner);
    EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
  }
}


TEST_F(CellTest2, bars_cutoff_corners) {
  size_t nbar_total = 0;
  for (int irep = 0; irep < NREP; ++irep) {
    // Test parameters:
    // - Random cell
    std::unique_ptr<cl::Cell> cell(create_random_cell(2*irep));
    // - Increasing cutoff
    double cutoff = (irep + 1)*0.1;
    // - Random center
    double center[3];
    fill_random_double(47332 + irep, center, 3, -2.0, 2.0);

    // Compute the bars.
    std::vector<int> bars;
    bars_cutoff(cell.get(), center, cutoff, &bars);
    EXPECT_LE(4, bars.size());
    nbar_total += bars.size();

    // Check the ranges
    double frac_corner[3] = {0, 0, 0};
    double cart_corner[3] = {0, 0, 0};
    cell->to_frac(center, frac_corner);
    int ibar = 2;
    for (int ifrac0 = bars[0]; ifrac0 < bars[1]; ++ifrac0) {
      // corner 0,0
      frac_corner[0] = ifrac0;
      frac_corner[1] = bars[ibar];
      cell->to_cart(frac_corner, cart_corner);
      EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
      // corner 1,0
      frac_corner[0] = ifrac0 + 1;
      frac_corner[1] = bars[ibar];
      cell->to_cart(frac_corner, cart_corner);
      EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
      // corner 0,1
      frac_corner[0] = ifrac0;
      frac_corner[1] = bars[ibar + 1];
      cell->to_cart(frac_corner, cart_corner);
      EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
      // corner 1,1
      frac_corner[0] = ifrac0 + 1;
      frac_corner[1] = bars[ibar + 1];
      cell->to_cart(frac_corner, cart_corner);
      EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
      ibar += 2;
    }
  }
  // Sufficiency check
  EXPECT_LE(NREP*10, nbar_total);
}


TEST_F(CellTest3, bars_cutoff_corners) {
  size_t nbar_total = 0;
  for (int irep = 0; irep < NREP; ++irep) {
    // Test parameters:
    // - Random cell
    std::unique_ptr<cl::Cell> cell(create_random_cell(2*irep));
    // - Increasing cutoff
    double cutoff = (irep + 1)*0.1;
    // - Random center
    double center[3];
    fill_random_double(47332 + irep, center, 3, -2.0, 2.0);

    // Compute the bars.
    std::vector<int> bars;
    bars_cutoff(cell.get(), center, cutoff, &bars);
    EXPECT_LE(6, bars.size());
    nbar_total += bars.size();

    // Check the ranges
    double frac_corner[3] = {0, 0, 0};
    double cart_corner[3] = {0, 0, 0};
    cell->to_frac(center, frac_corner);
    int ibar = 2;
    for (int ifrac0 = bars[0]; ifrac0 < bars[1]; ++ifrac0) {
      int begin1 = bars[ibar];
      int end1 = bars[ibar+1];
      ibar += 2;
      for (int ifrac1 = begin1; ifrac1 < end1; ++ifrac1) {
        // corner 0,0,0
        frac_corner[0] = ifrac0;
        frac_corner[1] = ifrac1;
        frac_corner[2] = bars[ibar];
        cell->to_cart(frac_corner, cart_corner);
        EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
        // corner 1,0,0
        frac_corner[0] = ifrac0 + 1;
        frac_corner[1] = ifrac1;
        frac_corner[2] = bars[ibar];
        cell->to_cart(frac_corner, cart_corner);
        EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
        // corner 0,1,0
        frac_corner[0] = ifrac0;
        frac_corner[1] = ifrac1 + 1;
        frac_corner[2] = bars[ibar];
        cell->to_cart(frac_corner, cart_corner);
        EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
        // corner 1,1,0
        frac_corner[0] = ifrac0 + 1;
        frac_corner[1] = ifrac1 + 1;
        frac_corner[2] = bars[ibar];
        cell->to_cart(frac_corner, cart_corner);
        EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
        // corner 0,0,1
        frac_corner[0] = ifrac0;
        frac_corner[1] = ifrac1;
        frac_corner[2] = bars[ibar + 1];
        cell->to_cart(frac_corner, cart_corner);
        EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
        // corner 1,0,1
        frac_corner[0] = ifrac0 + 1;
        frac_corner[1] = ifrac1;
        frac_corner[2] = bars[ibar + 1];
        cell->to_cart(frac_corner, cart_corner);
        EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
        // corner 0,1,1
        frac_corner[0] = ifrac0;
        frac_corner[1] = ifrac1 + 1;
        frac_corner[2] = bars[ibar + 1];
        cell->to_cart(frac_corner, cart_corner);
        EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
        // corner 1,1,1
        frac_corner[0] = ifrac0 + 1;
        frac_corner[1] = ifrac1 + 1;
        frac_corner[2] = bars[ibar + 1];
        cell->to_cart(frac_corner, cart_corner);
        EXPECT_LE(cutoff, vec3::distance(cart_corner, center));
        ibar += 2;
      }
    }
  }
  // Sufficiency check
  EXPECT_LE(NREP*14, nbar_total);
}


// BarIterator
// -----------

class BarIteratorTestP : public ::testing::TestWithParam<int> {
 public:
  virtual void SetUp() {
    nvec = GetParam();
  }

  std::unique_ptr<cl::Cell> create_random_cell(const unsigned int seed,
      const double scale = 1.0, const double ratio = 0.1, const bool cuboid = false) {
    return std::unique_ptr<cl::Cell>(cl::create_random_cell(seed, nvec, scale, ratio, cuboid));
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
    bars_cutoff(subcell.get(), center, cutoff, &bars);
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
