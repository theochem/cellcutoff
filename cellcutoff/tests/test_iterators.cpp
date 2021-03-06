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
// ~~~~~~~~~~~~~~~~


class RangesCutoffTestP : public CellTestP {};
class BarsCutoffTestP : public CellTestP {};
class CutoffTest1 : public CellTest1 {};
class CutoffTest2 : public CellTest2 {};
class CutoffTest3 : public CellTest3 {};


TEST_F(CutoffTest1, cutoff_ranges_example) {
  double center[3] = {6.3, 0.2, -0.8};
  int ranges_begin[1];
  int ranges_end[1];
  size_t ncell = 0;
  ncell = cutoff_ranges(mycell.get(), center, 1.0, ranges_begin, ranges_end);
  EXPECT_EQ(2, ncell);
  EXPECT_EQ(2, ranges_begin[0]);
  EXPECT_EQ(4, ranges_end[0]);
  ncell = cutoff_ranges(mycell.get(), center, 2.0, ranges_begin, ranges_end);
  EXPECT_EQ(3, ncell);
  EXPECT_EQ(2, ranges_begin[0]);
  EXPECT_EQ(5, ranges_end[0]);
  ncell = cutoff_ranges(mycell.get(), center, 3.0, ranges_begin, ranges_end);
  EXPECT_EQ(4, ncell);
  EXPECT_EQ(1, ranges_begin[0]);
  EXPECT_EQ(5, ranges_end[0]);
}


TEST_F(CutoffTest1, cutoff_ranges_edge) {
  double center[3] = {2.0, 0.2, -0.8};
  int ranges_begin[1];
  int ranges_end[1];
  size_t ncell = 0;
  ncell = cutoff_ranges(mycell.get(), center, 1.0, ranges_begin, ranges_end);
  EXPECT_EQ(2, ncell);
  EXPECT_EQ(0, ranges_begin[0]);
  EXPECT_EQ(2, ranges_end[0]);
  ncell = cutoff_ranges(mycell.get(), center, 2.0, ranges_begin, ranges_end);
  EXPECT_EQ(2, ncell);
  EXPECT_EQ(0, ranges_begin[0]);
  EXPECT_EQ(2, ranges_end[0]);
  ncell = cutoff_ranges(mycell.get(), center, 3.0, ranges_begin, ranges_end);
  EXPECT_EQ(4, ncell);
  EXPECT_EQ(-1, ranges_begin[0]);
  EXPECT_EQ(3, ranges_end[0]);
}


TEST_F(CutoffTest2, cutoff_ranges_example) {
  double center[3] = {6.3, 0.2, -5.0};
  int ranges_begin[2];
  int ranges_end[2];
  size_t ncell = 0;
  ncell = cutoff_ranges(mycell.get(), center, 1.1, ranges_begin, ranges_end);
  EXPECT_EQ(2*2, ncell);
  EXPECT_EQ(2, ranges_begin[0]);
  EXPECT_EQ(-2, ranges_begin[1]);
  EXPECT_EQ(4, ranges_end[0]);
  EXPECT_EQ(0, ranges_end[1]);
}


TEST_F(CutoffTest2, cutoff_ranges_edge) {
  double center[3] = {4.0, 0.2, -2.0};
  int ranges_begin[2];
  int ranges_end[2];
  size_t ncell = 0;
  ncell = cutoff_ranges(mycell.get(), center, 2.0, ranges_begin, ranges_end);
  EXPECT_EQ(2, ncell);
  EXPECT_EQ(1, ranges_begin[0]);
  EXPECT_EQ(-1, ranges_begin[1]);
  EXPECT_EQ(3, ranges_end[0]);
  EXPECT_EQ(0, ranges_end[1]);
}


TEST_F(CutoffTest3, cutoff_ranges_example) {
  double center[3] = {6.3, 2.2, -5.8};
  int ranges_begin[3];
  int ranges_end[3];
  size_t ncell = 0;
  ncell = cutoff_ranges(mycell.get(), center, 1.0, ranges_begin, ranges_end);
  EXPECT_EQ(2*3*1, ncell);
  EXPECT_EQ(2, ranges_begin[0]);
  EXPECT_EQ(1, ranges_begin[1]);
  EXPECT_EQ(-2, ranges_begin[2]);
  EXPECT_EQ(4, ranges_end[0]);
  EXPECT_EQ(4, ranges_end[1]);
  EXPECT_EQ(-1, ranges_end[2]);
}


TEST_F(CutoffTest3, cutoff_ranges_edge) {
  double center[3] = {10.0, -2.0, -6.0};
  int ranges_begin[3];
  int ranges_end[3];
  size_t ncell = 0;
  ncell = cutoff_ranges(mycell.get(), center, 2.0, ranges_begin, ranges_end);
  EXPECT_EQ(2*4*1, ncell);
  EXPECT_EQ(4, ranges_begin[0]);
  EXPECT_EQ(-4, ranges_begin[1]);
  EXPECT_EQ(-2, ranges_begin[2]);
  EXPECT_EQ(6, ranges_end[0]);
  EXPECT_EQ(0, ranges_end[1]);
  EXPECT_EQ(-1, ranges_end[2]);
}


TEST_P(RangesCutoffTestP, cutoff_ranges_domain) {
  double center[3] = {6.3, 2.2, -5.8};
  int ranges_begin[3];
  int ranges_end[3];
  EXPECT_THROW(cutoff_ranges(mycell.get(), center, -1.0, ranges_begin, ranges_end),
               std::domain_error);
  EXPECT_THROW(cutoff_ranges(mycell.get(), center, 0.0, ranges_begin, ranges_end),
               std::domain_error);
}


TEST_P(RangesCutoffTestP, cutoff_ranges_random) {
  int npoint_total = 0;
  for (int irep = 0; irep < NREP; ++irep) {
    std::unique_ptr<cl::Cell> cell(create_random_cell(irep));
    double center[3];
    int ranges_begin[3];
    int ranges_end[3];
    double cutoff = 0.3*(irep + 1);
    fill_random_double(irep + 2, center, 3, -5.0, 5.0);
    cutoff_ranges(cell.get(), center, cutoff, ranges_begin, ranges_end);
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
    // Check of the old API, will be removed from 1.0.
    int ranges_begin_bis[3];
    int ranges_end_bis[3];
    cell->ranges_cutoff(center, cutoff, ranges_begin_bis, ranges_end_bis);
    for (int ivec=0; ivec < nvec; ++ivec) {
      EXPECT_EQ(ranges_begin[ivec], ranges_begin_bis[ivec]);
      EXPECT_EQ(ranges_end[ivec], ranges_end_bis[ivec]);
    }
  }
  // Check sufficiency
  EXPECT_LT((NREP*NPOINT)/3, npoint_total);
}


// cutoff_bars
// ~~~~~~~~~~~

TEST_P(BarsCutoffTestP, cutoff_bars_domain) {
  double center[3] = {2.5, 3.4, -0.6};
  std::vector<int> bars;
  EXPECT_THROW(cutoff_bars(mycell.get(), center, 0.0, &bars), std::domain_error);
  EXPECT_THROW(cutoff_bars(mycell.get(), center, -1.0, &bars), std::domain_error);
  cl::Cell zero_cell(nullptr, 0);
  EXPECT_THROW(cutoff_bars(&zero_cell, center, 1.0, &bars), std::domain_error);
}


TEST_F(CutoffTest1, cutoff_bars_example) {
  // All the parameters
  double cutoff = 5.0;
  double center[3] = {2.5, 3.4, -0.6};

  // Call
  std::vector<int> bars;
  cutoff_bars(mycell.get(), center, cutoff, &bars);
  EXPECT_EQ(2, bars.size());

  // Check results
  // lower end: -2.5 (-2 #> 8)
  // upper end:  7.5 (8 #> 4) non-inclusive
  EXPECT_EQ(-2, bars[0]);
  EXPECT_EQ(4, bars[1]);
}


TEST_F(CutoffTest2, cutoff_bars_example) {
  // All the parameters
  double cutoff = 5.0;
  double center[3] = {2.5, 3.4, -0.6};

  // Call
  std::vector<int> bars;
  cutoff_bars(mycell.get(), center, cutoff, &bars);
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


TEST_F(CutoffTest3, cutoff_bars_example) {
  // All the parameters
  double cutoff = 1.9;
  double center[3] = {2.0, 2.0, 2.0};

  // Call
  std::vector<int> bars;
  cutoff_bars(mycell.get(), center, cutoff, &bars);
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


TEST_P(BarsCutoffTestP, cutoff_bars_random) {
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
    cutoff_bars(cell.get(), center, cutoff, &bars);

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

    // Check of the old API, will be removed from 1.0.
    std::vector<int> bars_bis;
    cell->bars_cutoff(center, cutoff, &bars_bis);
    EXPECT_EQ(bars.size(), bars_bis.size());
    EXPECT_EQ(bars, bars_bis);
  }
  // Sufficiency check
  EXPECT_LE(NREP*((nvec - 1)*3 + 1), ncell_total);
}


TEST_F(CutoffTest1, cutoff_bars_corners) {
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
    cutoff_bars(cell.get(), center, cutoff, &bars);
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


TEST_F(CutoffTest2, cutoff_bars_corners) {
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
    cutoff_bars(cell.get(), center, cutoff, &bars);
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


TEST_F(CutoffTest3, cutoff_bars_corners) {
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
    cutoff_bars(cell.get(), center, cutoff, &bars);
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
    cutoff_bars(subcell.get(), center, cutoff, &bars);
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
  // Third iteration
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
  // Second iteration
  EXPECT_NEAR(3.0, dit2.delta()[0], EPS);
  EXPECT_NEAR(4.1, dit2.delta()[1], EPS);
  EXPECT_NEAR(-4.0, dit2.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0*3.0 + 4.1*4.1 + 4.0*4.0), dit2.distance(), EPS);
  EXPECT_EQ(2, dit2.ipoint());
  ++dit2;
  EXPECT_TRUE(dit2.busy());
  // Third iteration
  EXPECT_NEAR(3.0, dit2.delta()[0], EPS);
  EXPECT_NEAR(4.2, dit2.delta()[1], EPS);
  EXPECT_NEAR(-4.0, dit2.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0*3.0 + 4.2*4.2 + 4.0*4.0), dit2.distance(), EPS);
  EXPECT_EQ(3, dit2.ipoint());
  ++dit2;
  EXPECT_FALSE(dit2.busy());
}


// serialize_icell
// ~~~~~~~~~~~~~~~

TEST(SerializeICell, examples) {
  EXPECT_EQ(0, cl::serialize_icell(0, 0, 0));
  EXPECT_EQ(1, cl::serialize_icell(0, 0, -1));
  EXPECT_EQ(2, cl::serialize_icell(0, -1, 0));
  EXPECT_EQ(3, cl::serialize_icell(0, -1, -1));
  EXPECT_EQ(4, cl::serialize_icell(-1, 0, 0));
  EXPECT_EQ(5, cl::serialize_icell(-1, 0, -1));
  EXPECT_EQ(6, cl::serialize_icell(-1, -1, 0));
  EXPECT_EQ(7, cl::serialize_icell(-1, -1, -1));
  EXPECT_EQ(8, cl::serialize_icell(0, 0, 1));
  EXPECT_EQ(16, cl::serialize_icell(0, 1, 0));
  EXPECT_EQ(24, cl::serialize_icell(1, 0, 0));
  EXPECT_EQ(25, cl::serialize_icell(1, 0, -1));
  EXPECT_EQ(26, cl::serialize_icell(1, -1, 0));
  EXPECT_EQ(28, cl::serialize_icell(-2, 0, 0));
  EXPECT_EQ(2696, cl::serialize_icell(6, 3, 2));
  EXPECT_EQ(2697, cl::serialize_icell(6, 3, -3));
  EXPECT_EQ(2698, cl::serialize_icell(6, -4, 2));
  EXPECT_EQ(2700, cl::serialize_icell(-7, 3, 2));
}

TEST(SerializeICell, unique) {
  std::set<size_t> s;
  for (int i0 = -10; i0 < 10; ++i0) {
    for (int i1 = -10; i1 < 10; ++i1) {
      for (int i2 = -10; i2 < 10; ++i2) {
        int icell[3]{i0, i1, i2};
        size_t serial1 = cl::serialize_icell(icell);
        size_t serial2 = cl::serialize_icell(i0, i1, i2);
        EXPECT_EQ(serial1, serial2);
        s.insert(serial1);
      }
    }
  }
  EXPECT_EQ(20*20*20, s.size());
}


// sensible_threshold
// ~~~~~~~~~~~~~~~~~~

TEST(SensibleThreshold, examples) {
  for (int irep = 0; irep < NREP; ++irep) {
    // Random points
    double points[3*NPOINT];
    fill_random_double(31 + irep, points, 3*NPOINT, -5.0, 5.0);

    // Non-periodic case
    cl::Cell cell0;
    double threshold0(cl::sensible_threshold(points, NPOINT, &cell0));
    EXPECT_LT(3.0, threshold0);
    EXPECT_GT(4.0, threshold0);
    int shape0[3];
    std::unique_ptr<cl::Cell> subcell0(cell0.create_subcell(threshold0, shape0));
    EXPECT_EQ(0, shape0[0]);
    EXPECT_EQ(0, shape0[1]);
    EXPECT_EQ(0, shape0[2]);
    EXPECT_TRUE(subcell0->cubic());
    EXPECT_LT(27.0, subcell0->volume());
    EXPECT_GT(64.0, subcell0->volume());

    // 1D periodic case
    std::unique_ptr<cl::Cell> cell1(cl::create_random_cell(irep + 50, 1, 10.0, 0.3));
    double threshold1(cl::sensible_threshold(points, NPOINT, cell1.get()));
    EXPECT_LT(2.0, threshold1);
    EXPECT_GT(5.0, threshold1);

    // 2D periodic case
    std::unique_ptr<cl::Cell> cell2(cl::create_random_cell(irep + 50, 2, 2.0, 0.3));
    double threshold2(cl::sensible_threshold(points, NPOINT, cell2.get()));
    EXPECT_LT(0.0, threshold2);
    EXPECT_GT(2.0, threshold2);

    // 3D periodic case
    std::unique_ptr<cl::Cell> cell3(cl::create_random_cell(irep + 50, 3, 5.0, 0.3));
    double threshold3(cl::sensible_threshold(nullptr, NPOINT, cell3.get()));
    EXPECT_LT(0.0, threshold3);
    EXPECT_GT(2.0, threshold3);
  }
}


// BoxSortedPoints
// ~~~~~~~~~~~~~~~

TEST(BoxSortedPointsTest, exception) {
  double* points(nullptr);
  cl::Cell cell(nullptr, 0);
  EXPECT_THROW(cl::BoxSortedPoints bsp(points, 0, &cell, 0.2), std::logic_error);
}


TEST(BoxSortedPointsTest, ranges_example_0) {
  double points[9]{3.1, -1.0, -0.5, 3.0, 2.9, 0.0, 0.7, -1.1, 0.1};
  cl::Cell cell(nullptr, 0);
  cl::BoxSortedPoints bsp(points, 3, &cell, 1.0);
  EXPECT_EQ(3, bsp.npoint());
  EXPECT_EQ(3, bsp.subcell()->nvec());
  EXPECT_NEAR(1.0, bsp.subcell()->vecs()[0], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[1], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[2], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[3], EPS);
  EXPECT_NEAR(1.0, bsp.subcell()->vecs()[4], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[5], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[6], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[7], EPS);
  EXPECT_NEAR(1.0, bsp.subcell()->vecs()[8], EPS);
  EXPECT_EQ(0, bsp.shape()[0]);
  EXPECT_EQ(0, bsp.shape()[1]);
  EXPECT_EQ(0, bsp.shape()[2]);
  EXPECT_TRUE(std::equal(points, points + 9, bsp.points()));
  EXPECT_EQ(2, bsp.ipoints()[0]);
  EXPECT_EQ(0, bsp.ipoints()[1]);
  EXPECT_EQ(1, bsp.ipoints()[2]);
  // for (const auto& n : *bsp.ranges())
  //  std::cout << "Key:[" << n.first << "] Value:[" << n.second[0] << "," << n.second[1] << "]\n";
  size_t serial0(cl::serialize_icell(0, -2, 0));
  EXPECT_EQ(0, bsp.ranges()->at(serial0)[0]);
  EXPECT_EQ(1, bsp.ranges()->at(serial0)[1]);
  size_t serial1(cl::serialize_icell(3, -1, -1));
  EXPECT_EQ(1, bsp.ranges()->at(serial1)[0]);
  EXPECT_EQ(2, bsp.ranges()->at(serial1)[1]);
  size_t serial2(cl::serialize_icell(3, 2, 0));
  EXPECT_EQ(2, bsp.ranges()->at(serial2)[0]);
  EXPECT_EQ(3, bsp.ranges()->at(serial2)[1]);
}


TEST(BoxSortedPointsTest, ranges_example_3) {
  double points[9]{3.1, -1.0, -0.5, 3.0, 2.9, 0.0, 0.7, -1.1, 0.1};
  double vecs[9]{2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0};
  cl::Cell cell(vecs, 3);
  EXPECT_EQ(3, cell.nvec());
  cl::BoxSortedPoints bsp(points, 3, &cell, 1.0);
  EXPECT_EQ(3, bsp.npoint());
  EXPECT_EQ(3, bsp.subcell()->nvec());
  EXPECT_NEAR(1.0, bsp.subcell()->vecs()[0], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[1], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[2], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[3], EPS);
  EXPECT_NEAR(1.0, bsp.subcell()->vecs()[4], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[5], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[6], EPS);
  EXPECT_NEAR(0.0, bsp.subcell()->vecs()[7], EPS);
  EXPECT_NEAR(1.0, bsp.subcell()->vecs()[8], EPS);
  EXPECT_EQ(2, bsp.shape()[0]);
  EXPECT_EQ(3, bsp.shape()[1]);
  EXPECT_EQ(4, bsp.shape()[2]);
  double wrapped_points[9] = {1.1, 2.0, 3.5, 1.0, 2.9, 0.0, 0.7, 1.9, 0.1};
  EXPECT_TRUE(std::equal(wrapped_points, wrapped_points + 9, bsp.points(),
                         [](double x1, double x2){ return fabs(x1 - x2) < EPS; }));
  EXPECT_EQ(2, bsp.ipoints()[0]);
  EXPECT_EQ(1, bsp.ipoints()[1]);
  EXPECT_EQ(0, bsp.ipoints()[2]);
  // for (const auto& n : *bsp.ranges())
  //   std::cout << "Key:[" << n.first << "] Value:[" << n.second[0] << "," << n.second[1] << "]\n";
  size_t serial0(cl::serialize_icell(0, 1, 0));
  EXPECT_EQ(0, bsp.ranges()->at(serial0)[0]);
  EXPECT_EQ(1, bsp.ranges()->at(serial0)[1]);
  size_t serial1(cl::serialize_icell(1, 2, 0));
  EXPECT_EQ(1, bsp.ranges()->at(serial1)[0]);
  EXPECT_EQ(2, bsp.ranges()->at(serial1)[1]);
  size_t serial2(cl::serialize_icell(1, 2, 3));
  EXPECT_EQ(2, bsp.ranges()->at(serial2)[0]);
  EXPECT_EQ(3, bsp.ranges()->at(serial2)[1]);
}


TEST(BoxSortedPointsTest, random_0) {
  for (int irep = 0; irep < NREP; ++irep) {
    // Non-periodic system.
    cl::Cell cell(nullptr, 0);

    // Generate random points, not yet wrapped in cell.
    double points[3*NPOINT];
    fill_random_double(3157 + irep, points, 3*NPOINT, -5.0, 5.0);

    // make sorted points.
    double threshold(irep % 2 == 0 ? 0.6 : -1.0);
    cl::BoxSortedPoints bsp(points, NPOINT, &cell, threshold);

    // Shape should be three zeros
    EXPECT_EQ(0, bsp.shape()[0]);
    EXPECT_EQ(0, bsp.shape()[1]);
    EXPECT_EQ(0, bsp.shape()[2]);

    // points should not have changed.
    EXPECT_TRUE(std::equal(points, points + 3*NPOINT, bsp.points()));
  }
}


TEST(BoxSortedPointsTest, random_1) {
  for (int irep = 0; irep < NREP; ++irep) {
    // Get a random 1D cell
    std::unique_ptr<cl::Cell> cell(cl::create_random_cell(irep*NPOINT, 1, 2.0));

    // Generate random points, not yet wrapped in cell.
    double points[3*NPOINT];
    fill_random_double(31 + irep, points, 3*NPOINT, -5.0, 5.0);

    // make sorted points.
    double threshold(irep % 2 == 0 ? 0.6 : -1.0);
    cl::BoxSortedPoints bsp(points, NPOINT, cell.get(), threshold);

    // Check shape
    EXPECT_LT(0, bsp.shape()[0]);
    EXPECT_EQ(0, bsp.shape()[1]);
    EXPECT_EQ(0, bsp.shape()[2]);

    // Check fractional coordinates of all  points
    for (size_t ipoint=0; ipoint < NPOINT; ++ipoint) {
      double frac0[3];
      bsp.subcell()->to_frac(bsp.points() + 3*ipoint, frac0);
      double frac1[3];
      cell->to_frac(bsp.points() + 3*ipoint, frac1);
      double frac2[3];
      cell->to_frac(points + 3*ipoint, frac2);
      EXPECT_LE(0.0, frac0[0]);
      EXPECT_GT(bsp.shape()[0], frac0[0]);
      EXPECT_LE(0.0, frac1[0]);
      EXPECT_GT(1.0, frac1[0]);
      double diff(frac1[0] - frac2[0]);
      EXPECT_NEAR(diff, round(diff), EPS);
      for (int ivec=1; ivec < 3; ++ivec) {
        EXPECT_NEAR(frac1[ivec], frac2[ivec], EPS);
      }
    }
  }
}


TEST(BoxSortedPointsTest, random_2) {
  for (int irep = 0; irep < NREP; ++irep) {
    // Get a random 2D cell
    std::unique_ptr<cl::Cell> cell(cl::create_random_cell(irep*NPOINT, 2, 2.0));

    // Generate random points, not yet wrapped in cell.
    double points[3*NPOINT];
    fill_random_double(31 + irep, points, 3*NPOINT, -5.0, 5.0);

    // make sorted points.
    double threshold(irep % 2 == 0 ? 0.6 : -1.0);
    cl::BoxSortedPoints bsp(points, NPOINT, cell.get(), threshold);

    // Check shape
    EXPECT_LT(0, bsp.shape()[0]);
    EXPECT_LT(0, bsp.shape()[1]);
    EXPECT_EQ(0, bsp.shape()[2]);

    // Check fractional coordinates of all  points
    for (size_t ipoint=0; ipoint < NPOINT; ++ipoint) {
      double frac0[3];
      bsp.subcell()->to_frac(bsp.points() + 3*ipoint, frac0);
      double frac1[3];
      cell->to_frac(bsp.points() + 3*ipoint, frac1);
      double frac2[3];
      cell->to_frac(points + 3*ipoint, frac2);
      for (int ivec=0; ivec < 2; ++ivec) {
        EXPECT_LE(0.0, frac0[ivec]);
        EXPECT_GT(bsp.shape()[ivec], frac0[ivec]);
        EXPECT_LE(0.0, frac1[ivec]);
        EXPECT_GT(1.0, frac1[ivec]);
        double diff(frac1[0] - frac2[0]);
        EXPECT_NEAR(diff, round(diff), EPS);
      }
      EXPECT_NEAR(frac1[2], frac2[2], EPS);
    }
  }
}


TEST(BoxSortedPointsTest, random_3) {
  for (int irep = 0; irep < NREP; ++irep) {
    // Get a random 3D cell
    std::unique_ptr<cl::Cell> cell(cl::create_random_cell(irep*NPOINT, 3, 2.0));

    // Generate random points, not yet wrapped in cell
    double points[3*NPOINT];
    fill_random_double(667 + irep, points, 3*NPOINT, -5.0, 5.0);

    // make sorted points
    double threshold(irep % 2 == 0 ? 0.2 : -1.0);
    cl::BoxSortedPoints bsp(points, NPOINT, cell.get(), threshold);

    // Check shape
    EXPECT_LT(0, bsp.shape()[0]);
    EXPECT_LT(0, bsp.shape()[1]);
    EXPECT_LT(0, bsp.shape()[2]);

    // Check subcell
    for (int ivec = 0; ivec < 3; ++ivec)
      EXPECT_NEAR(cell->spacings()[ivec],
                  bsp.subcell()->spacings()[ivec]*bsp.shape()[ivec],
                  EPS);

    // Check fractional coordinates of all  points
    for (size_t ipoint=0; ipoint < NPOINT; ++ipoint) {
      double frac0[3];
      bsp.subcell()->to_frac(bsp.points() + 3*ipoint, frac0);
      double frac1[3];
      cell->to_frac(bsp.points() + 3*ipoint, frac1);
      double frac2[3];
      cell->to_frac(points + 3*ipoint, frac2);
      for (int ivec=0; ivec < 3; ++ivec) {
        EXPECT_LE(0.0, frac0[ivec]);
        EXPECT_GT(bsp.shape()[ivec], frac0[ivec]);
        EXPECT_LE(0.0, frac1[ivec]);
        EXPECT_GT(1.0, frac1[ivec]);
        double diff(frac1[ivec] - frac2[ivec]);
        EXPECT_NEAR(diff, round(diff), EPS);
      }
    }

    // Test if all serials are from the list of allowed values
    std::vector<size_t> allowed_serials;
    for (int icell0 = 0; icell0 < bsp.shape()[0]; ++icell0) {
      for (int icell1 = 0; icell1 < bsp.shape()[1]; ++icell1) {
        for (int icell2 = 0; icell2 < bsp.shape()[2]; ++icell2) {
          allowed_serials.push_back(cl::serialize_icell(icell0, icell1, icell2));
        }
      }
    }
    for (const auto& iter : *bsp.ranges())
      EXPECT_NE(std::find(allowed_serials.begin(), allowed_serials.end(), iter.first),
                allowed_serials.end());
  }
}


TEST(BoxSortedPointsTest, ranges_example2_3) {
  double vecs[9]{3.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 3.0};
  cl::Cell cell(vecs, 3);
  for (int ivec = 0; ivec < 3; ++ivec) {
    double points[9]{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    points[ivec] = 2.5;
    cl::BoxSortedPoints bsp(points, 3, &cell, 1.0);
    /*
    for (const auto& n : *bsp.ranges()) {
      std::cout << "Key:[" << n.first << "] Value:[" << n.second[0] << "," << n.second[1] << "]\n";
      for (size_t ipoint = n.second[0]; ipoint < n.second[1]; ++ipoint) {
        std::cout << ipoint << " " << bsp.ipoints()[ipoint] << std::endl;
        for (int i = 0; i < 3; i++) {
          std::cout << "  " << points[3*bsp.ipoints()[ipoint] + i];
          std::cout << "  " << bsp.points()[3*bsp.ipoints()[ipoint] + i] << std::endl;
        }
      }
    }
    */
    EXPECT_EQ(2, bsp.ranges()->size());
    // icell0
    std::array<size_t, 2> range0(bsp.ranges()->at(cl::serialize_icell(0, 0, 0)));
    EXPECT_EQ(0, range0[0]);
    EXPECT_EQ(2, range0[1]);
    // icell1
    int icell[3]{0, 0, 0};
    icell[ivec] = 2;
    std::array<size_t, 2> range1(bsp.ranges()->at(cl::serialize_icell(icell)));
    EXPECT_EQ(2, range1[0]);
    EXPECT_EQ(3, range1[1]);
  }
}


// BoxCutoffIterator
// ~~~~~~~~~~~~~~~~~

TEST(BoxCutoffIteratorTest, exception_radius) {
  double vecs[9]{2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell cell(vecs, 3);
  const double points[6]{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  cl::BoxSortedPoints bsp(points, 2, &cell);
  const double center[3]{0.0, 0.0, 1.0};
  EXPECT_THROW(new cl::BoxCutoffIterator(&bsp, center, -1.0), std::domain_error);
}


TEST(BoxCutoffIteratorTest, exception_increment) {
  double vecs[9]{2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell cell(vecs, 3);
  const double points[6]{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  cl::BoxSortedPoints bsp(points, 2, &cell, 1.0);
  const double center[3]{0.0, 0.0, 1.0};
  cl::BoxCutoffIterator bci(&bsp, center, 1e-15);
  EXPECT_TRUE(bci.busy());
  EXPECT_THROW(bci++, std::logic_error);
}


TEST(BoxCutoffIteratorTest, example1) {
  // Problem definition: points, cell, center and cutoff.
  // 1) points
  double points[12]{9.0, 9.0, 99.0, 22.5, 0.0, 0.0, 4.0, 5.1, -3.0, 4.0, 5.2, -3.0};
  // 2) cell
  double vecs[9]{20.0, 0.0, 0.0, 0.0, 20.0, 0.0, 0.0, 0.0, 20.0};
  cl::Cell cell(vecs, 3);
  // 3) center and cutoff
  double center[3]{1.0, 1.0, 1.0};
  double cutoff = 10.0;

  // Construct the iterator.
  cl::BoxSortedPoints bsp(points, 4, &cell, 2.0);
  cl::BoxCutoffIterator bci(&bsp, center, cutoff);
  EXPECT_TRUE(bci.busy());

  // Iteration
  EXPECT_EQ(1, bci.ipoint());
  EXPECT_NEAR(1.5, bci.delta()[0], EPS);
  EXPECT_NEAR(-1.0, bci.delta()[1], EPS);
  EXPECT_NEAR(-1.0, bci.delta()[2], EPS);
  EXPECT_NEAR(sqrt(1.5*1.5 + 1.0 + 1.0), bci.distance(), EPS);
  ++bci;
  EXPECT_TRUE(bci.busy());
  // Iteration
  EXPECT_EQ(2, bci.ipoint());
  EXPECT_NEAR(3.0, bci.delta()[0], EPS);
  EXPECT_NEAR(4.1, bci.delta()[1], EPS);
  EXPECT_NEAR(-4.0, bci.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0*3.0 + 4.1*4.1 + 4.0*4.0), bci.distance(), EPS);
  ++bci;
  EXPECT_TRUE(bci.busy());
  // Iteration
  EXPECT_EQ(3, bci.ipoint());
  EXPECT_NEAR(3.0, bci.delta()[0], EPS);
  EXPECT_NEAR(4.2, bci.delta()[1], EPS);
  EXPECT_NEAR(-4.0, bci.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0*3.0 + 4.2*4.2 + 4.0*4.0), bci.distance(), EPS);
  ++bci;
  EXPECT_FALSE(bci.busy());
}


TEST(BoxCutoffIteratorTest, example2) {
  // Problem definition: points, cell, center and cutoff.
  // 1) points
  double points[12]{0.0, 0.0, 0.0, 22.5, 0.0, 0.0, 4.0, 5.1, -3.0, 4.0, 5.2, -3.0};
  // 2) cell
  cl::Cell cell(nullptr, 0);
  // 3) center and cutoff
  double center[3]{1.0, 1.0, 1.0};
  double cutoff = 10.0;

  // Construct the iterator.
  cl::BoxSortedPoints bsp(points, 4, &cell, 1.0);
  cl::BoxCutoffIterator bci(&bsp, center, cutoff);
  EXPECT_TRUE(bci.busy());

  // First iteration
  EXPECT_NEAR(-1.0, bci.delta()[0], EPS);
  EXPECT_NEAR(-1.0, bci.delta()[1], EPS);
  EXPECT_NEAR(-1.0, bci.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0), bci.distance(), EPS);
  EXPECT_EQ(0, bci.ipoint());
  ++bci;
  EXPECT_TRUE(bci.busy());
  // Second iteration
  EXPECT_NEAR(3.0, bci.delta()[0], EPS);
  EXPECT_NEAR(4.1, bci.delta()[1], EPS);
  EXPECT_NEAR(-4.0, bci.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0*3.0 + 4.1*4.1 + 4.0*4.0), bci.distance(), EPS);
  EXPECT_EQ(2, bci.ipoint());
  ++bci;
  EXPECT_TRUE(bci.busy());
  // Third iteration
  EXPECT_NEAR(3.0, bci.delta()[0], EPS);
  EXPECT_NEAR(4.2, bci.delta()[1], EPS);
  EXPECT_NEAR(-4.0, bci.delta()[2], EPS);
  EXPECT_NEAR(sqrt(3.0*3.0 + 4.2*4.2 + 4.0*4.0), bci.distance(), EPS);
  EXPECT_EQ(3, bci.ipoint());
  ++bci;
  EXPECT_FALSE(bci.busy());
}


// Instantiation of parameterized tests
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

INSTANTIATE_TEST_CASE_P(RangesCutoffTest0123, RangesCutoffTestP, ::testing::Range(0, 4));
INSTANTIATE_TEST_CASE_P(BarsCutoffTest0123, BarsCutoffTestP, ::testing::Range(1, 4));
INSTANTIATE_TEST_CASE_P(BarIteratorTest0123, BarIteratorTestP, ::testing::Range(0, 4));


// vim: textwidth=90 et ts=2 sw=2
