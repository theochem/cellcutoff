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


#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

#include <cellcutoff/decomposition.h>

#include "common.h"


namespace cl = cellcutoff;


// Point class
// ~~~~~~~~~~~

TEST(PointTest, constructor1) {
  double cart[3];
  fill_random_double(1475, cart, 3);
  cl::Point point(cart);
  // Make sure the contents is copied, not the pointer
  EXPECT_NE(cart, point.cart_);
  EXPECT_EQ(cart[0], point.cart_[0]);
  EXPECT_EQ(cart[1], point.cart_[1]);
  EXPECT_EQ(cart[2], point.cart_[2]);
  // When not given, icell must be zero
  EXPECT_EQ(0, point.icell_[0]);
  EXPECT_EQ(0, point.icell_[1]);
  EXPECT_EQ(0, point.icell_[2]);
}


TEST(PointTest, constructor2) {
  double cart[3];
  int icell[3];
  fill_random_double(1879, cart, 3);
  fill_random_int(5849, icell, 3, 0, 10);
  cl::Point point(cart, icell);
  // Make sure the contents is copied, not the pointer
  EXPECT_NE(cart, point.cart_);
  EXPECT_EQ(cart[0], point.cart_[0]);
  EXPECT_EQ(cart[1], point.cart_[1]);
  EXPECT_EQ(cart[2], point.cart_[2]);
  EXPECT_NE(icell, point.icell_);
  EXPECT_EQ(icell[0], point.icell_[0]);
  EXPECT_EQ(icell[1], point.icell_[1]);
  EXPECT_EQ(icell[2], point.icell_[2]);
}


TEST(PointTest, less_than) {
  double cart[3]{1.0, 2.0, 3.0};
  int icell[3]{0, 2, 3};
  cl::Point a(cart, icell);
  icell[0] = -1;
  cl::Point b(cart, icell);
  EXPECT_LT(b, a);
  icell[0] = 1;
  cl::Point c(cart, icell);
  EXPECT_LT(a, c);
  icell[0] = 0;
  icell[1] = 1;
  cl::Point d(cart, icell);
  EXPECT_LT(d, a);
  icell[1] = 3;
  cl::Point e(cart, icell);
  EXPECT_LT(a, e);
  icell[1] = 2;
  icell[2] = 2;
  cl::Point f(cart, icell);
  EXPECT_LT(f, a);
  icell[2] = 4;
  cl::Point g(cart, icell);
  EXPECT_LT(a, g);
  EXPECT_FALSE(g < a);
}


// Decomposition
// ~~~~~~~~~~~~~

TEST(DecompositionTest, assign_icell_domain) {
  double vecs[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell subcell0(vecs, 0);
  EXPECT_THROW(cl::assign_icell(subcell0, nullptr, 0, 0), std::domain_error);
  EXPECT_THROW(cl::assign_icell(subcell0, nullptr, nullptr, 0, 0), std::domain_error);
  cl::Cell subcell1(vecs, 1);
  EXPECT_THROW(cl::assign_icell(subcell1, nullptr, 0, 0), std::domain_error);
  EXPECT_THROW(cl::assign_icell(subcell1, nullptr, nullptr, 0, 0), std::domain_error);
  cl::Cell subcell2(vecs, 2);
  EXPECT_THROW(cl::assign_icell(subcell2, nullptr, 0, 0), std::domain_error);
  EXPECT_THROW(cl::assign_icell(subcell2, nullptr, nullptr, 0, 0), std::domain_error);
}


TEST(DecompositionTest, assign_icell_example) {
  std::vector<cl::Point> points;
  double cart0[3]{3.1, -1.0, -0.5};
  double cart1[3]{3.0, 2.9, 0.0};
  double cart2[3]{0.7, -1.1, 0.1};
  points.push_back(cl::Point(cart0));
  points.push_back(cl::Point(cart1));
  points.push_back(cl::Point(cart2));
  double vecs[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell subcell(vecs, 3);
  cl::assign_icell(subcell, points.data(), points.size(), sizeof(cl::Point));
  EXPECT_EQ(3, points[0].icell_[0]);
  EXPECT_EQ(-1, points[0].icell_[1]);
  EXPECT_EQ(-1, points[0].icell_[2]);
  EXPECT_EQ(3, points[1].icell_[0]);
  EXPECT_EQ(2, points[1].icell_[1]);
  EXPECT_EQ(0, points[1].icell_[2]);
  EXPECT_EQ(0, points[2].icell_[0]);
  EXPECT_EQ(-2, points[2].icell_[1]);
  EXPECT_EQ(0, points[2].icell_[2]);
}


TEST(DecompositionTest, assign_icell_example_shape) {
  std::vector<cl::Point> points;
  double cart0[3]{3.1, -1.0, -0.5};
  double cart1[3]{2.0, 2.9, 0.0};
  double cart2[3]{0.7, -1.1, 0.1};
  points.push_back(cl::Point(cart0));
  points.push_back(cl::Point(cart1));
  points.push_back(cl::Point(cart2));
  double vecs[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell subcell(vecs, 3);
  const int shape[3]{2, 3, 4};
  cl::assign_icell(subcell, shape, points.data(), points.size(), sizeof(cl::Point));
  EXPECT_EQ(1, points[0].icell_[0]);
  EXPECT_EQ(2, points[0].icell_[1]);
  EXPECT_EQ(3, points[0].icell_[2]);
  EXPECT_EQ(0, points[1].icell_[0]);
  EXPECT_EQ(2, points[1].icell_[1]);
  EXPECT_EQ(0, points[1].icell_[2]);
  EXPECT_EQ(0, points[2].icell_[0]);
  EXPECT_EQ(1, points[2].icell_[1]);
  EXPECT_EQ(0, points[2].icell_[2]);
  EXPECT_NEAR(1.1, points[0].cart_[0], EPS);
  EXPECT_NEAR(2.0, points[0].cart_[1], EPS);
  EXPECT_NEAR(3.5, points[0].cart_[2], EPS);
  EXPECT_NEAR(0.0, points[1].cart_[0], EPS);
  EXPECT_NEAR(2.9, points[1].cart_[1], EPS);
  EXPECT_NEAR(0.0, points[1].cart_[2], EPS);
  EXPECT_NEAR(0.7, points[2].cart_[0], EPS);
  EXPECT_NEAR(1.9, points[2].cart_[1], EPS);
  EXPECT_NEAR(0.1, points[2].cart_[2], EPS);
}


TEST(DecompositionTest, assign_icell_random_wrap) {
  for (int irep = 0; irep < NREP; ++irep) {
    // Get a random 3D cell
    std::unique_ptr<cl::Cell> cell(create_random_cell_nvec(irep*NPOINT, 3, 2));
    // Get a subcell
    double threshold = 0.2;
    int shape[3] = {-1, -1, -1};
    std::unique_ptr<cl::Cell> subcell(cell->create_subcell(threshold, shape));
    for (int ivec = 0; ivec < 3; ++ivec)
      EXPECT_NEAR(cell->spacings()[0], subcell->spacings()[0]*shape[0], EPS);

    // Generate random points, not yet wrapped in cell
    std::vector<cl::Point> points;
    for (int ipoint = 0; ipoint < NPOINT; ++ipoint) {
      double cart[3];
      fill_random_double(ipoint+3157, cart, 3, -5.0, 5.0);
      points.push_back(cl::Point(cart));
    }

    // Actual calculation of interest
    cl::assign_icell(*subcell, shape, points.data(), points.size(), sizeof(cl::Point));

    // Check all point and icell fields
    for (const auto& point : points) {
      double frac[3];
      subcell->to_frac(point.cart_, frac);
      EXPECT_LE(0.0, frac[0]);
      EXPECT_LE(0.0, frac[1]);
      EXPECT_LE(0.0, frac[2]);
      EXPECT_GT(shape[0], frac[0]);
      EXPECT_GT(shape[1], frac[1]);
      EXPECT_GT(shape[2], frac[2]);
      EXPECT_LE(0, point.icell_[0]);
      EXPECT_GT(shape[0], point.icell_[0]);
      EXPECT_LE(0, point.icell_[1]);
      EXPECT_GT(shape[1], point.icell_[1]);
      EXPECT_LE(0, point.icell_[2]);
      EXPECT_GT(shape[2], point.icell_[2]);
    }
  }
}


TEST(DecompositionTest, cell_map_example) {
  double vecs[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell cell(vecs, 3);
  for (int ivec = 0; ivec < 3; ++ivec) {
    double cart0[3] = {0.0, 0.0, 0.0};
    double cart1[3] = {0.0, 0.0, 0.0};
    double cart2[3] = {0.0, 0.0, 0.0};
    cart0[ivec] = 2.5;
    std::vector<cl::Point> points;
    points.push_back(cl::Point(cart0));
    points.push_back(cl::Point(cart1));
    points.push_back(cl::Point(cart2));
    cl::assign_icell(cell, points.data(), points.size(), sizeof(cl::Point));
    std::unique_ptr<cl::CellMap> cell_map(cl::create_cell_map(points.data(), points.size(), sizeof(cl::Point)));
    EXPECT_EQ(2, cell_map->size());
    // icell0
    std::array<int, 3> icell0{0, 0, 0};
    std::array<size_t, 2> range0(cell_map->at(icell0));
    EXPECT_EQ(1, range0[0]);
    EXPECT_EQ(3, range0[1]);
    // icell1
    std::array<int, 3> icell1{2*(ivec == 0), 2*(ivec == 1), 2*(ivec == 2)};
    std::array<size_t, 2> range1(cell_map->at(icell1));
    EXPECT_EQ(0, range1[0]);
    EXPECT_EQ(1, range1[1]);
  }
}


TEST(DecompositionTest, cell_map_points_not_ordered) {
  double vecs[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell cell(vecs, 3);
  for (int ivec = 0; ivec < 3; ++ivec) {
    double cart0[3] = {0.0, 0.0, 0.0};
    double cart1[3] = {0.0, 0.0, 0.0};
    double cart2[3] = {0.0, 0.0, 0.0};
    cart1[ivec] = 2.0;
    cart2[ivec] = 0.1;
    std::vector<cl::Point> points;
    points.push_back(cl::Point(cart0));
    points.push_back(cl::Point(cart1));
    points.push_back(cl::Point(cart2));
    cl::assign_icell(cell, points.data(), points.size(), sizeof(cl::Point));
    EXPECT_THROW(cl::create_cell_map(points.data(), points.size(), sizeof(cl::Point)), cl::points_not_grouped);
  }
}


TEST(DecompositionTest, random_cell_map) {
  for (int irep = 0; irep < NREP; ++irep) {
    // Apply entire machinery on random data
    std::vector<cl::Point> points;
    for (int ipoint = 0; ipoint < NPOINT; ++ipoint) {
      double cart[3];
      fill_random_double(ipoint+3157, cart, 3, -5.0, 5.0);
      points.push_back(cl::Point(cart));
    }
    std::unique_ptr<cl::Cell> subcell(create_random_cell_nvec(irep*NPOINT, 3));
    cl::assign_icell(*subcell, points.data(), points.size(), sizeof(cl::Point));
    cl::sort_by_icell(points.data(), points.size(), sizeof(cl::Point));
    std::unique_ptr<cl::CellMap> cell_map(cl::create_cell_map(points.data(), points.size(), sizeof(cl::Point)));
    // Check consistency of results: loop over map
    for (const auto& kv : *cell_map) {
      const size_t begin = kv.second[0];
      const size_t end = kv.second[1];
      for (size_t ipoint = begin; ipoint < end; ++ipoint) {
        EXPECT_EQ(kv.first[0], points.at(ipoint).icell_[0]);
        EXPECT_EQ(kv.first[1], points.at(ipoint).icell_[1]);
        EXPECT_EQ(kv.first[2], points.at(ipoint).icell_[2]);
      }
    }
    // Check consistency of results: loop over points
    for (size_t ipoint = 0; ipoint < points.size(); ++ipoint) {
      const int* icell_ = points.at(ipoint).icell_;
      const std::array<int, 3> icell{icell_[0], icell_[1], icell_[2]};
      EXPECT_GE(ipoint, cell_map->at(icell)[0]);
      EXPECT_LT(ipoint, cell_map->at(icell)[1]);
    }
  }
}


// robust_wrap
// ~~~~~~~~~~~

TEST(SmartWrap, examples) {
  EXPECT_EQ(0, cl::robust_wrap(-15, 5));
  EXPECT_EQ(0, cl::robust_wrap(-5, 5));
  EXPECT_EQ(2, cl::robust_wrap(-3, 5));
  EXPECT_EQ(4, cl::robust_wrap(-1, 5));
  EXPECT_EQ(0, cl::robust_wrap(0, 5));
  EXPECT_EQ(3, cl::robust_wrap(3, 5));
  EXPECT_EQ(0, cl::robust_wrap(5, 5));
  EXPECT_EQ(1, cl::robust_wrap(6, 5));
  EXPECT_EQ(0, cl::robust_wrap(10, 5));
  EXPECT_EQ(2, cl::robust_wrap(12, 5));
}


TEST(SmartWrap, examples_division) {
  int base = 5;
  for (int i = -20; i < 20; ++i) {
    int div = 10000;
    int m1 = cl::robust_wrap(i, base);
    int m2 = cl::robust_wrap(i, base, &div);
    EXPECT_EQ(m1, m2);
    EXPECT_EQ(i, div*base + m1);
  }
}


// icell_hash
// ~~~~~~~~~~

TEST(ICellHash, examples) {
  cl::icell_hash fn;
  EXPECT_EQ(0, fn(std::array<int, 3>{0, 0, 0}));
  EXPECT_EQ(1, fn(std::array<int, 3>{0, 0, -1}));
  EXPECT_EQ(2, fn(std::array<int, 3>{0, -1, 0}));
  EXPECT_EQ(3, fn(std::array<int, 3>{0, -1, -1}));
  EXPECT_EQ(4, fn(std::array<int, 3>{-1, 0, 0}));
  EXPECT_EQ(5, fn(std::array<int, 3>{-1, 0, -1}));
  EXPECT_EQ(6, fn(std::array<int, 3>{-1, -1, 0}));
  EXPECT_EQ(7, fn(std::array<int, 3>{-1, -1, -1}));
  EXPECT_EQ(8, fn(std::array<int, 3>{0, 0, 1}));
  EXPECT_EQ(16, fn(std::array<int, 3>{0, 1, 0}));
  EXPECT_EQ(24, fn(std::array<int, 3>{1, 0, 0}));
  EXPECT_EQ(25, fn(std::array<int, 3>{1, 0, -1}));
  EXPECT_EQ(26, fn(std::array<int, 3>{1, -1, 0}));
  EXPECT_EQ(28, fn(std::array<int, 3>{-2, 0, 0}));
  EXPECT_EQ(2696, fn(std::array<int, 3>{6, 3, 2}));
  EXPECT_EQ(2697, fn(std::array<int, 3>{6, 3, -3}));
  EXPECT_EQ(2698, fn(std::array<int, 3>{6, -4, 2}));
  EXPECT_EQ(2700, fn(std::array<int, 3>{-7, 3, 2}));
}

TEST(ICellHash, unique) {
  std::set<size_t> s;
  cl::icell_hash fn;
  for (int i0 = -10; i0 < 10; ++i0) {
    for (int i1 = -10; i1 < 10; ++i1) {
      for (int i2 = -10; i2 < 10; ++i2) {
        size_t key = fn(std::array<int, 3>{i0, i1, i2});
        s.insert(key);
      }
    }
  }
  EXPECT_EQ(20*20*20, s.size());
}

// vim: textwidth=90 et ts=2 sw=2
