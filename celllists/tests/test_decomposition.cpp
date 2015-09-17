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


#include <algorithm>

#include <gtest/gtest.h>

#include <celllists/decomposition.h>

#include "common.h"


namespace cl = celllists;


// Point class
// ~~~~~~~~~~~

TEST(PointTest, constructor1) {
  double cart[3];
  fill_random_double(1475, cart, 3);
  cl::Point point(5, cart);
  EXPECT_EQ(5, point.index);
  // Make sure the contents is copied, not the pointer
  EXPECT_NE(cart, point.cart.data());
  EXPECT_EQ(cart[0], point.cart[0]);
  EXPECT_EQ(cart[1], point.cart[1]);
  EXPECT_EQ(cart[2], point.cart[2]);
  // When not given, icell must be zero
  EXPECT_EQ(0, point.icell[0]);
  EXPECT_EQ(0, point.icell[1]);
  EXPECT_EQ(0, point.icell[2]);
}


TEST(PointTest, constructor2) {
  double cart[3];
  int icell[3];
  fill_random_double(1879, cart, 3);
  fill_random_int(5849, icell, 3, 0, 10);
  cl::Point point(4, cart, icell);
  EXPECT_EQ(4, point.index);
  // Make sure the contents is copied, not the pointer
  EXPECT_NE(cart, point.cart.data());
  EXPECT_EQ(cart[0], point.cart[0]);
  EXPECT_EQ(cart[1], point.cart[1]);
  EXPECT_EQ(cart[2], point.cart[2]);
  EXPECT_NE(icell, point.icell.data());
  EXPECT_EQ(icell[0], point.icell[0]);
  EXPECT_EQ(icell[1], point.icell[1]);
  EXPECT_EQ(icell[2], point.icell[2]);
}


TEST(PointTest, less_than) {
  double cart[3]{1.0, 2.0, 3.0};
  int icell[3]{0, 2, 3};
  cl::Point a(5, cart, icell);
  icell[0] = -1;
  cl::Point b(5, cart, icell);
  EXPECT_LT(b, a);
  icell[0] = 1;
  cl::Point c(5, cart, icell);
  EXPECT_LT(a, c);
  icell[0] = 0;
  icell[1] = 1;
  cl::Point d(5, cart, icell);
  EXPECT_LT(d, a);
  icell[1] = 3;
  cl::Point e(5, cart, icell);
  EXPECT_LT(a, e);
  icell[1] = 2;
  icell[2] = 2;
  cl::Point f(5, cart, icell);
  EXPECT_LT(f, a);
  icell[2] = 4;
  cl::Point g(5, cart, icell);
  EXPECT_LT(a, g);
  EXPECT_FALSE(g < a);
}


// Decomposition
// ~~~~~~~~~~~~~

TEST(DecompositionTest, assign_icell_domain) {
  double vecs[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell subcell0(vecs, 0);
  EXPECT_THROW(cl::assign_icell(subcell0, nullptr), std::domain_error);
  EXPECT_THROW(cl::assign_icell(subcell0, nullptr, nullptr), std::domain_error);
  cl::Cell subcell1(vecs, 1);
  EXPECT_THROW(cl::assign_icell(subcell1, nullptr), std::domain_error);
  EXPECT_THROW(cl::assign_icell(subcell1, nullptr, nullptr), std::domain_error);
  cl::Cell subcell2(vecs, 2);
  EXPECT_THROW(cl::assign_icell(subcell2, nullptr), std::domain_error);
  EXPECT_THROW(cl::assign_icell(subcell2, nullptr, nullptr), std::domain_error);
}


TEST(DecompositionTest, assign_icell_example) {
  std::vector<cl::Point> points;
  double cart0[3]{3.1, -1.0, -0.5};
  double cart1[3]{3.0, 2.9, 0.0};
  double cart2[3]{0.7, -1.1, 0.1};
  points.push_back(cl::Point(0, cart0));
  points.push_back(cl::Point(1, cart1));
  points.push_back(cl::Point(2, cart2));
  double vecs[9]{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  cl::Cell subcell(vecs, 3);
  cl::assign_icell(subcell, &points);
  EXPECT_EQ(3, points[0].icell[0]);
  EXPECT_EQ(-1, points[0].icell[1]);
  EXPECT_EQ(-1, points[0].icell[2]);
  EXPECT_EQ(3, points[1].icell[0]);
  EXPECT_EQ(2, points[1].icell[1]);
  EXPECT_EQ(0, points[1].icell[2]);
  EXPECT_EQ(0, points[2].icell[0]);
  EXPECT_EQ(-2, points[2].icell[1]);
  EXPECT_EQ(0, points[2].icell[2]);
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
      EXPECT_NEAR(cell->spacings()[0], subcell->spacings()[0]*shape[0], 1e-10);

    // Generate random points, wrapped in cell
    std::vector<cl::Point> points;
    for (int ipoint = 0; ipoint < NPOINT; ++ipoint) {
      double cart[3];
      fill_random_double(ipoint+3157, cart, 3, -5.0, 5.0);
      points.push_back(cl::Point(ipoint, cart));
    }

    // Actual calculation of interest
    cl::assign_icell(*subcell, &points, shape);

    // Check all icell fields, should be in range defined by shape
    for (const auto& point : points) {
      EXPECT_LE(0, point.icell[0]);
      EXPECT_GT(shape[0], point.icell[0]);
      EXPECT_LE(0, point.icell[1]);
      EXPECT_GT(shape[1], point.icell[1]);
      EXPECT_LE(0, point.icell[2]);
      EXPECT_GT(shape[2], point.icell[2]);
    }
  }
}


TEST(DecompositionTest, random_cell_map) {
  for (int irep = 0; irep < NREP; ++irep) {
    // Apply entire machinery on random data
    std::vector<cl::Point> points;
    for (int ipoint = 0; ipoint < NPOINT; ++ipoint) {
      double cart[3];
      fill_random_double(ipoint+3157, cart, 3, -5.0, 5.0);
      points.push_back(cl::Point(ipoint, cart));
    }
    std::unique_ptr<cl::Cell> subcell(create_random_cell_nvec(irep*NPOINT, 3));
    cl::assign_icell(*subcell, &points);
    std::sort(points.begin(), points.end());
    std::unique_ptr<cl::CellMap> cell_map(cl::create_cell_map(points));
    // Check consistency of results: loop over map
    for (const auto& kv : *cell_map) {
      const int begin = kv.second[0];
      const int end = kv.second[1];
      for (int ipoint = begin; ipoint < end; ++ipoint) {
        EXPECT_EQ(kv.first[0], points.at(ipoint).icell[0]);
        EXPECT_EQ(kv.first[1], points.at(ipoint).icell[1]);
        EXPECT_EQ(kv.first[2], points.at(ipoint).icell[2]);
      }
    }
    // Check consistency of results: loop over points
    int ipoint = 0;
    for (const auto& point : points) {
      EXPECT_GE(ipoint, cell_map->at(point.icell)[0]);
      EXPECT_LT(ipoint, cell_map->at(point.icell)[1]);
      ++ipoint;
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

// vim: textwidth=90 et ts=2 sw=2
