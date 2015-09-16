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
// aint with this program; if not, see <http://www.gnu.org/licenses/>
//
//--

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
    // Make sure the contents is copied, not the pointer
    EXPECT_NE(cart, point.cart);
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
    cl::Point point(5, cart, icell);
    // Make sure the contents is copied, not the pointer
    EXPECT_NE(cart, point.cart);
    EXPECT_EQ(cart[0], point.cart[0]);
    EXPECT_EQ(cart[1], point.cart[1]);
    EXPECT_EQ(cart[2], point.cart[2]);
    EXPECT_NE(icell, point.icell);
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
