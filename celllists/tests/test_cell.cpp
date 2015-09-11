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


#include <stdexcept>
#include <gtest/gtest.h>
#include <cstdlib>
#include <cmath>
#include "celllists/cell.h"

// Helper functions
// ================

//! Fills an array of doubles with random numbers in range ]-0.5*scale, 0.5*scale]
void fill_random_double(double* array, size_t size, unsigned int seed, double scale=1) {
    srand(seed);
    for (int i=0; i<size; i++)
        array[i] = (rand() + 1.0)/(RAND_MAX + 1.0) - 0.5;
}

//! Fills an array of int with random numbers in range [-range, range]
void fill_random_int(int* array, size_t size, unsigned int seed, int range) {
    srand(seed);
    for (int i=0; i<size; i++)
        array[i] = (rand() % (2*range+1))-range;
}

Cell create_random_cell(int nvec, unsigned int seed, double scale=1) {
    double rvecs[nvec*3];
    while (true) {
        try {
            fill_random_double(rvecs, nvec*3, seed, scale);
            return Cell(rvecs, nvec);
        }
        catch (singular_cell_vectors) {}
    }
}


// Fixtures
// ========

class CellTest1 : public ::testing::Test {
    public:
        Cell* mycell;
        virtual void SetUp() {
            double rvecs[3] = {2, 0, 0};
            mycell = new Cell(rvecs, 1);
        }
        virtual void TearDown() {
            delete mycell;
        }
};


class CellTest2 : public ::testing::Test {
    public:
        Cell* mycell;
        virtual void SetUp() {
            double rvecs[6] = {2, 0, 0, 0, 0, 4};
            mycell = new Cell(rvecs, 2);
        }
        virtual void TearDown() {
            delete mycell;
        }
};


class CellTest3 : public ::testing::Test {
    public:
        Cell* mycell;
        virtual void SetUp() {
            double rvecs[9] = {2, 0, 0, 0, 1, 0, 0, 0, 4};
            mycell = new Cell(rvecs, 3);
        }
        virtual void TearDown() {
            delete mycell;
        }
};


// Tests grouped by main method being tested
// =========================================

// Constructor
// -----------

TEST(CellTest, constructor_singular1) {
    double rvecs[3] = {0, 0, 0};
    EXPECT_THROW(Cell cell(rvecs, 1), singular_cell_vectors);
}

TEST(CellTest, constructor_singular2) {
    double rvecs[6] = {1, 0, 0, 1, 0, 0};
    EXPECT_THROW(Cell cell(rvecs, 2), singular_cell_vectors);
}

TEST(CellTest, constructor_singular3) {
    double rvecs[9] = {1, 0, 0, 0, 1, 0, 0.5, 0.5, 0};
    EXPECT_THROW(Cell cell(rvecs, 3), singular_cell_vectors);
}

TEST(CellTest, constructor_nvec_negative) {
    double rvecs[9] = {1, 0, 1, 0, 1, 0, 0.5, 0.5, 0};
    EXPECT_THROW(Cell cell(rvecs, -1), std::domain_error);
}

TEST(CellTest, constructor_nvec_too_large) {
    double rvecs[9] = {1, 0, 1, 0, 1, 0, 0.5, 0.5, 0};
    EXPECT_THROW(Cell cell(rvecs, 4), std::domain_error);
}

TEST(CellTest, constructor_simple) {
    double rvecs[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    for(int nvec=0; nvec <=3; nvec++) {
        Cell cell(rvecs, nvec);
        SCOPED_TRACE(nvec);
        EXPECT_EQ(cell.get_nvec(), nvec);
        EXPECT_EQ(cell.get_rvec(0, 0), 1.0);
        EXPECT_EQ(cell.get_rvec(0, 1), 0.0);
        EXPECT_EQ(cell.get_rvec(0, 1), 0.0);
        EXPECT_EQ(cell.get_rvec(1, 0), 0.0);
        EXPECT_EQ(cell.get_rvec(1, 1), 1.0);
        EXPECT_EQ(cell.get_rvec(1, 2), 0.0);
        EXPECT_EQ(cell.get_rvec(2, 0), 0.0);
        EXPECT_EQ(cell.get_rvec(2, 1), 0.0);
        EXPECT_EQ(cell.get_rvec(2, 2), 1.0);
        EXPECT_EQ(cell.get_gvec(0, 0), 1.0);
        EXPECT_EQ(cell.get_gvec(0, 1), 0.0);
        EXPECT_EQ(cell.get_gvec(0, 2), 0.0);
        EXPECT_EQ(cell.get_gvec(1, 0), 0.0);
        EXPECT_EQ(cell.get_gvec(1, 1), 1.0);
        EXPECT_EQ(cell.get_gvec(1, 2), 0.0);
        EXPECT_EQ(cell.get_gvec(2, 0), 0.0);
        EXPECT_EQ(cell.get_gvec(2, 1), 0.0);
        EXPECT_EQ(cell.get_gvec(2, 2), 1.0);
        EXPECT_EQ(cell.get_volume(), nvec > 0);
        EXPECT_EQ(cell.get_rspacing(0), 1.0);
        EXPECT_EQ(cell.get_rspacing(1), 1.0);
        EXPECT_EQ(cell.get_rspacing(2), 1.0);
        EXPECT_EQ(cell.get_gspacing(0), 1.0);
        EXPECT_EQ(cell.get_gspacing(1), 1.0);
        EXPECT_EQ(cell.get_gspacing(2), 1.0);
        EXPECT_EQ(cell.is_cubic(), true);
        EXPECT_EQ(cell.is_cuboid(), true);
    }
}


// wrap
// ----

TEST_F(CellTest1, wrap_example) {
    double delta[3] = {2.5, 4.3, 3.0};
    mycell->wrap(delta);
    EXPECT_EQ(delta[0], 0.5);
    EXPECT_EQ(delta[1], 4.3);
    EXPECT_EQ(delta[2], 3.0);
}

TEST_F(CellTest2, wrap_example2) {
    double delta[3] = {2.0, 5.3, 3.0};
    mycell->wrap(delta);
    EXPECT_EQ(delta[0], 0.0);
    EXPECT_EQ(delta[1], 5.3);
    EXPECT_EQ(delta[2], -1.0);
}

TEST_F(CellTest3, wrap_example3) {
    double delta[3] = {2.0, 0.3, 3.0};
    mycell->wrap(delta);
    EXPECT_EQ(delta[0], 0.0);
    EXPECT_EQ(delta[1], 0.3);
    EXPECT_EQ(delta[2], -1.0);
}

TEST_F(CellTest1, wrap_edges) {
    double delta[3] = {-1.0, -0.5, -2.0};
    mycell->wrap(delta);
    EXPECT_EQ(delta[0],  1.0);
    EXPECT_EQ(delta[1], -0.5);
    EXPECT_EQ(delta[2], -2.0);
}

TEST_F(CellTest2, wrap_edges) {
    double delta[3] = {-1.0, -0.5, -2.0};
    mycell->wrap(delta);
    EXPECT_EQ(delta[0],  1.0);
    EXPECT_EQ(delta[1], -0.5);
    EXPECT_EQ(delta[2],  2.0);
}

TEST_F(CellTest3, wrap_edges) {
    double delta[3] = {-1.0, -0.5, -2.0};
    mycell->wrap(delta);
    EXPECT_EQ(delta[0], 1.0);
    EXPECT_EQ(delta[1], 0.5);
    EXPECT_EQ(delta[2], 2.0);
}

TEST(CellTest, wrap_random1) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(1, irep);
        double delta[3];
        fill_random_double(delta, 3, 11, 5.0);
        double frac[3];
        cell.wrap(delta);
        cell.to_frac(delta, frac);
        EXPECT_LT(frac[0], 0.5);
        EXPECT_GE(frac[0], -0.5);
    }
}

TEST(CellTest, wrap_random2) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(2, irep);
        double delta[3];
        fill_random_double(delta, 3, 11, 5.0);
        double frac[3];
        cell.wrap(delta);
        cell.to_frac(delta, frac);
        EXPECT_LT(frac[0], 0.5);
        EXPECT_LT(frac[1], 0.5);
        EXPECT_GE(frac[0], -0.5);
        EXPECT_GE(frac[1], -0.5);
    }
}

TEST(CellTest, wrap_random3) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(3, irep);
        double delta[3];
        fill_random_double(delta, 3, 11, 5.0);
        double frac[3];
        cell.wrap(delta);
        cell.to_frac(delta, frac);
        EXPECT_LT(frac[0], 0.5);
        EXPECT_LT(frac[1], 0.5);
        EXPECT_LT(frac[2], 0.5);
        EXPECT_GE(frac[0], -0.5);
        EXPECT_GE(frac[1], -0.5);
        EXPECT_GE(frac[2], -0.5);
    }
}

TEST(CellTest, wrap_consistency1) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(1, irep);
        int coeffs[1];
        fill_random_int(coeffs, 1, irep, 5);
        double frac[3];
        double cart1[3];
        double cart2[3];
        fill_random_double(frac, 3, 11, 1.0);
        cell.to_cart(frac, cart1);
        cell.to_cart(frac, cart2);
        cell.add_rvec(cart2, coeffs);
        cell.wrap(cart2);
        EXPECT_NEAR(cart1[0], cart2[0], 1e-10);
        EXPECT_NEAR(cart1[1], cart2[1], 1e-10);
        EXPECT_NEAR(cart1[2], cart2[2], 1e-10);
    }
}

TEST(CellTest, wrap_consistency2) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(2, irep);
        int coeffs[2];
        fill_random_int(coeffs, 2, irep, 5);
        double frac[3];
        double cart1[3];
        double cart2[3];
        fill_random_double(frac, 3, 11, 1.0);
        cell.to_cart(frac, cart1);
        cell.to_cart(frac, cart2);
        cell.add_rvec(cart2, coeffs);
        cell.wrap(cart2);
        EXPECT_NEAR(cart1[0], cart2[0], 1e-10);
        EXPECT_NEAR(cart1[1], cart2[1], 1e-10);
        EXPECT_NEAR(cart1[2], cart2[2], 1e-10);
    }
}

TEST(CellTest, wrap_consistency3) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(3, irep);
        int coeffs[3];
        fill_random_int(coeffs, 3, irep, 5);
        double frac[3];
        double cart1[3];
        double cart2[3];
        fill_random_double(frac, 3, 11, 1.0);
        cell.to_cart(frac, cart1);
        cell.to_cart(frac, cart2);
        cell.add_rvec(cart2, coeffs);
        cell.wrap(cart2);
        EXPECT_NEAR(cart1[0], cart2[0], 1e-10);
        EXPECT_NEAR(cart1[1], cart2[1], 1e-10);
        EXPECT_NEAR(cart1[2], cart2[2], 1e-10);
    }
}

// to_frac and to_cart
// -------------------

TEST_F(CellTest1, to_frac_example) {
    double cart[3] = {2.5, 4.3, 3.0};
    double frac[3];
    mycell->to_frac(cart, frac);
    EXPECT_NEAR(frac[0], 1.25, 1e-10);
    EXPECT_NEAR(frac[1], 4.3, 1e-10);
    EXPECT_NEAR(frac[2], 3.0, 1e-10);
}

TEST_F(CellTest2, to_frac_example) {
    double cart[3] = {2.5, 4.3, 3.0};
    double frac[3];
    mycell->to_frac(cart, frac);
    EXPECT_NEAR(frac[0], 1.25, 1e-10);
    EXPECT_NEAR(frac[1], 0.75, 1e-10);
    EXPECT_NEAR(frac[2], -4.3, 1e-10);
}

TEST_F(CellTest3, to_frac_example) {
    double cart[3] = {2.5, 4.3, 3.0};
    double frac[3];
    mycell->to_frac(cart, frac);
    EXPECT_NEAR(frac[0], 1.25, 1e-10);
    EXPECT_NEAR(frac[1], 4.3, 1e-10);
    EXPECT_NEAR(frac[2], 0.75, 1e-10);
}

TEST_F(CellTest1, to_cart_example) {
    double frac[3] = {0.5, 0.2, -1.5};
    double cart[3];
    mycell->to_cart(frac, cart);
    EXPECT_NEAR(cart[0], 1.0, 1e-10);
    EXPECT_NEAR(cart[1], 0.2, 1e-10);
    EXPECT_NEAR(cart[2], -1.5, 1e-10);
}

TEST_F(CellTest2, to_cart_example) {
    double frac[3] = {0.5, 0.2, -1.5};
    double cart[3];
    mycell->to_cart(frac, cart);
    EXPECT_NEAR(cart[0], 1.0, 1e-10);
    EXPECT_NEAR(cart[1], 1.5, 1e-10);
    EXPECT_NEAR(cart[2], 0.8, 1e-10);
}

TEST_F(CellTest3, to_cart_example) {
    double frac[3] = {0.5, 0.2, -1.5};
    double cart[3];
    mycell->to_cart(frac, cart);
    EXPECT_NEAR(cart[0], 1.0, 1e-10);
    EXPECT_NEAR(cart[1], 0.2, 1e-10);
    EXPECT_NEAR(cart[2], -6.0, 1e-10);
}

TEST(CellTest, to_cart_to_frac_consistency1) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(1, irep);
        double frac[3];
        double cart1[3];
        double cart2[3];
        fill_random_double(cart1, 3, 11, 5.0);
        cell.to_frac(cart1, frac);
        cell.to_cart(frac, cart2);
        EXPECT_NEAR(cart1[0], cart2[0], 1e-10);
        EXPECT_NEAR(cart1[1], cart2[1], 1e-10);
        EXPECT_NEAR(cart1[2], cart2[2], 1e-10);
    }
}

TEST(CellTest, to_cart_to_frac_consistency2) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(2, irep);
        double frac[3];
        double cart1[3];
        double cart2[3];
        fill_random_double(cart1, 3, 11, 5.0);
        cell.to_frac(cart1, frac);
        cell.to_cart(frac, cart2);
        EXPECT_NEAR(cart1[0], cart2[0], 1e-10);
        EXPECT_NEAR(cart1[1], cart2[1], 1e-10);
        EXPECT_NEAR(cart1[2], cart2[2], 1e-10);
    }
}

TEST(CellTest, to_cart_to_frac_consistency3) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(3, irep);
        double frac[3];
        double cart1[3];
        double cart2[3];
        fill_random_double(cart1, 3, 11, 5.0);
        cell.to_frac(cart1, frac);
        cell.to_cart(frac, cart2);
        EXPECT_NEAR(cart1[0], cart2[0], 1e-10);
        EXPECT_NEAR(cart1[1], cart2[1], 1e-10);
        EXPECT_NEAR(cart1[2], cart2[2], 1e-10);
    }
}

// g_lincomb and dot_rvecs
// -----------------------

TEST_F(CellTest1, g_lincomb_example) {
    double coeffs[3] = {2.5, 4.3, 3.0};
    double gvec[3];
    mycell->g_lincomb(coeffs, gvec);
    EXPECT_NEAR(gvec[0], 1.25, 1e-10);
    EXPECT_NEAR(gvec[1], 4.3, 1e-10);
    EXPECT_NEAR(gvec[2], 3.0, 1e-10);
}

TEST_F(CellTest2, g_lincomb_example) {
    double coeffs[3] = {2.5, 4.3, 3.0};
    double gvec[3];
    mycell->g_lincomb(coeffs, gvec);
    EXPECT_NEAR(gvec[0],  1.25, 1e-10);
    EXPECT_NEAR(gvec[1], -3.0, 1e-10);
    EXPECT_NEAR(gvec[2],  1.075, 1e-10);
}

TEST_F(CellTest3, g_lincomb_example) {
    double coeffs[3] = {2.5, 4.3, 3.0};
    double gvec[3];
    mycell->g_lincomb(coeffs, gvec);
    EXPECT_NEAR(gvec[0], 1.25, 1e-10);
    EXPECT_NEAR(gvec[1], 4.3, 1e-10);
    EXPECT_NEAR(gvec[2], 0.75, 1e-10);
}

TEST_F(CellTest1, dot_rvecs_example) {
    double gvec[3] = {0.5, 0.2, -1.5};
    double dots[3];
    mycell->dot_rvecs(gvec, dots);
    EXPECT_NEAR(dots[0],  1.0, 1e-10);
    EXPECT_NEAR(dots[1],  0.2, 1e-10);
    EXPECT_NEAR(dots[2], -1.5, 1e-10);
}

TEST_F(CellTest2, dot_rvecs_example) {
    double gvec[3] = {0.5, 0.2, -1.5};
    double dots[3];
    mycell->dot_rvecs(gvec, dots);
    EXPECT_NEAR(dots[0],  1.0, 1e-10);
    EXPECT_NEAR(dots[1], -6.0, 1e-10);
    EXPECT_NEAR(dots[2], -0.2, 1e-10);
}

TEST_F(CellTest3, dot_rvecs_example) {
    double gvec[3] = {0.5, 0.2, -1.5};
    double dots[3];
    mycell->dot_rvecs(gvec, dots);
    EXPECT_NEAR(dots[0], 1.0, 1e-10);
    EXPECT_NEAR(dots[1], 0.2, 1e-10);
    EXPECT_NEAR(dots[2], -6.0, 1e-10);
}

TEST(CellTest, g_lincomb_dot_rvecs_consistency1) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(1, irep);
        double coeffs[3];
        double gvec[3];
        double dots[3];
        fill_random_double(coeffs, 3, 11, 5.0);
        cell.g_lincomb(coeffs, gvec);
        cell.dot_rvecs(gvec, dots);
        EXPECT_NEAR(coeffs[0], dots[0], 1e-10);
        EXPECT_NEAR(coeffs[1], dots[1], 1e-10);
        EXPECT_NEAR(coeffs[2], dots[2], 1e-10);
    }
}

TEST(CellTest, g_lincomb_dot_rvecs_consistency2) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(2, irep);
        double coeffs[3];
        double gvec[3];
        double dots[3];
        fill_random_double(coeffs, 3, 11, 5.0);
        cell.g_lincomb(coeffs, gvec);
        cell.dot_rvecs(gvec, dots);
        EXPECT_NEAR(coeffs[0], dots[0], 1e-10);
        EXPECT_NEAR(coeffs[1], dots[1], 1e-10);
        EXPECT_NEAR(coeffs[2], dots[2], 1e-10);
    }
}

TEST(CellTest, g_lincomb_dot_rvecs_consistency3) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(3, irep);
        double coeffs[3];
        double gvec[3];
        double dots[3];
        fill_random_double(coeffs, 3, 11, 5.0);
        cell.g_lincomb(coeffs, gvec);
        cell.dot_rvecs(gvec, dots);
        EXPECT_NEAR(coeffs[0], dots[0], 1e-10);
        EXPECT_NEAR(coeffs[1], dots[1], 1e-10);
        EXPECT_NEAR(coeffs[2], dots[2], 1e-10);
    }
}


// add_rvec
// --------

TEST(CellTest, add_rvec_consistency1) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(1, irep);
        int coeffs[1];
        fill_random_int(coeffs, 1, irep, 5);
        double cart1[3];
        double cart2[3];
        double frac1[3];
        double frac2[3];
        fill_random_double(cart1, 3, 11, 10.0);
        cart2[0] = cart1[0];
        cart2[1] = cart1[1];
        cart2[2] = cart1[2];
        cell.add_rvec(cart2, coeffs);
        cell.to_frac(cart1, frac1);
        cell.to_frac(cart2, frac2);
        EXPECT_NEAR(frac2[0] - frac1[0], (double)coeffs[0], 1e-10);
        EXPECT_NEAR(frac2[1] - frac1[1], 0.0, 1e-10);
        EXPECT_NEAR(frac2[2] - frac1[2], 0.0, 1e-10);
    }
}

TEST(CellTest, add_rvec_consistency2) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(2, irep);
        int coeffs[2];
        fill_random_int(coeffs, 2, irep, 5);
        double cart1[3];
        double cart2[3];
        double frac1[3];
        double frac2[3];
        fill_random_double(cart1, 3, 11, 10.0);
        cart2[0] = cart1[0];
        cart2[1] = cart1[1];
        cart2[2] = cart1[2];
        cell.add_rvec(cart2, coeffs);
        cell.to_frac(cart1, frac1);
        cell.to_frac(cart2, frac2);
        EXPECT_NEAR(frac2[0] - frac1[0], (double)coeffs[0], 1e-10);
        EXPECT_NEAR(frac2[1] - frac1[1], (double)coeffs[1], 1e-10);
        EXPECT_NEAR(frac2[2] - frac1[2], 0.0, 1e-10);
    }
}

TEST(CellTest, add_rvec_consistency3) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(3, irep);
        int coeffs[3];
        fill_random_int(coeffs, 3, irep, 5);
        double cart1[3];
        double cart2[3];
        double frac1[3];
        double frac2[3];
        fill_random_double(cart1, 3, 11, 10.0);
        cart2[0] = cart1[0];
        cart2[1] = cart1[1];
        cart2[2] = cart1[2];
        cell.add_rvec(cart2, coeffs);
        cell.to_frac(cart1, frac1);
        cell.to_frac(cart2, frac2);
        EXPECT_NEAR(frac2[0] - frac1[0], (double)coeffs[0], 1e-10);
        EXPECT_NEAR(frac2[1] - frac1[1], (double)coeffs[1], 1e-10);
        EXPECT_NEAR(frac2[2] - frac1[2], (double)coeffs[2], 1e-10);
    }
}


// The getters
// -----------

// get_nvec() is already tested above

TEST(CellTest, get_rvec) {
    for (int nvec=1; nvec<=3; nvec++) {
        double rvecs[nvec*3];
        SCOPED_TRACE(nvec);
        Cell* cell = NULL;
        fill_random_double(rvecs, nvec*3, 1487, 2.0);
        cell = new Cell(rvecs, nvec);
        EXPECT_EQ(mycell->get_rvec(0, 0), rvecs[0]);
        EXPECT_EQ(mycell->get_rvec(0, 1), rvecs[1]);
        EXPECT_EQ(mycell->get_rvec(0, 2), rvecs[2]);
        if (nvec < 2) continue;
        EXPECT_EQ(mycell->get_rvec(1, 0), rvecs[3]);
        EXPECT_EQ(mycell->get_rvec(1, 1), rvecs[4]);
        EXPECT_EQ(mycell->get_rvec(1, 2), rvecs[5]);
        if (nvec < 3) continue;
        EXPECT_EQ(mycell->get_rvec(2, 0), rvecs[6]);
        EXPECT_EQ(mycell->get_rvec(2, 1), rvecs[7]);
        EXPECT_EQ(mycell->get_rvec(2, 2), rvecs[8]);
        delete cell;
    }
}

TEST(CellTest, get_domain) {
    for (int nvec=1; nvec<=3; nvec++) {
        double rvecs[nvec*3];
        SCOPED_TRACE(nvec);
        Cell* cell = NULL;
        fill_random_double(rvecs, nvec*3, 1487, 2.0);
        cell = new Cell(rvecs, nvec);
        EXPECT_THROW(mycell->get_rvec(-1, 0), std::domain_error);
        EXPECT_THROW(mycell->get_rvec(3, 0), std::domain_error);
        EXPECT_THROW(mycell->get_rvec(0, -1), std::domain_error);
        EXPECT_THROW(mycell->get_rvec(0, 3), std::domain_error);
        EXPECT_THROW(mycell->get_gvec(-1, 0), std::domain_error);
        EXPECT_THROW(mycell->get_gvec(3, 0), std::domain_error);
        EXPECT_THROW(mycell->get_gvec(0, -1), std::domain_error);
        EXPECT_THROW(mycell->get_gvec(0, 3), std::domain_error);
        EXPECT_THROW(mycell->get_rlength(-1), std::domain_error);
        EXPECT_THROW(mycell->get_rlength(3), std::domain_error);
        EXPECT_THROW(mycell->get_glength(-1), std::domain_error);
        EXPECT_THROW(mycell->get_glength(3), std::domain_error);
        EXPECT_THROW(mycell->get_rspacing(-1), std::domain_error);
        EXPECT_THROW(mycell->get_rspacing(3), std::domain_error);
        EXPECT_THROW(mycell->get_gspacing(-1), std::domain_error);
        EXPECT_THROW(mycell->get_gspacing(3), std::domain_error);
        delete cell;
    }
}

TEST_F(CellTest1, get_example) {
    EXPECT_NEAR(mycell->get_gvec(0, 0), 0.5, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(0, 1), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(0, 2), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(1, 0), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(1, 1), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(1, 2), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(2, 0), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(2, 1), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(2, 2), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_volume(), 2.0, 1e-10);
    EXPECT_NEAR(mycell->get_rlength(0), 2.0, 1e-10);
    EXPECT_NEAR(mycell->get_rlength(1), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_rlength(2), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_rspacing(0), 2.0, 1e-10);
    EXPECT_NEAR(mycell->get_rspacing(1), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_rspacing(2), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_glength(0), 0.5, 1e-10);
    EXPECT_NEAR(mycell->get_glength(1), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_glength(2), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_gspacing(0), 0.5, 1e-10);
    EXPECT_NEAR(mycell->get_gspacing(1), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_gspacing(2), 1.0, 1e-10);
}

TEST_F(CellTest2, get_example) {
    EXPECT_NEAR(mycell->get_gvec(0, 0),  0.5, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(0, 1),  0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(0, 2),  0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(1, 0),  0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(1, 1),  0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(1, 2),  0.25, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(2, 0),  0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(2, 1), -1.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(2, 2),  0.0, 1e-10);
    EXPECT_NEAR(mycell->get_volume(), 8.0, 1e-10);
    EXPECT_NEAR(mycell->get_rlength(0), 2.0, 1e-10);
    EXPECT_NEAR(mycell->get_rlength(1), 4.0, 1e-10);
    EXPECT_NEAR(mycell->get_rlength(2), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_rspacing(0), 2.0, 1e-10);
    EXPECT_NEAR(mycell->get_rspacing(1), 4.0, 1e-10);
    EXPECT_NEAR(mycell->get_rspacing(2), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_glength(0), 0.5, 1e-10);
    EXPECT_NEAR(mycell->get_glength(1), 0.25, 1e-10);
    EXPECT_NEAR(mycell->get_glength(2), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_gspacing(0), 0.5, 1e-10);
    EXPECT_NEAR(mycell->get_gspacing(1), 0.25, 1e-10);
    EXPECT_NEAR(mycell->get_gspacing(2), 1.0, 1e-10);
}

TEST_F(CellTest3, get_example) {
    EXPECT_NEAR(mycell->get_gvec(0, 0), 0.5, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(0, 1), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(0, 2), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(1, 0), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(1, 1), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(1, 2), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(2, 0), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(2, 1), 0.0, 1e-10);
    EXPECT_NEAR(mycell->get_gvec(2, 2), 0.25, 1e-10);
    EXPECT_NEAR(mycell->get_volume(), 8.0, 1e-10);
    EXPECT_NEAR(mycell->get_rlength(0), 2.0, 1e-10);
    EXPECT_NEAR(mycell->get_rlength(1), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_rlength(2), 4.0, 1e-10);
    EXPECT_NEAR(mycell->get_rspacing(0), 2.0, 1e-10);
    EXPECT_NEAR(mycell->get_rspacing(1), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_rspacing(2), 4.0, 1e-10);
    EXPECT_NEAR(mycell->get_glength(0), 0.5, 1e-10);
    EXPECT_NEAR(mycell->get_glength(1), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_glength(2), 0.25, 1e-10);
    EXPECT_NEAR(mycell->get_gspacing(0), 0.5, 1e-10);
    EXPECT_NEAR(mycell->get_gspacing(1), 1.0, 1e-10);
    EXPECT_NEAR(mycell->get_gspacing(2), 0.25, 1e-10);
}

// is_cubic and is_cuboid
// ----------------------

TEST(CellTest, cubic_cuboid_random) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(1, irep);
        EXPECT_FALSE(cell.is_cubic());
        EXPECT_FALSE(cell.is_cuboid());
    }
}

TEST_F(CellTest1, cubic_cuboid_example) {
    EXPECT_TRUE(mycell->is_cubic());
    EXPECT_TRUE(mycell->is_cuboid());
}

TEST_F(CellTest2, cubic_cuboid_example) {
    EXPECT_FALSE(mycell->is_cubic());
    EXPECT_FALSE(mycell->is_cuboid());
}

TEST_F(CellTest3, cubic_cuboid_example) {
    EXPECT_FALSE(mycell->is_cubic());
    EXPECT_TRUE(mycell->is_cuboid());
}

// set_ranges_rcut
// ---------------

TEST_F(CellTest1, set_ranges_rcut_example) {
    double center[3] = {6.3, 0.2, -0.8};
    int ranges_begin[1];
    int ranges_end[1];
    mycell->set_ranges_rcut(center, 1.0, ranges_begin, ranges_end);
    EXPECT_EQ(ranges_begin[0], 2);
    EXPECT_EQ(ranges_end[0], 4);
    mycell->set_ranges_rcut(center, 2.0, ranges_begin, ranges_end);
    EXPECT_EQ(ranges_begin[0], 2);
    EXPECT_EQ(ranges_end[0], 5);
    mycell->set_ranges_rcut(center, 3.0, ranges_begin, ranges_end);
    EXPECT_EQ(ranges_begin[0], 1);
    EXPECT_EQ(ranges_end[0], 5);
}

TEST_F(CellTest1, set_ranges_rcut_edge) {
    double center[3] = {2.0, 0.2, -0.8};
    int ranges_begin[1];
    int ranges_end[1];
    mycell->set_ranges_rcut(center, 1.0, ranges_begin, ranges_end);
    EXPECT_EQ(ranges_begin[0], 0);
    EXPECT_EQ(ranges_end[0], 2);
    mycell->set_ranges_rcut(center, 2.0, ranges_begin, ranges_end);
    EXPECT_EQ(ranges_begin[0], 0);
    EXPECT_EQ(ranges_end[0], 2);
    mycell->set_ranges_rcut(center, 3.0, ranges_begin, ranges_end);
    EXPECT_EQ(ranges_begin[0], -1);
    EXPECT_EQ(ranges_end[0], 3);
}

TEST_F(CellTest2, set_ranges_rcut_example) {
    double center[3] = {6.3, 0.2, -5.0};
    int ranges_begin[2];
    int ranges_end[2];
    mycell->set_ranges_rcut(center, 1.1, ranges_begin, ranges_end);
    EXPECT_EQ(ranges_begin[0], 2);
    EXPECT_EQ(ranges_begin[1], -2);
    EXPECT_EQ(ranges_end[0], 4);
    EXPECT_EQ(ranges_end[1], 0);
}

TEST_F(CellTest2, set_ranges_rcut_edge) {
    double center[3] = {4.0, 0.2, -2.0};
    int ranges_begin[2];
    int ranges_end[2];
    mycell->set_ranges_rcut(center, 2.0, ranges_begin, ranges_end);
    EXPECT_EQ(ranges_begin[0], 1);
    EXPECT_EQ(ranges_begin[1], -1);
    EXPECT_EQ(ranges_end[0], 3);
    EXPECT_EQ(ranges_end[1], 0);
}

TEST_F(CellTest3, set_ranges_rcut_example) {
    double center[3] = {6.3, 2.2, -5.8};
    int ranges_begin[3];
    int ranges_end[3];
    mycell->set_ranges_rcut(center, 1.0, ranges_begin, ranges_end);
    EXPECT_EQ(ranges_begin[0], 2);
    EXPECT_EQ(ranges_begin[1], 1);
    EXPECT_EQ(ranges_begin[2], -2);
    EXPECT_EQ(ranges_end[0], 4);
    EXPECT_EQ(ranges_end[1], 4);
    EXPECT_EQ(ranges_end[2], -1);
}

TEST_F(CellTest3, set_ranges_rcut_edge) {
    double center[3] = {10.0, -2.0, -6.0};
    int ranges_begin[3];
    int ranges_end[3];
    mycell->set_ranges_rcut(center, 2.0, ranges_begin, ranges_end);
    EXPECT_EQ(ranges_begin[0], 4);
    EXPECT_EQ(ranges_begin[1], -4);
    EXPECT_EQ(ranges_begin[2], -2);
    EXPECT_EQ(ranges_end[0], 6);
    EXPECT_EQ(ranges_end[1], 0);
    EXPECT_EQ(ranges_end[2], -1);
}

TEST_F(CellTest3, set_ranges_rcut_domain) {
    double center[3] = {6.3, 2.2, -5.8};
    int ranges_begin[3];
    int ranges_end[3];
    EXPECT_THROW(mycell->set_ranges_rcut(center, -1.0, ranges_begin, ranges_end), std::domain_error);
    EXPECT_THROW(mycell->set_ranges_rcut(center, 0.0, ranges_begin, ranges_end), std::domain_error);
}

TEST(CellTest, set_ranges_rcut_random1) {
    for (int icell=0; icell < 100; icell++) {
        Cell cell = create_random_cell(1, icell);
        double center[3];
        int ranges_begin[1];
        int ranges_end[1];
        double rcut = 0.3*(icell+1);
        fill_random_double(center, 3, icell+2, 5.0);
        cell.set_ranges_rcut(center, rcut, ranges_begin, ranges_end);
        for (int ipoint=0; ipoint < 1000; ipoint++) {
            // Make a random point that is uniformly distributed within a sphere, but
            // that has a high probability of setting on the edge.
            double point[3];
            double frac[3];
            double norm;
            fill_random_double(point, 3, ipoint+icell*1000, 1.0);
            norm = sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2]);
            if (norm > 1) {
                point[0] /= norm;
                point[1] /= norm;
                point[2] /= norm;
            }
            point[0] = center[0] + rcut*point[0];
            point[1] = center[1] + rcut*point[1];
            point[2] = center[2] + rcut*point[2];
            cell.to_frac(point, frac);
            EXPECT_LE(ranges_begin[0], frac[0]);
            EXPECT_GE(ranges_end[0], frac[0]);
        }
    }
}

TEST(CellTest, set_ranges_rcut_random2) {
    for (int icell=0; icell < 100; icell++) {
        Cell cell = create_random_cell(2, icell);
        double center[3];
        int ranges_begin[2];
        int ranges_end[2];
        double rcut = 0.3*(icell+1);
        fill_random_double(center, 3, icell+3, 5.0);
        cell.set_ranges_rcut(center, rcut, ranges_begin, ranges_end);
        for (int ipoint=0; ipoint < 1000; ipoint++) {
            // Make a random point that is uniformly distributed within a sphere, but
            // that has a high probability of setting on the edge.
            double point[3];
            double frac[3];
            double norm;
            fill_random_double(point, 3, ipoint+icell*1000, 1.0);
            norm = sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2]);
            if (norm > 1) {
                point[0] /= norm;
                point[1] /= norm;
                point[2] /= norm;
            }
            point[0] = center[0] + rcut*point[0];
            point[1] = center[1] + rcut*point[1];
            point[2] = center[2] + rcut*point[2];
            cell.to_frac(point, frac);
            EXPECT_LE(ranges_begin[0], frac[0]);
            EXPECT_LE(ranges_begin[1], frac[1]);
            EXPECT_GE(ranges_end[0], frac[0]);
            EXPECT_GE(ranges_end[1], frac[1]);
        }
    }
}

TEST(CellTest, set_ranges_rcut_random3) {
    for (int icell=0; icell < 100; icell++) {
        Cell cell = create_random_cell(3, icell);
        double center[3];
        int ranges_begin[3];
        int ranges_end[3];
        double rcut = 0.3*(icell+1);
        fill_random_double(center, 3, icell+5, 5.0);
        cell.set_ranges_rcut(center, rcut, ranges_begin, ranges_end);
        for (int ipoint=0; ipoint < 1000; ipoint++) {
            // Make a random point that is uniformly distributed within a sphere, but
            // that has a high probability of setting on the edge.
            double point[3];
            double frac[3];
            double norm;
            fill_random_double(point, 3, ipoint+icell*1000, 1.0);
            norm = sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2]);
            if (norm > 1) {
                point[0] /= norm;
                point[1] /= norm;
                point[2] /= norm;
            }
            point[0] = center[0] + rcut*point[0];
            point[1] = center[1] + rcut*point[1];
            point[2] = center[2] + rcut*point[2];
            cell.to_frac(point, frac);
            EXPECT_LE(ranges_begin[0], frac[0]);
            EXPECT_LE(ranges_begin[1], frac[1]);
            EXPECT_LE(ranges_begin[2], frac[2]);
            EXPECT_GE(ranges_end[0], frac[0]);
            EXPECT_GE(ranges_end[1], frac[1]);
            EXPECT_GE(ranges_end[2], frac[2]);
        }
    }
}

// select_inside
// -------------

/* TODO

- Tests with random cells should use fixtures as well: provide iterator over
  random cells or something similar.
- Allocate a bunch of 3D vectors in the fixtures that are commonly used:
  frac, frac1, frac2, cart, cart1, cart2, delta, coeffs, ... Even if most tests
  don't use all of them, it is still a major convenience.
- Let set_ranges_rcut return the total number of cells contained in the computed
  ranges => convenient for memory allocation.

- select_inside would become easier if it calls set_ranges_rcut itself. Is this
  a good idea, i.e. it forces us to do allocation of indices inside the
  function. This may not be desirable in all cases. Better solution: provide
  convenience function that does everything in one go but also provide API
  for the separate parts.
- Rename: select_inside => select_inside_rcut

TEST_F(CellTest1, select_inside_example) {
    double rcut = 5.0
    double center[3] = {2.5, 3.4, -0.6}
    int ranges_begin[1];
    int ranges_end[1];
    mycell->set_ranges_rcut(center, rcut, ranges_begin, ranges_end);

    mycell->select_inside(center, rcut, ranges_begin, ranges_end)
}

*/


// smart_wrap
// ----------

TEST(SmartWrap, examples) {
    EXPECT_EQ(smart_wrap(-15, 5, true), 0);
    EXPECT_EQ(smart_wrap( -5, 5, true), 0);
    EXPECT_EQ(smart_wrap( -3, 5, true), 2);
    EXPECT_EQ(smart_wrap( -1, 5, true), 4);
    EXPECT_EQ(smart_wrap(  0, 5, true), 0);
    EXPECT_EQ(smart_wrap(  3, 5, true), 3);
    EXPECT_EQ(smart_wrap(  5, 5, true), 0);
    EXPECT_EQ(smart_wrap(  6, 5, true), 1);
    EXPECT_EQ(smart_wrap( 10, 5, true), 0);
    EXPECT_EQ(smart_wrap( 12, 5, true), 2);
    EXPECT_EQ(smart_wrap(-15, 5, false), -1);
    EXPECT_EQ(smart_wrap( -5, 5, false), -1);
    EXPECT_EQ(smart_wrap( -3, 5, false), -1);
    EXPECT_EQ(smart_wrap( -1, 5, false), -1);
    EXPECT_EQ(smart_wrap(  0, 5, false),  0);
    EXPECT_EQ(smart_wrap(  3, 5, false),  3);
    EXPECT_EQ(smart_wrap(  4, 5, false),  4);
    EXPECT_EQ(smart_wrap(  5, 5, false), -1);
    EXPECT_EQ(smart_wrap(  6, 5, false), -1);
    EXPECT_EQ(smart_wrap( 10, 5, false), -1);
    EXPECT_EQ(smart_wrap( 12, 5, false), -1);
}
