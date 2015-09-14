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
#include <cmath>
#include "celllists/cell.h"
#include "celllists/vec3.h"
#include "common.h"


// Fixtures
// ========

class CellTest : public ::testing::Test {
    public:
        int nvec;
        Cell* mycell;
        double myrvecs[9];
        double singrvecs[9];

        Cell* create_random_cell(unsigned int seed, double scale=1.0, bool cuboid=false) {
            return create_random_cell_nvec(seed, nvec, scale, cuboid);
        }

        virtual void SetUp() = 0;

        void set_up_data() {
            // Example tests
            std::fill(myrvecs, myrvecs+9, 0);
            myrvecs[0] = 2;
            myrvecs[4] = 1;
            myrvecs[8] = 4;
            if (nvec == 2) {
                myrvecs[4] = 0.0;
                myrvecs[5] = 4.0;
            }
            mycell = new Cell(myrvecs, nvec);
            // Singular cell vectors
            std::fill(singrvecs, singrvecs+9, 0);
            if (nvec > 1) {
                singrvecs[0] = 1.0;
                singrvecs[3] = 0.5;
            }
            if (nvec == 3) {
                singrvecs[3] = 0.0;
                singrvecs[4] = 2.0;
                singrvecs[6] = 0.5;
                singrvecs[7] = 0.8;
            }
        }
        virtual void TearDown() {
            delete mycell;
        }
};

class CellTestP : public CellTest,
                  public ::testing::WithParamInterface<int> {
    public:
        virtual void SetUp() {
            nvec = GetParam();
            set_up_data();
        }
};

class CellTest1 : public CellTest {
    public:
        virtual void SetUp() {
            nvec = 1;
            set_up_data();
        }
};


class CellTest2 : public CellTest {
    public:
        virtual void SetUp() {
            nvec = 2;
            set_up_data();
        }
};


class CellTest3 : public CellTest {
    public:
        virtual void SetUp() {
            nvec = 3;
            set_up_data();
        }
};


// Tests grouped by main method being tested
// =========================================

// Constructor
// -----------

TEST_P(CellTestP, constructor_singular) {
    EXPECT_THROW(Cell cell(singrvecs, GetParam()), singular_cell_vectors);
}

TEST_F(CellTest3, constructor_nvec_negative) {
    double rvecs[9] = {1, 0, 1, 0, 1, 0, 0.5, 0.5, 0};
    EXPECT_THROW(Cell cell(rvecs, -1), std::domain_error);
}

TEST_F(CellTest3, constructor_nvec_too_large) {
    double rvecs[9] = {1, 0, 1, 0, 1, 0, 0.5, 0.5, 0};
    EXPECT_THROW(Cell cell(rvecs, 4), std::domain_error);
}

TEST_P(CellTestP, constructor_simple) {
    double rvecs[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    Cell cell(rvecs, nvec);
    EXPECT_EQ(nvec, cell.get_nvec());
    EXPECT_EQ(1.0, cell.get_rvec(0, 0));
    EXPECT_EQ(0.0, cell.get_rvec(0, 1));
    EXPECT_EQ(0.0, cell.get_rvec(0, 1));
    EXPECT_EQ(0.0, cell.get_rvec(1, 0));
    EXPECT_EQ(1.0, cell.get_rvec(1, 1));
    EXPECT_EQ(0.0, cell.get_rvec(1, 2));
    EXPECT_EQ(0.0, cell.get_rvec(2, 0));
    EXPECT_EQ(0.0, cell.get_rvec(2, 1));
    EXPECT_EQ(1.0, cell.get_rvec(2, 2));
    EXPECT_EQ(1.0, cell.get_gvec(0, 0));
    EXPECT_EQ(0.0, cell.get_gvec(0, 1));
    EXPECT_EQ(0.0, cell.get_gvec(0, 2));
    EXPECT_EQ(0.0, cell.get_gvec(1, 0));
    EXPECT_EQ(1.0, cell.get_gvec(1, 1));
    EXPECT_EQ(0.0, cell.get_gvec(1, 2));
    EXPECT_EQ(0.0, cell.get_gvec(2, 0));
    EXPECT_EQ(0.0, cell.get_gvec(2, 1));
    EXPECT_EQ(1.0, cell.get_gvec(2, 2));
    EXPECT_EQ(nvec > 0, cell.get_volume());
    EXPECT_EQ(1.0, cell.get_rspacing(0));
    EXPECT_EQ(1.0, cell.get_rspacing(1));
    EXPECT_EQ(1.0, cell.get_rspacing(2));
    EXPECT_EQ(1.0, cell.get_gspacing(0));
    EXPECT_EQ(1.0, cell.get_gspacing(1));
    EXPECT_EQ(1.0, cell.get_gspacing(2));
    EXPECT_EQ(true, cell.is_cubic());
    EXPECT_EQ(true, cell.is_cuboid());
}


// wrap
// ----

TEST_F(CellTest1, wrap_example) {
    double delta[3] = {2.5, 4.3, 3.0};
    mycell->wrap(delta);
    EXPECT_EQ(0.5, delta[0]);
    EXPECT_EQ(4.3, delta[1]);
    EXPECT_EQ(3.0, delta[2]);
}

TEST_F(CellTest2, wrap_example) {
    double delta[3] = {2.0, 5.3, 3.0};
    mycell->wrap(delta);
    EXPECT_EQ(0.0, delta[0]);
    EXPECT_EQ(5.3, delta[1]);
    EXPECT_EQ(-1.0, delta[2]);
}

TEST_F(CellTest3, wrap_example) {
    double delta[3] = {2.0, 0.3, 3.0};
    mycell->wrap(delta);
    EXPECT_EQ(0.0, delta[0]);
    EXPECT_EQ(0.3, delta[1]);
    EXPECT_EQ(-1.0, delta[2]);
}

TEST_F(CellTest1, wrap_edges) {
    double delta[3] = {-1.0, -0.5, -2.0};
    mycell->wrap(delta);
    EXPECT_EQ(1.0, delta[0]);
    EXPECT_EQ(-0.5, delta[1]);
    EXPECT_EQ(-2.0, delta[2]);
}

TEST_F(CellTest2, wrap_edges) {
    double delta[3] = {-1.0, -0.5, -2.0};
    mycell->wrap(delta);
    EXPECT_EQ(1.0, delta[0]);
    EXPECT_EQ(-0.5, delta[1]);
    EXPECT_EQ(2.0, delta[2]);
}

TEST_F(CellTest3, wrap_edges) {
    double delta[3] = {-1.0, -0.5, -2.0};
    mycell->wrap(delta);
    EXPECT_EQ(1.0, delta[0]);
    EXPECT_EQ(0.5, delta[1]);
    EXPECT_EQ(2.0, delta[2]);
}

TEST_P(CellTestP, wrap_random) {
    int num_wrapped = 0;
    for (int irep=0; irep < NREP; irep++) {
        Cell* cell = create_random_cell(irep);
        double delta[3];
        fill_random_double(irep+NREP, delta, 3, -6.0, 6.0);
        double frac[3];

        // For test sufficiency check
        cell->to_rfrac(delta, frac);
        for (int ivec=0; ivec < nvec; ivec++) {
            if (fabs(frac[ivec] > 0.5)) num_wrapped++;
        }

        // Actual test
        cell->wrap(delta);
        cell->to_rfrac(delta, frac);
        for (int ivec=0; ivec < nvec; ivec++) {
            EXPECT_LT(frac[ivec], 0.5);
            EXPECT_GE(frac[ivec], -0.5);
        }
        delete cell;
    }
    // Check whether the test is sufficient.
    EXPECT_LT((NREP*nvec)/3, num_wrapped);
}

TEST_P(CellTestP, wrap_consistency) {
    for (int irep=0; irep < NREP; irep++) {
        Cell* cell = create_random_cell(irep);
        int coeffs[nvec];
        fill_random_int(irep, coeffs, nvec, -5, 5);
        double frac[3];
        double cart1[3];
        double cart2[3];
        fill_random_double(irep, frac, 3);
        cell->to_rcart(frac, cart1);
        cell->to_rcart(frac, cart2);
        cell->add_rvec(cart2, coeffs);
        cell->wrap(cart2);
        EXPECT_NEAR(cart2[0], cart1[0], 1e-10);
        EXPECT_NEAR(cart2[1], cart1[1], 1e-10);
        EXPECT_NEAR(cart2[2], cart1[2], 1e-10);
        delete cell;
    }
}

// to_rfrac and to_rcart
// ---------------------

TEST_F(CellTest1, to_rfrac_example) {
    double rcart[3] = {2.5, 4.3, 3.0};
    double rfrac[3];
    mycell->to_rfrac(rcart, rfrac);
    EXPECT_NEAR(1.25, rfrac[0], 1e-10);
    EXPECT_NEAR(4.3, rfrac[1], 1e-10);
    EXPECT_NEAR(3.0, rfrac[2], 1e-10);
}

TEST_F(CellTest2, to_rfrac_example) {
    double rcart[3] = {2.5, 4.3, 3.0};
    double rfrac[3];
    mycell->to_rfrac(rcart, rfrac);
    EXPECT_NEAR(1.25, rfrac[0], 1e-10);
    EXPECT_NEAR(0.75, rfrac[1], 1e-10);
    EXPECT_NEAR(-4.3, rfrac[2], 1e-10);
}

TEST_F(CellTest3, to_rfrac_example) {
    double rcart[3] = {2.5, 4.3, 3.0};
    double rfrac[3];
    mycell->to_rfrac(rcart, rfrac);
    EXPECT_NEAR(1.25, rfrac[0], 1e-10);
    EXPECT_NEAR(4.3, rfrac[1], 1e-10);
    EXPECT_NEAR(0.75, rfrac[2], 1e-10);
}

TEST_F(CellTest1, to_rcart_example) {
    double rfrac[3] = {0.5, 0.2, -1.5};
    double rcart[3];
    mycell->to_rcart(rfrac, rcart);
    EXPECT_NEAR(1.0, rcart[0], 1e-10);
    EXPECT_NEAR(0.2, rcart[1], 1e-10);
    EXPECT_NEAR(-1.5, rcart[2], 1e-10);
}

TEST_F(CellTest2, to_rcart_example) {
    double rfrac[3] = {0.5, 0.2, -1.5};
    double rcart[3];
    mycell->to_rcart(rfrac, rcart);
    EXPECT_NEAR(1.0, rcart[0], 1e-10);
    EXPECT_NEAR(1.5, rcart[1], 1e-10);
    EXPECT_NEAR(0.8, rcart[2], 1e-10);
}

TEST_F(CellTest3, to_rcart_example) {
    double rfrac[3] = {0.5, 0.2, -1.5};
    double rcart[3];
    mycell->to_rcart(rfrac, rcart);
    EXPECT_NEAR(1.0, rcart[0], 1e-10);
    EXPECT_NEAR(0.2, rcart[1], 1e-10);
    EXPECT_NEAR(-6.0, rcart[2], 1e-10);
}

TEST_P(CellTestP, to_rcart_to_rfrac_consistency) {
    for (int irep=0; irep < NREP; irep++) {
        Cell* cell = create_random_cell(irep);
        double rfrac[3];
        double rcart1[3];
        double rcart2[3];
        fill_random_double(irep, rcart1, 3, -5.0, 5.0);
        cell->to_rfrac(rcart1, rfrac);
        cell->to_rcart(rfrac, rcart2);
        EXPECT_NEAR(rcart2[0], rcart1[0], 1e-10);
        EXPECT_NEAR(rcart2[1], rcart1[1], 1e-10);
        EXPECT_NEAR(rcart2[2], rcart1[2], 1e-10);
        delete cell;
    }
}

// to_gcart and to_gfrac
// ---------------------

TEST_F(CellTest1, to_gcart_example) {
    double gfrac[3] = {2.5, 4.3, 3.0};
    double gcart[3];
    mycell->to_gcart(gfrac, gcart);
    EXPECT_NEAR(1.25, gcart[0], 1e-10);
    EXPECT_NEAR(4.3, gcart[1], 1e-10);
    EXPECT_NEAR(3.0, gcart[2], 1e-10);
}

TEST_F(CellTest2, to_gcart_example) {
    double gfrac[3] = {2.5, 4.3, 3.0};
    double gcart[3];
    mycell->to_gcart(gfrac, gcart);
    EXPECT_NEAR(1.25, gcart[0], 1e-10);
    EXPECT_NEAR(-3.0, gcart[1], 1e-10);
    EXPECT_NEAR(1.075, gcart[2], 1e-10);
}

TEST_F(CellTest3, to_gcart_example) {
    double gfrac[3] = {2.5, 4.3, 3.0};
    double gcart[3];
    mycell->to_gcart(gfrac, gcart);
    EXPECT_NEAR(1.25, gcart[0], 1e-10);
    EXPECT_NEAR(4.3, gcart[1], 1e-10);
    EXPECT_NEAR(0.75, gcart[2], 1e-10);
}

TEST_F(CellTest1, to_gfrac_example) {
    double gcart[3] = {0.5, 0.2, -1.5};
    double gfrac[3];
    mycell->to_gfrac(gcart, gfrac);
    EXPECT_NEAR(1.0, gfrac[0], 1e-10);
    EXPECT_NEAR(0.2, gfrac[1], 1e-10);
    EXPECT_NEAR(-1.5, gfrac[2], 1e-10);
}

TEST_F(CellTest2, to_gfrac_example) {
    double gcart[3] = {0.5, 0.2, -1.5};
    double gfrac[3];
    mycell->to_gfrac(gcart, gfrac);
    EXPECT_NEAR(1.0, gfrac[0], 1e-10);
    EXPECT_NEAR(-6.0, gfrac[1], 1e-10);
    EXPECT_NEAR(-0.2, gfrac[2], 1e-10);
}

TEST_F(CellTest3, to_gfrac_example) {
    double gcart[3] = {0.5, 0.2, -1.5};
    double gfrac[3];
    mycell->to_gfrac(gcart, gfrac);
    EXPECT_NEAR(1.0, gfrac[0], 1e-10);
    EXPECT_NEAR(0.2, gfrac[1], 1e-10);
    EXPECT_NEAR(-6.0, gfrac[2], 1e-10);
}

TEST_P(CellTestP, to_gcart_to_gfrac_consistency) {
    for (int irep=0; irep < NREP; irep++) {
        Cell* cell = create_random_cell(irep);
        double gfrac1[3];
        double gcart[3];
        double gfrac2[3];
        fill_random_double(irep, gfrac1, 3, -5.0, 5.0);
        cell->to_gcart(gfrac1, gcart);
        cell->to_gfrac(gcart, gfrac2);
        EXPECT_NEAR(gfrac2[0], gfrac1[0], 1e-10);
        EXPECT_NEAR(gfrac2[1], gfrac1[1], 1e-10);
        EXPECT_NEAR(gfrac2[2], gfrac1[2], 1e-10);
        delete cell;
    }
}

// add_rvec
// --------

TEST_P(CellTestP, add_rvec_consistency) {
    for (int irep=0; irep < NREP; irep++) {
        Cell* cell = create_random_cell(irep);
        int coeffs[nvec];
        fill_random_int(irep, coeffs, nvec, -5, 5);
        double cart1[3];
        double cart2[3];
        double frac1[3];
        double frac2[3];
        fill_random_double(irep, cart1, 3, -10.0, 10.0);
        cart2[0] = cart1[0];
        cart2[1] = cart1[1];
        cart2[2] = cart1[2];
        cell->add_rvec(cart2, coeffs);
        cell->to_rfrac(cart1, frac1);
        cell->to_rfrac(cart2, frac2);
        for (int ivec=0; ivec < nvec; ivec++) {
            EXPECT_NEAR(coeffs[ivec], frac2[ivec] - frac1[ivec], 1e-10);
        }
        for (int ivec=nvec; ivec < 3; ivec++) {
            EXPECT_NEAR(0.0, frac2[ivec] - frac1[ivec], 1e-10);
        }
        delete cell;
    }
}

// The getters
// -----------

// get_nvec() is already tested above

TEST_P(CellTestP, get_rvec) {
    double rvecs[nvec*3];
    Cell* cell = NULL;
    while (true) {
        try {
            fill_random_double(1487, rvecs, nvec*3, -2.0, 2.0);
            cell = new Cell(rvecs, nvec);
            break;
        } catch (singular_cell_vectors) {}
    }
    for (int ivec=0; ivec < nvec; ivec++) {
        EXPECT_EQ(rvecs[3*ivec+0], cell->get_rvec(ivec, 0));
        EXPECT_EQ(rvecs[3*ivec+1], cell->get_rvec(ivec, 1));
        EXPECT_EQ(rvecs[3*ivec+2], cell->get_rvec(ivec, 2));
    }
    delete cell;
}

TEST_P(CellTestP, get_domain) {
    double rvecs[nvec*3];
    Cell* cell = NULL;
    fill_random_double(1487, rvecs, nvec*3, -2.0, 2.0);
    cell = new Cell(rvecs, nvec);
    EXPECT_THROW(cell->get_rvec(-1, 0), std::domain_error);
    EXPECT_THROW(cell->get_rvec(3, 0), std::domain_error);
    EXPECT_THROW(cell->get_rvec(0, -1), std::domain_error);
    EXPECT_THROW(cell->get_rvec(0, 3), std::domain_error);
    EXPECT_THROW(cell->get_gvec(-1, 0), std::domain_error);
    EXPECT_THROW(cell->get_gvec(3, 0), std::domain_error);
    EXPECT_THROW(cell->get_gvec(0, -1), std::domain_error);
    EXPECT_THROW(cell->get_gvec(0, 3), std::domain_error);
    EXPECT_THROW(cell->get_rlength(-1), std::domain_error);
    EXPECT_THROW(cell->get_rlength(3), std::domain_error);
    EXPECT_THROW(cell->get_glength(-1), std::domain_error);
    EXPECT_THROW(cell->get_glength(3), std::domain_error);
    EXPECT_THROW(cell->get_rspacing(-1), std::domain_error);
    EXPECT_THROW(cell->get_rspacing(3), std::domain_error);
    EXPECT_THROW(cell->get_gspacing(-1), std::domain_error);
    EXPECT_THROW(cell->get_gspacing(3), std::domain_error);
    delete cell;
}

TEST_F(CellTest1, get_example) {
    EXPECT_NEAR(0.5, mycell->get_gvec(0, 0), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(0, 1), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(0, 2), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(1, 0), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_gvec(1, 1), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(1, 2), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(2, 0), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(2, 1), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_gvec(2, 2), 1e-10);
    EXPECT_NEAR(2.0, mycell->get_volume(), 1e-10);
    EXPECT_NEAR(2.0, mycell->get_rlength(0), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_rlength(1), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_rlength(2), 1e-10);
    EXPECT_NEAR(2.0, mycell->get_rspacing(0), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_rspacing(1), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_rspacing(2), 1e-10);
    EXPECT_NEAR(0.5, mycell->get_glength(0), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_glength(1), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_glength(2), 1e-10);
    EXPECT_NEAR(0.5, mycell->get_gspacing(0), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_gspacing(1), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_gspacing(2), 1e-10);
}

TEST_F(CellTest2, get_example) {
    EXPECT_NEAR(0.5, mycell->get_gvec(0, 0), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(0, 1), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(0, 2), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(1, 0), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(1, 1), 1e-10);
    EXPECT_NEAR(0.25, mycell->get_gvec(1, 2), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(2, 0), 1e-10);
    EXPECT_NEAR(-1.0, mycell->get_gvec(2, 1), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(2, 2), 1e-10);
    EXPECT_NEAR(8.0, mycell->get_volume(), 1e-10);
    EXPECT_NEAR(2.0, mycell->get_rlength(0), 1e-10);
    EXPECT_NEAR(4.0, mycell->get_rlength(1), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_rlength(2), 1e-10);
    EXPECT_NEAR(2.0, mycell->get_rspacing(0), 1e-10);
    EXPECT_NEAR(4.0, mycell->get_rspacing(1), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_rspacing(2), 1e-10);
    EXPECT_NEAR(0.5, mycell->get_glength(0), 1e-10);
    EXPECT_NEAR(0.25, mycell->get_glength(1), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_glength(2), 1e-10);
    EXPECT_NEAR(0.5, mycell->get_gspacing(0), 1e-10);
    EXPECT_NEAR(0.25, mycell->get_gspacing(1), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_gspacing(2), 1e-10);
}

TEST_F(CellTest3, get_example) {
    EXPECT_NEAR(0.5, mycell->get_gvec(0, 0), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(0, 1), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(0, 2), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(1, 0), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_gvec(1, 1), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(1, 2), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(2, 0), 1e-10);
    EXPECT_NEAR(0.0, mycell->get_gvec(2, 1), 1e-10);
    EXPECT_NEAR(0.25, mycell->get_gvec(2, 2), 1e-10);
    EXPECT_NEAR(8.0, mycell->get_volume(), 1e-10);
    EXPECT_NEAR(2.0, mycell->get_rlength(0), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_rlength(1), 1e-10);
    EXPECT_NEAR(4.0, mycell->get_rlength(2), 1e-10);
    EXPECT_NEAR(2.0, mycell->get_rspacing(0), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_rspacing(1), 1e-10);
    EXPECT_NEAR(4.0, mycell->get_rspacing(2), 1e-10);
    EXPECT_NEAR(0.5, mycell->get_glength(0), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_glength(1), 1e-10);
    EXPECT_NEAR(0.25, mycell->get_glength(2), 1e-10);
    EXPECT_NEAR(0.5, mycell->get_gspacing(0), 1e-10);
    EXPECT_NEAR(1.0, mycell->get_gspacing(1), 1e-10);
    EXPECT_NEAR(0.25, mycell->get_gspacing(2), 1e-10);
}

// is_cubic and is_cuboid
// ----------------------

TEST_P(CellTestP, cubic_cuboid_random) {
    for (int irep=0; irep < NREP; irep++) {
        Cell* cell = create_random_cell(irep);
        EXPECT_FALSE(cell->is_cubic());
        EXPECT_FALSE(cell->is_cuboid());
        delete cell;
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
    int ncell = 0;
    ncell = mycell->set_ranges_rcut(center, 1.0, ranges_begin, ranges_end);
    EXPECT_EQ(2, ncell);
    EXPECT_EQ(2, ranges_begin[0]);
    EXPECT_EQ(4, ranges_end[0]);
    ncell = mycell->set_ranges_rcut(center, 2.0, ranges_begin, ranges_end);
    EXPECT_EQ(3, ncell);
    EXPECT_EQ(2, ranges_begin[0]);
    EXPECT_EQ(5, ranges_end[0]);
    ncell = mycell->set_ranges_rcut(center, 3.0, ranges_begin, ranges_end);
    EXPECT_EQ(4, ncell);
    EXPECT_EQ(1, ranges_begin[0]);
    EXPECT_EQ(5, ranges_end[0]);
}

TEST_F(CellTest1, set_ranges_rcut_edge) {
    double center[3] = {2.0, 0.2, -0.8};
    int ranges_begin[1];
    int ranges_end[1];
    int ncell = 0;
    ncell = mycell->set_ranges_rcut(center, 1.0, ranges_begin, ranges_end);
    EXPECT_EQ(2, ncell);
    EXPECT_EQ(0, ranges_begin[0]);
    EXPECT_EQ(2, ranges_end[0]);
    ncell = mycell->set_ranges_rcut(center, 2.0, ranges_begin, ranges_end);
    EXPECT_EQ(2, ncell);
    EXPECT_EQ(0, ranges_begin[0]);
    EXPECT_EQ(2, ranges_end[0]);
    ncell = mycell->set_ranges_rcut(center, 3.0, ranges_begin, ranges_end);
    EXPECT_EQ(4, ncell);
    EXPECT_EQ(-1, ranges_begin[0]);
    EXPECT_EQ(3, ranges_end[0]);
}

TEST_F(CellTest2, set_ranges_rcut_example) {
    double center[3] = {6.3, 0.2, -5.0};
    int ranges_begin[2];
    int ranges_end[2];
    int ncell = 0;
    ncell = mycell->set_ranges_rcut(center, 1.1, ranges_begin, ranges_end);
    EXPECT_EQ(2*2, ncell);
    EXPECT_EQ(2, ranges_begin[0]);
    EXPECT_EQ(-2, ranges_begin[1]);
    EXPECT_EQ(4, ranges_end[0]);
    EXPECT_EQ(0, ranges_end[1]);
}

TEST_F(CellTest2, set_ranges_rcut_edge) {
    double center[3] = {4.0, 0.2, -2.0};
    int ranges_begin[2];
    int ranges_end[2];
    int ncell = 0;
    ncell = mycell->set_ranges_rcut(center, 2.0, ranges_begin, ranges_end);
    EXPECT_EQ(2, ncell);
    EXPECT_EQ(1, ranges_begin[0]);
    EXPECT_EQ(-1, ranges_begin[1]);
    EXPECT_EQ(3, ranges_end[0]);
    EXPECT_EQ(0, ranges_end[1]);
}

TEST_F(CellTest3, set_ranges_rcut_example) {
    double center[3] = {6.3, 2.2, -5.8};
    int ranges_begin[3];
    int ranges_end[3];
    int ncell = 0;
    ncell = mycell->set_ranges_rcut(center, 1.0, ranges_begin, ranges_end);
    EXPECT_EQ(2*3*1, ncell);
    EXPECT_EQ(2, ranges_begin[0]);
    EXPECT_EQ(1, ranges_begin[1]);
    EXPECT_EQ(-2, ranges_begin[2]);
    EXPECT_EQ(4, ranges_end[0]);
    EXPECT_EQ(4, ranges_end[1]);
    EXPECT_EQ(-1, ranges_end[2]);
}

TEST_F(CellTest3, set_ranges_rcut_edge) {
    double center[3] = {10.0, -2.0, -6.0};
    int ranges_begin[3];
    int ranges_end[3];
    int ncell = 0;
    ncell = mycell->set_ranges_rcut(center, 2.0, ranges_begin, ranges_end);
    EXPECT_EQ(2*4*1, ncell);
    EXPECT_EQ(4, ranges_begin[0]);
    EXPECT_EQ(-4, ranges_begin[1]);
    EXPECT_EQ(-2, ranges_begin[2]);
    EXPECT_EQ(6, ranges_end[0]);
    EXPECT_EQ(0, ranges_end[1]);
    EXPECT_EQ(-1, ranges_end[2]);
}

TEST_P(CellTestP, set_ranges_rcut_domain) {
    double center[3] = {6.3, 2.2, -5.8};
    int ranges_begin[nvec];
    int ranges_end[nvec];
    EXPECT_THROW(mycell->set_ranges_rcut(center, -1.0, ranges_begin, ranges_end), std::domain_error);
    EXPECT_THROW(mycell->set_ranges_rcut(center, 0.0, ranges_begin, ranges_end), std::domain_error);
}

TEST_P(CellTestP, set_ranges_rcut_random) {
    int npoint_total = 0;
    for (int icell=0; icell < NREP; icell++) {
        Cell* cell = create_random_cell(icell);
        double center[3];
        int ranges_begin[nvec];
        int ranges_end[nvec];
        double rcut = 0.3*(icell+1);
        fill_random_double(icell+2, center, 3, -5.0, 5.0);
        cell->set_ranges_rcut(center, rcut, ranges_begin, ranges_end);
        for (int ipoint=0; ipoint < NPOINT; ipoint++) {
            double point[3];
            double norm;
            random_point(ipoint+icell*NPOINT, point, rcut, center, norm);
            if (norm <= rcut) {
                double frac[3];
                cell->to_rfrac(point, frac);
                for (int ivec=0; ivec < nvec; ivec++) {
                    EXPECT_LE(ranges_begin[ivec], frac[ivec]);
                    EXPECT_GE(ranges_end[ivec], frac[ivec]);
                }
                npoint_total++;
            }
        }
        delete cell;
    }
    // Check sufficiency
    EXPECT_LT((NREP*NPOINT)/3, npoint_total);
}

// select_inside
// -------------

TEST_P(CellTestP, select_inside_rcut_domain) {
    double center[3] = {2.5, 3.4, -0.6};
    int shape[3] = {10, 10, 10};
    bool pbc[3] = {true, true, true};
    EXPECT_THROW(mycell->select_inside_rcut(center, 0.0, shape, pbc, NULL), std::domain_error);
    EXPECT_THROW(mycell->select_inside_rcut(center, -1.0, shape, pbc, NULL), std::domain_error);
    Cell zero_cell = Cell(NULL, 0);
    EXPECT_THROW(zero_cell.select_inside_rcut(center, 1.0, shape, pbc, NULL), std::domain_error);
}

TEST_F(CellTest1, select_inside_rcut_example) {
    // All the parameters
    double rcut = 5.0;
    double center[3] = {2.5, 3.4, -0.6};
    int shape[1] = {10};
    bool pbc[1] = {true};

    // Call
    int nbar = mycell->select_inside_rcut(center, rcut, shape, pbc, NULL);
    ASSERT_EQ(1, nbar);
    int bars[2*nbar];
    nbar = mycell->select_inside_rcut(center, rcut, shape, pbc, bars);

    // Check results
    // lower end: -2.5 (-2 #> 8)
    // upper end:  7.5 (8 #> 4) non-inclusive
    EXPECT_EQ(1, nbar);
    EXPECT_EQ(-2, bars[0]);
    EXPECT_EQ(4, bars[1]);
}


TEST_F(CellTest2, select_inside_rcut_example) {
    // All the parameters
    double rcut = 5.0;
    double center[3] = {2.5, 3.4, -0.6};
    int shape[2] = {10, 5};
    bool pbc[2] = {true, false};

    // Call
    int nbar = mycell->select_inside_rcut(center, rcut, shape, pbc, NULL);
    ASSERT_EQ(6, nbar);
    int bars[3*nbar];
    nbar = mycell->select_inside_rcut(center, rcut, shape, pbc, bars);
    EXPECT_EQ(6, nbar);

    // Test
    for (int ibar=0; ibar<nbar; ibar++) {
        EXPECT_EQ(ibar-2, bars[3*ibar]);
        EXPECT_EQ(0, bars[3*ibar+1]);
        if ((ibar == 0) || (ibar == 1) || (ibar == 5)) {
            EXPECT_EQ(1, bars[3*ibar+2]);
        } else {
            EXPECT_EQ(2, bars[3*ibar+2]);
        }
    }
}

TEST_F(CellTest3, select_inside_rcut_example) {
    // All the parameters
    double rcut = 1.9;
    double center[3] = {2.0, 2.0, 2.0};
    int shape[3] = {10, 5, 7};
    bool pbc[3] = {true, true, true};

    // Call
    int nbar = mycell->select_inside_rcut(center, rcut, shape, pbc, NULL);
    EXPECT_EQ(8, nbar);
    int bars[4*nbar];
    nbar = mycell->select_inside_rcut(center, rcut, shape, pbc, bars);
    EXPECT_EQ(8, nbar);

    // Test
    for (int ibar=0; ibar<nbar; ibar++) {
        int* bar = bars + ibar*4;
        EXPECT_EQ(bar[0], ibar/4);
        EXPECT_EQ(bar[1], ibar%4);
        EXPECT_EQ(bar[2], 0);
        EXPECT_EQ(bar[3], 1);
    }
}

TEST_P(CellTestP, select_inside_rcut_random) {
    int nbar_total = 0;
    for (int irep=0; irep < NREP; irep++) {
        // Test parameters:
        // - Random cell
        Cell* cell = create_random_cell(2*irep, 1.0, false);
        // - Increasing rcut
        double rcut = (irep+1)*0.1;
        // - Random center
        double center[3];
        fill_random_double(47332+irep, center, 3, -1.0, 1.0);
        // - Alternating values for shape and pbc
        int shape[nvec];
        bool pbc[nvec];
        for (int ivec=0; ivec < nvec; ivec++) {
            shape[ivec] = ((irep*(ivec+1)) % 5) + 1;
            pbc[ivec] = (irep >> ivec) % 2;
        }

        // Compute the bars.
        int nbar1 = cell->select_inside_rcut(center, rcut, shape, pbc, NULL);
        int bars[(nvec+2)*nbar1];
        int nbar2 = cell->select_inside_rcut(center, rcut, shape, pbc, bars);
        EXPECT_EQ(nbar1, nbar2);
        nbar_total += nbar2;

        // Construct a random vector in a cubic box around the cutoff sphere.
        double cart[3];
        fill_random_double(123+irep, cart, 3, -rcut*1.1, rcut*1.1);
        double norm = vec3::norm(cart);
        // Center of the box must coincide with center of the sphere.
        cart[0] += center[0];
        cart[1] += center[1];
        cart[2] += center[2];
        // For the rest of the test, we need this random vector in fractional coordinates.
        double frac[3];
        cell->to_rfrac(cart, frac);

        // Does the fractional coordinate fit in one of the bars?
        int index[3] = {
            static_cast<int>(floor(frac[0])),
            static_cast<int>(floor(frac[1])),
            static_cast<int>(floor(frac[2]))
        };
        bool in_bar = false;
        for (int ibar=0; ibar < nbar2; ibar++) {
            int* bar = bars + (nvec+1)*ibar;
            if (nvec > 1) {
                if (bar[0] != index[0])
                    continue;
            }
            if (nvec > 2) {
                if (bar[1] != index[1])
                    continue;
            }
            if (index[nvec-1] < bar[nvec-1])
                continue;
            if (index[nvec-1] >= bar[nvec])
                continue;
            in_bar = true;
            break;
        }

        // Does the relative vector sit in the cutoff sphere, taking into account
        // non-periodic boundaries that truncate the cutoff sphere.
        bool in_sphere = (norm < rcut);
        if (in_sphere) {
            for (int ivec=0; ivec < nvec; ivec++) {
                if (!pbc[ivec]) {
                    if (frac[ivec] < 0) {
                        in_sphere = false;
                        break;
                    }
                    if (frac[ivec] >= shape[ivec]) {
                        in_sphere = false;
                        break;
                    }
                }
            }
        }

        // First test: if the vector is in the cutoff and non-periodic boundaries, the
        //             point must be in a bar.
        if (in_sphere) {
            EXPECT_TRUE(in_bar);
        }
        // Second test: if not in a bar, the norm must be larger than the cutoff, or the
        //              point is outside a non-periodic boundary.
        else if (!in_bar) {
            EXPECT_FALSE(in_sphere);
        }
        // Clean up
        delete cell;
    }
    // Sufficiency check
    EXPECT_LE(NREP*((nvec-1)*3 + 1), nbar_total);
}


TEST_P(CellTestP, select_inside_rcut_corners) {
    int nbar_total = 0;
    for (int irep=0; irep < NREP; irep++) {
        // Test parameters:
        // - Random cell
        Cell* cell = create_random_cell(2*irep, 1.0, true);
        // - Increasing rcut
        double rcut = (irep+1)*0.1;
        // - Random center
        double center[3];
        fill_random_double(47332+irep, center, 3, -2.0, 2.0);
        // - Alternating values for shape and pbc
        int shape[nvec];
        bool pbc[nvec];
        for (int ivec=0; ivec < nvec; ivec++) {
            shape[ivec] = ((irep*(ivec+1)) % 5) + 1;
            pbc[ivec] = true;
        }

        // Compute the bars.
        int nbar1 = cell->select_inside_rcut(center, rcut, shape, pbc, NULL);
        //if (nbar1 > 1000) {delete cell; continue;}
        int bars[(nvec+2)*nbar1];
        int nbar2 = cell->select_inside_rcut(center, rcut, shape, pbc, bars);
        EXPECT_EQ(nbar1, nbar2);
        nbar_total += nbar2;

        // Test if the corners of each bar fall outside of the sphere
        for (int ibar=0; ibar < nbar2; ibar++) {
            int* bar = bars + (nvec+1)*ibar;

            double frac_corner[3] = {0, 0, 0};
            double cart_corner[3] = {0, 0, 0};
            double dist;

            // loop of begin and end of the bar (last two integers in the bar).
            cell->to_rfrac(center, frac_corner);
            for (int ilast=0; ilast < 2; ilast++) {
                frac_corner[nvec-1] = bar[nvec-ilast];
                if (nvec==1) {
                    cell->to_rcart(frac_corner, cart_corner);
                    dist = vec3::distance(cart_corner, center);
                    EXPECT_GT(dist, rcut);
                } else if (nvec==2) {
                    //
                    frac_corner[0] = bar[0];
                    cell->to_rcart(frac_corner, cart_corner);
                    dist = vec3::distance(cart_corner, center);
                    EXPECT_GT(dist, rcut);
                    //
                    frac_corner[0] = bar[0]+1;
                    cell->to_rcart(frac_corner, cart_corner);
                    dist = vec3::distance(cart_corner, center);
                    EXPECT_GT(dist, rcut);
                } else if (nvec==3) {
                    //
                    frac_corner[0] = bar[0];
                    frac_corner[1] = bar[1];
                    cell->to_rcart(frac_corner, cart_corner);
                    dist = vec3::distance(cart_corner, center);
                    EXPECT_GT(dist, rcut);
                    //
                    frac_corner[0] = bar[0]+1;
                    frac_corner[1] = bar[1];
                    cell->to_rcart(frac_corner, cart_corner);
                    dist = vec3::distance(cart_corner, center);
                    EXPECT_GT(dist, rcut);
                    //
                    frac_corner[0] = bar[0];
                    frac_corner[1] = bar[1]+1;
                    cell->to_rcart(frac_corner, cart_corner);
                    dist = vec3::distance(cart_corner, center);
                    EXPECT_GT(dist, rcut);
                    //
                    frac_corner[0] = bar[0]+1;
                    frac_corner[1] = bar[1]+1;
                    cell->to_rcart(frac_corner, cart_corner);
                    dist = vec3::distance(cart_corner, center);
                    EXPECT_GT(dist, rcut);
                }
            }
        }
    }
    // Sufficiency check
    EXPECT_LE(NREP*((nvec-1)*3 + 1), nbar_total);
}

// Instantiation of parameterized tests
// ------------------------------------

INSTANTIATE_TEST_CASE_P(CellTest123, CellTestP, ::testing::Range(1, 4));


// smart_wrap
// ----------

TEST(SmartWrap, examples) {
    EXPECT_EQ(0, smart_wrap(-15, 5, true));
    EXPECT_EQ(0, smart_wrap(-5, 5, true));
    EXPECT_EQ(2, smart_wrap(-3, 5, true));
    EXPECT_EQ(4, smart_wrap(-1, 5, true));
    EXPECT_EQ(0, smart_wrap(0, 5, true));
    EXPECT_EQ(3, smart_wrap(3, 5, true));
    EXPECT_EQ(0, smart_wrap(5, 5, true));
    EXPECT_EQ(1, smart_wrap(6, 5, true));
    EXPECT_EQ(0, smart_wrap(10, 5, true));
    EXPECT_EQ(2, smart_wrap(12, 5, true));
    EXPECT_EQ(-1, smart_wrap(-15, 5, false));
    EXPECT_EQ(-1, smart_wrap(-5, 5, false));
    EXPECT_EQ(-1, smart_wrap(-3, 5, false));
    EXPECT_EQ(-1, smart_wrap(-1, 5, false));
    EXPECT_EQ(0, smart_wrap(0, 5, false));
    EXPECT_EQ(3, smart_wrap(3, 5, false));
    EXPECT_EQ(4, smart_wrap(4, 5, false));
    EXPECT_EQ(-1, smart_wrap(5, 5, false));
    EXPECT_EQ(-1, smart_wrap(6, 5, false));
    EXPECT_EQ(-1, smart_wrap(10, 5, false));
    EXPECT_EQ(-1, smart_wrap(12, 5, false));
}
