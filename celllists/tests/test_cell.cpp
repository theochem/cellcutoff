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
#include "celllists/cell.h"

// Helper functions

Cell create_random_cell(long nvec, unsigned int seed) {
    srand(seed);
    double rvecs[nvec*3];
    while (true) {
        try {
            for (long ivec=0; ivec<nvec*3; ivec++) {
                rvecs[ivec] = rand()/(RAND_MAX + 1.0);
            }
            return Cell(rvecs, nvec);
        }
        catch (singular_cell_vectors) {}
    }
}

// Actual tests

TEST(cell_test, constructor_singular1) {
    double rvecs[3] = {0, 0, 0};
    ASSERT_THROW(Cell cell(rvecs, 1), singular_cell_vectors);
}

TEST(cell_test, constructor_singular2) {
    double rvecs[6] = {1, 0, 0, 1, 0, 0};
    ASSERT_THROW(Cell cell(rvecs, 2), singular_cell_vectors);
}

TEST(cell_test, constructor_singular3) {
    double rvecs[9] = {1, 0, 0, 0, 1, 0, 0.5, 0.5, 0};
    ASSERT_THROW(Cell cell(rvecs, 3), singular_cell_vectors);
}

TEST(cell_test, constructor_nvec_negative) {
    double rvecs[9] = {1, 0, 1, 0, 1, 0, 0.5, 0.5, 0};
    ASSERT_THROW(Cell cell(rvecs, -1), std::domain_error);
}

TEST(cell_test, constructor_nvec_too_large) {
    double rvecs[9] = {1, 0, 1, 0, 1, 0, 0.5, 0.5, 0};
    ASSERT_THROW(Cell cell(rvecs, 4), std::domain_error);
}

TEST(cell_test, constructor_simple) {
    double rvecs[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    for(int nvec=0; nvec <=3; nvec++) {
        Cell cell(rvecs, nvec);
        SCOPED_TRACE(nvec);
        ASSERT_EQ(cell.get_nvec(), nvec);
        ASSERT_EQ(cell.get_rvec(0, 0), 1.0);
        ASSERT_EQ(cell.get_rvec(0, 1), 0.0);
        ASSERT_EQ(cell.get_rvec(0, 1), 0.0);
        ASSERT_EQ(cell.get_rvec(1, 0), 0.0);
        ASSERT_EQ(cell.get_rvec(1, 1), 1.0);
        ASSERT_EQ(cell.get_rvec(1, 2), 0.0);
        ASSERT_EQ(cell.get_rvec(2, 0), 0.0);
        ASSERT_EQ(cell.get_rvec(2, 1), 0.0);
        ASSERT_EQ(cell.get_rvec(2, 2), 1.0);
        ASSERT_EQ(cell.get_gvec(0, 0), 1.0);
        ASSERT_EQ(cell.get_gvec(0, 1), 0.0);
        ASSERT_EQ(cell.get_gvec(0, 2), 0.0);
        ASSERT_EQ(cell.get_gvec(1, 0), 0.0);
        ASSERT_EQ(cell.get_gvec(1, 1), 1.0);
        ASSERT_EQ(cell.get_gvec(1, 2), 0.0);
        ASSERT_EQ(cell.get_gvec(2, 0), 0.0);
        ASSERT_EQ(cell.get_gvec(2, 1), 0.0);
        ASSERT_EQ(cell.get_gvec(2, 2), 1.0);
        ASSERT_EQ(cell.get_volume(), nvec > 0);
        ASSERT_EQ(cell.get_rspacing(0), 1.0);
        ASSERT_EQ(cell.get_rspacing(1), 1.0);
        ASSERT_EQ(cell.get_rspacing(2), 1.0);
        ASSERT_EQ(cell.get_gspacing(0), 1.0);
        ASSERT_EQ(cell.get_gspacing(1), 1.0);
        ASSERT_EQ(cell.get_gspacing(2), 1.0);
        ASSERT_EQ(cell.is_cubic(), true);
        ASSERT_EQ(cell.is_cuboid(), true);
    }
}

TEST(cell_test, wrap_edges) {
    double rvecs[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    Cell cell(rvecs, 3);
    double delta[3] = {-0.5, -0.5, -0.5};
    cell.wrap(delta);
    ASSERT_EQ(delta[0], 0.5);
    ASSERT_EQ(delta[1], 0.5);
    ASSERT_EQ(delta[2], 0.5);
}

TEST(cell_test, wrap_random) {
    for (int irep=0; irep < 100; irep++) {
        Cell cell = create_random_cell(3, irep);
        double delta[3] = {0.5*irep, 0.9*irep, 1.3*irep};
        double frac[3];
        cell.wrap(delta);
        cell.to_frac(delta, frac);
        ASSERT_LT(frac[0], 0.5);
        ASSERT_LT(frac[1], 0.5);
        ASSERT_LT(frac[2], 0.5);
        ASSERT_GE(frac[0], -0.5);
        ASSERT_GE(frac[1], -0.5);
        ASSERT_GE(frac[2], -0.5);
    }
}
