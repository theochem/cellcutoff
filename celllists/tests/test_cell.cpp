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
#include "celllists/cell.h"

TEST(general_cell_test, constructor_singular1) {
    double rvecs[3] = {0, 0, 0};
    ASSERT_THROW(GeneralCell cell(rvecs, 1), std::domain_error);
}

TEST(general_cell_test, constructor_singular2) {
    double rvecs[6] = {1, 0, 0, 1, 0, 0};
    ASSERT_THROW(GeneralCell cell(rvecs, 2), std::domain_error);
}

TEST(general_cell_test, constructor_singular3) {
    double rvecs[9] = {1, 0, 0, 0, 1, 0, 0.5, 0.5, 0};
    ASSERT_THROW(GeneralCell cell(rvecs, 3), std::domain_error);
}

TEST(general_cell_test, constructor_nvec_negative) {
    double rvecs[9] = {1, 0, 1, 0, 1, 0, 0.5, 0.5, 0};
    ASSERT_THROW(GeneralCell cell(rvecs, -1), std::domain_error);
}

TEST(general_cell_test, constructor_nvec_too_large) {
    double rvecs[9] = {1, 0, 1, 0, 1, 0, 0.5, 0.5, 0};
    ASSERT_THROW(GeneralCell cell(rvecs, 4), std::domain_error);
}

TEST(general_cell_test, volume) {
    double rvecs[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    GeneralCell cell(rvecs, 3);
    EXPECT_EQ(cell.get_volume(), 1.0);
}
