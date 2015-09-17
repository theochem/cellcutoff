# CellList is a 3D domain decomposition library.
# Copyright (C) 2011-2015 The CellList Development Team
#
# This file is part of CellList.
#
# CellList is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# CellList is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--

# Make use of the GTest CMake module
find_package(GTest REQUIRED)

# Make sure the header files are found
include_directories(${GTEST_INCLUDE_DIRS})

# Define test source files
set(TEST_SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/test_cell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_decomposition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/test_sphere_slice.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/common.cpp
)

# The testing executable
add_executable(test_celllists EXCLUDE_FROM_ALL ${TEST_SOURCE_FILES})
set_property(TARGET test_celllists PROPERTY CXX_STANDARD 11)

# How the testing executable has to be linked...
target_link_libraries(test_celllists celllists ${GTEST_BOTH_LIBRARIES} gcov pthread)

# Use our own test target because ctest hides relevant output.
IF(CELLLISTS_DO_COVERAGE)
    add_custom_target(clean_coverage
                      COMMAND ${PROJECT_SOURCE_DIR}/tools/coverage_clean.sh
                      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                      COMMENT "Removing old coverage statistics")

    add_custom_target(check
                      COMMAND test_celllists
                      DEPENDS test_celllists clean_coverage
                      COMMENT "Running tests")

    add_custom_target(coverage
                      COMMAND ${PROJECT_SOURCE_DIR}/tools/coverage_analyze.py
                      DEPENDS check
                      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
                      COMMENT "Coverage analysis")

ELSE()
    add_custom_target(check
                      COMMAND test_celllists
                      DEPENDS test_celllists
                      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      COMMENT "Running tests")


ENDIF()
