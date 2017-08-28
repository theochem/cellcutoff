# -*- coding: utf-8 -*-
# CellCutoff is a library for periodic boundary conditions and real-space cutoff calculations.
# Copyright (C) 2017 The CellCutoff Development Team
#
# This file is part of CellCutoff.
#
# CellCutoff is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# CellCutoff is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --

cmake_minimum_required(VERSION 3.1)
project(CellCutoff CXX)

# Get the version from an include file:
include(CMakeListsVersion.txt.in)
set(CELLCUTOFF_VERSION ${GIT_TAG_VERSION})
set(CELLCUTOFF_SOVERSION ${GIT_TAG_SOVERSION})
set(PROJECT_BRIEF "Library for periodic boundary conditions and real-space cutoff calculations")

# The googletest snippets were taken from:
#    http://crascit.com/2015/07/25/cmake-gtest/
# and modified to work with googletest-1.8.0 instead of the latest git revision.

# Download and unpack googletest at configure time
configure_file(CMakeListsGTest.txt.in
               googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)
execute_process(COMMAND ${CMAKE_COMMAND} --build .
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/googletest-download)

# Add googletest directly to our build. This adds the following targets: gtest,
# gtest_main, gmock and gmock_main
add_subdirectory(${CMAKE_BINARY_DIR}/googletest-src
                 ${CMAKE_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL )

# NOTE: A change request has been submitted to the GoogleTest
#       project to remove the need for the following lines.
#       They add the gtest and gmock include directories as
#       target dependencies, since they are not added by
#       googletest (yet).
target_include_directories(gtest INTERFACE
    "${CMAKE_BINARY_DIR}/googletest-src/include")
target_include_directories(gtest_main INTERFACE
    "${CMAKE_BINARY_DIR}/googletest-src/include")

# # Only switch this one when done with googletest.
# set(CMAKE_VERBOSE_MAKEFILE "ON")

# General compiler flags
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wconversion -Wextra")

# Additional compiler settings when doing a debug build
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # When coding (debug mode), we want to know when there might be trouble ...
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_CXX_FLAGS} -coverage")
        set(CELLCUTOFF_DO_COVERAGE "ON")
        # Google test uses comparisons between signed and unsigned ...
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")
    endif()
endif()

# Compilation of cellcutoff according to cellcutoff/CMakeLists.txt
add_subdirectory(cellcutoff)

# Documentation with doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${PROJECT_SOURCE_DIR}/doxygen.conf.in ${PROJECT_BINARY_DIR}/doxygen.conf @ONLY)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doxygen.conf
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

# Source package with CPack
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME
  "cellcutoff-${CELLCUTOFF_VERSION}")
set(CPACK_SOURCE_IGNORE_FILES
  "/build/;/release/;/debug/;/.git/;/python-cellcutoff/;/HEADER/;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)
add_custom_target(sdist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
