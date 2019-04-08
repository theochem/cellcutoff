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

# General compiler flags
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wconversion -Wextra")

# Additional compiler settings when doing a debug build
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    # SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
    # SET(CMAKE_INSTALL_RPATH $ENV{LIBRARY_PATH})
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
