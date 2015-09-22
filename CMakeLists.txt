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

cmake_minimum_required(VERSION 3.0)
project(CellLists CXX)

# Version information
set(CELLLISTS_MAJOR 0)
set(CELLLISTS_MINOR 0)
set(CELLLISTS_PATCH 0)

set(CELLLISTS_VERSION ${CELLLISTS_MAJOR}.${CELLLISTS_MINOR}.${CELLLISTS_PATCH})
set(CELLLISTS_SOVERSION ${CELLLISTS_MAJOR}.${CELLLISTS_MINOR})
set(PROJECT_BRIEF "A 3D domain decomposition library")

# Configure the compiler
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wconversion -Wextra -Werror -pedantic")

# Additional compiler settings when doing a debug build
string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE)
if(CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    #SET (CMAKE_VERBOSE_MAKEFILE "ON")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        # When coding (debug mode), we want to know when there might be trouble ...
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_CXX_FLAGS} -coverage")
        set(CELLLISTS_DO_COVERAGE "ON")
    endif()
endif()


# Compilation of celllists according to celllists/CMakeLists.txt
add_subdirectory(celllists)

# Documentation with doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${PROJECT_SOURCE_DIR}/doxygen.conf.in ${PROJECT_BINARY_DIR}/doxygen.conf @ONLY)
    add_custom_target(doc ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/doxygen.conf
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)
