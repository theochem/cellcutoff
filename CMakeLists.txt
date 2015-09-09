# -*- coding: utf-8 -*-
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
set (CELLLISTS_MAJOR 0)
set (CELLLISTS_MINOR 0)
set (CELLLISTS_PATCH 0)

set (CELLLISTS_VERSION ${CELLLISTS_MAJOR}.${CELLLISTS_MINOR}.${CELLLISTS_PATCH})
set (CELLLISTS_SOVERSION ${CELLLISTS_MAJOR}.${CELLLISTS_MINOR})

# Configure a header file to pass some of the CMake settings
# to the source code.
configure_file (
  ${PROJECT_SOURCE_DIR}/celllists/config.h.in
  ${PROJECT_BINARY_DIR}/celllists/config.h
)

# Explain that the current directory has to be in the source path/
# This allows use to use include "celllists/foo.h", which is consistent
# with the way external programs would compile against celllists.
include_directories(${PROJECT_SOURCE_DIR})

# Define source files
set(SOURCE_FILES
  ${PROJECT_SOURCE_DIR}/celllists/cell.cpp
)

# Define header files
set(HEADER_FILES
  ${PROJECT_SOURCE_DIR}/celllists/cell.h
  ${PROJECT_BINARY_DIR}/celllists/config.h
)

# Define the shared library
add_library(celllists SHARED ${SOURCE_FILES})
set_property(TARGET celllists PROPERTY VERSION ${CELLLISTS_VERSION})
set_property(TARGET celllists PROPERTY SOVERSION ${CELLLISTS_SOVERSION})

# Install include files
install (FILES ${HEADER_FILES} DESTINATION include/cellists)

# Install library
install (TARGETS celllists DESTINATION lib)

# Unit Tests with Google testing framework
enable_testing()
add_subdirectory(celllists/tests)
