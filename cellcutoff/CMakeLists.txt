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

# Configure a header file to pass some of the CMake settings
# to the source code.
configure_file (
  ${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

# Explain that the current directory has to be in the source path/
# This allows use to use include "cellcutoff/foo.h", which is consistent
# with the way external programs would compile against cellcutoff.
include_directories(${PROJECT_SOURCE_DIR})

# Define source files
set(SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/cell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/decomposition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/iterators.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/sphere_slice.cpp
)

# Define header files
set(HEADER_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/cell.h
  ${CMAKE_CURRENT_SOURCE_DIR}/decomposition.h
  ${CMAKE_CURRENT_SOURCE_DIR}/iterators.h
  ${CMAKE_CURRENT_SOURCE_DIR}/sphere_slice.h
  ${CMAKE_CURRENT_SOURCE_DIR}/vec3.h
  ${CMAKE_CURRENT_BINARY_DIR}/config.h
)

# Define the shared library
add_library(cellcutoff SHARED ${SOURCE_FILES})
set_property(TARGET cellcutoff PROPERTY VERSION ${CELLCUTOFF_VERSION})
set_property(TARGET cellcutoff PROPERTY SOVERSION ${CELLCUTOFF_SOVERSION})

# Install include files
install (FILES ${HEADER_FILES} DESTINATION include/cellcutoff)

# Install library
install (TARGETS cellcutoff DESTINATION lib)

# Unit Tests with Google testing framework
add_subdirectory(tests)
