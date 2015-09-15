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
// aint with this program; if not, see <http://www.gnu.org/licenses/>
//
//--


#include "common.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

#include <gtest/gtest.h>

#include <celllists/vec3.h>


namespace cl = celllists;
namespace vec3 = celllists::vec3;


unsigned int get_next_seed(std::minstd_rand gen) {
  std::uniform_int_distribution<unsigned int>
    dis_seed(0, std::numeric_limits<unsigned int>::max());
  return dis_seed(gen);
}


//! Fills an array of doubles with random numbers in range ]-0.5*scale, 0.5*scale]
unsigned int fill_random_double(const unsigned int seed, double* array, const int size,
  const double low, const double high) {

  // Parameter check
  if (size <= 0)
    throw std::domain_error("Array size must be strictly positive.");

  // Fill the array with random data for given seed.
  std::minstd_rand gen(seed);
  std::uniform_real_distribution<double> dis(low, high);
  for (int i=0; i < size; ++i)
    array[i] = dis(gen);

  // Generate a different seed for the next call
  return get_next_seed(gen);
}

//! Fills an array of int with random numbers in range [-range, range]
unsigned int fill_random_int(const unsigned int seed, int* array, const int size,
  const int begin, const int end) {
  // Parameter check
  if (size <= 0)
    throw std::domain_error("Array size must be strictly positive.");
  if (begin > end)
    throw std::domain_error("Begin cannot be larger than end.");

  // Fill the array with random data for given seed.
  std::minstd_rand gen(seed);
  std::uniform_int_distribution<int> dis(begin, end);
  for (int i=0; i < size; ++i)
    array[i] = dis(gen);

  // Generate a different seed for the next call
  return get_next_seed(gen);
}

//! Fills and array of int with a random permutation
unsigned int fill_random_permutation(const unsigned int seed, int* array,
  const int size) {
  // Parameter check
  if (size <= 0)
    throw std::domain_error("Array size must be strictly positive.");

  // Fill the array with integers in order.
  for (int i=0; i < size; ++i)
    array[i] = static_cast<int>(i);

  // Make a random permutation
  std::minstd_rand gen(seed);
  std::shuffle(array, array + size, gen);

  // Generate a different seed for the next call
  return get_next_seed(gen);
}

//! Random cell with a volume larger than (0.1*scale)**nvec
std::unique_ptr<cl::Cell> create_random_cell_nvec(unsigned int seed, const int nvec,
  const double scale, const bool cuboid) {
  // Range check
  if ((nvec <= 0) || (nvec > 3)) {
    throw std::domain_error("A random cell must be 1D, 2D or 2D periodic.");
  }
  // Randomly construct a cell till a decent one (sufficient volume) is found.
  double rvecs[9];
  while (true) {
    seed = fill_random_double(seed, rvecs, 9, -scale, scale);
    if (cuboid) {
      rvecs[1] = 0.0;
      rvecs[2] = 0.0;
      if (nvec > 1) {
        rvecs[3] = 0.0;
        rvecs[5] = 0.0;
      }
      if (nvec > 2) {
        rvecs[6] = 0.0;
        rvecs[7] = 0.0;
      }
    }
    try {
      std::unique_ptr<cl::Cell> cell(new cl::Cell(rvecs, nvec));
      if (cell->get_volume() > pow(0.1*scale, nvec))
        return cell;
    } catch (cl::singular_cell_vectors) { }
  }
}

//! Compute a random point in a cubic box centered around center. Also computes distance.
unsigned int random_point(unsigned int seed,  const double* center,
  const double rcut, double* point, double* norm) {
  seed = fill_random_double(seed, point, 3, -rcut, rcut);
  *norm = vec3::norm(point);
  vec3::iadd(point, center);
  return seed;
}


TEST(CommonTest, domain) {
  EXPECT_THROW(fill_random_double(0, nullptr, 0, 0.0, 1.0), std::domain_error);
  EXPECT_THROW(fill_random_double(0, nullptr, -1, 0.0, 1.0), std::domain_error);
  EXPECT_THROW(fill_random_int(0, nullptr, 0, 0, 1), std::domain_error);
  EXPECT_THROW(fill_random_int(0, nullptr, -1, 0, 1), std::domain_error);
  EXPECT_THROW(fill_random_int(0, nullptr, 1, 1, 0), std::domain_error);
  EXPECT_THROW(fill_random_permutation(0, nullptr, 0), std::domain_error);
  EXPECT_THROW(fill_random_permutation(0, nullptr, -1), std::domain_error);
  EXPECT_THROW(create_random_cell_nvec(-1, 0, 1, false), std::domain_error);
  EXPECT_THROW(create_random_cell_nvec(0, 0, 1, false), std::domain_error);
  EXPECT_THROW(create_random_cell_nvec(4, 0, 1, false), std::domain_error);
}
