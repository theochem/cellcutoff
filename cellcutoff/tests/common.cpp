// CellCutoff is a library for periodic boundary conditions and real-space cutoff calculations.
// Copyright (C) 2017 The CellCutoff Development Team
//
// This file is part of CellCutoff.
//
// CellCutoff is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// CellCutoff is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
// --


#include "cellcutoff/tests/common.h"

#include <cmath>

#include <algorithm>
#include <limits>
#include <random>
#include <stdexcept>

#include <cellcutoff/vec3.h>


namespace cl = cellcutoff;
namespace vec3 = cellcutoff::vec3;


//! Internal helper that just makes a useful random seed for a random generates.
unsigned int get_next_seed(std::minstd_rand gen) {
  std::uniform_int_distribution<unsigned int>
      dis_seed(0, std::numeric_limits<unsigned int>::max());
  return dis_seed(gen);
}


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


unsigned int random_point(unsigned int seed,  const double* center,
    const double cutoff, double* point, double* norm) {
  seed = fill_random_double(seed, point, 3, -cutoff, cutoff);
  *norm = vec3::norm(point);
  vec3::iadd(point, center);
  return seed;
}


// Fixtures
// --------

void CellTest::set_up_data() {
  // Example tests
  std::fill(myvecs, myvecs + 9, 0.0);
  myvecs[0] = 2;
  myvecs[4] = 1;
  myvecs[8] = 4;
  if (nvec == 2) {
    myvecs[4] = 0.0;
    myvecs[5] = 4.0;
  }
  mycell.reset(new cl::Cell(myvecs, nvec));
  // Singular cell vectors
  std::fill(singvecs, singvecs + 9, 0.0);
  if (nvec > 1) {
    singvecs[0] = 1.0;
    singvecs[3] = 0.5;
  }
  if (nvec == 3) {
    singvecs[3] = 0.0;
    singvecs[4] = 2.0;
    singvecs[6] = 0.5;
    singvecs[7] = 0.8;
  }
}


std::unique_ptr<cl::Cell> CellTest::create_random_cell(const unsigned int seed,
    const double scale, const double ratio, const bool cuboid) {
  return std::unique_ptr<cl::Cell>(cl::create_random_cell(seed, nvec, scale, ratio, cuboid));
}


TEST(CommonTest, domain) {
  EXPECT_THROW(fill_random_double(0, nullptr, 0, 0.0, 1.0), std::domain_error);
  EXPECT_THROW(fill_random_double(0, nullptr, -1, 0.0, 1.0), std::domain_error);
  EXPECT_THROW(fill_random_int(0, nullptr, 0, 0, 1), std::domain_error);
  EXPECT_THROW(fill_random_int(0, nullptr, -1, 0, 1), std::domain_error);
  EXPECT_THROW(fill_random_int(0, nullptr, 1, 1, 0), std::domain_error);
  EXPECT_THROW(fill_random_permutation(0, nullptr, 0), std::domain_error);
  EXPECT_THROW(fill_random_permutation(0, nullptr, -1), std::domain_error);
}

// vim: textwidth=90 et ts=2 sw=2
