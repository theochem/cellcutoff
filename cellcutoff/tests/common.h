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


#ifndef CELLCUTOFF_TESTS_COMMON_H_
#define CELLCUTOFF_TESTS_COMMON_H_

#include <memory>

#include <gtest/gtest.h>

#include "cellcutoff/cell.h"


namespace cl = cellcutoff;


// Number of repetitions on a test with random data
#define NREP 100
// Number of random Cartesian points to sample in a test
#define NPOINT 1000
// Allowed tolerance when comparing doubles (EPS stands for epsilon.)
#define EPS 1e-10


/* Some notes:

    - Usually, the output argument is last. The exception in this module is the size
      of the output argument, which comes after the actual output argument.

*/


//! Fills an array of doubles with random numbers in range ]-0.5*scale, 0.5*scale]
unsigned int fill_random_double(const unsigned int seed, double* array, const int size,
    const double low = -0.5, const double high = 0.5);

//! Fills an array of int with random numbers in range [-range, range]
unsigned int fill_random_int(const unsigned int seed, int* array, const int size,
    const int begin, const int end);

//! Fills and array of int with a random permutation
unsigned int fill_random_permutation(const unsigned int seed, int* array, const int size);

//! Compute a random point in a cubic box centered around center. Also computes distance.
unsigned int random_point(const unsigned int seed,  const double* center,
    const double cutoff, double* point, double* norm);


// Fixtures
// ========

//! Abstract base class for CellTest fixtures.
class CellTest : public ::testing::Test {
 public:
  //! Fix number of cell vectors, should be defined in subclass.
  virtual void SetUp() = 0;
  //! Define some example arrays to work with in the tests.
  void set_up_data();
  //! Create a random cell with the right number of cell vectors.
  std::unique_ptr<cl::Cell> create_random_cell(const unsigned int seed,
      const double scale = 1.0, const double ratio = 0.1, const bool cuboid = false);

  //! Number of cell vectors (fixed by setup).
  int nvec;
  //! Cell vectors for a trivial cuboid cell.
  double myvecs[9];
  //! The cell object derived from myvecs.
  std::unique_ptr<cl::Cell> mycell;
  //! A singular set of cell vectors.
  double singvecs[9];
};


//! Parameterized CellTest fixtures. (Tests get repeated for nvec = 0, 1, 2 and 3.)
class CellTestP : public CellTest, public ::testing::WithParamInterface<int> {
 public:
  //! Fix the number of cell vectors (paramerized).
  virtual void SetUp() {
    nvec = GetParam();
    set_up_data();
  }
};


//! CellTest fixture with nvec = 0.
class CellTest0 : public CellTest {
 public:
  //! Set the number of cell vectors to zero.
  virtual void SetUp() {
    nvec = 0;
    set_up_data();
  }
};


//! CellTest fixture with nvec = 1.
class CellTest1 : public CellTest {
 public:
  //! Set the number of cell vectors to one.
  virtual void SetUp() {
    nvec = 1;
    set_up_data();
  }
};


//! CellTest fixture with nvec = 2.
class CellTest2 : public CellTest {
 public:
  //! Set the number of cell vectors to two.
  virtual void SetUp() {
    nvec = 2;
    set_up_data();
  }
};


//! CellTest fixture with nvec = 3.
class CellTest3 : public CellTest {
 public:
  //! Set the number of cell vectors to three.
  virtual void SetUp() {
    nvec = 3;
    set_up_data();
  }
};


#endif  // CELLCUTOFF_TESTS_COMMON_H_

// vim: textwidth=90 et ts=2 sw=2
