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


/* Usage tests in this file combine different features from cellcutoff, comparable to how
   this library is to be used by others. Strictly speaking, these aren't unit tests, but
   we want these tests to pass as well! They check non-trivial outcomes that are more a
   validation of the algorithm choices, rather than algorithm implementations. */


#include <cmath>

#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

#include <cellcutoff/cell.h>
#include <cellcutoff/decomposition.h>
#include <cellcutoff/iterators.h>
#include <cellcutoff/vec3.h>

#include "cellcutoff/tests/common.h"


namespace cl = cellcutoff;
namespace vec3 = cellcutoff::vec3;


class UsageTestP : public ::testing::TestWithParam<int> {
 public:
  virtual void SetUp() {
    nvec = GetParam();
  }

  std::unique_ptr<cl::Cell> create_random_cell(const unsigned int seed,
      const double scale = 1.0, const double ratio = 0.1, const bool cuboid = false) {
    return std::unique_ptr<cl::Cell>(cl::create_random_cell(seed, nvec, scale, ratio, cuboid));
  }

  int nvec;
};


TEST_P(UsageTestP, points_within_cutoff) {
  size_t npoint_total = 0;
  for (int irep = 0; irep < NREP; ++irep) {
    // Select a random center and a cutoff
    double cutoff = 1.0 + static_cast<double>(irep)/NREP;
    double center[3];
    unsigned int seed  = fill_random_double(irep, center, 3, -cutoff, cutoff);

    // Generate a set of random points (in and beyond cutoff)
    std::vector<cl::Point> points;
    for (int ipoint = 0; ipoint < NPOINT; ++ipoint) {
      double cart[3];
      seed = fill_random_double(seed, cart, 3, -cutoff, cutoff);
      points.push_back(cl::Point(cart));
    }

    // Make a reasonable cell and subcell
    std::unique_ptr<cl::Cell> cell(create_random_cell(seed, cutoff, 0.6, false));
    cell->iwrap_box(center);
    int shape[3] = {-1, -1, -1};
    std::unique_ptr<cl::Cell> subcell(cell->create_subcell(cutoff*0.2, shape));
    EXPECT_LT(-1*(nvec <= 0), shape[0]);
    EXPECT_LT(-1*(nvec <= 1), shape[1]);
    EXPECT_LT(-1*(nvec <= 2), shape[2]);

    // Assign icells to points, sort and make a cell_map.
    cl::assign_icell(*subcell, shape, points.data(), points.size(), sizeof(cl::Point));
    std::sort(points.begin(), points.end());
    std::unique_ptr<cl::CellMap>
      cell_map(cl::create_cell_map(points.data(), points.size(), sizeof(cl::Point)));

    // Compute the points within the cutoff in the most efficient way, i.e. using
    // the DeltaIterator.
    std::vector<size_t> ipoints_bars;
    std::set<std::array<double, 3>> deltas;
    for (cl::DeltaIterator dit(*subcell, shape, center, cutoff, points.data(), points.size(),
         sizeof(cl::Point), *cell_map); dit.busy(); ++dit) {
      ++npoint_total;
      ipoints_bars.push_back(dit.ipoint());
      EXPECT_NEAR(dit.distance(), vec3::norm(dit.delta()), EPS);
      // No duplicates allowed!
      std::array<double, 3> delta;
      std::copy(dit.delta(), dit.delta() + 3, delta.data());
      EXPECT_EQ(0, deltas.count(delta));
      deltas.insert(delta);
    }
    std::sort(ipoints_bars.begin(), ipoints_bars.end());

    // Compute the points within the cutoff in a less efficient way, i.e. using
    // Cell::ranges_cutoff. Results should match.
    int ranges_begin[3];
    int ranges_end[3];
    size_t ncell = subcell->ranges_cutoff(center, cutoff, ranges_begin, ranges_end);
    EXPECT_LT(0, ncell);
    std::vector<size_t> ipoints_ranges;
    for (int icell0 = ranges_begin[0]; icell0 < ranges_end[0]; ++icell0) {
      for (int icell1 = ranges_begin[1]; icell1 < ranges_end[1]; ++icell1) {
        for (int icell2 = ranges_begin[2]; icell2 < ranges_end[2]; ++icell2) {
          std::array<int, 3> icell;
          int coeffs[3];
          icell[0] = cl::robust_wrap(icell0, shape[0], &coeffs[0]);
          icell[1] = cl::robust_wrap(icell1, shape[1], &coeffs[1]);
          icell[2] = cl::robust_wrap(icell2, shape[2], &coeffs[2]);
          auto it = cell_map->find(icell);
          if (it != cell_map->end()) {
            for (size_t ipoint = it->second[0]; ipoint < it->second[1]; ++ipoint) {
              double cart[3];
              std::copy(points[ipoint].cart_, points[ipoint].cart_ + 3, cart);
              cell->iadd_vec(cart, coeffs);
              double d = vec3::distance(center, cart);
              if (d < cutoff) {
                ipoints_ranges.push_back(ipoint);
              }
            }
          }
        }
      }
    }
    std::sort(ipoints_ranges.begin(), ipoints_ranges.end());

    // Compare the two so far;
    EXPECT_EQ(ipoints_bars.size(), ipoints_ranges.size());
    for (size_t i = 0; i < ipoints_bars.size(); ++i) {
      EXPECT_EQ(ipoints_bars.at(i), ipoints_ranges.at(i));
    }

    // If aperiodic, compute the points within the cutoff in a dumb way: just try them
    // all! Results should match.
    if (nvec == 0) {
      std::vector<size_t> ipoints_dumb;
      for (size_t ipoint = 0; ipoint < points.size(); ++ipoint) {
        double d = vec3::distance(center, points[ipoint].cart_);
        if (d < cutoff) {
          ipoints_dumb.push_back(ipoint);
        }
      }

      // Compare
      EXPECT_EQ(ipoints_bars.size(), ipoints_dumb.size());
      for (size_t i = 0; i < ipoints_bars.size(); ++i) {
        EXPECT_EQ(ipoints_bars.at(i), ipoints_dumb.at(i));
      }
    }
  }
  // Sufficiency check
  EXPECT_LE((NREP*NPOINT)/10, npoint_total);
}


// Instantiation of parameterized tests
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

INSTANTIATE_TEST_CASE_P(UsageTest0123, UsageTestP, ::testing::Range(0, 4));

// vim: textwidth=90 et ts=2 sw=2
