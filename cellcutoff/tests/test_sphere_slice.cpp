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


#include <cmath>
#include <cstdlib>

#include <memory>
#include <stdexcept>
#include <utility>

#include <gtest/gtest.h>

#include <cellcutoff/sphere_slice.h>
#include <cellcutoff/vec3.h>

#include "cellcutoff/tests/common.h"


namespace cl = cellcutoff;
namespace vec3 = cellcutoff::vec3;


class SphereSliceTest : public ::testing::Test {
 public:
  void SetUp() {
    my_center[0] = 0.4;
    my_center[1] = -2.0;
    my_center[2] = 1.0;
    std::fill(easy_normals, easy_normals + 9, 0.0);
    easy_normals[0] = 1.0;
    easy_normals[4] = 1.0;
    easy_normals[8] = 1.0;
  }

  std::unique_ptr<cl::SphereSlice> create_random_problem(unsigned int seed,
      const double radius, double* center, double* normals) {
    // get random double and keep seed of more random data is needed
    seed = fill_random_double(seed, center, 3);
    double vol;
    do {
      seed = fill_random_double(seed, normals, 9);
      vol = fabs(vec3::triple(normals, normals + 3, normals + 6));
    } while (vol < 0.001);
    return std::unique_ptr<cl::SphereSlice>(new cl::SphereSlice(center, normals,
      radius));
  }

  void random_cut(const unsigned int seed, const cl::SphereSlice &slice,
      const int id_cut, double *cut, double *cut_min, double *cut_max) {
    // Do a solve_full_low, to know over which range we can cut
    slice.solve_full_low(id_cut, cut_min, cut_max, nullptr, nullptr);
    fill_random_double(seed, cut, 1, *cut_min, *cut_max);
  }

  void random_slice(const unsigned int seed, const cl::SphereSlice &slice,
      const int id_cut, double *cut_begin, double *cut_end,
    double *cut_min, double *cut_max) {
    // Do a solve_full_low, to know over which range we can make a
    // disc-like slice that still intersects with the sphere.
    slice.solve_full_low(id_cut, cut_min, cut_max, nullptr, nullptr);
    EXPECT_LT(*cut_min, *cut_max);
    // Select two cut positions between cut_min and cut_max
    double cuts[2];
    fill_random_double(seed, cuts, 2, *cut_min, *cut_max);
    if (cuts[1] < cuts[0]) std::swap(cuts[0], cuts[1]);
    *cut_begin = cuts[0];
    *cut_end = cuts[1];
  }

  double my_center[3];
  double easy_normals[9];
};


TEST_F(SphereSliceTest, domain) {
  EXPECT_THROW(cl::SphereSlice(my_center, easy_normals, 0.0), std::domain_error);
  EXPECT_THROW(cl::SphereSlice(my_center, easy_normals, -1.0), std::domain_error);
  const double degenerate_normals[9] = {1.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0};
  EXPECT_THROW(cl::SphereSlice(my_center, degenerate_normals, 5.0), std::domain_error);
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  EXPECT_THROW(slice.solve_range(-1, &begin, &end), std::domain_error);
  EXPECT_THROW(slice.solve_range(3, &begin, &end), std::domain_error);
  EXPECT_THROW(slice.set_cut_begin_end(-1, 1, 2), std::domain_error);
  EXPECT_THROW(slice.set_cut_begin_end(2, 1, 2), std::domain_error);
  EXPECT_THROW(slice.set_cut_begin_end(0, 1, -2), std::domain_error);
  EXPECT_THROW(slice.set_cut_begin_end(1, 2, 1), std::domain_error);
  EXPECT_THROW(slice.solve_full_low(-1, &begin, &end, nullptr, nullptr), std::domain_error);
  EXPECT_THROW(slice.solve_full_low(3, &begin, &end, nullptr, nullptr), std::domain_error);
  EXPECT_THROW(slice.solve_plane_low(-1, 0, 0.0, &begin, &end, nullptr, nullptr),
               std::domain_error);
  EXPECT_THROW(slice.solve_plane_low(3, 0, 0.0, &begin, &end, nullptr, nullptr),
               std::domain_error);
  EXPECT_THROW(slice.solve_plane_low(0, -1, 0.0, &begin, &end, nullptr, nullptr),
               std::domain_error);
  EXPECT_THROW(slice.solve_plane_low(0, 3, 0.0, &begin, &end, nullptr, nullptr),
               std::domain_error);
  EXPECT_THROW(slice.solve_line_low(-1, 0, 1, 0.0, 0.0, &begin, &end, nullptr, nullptr),
               std::domain_error);
  EXPECT_THROW(slice.solve_line_low(3, 0, 1, 0.0, 0.0, &begin, &end, nullptr, nullptr),
              std::domain_error);
  EXPECT_THROW(slice.solve_line_low(0, -1, 1, 0.0, 0.0, &begin, &end, nullptr, nullptr),
               std::domain_error);
  EXPECT_THROW(slice.solve_line_low(0, 3, 1, 0.0, 0.0, &begin, &end, nullptr, nullptr),
               std::domain_error);
  EXPECT_THROW(slice.solve_line_low(0, 1, -1, 0.0, 0.0, &begin, &end, nullptr, nullptr),
               std::domain_error);
  EXPECT_THROW(slice.solve_line_low(0, 1, 3, 0.0, 0.0, &begin, &end, nullptr, nullptr),
               std::domain_error);
  EXPECT_THROW(slice.compute_plane_intersection(-1, 0, 0.0, 0.0, nullptr), std::domain_error);
  EXPECT_THROW(slice.compute_plane_intersection(3, 0, 0.0, 0.0, nullptr), std::domain_error);
  EXPECT_THROW(slice.compute_plane_intersection(0, -1, 0.0, 0.0, nullptr), std::domain_error);
  EXPECT_THROW(slice.compute_plane_intersection(0, 3, 0.0, 0.0, nullptr), std::domain_error);
}


TEST_F(SphereSliceTest, solve_full_low_example1) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.solve_range(0, &begin, &end);
  EXPECT_DOUBLE_EQ(-4.6, begin);
  EXPECT_DOUBLE_EQ(5.4, end);
}


TEST_F(SphereSliceTest, solve_plane_low_example_ortho) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.solve_plane_low(1, 0, 3.4, &begin, &end, nullptr, nullptr);
  EXPECT_DOUBLE_EQ(-6.0, begin);
  EXPECT_DOUBLE_EQ(2.0, end);
}


TEST_F(SphereSliceTest, solve_plane_low_example_angle) {
  easy_normals[3] = 1.0;
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.solve_plane_low(1, 0, 3.4, &begin, &end, nullptr, nullptr);
  EXPECT_DOUBLE_EQ(-2.6, begin);
  EXPECT_DOUBLE_EQ(5.4, end);
}


TEST_F(SphereSliceTest, solve_plane_low_example_nonexisting) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.solve_plane_low(1, 0, 50.0, &begin, &end, nullptr, nullptr);
  EXPECT_TRUE(std::isnan(begin));
  EXPECT_TRUE(std::isnan(end));
}


TEST_F(SphereSliceTest, solve_range_1_example1) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.set_cut_begin_end(0, 0.0, 1.0);
  slice.solve_range(1, &begin, &end);
  EXPECT_DOUBLE_EQ(-7.0, begin);
  EXPECT_DOUBLE_EQ(3.0, end);
}


TEST_F(SphereSliceTest, solve_range_1_example2) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.set_cut_begin_end(0, 3.4, 4.4);
  slice.solve_range(1, &begin, &end);
  EXPECT_DOUBLE_EQ(-6.0, begin);
  EXPECT_DOUBLE_EQ(2.0, end);
}


TEST_F(SphereSliceTest, solve_range_1_example3) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.set_cut_begin_end(0, -3.6, -2.6);
  slice.solve_range(1, &begin, &end);
  EXPECT_DOUBLE_EQ(-6.0, begin);
  EXPECT_DOUBLE_EQ(2.0, end);
}


TEST_F(SphereSliceTest, solve_range_1_example_nofound) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.set_cut_begin_end(0, -10, -9);
  EXPECT_THROW(slice.solve_range(1, &begin, &end), cl::no_solution_found);
}


TEST_F(SphereSliceTest, solve_line_low_example_ortho) {
  cl::SphereSlice slice(my_center, easy_normals, sqrt(14));
  double begin, end;
  slice.solve_line_low(2, 0, 1, 3.4, -1.0, &begin, &end, nullptr, nullptr);
  EXPECT_DOUBLE_EQ(-1.0, begin);
  EXPECT_DOUBLE_EQ(3.0, end);
}


TEST_F(SphereSliceTest, solve_line_low_example_angle) {
  easy_normals[1] = 1.0;
  cl::SphereSlice slice(my_center, easy_normals, sqrt(14));
  double begin, end;
  slice.solve_line_low(2, 0, 1, 1.4, 0.0, &begin, &end, nullptr, nullptr);
  EXPECT_DOUBLE_EQ(-2.0, begin);
  EXPECT_DOUBLE_EQ(4.0, end);
}


TEST_F(SphereSliceTest, solve_line_low_example_nonexisting) {
  cl::SphereSlice slice(my_center, easy_normals, sqrt(14));
  double begin, end;
  slice.solve_line_low(2, 0, 1, 50.0, 30.0, &begin, &end, nullptr, nullptr);
  EXPECT_TRUE(std::isnan(begin));
  EXPECT_TRUE(std::isnan(end));
}


TEST_F(SphereSliceTest, solve_range_2_example1) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.set_cut_begin_end(0, 0.0, 1.0);
  slice.set_cut_begin_end(1, -2.2, 1.0);
  slice.solve_range(2, &begin, &end);
  EXPECT_DOUBLE_EQ(-4.0, begin);
  EXPECT_DOUBLE_EQ(6.0, end);
}


TEST_F(SphereSliceTest, solve_range_2_example2) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.set_cut_begin_end(0, 3.4, 4.0);
  slice.set_cut_begin_end(1, -2.2, 1.0);
  slice.solve_range(2, &begin, &end);
  EXPECT_DOUBLE_EQ(-3.0, begin);
  EXPECT_DOUBLE_EQ(5.0, end);
}


TEST_F(SphereSliceTest, solve_range_2_example3) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.set_cut_begin_end(0, -5.0, -2.6);
  slice.set_cut_begin_end(1, -2.2, 1.0);
  slice.solve_range(2, &begin, &end);
  EXPECT_DOUBLE_EQ(-3.0, begin);
  EXPECT_DOUBLE_EQ(5.0, end);
}


TEST_F(SphereSliceTest, solve_range_2_example4) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  slice.set_cut_begin_end(0, 0.0, 1.0);
  slice.set_cut_begin_end(1, -6.2, -5.0);
  slice.solve_range(2, &begin, &end);
  EXPECT_DOUBLE_EQ(-3.0, begin);
  EXPECT_DOUBLE_EQ(5.0, end);
}


TEST_F(SphereSliceTest, solve_range_2_example5) {
  cl::SphereSlice slice(my_center, easy_normals, 5.0);
  double begin, end;
  // 0.4, -2.0, 1.0
  slice.set_cut_begin_end(0, 0.0, 1.0);
  slice.set_cut_begin_end(1, 1.0, 2.0);
  slice.solve_range(2, &begin, &end);
  EXPECT_DOUBLE_EQ(-3.0, begin);
  EXPECT_DOUBLE_EQ(5.0, end);
}


TEST_F(SphereSliceTest, solve_range_2_example6) {
  cl::SphereSlice slice(my_center, easy_normals, sqrt(14.0));
  double begin, end;
  // 0.4, -2.0, 1.0
  slice.set_cut_begin_end(0, 1.4, 2.0);
  slice.set_cut_begin_end(1, 0.0, 1.0);
  slice.solve_range(2, &begin, &end);
  EXPECT_DOUBLE_EQ(-2.0, begin);
  EXPECT_DOUBLE_EQ(4.0, end);
}


TEST_F(SphereSliceTest, solve_range_2_example7) {
  cl::SphereSlice slice(my_center, easy_normals, sqrt(14.0));
  double begin, end;
  slice.set_cut_begin_end(0, 1.4, 2.0);
  slice.set_cut_begin_end(1, -5.2, -4.0);
  slice.solve_range(2, &begin, &end);
  EXPECT_DOUBLE_EQ(-2.0, begin);
  EXPECT_DOUBLE_EQ(4.0, end);
}


TEST_F(SphereSliceTest, solve_range_2_example8) {
  cl::SphereSlice slice(my_center, easy_normals, sqrt(14.0));
  double begin, end;
  slice.set_cut_begin_end(0, -2.0, -0.6);
  slice.set_cut_begin_end(1, 0.0, 1.0);
  slice.solve_range(2, &begin, &end);
  EXPECT_DOUBLE_EQ(-2.0, begin);
  EXPECT_DOUBLE_EQ(4.0, end);
}


TEST_F(SphereSliceTest, solve_range_2_example9) {
  cl::SphereSlice slice(my_center, easy_normals, sqrt(14.0));
  double begin, end;
  slice.set_cut_begin_end(0, -2.0, -0.6);
  slice.set_cut_begin_end(1, -5.2, -4.0);
  slice.solve_range(2, &begin, &end);
  EXPECT_DOUBLE_EQ(-2.0, begin);
  EXPECT_DOUBLE_EQ(4.0, end);
}


TEST_F(SphereSliceTest, solve_range_2_example_nofound) {
  cl::SphereSlice slice(my_center, easy_normals, 1.0);
  double begin, end;
  slice.set_cut_begin_end(0, 20.0, 30.0);
  slice.set_cut_begin_end(1, -5.2, -4.0);
  EXPECT_THROW(slice.solve_range(2, &begin, &end), cl::no_solution_found);
}


TEST_F(SphereSliceTest, compute_plane_intersection_example1) {
  cl::SphereSlice slice(my_center, easy_normals, sqrt(14.0));

  double line_center[3];
  double dist_sq = slice.compute_plane_intersection(0, 1, 3.0, 4.0, line_center);
  EXPECT_DOUBLE_EQ(dist_sq, 25.0);
  EXPECT_DOUBLE_EQ(line_center[0], 3.0);
  EXPECT_DOUBLE_EQ(line_center[1], 4.0);
  EXPECT_DOUBLE_EQ(line_center[2], 0.0);
}


TEST_F(SphereSliceTest, solve_full_low_random) {
  int num_inside = 0;
  for (int irep=0; irep < NREP; ++irep) {
    // Test parameters
    double radius = (irep + 1)*0.1;
    double center[3];
    double normals[9];
    std::unique_ptr<cl::SphereSlice>
        slice(create_random_problem(irep, radius, center, normals));

    // Do a solve_full_low
    double begin, end;
    double point_begin[3];
    double point_end[3];
    slice->solve_full_low(0, &begin, &end, point_begin, point_end);

    // order of begin and end must be right
    EXPECT_LE(begin, end);

    // Check consistency begin, point_begin
    EXPECT_NEAR(begin, vec3::dot(point_begin, normals), EPS);

    // Check consistency begin, point_end
    EXPECT_NEAR(end, vec3::dot(point_end, normals), EPS);

    // Check consistency when not using point_begin, point_end.
    double begin_bis, end_bis;
    slice->solve_full_low(0, &begin_bis, &end_bis, nullptr, nullptr);
    EXPECT_DOUBLE_EQ(begin, begin_bis);
    EXPECT_DOUBLE_EQ(end, end_bis);

    for (int ipoint=0; ipoint < NPOINT; ++ipoint) {
      // Random point
      double point[3];
      double norm;
      random_point(ipoint, center, radius, point, &norm);

      // If the point is in the sphere, test if reduced coordinate falls
      // in the range [begin,end].
      double proj = vec3::dot(point, normals);
      if (norm < radius) {
        EXPECT_LE(begin, proj);
        EXPECT_GE(end, proj);
        ++num_inside;
      }
    }
  }
  // Check sufficiency
  EXPECT_LT((NREP*NPOINT)/3, num_inside);
}


TEST_F(SphereSliceTest, solve_range_0_random) {
  for (int irep=0; irep < NREP; ++irep) {
    // Test parameters
    double radius = (irep + 1)*0.1;
    double center[3];
    double normals[9];
    std::unique_ptr<cl::SphereSlice> slice(create_random_problem(irep, radius, center, normals));

    // Do a solve_full_low
    double begin, end;
    slice->solve_full_low(0, &begin, &end, nullptr, nullptr);

    // order of begin and end must be right
    EXPECT_LE(begin, end);

    // Do a solve_range_zero
    double begin_bis, end_bis;
    slice->solve_range_0(&begin_bis, &end_bis);

    // Should be the same
    EXPECT_EQ(begin, begin_bis);
    EXPECT_EQ(end, end_bis);
  }
}


TEST_F(SphereSliceTest, solve_plane_low_random) {
  for (int irep=0; irep < NREP; ++irep) {
    // Test parameters
    double radius = (irep + 1)*0.1;
    double center[3];
    double normals[9];
    std::unique_ptr<cl::SphereSlice> slice(create_random_problem(irep, radius, center, normals));

    // Select the random ids for cut_normal and axis
    int permutation[3];
    fill_random_permutation(irep*5 + 123, permutation, 3);
    int id_cut = permutation[0];
    int id_axis = permutation[1];

    // Name some normals for convenience
    double* cut_normal = normals + 3*id_cut;
    double* axis = normals + 3*id_axis;

    // Select randomized place to cut the sphere
    double cut, cut_min, cut_max;
    random_cut(irep + 12345, *slice, id_cut, &cut, &cut_min, &cut_max);

    // Actual computation
    double begin, end;
    double point_begin[3];
    double point_end[3];
    slice->solve_plane_low(id_axis, id_cut, cut, &begin, &end, point_begin, point_end);

    // It should have worked...
    EXPECT_FALSE(std::isnan(begin));
    EXPECT_FALSE(std::isnan(end));

    // order of begin and end must be right
    EXPECT_LE(begin, end);

    // Points must be on sphere...
    EXPECT_NEAR(radius, vec3::distance(center, point_begin), EPS);
    EXPECT_NEAR(radius, vec3::distance(center, point_end), EPS);

    // Points should be consistent with begin and end
    EXPECT_NEAR(begin, vec3::dot(point_begin, axis), EPS);
    EXPECT_NEAR(end, vec3::dot(point_end, axis), EPS);

    // Points should be on the cut
    EXPECT_NEAR(cut, vec3::dot(point_begin, cut_normal), EPS);
    EXPECT_NEAR(cut, vec3::dot(point_end, cut_normal), EPS);

    // Check that cut_normal, axis and (point_end - point_begin) lie in
    // the same plane.
    double delta[3];
    vec3::delta(point_end, point_begin, delta);
    EXPECT_NEAR(0.0, vec3::triple(cut_normal, axis, delta), EPS);

    // Check that cut_normal is orthogonal to delta.
    EXPECT_NEAR(0.0, vec3::dot(cut_normal, delta), EPS);

    // Parameterize the circle
    double circle_center[3];
    vec3::copy(point_begin, circle_center);
    vec3::iadd(circle_center, point_end);
    vec3::iscale(circle_center, 0.5);
    double circle_basis0[3];
    vec3::copy(delta, circle_basis0);
    vec3::iscale(circle_basis0, 0.5);
    double circle_radius = vec3::norm(circle_basis0);
    double circle_basis1[3];
    vec3::cross(cut_normal, circle_basis0, circle_basis1);
    vec3::iscale(circle_basis1, circle_radius/vec3::norm(circle_basis1));

    // Check the circle center
    EXPECT_NEAR(cut, vec3::dot(cut_normal, circle_center), EPS);

    // Scan the circle close to the solutions, where optimality is easily
    // tested.
    for (double angle = -0.01; angle < 0.01; angle += 0.001) {
      double other[3];
      double proj_axis;
      // Check end
      vec3::copy(circle_center, other);
      vec3::iadd(other, circle_basis0, cos(angle));
      vec3::iadd(other, circle_basis1, sin(angle));
      proj_axis = vec3::dot(other, axis);
      EXPECT_GE(end + EPS, proj_axis);
      // Check begin
      vec3::copy(circle_center, other);
      vec3::iadd(other, circle_basis0, -cos(angle));
      vec3::iadd(other, circle_basis1, -sin(angle));
      proj_axis = vec3::dot(other, axis);
      EXPECT_LE(begin-EPS, proj_axis);
    }

    // Call without point_* arguments
    double begin_bis, end_bis;
    slice->solve_plane_low(id_axis, id_cut, cut, &begin_bis, &end_bis, nullptr, nullptr);

    // Result should be the same
    EXPECT_NEAR(begin, begin_bis, EPS);
    EXPECT_NEAR(end, end_bis, EPS);
  }
}


TEST_F(SphereSliceTest, solve_range_1_random) {
  int num_inside = 0;
  for (int irep=0; irep < NREP; ++irep) {
    // Test parameters
    double radius = (irep + 1)*0.1;
    double center[3];
    double normals[9];
    std::unique_ptr<cl::SphereSlice> slice(create_random_problem(irep, radius, center, normals));

    // Name some normals for convenience
    double* cut_normal = normals;
    double* axis = normals + 3;

    double cut_begin, cut_end;
    double cut_min, cut_max;
    random_slice(irep*2 + 1, *slice, 0, &cut_begin, &cut_end, &cut_min, &cut_max);

    // Do the computation, this should always work!
    double axis_begin, axis_end;
    slice->set_cut_begin_end(0, cut_begin, cut_end);
    slice->solve_range_1(&axis_begin, &axis_end);

    // order of begin and end must be right
    EXPECT_LE(axis_begin, axis_end);

    // Check if the solution of solve_range_1 is better or as good as the
    // solutions of the separate parts
    double axis_begin0, axis_end0;
    slice->solve_plane_low(1, 0, cut_begin, &axis_begin0, &axis_end0, nullptr, nullptr);
    EXPECT_FALSE(std::isnan(axis_begin0));
    EXPECT_FALSE(std::isnan(axis_end0));
    EXPECT_LE(axis_begin0, axis_end0);
    EXPECT_LE(axis_begin, axis_begin0);
    EXPECT_GE(axis_end, axis_end0);
    double axis_begin1, axis_end1;
    slice->solve_plane_low(1, 0, cut_end, &axis_begin1, &axis_end1, nullptr, nullptr);
    EXPECT_FALSE(std::isnan(axis_begin1));
    EXPECT_FALSE(std::isnan(axis_end1));
    EXPECT_LE(axis_begin1, axis_end1);
    EXPECT_LE(axis_begin, axis_begin1);
    EXPECT_GE(axis_end, axis_end1);
    double point_begin[3];
    double point_end[3];
    // If the sphere solution is in the proper range, it is the solution
    double axis_begin_sphere, axis_end_sphere;
    slice->solve_full_low(1, &axis_begin_sphere, &axis_end_sphere, point_begin, point_end);
    EXPECT_LE(axis_begin_sphere, axis_end_sphere);
    double proj_begin = vec3::dot(cut_normal, point_begin);
    if ((proj_begin > cut_begin) && (proj_begin < cut_end)) {
      EXPECT_EQ(axis_begin, axis_begin_sphere);
    }
    double proj_end = vec3::dot(cut_normal, point_end);
    if ((proj_end > cut_begin) && (proj_end < cut_end)) {
      EXPECT_EQ(axis_end, axis_end_sphere);
    }

    // Test by generating lots of random points. If a point is in the disc,
    // its fractional coordinate should be between begin and end.
    for (int ipoint=0; ipoint < NPOINT; ++ipoint) {
      // Random point
      double point[3];
      double norm;
      random_point(ipoint, center, radius, point, &norm);
      ASSERT_NEAR(norm, vec3::distance(point, center), EPS);
      if (norm < radius) {
        // Projection on axis should always be in the "sphere range"
        double proj2 = vec3::dot(point, axis);
        EXPECT_LE(axis_begin_sphere, proj2);
        EXPECT_GE(axis_end_sphere, proj2);
        // If in the slice, test if proj2 is in bounds given by solve_range_1
        double proj = vec3::dot(point, cut_normal);
        EXPECT_LE(cut_min, proj);
        EXPECT_GE(cut_max, proj);
        if ((proj > cut_begin) && (proj < cut_end)) {
          EXPECT_LE(axis_begin, proj2);
          EXPECT_GE(axis_end, proj2);
          ++num_inside;
        }
      }
    }
  }
  // Sufficiency check
  EXPECT_LT((NREP*NPOINT)/10, num_inside);
}


TEST_F(SphereSliceTest, solve_line_low_random) {
  int num_inside = 0;
  for (int irep=0; irep < NREP; ++irep) {
    // Test parameters
    double radius = (irep + 1)*0.1;
    double center[3];
    double normals[9];
    std::unique_ptr<cl::SphereSlice> slice(create_random_problem(irep, radius, center, normals));

    // Select the random ids for cut_normal and axis
    int permutation[3];
    fill_random_permutation(irep*5 + 123, permutation, 3);
    int id_cut0 = permutation[0];
    int id_cut1 = permutation[1];
    int id_axis = permutation[2];

    // Name some normals for convenience
    double* cut0_normal = normals + 3*id_cut0;
    double* cut1_normal = normals + 3*id_cut1;
    double* axis = normals + 3*id_axis;

    // Select randomized places to cut the sphere
    double cut0, cut0_min, cut0_max;
    double cut1, cut1_min, cut1_max;
    random_cut(irep + 12345, *slice, id_cut0, &cut0, &cut0_min, &cut0_max);
    random_cut(irep*2 + 114, *slice, id_cut1, &cut1, &cut1_min, &cut1_max);

    // Compute the distance from the origin to the intersection
    double delta_cut0 = cut0 - vec3::dot(center, cut0_normal);
    double delta_cut1 = cut1 - vec3::dot(center, cut1_normal);
    double dist_line_center = sqrt(slice->compute_plane_intersection(id_cut0, id_cut1,
      delta_cut0, delta_cut1, nullptr));

    // Actual computation, may not work in some cases.
    double begin, end;
    double point_begin[3];
    double point_end[3];
    slice->solve_line_low(id_axis, id_cut0, id_cut1, cut0, cut1, &begin, &end,
      point_begin, point_end);

    if (dist_line_center > radius) {
      // It should not have worked ...
      EXPECT_TRUE(std::isnan(begin));
      EXPECT_TRUE(std::isnan(end));
      continue;
    } else {
      // It should have worked...
      EXPECT_FALSE(std::isnan(begin));
      EXPECT_FALSE(std::isnan(end));
      ++num_inside;
    }

    // order of begin and end must be right
    EXPECT_LE(begin, end);

    // Points must be on sphere...
    EXPECT_NEAR(radius, vec3::distance(center, point_begin), EPS);
    EXPECT_NEAR(radius, vec3::distance(center, point_end), EPS);

    // Points should be consistent with begin and end
    EXPECT_NEAR(begin, vec3::dot(point_begin, axis), EPS);
    EXPECT_NEAR(end, vec3::dot(point_end, axis), EPS);

    // Points should be on the cut
    EXPECT_NEAR(cut0, vec3::dot(point_begin, cut0_normal), EPS);
    EXPECT_NEAR(cut0, vec3::dot(point_end, cut0_normal), EPS);
    EXPECT_NEAR(cut1, vec3::dot(point_begin, cut1_normal), EPS);
    EXPECT_NEAR(cut1, vec3::dot(point_end, cut1_normal), EPS);

    // Check that cut0_normal and cut1_normal are orthogonal to delta.
    double delta[3];
    vec3::delta(point_end, point_begin, delta);
    EXPECT_NEAR(0.0, vec3::dot(cut0_normal, delta), EPS);
    EXPECT_NEAR(0.0, vec3::dot(cut1_normal, delta), EPS);

    // Get the line center and radius
    double line_center[3];
    vec3::copy(point_begin, line_center);
    vec3::iadd(line_center, point_end);
    vec3::iscale(line_center, 0.5);

    // Check the line center
    EXPECT_NEAR(cut0, vec3::dot(cut0_normal, line_center), EPS);
    EXPECT_NEAR(cut1, vec3::dot(cut1_normal, line_center), EPS);

    // Call without point_* arguments
    double begin_bis, end_bis;
    slice->solve_line_low(id_axis, id_cut0, id_cut1, cut0, cut1, &begin_bis,
      &end_bis, nullptr, nullptr);

    // Result should be the same
    EXPECT_NEAR(begin, begin_bis, EPS);
    EXPECT_NEAR(end, end_bis, EPS);
  }
  // Sufficiency check
  EXPECT_LT(NREP/3, num_inside);
}


TEST_F(SphereSliceTest, solve_range_2_random) {
  int num_inside = 0;
  for (int irep=0; irep < NREP; ++irep) {
    // Test parameters
    double radius = (irep + 1)*0.1;
    double center[3];
    double normals[9];
    std::unique_ptr<cl::SphereSlice> slice(
      create_random_problem(irep, radius, center, normals));

    // Name some normals for convenience
    double* cut0_normal = normals;
    double* cut1_normal = normals + 3;
    double* axis = normals + 6;

    // Take two random slices to form a random bar
    double cut0_begin, cut0_end, cut0_min, cut0_max;
    random_slice(irep*3 + 1, *slice, 0, &cut0_begin, &cut0_end, &cut0_min, &cut0_max);
    double cut1_begin, cut1_end, cut1_min, cut1_max;
    random_slice(irep*3 + 2, *slice, 1, &cut1_begin, &cut1_end, &cut1_min, &cut1_max);

    // Do the computation
    double axis_begin, axis_end;
    slice->set_cut_begin_end(0, cut0_begin, cut0_end);
    slice->set_cut_begin_end(1, cut1_begin, cut1_end);
    try {
      slice->solve_range_2(&axis_begin, &axis_end);
    } catch (cl::no_solution_found) {
      // The rest of this iteration of the test makes no sense.
      continue;
    }

    // order of begin and end must be right
    EXPECT_LE(axis_begin, axis_end);

    // Check if the solution of solve_range_2 is better or as good as the
    // solutions of the 9 separate parts (A-I).
    bool exists;

    // * case A: cut0_begin  cut1_begin
    double axis_begin_a, axis_end_a;
    double point_begin_a[3], point_end_a[3];
    slice->solve_line_low(2, 0, 1, cut0_begin, cut1_begin, &axis_begin_a, &axis_end_a,
                          point_begin_a, point_end_a);
    exists = std::isfinite(axis_begin_a) && std::isfinite(axis_end_a);
    if (exists) {
      EXPECT_LE(axis_begin_a, axis_end_a + EPS);
      EXPECT_LE(axis_begin, axis_begin_a + EPS);
      EXPECT_GE(axis_end, axis_end_a - EPS);
    }

    // * case B: cut0_begin  cut1_end
    double axis_begin_b, axis_end_b;
    double point_begin_b[3], point_end_b[3];
    slice->solve_line_low(2, 0, 1, cut0_begin, cut1_end, &axis_begin_b, &axis_end_b,
                          point_begin_b, point_end_b);
    exists = std::isfinite(axis_begin_b) && std::isfinite(axis_end_b);
    if (exists) {
      EXPECT_LE(axis_begin_b, axis_end_b + EPS);
      EXPECT_LE(axis_begin, axis_begin_b + EPS);
      EXPECT_GE(axis_end, axis_end_b - EPS);
    }

    // * case C: cut0_end    cut1_begin
    double axis_begin_c, axis_end_c;
    double point_begin_c[3], point_end_c[3];
    slice->solve_line_low(2, 0, 1, cut0_end, cut1_begin, &axis_begin_c, &axis_end_c,
                          point_begin_c, point_end_c);
    exists = std::isfinite(axis_begin_c) && std::isfinite(axis_end_c);
    if (exists) {
      EXPECT_LE(axis_begin_c, axis_end_c + EPS);
      EXPECT_LE(axis_begin, axis_begin_c + EPS);
      EXPECT_GE(axis_end, axis_end_c - EPS);
    }

    // * case D: cut0_end    cut1_end
    double axis_begin_d, axis_end_d;
    double point_begin_d[3], point_end_d[3];
    slice->solve_line_low(2, 0, 1, cut0_end, cut1_end, &axis_begin_d, &axis_end_d,
                          point_begin_d, point_end_d);
    exists = std::isfinite(axis_begin_d) && std::isfinite(axis_end_d);
    if (exists) {
      EXPECT_LE(axis_begin_d, axis_end_d + EPS);
      EXPECT_LE(axis_begin, axis_begin_d + EPS);
      EXPECT_GE(axis_end, axis_end_d - EPS);
    }

    // * case E: cut0_begin
    double axis_begin_e, axis_end_e;
    double frac1_begin_e, frac1_end_e;
    double point_begin_e[3], point_end_e[3];
    slice->solve_plane_low(2, 0, cut0_begin, &axis_begin_e, &axis_end_e,
                           point_begin_e, point_end_e);
    exists = std::isfinite(axis_begin_e) && std::isfinite(axis_end_e);
    if (exists) {
      EXPECT_LE(axis_begin_e, axis_end_e);
      frac1_begin_e = vec3::dot(cut1_normal, point_begin_e);
      if ((frac1_begin_e > cut1_begin) && (frac1_begin_e < cut1_end)) {
        EXPECT_LE(axis_begin, axis_begin_e);
      }
      frac1_end_e = vec3::dot(cut1_normal, point_end_e);
      if ((frac1_end_e > cut1_begin) && (frac1_end_e < cut1_end)) {
        EXPECT_GE(axis_end, axis_end_e);
      }
    }

    // * case F: cut0_end
    double axis_begin_f, axis_end_f;
    double frac1_begin_f, frac1_end_f;
    double point_begin_f[3], point_end_f[3];
    slice->solve_plane_low(2, 0, cut0_end, &axis_begin_f, &axis_end_f, point_begin_f, point_end_f);
    exists = std::isfinite(axis_begin_f) && std::isfinite(axis_end_f);
    if (exists) {
      EXPECT_LE(axis_begin_f, axis_end_f);
      frac1_begin_f = vec3::dot(cut1_normal, point_begin_f);
      if ((frac1_begin_f > cut1_begin) && (frac1_begin_f < cut1_end)) {
        EXPECT_LE(axis_begin, axis_begin_f);
      }
      frac1_end_f = vec3::dot(cut1_normal, point_end_f);
      if ((frac1_end_f > cut1_begin) && (frac1_end_f < cut1_end)) {
        EXPECT_GE(axis_end, axis_end_f);
      }
    }

    // * case G: cut1_begin
    double axis_begin_g, axis_end_g;
    double frac0_begin_g, frac0_end_g;
    double point_begin_g[3], point_end_g[3];
    slice->solve_plane_low(2, 1, cut1_begin, &axis_begin_g, &axis_end_g,
                           point_begin_g, point_end_g);
    exists = std::isfinite(axis_begin_g) && std::isfinite(axis_end_g);
    if (exists) {
      EXPECT_LE(axis_begin_g, axis_end_g);
      frac0_begin_g = vec3::dot(cut0_normal, point_begin_g);
      if ((frac0_begin_g > cut0_begin) && (frac0_begin_g < cut0_end)) {
        EXPECT_LE(axis_begin, axis_begin_g);
      }
      frac0_end_g = vec3::dot(cut0_normal, point_end_g);
      if ((frac0_end_g > cut0_begin) && (frac0_end_g < cut0_end)) {
        EXPECT_GE(axis_end, axis_end_g);
      }
    }

    // * case H: cut1_end
    double axis_begin_h, axis_end_h;
    double frac0_begin_h, frac0_end_h;
    double point_begin_h[3], point_end_h[3];
    slice->solve_plane_low(2, 1, cut1_end, &axis_begin_h, &axis_end_h, point_begin_h, point_end_h);
    exists = std::isfinite(axis_begin_h) && std::isfinite(axis_end_h);
    if (exists) {
      EXPECT_LE(axis_begin_h, axis_end_h);
      frac0_begin_h = vec3::dot(cut0_normal, point_begin_h);
      if ((frac0_begin_h > cut0_begin) && (frac0_begin_h < cut0_end)) {
        EXPECT_LE(axis_begin, axis_begin_h);
      }
      frac0_end_h = vec3::dot(cut0_normal, point_end_h);
      if ((frac0_end_h > cut0_begin) && (frac0_end_h < cut0_end)) {
        EXPECT_GE(axis_end, axis_end_h);
      }
    }

    // * case I: (none)
    //   If the sphere solution is in the proper range, it is the solution.
    double axis_begin_i, axis_end_i;
    double frac0_begin_i, frac0_end_i;
    double frac1_begin_i, frac1_end_i;
    double point_begin_i[3], point_end_i[3];

    slice->solve_full_low(2, &axis_begin_i, &axis_end_i, point_begin_i, point_end_i);
    EXPECT_FALSE(std::isnan(axis_begin_i));
    EXPECT_FALSE(std::isnan(axis_end_i));
    EXPECT_LE(axis_begin_i, axis_end_i);
    frac0_begin_i = vec3::dot(cut0_normal, point_begin_i);
    frac1_begin_i = vec3::dot(cut1_normal, point_begin_i);
    if ((frac0_begin_i > cut0_begin) && (frac0_begin_i < cut0_end) &&
      (frac1_begin_i > cut1_begin) && (frac1_begin_i < cut1_end)) {
      EXPECT_EQ(axis_begin, axis_begin_i);
    }
    frac0_end_i = vec3::dot(cut0_normal, point_end_i);
    frac1_end_i = vec3::dot(cut1_normal, point_end_i);
    if ((frac0_end_i > cut0_begin) && (frac0_end_i < cut0_end) &&
      (frac1_end_i > cut1_begin) && (frac1_end_i < cut1_end)) {
      EXPECT_EQ(axis_end, axis_end_i);
    }

    // Test by generating lots of random points. If a point is in the bar,
    // its fractional coordinate should be between begin and end.
    for (int ipoint=0; ipoint < NPOINT; ++ipoint) {
      // Random point
      double point[3];
      double norm;
      random_point(ipoint, center, radius, point, &norm);
      ASSERT_NEAR(norm, vec3::distance(point, center), EPS);
      if (norm < radius) {
        // Projection on axis should always be in the "sphere range"
        double frac_axis = vec3::dot(point, axis);
        EXPECT_LE(axis_begin_i, frac_axis);
        EXPECT_GE(axis_end_i, frac_axis);
        // If in the slice, test if proj2 is in bounds given by solve_range_1
        double frac_cut0 = vec3::dot(point, cut0_normal);
        double frac_cut1 = vec3::dot(point, cut1_normal);
        EXPECT_LE(cut0_min, frac_cut0);
        EXPECT_GE(cut0_max, frac_cut0);
        EXPECT_LE(cut1_min, frac_cut1);
        EXPECT_GE(cut1_max, frac_cut1);
        if ((frac_cut0 > cut0_begin) && (frac_cut0 < cut0_end) &&
          (frac_cut1 > cut1_begin) && (frac_cut1 < cut1_end)) {
          EXPECT_LE(axis_begin, frac_axis);
          EXPECT_GE(axis_end, frac_axis);
          ++num_inside;
        }
      }
    }
  }
  // Sufficiency check
  EXPECT_LT((NREP*NPOINT)/20, num_inside);
}


TEST(ComputeBeginEndTest, example1) {
  // Test parameters
  const double center[3] = {1.0, 0.5, -2.0};
  const double ortho[3] = {-0.5, 2.2, 1.5};
  const double axis[3] = {1.0, 2.0, 3.0};

  // Check the output of a call with all arguments.
  double begin, end;
  double point_begin[3];
  double point_end[3];
  cl::compute_begin_end(center, ortho, axis, &begin, &end, point_begin, point_end);
  EXPECT_DOUBLE_EQ(1.5, point_begin[0]);
  EXPECT_DOUBLE_EQ(-1.7, point_begin[1]);
  EXPECT_DOUBLE_EQ(-3.5, point_begin[2]);
  EXPECT_DOUBLE_EQ(0.5, point_end[0]);
  EXPECT_DOUBLE_EQ(2.7, point_end[1]);
  EXPECT_DOUBLE_EQ(-0.5, point_end[2]);
  EXPECT_DOUBLE_EQ(-12.4, begin);
  EXPECT_DOUBLE_EQ(4.4, end);

  // Also compare to the output of a call without point_* arguments.
  double begin_bis, end_bis;
  cl::compute_begin_end(center, ortho, axis, &begin_bis, &end_bis, nullptr, nullptr);

  EXPECT_DOUBLE_EQ(begin, begin_bis);
  EXPECT_DOUBLE_EQ(end, end_bis);
}

// vim: textwidth=90 et ts=2 sw=2
