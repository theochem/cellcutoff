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
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--


#include <stdexcept>
#include <cstdlib>
#include <gtest/gtest.h>
#include <cmath>
#include "celllists/sphere_slice.h"
#include "celllists/vec3.h"
#include "common.h"


/*
    TODO
    - unsigned int seed always first in random functions
    - min and max arguments for fill_random_double, fill_random_long
    - randomize arguments of solve_circle(1, 0, ...) and solve_line(2, 0, 1,
      ...) in tests.
 */

class SphereSliceTest : public ::testing::Test {
    public:
        double my_center[3];
        double easy_normals[9];

        void SetUp() {
            my_center[0] = 0.4;
            my_center[1] = -2.0;
            my_center[2] = 1.0;
            std::fill(easy_normals, easy_normals+9, 0.0);
            easy_normals[0] = 1.0;
            easy_normals[4] = 1.0;
            easy_normals[8] = 1.0;
        }

        SphereSlice* create_random_problem(unsigned int seed, double rcut,
            double* center, double* normals)
        {
            fill_random_double(center, 3, seed);
            fill_random_double(normals, 9, 1342+seed);
            return new SphereSlice(center, normals, rcut);
        }

        void random_point(unsigned int seed, double rcut, const double* center,
            double* point, double &norm)
        {
            fill_random_double(point, 3, seed, rcut);
            norm = vec3::norm(point);
            point[0] += center[0];
            point[1] += center[1];
            point[2] += center[2];
        }

        void random_cut(unsigned int seed, SphereSlice* slice, int id_cut,
            double &cut, double &cut_min, double &cut_max) {
            // Do a solve_sphere, to know over which range we can cut
            slice->solve_sphere(id_cut, cut_min, cut_max, NULL, NULL);
            double x;
            fill_random_double(&x, 1, seed);
            cut = cut_min + (x+0.5)*(cut_max - cut_min);
        }

        void random_slice(unsigned int seed, SphereSlice* slice, int id_cut,
            double &cut_begin, double &cut_end, double &cut_min, double &cut_max) {
            // Do a solve_sphere, to know over which range we can make a
            // disc-like slice that still intersects with the sphere.
            slice->solve_sphere(id_cut, cut_min, cut_max, NULL, NULL);
            EXPECT_LT(cut_min, cut_max);
            // Select two cut positions between cut_min and cut_max
            double cuts[2];
            fill_random_double(cuts, 2, seed);
            if (cuts[1] < cuts[0]) {
                double tmp = cuts[0];
                cuts[0] = cuts[1];
                cuts[1] = tmp;
            }
            cut_begin = (cuts[0] + 0.5)*(cut_max - cut_min) + cut_min;
            cut_end = (cuts[1] + 0.5)*(cut_max - cut_min) + cut_min;
        }

};


TEST_F(SphereSliceTest, domain) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    EXPECT_THROW(slice.solve_range(-1, begin, end), std::domain_error);
    EXPECT_THROW(slice.solve_range(3, begin, end), std::domain_error);
    EXPECT_THROW(slice.set_cut_begin_end(-1, 1, 2), std::domain_error);
    EXPECT_THROW(slice.set_cut_begin_end(2, 1, 2), std::domain_error);
    EXPECT_THROW(slice.set_cut_begin_end(0, 1, -2), std::domain_error);
    EXPECT_THROW(slice.set_cut_begin_end(1, 2, 1), std::domain_error);
}

TEST_F(SphereSliceTest, sphere_example1) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    slice.solve_range(0, begin, end);
    EXPECT_DOUBLE_EQ(-4.6, begin);
    EXPECT_DOUBLE_EQ(5.4, end);
}

TEST_F(SphereSliceTest, disc_slice_example1) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    slice.set_cut_begin_end(0, 0.0, 1.0);
    slice.solve_range(1, begin, end);
    EXPECT_DOUBLE_EQ(-7.0, begin);
    EXPECT_DOUBLE_EQ(3.0, end);
}

TEST_F(SphereSliceTest, disc_slice_example2) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    slice.set_cut_begin_end(0, 3.4, 4.4);
    slice.solve_range(1, begin, end);
    EXPECT_DOUBLE_EQ(-6.0, begin);
    EXPECT_DOUBLE_EQ(2.0, end);
}

TEST_F(SphereSliceTest, disc_slice_example3) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    slice.set_cut_begin_end(0, -3.6, -2.6);
    slice.solve_range(1, begin, end);
    EXPECT_DOUBLE_EQ(-6.0, begin);
    EXPECT_DOUBLE_EQ(2.0, end);
}

TEST_F(SphereSliceTest, solve_circle_ortho) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    slice.solve_circle(1, 0, 3.4, begin, end, NULL, NULL);
    EXPECT_DOUBLE_EQ(-6.0, begin);
    EXPECT_DOUBLE_EQ(2.0, end);
}

TEST_F(SphereSliceTest, solve_circle_angle) {
    easy_normals[3] = 1.0;
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    slice.solve_circle(1, 0, 3.4, begin, end, NULL, NULL);
    EXPECT_DOUBLE_EQ(-2.6, begin);
    EXPECT_DOUBLE_EQ(5.4, end);
}

TEST_F(SphereSliceTest, bar_slice_example1) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    slice.set_cut_begin_end(0, 0.0, 1.0);
    slice.set_cut_begin_end(1, -2.2, 1.0);
    slice.solve_range(2, begin, end);
    EXPECT_DOUBLE_EQ(6.0, begin);
    EXPECT_DOUBLE_EQ(-4.0, end);
}

TEST_F(SphereSliceTest, bar_slice_example2) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    slice.set_cut_begin_end(0, 3.4, 4.0);
    slice.set_cut_begin_end(1, -2.2, 1.0);
    slice.solve_range(2, begin, end);
    EXPECT_DOUBLE_EQ(5.0, begin);
    EXPECT_DOUBLE_EQ(-3.0, end);
}

TEST_F(SphereSliceTest, bar_slice_example3) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    slice.set_cut_begin_end(0, -5.0, -3.0);
    slice.set_cut_begin_end(1, -2.2, 1.0);
    slice.solve_range(2, begin, end);
    EXPECT_DOUBLE_EQ(5.0, begin);
    EXPECT_DOUBLE_EQ(-3.0, end);
}

TEST_F(SphereSliceTest, bar_slice_example4) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    slice.set_cut_begin_end(0, 0.0, 1.0);
    slice.set_cut_begin_end(1, -6.2, -5.0);
    slice.solve_range(2, begin, end);
    EXPECT_DOUBLE_EQ(5.0, begin);
    EXPECT_DOUBLE_EQ(-3.0, end);
}

TEST_F(SphereSliceTest, bar_slice_example5) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, 5.0);
    double begin, end;
    // 0.4, -2.0, 1.0
    slice.set_cut_begin_end(0, 0.0, 1.0);
    slice.set_cut_begin_end(1, 1.0, 2.0);
    slice.solve_range(2, begin, end);
    EXPECT_DOUBLE_EQ(5.0, begin);
    EXPECT_DOUBLE_EQ(-3.0, end);
}

TEST_F(SphereSliceTest, bar_slice_example6) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, sqrt(14.0));
    double begin, end;
    // 0.4, -2.0, 1.0
    slice.set_cut_begin_end(0, 1.4, 2.0);
    slice.set_cut_begin_end(1, 0.0, 1.0);
    slice.solve_range(2, begin, end);
    EXPECT_DOUBLE_EQ(2.0, begin);
    EXPECT_DOUBLE_EQ(0.0, end);
}

TEST_F(SphereSliceTest, bar_slice_example7) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, sqrt(14.0));
    double begin, end;
    slice.set_cut_begin_end(0, 1.4, 2.0);
    slice.set_cut_begin_end(1, -5.2, -4.0);
    slice.solve_range(2, begin, end);
    EXPECT_DOUBLE_EQ(2.0, begin);
    EXPECT_DOUBLE_EQ(0.0, end);
}

TEST_F(SphereSliceTest, bar_slice_example8) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, sqrt(14.0));
    double begin, end;
    slice.set_cut_begin_end(0, -2.0, -0.6);
    slice.set_cut_begin_end(1, 0.0, 1.0);
    slice.solve_range(2, begin, end);
    EXPECT_DOUBLE_EQ(2.0, begin);
    EXPECT_DOUBLE_EQ(0.0, end);
}

TEST_F(SphereSliceTest, bar_slice_example9) {
    SphereSlice slice = SphereSlice(my_center, easy_normals, sqrt(14.0));
    double begin, end;
    slice.set_cut_begin_end(0, -2.0, -0.6);
    slice.set_cut_begin_end(1, -5.2, -4.0);
    slice.solve_range(2, begin, end);
    EXPECT_DOUBLE_EQ(2.0, begin);
    EXPECT_DOUBLE_EQ(0.0, end);
}

TEST_F(SphereSliceTest, solve_sphere_random) {
    for (int irep=0; irep < 100; irep++) {
        // Test parameters
        double rcut = (irep+1)*0.1;
        double center[3];
        double normals[9];
        SphereSlice* slice = create_random_problem(irep, rcut, center, normals);

        // Do a solve_sphere
        double begin, end;
        double point_begin[3];
        double point_end[3];
        slice->solve_sphere(0, begin, end, point_begin, point_end);

        // Check consistency begin, point_begin
        EXPECT_NEAR(begin, vec3::dot(point_begin, normals), 1e-10);

        // Check consistency begin, point_end
        EXPECT_NEAR(end, vec3::dot(point_end, normals), 1e-10);

        // Check consistency when not using point_begin, point_end.
        double begin_bis, end_bis;
        slice->solve_sphere(0, begin_bis, end_bis, NULL, NULL);
        EXPECT_DOUBLE_EQ(begin, begin_bis);
        EXPECT_DOUBLE_EQ(end, end_bis);

        for (int ipoint=0; ipoint < 1000; ipoint++) {
            // Random point
            double point[3];
            double norm;
            random_point(ipoint, rcut, center, point, norm);

            // If the point is in the sphere, test if reduced coordinate falls
            // in the range [begin,end].
            double proj = vec3::dot(point, normals);
            if (norm < rcut) {
                EXPECT_LE(begin, proj);
                EXPECT_GE(end, proj);
            }
        }
    }
}

TEST_F(SphereSliceTest, solve_range_0_random) {
    for (int irep=0; irep < 100; irep++) {
        // Test parameters
        double rcut = (irep+1)*0.1;
        double center[3];
        double normals[9];
        SphereSlice* slice = create_random_problem(irep, rcut, center, normals);

        // Do a solve_sphere
        double begin, end;
        slice->solve_sphere(0, begin, end, NULL, NULL);

        // Do a solve_range_zero
        double begin_bis, end_bis;
        slice->solve_range_0(begin_bis, end_bis);

        // Should be the same
        EXPECT_EQ(begin, begin_bis);
        EXPECT_EQ(end, end_bis);
    }
}

TEST_F(SphereSliceTest, solve_circle_random) {
    for (int irep=0; irep < 100; irep++) {
        // Test parameters
        double rcut = (irep+1)*0.1;
        double center[3];
        double normals[9];
        SphereSlice* slice = create_random_problem(irep, rcut, center, normals);

        // Name some normals for convenience
        double* cut_normal = normals;
        double* axis = normals + 3;

        // Select randomized place to cut the sphere
        double cut, cut_min, cut_max;
        random_cut(irep+12345, slice, 0, cut, cut_min, cut_max);

        // Actual computation
        double begin, end;
        double point_begin[3];
        double point_end[3];
        bool exists = slice->solve_circle(1, 0, cut, begin, end, point_begin, point_end);

        // It should have worked...
        EXPECT_TRUE(exists);

        // Points must be on sphere...
        EXPECT_NEAR(rcut, vec3::distance(center, point_begin), 1e-10);
        EXPECT_NEAR(rcut, vec3::distance(center, point_end), 1e-10);

        // Points should be consistent with begin and end
        EXPECT_NEAR(begin, vec3::dot(point_begin, axis), 1e-10);
        EXPECT_NEAR(end, vec3::dot(point_end, axis), 1e-10);

        // Points should be on the cut
        EXPECT_NEAR(cut, vec3::dot(point_begin, cut_normal), 1e-10);
        EXPECT_NEAR(cut, vec3::dot(point_end, cut_normal), 1e-10);

        // Check that cut_normal, axis and (point_end - point_begin) lie in
        // the same plane.
        double delta[3];
        vec3::delta(point_end, point_begin, delta);
        EXPECT_NEAR(0.0, vec3::triple_product(cut_normal, axis, delta), 1e-10);

        // Check that cut_normal is orthogonal to delta.
        EXPECT_NEAR(0.0, vec3::dot(cut_normal, delta), 1e-10);

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
        EXPECT_NEAR(cut, vec3::dot(cut_normal, circle_center), 1e-10);

        // Scan the circle close to the solutions, where optimality is easily
        // tested.
        for (double angle=-0.01; angle < 0.01; angle += 0.001) {
            double other[3];
            double proj_axis;
            // Check end
            vec3::copy(circle_center, other);
            vec3::iadd(other, circle_basis0, cos(angle));
            vec3::iadd(other, circle_basis1, sin(angle));
            proj_axis = vec3::dot(other, axis);
            EXPECT_GE(end+1e-10, proj_axis);
            // Check begin
            vec3::copy(circle_center, other);
            vec3::iadd(other, circle_basis0, -cos(angle));
            vec3::iadd(other, circle_basis1, -sin(angle));
            proj_axis = vec3::dot(other, axis);
            EXPECT_LE(begin-1e-10, proj_axis);
        }

        // Call without point_* arguments
        double begin_bis, end_bis;
        exists = slice->solve_circle(1, 0, cut, begin_bis, end_bis, NULL, NULL);

        // It should have worked...
        EXPECT_TRUE(exists);

        // Result should be the same
        EXPECT_DOUBLE_EQ(begin, begin_bis);
        EXPECT_DOUBLE_EQ(end, end_bis);

        // Clean up
        delete slice;
    }
}

TEST_F(SphereSliceTest, solve_range_1_random) {
    for (int irep=0; irep < 100; irep++) {
        // Test parameters
        double rcut = (irep+1)*0.1;
        double center[3];
        double normals[9];
        SphereSlice* slice = create_random_problem(irep, rcut, center, normals);

        // Name some normals for convenience
        double* cut_normal = normals;
        double* axis = normals + 3;

        double cut_begin, cut_end;
        double cut_min, cut_max;
        random_slice(irep*2+1, slice, 0, cut_begin, cut_end, cut_min, cut_max);

        // Do the computation
        double axis_begin, axis_end;
        slice->set_cut_begin_end(0, cut_begin, cut_end);
        slice->solve_range_1(axis_begin, axis_end);

        // Check if the solution of solve_range_1 is better or as good as the
        // solutions of the separate parts
        double axis_begin0, axis_end0;
        double exists;
        exists = slice->solve_circle(1, 0, cut_begin, axis_begin0, axis_end0, NULL, NULL);
        EXPECT_TRUE(exists);
        EXPECT_LE(axis_begin, axis_begin0);
        EXPECT_GE(axis_end, axis_end0);
        double axis_begin1, axis_end1;
        exists = slice->solve_circle(1, 0, cut_end, axis_begin1, axis_end1, NULL, NULL);
        EXPECT_TRUE(exists);
        EXPECT_LE(axis_begin, axis_begin1);
        EXPECT_GE(axis_end, axis_end1);
        double point_begin[3];
        double point_end[3];
        // If the sphere solution is in the proper range, it is the solution
        double axis_begin_sphere, axis_end_sphere;
        slice->solve_sphere(1, axis_begin_sphere, axis_end_sphere, point_begin, point_end);
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
        for (int ipoint=0; ipoint < 100; ipoint++) {
            // Random point
            double point[3];
            double norm;
            random_point(ipoint, rcut, center, point, norm);
            ASSERT_DOUBLE_EQ(norm, vec3::distance(point, center));
            if (norm < rcut) {
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
                }
            }
        }
        delete slice;
    }
}


TEST_F(SphereSliceTest, solve_line_random) {
    for (int irep=0; irep < 100; irep++) {
        // Test parameters
        double rcut = (irep+1)*0.1;
        double center[3];
        double normals[9];
        SphereSlice* slice = create_random_problem(irep, rcut, center, normals);

        // Name some normals for convenience
        double* cut0_normal = normals;
        double* cut1_normal = normals + 3;
        double* axis = normals + 6;

        // Select randomized places to cut the sphere
        double cut0, cut0_min, cut0_max;
        double cut1, cut1_min, cut1_max;
        random_cut(irep+12345, slice, 0, cut0, cut0_min, cut0_max);
        random_cut(irep*2+114, slice, 1, cut1, cut1_min, cut1_max);

        // Actual computation
        double begin, end;
        double point_begin[3];
        double point_end[3];
        bool exists = slice->solve_line(2, 0, 1, cut0, cut1, begin, end, point_begin, point_end);

        // It should have worked...
        EXPECT_TRUE(exists);

        // Points must be on sphere...
        EXPECT_NEAR(rcut, vec3::distance(center, point_begin), 1e-10);
        EXPECT_NEAR(rcut, vec3::distance(center, point_end), 1e-10);

        // Points should be consistent with begin and end
        EXPECT_NEAR(begin, vec3::dot(point_begin, axis), 1e-10);
        EXPECT_NEAR(end, vec3::dot(point_end, axis), 1e-10);

        // Points should be on the cut
        EXPECT_NEAR(cut0, vec3::dot(point_begin, cut0_normal), 1e-10);
        EXPECT_NEAR(cut0, vec3::dot(point_end, cut0_normal), 1e-10);
        EXPECT_NEAR(cut1, vec3::dot(point_begin, cut1_normal), 1e-10);
        EXPECT_NEAR(cut1, vec3::dot(point_end, cut1_normal), 1e-10);

        // Check that cut0_normal and cut1_normal are orthogonal to delta.
        double delta[3];
        vec3::delta(point_end, point_begin, delta);
        EXPECT_NEAR(0.0, vec3::dot(cut0_normal, delta), 1e-10);
        EXPECT_NEAR(0.0, vec3::dot(cut1_normal, delta), 1e-10);

        // Get the line center and radius
        double line_center[3];
        vec3::copy(point_begin, line_center);
        vec3::iadd(line_center, point_end);
        vec3::iscale(line_center, 0.5);
        // TODO: remove: double line_radius = vec3::distance(point_begin, point_end);

        // Check the line center
        EXPECT_NEAR(cut0, vec3::dot(cut0_normal, line_center), 1e-10);
        EXPECT_NEAR(cut1, vec3::dot(cut1_normal, line_center), 1e-10);

        // Randomly sample the sphere surface close to the solution. These
        // random samples should all be worse.
        for (int isample=0; isample<100; isample++) {
            double other[3];
            double frac_axis;
            // Check end
            fill_random_double(other, 3, isample+irep, rcut*1e-4);
            vec3::iadd(other, point_end);
            vec3::iadd(other, center, -1);
            vec3::iscale(other, 1.0/vec3::norm(other));
            vec3::iadd(other, center);
            frac_axis = vec3::dot(other, axis);
            EXPECT_GE(end+1e-10, frac_axis);
            // Check begin
            fill_random_double(other, 3, isample+irep, rcut*1e-4);
            vec3::iadd(other, point_begin);
            vec3::iadd(other, center, -1);
            vec3::iscale(other, 1.0/vec3::norm(other));
            vec3::iadd(other, center);
            frac_axis = vec3::dot(other, axis);
            EXPECT_LE(begin-1e-10, frac_axis);
        }

        // Call without point_* arguments
        double begin_bis, end_bis;
        exists = slice->solve_line(2, 0, 1, cut0, cut1, begin_bis, end_bis, NULL, NULL);

        // It should have worked...
        EXPECT_TRUE(exists);

        // Result should be the same
        EXPECT_DOUBLE_EQ(begin, begin_bis);
        EXPECT_DOUBLE_EQ(end, end_bis);

        // Clean up
        delete slice;
    }
}


TEST_F(SphereSliceTest, solve_range_2_random) {
    for (int irep=0; irep < 100; irep++) {
        // Test parameters
        double rcut = (irep+1)*0.1;
        double center[3];
        double normals[9];
        SphereSlice* slice = create_random_problem(irep, rcut, center, normals);

        // Name some normals for convenience
        double* cut0_normal = normals;
        double* cut1_normal = normals + 3;
        double* axis = normals + 6;

        // Take two random slices to form a random bar
        double cut0_begin, cut0_end, cut0_min, cut0_max;
        random_slice(irep*3+1, slice, 0, cut0_begin, cut0_end, cut0_min, cut0_max);
        double cut1_begin, cut1_end, cut1_min, cut1_max;
        random_slice(irep*3+2, slice, 1, cut1_begin, cut1_end, cut1_min, cut1_max);

        // Do the computation
        double axis_begin, axis_end;
        slice->set_cut_begin_end(0, cut0_begin, cut0_end);
        slice->set_cut_begin_end(1, cut1_begin, cut1_end);
        slice->solve_range_2(axis_begin, axis_end);

        // Check if the solution of solve_range_2 is better or as good as the
        // solutions of the 9 separate parts (A-I).
        double exists;

        // * case A: cut0_begin  cut1_begin
        double axis_begin_a, axis_end_a;
        double frac1_begin_a, frac1_end_a;
        double point_begin_a[3], point_end_a[3];
        exists = slice->solve_line(2, 0, 1, cut0_begin, cut1_begin, axis_begin_a, axis_end_a, point_begin_a, point_end_a);
        EXPECT_TRUE(exists);
        EXPECT_LE(axis_begin, axis_begin_a);
        EXPECT_GE(axis_end, axis_end_a);

        // * case B: cut0_begin  cut1_end
        double axis_begin_b, axis_end_b;
        double frac1_begin_b, frac1_end_b;
        double point_begin_b[3], point_end_b[3];
        exists = slice->solve_line(2, 0, 1, cut0_begin, cut1_end, axis_begin_b, axis_end_b, point_begin_b, point_end_b);
        EXPECT_TRUE(exists);
        EXPECT_LE(axis_begin, axis_begin_b);
        EXPECT_GE(axis_end, axis_end_b);

        // * case C: cut0_end    cut1_begin
        double axis_begin_c, axis_end_c;
        double frac1_begin_c, frac1_end_c;
        double point_begin_c[3], point_end_c[3];
        exists = slice->solve_line(2, 0, 1, cut0_end, cut1_begin, axis_begin_c, axis_end_c, point_begin_c, point_end_c);
        EXPECT_TRUE(exists);
        EXPECT_LE(axis_begin, axis_begin_c);
        EXPECT_GE(axis_end, axis_end_c);

        // * case D: cut0_end    cut1_end
        double axis_begin_d, axis_end_d;
        double frac1_begin_d, frac1_end_d;
        double point_begin_d[3], point_end_d[3];
        exists = slice->solve_line(2, 0, 1, cut0_end, cut1_end, axis_begin_d, axis_end_d, point_begin_d, point_end_d);
        EXPECT_TRUE(exists);
        EXPECT_LE(axis_begin, axis_begin_d);
        EXPECT_GE(axis_end, axis_end_d);

        // * case E: cut0_begin
        double axis_begin_e, axis_end_e;
        double frac1_begin_e, frac1_end_e;
        double point_begin_e[3], point_end_e[3];
        exists = slice->solve_circle(2, 0, cut0_begin, axis_begin_e, axis_end_e, point_begin_e, point_end_e);
        EXPECT_TRUE(exists);
        frac1_begin_e = vec3::dot(cut1_normal, point_begin_e);
        if ((frac1_begin_e > cut1_begin) && (frac1_begin_e < cut1_end))
            EXPECT_LE(axis_begin, axis_begin_e);
        frac1_end_e = vec3::dot(cut1_normal, point_end_e);
        if ((frac1_end_e > cut1_begin) && (frac1_end_e < cut1_end))
            EXPECT_GE(axis_end, axis_end_e);

        // * case F: cut0_end
        double axis_begin_f, axis_end_f;
        double frac1_begin_f, frac1_end_f;
        double point_begin_f[3], point_end_f[3];
        exists = slice->solve_circle(2, 0, cut0_end, axis_begin_f, axis_end_f, point_begin_f, point_end_f);
        EXPECT_TRUE(exists);
        frac1_begin_f = vec3::dot(cut1_normal, point_begin_f);
        if ((frac1_begin_f > cut1_begin) && (frac1_begin_f < cut1_end))
            EXPECT_LE(axis_begin, axis_begin_f);
        frac1_end_f = vec3::dot(cut1_normal, point_end_f);
        if ((frac1_end_f > cut1_begin) && (frac1_end_f < cut1_end))
            EXPECT_GE(axis_end, axis_end_f);

        // * case G: cut1_begin
        double axis_begin_g, axis_end_g;
        double frac0_begin_g, frac0_end_g;
        double point_begin_g[3], point_end_g[3];
        exists = slice->solve_circle(2, 0, cut1_begin, axis_begin_g, axis_end_g, point_begin_g, point_end_g);
        EXPECT_TRUE(exists);
        frac0_begin_g = vec3::dot(cut0_normal, point_begin_g);
        if ((frac0_begin_g > cut0_begin) && (frac0_begin_g < cut0_end))
            EXPECT_LE(axis_begin, axis_begin_g);
        frac0_end_g = vec3::dot(cut0_normal, point_end_g);
        if ((frac0_end_g > cut0_begin) && (frac0_end_g < cut0_end))
            EXPECT_GE(axis_end, axis_end_g);

        // * case H: cut1_end
        double axis_begin_h, axis_end_h;
        double frac0_begin_h, frac0_end_h;
        double point_begin_h[3], point_end_h[3];
        exists = slice->solve_circle(2, 0, cut1_end, axis_begin_h, axis_end_h, point_begin_h, point_end_h);
        EXPECT_TRUE(exists);
        frac0_begin_h = vec3::dot(cut0_normal, point_begin_h);
        if ((frac0_begin_h > cut0_begin) && (frac0_begin_h < cut0_end))
            EXPECT_LE(axis_begin, axis_begin_h);
        frac0_end_h = vec3::dot(cut0_normal, point_end_h);
        if ((frac0_end_h > cut0_begin) && (frac0_end_h < cut0_end))
            EXPECT_GE(axis_end, axis_end_h);

        // * case I: (none)
        //   If the sphere solution is in the proper range, it is the solution.
        double axis_begin_i, axis_end_i;
        double frac0_begin_i, frac0_end_i;
        double frac1_begin_i, frac1_end_i;
        double point_begin_i[3], point_end_i[3];

        slice->solve_sphere(2, axis_begin_i, axis_end_i, point_begin_i, point_end_i);
        frac0_begin_i = vec3::dot(cut0_normal, point_begin_i);
        frac1_begin_i = vec3::dot(cut1_normal, point_begin_i);
        if ((frac0_begin_i > cut0_begin) && (frac0_begin_i < cut0_end) &&
            (frac1_begin_i > cut1_begin) && (frac1_begin_i < cut1_end))
            EXPECT_EQ(axis_begin, axis_begin_i);
        frac0_end_i = vec3::dot(cut0_normal, point_end_i);
        frac1_end_i = vec3::dot(cut1_normal, point_end_i);
        if ((frac0_end_i > cut0_begin) && (frac0_end_i < cut0_end) &&
            (frac1_end_i > cut1_begin) && (frac1_end_i < cut1_end))
            EXPECT_EQ(axis_end, axis_end_i);

        // Test by generating lots of random points. If a point is in the bar,
        // its fractional coordinate should be between begin and end.
        for (int ipoint=0; ipoint < 100; ipoint++) {
            // Random point
            double point[3];
            double norm;
            random_point(ipoint, rcut, center, point, norm);
            ASSERT_DOUBLE_EQ(norm, vec3::distance(point, center));
            if (norm < rcut) {
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
                }
            }
        }
        delete slice;
    }
}
