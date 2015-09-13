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

        // Do a solve_range_0, to know over which range we can cut
        double cut_min, cut_max;
        slice->solve_range_0(cut_min, cut_max);
        // TODO: use seed
        double cut = cut_min + (rand()/(RAND_MAX + 1.0))*(cut_max - cut_min);

        // Actual computation
        double begin, end;
        double point_begin[3];
        double point_end[3];
        bool exists = slice->solve_circle(1, 0, cut, begin, end, point_begin, point_end);

        // It should have worked...
        ASSERT_TRUE(exists);

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
        //double delta_center[3];
        //vec3::delta(center, circle_center, delta_center);
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
        ASSERT_TRUE(exists);

        // Result should be the same
        EXPECT_DOUBLE_EQ(begin, begin_bis);
        EXPECT_DOUBLE_EQ(end, end_bis);
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

        // Do a solve_range_0, to know over which range we can make a disc-like
        // slice that still intersects with the sphere.
        double cut_min, cut_max;
        slice->solve_range_0(cut_min, cut_max);
        EXPECT_LT(cut_min, cut_max);
        // Select two cut positions between cut_min and cut_max
        double cuts[2];
        fill_random_double(cuts, 2, irep*2+1);
        if (cuts[1] < cuts[0]) {
            double tmp = cuts[0];
            cuts[0] = cuts[1];
            cuts[1] = tmp;
        }
        double cut_begin = (cuts[0] + 0.5)*(cut_max - cut_min) + cut_min;
        double cut_end = (cuts[1] + 0.5)*(cut_max - cut_min) + cut_min;

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
                    //EXPECT_GE(axis_end, proj2);
                }
            }
        }

        delete slice;
    }
}
