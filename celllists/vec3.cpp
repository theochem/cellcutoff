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


#include <cmath>
#include "celllists/vec3.h"

namespace vec3 {

//! Computes the norm of a vector
double norm(const double* vec) {
    return sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

//! Computes the squared norm of a vector
double normsq(const double* vec) {
    return vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
}

//! Computes the dot product of two vectors
double dot(const double* vec1, const double* vec2) {
    return vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2];
}

//! Computes the distance between two vectors
double distance(const double* vec1, const double* vec2) {
    double d, sum;
    d = vec1[0] - vec2[0];
    sum = d*d;
    d = vec1[1] - vec2[1];
    sum += d*d;
    d = vec1[2] - vec2[2];
    sum += d*d;
    return sqrt(sum);
}

//! Computes the triple product of the three given vectors
double triple(const double* vec1, const double* vec2, const double* vec3) {
    return vec1[0]*(vec2[1]*vec3[2] - vec2[2]*vec3[1]) +
           vec1[1]*(vec2[2]*vec3[0] - vec2[0]*vec3[2]) +
           vec1[2]*(vec2[0]*vec3[1] - vec2[1]*vec3[0]);
}

//! Computes the cross product of the first two vectors and assigns result to the third
void cross(const double* vec1, const double* vec2, double* vec3) {
    vec3[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
    vec3[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
    vec3[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];
}

//! Scale the vector in-place
void iscale(double* vec, double scale) {
    vec[0] *= scale;
    vec[1] *= scale;
    vec[2] *= scale;
}

//! Copy contents of first argument to the second
void copy(const double* source, double* dest) {
    dest[0] = source[0];
    dest[1] = source[1];
    dest[2] = source[2];
}

//! Add to the first argument a rescaled second argument
void iadd(double* output, const double* term, double scale) {
    output[0] += scale*term[0];
    output[1] += scale*term[1];
    output[2] += scale*term[2];
}

//! Construct a relative vector
void delta(const double* begin, const double* end, double* output) {
    output[0] = end[0] - begin[0];
    output[1] = end[1] - begin[1];
    output[2] = end[2] - begin[2];
}

//! Compute a matrix-vector product
void matvec(const double* mat, const double* vec, double* output) {
    output[0] = mat[0]*vec[0] + mat[1]*vec[1] + mat[2]*vec[2];
    output[1] = mat[3]*vec[0] + mat[4]*vec[1] + mat[5]*vec[2];
    output[2] = mat[6]*vec[0] + mat[7]*vec[1] + mat[8]*vec[2];

}

//! Compute a matrix^T-vector product
void tmatvec(const double* mat, const double* vec, double* output) {
    output[0] = mat[0]*vec[0] + mat[3]*vec[1] + mat[6]*vec[2];
    output[1] = mat[1]*vec[0] + mat[4]*vec[1] + mat[7]*vec[2];
    output[2] = mat[2]*vec[0] + mat[5]*vec[1] + mat[8]*vec[2];

}

}
