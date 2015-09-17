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


#ifndef CELLLISTS_VEC3_H_
#define CELLLISTS_VEC3_H_


namespace celllists {
namespace vec3 {


/* Some notes:

    - Usually, the output argument is last. The in-place functions do the opposite with
      the input-output argument.

    - This is intentionally not made object-oriented to allow for simple pointer
      arithmetics, e.g. norm(rvecs+3);

*/


//! Computes the norm of a vector
double norm(const double* vec);

//! Computes the squared norm of a vector
double normsq(const double* vec);

//! Computes the dot product of two vectors
double dot(const double* vec1, const double* vec2);

//! Computes the distance between two vectors
double distance(const double* vec1, const double* vec2);

//! Computes the triple product of the three given vectors
double triple(const double* vec1, const double* vec2, const double* vec3);

//! Computes the cross product of the first two vectors and assigns result to the third
void cross(const double* vec1, const double* vec2, double* vec3);

//! Scale the vector in-place
void iscale(double* vec, double scale);

//! Copy contents of first argument to the second
void copy(const double* source, double* dest);

//! Add to the first argument the second argument
void iadd(double* output, const double* term);

//! Add to the first argument a rescaled second argument
void iadd(double* output, const double* term, double scale);

//! Construct a relative vector
void delta(const double* begin, const double* end, double* output);

//! Compute a matrix-vector product
void matvec(const double* mat, const double* vec, double* output);

//! Compute a matrix^T-vector product
void tmatvec(const double* mat, const double* vec, double* output);


#endif  // CELLLISTS_VEC3_H_


}  // namespace vec3
}  // namespace celllists


// vim: textwidth=90 et ts=2 sw=2
