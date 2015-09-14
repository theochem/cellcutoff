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


#ifndef CELLLIST_VEC3_H_
#define CELLLIST_VEC3_H_

namespace vec3 {

double norm(const double* vec);
double normsq(const double* vec);
double dot(const double* vec1, const double* vec2);
double distance(const double* vec1, const double* vec2);
double triple(const double* vec1, const double* vec2, const double* vec3);
void cross(const double* vec1, const double* vec2, double* vec3);
void iscale(double* vec, double scale);
void copy(const double* source, double* dest);
void iadd(double* output, const double* term, double scale=1);
void delta(const double* begin, const double* end, double* output);

}

#endif
