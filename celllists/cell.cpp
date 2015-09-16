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


#include "celllists/cell.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "celllists/vec3.h"


namespace celllists {


Cell::Cell(const double* vecs, const int nvec): nvec_(nvec) {
  // check if nvec is sensible
  if ((nvec_ < 0) || (nvec_ > 3))
    throw std::domain_error("The number of cell vectors must be 0, 1, 2 or 3.");

  // copy the given _vecs and _nvec:
  std::copy(vecs, vecs + nvec_*3, vecs_);

  // compute the volume
  switch (nvec_) {
    case 0: {
      volume_ = NAN;
      break;
    }
    case 1: {
      volume_ = vec3::norm(vecs_);
      break;
    }
    case 2: {
      double tmp = vec3::dot(vecs_, vecs_ + 3);
      tmp = vec3::normsq(vecs_)*vec3::normsq(vecs_ + 3) - tmp*tmp;
      volume_ = (tmp > 0.0) ? sqrt(tmp) : 0.0;
      break;
    }
    case 3: {
      volume_ = fabs(vec3::triple(vecs_, vecs_ + 3, vecs_ + 6));
      break;
    }
  }

  // If the volume is zero and nvec > 0, raise an error. In this case, the
  // reciprocal cell vectors can not be computed.
  if (volume_ == 0.0)
    throw singular_cell_vectors("The cell vectors are degenerate");
  gvolume_ = 1.0/volume_;

  // complete the list of vecs in case nvec < 3
  switch (nvec_) {
    case 0: {
      // Just put in the identity matrix.
      std::fill(vecs_, vecs_ + 9, 0.0);
      vecs_[0] = 1.0;
      vecs_[4] = 1.0;
      vecs_[8] = 1.0;
      break;
    }
    case 1: {
      // Add two vecs that are orthogonal to the given vec, orthogonal
      // to each other and normalized. The three vectors will be
      // right-handed.
      // 1) find the component of the given vector with the smallest
      // absolute value
      int ismall = 0;
      if (fabs(vecs_[1]) < fabs(vecs_[0])) {
        ismall = 1;
        if (fabs(vecs_[2]) <= fabs(vecs_[1])) {
          ismall = 2;
        }
      } else if (fabs(vecs_[2]) < fabs(vecs_[0])) {
        ismall = 2;
      }
      // 2) store a temporary vector in position 3
      std::fill(vecs_ + 6, vecs_ + 9, 0.0);
      vecs_[ismall + 6] = 1.0;
      // 3) compute the cross product of vector 3 and 1
      vec3::cross(vecs_ + 6, vecs_, vecs_ + 3);
      // 4) normalize
      double norm = vec3::norm(vecs_ + 3);
      vec3::iscale(vecs_ + 3, 1.0/norm);
      // the rest is done in case 2, so no break here!
    }
    case 2: {
      // Add one vec that is normalized and orthogonal to the two given
      // vecs. The three vectors will be right-handed.
      // 1) compute the cross product of vector 1 and 2
      vec3::cross(vecs_, vecs_ + 3, vecs_ + 6);
      // 2) normalize
      vec3::iscale(vecs_ + 6, 1.0/vec3::norm(vecs_ + 6));
    }
  }

  // Now we assume that vecs contains a set of three well-behaved
  // non-degenerate vectors. Cramer's rule is used to compute the reciprocal
  // space vectors. This is fairly ugly in terms of numerical stability but
  // it keeps things simple.
  vec3::cross(vecs_ + 3, vecs_ + 6, gvecs_);
  vec3::cross(vecs_ + 6, vecs_, gvecs_ + 3);
  vec3::cross(vecs_, vecs_ + 3, gvecs_ + 6);
  // inverse of determinant
  double denom = 1.0/vec3::dot(gvecs_, vecs_);
  // inverse
  vec3::iscale(gvecs_, denom);
  vec3::iscale(gvecs_ + 3, denom);
  vec3::iscale(gvecs_ + 6, denom);

  // compute the spacings and the lengths of the cell vectors
  for (int ivec = 2; ivec >= 0; --ivec) {
    lengths_[ivec] = vec3::norm(vecs_ + 3*ivec);
    glengths_[ivec] = vec3::norm(gvecs_ + 3*ivec);
    spacings_[ivec] = 1.0/glengths_[ivec];
    gspacings_[ivec] = 1.0/lengths_[ivec];
  }
}


Cell* Cell::create_reciprocal() const {
  return new Cell(gvecs_, nvec_, vecs_, gvolume_, volume_, glengths_, lengths_,
                  gspacings_, spacings_);
}


Cell* Cell::create_subcell(const int* shape, const double* spacings, bool* pbc) {
  // Start by copying all three cell vectors, active or not.
  double new_vecs[9];
  std::copy(vecs_, vecs_ + 9, new_vecs);
  // Divide the 'active' vectors by the corresponding shape value.
  for (int ivec = 0; ivec < nvec_; ++ivec) {
    vec3::iscale(new_vecs + 3*ivec, 1.0/shape[ivec]);
    pbc[ivec] = true;
  }
  // Multiply the 'inactive' vectors by the corresponding spacing.
  for (int ivec = nvec_; ivec < 3; ++ivec) {
    vec3::iscale(new_vecs + 3*ivec, spacings[ivec - nvec_]);
    pbc[ivec] = false;
  }
  // Return the subcell, always 3D periodic
  return new Cell(new_vecs, 3);
}




const double* Cell::vec(const int ivec) const {
  if ((ivec < 0) || (ivec >= 3)) {
    throw std::domain_error("ivec must be 0, 1 or 2.");
  }
  return vecs_ + 3*ivec;
}


const double* Cell::gvec(const int ivec) const {
  if ((ivec < 0) || (ivec >= 3)) {
    throw std::domain_error("ivec must be 0, 1 or 2.");
  }
  return gvecs_ + 3*ivec;
}


void Cell::to_frac(const double* rcart, double* rfrac) const {
  // Transfrom to real-space fractional coordinates
  vec3::matvec(gvecs_, rcart, rfrac);
}


void Cell::to_cart(const double* rfrac, double* rcart) const {
  // Transfrom to real-space Cartesian coordinates
  vec3::tmatvec(vecs_, rfrac, rcart);
}


void Cell::iwrap(double* delta) const {
  double x;
  if (nvec_ == 0) return;
  // Compute the first fractional coordinates, subtract one half and ceil. The round
  // founction is intentionally not used here! The half-ways case is always up instead
  // of away from zero.
  x = ceil(vec3::dot(gvecs_, delta) - 0.5);
  vec3::iadd(delta, vecs_, -x);
  if (nvec_ == 1) return;
  // Compute the second fractional coordinates, subtract one half and ceil.
  x = ceil(vec3::dot(gvecs_ + 3, delta) - 0.5);
  vec3::iadd(delta, vecs_ + 3, -x);
  if (nvec_ == 2) return;
  // Compute the third fractional coordinates, subtract one half and ceil.
  x = ceil(vec3::dot(gvecs_ + 6, delta) - 0.5);
  vec3::iadd(delta, vecs_ + 6, -x);
}


void Cell::iadd_vec(double* delta, const int* coeffs) const {
  // Simply adds an linear combination of real cell vectors to delta.
  if (nvec_ == 0) return;
  vec3::iadd(delta, vecs_, coeffs[0]);
  if (nvec_ == 1) return;
  vec3::iadd(delta, vecs_ + 3, coeffs[1]);
  if (nvec_ == 2) return;
  vec3::iadd(delta, vecs_ + 6, coeffs[2]);
}


bool Cell::cubic() const {
  if (!cuboid()) return false;
  if (nvec_ < 2) return true;
  if (vecs_[0] != vecs_[4]) return false;
  if (nvec_ < 3) return true;
  if (vecs_[0] != vecs_[8]) return false;
  return true;
}


bool Cell::cuboid() const {
  if (nvec_ < 1) return true;
  if (vecs_[1] != 0.0) return false;
  if (vecs_[2] != 0.0) return false;
  if (nvec_ < 2) return true;
  if (vecs_[3] != 0.0) return false;
  if (vecs_[5] != 0.0) return false;
  if (nvec_ < 3) return true;
  if (vecs_[6] != 0.0) return false;
  if (vecs_[7] != 0.0) return false;
  return true;
}


int Cell::ranges_cutoff(const double* center, const double cutoff, int* ranges_begin,
    int* ranges_end) const {
  if (cutoff <= 0) {
    throw std::domain_error("cutoff must be strictly positive.");
  }
  double frac[3];
  int ncell = 1;
  to_frac(center, frac);
  for (int ivec = nvec_ - 1; ivec >= 0; --ivec) {
    // Use spacings between planes to find first plane before cutoff sphere and last
    // plane after cutoff sphere. To this end, we must divide cutoff by the spacing
    // between planes.
    double frac_cutoff = cutoff/spacings_[ivec];
    ranges_begin[ivec] = static_cast<int>(floor(frac[ivec] - frac_cutoff));
    ranges_end[ivec] = static_cast<int>(ceil(frac[ivec] + frac_cutoff));
    ncell *= (ranges_end[ivec] - ranges_begin[ivec]);
  }
  return ncell;
}


size_t Cell::bars_cutoff(const double* center, const double cutoff,
    const int* shape, const bool* pbc, std::vector<int>* bars) const {
  // Check arguments
  if (nvec_ == 0) {
    throw std::domain_error("The cell must be at least 1D periodic.");
  }
  if (cutoff <= 0) {
    throw std::domain_error("cutoff must be strictly positive.");
  }
  // For all the heavy work, a SphereSlice object is used that precomputes a lot.
  SphereSlice sphere_slice(center, gvecs_, cutoff);
  // Prefix is used to keep track of current bar indices while going into recursion.
  std::vector<int> prefix;
  // Compute bars and return the number of bars
  bars_cutoff_low(&sphere_slice, shape, pbc, &prefix, bars);
  return bars->size()/(nvec_ + 1);
}


Cell::Cell(const double* vecs, const int nvec, const double* gvecs,
    const double volume, const double gvolume,
    const double* lengths, const double* glengths,
    const double* spacings, const double* gspacings)
    : nvec_(nvec), volume_(volume), gvolume_(gvolume) {
  // Just copy everything to data members;
  std::copy(vecs, vecs + nvec_*3, vecs_);
  std::copy(gvecs, gvecs + nvec_*3, gvecs_);
  std::copy(lengths, lengths + nvec_, lengths_);
  std::copy(glengths, glengths + nvec_, glengths_);
  std::copy(spacings, spacings + nvec_, spacings_);
  std::copy(gspacings, gspacings + nvec_, gspacings_);
}


void Cell::bars_cutoff_low(SphereSlice* slice, const int* shape,
    const bool* pbc, std::vector<int>* prefix, std::vector<int>* bars) const {
  // Get the vector index for which the range is currently searched
  int ivec = static_cast<int>(prefix->size());
  // Use SphereSlice object to solve the hard of problem of finding begin and end.
  double begin_exact = 0.0;
  double end_exact = 0.0;
  slice->solve_range(ivec, &begin_exact, &end_exact);
  int begin = static_cast<int>(floor(begin_exact));
  int end = static_cast<int>(ceil(end_exact));
  // Truncate the begin-end range if there are non-periodic bounds.
  if (!pbc[ivec]) {
    if (begin < 0) begin = 0;
    if (end > shape[ivec]) end = shape[ivec];
  }

  if (ivec == nvec_ - 1) {
    // If we are dealing with the last recursion, just store the bar.
    for (const auto& i : *prefix)
      bars->push_back(i);
    bars->push_back(begin);
    bars->push_back(end);
  } else {
    // If this is not yet the last recursion, iterate over the range of integer
    // fractional coordinates, and go one recursion deeper in each iteration.
    for (int i = begin; i < end; ++i) {
      // Make sure the following recursion knows the indices of the current bar.
      prefix->push_back(i);
      // Make a new cut in the sphere slice.
      slice->set_cut_begin_end(ivec, i, i + 1);
      // Recursive call in which the remaining details of the bar/bars is/are solved.
      bars_cutoff_low(slice, shape, pbc, prefix, bars);
      // Remove the last element of prefix as it is no longer applicable.
      prefix->pop_back();
    }
  }
}


int smart_wrap(int i, const int shape, const bool pbc) {
  if (pbc) {
    i %= shape;
    if (i < 0) i += shape;
    return i;
  } else if (i < 0) {
    return -1;
  } else if (i >= shape) {
    return -1;
  } else {
    return i;
  }
}


}  // namespace celllists

// vim: textwidth=90 et ts=2 sw=2
