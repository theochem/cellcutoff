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

/** @file */

#ifndef CELLLIST_CELL_H_
#define CELLLIST_CELL_H_


/** @brief
        Abstract base class for cells

    Upon construction, an object of this class acts as a read-only representation of the
    cell. Reciprocal cell vectors, vector lengths and spacings between planes are computed
    immediately. All sorts of manipulations of fractional or Cartesian coordinates are
    supported.

    Note that this class is specific for 3D systems, even though lower-dimensional
    periodic boundary conditions may be supported by some subclasses. In case of 1D or 2D
    PBC, the cell vectors are internally extended with orthonormal basis vectors to
    guarantee an invertible transformation between Cartesian and fractional coordinates.
    In that case, the fractional coordinates are actually also Cartesian coordinates in
    directions orthogonal to the available cell vectors.
 */
class BaseCell {
    public:
        /** @brief
                Wrap a (relative) vector back into the cell [-0.5, 0.5[.

            @param delta
                A pointer to 3 doubles with the (relative) vector. It will be
                modified in-place.

            After calling the wrap method, the fractional coordinates of delta will be
            in the range [-0.5, 0.5[.

            This is an approximate implementation of the minimum image convention that
            sometimes fails in very skewed cells, i.e. the wrapped vector is not always
            the shortest relative vector between a reference point and all of its periodic
            images. For more details see:
            http://scicomp.stackexchange.com/questions/3107/minimum-image-convention-for-triclinic-unit-cell
        */
        virtual void wrap(double* delta) const = 0;

        /** @brief
                Convert Cartesian real-space coordinates to fractional coordinates.

            @param cart
                A pointer to 3 doubles containing the input real-space vector.

            @param frac
                A pointer to 3 doubles in which the output is written.

            This effectively computes the dot products of the Cartesian vector with the
            reciprocal cell vectors.
         */
        virtual void to_frac(const double* cart, double* frac) const = 0;

        /** @brief
                Convert fractional coordinates to Cartesian real-space coordinates.

            @param frac
                A pointer to 3 doubles containing the input fractional
                coordinates.

            @param cart
                A pointer to 3 doubles to which the output is written

            This effectively computes a linear combination of real-space cell vectors.
         */
        virtual void to_cart(const double* frac, double* cart) const = 0;

        /** @brief
                Construct a linear combination of reciprocal cell vectors.

            @param coeffs
                A pointer to 3 doubles containing the coefficients for the
                linear combination.

            @param gvec
                A pointer to 3 doubles to which the output is written
         */
        virtual void g_lincomb(const double* coeffs, double* gvec) const = 0;

        /** @brief
                Compute the dot products of fractional coordinates with each cell vector.

            @param frac
                A pointer to 3 doubles containing the input fractional
                coordinates.

            @param dots
                A pointer to 3 doubles to which the output is written.
         */
        virtual void dot_rvecs(const double* frac, double* dots) const = 0;

        /** @brief
                Add a linear combination of cell vectors to delta.

            @param delta
                A pointer to 3 doubles for the real-space vector to which the
                linear combination is added in-place.

            @param coeffs
                A pointer to 3 doubles with the coefficients of the linear
                combination.
         */
        virtual void add_rvec(double* delta, const long* coeffs) const = 0;


        //! Returns the number of periodic dimensions.
        virtual int get_nvec() const = 0;
        //! Returns the volume (or area or length) of the cell.
        virtual double get_volume() const = 0;
        //! Returns the spacing between the i-th real-space crystal plane
        virtual double get_rspacing(int i) const = 0;
        //! Returns the spacing between the i-th reciprocal crystal plane
        virtual double get_gspacing(int i) const = 0;
        //! Returns the length of the i-th real-space cell vector
        virtual double get_rlength(int i) const = 0;
        //! Returns the length of the i-th reciprocal cell vector
        virtual double get_glength(int i) const = 0;

        /** @brief
                Get the ranges of cells within a cutoff radius.

            @param center
                A pointer to 3 doubles that specify the center of the cutoff
                sphere in real-space.

            @param rcut
                The cutoff radius.

            @param ranges_begin
                A pointer to `nvec` longs to which the begin of each range of
                periodic images along a periodic boundary condition is written.

            @param ranges_end
                A pointer to `nvec` longs to which the end of each range of
                periodic images along a periodic boundary condition is written.
                Then end values are non-inclusive as in Python ranges.

            This function effectively defines a supercell that is guaranteed to
            enclose the cutoff sphere.
         */
        virtual void set_ranges_rcut(const double* center, double rcut, long*
            ranges_begin, long* ranges_end) const = 0;

        /** @brief
                Selects a list of cells inside a cutoff sphere.

            @return
                The number of periodic images inside the cutoff sphere.

            @param origin
                A pointer of 3 doubles with the origin of a supercell.

            @param center
                A pointer of 3 doubles with the center of the cutoff sphere.

            @param rcut
                The cutoff radius.

            @param ranges_begin
                As obtained with set_ranges_rcut().

            @param ranges_end
                As obtained with set_ranges_rcut().

            @param shape
                A pointer of 3 longs with the shape of the supercell.

            @param pbc
                A pointer to integer flags indicating the periodicity of the
                supercell along each cell vector.

            @param indices
                A sufficiently large pre-allocated output array to which the
                indexes of the selected periodic images are written. The number
                of rows is the product of the lengths of the ranges specified by
                ranges_begin and ranges_end. The number of columns equals `nvec`.
                The elements are stored in row-major order.
          */
        virtual long select_inside(const double* origin, const double* center,
            double rcut, const long* ranges_begin, const long* ranges_end, const long*
            shape, const long* pbc, long* indices) const = 0;
};


/** @brief
        3D/2D/1D (Tri)clinic cell and derived quantities.
 */
class GeneralCell : public BaseCell {
    private:
        double rvecs[9];       //!< real-space vectors,       one per row, row-major
        double gvecs[9];       //!< reciprocal-space vectors, one per row, row-major
        double rlengths[3];    //!< real-space vector lengths
        double glengths[3];    //!< reciprocal-space vector lengths
        double rspacings[3];   //!< spacing between real-space crystal planes
        double gspacings[3];   //!< spacing between reciprocal-space crystal planes
        double volume;         //!< volume (or area or length) of the cell
        const int nvec;        //!< number of defined cell vectors
    public:
        /** @brief
                Construct a Cell object.

            @param _rvecs
                A pointer to `3*nvec` doubles that represent the real-space vectors in
                row-major ordering. The vectors may not have a linear dependency.

            @param _nvec
                The number of cell vectors. This corresponds to the number of periodic
                dimensions of a unit cell. `nvec` must be 0, 1, 2 or 3.
        */
        GeneralCell(const double* _rvecs, int _nvec);

        void wrap(double* delta) const;
        void to_frac(const double* cart, double* frac) const;
        void to_cart(const double* frac, double* cart) const;
        void g_lincomb(const double* coeffs, double* gvec) const;
        void dot_rvecs(const double* frac, double* dots) const;
        void add_rvec(double* delta, const long* coeffs) const;

        int get_nvec() const {return nvec;};
        double get_volume() const {return volume;};
        double get_rspacing(int i) const;
        double get_gspacing(int i) const;
        double get_rlength(int i) const;
        double get_glength(int i) const;

        void set_ranges_rcut(const double* center, double rcut, long* ranges_begin,
            long* ranges_end) const;
        long select_inside(const double* origin, const double* center, double rcut,
            const long* ranges_begin, const long* ranges_end, const long* shape,
            const long* pbc, long* indices) const;
};

/**
    @brief
        A standardized modulo operation that works across all compilers.

    @param i
        The numerator of the integer division.

    @param shape
        The denominator of the integer division.

    @param pbc
        Whether periodic boundary conditions apply.

    @return
        `i % shape` (guaranteed to be positive) if `pbc` is non-zero. If `pbc` is zero,
        `i` is returned if it lies in `[0,shape[`, `-1` is returned otherwise.
 */
long smart_wrap(long i, long shape, long pbc);

#endif
