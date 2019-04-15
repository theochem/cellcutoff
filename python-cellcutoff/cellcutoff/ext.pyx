# CellCutoff is a library for periodic boundary conditions and real-space cutoff calculations.
# Copyright (C) 2017 The CellCutoff Development Team
#
# This file is part of CellCutoff.
#
# CellCutoff is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# CellCutoff is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
# --
# cython: linetrace=True, embedsignature=True, language_level=3
'''Python wrapper for the CellCutoff library'''


import numpy as np
cimport numpy as np
np.import_array()

from libc.string cimport memcpy
from libcpp cimport bool

cimport cellcutoff.cell as cell
cimport cellcutoff.iterators


__all__ = ['Cell', 'ranges_cutoff', 'create_random_cell']


def check_array_arg(name, arg, expected_shape):
    if not arg.flags['C_CONTIGUOUS']:
        raise TypeError('Argument %s must be C_CONTIHUOUS.' % arg)
    for i, n in enumerate(expected_shape):
        if n >= 0 and arg.shape[i] != n:
            raise TypeError(('Axis %i of argument %s has length %i but while '
                             'expecting %i') % (i, name, arg.shape[i], n))


cdef class Cell(object):
    def __cinit__(self, np.ndarray[double, ndim=2] vecs=None, initvoid=False):
        """C-level constructor, allocates low-level object unless initvoid==True."""
        if initvoid:
            self._this = NULL
        elif vecs is None or (vecs.shape[0] == 0 and vecs.shape[1] == 3):
            self._this = new cell.Cell()
        else:
            check_array_arg('vecs', vecs, (-1, 3))
            if vecs.shape[0] > 3:
                raise TypeError("At most three cell vectors allows.")
            nvec = vecs.shape[0]
            self._this = new cell.Cell(&vecs[0,0], nvec)

    def __init__(self, np.ndarray[double, ndim=2] vecs=None):
        """Initialize a cell object.

        Parameters
        ----------
        vecs
            The cell vectors, each vector being one row. At most three rows
            are allowed and exactly three columbs are expected. None is
            equivalent np.zeros((0, 3)).

        """
        pass

    def __dealloc__(self):
        """Dealocate the cell object, also dealocates low-level object."""
        if self._this != NULL:
            del self._this

    def subcell(self, double threshold):
        """Generate an integer division of the current cell.

        Parameters
        ----------
        threshold
            The maximum distance between the subcell faces.

        Returns
        -------
        subcell
            The subcell.
        shape
            The number of times the subcell can be repeated along each cell
            vector to obtain the current cell.

        """
        cdef np.ndarray[int, ndim=1] shape = np.zeros(3, np.intc)
        cdef cell.Cell* cpp_cell = self._this.create_subcell(threshold, &shape[0])
        cdef Cell subcell = Cell.__new__(Cell, None, initvoid=True)
        subcell._this = cpp_cell
        return subcell, shape

    def reciprocal(self):
        """Return the reciprocal cell."""
        cdef cell.Cell* cpp_cell = self._this.create_reciprocal()
        cdef Cell cell = Cell.__new__(Cell, None, initvoid=True)
        cell._this = cpp_cell
        return cell

    property nvec:
        def __get__(self):
            """The dimensionality of the cell."""
            return self._this.nvec()

    property vecs:
        def __get__(self):
            """The cell vectors."""
            cdef np.ndarray[double, ndim=2] vecs = np.zeros((3, 3), float)
            memcpy(&vecs[0, 0], self._this.vecs(), sizeof(double)*9);
            return vecs

    property gvecs:
        def __get__(self):
            """The reciprocal cell vectors."""
            cdef np.ndarray[double, ndim=2] gvecs = np.zeros((3, 3), float)
            memcpy(&gvecs[0, 0], self._this.gvecs(), sizeof(double)*9);
            return gvecs

    property volume:
        def __get__(self):
            """The cell volume."""
            return self._this.volume()

    property gvolume:
        def __get__(self):
            """The reciprocal cell volume."""
            return self._this.gvolume()

    property lengths:
        def __get__(self):
            """The lengths of the cell vectors."""
            cdef np.ndarray[double, ndim=1] lengths = np.zeros(3, float)
            memcpy(&lengths[0], self._this.lengths(), sizeof(double)*3);
            return lengths

    property glengths:
        def __get__(self):
            """The lengths of the reciprocal cell vectors."""
            cdef np.ndarray[double, ndim=1] glengths = np.zeros(3, float)
            memcpy(&glengths[0], self._this.glengths(), sizeof(double)*3);
            return glengths

    property spacings:
        def __get__(self):
            """The distances between the cell planes."""
            cdef np.ndarray[double, ndim=1] spacings = np.zeros(3, float)
            memcpy(&spacings[0], self._this.spacings(), sizeof(double)*3);
            return spacings

    property gspacings:
        def __get__(self):
            """The distances between the reciprocal cell planes."""
            cdef np.ndarray[double, ndim=1] gspacings = np.zeros(3, float)
            memcpy(&gspacings[0], self._this.gspacings(), sizeof(double)*3);
            return gspacings

    property cubic:
        def __get__(self):
            """Is the cell cubic?"""
            return self._this.cubic()

    property cuboid:
        def __get__(self):
            """Is the cell cuboid?"""
            return self._this.cuboid()

    def to_frac(self, np.ndarray[double, ndim=1] cart not None):
        """Convert Cartesian cell vectors to fractional ones."""
        check_array_arg('cart', cart, (3,))
        cdef np.ndarray[double, ndim=1] frac = np.zeros(3, float)
        self._this.to_frac(&cart[0], &frac[0])
        return frac

    def to_cart(self, np.ndarray[double, ndim=1] frac not None):
        """Convert fractional cell vectors to Cartesian ones."""
        check_array_arg('frac', frac, (3,))
        cdef np.ndarray[double, ndim=1] cart = np.zeros(3, float)
        self._this.to_cart(&frac[0], &cart[0])
        return cart

    def iwrap_mic(self, delta):
        """Wrap Cartesian vectors in-place to fractional range [-0.5, 0.5[."""
        if isinstance(delta, np.ndarray):
            if len(delta.shape) == 1:
                self._iwrap_mic_one(delta)
            elif len(delta.shape) == 2:
                self._iwrap_mic_many(delta)
            else:
                raise TypeError('The argument delta must be a one- or two-dimensional numpy array.')
        else:
            raise TypeError('The argument delta must be a numpy array.')

    cdef _iwrap_mic_one(self, np.ndarray[double, ndim=1] delta):
        check_array_arg('delta', delta, (3,))
        self._this.iwrap_mic(&delta[0])

    cdef _iwrap_mic_many(self, np.ndarray[double, ndim=2] deltas):
        check_array_arg('deltas', deltas, (-1, 3))
        cdef int i
        cdef int n = deltas.shape[0]
        for i in range(n):
            self._this.iwrap_mic(&deltas[i, 0])

    def iwrap_box(self, delta):
        """Wrap Cartesian vectors in-place to fractional range [-0, 1[."""
        if isinstance(delta, np.ndarray):
            if len(delta.shape) == 1:
                self._iwrap_box_one(delta)
            elif len(delta.shape) == 2:
                self._iwrap_box_many(delta)
            else:
                raise TypeError('The argument delta must be a one- or two-dimensional numpy array.')
        else:
            raise TypeError('The argument delta must be a numpy array.')

    cdef _iwrap_box_one(self, np.ndarray[double, ndim=1] delta):
        check_array_arg('delta', delta, (3,))
        self._this.iwrap_box(&delta[0])

    cdef _iwrap_box_many(self, np.ndarray[double, ndim=2] deltas):
        check_array_arg('deltas', deltas, (-1, 3))
        cdef int i
        cdef int n = deltas.shape[0]
        for i in range(n):
            self._this.iwrap_box(&deltas[i, 0])


def ranges_cutoff(Cell cell, np.ndarray[double, ndim=1] center not None, double cutoff):
    """Get the ranges of periodic images which contain a given cutoff sphere.

    This function assumes the space is divided into boxes by crystal planes
    at integer indexes. For example, these planes for the first cell vector
    are parallel to the second and third cell vector and have one point in
    (first vector)*index where index is the integer index for these planes.
    Similar definitions are used for the other two cell vectors. The
    returned ranges are arrays referring to the integer indexes that
    demarcate the cutoff sphere. One could interpret the result as a
    supercell that contains the entire cutoff sphere.

    Parameters
    ----------
    center
        The center of the cutoff sphere.
    cutoff
        The radius of the cutoff sphere.

    Returns
    -------
    ranges_begin
        The lower bounds of the intervals containing the cutoff sphere.
        These integers are the highest indices of the crystal planes below
        the cutoff sphere.
    ranges_end
        The (non-inclusive) upper bounds of the intervals containing the
        cutoff sphere. These integers are the lowest indices of the crystal
        planes above the cutoff sphere.

    """
    check_array_arg('center', center, (3,))
    cdef np.ndarray[int, ndim=1] ranges_begin = np.zeros(3, np.intc)
    cdef np.ndarray[int, ndim=1] ranges_end = np.zeros(3, np.intc)
    cellcutoff.iterators.ranges_cutoff(cell._this, &center[0], cutoff,
                                       &ranges_begin[0], &ranges_end[0])
    return ranges_begin, ranges_end


def create_random_cell(int seed, int nvec, double scale=10.0,
                       double ratio=0.1, bool cuboid=False) -> Cell:
    """Return a sensible random cell.

    Parameters
    ----------
    seed
        A random seed, for reproducile testing.
    nvec
        The dimensionality of the cell.
    scale
        Determines the overall size of the random cell. The components of the
        cell vectors are sampled from a uniform random distribution over the
        interval [-scale, scale].
    ratio
        Controls the minimal volume of the cell. Random cell vectors are tried
        until the volume is above (ratio*scale)**nvec. When a suitable volume is
        found, these cell vectors are used. The closer to 2, the more cubic the
        cell. High values of this parameter may make this function slow because
        many trials will be needed.
    cuboid
        When True, a random orthorhombic cell is generated.

    Returns
    -------
    cell
        The random cell object.

    """
    cdef cell.Cell* cpp_cell = cell.create_random_cell(seed, nvec, scale, ratio, cuboid)
    cdef Cell result = Cell.__new__(Cell, None, initvoid=True)
    result._this = cpp_cell
    return result
