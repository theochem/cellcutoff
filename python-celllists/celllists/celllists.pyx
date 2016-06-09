#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CellList is a 3D domain decomposition library.
# Copyright (C) 2011-2015 The CellList Development Team
#
# This file is part of CellList.
#
# CellList is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# CellList is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
#--
'''Python wrapper for the CellLists library'''


import numpy as np
cimport numpy as np
np.import_array()

from libc.string cimport memcpy

cimport cell


__all__ = ['Cell']


def check_array_arg(name, arg, expected_shape):
    if not arg.flags['C_CONTIGUOUS']:
        raise TypeError('Argument %s must be C_CONTIHUOUS.' % arg)
    for i, n in enumerate(expected_shape):
        if n >= 0 and arg.shape[i] != n:
            raise TypeError(('Axis %i of argument %s has length %i but while '
                             'expecting %i') % (i, name, arg.shape[i], n))


cdef class Cell(object):
    def __cinit__(self, np.ndarray[double, ndim=2] vecs=None, initvoid=False):
        if initvoid:
            self._this = NULL
        elif vecs is None:
            self._this = new cell.Cell()
        else:
            check_array_arg('vecs', vecs, (-1, 3))
            assert vecs.shape[0] <= 3
            nvec = vecs.shape[0]
            self._this = new cell.Cell(&vecs[0,0], nvec)

    def __init__(self, np.ndarray[double, ndim=2] vecs=None):
        pass

    def __dealloc__(self):
        if self._this != NULL:
            del self._this

    def subcell(self, double threshold):
        cdef np.ndarray[int, ndim=1] shape = np.zeros(3, np.intc)
        cdef cell.Cell* cpp_cell = self._this.create_subcell(threshold, &shape[0])
        cdef Cell cell = Cell.__new__(Cell, None, initvoid=True)
        cell._this = cpp_cell
        return cell, shape

    def reciprocal(self):
        cdef cell.Cell* cpp_cell = self._this.create_reciprocal()
        cdef Cell cell = Cell.__new__(Cell, None, initvoid=True)
        cell._this = cpp_cell
        return cell

    property nvec:
        def __get__(self):
            return self._this.nvec()

    property vecs:
        def __get__(self):
            cdef np.ndarray[double, ndim=2] vecs = np.zeros((3, 3), float)
            memcpy(&vecs[0, 0], self._this.vecs(), sizeof(double)*9);
            return vecs

    property gvecs:
        def __get__(self):
            cdef np.ndarray[double, ndim=2] gvecs = np.zeros((3, 3), float)
            memcpy(&gvecs[0, 0], self._this.gvecs(), sizeof(double)*9);
            return gvecs

    property volume:
        def __get__(self):
            return self._this.volume()

    property gvolume:
        def __get__(self):
            return self._this.gvolume()

    property lengths:
        def __get__(self):
            cdef np.ndarray[double, ndim=1] lengths = np.zeros(3, float)
            memcpy(&lengths[0], self._this.lengths(), sizeof(double)*3);
            return lengths

    property glengths:
        def __get__(self):
            cdef np.ndarray[double, ndim=1] glengths = np.zeros(3, float)
            memcpy(&glengths[0], self._this.glengths(), sizeof(double)*3);
            return glengths

    property spacings:
        def __get__(self):
            cdef np.ndarray[double, ndim=1] spacings = np.zeros(3, float)
            memcpy(&spacings[0], self._this.spacings(), sizeof(double)*3);
            return spacings

    property gspacings:
        def __get__(self):
            cdef np.ndarray[double, ndim=1] gspacings = np.zeros(3, float)
            memcpy(&gspacings[0], self._this.gspacings(), sizeof(double)*3);
            return gspacings

    property cubic:
        def __get__(self):
            return self._this.cubic()

    property cuboid:
        def __get__(self):
            return self._this.cuboid()

    def to_frac(self, np.ndarray[double, ndim=1] cart not None):
        check_array_arg('cart', cart, (3,))
        cdef np.ndarray[double, ndim=1] frac = np.zeros(3, float)
        self._this.to_frac(&cart[0], &frac[0])
        return frac

    def to_cart(self, np.ndarray[double, ndim=1] frac not None):
        check_array_arg('frac', frac, (3,))
        cdef np.ndarray[double, ndim=1] cart = np.zeros(3, float)
        self._this.to_cart(&frac[0], &cart[0])
        return cart

    def iwrap_mic(self, delta):
        if isinstance(delta, np.ndarray):
            if len(delta.shape) == 1:
                self._iwrap_mic_one(delta)
            elif len(delta.shape) == 2:
                self._iwrap_mic_many(delta)
            else:
                raise TypeError('The argument delta must be a one- or two-dimensional numpy array.')
        else:
            raise TypeError('The argument delta must be a numpy array.')

    def _iwrap_mic_one(self, np.ndarray[double, ndim=1] delta not None):
        check_array_arg('delta', delta, (3,))
        self._this.iwrap_mic(&delta[0])

    def _iwrap_mic_many(self, np.ndarray[double, ndim=2] deltas not None):
        check_array_arg('deltas', deltas, (-1, 3))
        cdef int i
        cdef int n = deltas.shape[0]
        for i in range(n):
            self._this.iwrap_mic(&deltas[i, 0])

    def iwrap_box(self, delta):
        if isinstance(delta, np.ndarray):
            if len(delta.shape) == 1:
                self._iwrap_box_one(delta)
            elif len(delta.shape) == 2:
                self._iwrap_box_many(delta)
            else:
                raise TypeError('The argument delta must be a one- or two-dimensional numpy array.')
        else:
            raise TypeError('The argument delta must be a numpy array.')

    def _iwrap_box_one(self, np.ndarray[double, ndim=1] delta not None):
        check_array_arg('delta', delta, (3,))
        self._this.iwrap_box(&delta[0])

    def _iwrap_box_many(self, np.ndarray[double, ndim=2] deltas not None):
        check_array_arg('deltas', deltas, (-1, 3))
        cdef int i
        cdef int n = deltas.shape[0]
        for i in range(n):
            self._this.iwrap_box(&deltas[i, 0])

    def ranges_cutoff(self, np.ndarray[double, ndim=1] center not None, double cutoff):
        check_array_arg('center', center, (3,))
        cdef np.ndarray[int, ndim=1] ranges_begin = np.zeros(3, np.intc)
        cdef np.ndarray[int, ndim=1] ranges_end = np.zeros(3, np.intc)
        self._this.ranges_cutoff(&center[0], cutoff, &ranges_begin[0], &ranges_end[0])
        return ranges_begin, ranges_end
