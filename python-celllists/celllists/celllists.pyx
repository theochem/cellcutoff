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


cdef class Cell(object):
    def __cinit__(self, np.ndarray[double, ndim=2] vecs=None, initvoid=False):
        if initvoid:
            self._this = NULL
        elif vecs is None:
            self._this = new cell.Cell()
        else:
            assert vecs.flags['C_CONTIGUOUS']
            assert vecs.shape[0] <= 3
            assert vecs.shape[1] == 3
            nvec = vecs.shape[0]
            self._this = new cell.Cell(&vecs[0,0], nvec)

    def __init__(self, np.ndarray[double, ndim=2] vecs=None):
        pass

    def __dealloc__(self):
        if self._this != NULL:
            del self._this

    def create_subcell(self, double threshold):
        cdef np.ndarray[int, ndim=1] shape = np.zeros(3, int)
        cdef cell.Cell* cpp_cell = self._this.create_subcell(threshold, &shape[0])
        cell = Cell(initvoid=True)
        cell._this = cpp_cell
        return cell, shape

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
