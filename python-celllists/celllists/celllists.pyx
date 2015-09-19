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


cdef class Cell:
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

    property vecs:
        def __get__(self):
            cdef np.ndarray[double, ndim=2] vecs = np.zeros((3, 3), float)
            memcpy(&vecs[0, 0], self._this.vecs(), sizeof(double)*9);
            return vecs

    property nvec:
        def __get__(self):
            return self._this.nvec()
