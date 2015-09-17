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

cimport cell


__all__ = ['Cell']


cdef class Cell:
    def __cinit__(self, np.ndarray[double, ndim=2] vecs=None):
        cdef int nvec = 0
        self._this = NULL
        if vecs is None:
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
