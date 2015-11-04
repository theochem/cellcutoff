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
'''Unit tests'''


from celllists import *

import numpy as np


def test_subcell():
    c = Cell(np.random.uniform(-10, 10, (3, 3)))
    sc, shape = c.subcell(0.1)
    assert sc.spacings.max() < 0.1


def test_reciprocal():
    c = Cell(np.random.uniform(-10, 10, (3, 3)))
    gc = c.reciprocal()
    assert abs(np.dot(c.vecs, gc.vecs.T) - np.identity(3)).max() < 1e-3
