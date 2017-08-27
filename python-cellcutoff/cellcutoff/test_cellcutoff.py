#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
#
# --
"""Unit tests."""


from cellcutoff import Cell

import numpy as np


def test_subcell():
    subcell, shape = Cell(np.random.uniform(-10, 10, (3, 3))).subcell(0.1)
    assert subcell.spacings.max() < 0.1
    assert min(shape) >= 1


def test_reciprocal():
    cell = Cell(np.random.uniform(-10, 10, (3, 3)))
    gcell = cell.reciprocal()
    assert abs(np.dot(cell.vecs, gcell.vecs.T) - np.identity(3)).max() < 1e-3
