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
"""Unit tests."""


import numpy as np
from numpy.testing import assert_equal, assert_allclose

from pytest import raises

from cellcutoff import Cell, create_random_cell


def test_subcell():
    for nvec in 1, 2, 3:
        cell = create_random_cell(1, nvec)
        subcell, shape = cell.subcell(0.1)
        assert subcell.spacings[:nvec].max() < 0.1
        assert min(shape[:nvec]) >= 1
        assert_equal(shape[nvec:], [0] * (3 - nvec))


def test_reciprocal():
    for nvec in 1, 2, 3:
        cell = create_random_cell(2, nvec)
        gcell = cell.reciprocal()
        assert_allclose(np.dot(cell.vecs, gcell.vecs.T), np.identity(3), atol=1e-8)


def test_errors():
    vecs = np.zeros((2, 3)).T
    with raises(TypeError):
        Cell(vecs)
    vecs = np.zeros(5)
    with raises(ValueError):
        Cell(vecs)
    vecs = np.zeros((5, 3))
    with raises(TypeError):
        Cell(vecs)
    vecs = np.zeros((3, 3))
    with raises(ValueError):
        Cell(vecs)
    vecs = np.zeros((3, 5))
    with raises(TypeError):
        Cell(vecs)


def test_nvec():
    vecs = np.identity(3)
    assert Cell(vecs).nvec == 3
    assert Cell(vecs[:2]).nvec == 2
    assert Cell(vecs[:1]).nvec == 1
    assert Cell(vecs[:0]).nvec == 0
    assert Cell(None).nvec == 0


def test_properties():
    for nvec in 1, 2, 3:
        cell = create_random_cell(3, nvec)
        assert_allclose(np.dot(cell.vecs, cell.gvecs.T), np.identity(3), atol=1e-8)
        assert_allclose(cell.volume, abs(np.linalg.det(cell.vecs)))
        assert_allclose(cell.gvolume, abs(np.linalg.det(cell.gvecs)))
        assert_allclose(cell.volume * cell.gvolume, 1.0)
        assert_allclose((cell.vecs**2).sum(axis=1)**0.5, cell.lengths)
        assert_allclose((cell.gvecs**2).sum(axis=1)**0.5, cell.glengths)
        assert_allclose(cell.spacings, 1.0/cell.glengths)
        assert_allclose(cell.gspacings, 1.0/cell.lengths)
        assert not cell.cuboid
        assert not cell.cubic


def test_cart_frac():
    for nvec in 1, 2, 3:
        cell = create_random_cell(4, nvec)
        cart1 = np.random.uniform(-20.0, 20.0, 3)
        frac1 = cell.to_frac(cart1)
        cart2 = cell.to_cart(frac1)
        frac2 = cell.to_frac(cart2)
        assert_allclose(cart1, cart2)
        assert_allclose(frac1, frac2)
        assert_allclose(cart1, np.dot(frac1, cell.vecs))


def test_iwrap_mic():
    for nvec in 1, 2, 3:
        for _ in range(50):
            cell = create_random_cell(5, nvec, 1.0)
            cart1 = np.random.uniform(-20.0, 20.0, 3)
            cart2 = cart1.copy()
            cell.iwrap_mic(cart2)
            frac1 = cell.to_frac(cart1)
            frac2 = cell.to_frac(cart2)
            assert (abs(frac2[:nvec]) < 0.5 + 1e-8).all()
            delta = (frac1 - frac2)[:nvec]
            assert_allclose(delta, np.round(delta), atol=1e-8)


def test_iwrap_mic_many():
    for nvec in 1, 2, 3:
        cell = create_random_cell(6, nvec, 1.0)
        cart1 = np.random.uniform(-20.0, 20.0, (10, 3))
        cart2 = cart1.copy()
        cell.iwrap_mic(cart2)
        frac1 = np.array([cell.to_frac(c1) for c1 in cart1])
        frac2 = np.array([cell.to_frac(c2) for c2 in cart2])
        assert (abs(frac2[:, :nvec]) < 0.5 + 1e-8).all()
        delta = (frac1 - frac2)[:, :nvec]
        assert_allclose(delta, np.round(delta), atol=1e-8)


def test_iwrap_box():
    for nvec in 1, 2, 3:
        for _ in range(50):
            cell = create_random_cell(7, nvec, 1.0)
            cart1 = np.random.uniform(-20.0, 20.0, 3)
            cart2 = cart1.copy()
            cell.iwrap_box(cart2)
            frac1 = cell.to_frac(cart1)
            frac2 = cell.to_frac(cart2)
            assert (frac2[:nvec] < 1.0 + 1e-8).all()
            assert (frac2[:nvec] >= 0.0 - 1e-8).all()
            delta = (frac1 - frac2)[:nvec]
            assert_allclose(delta, np.round(delta), atol=1e-8)


def test_iwrap_box_many():
    for nvec in 1, 2, 3:
        cell = create_random_cell(8, nvec, 1.0)
        cart1 = np.random.uniform(-20.0, 20.0, (10, 3))
        cart2 = cart1.copy()
        cell.iwrap_box(cart2)
        frac1 = np.array([cell.to_frac(c1) for c1 in cart1])
        frac2 = np.array([cell.to_frac(c2) for c2 in cart2])
        assert (frac2[:, :nvec] < 1.0 + 1e-8).all()
        assert (frac2[:, :nvec] >= 0.0 - 1e-8).all()
        delta = (frac1 - frac2)[:, :nvec]
        assert_allclose(delta, np.round(delta), atol=1e-8)


def test_ranges_cutoff_simple():
    cell = Cell(np.identity(3) * 5.0)
    ranges_begin, ranges_end = cell.ranges_cutoff(np.array([0.0, 0.0, 0.0]), 4.0)
    assert_equal(ranges_begin, [-1, -1, -1])
    assert_equal(ranges_end, [1, 1, 1])
    ranges_begin, ranges_end = cell.ranges_cutoff(np.array([2.5, 2.5, 2.5]), 5.0)
    assert_equal(ranges_begin, [-1, -1, -1])
    assert_equal(ranges_end, [2, 2, 2])
