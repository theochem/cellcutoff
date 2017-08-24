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


import numpy as np
from setuptools import setup, Extension
import Cython.Build


setup(
    name='python-cellcutoff',
    version='0.0.0',
    description='CellCutoff is a ibrary for periodic boundary conditions and real-space cutoff calculations.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    cmdclass = {'build_ext': Cython.Build.build_ext},
    packages = ['cellcutoff'],
    package_data = {
        'cellcutoff': ['cellcutoff.pxd', 'cell.pxd'],
    },
    zip_safe=False,
    ext_modules=[
        Extension(
            "cellcutoff.cellcutoff",
            sources=['cellcutoff/cellcutoff.pyx'],
            depends=['cellcutoff/cellcutoff.pxd', 'cellcutoff/cell.pxd'],
            libraries=['cellcutoff'],
            include_dirs=[np.get_include()],
            extra_compile_args=['-std=c++11', '-Wall', '-pedantic'],
            language="c++"),
    ],
)
