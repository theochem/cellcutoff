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
"""Package build and install script."""


import numpy as np
from setuptools import setup, Extension
import Cython.Build


# Load the version from version.py, without importing it. This assumes that the last
# line in the file contains a variable defining the version string with single quotes.
with open('cellcutoff/version.py', 'r') as f:
    version = f.read().split('=')[-1].replace('\'', '').strip()


setup(
    name='python-cellcutoff',
    version=version,
    description='CellCutoff is a ibrary for periodic boundary conditions '
                'and real-space cutoff calculations.',
    author='Toon Verstraelen',
    url='https://github.com/theochem/cellcutoff',
    author_email='Toon.Verstraelen@UGent.be',
    cmdclass={'build_ext': Cython.Build.build_ext},
    packages=['cellcutoff'],
    package_data={
        'cellcutoff': ['cellcutoff.pxd', 'cell.pxd'],
    },
    zip_safe=False,
    ext_modules=[Extension(
        "cellcutoff.ext",
        sources=['cellcutoff/ext.pyx'],
        depends=['cellcutoff/cellcutoff.pxd', 'cellcutoff/cell.pxd'],
        libraries=['cellcutoff'],
        include_dirs=[np.get_include()],
        extra_compile_args=['-std=c++11', '-Wall', '-pedantic'],
        language="c++",
    )],
    classifiers=[
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Intended Audience :: Science/Research',
    ],)
