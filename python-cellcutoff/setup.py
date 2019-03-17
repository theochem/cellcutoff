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

import os
import sys

from Cython.Distutils.build_ext import new_build_ext as build_ext
import numpy as np
from setuptools import setup, Extension


NAME = 'cellcutoff'


def get_version():
    """Read __version__ from version.py, with exec to avoid importing it."""
    try:
        with open(os.path.join(NAME, 'version.py'), 'r') as f:
            myglobals = {}
            exec(f.read(), myglobals)  # pylint: disable=exec-used
        return myglobals['__version__']
    except IOError:
        return "0.0.0.post0"


class DebugBuildExt(build_ext):
    """Fix rpath bug in distutils for OSX."""

    def finalize_options(self):
        super().finalize_options()
        # Special treatment of rpath in case of OSX, to work around python
        # distutils bug. This constructs proper rpath arguments for clang.
        if sys.platform[:6] == "darwin":
            for path in self.rpath:
                for ext in self.extensions:
                    ext.extra_link_args.append("-Wl,-rpath," + path)
            self.rpath[:] = []


setup(
    name=NAME,
    version=get_version(),
    package_dir={NAME: NAME},
    packages=[NAME],
    description='CellCutoff is a library for periodic boundary conditions '
                'and real-space cutoff calculations.',
    author='The CellCutoff development team',
    url='https://github.com/theochem/cellcutoff',
    cmdclass={'build_ext': DebugBuildExt},
    package_data={
        'cellcutoff': ['ext.pxd', 'cell.pxd'],
    },
    ext_modules=[Extension(
        "cellcutoff.ext",
        sources=['cellcutoff/ext.pyx'],
        depends=['cellcutoff/ext.pxd', 'cellcutoff/cell.pxd'],
        libraries=['cellcutoff'],
        include_dirs=[np.get_include()],
        extra_compile_args=['-std=c++11', '-Wall', '-pedantic'],
        language="c++",
    )],
    classifiers=[
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Intended Audience :: Science/Research',
    ])
