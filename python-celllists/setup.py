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

import numpy as np
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


def parse_cpath():
    import os
    return [s for s in os.getenv('CPATH').split(':') if len(s) > 0]


setup(
    name='python-celllists',
    version='0.0.0',
    description='CellList is a 3D domain decomposition library.',
    author='Toon Verstraelen',
    author_email='Toon.Verstraelen@UGent.be',
    cmdclass = {'build_ext': build_ext},
    packages = ['celllists'],
    package_data = {
        'celllists': ['celllists.pxd', 'cell.pxd'],
    },
    ext_modules=[
        Extension("celllists.celllists",
            sources=['celllists/celllists.pyx'],
            depends=['celllists/celllists.pxd', 'celllists/cell.pxd'],
            libraries=['celllists'],
            include_dirs=[np.get_include()] + parse_cpath(),
            extra_compile_args=['-std=c++11', '-Wall', '-pedantic'],
            language="c++")],
)