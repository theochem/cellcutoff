|Build Status| |codecov|

CellCutoff is a library for periodic boundary conditions and real-space
cutoff calculations.

Dependencies
============

-  C++11 compiler (tested: GNU and Intel)
-  Python, >=2.7.x, <3.x
-  Numpy, >1.9
-  Cython, >= 0.24.1

Runtime environment configuration
=================================

The instructions below explain how to install everything in
``${HOME}/.local``, such that you do not need root permissions to
install all the software. To make this work, some environment variables
must be set, e.g. in your ``~/.bashrc``.

::

    export PATH=${HOME}/.local/bin:${PATH}
    export LD_LIBRARY_PATH=${HOME}/.local/lib:${HOME}/.local/lib64:${LD_LIBRARY_PATH}

Installation of ``cellcutoff``
==============================

Build (for installation in home directory):

::

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=${HOME}/.local
    make

Testing (in the build directory):

::

    make check

Install:

::

    make install

Installation of ``python-cellcutoff``
=====================================

In-place build and test

::

    cd python-cellcutoff
    ./setup.py build_ext -i -L${LD_LIBRARY_PATH}
    nosetests -v

Build and install (into home directory):

::

    ./setup.py build_ext -L${LD_LIBRARY_PATH}
    ./setup.py install --user

.. |Build Status| image:: https://travis-ci.org/theochem/cellcutoff.svg?branch=master
   :target: https://travis-ci.org/theochem/cellcutoff
.. |codecov| image:: https://codecov.io/gh/theochem/cellcutoff/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/theochem/cellcutoff
