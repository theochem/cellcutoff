.. image:: https://travis-ci.com/theochem/cellcutoff.svg?branch=master
    :target: https://travis-ci.com/theochem/cellcutoff/builds
.. image:: https://img.shields.io/codecov/c/github/theochem/cellcutoff/master.svg
    :target: https://codecov.io/gh/theochem/cellcutoff
.. image:: https://img.shields.io/conda/v/theochem/cellcutoff.svg
    :target: https://anaconda.org/theochem/cellcutoff
.. image:: https://img.shields.io/conda/vn/theochem/python-cellcutoff.svg
    :target: https://anaconda.org/theochem/python-cellcutoff
.. image:: https://img.shields.io/conda/pn/theochem/cellcutoff.svg
    :target: https://anaconda.org/theochem/cellcutoff
.. image:: https://img.shields.io/github/release/theochem/cellcutoff.svg
    :target: https://github.com/theochem/cellcutoff/releases

CellCutoff is a library for periodic boundary conditions and real-space
cutoff calculations.


Installation
============

When you are interested in using cellcutoff (without needing to modify it), you
can install cellcutoff with conda. After installing and activating a miniconda
environment, run:

.. code-block:: bash

  conda install -c theochem cellcutoff python-cellcutoff

If you are interesed in working on the development of cellcutoff, you first need
to check out the latest version from the git repository

.. code-block:: bash

  git clone git@github.com:theochem/cellcutoff.git
  cd cellcutoff

Then install Roberto and run it in the root of the repository:

.. code-block:: bash

  pip install --user --upgrade 'roberto<2.0.0'
  rob quality

This will build cellcutoff in-place and run all tests. More details for
potential contributors are given in `CONTRIBUTING.rst <CONTRIBUTING.rst>`_.
