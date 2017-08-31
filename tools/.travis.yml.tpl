env:
  matrix:
  - MYCONDAPY=2.7
  #- MYCONDAPY=3.5
  #- MYCONDAPY=3.6
  global:
    - secure: "${TPL_ANACONDA_TOKEN}"
    - secure: "${TPL_GITHUB_TOKEN}"

# Do not use Travis Python to save some time.
language: generic
os:
- linux
- osx
osx_image: xcode6.4
dist: trusty
sudo: false

branches:
  only:
    - master
    - /^[0-9]+\.[0-9]+(\.[0-9]+)?([ab][0-9]+)?$/


install:
# Get miniconda. Take the right version, so re-installing python is only needed for 3.5.
- if [[ "$MYCONDAPY" == "2.7" ]]; then
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
      wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh;
    else
      wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda.sh;
    fi;
  else
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
      wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
    else
      wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
    fi;
  fi
- bash miniconda.sh -b -p $HOME/miniconda
- source $HOME/miniconda/bin/activate
- hash -r

# Configure conda and get a few essentials
- conda config --set always_yes yes
- conda config --add channels theochem
- conda update -q conda
# Get the right python version for building. This only does something for 3.5.
# Install extra package needed to make things work. Most things can be listed as
# dependencies on metal.yaml and setup.py, unless setup.py already imports them.
# Install conda tools for packaging and uploading
- conda install python=${MYCONDAPY} numpy cython cppcheck doxygen nose conda-build anaconda-client gtest cmake
# Install more recent stuff with pip
- pip install pylint codecov coverage pycodestyle pydocstyle flake8
# Get the latest cpplint
- wget https://raw.githubusercontent.com/google/styleguide/gh-pages/cpplint/cpplint.py
- chmod +x cpplint.py
# Show conda info for debugging
- conda info -a

# Install the latest cardboardlinter
- if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
    pip install --upgrade git+https://github.com/theochem/cardboardlint.git@master#egg=cardboardlint;
  fi

# Set the version info from the git tag
- git fetch origin --tags
- export PROJECT_VERSION=$(python tools/gitversion.py)
- python tools/gitversion.py python > ${PYDIR}/${PROJECT_NAME}/version.py
- python tools/gitversion.py cmake > CMakeListsVersion.txt.in

script:
# Static linting
# ~~~~~~~~~~~~~~
- if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
    cardboardlinter --refspec $TRAVIS_BRANCH -f static &&
    (cd ${PYDIR}; cardboardlinter --refspec $TRAVIS_BRANCH -f static);
  fi

# Conda-based tests
# ~~~~~~~~~~~~~~~~~

# This is needed to simulate the end-user situation, and it also creates packages for a
# conda release.

# Build the conda packages
- conda build -q tools/conda.recipe
- (cd ${PYDIR}; ./setup.py; conda build -q tools/conda.recipe)
# Install Conda packages
- conda install --use-local ${CONDA_PKG_NAME_CPP} ${CONDA_PKG_NAME_PY}
# Run the unittests, using the installed package
- ${CONDA_PREFIX}/libexec/${PROJECT_NAME}/test_${PROJECT_NAME}
- (cd; nosetests ${PROJECT_NAME} -v --detailed-errors)

# In-place building and testing
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This is allows for coverage analysis and dynamic linting. The compiler settings used
# here are not suitable for releases, so we need to recompile and rerun the tests.

# Install GCC compilers for in-place builds, even on OSX because clang does not manage to
# compiler our C++ code.
- conda install gcc libgcc;

# PY project
# Uninstall conda package, to be sure. The conda cpp package is still used.
- conda uninstall ${CONDA_PKG_NAME_PY}
- (cd ${PYDIR}; python setup.py build_ext -i --define CYTHON_TRACE_NOGIL)
# Run nosetests without coverage.xml output. That file is broken by nosetests (pyx files
# not include) and gets priority over .coverage, which contains everything.
- (cd ${PYDIR}; nosetests ${PROJECT_NAME} -v --detailed-errors --with-coverage --cover-package=${PROJECT_NAME} --cover-tests --cover-erase --cover-inclusive --cover-branches; coverage xml -i)

- if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
    (cd ${PYDIR}; cardboardlinter --refspec $TRAVIS_BRANCH -f 'dynamic');
  fi

# CPP project
# Uninstall conda package, to be sure.
- conda uninstall ${CONDA_PKG_NAME_CPP}

# Install dependencies
- ${TPL_DEPENDENCIES}
- (mkdir debug; cd debug; cmake cmake -DCMAKE_BUILD_TYPE=debug ..; make)
# DYLD_LIBRARY_PATH is needed on OSX to let it find the shared libraries in
# ${CONDA_PREFIX}/lib. It should not have any effect on Linux.
- (cd debug; DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${CONDA_PREFIX}/lib/ ./${PROJECT_NAME}/tests/test_${PROJECT_NAME})

- if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
    cardboardlinter --refspec $TRAVIS_BRANCH -f 'dynamic';
  fi

# Some other stuff
# ----------------

# Make CPP source package for github deployment
- (cd debug; make sdist)
# Build PY source package for github deployment.
- (cd ${PYDIR}; python setup.py sdist)

after_success:
# Upload the coverage analysis
- codecov --file ${PYDIR}/coverage.xml
- codecov

# In deployment, the env var TRAVIS_TAG contains the name of the current tag, if any.
deploy:
- provider: releases
  skip_cleanup: true
  api_key: ${GITHUB_TOKEN}
  file:
  - debug/${CONDA_PKG_NAME_CPP}-${TRAVIS_TAG}.tar.gz
  - ${PYDIR}/dist/${CONDA_PKG_NAME_PY}-${TRAVIS_TAG}.tar.gz
  on:
    repo: ${GITHUB_REPO_NAME}
    tags: true
    condition: "$TRAVIS_TAG != *[ab]* && $MYCONDAPY == 2.7 && $TRAVIS_OS_NAME == linux"
  prerelease: false

- provider: releases
  skip_cleanup: true
  api_key: ${GITHUB_TOKEN}
  file:
  - debug/${CONDA_PKG_NAME_CPP}-${TRAVIS_TAG}.tar.gz
  - ${PYDIR}/dist/${CONDA_PKG_NAME_PY}-${TRAVIS_TAG}.tar.gz
  on:
    repo: ${GITHUB_REPO_NAME}
    tags: true
    condition: "$TRAVIS_TAG == *[ab]* && $MYCONDAPY == 2.7 && $TRAVIS_OS_NAME == linux"
  prerelease: true

- provider: script
  skip_cleanup: true
  script: anaconda -t $ANACONDA_TOKEN upload --force -l alpha ${HOME}/miniconda/conda-bld/*/${CONDA_PKG_NAME_CPP}-*.tar.bz2 ${HOME}/miniconda/conda-bld/*/${CONDA_PKG_NAME_PY}-*.tar.bz2
  on:
    repo: ${GITHUB_REPO_NAME}
    tags: true
    condition: "$TRAVIS_TAG == *a*"

- provider: script
  skip_cleanup: true
  script: anaconda -t $ANACONDA_TOKEN upload --force -l beta ${HOME}/miniconda/conda-bld/*/${CONDA_PKG_NAME_CPP}-*.tar.bz2 ${HOME}/miniconda/conda-bld/*/${CONDA_PKG_NAME_PY}-*.tar.bz2
  on:
    repo: ${GITHUB_REPO_NAME}
    tags: true
    condition: "$TRAVIS_TAG == *b*"

- provider: script
  skip_cleanup: true
  script: anaconda -t $ANACONDA_TOKEN upload --force -l main ${HOME}/miniconda/conda-bld/*/${CONDA_PKG_NAME_CPP}-*.tar.bz2 ${HOME}/miniconda/conda-bld/*/${CONDA_PKG_NAME_PY}-*.tar.bz2
  on:
    repo: ${GITHUB_REPO_NAME}
    tags: true
    condition: "$TRAVIS_TAG != *[ab]*"
