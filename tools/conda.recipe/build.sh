#!/usr/bin/env bash

rm -rf build
mkdir build
cd build
cmake \
  -DCMAKE_INSTALL_PREFIX=${PREFIX} \
  -DCMAKE_BUILD_TYPE=release \
  ..
VERBOSE=1 make install
