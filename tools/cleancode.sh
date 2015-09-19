#!/bin/bash
echo "Cleaning code in ${PWD} and subdirectories."
# Remove files that "contaminate" the source tree. The cmake build dir should be included
# here.
rm -vf python-celllists/celllists/celllists.cpp
rm -vf python-celllists/celllists/*.pyc
rm -vfr python-celllists/build
# split output of find at newlines.
IFS=$'\n'
# send all relevant files to the code cleaner
find celllists *.* tools python-celllists | \
    egrep "(\.cpp$)|(\.h$)|(\.in$)|(\.sh$)|(\.py$)|(\.pyx$)|(\.pxd$)|(\.txt$)|(\.conf$)|(.gitignore$)" | \
    xargs ./tools/codecleaner.py
