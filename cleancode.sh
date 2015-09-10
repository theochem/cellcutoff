#!/bin/bash
echo "Cleaning code in ${PWD} and subdirectories."
# split output of find at newlines.
IFS=$'\n'
# send all relevant files to the code cleaner
find celllists *.* tools | \
    egrep "(\.cpp$)|(\.h$)|(\.in$)|(\.sh$)|(\.py$)|(\.txt$)|(\.conf$)|(.gitignore$)" | \
    xargs ./tools/codecleaner.py
