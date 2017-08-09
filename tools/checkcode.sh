#!/usr/bin/env bash
./tools/cpplint.py --linelength=90 cellcutoff/*.cpp cellcutoff/*.h
./tools/cpplint.py --linelength=120 cellcutoff/tests/*.cpp cellcutoff/tests/*.h
