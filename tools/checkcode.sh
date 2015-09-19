#!/usr/bin/env bash
./tools/cpplint.py --linelength=90 celllists/*.cpp celllists/*.h
./tools/cpplint.py --linelength=120 celllists/tests/*.cpp celllists/tests/*.h
