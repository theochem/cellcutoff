#!/usr/bin/env bash

rm -vr gcov
find -type f | grep "gcda$" | xargs rm -v
