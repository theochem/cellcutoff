#!/usr/bin/env bash

NAME=cellcutoff
VERSION=$(git describe --tags)
git archive --prefix=${NAME}-${VERSION}/ HEAD | bzip2 > ${NAME}-${VERSION}.tar.bz2
