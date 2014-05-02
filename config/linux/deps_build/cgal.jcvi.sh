#!/bin/bash
export CC=`which gcc`
export CXX=`which g++`
cmake -DCMAKE_BUILD_TYPE=Release \
    -D CMAKE_TOOLCHAIN_FILE=$PRODDL_SRC/config/linux/toolchain.cmake \
    -D CMAKE_INSTALL_PREFIX=$PRODDL_DEPS \
    $PRODDL_DEPS_SRC/CGAL-4.4

