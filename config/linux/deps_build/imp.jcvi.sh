#!/bin/bash
#export PATH=/usr/local/packages/gsl/bin:$PATH
export CC=`which gcc`
export CXX=`which g++`
cmake -DCMAKE_BUILD_TYPE=Release \
    -D CMAKE_TOOLCHAIN_FILE=$PRODDL_SRC/config/linux/toolchain.cmake \
    -D CMAKE_INSTALL_PREFIX=$IMP_HOME \
    $PRODDL_DEPS_SRC/imp

