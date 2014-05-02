#!/bin/bash
export CC=`which gcc`
export CXX=`which g++`
export EXTRA_CFLAGS=$CFLAGS
export CFLAGS=
export EXTRA_CXXFLAGS=$CXXFLAGS
export CXXFLAGS=
cmake \
    -D DEVEL_MODE=1 \
    -D PRODDL_TARGET_ENV=linux \
    -D CMAKE_TOOLCHAIN_FILE=$PRODDL_SRC/config/linux/toolchain.cmake \
    -D CMAKE_INSTALL_PREFIX=$PRODDL_HOME \
    -D CMAKE_BUILD_TYPE=Release \
    $PRODDL_SRC

#cmake --build .
#make install

