cmake -DCMAKE_TOOLCHAIN_FILE=../mingw-cross-test/mingw-w64-toolchain \
-DCMAKE_PREFIX_PATH=/mnt/hgfs/proddl/deps/HDF5_1.8.10/cmake/hdf5 \
-DPYTHON_INCLUDE_DIRS=/mnt/hgfs/proddl/deps-cross/include/python2.7 \
-DSWIG_DIR=/mnt/hgfs/proddl/deps_build/swigwin-2.0.9 \
../imp
