cmake_minimum_required(VERSION 2.8)
FIND_PACKAGE(HDF5 NO_MODULE)
ADD_EXECUTABLE(test_cross test_cross.cpp)
