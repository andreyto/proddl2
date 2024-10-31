project(BOB_IO CXX C)

include_directories(${PROJECT_SOURCE_DIR}/include)

set(bob_io_sources 
src/File.cc
src/HDF5ArrayFile.cc
src/HDF5Attribute.cc
src/HDF5Dataset.cc
src/HDF5File.cc
src/HDF5Group.cc
src/HDF5Types.cc
src/HDF5Utils.cc
src/array.cc
src/array_type.cc
src/blitz_array.cc
src/cast.cc
src/logging.cc
)

MESSAGE("PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
MESSAGE("PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}")

find_package( HDF5 COMPONENTS C REQUIRED )
#somehow the variable below might not be
#defined even when HDF5_LIBRARIES contains the 
#corresponding entries
if(NOT DEFINED HDF5_C_LIBRARIES)
    set(HDF5_C_LIBRARIES ${HDF5_LIBRARIES})
endif()
find_package( Blitz++ REQUIRED )

find_package( Boost REQUIRED COMPONENTS filesystem system iostreams program_options unit_test_framework)

set(libs ${Boost_LIBRARIES} ${HDF5_LIBRARIES} ${Blitz_LIBRARIES})

add_library(bob_io ${bob_io_sources})
target_link_libraries(bob_io ${libs})
target_include_directories(bob_io PRIVATE ${HDF5_INCLUDE_DIRS})

add_executable(test_hdf5 src/test/hdf5.cc)
target_link_libraries(test_hdf5 bob_io ${libs})
target_include_directories(test_hdf5 PRIVATE ${HDF5_INCLUDE_DIRS})

MESSAGE("HDF5_C_LIBRARIES=${HDF5_C_LIBRARIES}")
MESSAGE("HDF5_FOUND=${HDF5_FOUND}")
MESSAGE("HDF5_LIBRARIES=${HDF5_LIBRARIES}")

add_test(NAME test_hdf5 COMMAND test_hdf5)

set(BOB_IO_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include ${HDF5_INCLUDE_DIRS} CACHE INTERNAL "Bob IO include dirs")
