include(Util)
include_directories(${PROJECT_SOURCE_DIR}/include ${BOB_IO_INCLUDE_DIRS} ${FFTW_INCLUDE_DIRS} ${Blitz_INCLUDE_DIR})

MESSAGE("PROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}")
MESSAGE("PROJECT_BINARY_DIR=${PROJECT_BINARY_DIR}")

#if(MSVC)
#set_property(DIRECTORY
#	PROPERTY LINK_FLAGS_DEBUG NODEFAULTLIB=LIBCMT NODEFAULTLIB=LIBCMTD
#	)
#endif()

#use GTest
include_directories(${PROJECT_SOURCE_DIR}/vendor/gtest/include)
#link_directories(${PROJECT_BINARY_DIR}/vendor/gtest)

# Alternative way to find Google Test's files if we have not built it
# from higher level CMakeLists.txt
#include_directories(${GTEST_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
# Where Google Test's libraries can be found.
#link_directories(${gtest_BINARY_DIR}/src)

set(GTEST_MAIN Testing/gtestc_main.cpp)
set(TEST_ARGS --test-data-dir ${TEST_DATA_DIR})

### Libraries

set(common_sources 
	Common/options_io_json.cpp 
	Common/logger.cpp 
	Common/string_util.cpp 
	External/jsoncpp/jsoncpp.cpp
	External/Dbg/dbg.cpp
	Geom/gdiam_simple.cpp
	External/Gdiam/gdiam.cpp
    External/Pdb++/pdb_all.cc
	)
add_library(proddl ${common_sources})

### Tests

add_test_gtest(test_gtestc_samp SOURCES Testing/test_gtestc_samp.cpp)

add_test_gtest(test_bob_hdf5 SOURCES IO/test_bob_hdf5.cpp LIBS ${Boost_LIBRARIES} bob_io)

add_test_gtest(test_bz_tiny SOURCES Common/test_bz_tiny.cpp LIBS ${Blitz_LIBRARIES})

add_test_gtest(test_options SOURCES Common/test_options.cpp LIBS proddl ${Boost_LIBRARIES})

add_test_gtest(test_dock_io_bin SOURCES IO/test_dock_io_bin.cpp LIBS proddl ${Boost_LIBRARIES})

add_test_gtest(test_common SOURCES 
	Common/test_logger.cpp 
	Common/test_queue.cpp
	Common/test_bz_cast.cpp
	Common/test_c_array.cpp
	Common/test_nd_index_iter.cpp
	Common/test_bz_ext.cpp
	LIBS proddl ${Boost_LIBRARIES} ${Blitz_LIBRARIES})

add_test_gtest(test_geom SOURCES 
	External/Gdiam/gdiam_test.cpp
	Geom/test_transformation.cpp 
#	Geom/test_rotational_grid.cpp
	Geom/test_pairdist.cpp
	Geom/test_gdiam_simple.cpp
	Geom/test_cluster.cpp
	Grid/test_grid.cpp
	LIBS proddl ${Boost_LIBRARIES} ${Blitz_LIBRARIES})

add_test_gtest(test_math_os SOURCES 
	Math/test_pfactors.cpp
	LIBS proddl ${Boost_LIBRARIES} ${Blitz_LIBRARIES})


### Programs

add_executable(${EXE_PREFIX}dock-fft docking_main.cpp)
target_link_libraries(${EXE_PREFIX}dock-fft proddl ${Boost_LIBRARIES} ${FFTW_LIBRARIES} bob_io)

add_executable(${EXE_PREFIX}export export_models.cpp)
target_link_libraries(${EXE_PREFIX}export proddl ${Boost_LIBRARIES})

### Install

install(TARGETS ${EXE_PREFIX}dock-fft ${EXE_PREFIX}export RUNTIME DESTINATION bin)
