cmake_minimum_required(VERSION 2.8.5)

PROJECT(PRODDL)

SET( ${PROJECT_NAME}_MAJOR_VERSION 2 )
SET( ${PROJECT_NAME}_MINOR_VERSION 0 )
SET( ${PROJECT_NAME}_PATCH_LEVEL 0 )

SET( ${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_MAJOR_VERSION}.${${PROJECT_NAME}_MINOR_VERSION}.${${PROJECT_NAME}_PATCH_LEVEL}" )

set(PACKAGE_VERSION ${${PROJECT_NAME}_VERSION})

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#MESSAGE("CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}")
message("CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}")

if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-g -std=c++11 ${CMAKE_CXX_FLAGS}")
    message("CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}")
    message("CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}")
endif()

set(DEVEL_MODE 0 CACHE INTEGER "Set to 1 for development mode")

set_property(DIRECTORY
	PROPERTY COMPILE_DEFINITIONS_DEBUG ATDEBUG ATDEBUG_LEVEL=6 ATLOG_LEVEL=6
	)
set_property(DIRECTORY
	PROPERTY COMPILE_DEFINITIONS_RELEASE ATDEBUG_LEVEL=1 ATLOG_LEVEL=1 NODEBUG
	)

if(DEVEL_MODE)
	#ATALWAYS() assertion failure and AT_THROW() break into a debugger
	set_property(DIRECTORY
		APPEND
		PROPERTY COMPILE_DEFINITIONS_DEBUG ATDEBUG_BREAK
		)
	#GTest variables below are aimed at invoking the debugger on assertions and letting
	#the exception to go unhandled by GTest so that the debugger can examine the call stack
	#maybe also GTEST_THROW_ON_FAILURE=1 ?
	set_property(TEST
		APPEND
		PROPERTY ENVIRONMENT GTEST_BREAK_ON_FAILURE=1 GTEST_CATCH_EXCEPTIONS=0 
		)
endif()

set(T_NUM "double" CACHE STRING "C datatype for real numbers")
add_definitions(-DPRODDL_T_NUM=${T_NUM})

find_package( HDF5 REQUIRED )
find_package( Blitz++ REQUIRED )
#something is wrong with Blitz lib on Windows - symbols reported
#as duplicate with C run-time lib. Lib is only needed if BZ_DEBUG is defined
set(Blitz_LIBRARIES "")
#if( "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" )
#	set(Blitz_LIBRARIES ${Blitz_Debug_LIBRARY})
#endif()
MESSAGE("Blitz_LIBRARIES=${Blitz_LIBRARIES}")
set(FFTW_USE_STATIC_LIBS TRUE)
find_package( FFTW REQUIRED )
# Sets PYTHONINTERP_FOUND and PYTHON_EXECUTABLE.
find_package(PythonInterp REQUIRED)
set(PRODDL_PYTHON_VERSION_MAJOR_MINOR "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")

# Search for Boost ###########################
#set(BOOST_ROOT "" CACHE FILEPATH "Optional root directory for Boost library search")
if (NOT Boost_FOUND)
    MESSAGE("Looking for Boost with find_package()...")
    #run "cmake --help-module FindBoost" for a full list of options
    find_package( Boost REQUIRED COMPONENTS filesystem iostreams program_options)
endif (NOT Boost_FOUND)

if (Boost_FOUND)
    MESSAGE("Boost information:")
    MESSAGE("BOOST_ROOT=${BOOST_ROOT}")
    MESSAGE("   Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
    MESSAGE("   Boost_LIBRARIES: ${Boost_LIBRARIES}")
    MESSAGE("   Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}")

    link_directories ( ${Boost_LIBRARY_DIRS} )
    include_directories ( ${Boost_INCLUDE_DIRS} )
else (Boost_FOUND)
    message(FATAL_ERROR "Boost libraries are not found")    
endif (Boost_FOUND)
###############################################

#find_library(PACK_LIB MSVCR100)
#find_file(PACK_LIB MSVCR100.DLL)
#find_file(PACK_LIB hdf5dll.dll)
#find_file(PACK_LIB boost_unit_test_framework-vc100-mt-1_44.dll)
#message(STATUS "PACK_LIB=${PACK_LIB}")

#disable autolinking of missing boost_zlib from boost iostreams
add_definitions(-DBOOST_IOSTREAMS_NO_LIB)
#In Boost 1.44 we need to define it to use version 3 of filesystem API
add_definitions(-DBOOST_FILESYSTEM_VERSION=3)
##TODO file_string() was replaced with string() at some point -
##figure out when and make it conditional
add_definitions(-DPRODDL_BOOST_FILE_STRING=string)

# For make-based builds, defines make target named test.
# For Visual Studio builds, defines Visual Studio project named RUN_TESTS.
enable_testing()

set(TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/test_data
    CACHE FILEPATH "Root directory for input data files used by the test cases")

# Override internal GTest option to get rid of multiple symbol linker errors (see link below)
# https://code.google.com/p/googletest/wiki/FAQ#I_am_building_my_project_with_Google_Test_in_Visual_Studio_and_a
SET(gtest_force_shared_crt ON CACHE BOOL "Use shared (DLL) run-time lib even when Google Test is built as static lib.")

#We want to be able to test from the build diretory,
#without having to run install every time. We unify
#the build output locations below.
# On Windows, things end up under Debug/ or Release/.
if (WIN32)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
else (WIN32)
  # If not on Windows, put them under more standard UNIX-like places.
  # Otherwise they will be under each sub-projects directory
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)
endif (WIN32)

set(PY_INSTALL_SUBDIR py)
set(PY_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/${PY_INSTALL_SUBDIR}")
set(PY_BUILD_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/${PY_INSTALL_SUBDIR}")
set(PY_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/python")
set(PY_PACKAGE_NAME "${PROJECT_NAME}")
set(EXE_PREFIX "proddl-")
set(PY_SOURCE_PACKAGE_SUBDIR "lib")
set(PY_SOURCE_PACKAGE_DIR "${PY_SOURCE_DIR}/${PY_SOURCE_PACKAGE_SUBDIR}")

include(ProddlGenerateProfiles)

add_subdirectory(vendor/bob)

add_subdirectory(vendor/gtest)

add_subdirectory(src/PRODDL)

add_subdirectory(${PY_SOURCE_DIR})

install(DIRECTORY ${PROJECT_SOURCE_DIR}/${DATA_DIR_BASE}/ DESTINATION ${DATA_DIR_BASE} USE_SOURCE_PERMISSIONS)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/${TEST_DATA_DIR_BASE}/ DESTINATION ${TEST_DATA_DIR_BASE} USE_SOURCE_PERMISSIONS)

include(ProddlInstaller)
