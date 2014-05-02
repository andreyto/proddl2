include(CMakeParseArguments)
include(Checks)

## Get delimiter for entries within PATH, PYTHONPATH etc
function(get_path_separator python var_path_sep)
    execute_process(COMMAND "${python}" -c "import os; print (os.pathsep),"
        RESULT_VARIABLE cmd_error_status
        OUTPUT_VARIABLE path_sep
	OUTPUT_STRIP_TRAILING_WHITESPACE
	)
    if( cmd_error_status )
        message(FATAL_ERROR execute_process failed)
    endif()
    set(${var_path_sep} "${path_sep}" PARENT_SCOPE)
endfunction()

function(string_join SEP VALUES OUTPUT)
  string (REGEX REPLACE "([^\\]|^);" "\\1${SEP}" _TMP_STR "${VALUES}")
  string (REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") #fixes escaping
  set (${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
endfunction()

## Function to build native string for environment variables like PATH
## PATHS must be a CMake list of paths in CMake format. It can be obtained
## like this: file(TO_CMAKE_PATH "$ENV{PYTHONPATH}" PYTHONPATH)
## PATHS_SEP is platform separator of elements in PATH and alike (";"
## on Windows, ":" on Linux - can be obtained with python os.pathsep)
## OUT name of output variable
## The reason for this function is that although TO_CMAKE_PATH knows PATHS_SEP
## and splits into list correctly, the reciprocal function TO_NATIVE_PATH does
## not replace the CMake ";"
function(to_native_path_list PATHS PATHS_SEP OUT)
    file(TO_NATIVE_PATH "${PATHS}" _tmp_var)
    string(REPLACE ";" "${PATHS_SEP}" _tmp_var "${_tmp_var}")
    set(${OUT} "${_tmp_var}" PARENT_SCOPE)
endfunction()


## Add a CTEST test compiled against GTest
## add_test_gtest(name SOURCES source [,source,...] [LIBS lib,...] [GTEST_NAME gtest_main.c] [TEST_ARGS ...])
## Where GTEST_NAME and TEST_ARGS will be taken from variables of the same name
## if not passed as named arguments
function(add_test_gtest name)
    set(multi_value_args SOURCES LIBS TEST_ARGS)
	set(one_value_args GTEST_MAIN)
    cmake_parse_arguments(arg "" "${one_value_args}" "${multi_value_args}" ${ARGN})
    if (NOT DEFINED arg_GTEST_MAIN)
        at_assert_defined(GTEST_MAIN)
        set(arg_GTEST_MAIN ${GTEST_MAIN})
    endif ()
	if(NOT DEFINED arg_TEST_ARGS)
        at_assert_defined(TEST_ARGS)
        set(arg_TEST_ARGS ${TEST_ARGS})
	endif ()
	at_assert_defined_prefixed(arg SOURCES)
	add_executable(${name} ${arg_SOURCES} ${arg_GTEST_MAIN})
	target_link_libraries(${name} ${arg_LIBS} gtest)
	add_test(NAME ${name} COMMAND ${name} ${arg_TEST_ARGS})
endfunction()

macro(arg_from_var var require)
    if (NOT DEFINED arg_${var})
		if(require)
			at_assert_defined(${var})
		endif ()
        set(arg_${var} ${${var}})
    endif ()
endmacro()
