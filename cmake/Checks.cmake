include(CMakeParseArguments)

#There is no eval() in CMake, and macro parameters behave like
#function parameters, not like C preprocessor parameters, so
#there is no way to pass an expression as 'test' parameter here
#to be evaluated. 'test' must expand to the final value that
#you want to assert for inside the if().
#TODO: implement eval() macro as described here:
#(their example is actually 'exec', there is no way
#to emulate an eval() because you need to use math() for
#numbers, string() for strings and if() for logicals.
#http://www.cmake.org/Wiki/CMake/Language_Syntax
#and then eval("if(${test}) set(ret TRUE) else() set(ret FALSE) endif()")
function(at_assert test)
    set(option_args NO_PRINT_TEST)
    set(one_value_args COMMENT_FAIL COMMENT_PASS ERROR_LEVEL)
    cmake_parse_arguments(arg "${option_args}" "${one_value_args}" "" ${ARGN})
    set(test ${test})
    if (test)
        if (NOT arg_COMMENT_PASS STREQUAL "")
            message(STATUS "${arg_COMMENT_PASS}")
        endif ()
    else ()
        if(NOT DEFINED arg_ERROR_LEVEL)
            set(arg_ERROR_LEVEL "FATAL_ERROR")
        endif()
        if(DEFINED arg_COMMENT_FAIL)
            set(msg "${arg_COMMENT_FAIL} ")
        endif()
        if (NOT arg_NO_PRINT_TEST)
            set(msg "${msg}\"${test}\"")
        endif()
        message(${arg_ERROR_LEVEL} "${msg}")
    endif ()
endfunction(at_assert)

#at_assert("1 STREQUAL 2" COMMENT_FAIL "Failed" COMMENT_PASS "Passed")
#at_assert(0 COMMENT_PASS "Passed" NO_PRINT_TEST COMMENT_FAIL "Failed")
#at_assert("1 EQUAL 2" COMMENT_FAIL "Failed" COMMENT_PASS "Passed")

function(at_assert_defined)
    foreach(var IN LISTS ARGN)
        if(NOT DEFINED "${var}")
            message(FATAL_ERROR "Variable ${var} must be defined")
        endif()
    endforeach()
endfunction()

function(at_assert_defined_prefixed prefix)
    foreach(var_base IN LISTS ARGN)
        set(var "${prefix}_${var_base}")
        if(NOT DEFINED "${var}")
            message(FATAL_ERROR "Variable ${var_base} must be defined")
        endif()
    endforeach()
endfunction()

