include(Checks)
include(CMakeParseArguments)

## Write CMake script_body to a file and return command line
## for running the script and the file name of the script
## These values can be used as arguments of custom_command():
## add_custom_command(COMMAND script_command DEPENDS script_file ...).
## If your script body is a multiline string that contains:
## "
## set(ENV{VAR} "${val}")
## execute_process(COMMAND whatever ...)
## "
## and 'val' is defined in your CMakeLists.txt,
## then this is the way to generate custom_command target that 
## is built under a modified shell environment.
## If you need to set variables like PATH that under Windows use
## the same semicolon separator as CMake uses as a list separator,
## then you have to be extra careful.
## You can provide extra arguments between SCRIPT_BODY and 
## named options. Those will be passed as additional arguments to
## CMake.
function(make_script_command SCRIPT_BODY) 
    set(oneValueArgs SCRIPT_NAME VAR_SCRIPT_COMMAND VAR_SCRIPT_FILE)
    cmake_parse_arguments(arg "" "${oneValueArgs}" "" ${ARGN})
    at_assert_defined_prefixed(arg ${oneValueArgs})
    set(scr_name ${CMAKE_CURRENT_BINARY_DIR}/${arg_SCRIPT_NAME}.cmake)
    file(WRITE ${scr_name} "${SCRIPT_BODY}")
    set(cmd ${CMAKE_COMMAND} ${arg_UNPARSED_ARGUMENTS} -P ${scr_name})
    set (${arg_VAR_SCRIPT_COMMAND} ${cmd} PARENT_SCOPE)
    set (${arg_VAR_SCRIPT_FILE} ${scr_name} PARENT_SCOPE)
endfunction()

