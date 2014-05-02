include(CustomCommands)
include(Util)
include(CMakeParseArguments)
include(Checks)

## CMD_BUILD_INSTALL_PRE and CMD_BUILD_INSTALL_POST are multiline strings that
## will be pasted into the end beginning and end of the build-install target
##CMAKE script.
## You should also define other variables needed by your specific setup.py.in
## substitutions. If necessary to make these variables local, wrap the call to
## this function in another function's body
function(install_python_package)
    set(multi_value_args DEPS)
	set(one_value_args 
		PYTHON_EXECUTABLE 
		PY_BUILD_INSTALL_PREFIX 
		PY_INSTALL_PREFIX
		PACKAGE_VERSION 
		CMD_BUILD_INSTALL_PRE
		CMD_BUILD_INSTALL_POST
		ALWAYS_COPY)
    cmake_parse_arguments(arg "" "${one_value_args}" "${multi_value_args}" ${ARGN})
    arg_from_var(PYTHON_EXECUTABLE TRUE)
	arg_from_var(PY_BUILD_INSTALL_PREFIX TRUE)
	arg_from_var(PY_INSTALL_PREFIX TRUE)
	arg_from_var(PACKAGE_VERSION TRUE)
	at_assert_defined_prefixed(arg DEPS PYTHON_EXECUTABLE PY_BUILD_INSTALL_PREFIX PY_INSTALL_PREFIX PACKAGE_VERSION)

	if(NOT DEFINED arg_ALWAYS_COPY)
		if (WIN32)
			set(arg_ALWAYS_COPY FALSE)
		else()
			set(arg_ALWAYS_COPY TRUE)
		endif()
	endif()

	set(PY_BUILD_INSTALL_PREFIX ${arg_PY_BUILD_INSTALL_PREFIX})
	set(PY_INSTALL_PREFIX ${arg_PY_INSTALL_PREFIX})
	set(PACKAGE_VERSION ${arg_PACKAGE_VERSION})
	set(PY_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
	set(PY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/py_build")
	set(SETUP_PY_IN "${PY_SOURCE_DIR}/setup.py.in")
	#easy_install needs the source dir as the last argument,
	#and expects that dir to contain setup.py,
	#plus it will download distribute tarball into it.
	#We can either copy all python source into cmake binary dir
	#and build in that copy, or build in the original source dir.
	#We choose the late so that the development cycle modify/test 
	#is easier.
	set(SETUP_PY    "${PY_SOURCE_DIR}/setup.py")
	set(DEPS        ${arg_DEPS})
	set(OUTPUT      "${PY_BUILD_DIR}/timestamp")

	configure_file(${SETUP_PY_IN} ${SETUP_PY} @ONLY)

	file(MAKE_DIRECTORY ${arg_PY_BUILD_INSTALL_PREFIX})
	file(MAKE_DIRECTORY ${PY_BUILD_DIR})

	get_path_separator(${arg_PYTHON_EXECUTABLE} PATH_SEP)

	#Doing just 'setup.py build' or 'setup.py develop'
	#does not generate
	#entry point scripts, thus making in-build testing
	#impossible. Thus, we do 'setup.py install' twice -
	#during build and during install, into different locations.

	#We use --site-dirs argument for easy_install to force installation
	#into a directory which is not in the PYTHONPATH 
	#(i.e. .pth files would not work where). easy_install still complains,
	#in a slightly differnt way, so we also modify the setup.py above.

	if(arg_ALWAYS_COPY)
		set(ei_always_copy "--always-copy")
	else()
		set(ei_always_copy "")
	endif()

	set(cmd_build_install "${arg_PYTHON_EXECUTABLE} ${SETUP_PY} easy_install 
						--install-dir ${arg_PY_BUILD_INSTALL_PREFIX}
						--site-dirs ${arg_PY_BUILD_INSTALL_PREFIX}
						--build-directory ${PY_BUILD_DIR} 
						${ei_always_copy}
						.")

	#lines in a string body should be deindented or else indent goes to the output and 
	#breaks Python import of the file
	set(sitecustomize_body "##This file causes .pth files to be processed if the directory
##of this file is on PYTHONPATH. Placing this file in such a directory
##is the way to create an additional site-packages like directory in which
##eggs can be installed without having to edit any system or user wide Python
##config files.
import site
import os,inspect
site.addsitedir(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
")

	## Create a script multiline string that modifies the environment and create
	## a build-time target out of that script. This is not currently needed here,
	## but we keep the infrastructure in case we will have to build a C++ extension,
	## in which case it might become necessary

	message("PYTHON_EXECUTABLE=${PYTHON_EXECUTABLE}") 

	file(TO_CMAKE_PATH "$ENV{PYTHONPATH}" PYTHONPATH)
	list(APPEND PYTHONPATH "${arg_PY_BUILD_INSTALL_PREFIX}")

	set(cmd "

	set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})
	include(Util)
	#set(ENV{CFLAGS} -O3 -DNDEBUG -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0)
	#set(ENV{CXXFLAGS} -O3 -DNDEBUG -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0)
	${arg_CMD_BUILD_INSTALL_PRE}
	#Make sure that Python installer does not complain about install path not being in
	#PYTHONPATH
	to_native_path_list(\"${PYTHONPATH}\" \"${PATH_SEP}\" ENV{PYTHONPATH})
	execute_process(
		COMMAND ${cmd_build_install}
		RESULT_VARIABLE cmd_error_status
		)
	if( cmd_error_status )
		message(FATAL_ERROR \"execute_process failed\")
	endif()
	file(WRITE ${arg_PY_BUILD_INSTALL_PREFIX}/sitecustomize.py \"${sitecustomize_body}\")
	${arg_CMD_BUILD_INSTALL_POST}
	")

	make_script_command("${cmd}"
		SCRIPT_NAME py_build 
		VAR_SCRIPT_COMMAND cmd_build_install_scr 
		VAR_SCRIPT_FILE build_install_scr)

	add_custom_command(OUTPUT ${OUTPUT}
						COMMAND ${cmd_build_install_scr}
						#COMMAND ${cmd_develop}
						#COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
						DEPENDS ${build_install_scr} ${DEPS}
						WORKING_DIRECTORY ${PY_SOURCE_DIR} 
						VERBATIM
						COMMENT "Command to build and install Python modules into build dir")

	add_custom_target(PY_BUILD ALL 
		DEPENDS ${OUTPUT} 
		COMMENT "Target to build and install Python modules into build dir")

	#trailing slash copies files but not the directory itself
	install(DIRECTORY ${arg_PY_BUILD_INSTALL_PREFIX}/ DESTINATION ${arg_PY_INSTALL_PREFIX} USE_SOURCE_PERMISSIONS)
endfunction()
