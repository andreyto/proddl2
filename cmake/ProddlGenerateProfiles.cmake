
set(config_dir "${CMAKE_INSTALL_PREFIX}/config")
set(bin_dir "${CMAKE_INSTALL_PREFIX}/bin")

set(DEPS_RC_BASE "proddl-deps.rc")
set(PRODDL_DEPS_RC "${config_dir}/${DEPS_RC_BASE}")
set(PRODDL_WRAPPER "${bin_dir}/proddl-wrapper")
set(WRAPPER "${PRODDL_WRAPPER}")
set(CONF_BASE proddl.json)
set(PRODDL_CONFIG "${config_dir}/${CONF_BASE}")
set(PRODDL_ROOT "${CMAKE_INSTALL_PREFIX}")
set(PRODDL_HOME "${PRODDL_INSTALL_PREFIX}")
set(DATA_DIR_BASE "data")
set(TEST_DATA_DIR_BASE "test_data")

if (WIN32)
	foreach(var_name PRODDL_WRAPPER WRAPPER PRODDL_DEPS_RC)
		set(${var_name} "${${var_name}}.bat")
	endforeach()
endif()

foreach(dir_name config bin)
    set(dir_out "${CMAKE_BINARY_DIR}/${dir_name}")
    file(MAKE_DIRECTORY "${dir_out}")
    file(GLOB arch_files_in ${PROJECT_SOURCE_DIR}/${dir_name}/${PRODDL_TARGET_ENV}/*)
    file(GLOB noarch_files_in ${PROJECT_SOURCE_DIR}/${dir_name}/noarch/*)
    foreach( file_in IN LISTS arch_files_in noarch_files_in)
        if(NOT file_in MATCHES ".*\\.cmake\$" AND NOT file_in MATCHES ".*build.*")  
            get_filename_component(file_out "${file_in}" NAME)
            #remove '.in' from file_out name if present
            #TODO: use COPY option to configure_file if there is no '.in'
            string(REGEX REPLACE ".in\$" "" file_out "${file_out}")
            set(file_out "${dir_out}/${file_out}")
            configure_file(${file_in} ${file_out} @ONLY)
        endif()
    endforeach()
endforeach()

install(DIRECTORY "${CMAKE_BINARY_DIR}/bin" 
    DESTINATION "${CMAKE_INSTALL_PREFIX}" 
    FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                GROUP_EXECUTE GROUP_READ)

install(DIRECTORY "${CMAKE_BINARY_DIR}/config" 
    DESTINATION "${CMAKE_INSTALL_PREFIX}")
