## Same as configure_file() but will raise error
## on unmatched substitution patterns (not finished yet)
function(configure_file_assert file_in file_out)
    set(rx "@[a-zA-Z0-9_]+@")
    configure_file(${file_in} ${file_out} @ONLY ${ARGN})
    file(STRINGS ${file_in} lines_var REGEX "${rx}")
    message(FATAL_ERROR "Not finished")
endfunction()
