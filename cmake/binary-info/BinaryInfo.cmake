include(ConfigGitRevision)

execute_process(
    COMMAND whoami
    TIMEOUT 1
    OUTPUT_VARIABLE USER_NAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

execute_process(
    COMMAND hostname
    TIMEOUT 1
    OUTPUT_VARIABLE HOST_NAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/binary-info/binary_info.py.in
    ${CMAKE_BINARY_DIR}/binary_info.py
    )
    
add_custom_target(
    generate_binary_info
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_BINARY_DIR}/binary_info.py > ${CMAKE_BINARY_DIR}/binary_info.F90
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_BINARY_DIR}/binary_info.py
    )

set_source_files_properties(${CMAKE_BINARY_DIR}/binary_info.F90 PROPERTIES GENERATED 1)
