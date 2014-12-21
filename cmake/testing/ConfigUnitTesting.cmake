# gather unit tests ==> create UnitTests.cmake
execute_process(
    COMMAND python ${CMAKE_SOURCE_DIR}/cmake/testing/unit-test-gather.py ${PROJECT_SOURCE_DIR}
    )

execute_process(COMMAND mkdir -p ${PROJECT_BINARY_DIR}/unit-tests)
execute_process(COMMAND cp ${PROJECT_SOURCE_DIR}/cmake/testing/unit-test-wrap.py ${PROJECT_BINARY_DIR}/unit-tests)

# set unit test libs
set(UNIT_TEST_LIBS
    dirac
    ${EXTERNAL_LIBS}
    )

macro(add_unit_test _path _name)
    add_custom_target(
        ${_name}_files_copied
        COMMAND cp ${_path}/* ${PROJECT_BINARY_DIR}/unit-tests/${_name}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        )
    add_executable(
        ${_name}.x
        ${_path}/test.F90
        )
    set_target_properties(
        ${_name}.x
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
        "unit-tests/${_name}"
        )
    target_link_libraries(
        ${_name}.x
        ${UNIT_TEST_LIBS}
        )
    set_property(
        TARGET ${_name}.x 
        PROPERTY 
        LINKER_LANGUAGE Fortran
        )
    add_dependencies(
        ${_name}.x
        ${_name}_files_copied
        )
    add_dependencies(${_name}.x UNIT_TEST_LIBS)
    add_test(
        ${_name}
        python ${PROJECT_BINARY_DIR}/unit-tests/unit-test-wrap.py --build-dir=${PROJECT_BINARY_DIR} --name=${_name}
        )
endmacro()

include(UnitTests)
