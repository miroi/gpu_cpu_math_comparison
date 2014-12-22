# set cdash buildname
set(BUILDNAME
    "${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}-${CMAKE_Fortran_COMPILER_ID}-${BLAS_TYPE}-${cmake_build_type_tolower}"
    CACHE STRING
    "Name of build on the dashboard"
    )

if(NOT ENABLE_BENCHMARKS)
    # set ctest own timeout
    set(DART_TESTING_TIMEOUT
        "1500"
        CACHE STRING
        "Set timeout in seconds for every single test"
        )
endif()

macro(dirac_test _name _label)
    add_test(${_name} ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/test/${_name}/test --binary-dir=${PROJECT_BINARY_DIR} --work-dir=${PROJECT_BINARY_DIR}/test/${_name} --verbose)
    if(NOT "${_label}" STREQUAL "")
        set_tests_properties(${_name} PROPERTIES LABELS "${_label}")
    endif()
endmacro()

if(RUN_UNIT_TESTS)
  include(ConfigUnitTesting)
endif()

#include(Tests)
include(CTest)
enable_testing()
