message("-- Who compiled             : ${USER_NAME}")
message("-- Host                     : ${HOST_NAME}")
message("-- System name              : ${CMAKE_SYSTEM_NAME}")
message("-- System                   : ${CMAKE_SYSTEM}")
message("-- CMake version            : ${CMAKE_VERSION}")
message("-- CMake generator          : ${CMAKE_GENERATOR}")
message("-- Python version           : ${PYTHON_VERSION_STRING}")
message("-- Processor                : ${CMAKE_SYSTEM_PROCESSOR}")
message("-- 64-bit integers          : ${ENABLE_64BIT_INTEGERS}")
message("-- MPI                      : ${ENABLE_MPI}")

message("-- Fortran compiler         : ${CMAKE_Fortran_COMPILER}")
message("-- Fortran compiler version : ${CMAKE_Fortran_COMPILER_ID} ${PYTHON_Fortran_VERSION}")
message("-- Fortran compiler flags   : ${CMAKE_Fortran_FLAGS}  ${CMAKE_Fortran_FLAGS_${cmake_build_type_toupper}}")

message("-- C compiler               : ${CMAKE_C_COMPILER}")
message("-- C compiler version       : ${CMAKE_C_COMPILER_ID} ${PYTHON_C_VERSION}")
message("-- C compiler flags         : ${CMAKE_C_FLAGS}  ${CMAKE_C_FLAGS_${cmake_build_type_toupper}}")

message("-- C++ compiler             : ${CMAKE_CXX_COMPILER}")
message("-- C++ compiler version     : ${CMAKE_CXX_COMPILER_ID} ${PYTHON_CXX_VERSION}")
message("-- C++ compiler flags       : ${CMAKE_CXX_FLAGS}  ${CMAKE_CXX_FLAGS_${cmake_build_type_toupper}}")

message("-- Definitions              : ${_list_of_definitions}")

message("-- BLAS                     : ${BLAS_LIBRARIES}${BLAS_LIBRARIES_INFO}")
message("-- LAPACK                   : ${LAPACK_LIBRARIES}${LAPACK_LIBRARIES_INFO}")
message("-- Libraries                : ${EXTERNAL_LIBS}")
message("-- Explicit libs            : ${EXPLICIT_LIBS}${EXPLICIT_LIBS_INFO}")
message("-- Static linking           : ${ENABLE_STATIC_LINKING}")
message("-- Git branch               : ${GIT_BRANCH}")
message("-- Git revision             : ${GIT_REVISION}")
message("-- Last commit author       : ${GIT_LAST_COMMIT_AUTHOR}")
message("-- Last commit date         : ${GIT_LAST_COMMIT_DATE}")

message("-- Configuration time       : ${CONFIGURATION_TIME} UTC")

# get size of static allocations
foreach(_binary ${STATIC_MEM_INFO_BINARIES})
    add_custom_target(
        static_mem_${_binary}
            COMMAND ${CMAKE_SOURCE_DIR}/cmake/binary-info/get_static_size.py ${_binary}.x
        )
endforeach()
