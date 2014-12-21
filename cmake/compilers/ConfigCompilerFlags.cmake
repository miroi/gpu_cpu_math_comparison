# This variables are not guaranteed to be defined for all compilers or languages
if(NOT DEFINED CMAKE_C_COMPILER_ID)
    message(FATAL_ERROR "CMAKE_C_COMPILER_ID variable is not defined! (CMake Error)")
endif()

if(NOT DEFINED CMAKE_CXX_COMPILER_ID)
    message(FATAL_ERROR "CMAKE_CXX_COMPILER_ID variable is not defined! (CMake Error)")
endif()

if(NOT DEFINED CMAKE_Fortran_COMPILER_ID)
    message(FATAL_ERROR "CMAKE_Fortran_COMPILER_ID variable is not defined! (CMake Error)")
endif()

# macro for saving flags to cache
include(SaveCompilerFlags)

if(CMAKE_C_COMPILER_WORKS AND CMAKE_CXX_COMPILER_WORKS AND CMAKE_Fortran_COMPILER_WORKS)

    # check whether we can use -xHost on this machine
    include(ConfigXHostFlag)

    # set flags
    include(CFlags)
    include(CXXFlags)
    include(FortranFlags)

    # check for correct compiler version with Python script and CMake variable (if defined)
    include(FindCompilerVersion)
endif()
