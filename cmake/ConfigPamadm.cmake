if(EXISTS ${CMAKE_SOURCE_DIR}/pamadm)
    # check whether all pamadm-generated files exist
    if(EXISTS "${PROJECT_BINARY_DIR}/pamadm-generated/Sources.cmake"      AND
       EXISTS "${PROJECT_BINARY_DIR}/pamadm-generated/Definitions.cmake"  AND
       EXISTS "${PROJECT_BINARY_DIR}/pamadm-generated/IncludeFiles.cmake" AND
       EXISTS "${PROJECT_BINARY_DIR}/pamadm-generated/Tests.cmake"        AND
       EXISTS "${PROJECT_BINARY_DIR}/pamadm-generated/pamadm-flags")
       # they exist, read previously used flags from "pamadm-flags"
       file(READ ${PROJECT_BINARY_DIR}/pamadm-generated/pamadm-flags PAMADM_FLAGS)
       # check for Windows
       if(PAMADM_FLAGS STREQUAL "")
           message(FATAL_ERROR "No pamadm flags in ${PROJECT_BINARY_DIR}/pamadm-generated/pamadm-flags")
       endif()
    else()
        # set flags to '-a'
        set(PAMADM_FLAGS "-a")
    endif()

    if("${PAMADM_FLAGS}" MATCHES "builddir")
        execute_process(
            COMMAND python pamadm ${PAMADM_FLAGS}
            TIMEOUT 1
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_QUIET
        )
    else()
        execute_process(
            COMMAND python pamadm ${PAMADM_FLAGS} --builddir=${PROJECT_BINARY_DIR}
            TIMEOUT 1
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
            OUTPUT_QUIET
        )
    endif()

    # tweak null device for Windows
    if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        set(NULL_DEVICE "NUL")
    else()
        set(NULL_DEVICE "/dev/null")
    endif()

    add_custom_target(
        run_pamadm
        COMMAND python pamadm ${PAMADM_FLAGS} > ${NULL_DEVICE}
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        )

    message("-- Packaging sources and tests using pamadm ${PAMADM_FLAGS} (to change flags, run pamadm manually with new flags)")
else()
    # pamadm does not exist
    # this must then be released code
    # in this case don't bother
endif()
