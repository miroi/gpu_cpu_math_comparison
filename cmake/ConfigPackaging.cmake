foreach(
    EXECUTABLE
    ${LIST_OF_EXECUTABLES}
    )
    install(
        TARGETS ${EXECUTABLE}
        DESTINATION share/dirac
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ             GROUP_EXECUTE
        WORLD_READ             WORLD_EXECUTE
        )
endforeach()

install(
    FILES ${PROJECT_BINARY_DIR}/pam
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )

# workaround to install pam-dirac symlink:
# 1) copy pam to pam-dirac
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${PROJECT_BINARY_DIR}/pam ${PROJECT_BINARY_DIR}/pam-dirac)")
# 2) install real file pam-dirac
install(
    FILES ${PROJECT_BINARY_DIR}/pam-dirac
    DESTINATION bin
    PERMISSIONS
    OWNER_READ OWNER_WRITE OWNER_EXECUTE
    GROUP_READ             GROUP_EXECUTE
    WORLD_READ             WORLD_EXECUTE
    )
# 3) remove real file
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_INSTALL_PREFIX}/bin/pam-dirac)")
# 4) create symlink
install(CODE "EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_INSTALL_PREFIX}/share/dirac/pam ${CMAKE_INSTALL_PREFIX}/bin/pam-dirac)")

# write git hash to build dir
file(WRITE ${PROJECT_BINARY_DIR}/GIT_HASH "${GIT_REVISION}")
# copy version info to install dir
install(
    FILES ${PROJECT_BINARY_DIR}/GIT_HASH ${PROJECT_SOURCE_DIR}/VERSION
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
    )

install(
    FILES ${PROJECT_SOURCE_DIR}/utils/wrapper.py
    DESTINATION share/dirac
    PERMISSIONS
    OWNER_READ OWNER_WRITE
    GROUP_READ
    WORLD_READ
    )

install(
    DIRECTORY
    ${PROJECT_SOURCE_DIR}/basis
    ${PROJECT_SOURCE_DIR}/basis_dalton
    ${PROJECT_SOURCE_DIR}/basis_ecp
    DESTINATION share/dirac
    PATTERN .git EXCLUDE
    )

include(InstallRequiredSystemLibraries)

set(CPACK_SOURCE_GENERATOR "TGZ")

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "DIRAC")
set(CPACK_PACKAGE_VENDOR "R. Bast, H. J. Aa. Jensen, T. Saue, and L. Visscher, and contributors")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "R. Bast")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/doc/installation/general.rst")
set(CPACK_RESOURCE_FILE_LICENSE    "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_FILE_NAME        "DIRAC")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "DIRAC-${DIRAC_VERSION}-Source")
set(CPACK_PACKAGE_INSTALL_DIRECTORY
    "DIRAC ${DIRAC_VERSION}"
    )

SET(EXPORT_DIR ${CMAKE_BINARY_DIR}/export)
SET(CPACK_SOURCE_INSTALLED_DIRECTORIES
    "${EXPORT_DIR};/"
    )

include(CPack)
include(IncludeFiles)

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_custom_target(release
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/fetch_external_sources ${CMAKE_SOURCE_DIR}
        COMMAND mkdir -p ${EXPORT_DIR}
        COMMAND tar cf ${EXPORT_DIR}/sources.tar -C ${CMAKE_SOURCE_DIR} ${DIRAC_SOURCES} ${INCLUDE_FILES}
        COMMAND tar xf ${EXPORT_DIR}/sources.tar -C ${EXPORT_DIR}
        COMMAND rm     ${EXPORT_DIR}/sources.tar
        COMMAND rm -rf ${EXPORT_DIR}/build
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/remove_unreleased_code  ${EXPORT_DIR} ${PROJECT_BINARY_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/basis                        ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/basis_dalton                 ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/basis_ecp                    ${EXPORT_DIR}
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/remove_basis_sets       ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/test                         ${EXPORT_DIR}
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/remove_unreleased_tests ${EXPORT_DIR} ${PROJECT_BINARY_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/CMakeLists.txt               ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/COPYING                      ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/CHANGELOG.rst                ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/CTestConfig.cmake            ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/atomic-start-x2c-pam         ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/cmake                        ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/external                     ${EXPORT_DIR}
        COMMAND rm -f  ${EXPORT_DIR}/external/pcmsolver/.git  # prevens "git init; git add ." inside tarball
        COMMAND rm -f  ${EXPORT_DIR}/external/unit-tests/.git # prevens "git init; git add ." inside tarball
        COMMAND cp -r  ${PROJECT_BINARY_DIR}/pamadm-generated           ${EXPORT_DIR}/cmake
        COMMAND cp -rL ${CMAKE_SOURCE_DIR}/doc                          ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/pam.in                       ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/setup                        ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/utils                        ${EXPORT_DIR}
        COMMAND cp -r  ${CMAKE_SOURCE_DIR}/src/xcfun                    ${EXPORT_DIR}/src
        COMMAND cp     ${CMAKE_SOURCE_DIR}/src/main/main.F90            ${EXPORT_DIR}/src/main
        COMMAND cp     ${CMAKE_SOURCE_DIR}/VERSION                      ${EXPORT_DIR}
        COMMAND cp     ${CMAKE_SOURCE_DIR}/.gitignore                   ${EXPORT_DIR}
        COMMAND ${CMAKE_SOURCE_DIR}/maintenance/remove_copyright_tags   ${EXPORT_DIR}
        COMMAND echo "${GIT_REVISION}"                                > ${EXPORT_DIR}/cmake/GIT_HASH
        COMMAND make package_source
        COMMENT "Packaging source files"
        VERBATIM
        )
endif()
