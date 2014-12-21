if(NOT DEFINED DEFAULT_Fortran_FLAGS_SET OR RESET_FLAGS)

# custom flags are defined by setup --custom-fc-flags
if(DEFINED CUSTOM_Fortran_FLAGS)
    # set custom compiler flags (for every build type)
    set(CMAKE_Fortran_FLAGS "${CUSTOM_Fortran_FLAGS}")
    # special flags for build types will be empty
    set(CMAKE_Fortran_FLAGS_DEBUG "")
    set(CMAKE_Fortran_FLAGS_RELEASE "")
    set(CMAKE_Fortran_FLAGS_PROFILE "")
else()
    # custom flags are not defined
    if(CMAKE_Fortran_COMPILER_ID MATCHES GNU) # this is gfortran
        set(CMAKE_Fortran_FLAGS         "-g -fcray-pointer -fbacktrace -DVAR_GFORTRAN -DVAR_MFDS -fno-range-check")
                                        # -fcray-pointer is for VAR_MPI2
        set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops -w")
        set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
        
        if(ENABLE_STATIC_LINKING)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} -static"
                )
        endif()
        
        if(ENABLE_64BIT_INTEGERS)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
                )
        endif()
        
        if(ENABLE_BOUNDS_CHECK)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} -fbounds-check"
                )
        endif()
        
        if(ENABLE_CODE_COVERAGE)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} -fprofile-arcs -ftest-coverage"
                )
        endif()
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
        # Intel has different flags on Windows and Linux + MacOS
        if(CMAKE_SYSTEM_NAME MATCHES "Windows")
            set(CMAKE_Fortran_FLAGS         "/w /assume:byterecl /DVAR_IFORT /Z7 /traceback")
            set(CMAKE_Fortran_FLAGS_DEBUG   "/O0")
            set(CMAKE_Fortran_FLAGS_RELEASE "/O3 /Qip")
            
            # intel compiler xHost flag
            if(XHOST_FLAG_AVAILABLE)
                set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} /QxHost")
            endif()
            
            # windows: pg - none
            set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} /debug:full")

            # to do - make a change in setup or Math.cmake, or here?
            if(DEFINED MKL_FLAG)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_FLAG}")
            endif()

            # windows: none
            if(ENABLE_STATIC_LINKING)
                set(CMAKE_Fortran_FLAGS
                    "${CMAKE_Fortran_FLAGS}"
                    )
            endif()
            
            if(ENABLE_64BIT_INTEGERS)
                set(CMAKE_Fortran_FLAGS
                    "${CMAKE_Fortran_FLAGS} /4I8"
                    )
            endif()
            
            if(ENABLE_BOUNDS_CHECK)
                # miro: removed -ftrapuv due to problems with floating point execptions
                # to do: -fpstkchk for windows
                set(CMAKE_Fortran_FLAGS
                    "${CMAKE_Fortran_FLAGS} /check:bounds /check:pointers /check:uninit /check:output_conversion"
                    )
            endif()
        else()
            set(CMAKE_Fortran_FLAGS         "-w -assume byterecl -DVAR_IFORT -g -traceback")
            set(CMAKE_Fortran_FLAGS_DEBUG   "-O0")
            set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip")
            
            # intel compiler xHost flag
            if(XHOST_FLAG_AVAILABLE)
                set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -xHost")
            endif()
            
            set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")

            if(DEFINED MKL_FLAG)
                set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MKL_FLAG}")
            endif()

            if(ENABLE_STATIC_LINKING)
                set(CMAKE_Fortran_FLAGS
                    "${CMAKE_Fortran_FLAGS} -static-libgcc -static-intel"
                    )
            endif()
            
            if(ENABLE_64BIT_INTEGERS)
                set(CMAKE_Fortran_FLAGS
                    "${CMAKE_Fortran_FLAGS} -i8"
                    )
            endif()
            
            if(ENABLE_BOUNDS_CHECK)
                # miro: removed -ftrapuv due to problems with floating point execptions
                set(CMAKE_Fortran_FLAGS
                    "${CMAKE_Fortran_FLAGS} -check bounds -fpstkchk -check pointers -check uninit -check output_conversion"
                    )
            endif()

            if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
                message("--Switch off warnings due to incompatibility XCode 4 and Intel 11 on OsX 10.6")
                set(CMAKE_Fortran_FLAGS
                    "${CMAKE_Fortran_FLAGS} -Qoption,ld,-w"
                    )
            endif()
        endif()
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES PGI)
        set(CMAKE_Fortran_FLAGS         "-DVAR_PGF90")
        set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
        set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
        
        if(ENABLE_STATIC_LINKING)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} -Bstatic"
                )
        endif()
        
        if(ENABLE_64BIT_INTEGERS)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} -i8"
                )
        endif()
        
        if(ENABLE_BOUNDS_CHECK)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} "
                )
        endif()
        
        if(ENABLE_CODE_COVERAGE)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} "
                )
        endif()
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES XL)
        set(CMAKE_Fortran_FLAGS         "-qzerosize -qextname -qsuppress=cmpmsg")
        set(CMAKE_Fortran_FLAGS_DEBUG   "-g")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
        set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
        
        if(ENABLE_64BIT_INTEGERS)
            set(CMAKE_Fortran_FLAGS
                "${CMAKE_Fortran_FLAGS} -qintsize=8 -q64"
                )
        endif()

        set_source_files_properties(${FREE_FORTRAN_SOURCES}
            PROPERTIES COMPILE_FLAGS
            "-qfree=f90"
            )
        set_source_files_properties(${FIXED_FORTRAN_SOURCES}
            PROPERTIES COMPILE_FLAGS
            "-qfixed"
            )
    else()
        message(FATAL_ERROR "Vendor of your Fortran compiler is not supported")
    endif()

    if(DEFINED EXTRA_Fortran_FLAGS)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_Fortran_FLAGS}")
    endif()
endif()
    
save_compiler_flags(Fortran)
endif()
