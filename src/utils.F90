subroutine fortran_compiler_info
! see  http://fortranwiki.org/fortran/show/Predefined+preprocessor+macros ...
write(*,'(A)') '--- Fortran compiler info --- '
#if defined VAR_GFORTRAN
! http://gcc.gnu.org/onlinedocs/cpp/Common-Predefined-Macros.html
write(*,'(A,3(I1,A))') 'GNU, gfortran, version  ',__GNUC__,'.',__GNUC_MINOR__,'.',__GNUC_PATCHLEVEL__
write(*,*) 'version:',__VERSION__
#endif
#if defined VAR_IFORT
!write(*,'(A)') 'Intel Fortran info :'
#if defined __ECC
write(*,*) __ECC
#endif
#if defined __INTEL_COMPILER
write(*,*) 'Intel Fortran, version:',__INTEL_COMPILER
#endif
write(*,*) 'intel compiler build date:', __INTEL_COMPILER_BUILD_DATE
!#include "mkl.fi"
!DEC$IF DEFINED INTEL_MKL_VERSION
!#pragma message "INTEL_MKL_VERSION=",INTEL_MKL_VERSION
!write(*,*) 'Intel MKL  version:',INTEL_MKL_VERSION
!!DEC$IF INTEL_MKL_VERSION .EQ. 100304
!!* 	Code to be conditionally compiled 
!!DEC$ENDIF
!DEC$ENDIF
#endif
end subroutine fortran_compiler_info

subroutine CulaAbouts
implicit none
#if defined USE_CULA
! get info about CULA ...
! see  nm -D /usr/local/cula/lib64/libcula_fortran.so | less
! CULA Programmer’s Guide,  www.culatools.com,  Release R14 (CUDA 4.1)
integer, external ::  CULA_GETVERSION, CULA_GETCUBLASRUNTIMEVERSION, CULA_GETCUDADRIVERVERSION
integer, external :: CULA_GETCUDARUNTIMEVERSION, CULA_GETDEVICEINFO, CULA_GETCUDAMINIMUMVERSION
write(*,'(A)') '--- CULA (and CUDA) GPU info ---'
  write(*,'(1X,A,I5)') 'CULA version=',CULA_GETVERSION()
  write(*,'(1X,A,I5)') 'CULA_CUBLAS runtime version :',CULA_GETCUBLASRUNTIMEVERSION()
  write(*,'(1X,A,I5)') 'CULA_CUDA   runtime version :',CULA_GETCUDARUNTIMEVERSION()
  write(*,'(1X,A,I5)') 'CULA_CUDA   driver version  :', CULA_GETCUDADRIVERVERSION()
 if (CULA_GETCUDARUNTIMEVERSION().LT.CULA_GETCUDAMINIMUMVERSION()) then
   write(*,*) ' CULA_GETCUDAMINIMUMVERSION=',CULA_GETCUDAMINIMUMVERSION()
   write(*,*) "warning: CUDA runtime version < cuda minumum required version !"
 endif
 if (CULA_GETCUBLASRUNTIMEVERSION().LT.CULA_GETCUDAMINIMUMVERSION()) then
   write(*,*) "warning: CUBLAS runtime version < cuda minumum required version !"
   write(*,*) '  cublasRuntimeVersion',cula_GetCublasRuntimeVersion()
 endif
 if (CULA_GETCUDADRIVERVERSION() .LT.CULA_GETCUDAMINIMUMVERSION()) then
   write(*,*) "CUDA driver version < cuda minumum required version !"
 endif
#endif
end subroutine CulaAbouts

subroutine CudaAbouts
! No fortran functions about CUDA version, obtainable only from CULA
! see less /usr/local/cuda/src/fortran.c
implicit none
external :: CudaAbouts_C
call CudaAbouts_C()
end subroutine CudaAbouts
