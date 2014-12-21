!---------------------------------------------------------------------------------------------
!
! Simple test of Fortran GPU CUDA & CULA libraries
! in the DIRAC environment.
!
! This program measures performances of few  GPU CUBLAS/CULA
! mathematical routines against corresponding CPU BLAS/LAPACK counterparts.
!
! Note:  CULA Dense Free Edition implements popular single precision real and complex functions from LAPACK,
! full commercial version of CULA meets all the LAPACK/BLAS needs of the DIRAC suite,
! see http://www.culatools.com/downloads/dense/
!
! Written by Miro Ilias, Banska Bystrica, March/April 2012
!
!-----------------------------------------------------------------------------------------------

Program GPU_math_tests
  use gpu_stuff
  use fortran_timing

  call fortran_compiler_info
  call c_compiler_info

  call CudaAbouts
  call CulaAbouts

  call print_time_info

  call read_arguments
  select case (method)
    case("sgemm_cuda")
      call sgemm_cuda_test
    case("sgemm_cuda_c")
      call  c_gpu_math_tests(n,method,verb)
    case("sgemm_cula")
      call sgemm_cula_test
    case("sgemm_cula_c")
      call c_gpu_math_tests(n,method,verb)
    case("dgemm_cuda")
      call dgemm_cuda_test
    case("dgemm_cula")
      call dgemm_cula_test
    case("sgemv_cuda")
      call sgemv_cuda_test
    case("sgemv_cuda_c")
      call  c_gpu_math_tests(n,method,verb)
    case("sgemv_cula")
      call sgemv_cula_test
    case("dgemv_cula")
      call dgemv_cula_test
    case("sgemv_cula_c")
      call  c_gpu_math_tests(n,method,verb)
    case("sgeqrf_cula")
      call sgeqrf_test
    case("sgeqrf_cula_c")
      call  c_gpu_math_tests(n,method,verb)
    case("sgesv_cula")
      call sgesv_test
    case("sgesv_cula_c")
      call  c_gpu_math_tests(n,method,verb)
    case default
      write(*,*) 'no method for GPU math testing selected'
      call usage
      stop
  end select
End Program GPU_math_tests
