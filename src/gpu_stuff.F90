module gpu_stuff

  use fortran_timing

implicit none

integer, parameter :: size_of_real  =4
integer, parameter :: size_of_double=8
integer :: AllocateStatus, stat, ios, info

! CUDA specifix library obejcts 
#if USE_CUBLAS
external ::  cublas_init, cublas_free, cublas_set_matrix
external ::  cublas_get_matrix,cublas_set_vector,cublas_get_vector
external ::  cublas_shutdown, cublas_alloc
external ::  cublas_sgemm, cublas_sgemv, cublas_dgemm
integer ::   cublas_alloc, cublas_set_matrix, cublas_get_matrix, cublas_free
integer ::   cublas_set_vector, cublas_get_vector, cublas_sgemm, cublas_sgemv,cublas_dgemm
#endif

! CULA specific library objects 
#if USE_CULA
integer, external  :: cula_initialize, cula_shutdown,cula_sgemv, cula_dgemv
integer, external  :: CULA_SGEQRF, cula_sgesv, cula_sgemm, cula_dgemm
#endif

real(kind=size_of_real),   external :: snrm2 ! for checking  the norm of a vector ..
real(kind=size_of_double), external :: dnrm2 ! for checking  the norm of a vector ..

CHARACTER*100  BUFFER

integer :: verb=0 ! default verbosity level

integer :: n,ii
integer :: i,j

! common allocatble objects - both for CUDA (BLAS) and CULA (BLAS+LAPACK) routines...
real(kind=size_of_real),   allocatable:: A_r(:,:),B_r(:,:),C_r(:,:),x_r(:),y_r(:)
real(kind=size_of_double), allocatable:: A_d(:,:),B_d(:,:),C_d(:,:),x_d(:),y_d(:)
REAL(kind=size_of_real),   ALLOCATABLE :: TAU_r(:), TAU1_r(:)
REAL(kind=size_of_real),   ALLOCATABLE :: WORK_r(:)
integer, allocatable :: ipiv(:), ipiv1(:)

! extra stuff for checking QR-decomposition from SGEQRF ouputs
real(kind=size_of_real),  allocatable:: H_r(:,:),Q_r(:,:),Q1_r(:,:),v_r(:), A1_r(:,:), A2_r(:,:)
real(kind=size_of_real),  allocatable:: b1_r(:),b2_r(:)

#if ARCH_64
! 8-byte pointers on 64-bit architectures for CUDA routines (CULA does not need them)
integer*8 :: devPtrA,devPtrB,devPtrC,devPtrx,devPtry
#else
integer   :: devPtrA,devPtrB,devPtrC,devPtrx,devPtry
#endif

integer, parameter  :: maxnu=20
integer :: no_rout=0
character*15 :: routine(maxnu), method

contains

subroutine usage
 integer :: no_rout1=0
routine(1)="sgemm_cuda"
routine(2)="sgemv_cuda"
routine(3)="dgemm_cuda" 
routine(4)="sgemm_cuda_c"
routine(5)="sgemv_cuda_c"
routine(6)="dgemm_cuda_c"
no_rout=6
#if USE_CULA
 ! cula lapack
 no_rout1=no_rout
routine(no_rout1+1)="sgemm_cula"
routine(no_rout1+2)="sgemv_cula"
routine(no_rout1+3)="sgesv_cula"
routine(no_rout1+4)="sgeqrf_cula"
routine(no_rout1+5)="sgemm_cula_c"
routine(no_rout1+6)="sgemv_cula_c"
routine(no_rout1+7)="dgemm_cula"
routine(no_rout1+8)="dgemv_cula"
routine(no_rout1+9)="sgesv_cula_c"
routine(no_rout1+10)="dgemm_cula_c"
  no_rout=no_rout1+10
#endif
write(*,*) 'usage: gpu_test.x  <routine>  <n> <verbosity level>'
! using own C-function trim 
write(*,*) 'routine: ',(trim(routine(i)),',',i=1,no_rout)
end subroutine usage

subroutine read_arguments
CALL GETARG(1,BUFFER); READ(BUFFER,*,IOSTAT=IOS) method
IF (IOS.NE.0) THEN
    call usage
    stop "error in reading online argument !"
ELSE
    continue
    !print *,'online argument read ... method=',method
ENDIF
CALL GETARG(2,BUFFER); READ(BUFFER,*,IOSTAT=IOS) n
IF (IOS.NE.0) THEN
  call usage
  stop "error in reading online argument !"
ELSE
  continue   
  !print *,'online argument read in ...n=',n
ENDIF
! read in verbosity level - optional argument
CALL GETARG(3,BUFFER); READ(BUFFER,*,IOSTAT=ios) verb
IF (ios.NE.0) THEN
  write(*,'(1x,a,i1)') 'using default verbosity level:',verb
ELSE
  write(*,'(1x,a,i2)') 'online argument read in - verbosity level:',verb
ENDIF
end subroutine read_arguments

subroutine sgemm_cuda_test
write(*,'(2x,a)') 'cuda_sgemm GPU vs CPU test'
allocate (A_r(n,n),B_r(n,n),C_r(n,n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
  stop "error in matrix allocations !"
else
  if (verb.ge.1) then
    write(*,*) "A_r,B_r,C_r matrixes allocated in CPU memory"
  endif
endif
do i=1,n
do j=1,n
  if (i.eq.j) then
    A_r(i,j)=1.0e0
  else
    A_r(i,j)=0.0e0
  endif
  B_r(i,j)=float(i+j) 
  C_r(i,j)=0.0e0
enddo
enddo
!call system_clock( c1, cr, cm )
call cublas_init
stat= cublas_Alloc(n*n,size_of_real, devPtrA)
if (stat .NE. 0) then
  write(*,*) "device memory allocation of devPtrA failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrA OK"
  endif
endif
stat=cublas_Alloc(n*n,size_of_real, devPtrB)
if (stat .NE. 0) then
  write(*,*) "device memory allocation of devPtrB failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrB OK"
  endif
endif
stat=cublas_Alloc(n*n,size_of_real, devPtrC)
if (stat .NE. 0) then
  write(*,*) "device memory allocation of devPtrC failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrC OK"
  endif
endif
! ... start the time measurement
call system_clock( c1, cr, cm )
stat=cublas_Set_Matrix(n,n,size_of_real,A_r,n,devPtrA,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrA)
  write(*,*) "data upload mtx A_r(CPU) -> devPtrA(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx A_r(CPU) -> devPtrA(GPU) : OK"
  endif
endif
stat=cublas_Set_Matrix(n,n,size_of_real,B_r,n,devPtrB,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrB)
  write(*,*) "data upload mtx B_r(CPU) -> devPtrB(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx B_r(CPU) -> devPtrB(GPU) : OK"
  endif
endif
stat=cublas_Set_Matrix(n,n,size_of_real,C_r,n,devPtrC,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrC)
  write(*,*) "data upload mtx C_r(CPU) -> devPtrC(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx C_r(CPU) -> devPtrC(GPU) : OK"
  endif
endif
stat=cublas_SGEMM('n','n', n,n,n,1.0e0,devPtrA,n,devPtrB,n,1.0e0,devPtrC,n)
if (stat .NE. 0) then
  stop "cublas_SGEMM failed!"
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx C_r(CPU) -> devPtrC(GPU) : OK"
  endif
endif
stat=cublas_Get_Matrix(n,n,size_of_real,devPtrC,n,C_r,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrC)
write(*,*) "data upload failed"
call cublas_shutdown
stop
endif
call system_clock( count=c2 )
call print_time("GPU CUDA SGEMM")
call mtx_check(C_r)
stat=cublas_free(devPtrA)
stat=cublas_free(devPtrB)
stat=cublas_free(devPtrC)
call cublas_shutdown
! now CPU sgemm...
call system_clock( c1, cr, cm )
call sgemm('n','n',n,n,n,1.0e0,A_r,n,B_r,n,-2.0e0,C_r,n)
call system_clock( count=c2 )
call print_time("CPU BLAS SGEMM")
call mtx_check(C_r)
call cublas_shutdown
end subroutine sgemm_cuda_test

subroutine sgemm_cula_test
#if defined USE_CULA
write(*,'(2x,a)') 'cula_sgemm GPU vs CPU test'
allocate (A_r(n,n),B_r(n,n),C_r(n,n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
   stop "error in main matrix allocations !"
else
   if (verb.ge.1) then
     write(*,*) "A_d,B_d,C_d matrixes allocated in CPU memory"
   endif
endif
do i=1,n
do j=1,n
 if (i.eq.j) then ! A_d is unit matrix
   A_r(i,j)=1.0e0
 else
   A_r(i,j)=0.0e0
 endif
 B_r(i,j)=dfloat(i+j) ! B_d is symmetric matrix
 C_r(i,j)=0.0e0
enddo
enddo
! start the time measurement
call system_clock( c1, cr, cm )
stat = CULA_INITIALIZE(); call check_cula_status(stat)
stat = CULA_SGEMM('n','n',n,n,n,1.0e0,A_r,n,B_r,n,1.0e0,C_r,n)
call system_clock( count=c2 )
call print_time("GPU CULA_SGEMM")
call check_cula_status(stat)
stat = CULA_SHUTDOWN(); call check_cula_status(stat)
! ---
call mtx_check(C_r)
! now CPU sgemm...
call system_clock( c1, cr, cm )
call sgemm('n','n',n,n,n,1.0e0,A_r,n,B_r,n,-2.0e0,C_r,n)
call system_clock( count=c2 )
call print_time("CPU BLAS SGEMM")
call mtx_check(C_r)
#else
  print *,'CULA library not active !'
#endif
end subroutine sgemm_cula_test

subroutine dgemm_cula_test
#if defined USE_CULA
write(*,'(2x,a)') 'cula_dgemm GPU vs CPU test'
allocate (A_d(n,n),B_d(n,n),C_d(n,n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
   stop "error in dgemm_cula_test allocations !"
else
   if (verb.ge.1) then
     write(*,*) "A_d,B_d,C_d matrixes allocated in CPU memory for dgemm_cula_test"
   endif
endif
do i=1,n
do j=1,n
 if (i.eq.j) then ! A_d is unit matrix
   A_d(i,j)=1.0e0
 else
   A_d(i,j)=0.0e0
 endif
 B_d(i,j)=dfloat(i+j) ! B_d is symmetric matrix
 C_d(i,j)=0.0e0
enddo
enddo
! start the time measurement
call system_clock( c1, cr, cm )
stat = CULA_INITIALIZE(); call check_cula_status(stat)
stat = CULA_DGEMM('n','n',n,n,n,1.0d0,A_d,n,B_d,n,1.0e0,C_d,n); call check_cula_status(stat)
call system_clock( count=c2 )
call print_time("GPU CULA_DGEMM")
stat = CULA_SHUTDOWN(); call check_cula_status(stat)
! ---
call mtx_check_d(C_d)
! now CPU sgemm...
call system_clock( c1, cr, cm )
call dgemm('n','n',n,n,n,1.0d0,A_d,n,B_d,n,-2.0e0,C_d,n)
call system_clock( count=c2 )
call print_time("CPU BLAS DGEMM")
call mtx_check_d(C_d)
#else
  print *,'CULA library not active !'
#endif
end subroutine dgemm_cula_test


subroutine dgemm_cuda_test
write(*,'(2x,a)') 'cuda_dgemm GPU vs CPU test'
allocate (A_d(n,n),B_d(n,n),C_d(n,n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
   stop "error in main matrix allocations !"
else
   if (verb.ge.1) then
     write(*,*) "A_d,B_d,C_d matrixes allocated in CPU memory"
   endif
endif
do i=1,n
do j=1,n
 if (i.eq.j) then ! A_d is unit matrix
   A_d(i,j)=1.0e0
 else
   A_d(i,j)=0.0e0
 endif
 B_d(i,j)=dfloat(i+j) ! B_d is symmetric matrix
 C_d(i,j)=0.0e0
enddo
enddo
!call system_clock( c1, cr, cm )
call cublas_init
stat= cublas_Alloc(n*n,size_of_double, devPtrA)
if (stat .NE. 0) then
write(*,*) "device memory allocation of devPtrA failed"
call cublas_shutdown
stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrA done"
  endif
endif
stat=cublas_Alloc(n*n,size_of_double,devPtrB)
if (stat .NE. 0) then
write(*,*) "device memory allocation of devPtrB failed"
call cublas_shutdown
stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrB done"
  endif
endif
stat=cublas_Alloc(n*n,size_of_double, devPtrC)
if (stat .NE. 0) then
write(*,*) "device memory allocation of devPtrC failed"
call cublas_shutdown
stop
else
  if (verb.ge.1) then
    write(*,*) "device memory allocation of devPtrC done"
  endif
endif
! start the time measurement
call system_clock( c1, cr, cm )
stat=cublas_Set_Matrix(n,n,size_of_double,A_d,n,devPtrA,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrA)
  write(*,*) "data upload mtx A_d(CPU) -> devPtrA(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx A_d(CPU) -> devPtrA(GPU) : OK"
  endif
endif
stat=cublas_Set_Matrix(n,n,size_of_double,B_d,n,devPtrB,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrB)
  write(*,*) "data upload mtx B_d(CPU) -> devPtrA(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx B_d(CPU) -> devPtrB(GPU) : OK"
  endif
endif
stat=cublas_Set_Matrix(n,n,size_of_double,C_d,n,devPtrC,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrC)
  write(*,*) "data upload mtx C_d(CPU) -> devPtrC(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx C_d(CPU) -> devPtrC(GPU) : OK"
  endif
endif
stat=cublas_DGEMM('n','n', n,n,n,1.0d0,devPtrA,n,devPtrB,n,1.0d0,devPtrC,n)
if (stat .NE. 0) then
  stop "cublas_DGEMM failed!"
else
  if (verb.ge.1) then
    write(*,*) "cublas_DGEMM done OK"
  endif
endif
stat=cublas_Get_Matrix(n,n,size_of_double,devPtrC,n,C_d,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrC)
  write(*,*) "data download GPU devPtrC -> CPU mtx C_d failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data download GPU devPtrC -> CPU mtx C_d OK"
  endif
endif
call system_clock( count=c2 )
call print_time("GPU CUDA DGEMM")
call mtx_check_d(C_d)
stat=cublas_free(devPtrA)
stat=cublas_free(devPtrB)
stat=cublas_free(devPtrC)
call cublas_shutdown
! now CPU sgemm...
call system_clock( c1, cr, cm )
call dgemm('n','n',n,n,n,1.0d0,A_d,n,B_d,n,-2.0e0,C_d,n)
call system_clock( count=c2 )
call print_time("CPU BLAS DGEMM")
call mtx_check_d(C_d)
call cublas_shutdown
end subroutine dgemm_cuda_test

subroutine sgemv_cuda_test
write(*,'(2x,a)') 'cuda_sgemv GPU vs CPU test'
allocate (A_r(n,n),x_r(n),y_r(n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
  stop "error in allocations !"
else
  if (verb.ge.1) then
    write(*,*) "A_r,x_r,y_r objects allocated in CPU memory"
  endif
endif
! fill allocated objects
do i=1,n
  do j=1,n
    if (i.eq.j) then
      A_r(i,j)=1.0e0
    else
     A_r(i,j)=0.0e0
    endif
  enddo
  x_r(i)=float(1)
  y_r(i)=0.0e0
enddo
!call system_clock( c1, cr, cm )
call cublas_init
stat=cublas_Alloc(n*n,size_of_real, devPtrA)
if (stat .NE. 0) then
  write(*,*) "GPU device memory allocation of devPtrA failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "GPU device memory allocation of devPtrA OK"
  endif
endif
stat= cublas_Alloc(n,size_of_real,devPtrx)
if (stat .NE. 0) then
  write(*,*) "GPU device memory allocation of devPtrx failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "GPU device memory allocation of devPtrx OK"
  endif
endif
stat=cublas_Alloc(n,size_of_real,devPtry)
if (stat .NE. 0) then
  write(*,*) "GPU device memory allocation of devPtry failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "GPU device memory allocation of devPtry OK"
  endif
endif
! start the time measurement
call system_clock( c1, cr, cm )
stat=cublas_Set_Matrix(n,n,size_of_real,A_r,n,devPtrA,n)
if (stat .NE. 0) then
  stat=cublas_free (devPtrA)
  write(*,*) "data upload mtx A_r(CPU) -> devPtrA(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload mtx A_r(CPU) -> devPtrA(GPU) : OK"
  endif
endif
stat=cublas_set_vector(n,size_of_real,x_r,1,devPtrx,1)
if (stat .NE. 0) then
  stat=cublas_free (devPtrx)
  write(*,*) "data upload vector x_r(CPU) -> devPtrx(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload vector x_r(CPU) -> devPtrx(GPU) : OK"
  endif
endif
stat=cublas_set_vector(n,size_of_real,y_r,1,devPtry,1)
if (stat .NE. 0) then
  stat=cublas_free (devPtry)
  write(*,*) "data upload vector y_r(CPU) -> devPtry(GPU) failed"
  call cublas_shutdown
  stop
else
  if (verb.ge.1) then
    write(*,*) "data upload vector y_r(CPU) -> devPtry(GPU) : OK"
  endif
endif
!y = beta*y+alpha*A*x
stat=cublas_sgemv('n',n,n,1.0e0,devPtrA,n,devPtrx,1,0.0e0,devPtry,1)
if (stat .NE. 0) then
  stop "cublas_sgemv failed!"
else
  if (verb.ge.1) then
    write(*,*) 'cublas_sgemv done'
  endif
endif
! get the vector from GPU (devPtry) into y_r...
stat=cublas_get_vector(n,size_of_real,devPtry,1,y_r,1)
call system_clock( count=c2 )
call print_time("GPU sgemv")
write(*,*) 'GPU cublas_sgemv vector control output: sum,norm=',sum(y_r),snrm2(n,y_r,1)
stat=cublas_free(devPtrA)
stat=cublas_free(devPtrx)
stat=cublas_free(devPtry)
call system_clock( c1, cr, cm )
!now perform CPU lib routine:  y = beta*y+alpha*A*x
call sgemv('n',n,n,-2.0e0,A_r,n,x_r,1,1.0e0,y_r,1)
call system_clock( count=c2 )
call print_time("CPU sgemv")
write(*,*) 'CPU sgemv vector control output: sum,norm=',sum(y_r),snrm2(n,y_r,1)
end subroutine sgemv_cuda_test

subroutine sgemv_cula_test
#if defined USE_CULA
write(*,'(2x,a)') 'cula_sgemv GPU vs CPU test'
allocate (A_r(n,n),x_r(n),y_r(n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
  stop "error in allocations !"
else
  if (verb.ge.1) then
    write(*,*) "A_r,x_r,y_r objects allocated in CPU memory"
  endif
endif
! fill allocated objects
do i=1,n
do j=1,n
 if (i.eq.j) then
   A_r(i,j)=1.0e0
 else
   A_r(i,j)=0.0e0
 endif
enddo
 x_r(i)=float(1)
 y_r(i)=0.0e0
enddo
call system_clock( c1, cr, cm )
stat = CULA_INITIALIZE(); call check_cula_status(stat)
stat=cula_sgemv('n',n,n,1.0e0,A_r,n,x_r,1,0.0e0,y_r,1)
call system_clock( count=c2 )
call print_time("GPU CULA SGEMV")
stat = CULA_SHUTDOWN(); call check_cula_status(stat)
write(*,*) 'GPU cula_sgemv vector control output: sum,norm=',sum(y_r),snrm2(n,y_r,1)
! ---
call system_clock( c1, cr, cm )
!now perform CPU lib routine:  y = beta*y+alpha*A*x
call sgemv('n',n,n,-2.0e0,A_r,n,x_r,1,1.0e0,y_r,1)
call system_clock( count=c2 )
call print_time("CPU BLAS SGEMV")
write(*,*) 'CPU sgemv vector control output: sum,norm=',sum(y_r),snrm2(n,y_r,1)
#else
 print *,'sorry, CULA library is not active'
#endif
end subroutine sgemv_cula_test

subroutine dgemv_cula_test
#if defined USE_CULA
write(*,'(2x,a)') 'cula_dgemv GPU vs CPU test'
allocate (A_d(n,n),x_d(n),y_d(n),STAT=AllocateStatus)
if (AllocateStatus.ne.0) then
  stop "error in matrix,vectors allocations in dgemv_cula_test !"
else
  if (verb.ge.1) then
    write(*,*) "A_d,x_d,y_d objects allocated in CPU memory"
  endif
endif
! fill allocated objects
do i=1,n
do j=1,n
 if (i.eq.j) then
   A_d(i,j)=1.0e0
 else
   A_d(i,j)=0.0e0
 endif
enddo
 x_d(i)=dfloat(1)
 y_d(i)=0.0d0
enddo
call system_clock( c1, cr, cm )
stat = CULA_INITIALIZE(); call check_cula_status(stat)
stat=cula_dgemv('n',n,n,1.0d0,A_d,n,x_d,1,0.0d0,y_d,1); call check_cula_status(stat)
call system_clock( count=c2 )
call print_time("GPU CULA SGEMV")
stat = CULA_SHUTDOWN(); call check_cula_status(stat)
write(*,*) 'GPU cula_dgemv vector control output: sum,norm=',sum(y_d),dnrm2(n,y_d,1)
! ---
call system_clock( c1, cr, cm )
!now perform CPU lib routine:  y = beta*y+alpha*A*x
call dgemv('n',n,n,-2.0d0,A_d,n,x_d,1,1.0d0,y_d,1)
call system_clock( count=c2 )
call print_time("CPU BLAS DGEMV")
write(*,*) 'CPU dgemv vector control output: sum,norm=',sum(y_d),dnrm2(n,y_d,1)
#else
 print *,'sorry, CULA library is not active'
#endif
end subroutine dgemv_cula_test


subroutine sgeqrf_test
#if defined USE_CULA
! test on matrix QR decomposition ....
  real(kind=size_of_real) :: diffs = 0.0e0
#if defined VAR_IFORT
  real(kind=size_of_real), external :: rand
#endif
  write(*,'(2x,a)') '=== Going to perform GPU cula_sgeqrf vs CPU lapack  test ==='
  allocate (A_r(n,n),B_r(n,n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in allocations of A_r,B_r !"
  else
    if (verb.ge.1) then
      write(*,*) "A_r, B_r matrixes allocated in CPU memory"
    endif
  endif
  allocate(WORK_r(n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in WORK_r array allocation"
  else
    if (verb.ge.1) then
      print *,"allocated array, WORK_r(LWORK)"
    endif
  endif
  allocate(TAU_r(n),TAU1_r(n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in TAU_r, TAU1_r array allocation"
  else
    if (verb.ge.1) then
      print *,"allocated auxiliary arrays, TAU_r(N), TAU1_r(N),N=",n
    endif
  endif
! allocate extra stuff for throught checking of result:
  if (verb.ge.1) then
     allocate(H_r(n,n),Q_r(n,n),Q1_r(n,n),v_r(n),A1_r(n,n),A2_r(n,n),STAT=AllocateStatus)
     if (AllocateStatus.ne.0) then
        stop "error in H_r(n,n),Q_r(n,n),Q1_r(n,n),v_r(n),A1_r(n,n),A2_r(n,n) allocations"
     else
        print *,"H_r(n,n),Q_r(n,n),Q1_r(n,n),v_r(n),A1_r(n,n),A2_r(n,n) allocations OK." 
     endif
  endif
! fill allocated objects A_r, B_r with (random?) data for later checking
  if (verb.ge.3) then
    print *,'    Entering matrix for QR decomposition:'
  endif
  write(*,'(2X,A)') '- matrix for QR decomposition elements are [i,j]=i+j'
  do j=1,n
  do i=1,n
     ! symmetric matrix ....
     ! A_r(i,j) = rand(0)
     A_r(i,j) = float(i)+float(j) 
     ! B the same as A
     B_r(i,j)=A_r(i,j)
     if (verb.ge.1) then ! save mtx for further tests of GPU/CPU data ...
        A1_r(i,j) = A_r(i,j) 
        A2_r(i,j) = A_r(i,j)
     endif
     if (verb.ge.3) then
       write(*,"(7X,A,i2,A,i2,a,f12.7)")  "A[",i,",",j,"]=",A_r(i,j)
     endif
  enddo
  enddo
  !-----------
  stat = CULA_INITIALIZE(); call check_cula_status(stat)
  ! ... start the time measurement
  call system_clock( c1, cr, cm )
  stat = CULA_SGEQRF(N, N, A_r, N, TAU_r)
  call system_clock( count=c2 )
  call check_cula_status(stat)
  call print_time("GPU CULA_SGEQRF")

  stat = CULA_SHUTDOWN(); call check_cula_status(stat)

  ! ... start again the time measurement
  call system_clock( c1, cr, cm )
! time for lapack SGEQRF
!see  http://www.netlib.org/lapack/single/sgeqrf.f
  CALL SGEQRF(N, N, B_r, N, TAU1_r, WORK_r, N , stat)
  call system_clock( count=c2 )
  IF (stat.ne.0) then
    stop "error in CPU LAPACK SGEQRF..."
  else
    if (verb.ge.0) then
      print *,"LAPACK SGEQRF finished successfully ...."
    endif
  endif
  call print_time("CPU LAPACK SGEQRF")
! check numbers in A_r, B_r - must be the same...
  if (verb.ge.3) then
     write(*,"(4X,A)") "---------------------------------------------------------------"
     write(*,"(4X,A)")  "QR results:     GPU-SGEQRF    CPU-SGEQRF   abs.diff(GPU-CPU)"
     write(*,"(4X,A)") "---------------------------------------------------------------"
  endif

  do i=1,n
  do j=i,n
    diffs = diffs + abs(A_r(i,j)-B_r(i,j))  
    if (verb.ge.3) then
       write(*,'(7X,A,I2,A,I2,A,F12.6,A,F12.6,A,F12.6)') "R[",i,",",j,"]=",A_r(i,j),"  ",B_r(i,j),":", abs(A_r(i,j)-B_r(i,j))
    endif
  enddo
  enddo
  if (verb.ge.3) then
    write(*,"(5X,A)") "---------  ONLY INFORMATIVE RESULTS !!!  ------------"
  endif
  do i=2,n
  do j=1,i
    diffs = diffs + abs(A_r(i,j)-B_r(i,j))  
    if (verb.ge.3.and.i.ne.j) then
       write(*,'(4X,A,I2,A,I2,A,F12.6,A,F12.6,A,F12.6)') "e.r.[",i,",",j,"]=",A_r(i,j),"  ",B_r(i,j),":", abs(A_r(i,j)-B_r(i,j))
    endif
  enddo
  enddo
  if (verb.ge.3) then
    write(*,"(5X,A)") "-----  ONLY INFORMATIVE RESULT !!!    ---------------"
  endif
  do i=1, n
    diffs =  diffs + abs(TAU_r(i)-TAU1_r(i))
    if (verb.ge.3) then
      write(*,'(8X,A,I2,A,F12.6,A,F12.6,A,F12.6)') "TAU[",i,"]=",TAU_r(i),"  ",TAU1_r(i),":",abs(TAU_r(i)-TAU1_r(i))
    endif
  enddo
  if (verb.ge.3) then
    write(*,"(5X,A)") "-------  ONLY INFORMATIVE RESULT !!!    -------------"
    write(*,*) 'CULA_SGEQRF vs lapack SGEQRF averaged absolute diffs are ',diffs/float(n*n)
  endif

  if (verb.ge.1) then
      write(*,"(1x,a)")  & 
"--- Going to do a thorough check of QR decomposed matrixes by using CPU data/routines (shall take some time ) ---"
      ! make Q(i) unit matrix
      do i = 1,n
      do j = 1,n
        if (i.eq.j) then
          Q_r(i,j)=1.0e0
        else
          Q_r(i,j)=0.0e0
        endif
      enddo
      enddo

      ! own computing H(i) = I - tau(i) * v[i] * v[i]^T and Q(i) = H(1)H(2)...
    do ii = 1,n
      ! compute H(i) = I - tau(i) * v * v^T, where v(1:i-1) = 0 and v(i) = 1; v(i+1:n)=A(i+1:n,i)
      do j = 1, ii-1
        v_r(j) = 0.0e0
      enddo
      v_r(ii) = 1.0e0
      do j = ii+1,n
        v_r(j)=B_r(j,ii)  ! CPU data
      enddo
      ! make H(i) unit matrix
      do i = 1,n
      do j = 1,n
        if (i.eq.j) then
          H_r(i,j)=1.0e0
        else
          H_r(i,j)=0.0e0
        endif
      enddo
      enddo
      ! own computing H(i) = I - tau(i) * v[i] * v[i]^T
      call sgemm('n','t',n,n,1,-TAU1_r(ii),v_r,n,v_r,n,+1.0e0,H_r,n)
      if (verb.ge.5) then
        write(*,'(3X,A,I2,A)') "current vector v[",ii,"]:"
        do i=1,n
          write(*,*) i,v_r(i)
        enddo
        ! ------------------ !
        write(*,'(3X,A,I2,A)') "current matrix H(",ii,") = I - tau(i) * v[i] * v[i]^T:"
        do i=1,n
        do j=1,n
          write(*,*) i,j,H_r(i,j)
        enddo
        enddo
      endif
      ! compute  Q = H(1) H(2) . . . H(n), or Q1 := Q*H(i) plus Q:=Q1
      call sgemm('n','n',n,n,n,+1.0e0,Q_r,n,H_r,n,+0.0e0,Q1_r,n)
      ! (n, x, incx, y, incy)
      call scopy(n*n,Q1_r,1,Q_r,1)
      if (verb.ge.5) then
        write(*,'(3X,A,I2,A)') "current matrix Q(",ii,"):"
        do i=1,n
        do j=1,n
          write(*,*) i,j,Q_r(i,j)
        enddo
        enddo
      endif
    enddo
    ! orthonormality check of Q^T* Q
    call sgemm('t','n',n,n,n,+1.0e0,Q_r,n,Q_r,n,+0.0e0,Q1_r,n)
    write(*,'(3X,A)') "orthoronormality check of Q^T*Q=I:"
    if (verb.ge.5) then
        do i=1,n
        do j=1,n
          write(*,*) i,j,Q1_r(i,j)
        enddo
        enddo
    endif
    call check_unit_mtx_r(Q1_r)
    ! orthonormality check of Q* Q^T
    call sgemm('n','t',n,n,n,+1.0e0,Q_r,n,Q_r,n,+0.0e0,Q1_r,n)
    write(*,'(3X,A)') "orthoronormality check of Q*Q^T=I:"
    if (verb.ge.5) then
        do i=1,n
        do j=1,n
          write(*,*) i,j,Q1_r(i,j)
        enddo
        enddo
    endif
    call check_unit_mtx_r(Q1_r)
    ! Finally check A = Q R
    ! frst zero lower triangle of R contained in B_r - CPU data
    do i=2,n
    do j=1,i
      if (i.ne.j) then
        B_r(i,j) = 0.0e0
      endif
    enddo
    enddo
    ! Finally check A = Q R, what is A1_r = Q_r * B_r - A1_r = 0
    write(*,'(3X,A)') "--- FINAL CHECK of CPU SGEQRF results through A-QR=0: ---"
    call sgemm('n','n',n,n,n,+1.0e0,Q_r,n,B_r,n,-1.0e0,A1_r,n)
    call mtx_check(A1_r)

    write(*,'(2X,A)') 'Going to do a thorough check A = QR by using GPU data/routines (shall take some time):'
      ! make Q(i) unit matrix
      do i = 1,n
      do j = 1,n
        if (i.eq.j) then
          Q_r(i,j)=1.0e0
        else
          Q_r(i,j)=0.0e0
        endif
      enddo
      enddo
    ! own computing H(i) = I - tau(i) * v[i] * v[i]^T and Q(i) = H(1)H(2) ...
    do ii = 1,n
      ! compute H(i) = I - tau(i) * v * v^T, where v(1:i-1) = 0 and v(i) = 1; v(i+1:n)=A(i+1:n,i)
      do j = 1, ii-1
        v_r(j) = 0.0e0
      enddo
      v_r(ii) = 1.0e0
      do j = ii+1,n
        v_r(j)=A_r(j,ii)  ! GPU data
      enddo
      ! initialize  H(i)  as  unit matrix
      do i = 1,n
      do j = 1,n
        if (i.eq.j) then
          H_r(i,j)=1.0e0
        else
          H_r(i,j)=0.0e0
        endif
      enddo
      enddo
      ! own computing H(i) = I - tau(i) * v[i] * v[i]^T
      call sgemm('n','t',n,n,1,-TAU_r(ii),v_r,n,v_r,n,+1.0e0,H_r,n) ! using GPU data - TAU_r
      if (verb.ge.5) then
        write(*,'(3X,A,I2,A)') "current vector v[",ii,"]:"
        do i=1,n
          write(*,*) i,v_r(i)
        enddo
        ! ------------------ !
        write(*,'(3X,A,I2,A)') "current matrix H(",ii,") = I - tau(i) * v[i] * v[i]^T:"
        do i=1,n
        do j=1,n
          write(*,*) i,j,H_r(i,j)
        enddo
        enddo
      endif
      ! compute  Q = H(1) H(2) . . . H(n), or Q1 := Q*H(i) plus Q:=Q1
      call sgemm('n','n',n,n,n,+1.0e0,Q_r,n,H_r,n,+0.0e0,Q1_r,n)
      ! (n, x, incx, y, incy)
      call scopy(n*n,Q1_r,1,Q_r,1)
      if (verb.ge.5) then
        write(*,'(3X,A,I2,A)') "current matrix Q(",ii,"):"
        do i=1,n
        do j=1,n
          write(*,*) i,j,Q_r(i,j)
        enddo
        enddo
      endif
    enddo
    ! orthonormality check of Q^T* Q
    stat = CULA_INITIALIZE(); call check_cula_status(stat)
    stat = cula_sgemm('t','n',n,n,n,+1.0e0,Q_r,n,Q_r,n,+0.0e0,Q1_r,n); call check_cula_status(stat)
    !stat = CULA_SHUTDOWN(); call check_cula_status(stat)
    write(*,'(3X,A)') "orthoronormality check of Q^T*Q=I:"
    if (verb.ge.5) then
        do i=1,n
        do j=1,n
          write(*,*) i,j,Q1_r(i,j)
        enddo
        enddo
    endif
    call check_unit_mtx_r(Q1_r)
    ! orthonormality check of Q* Q^T
    stat = cula_sgemm('n','t',n,n,n,+1.0e0,Q_r,n,Q_r,n,+0.0e0,Q1_r,n); call check_cula_status(stat)
    write(*,'(3X,A)') "orthoronormality check of Q*Q^T=I:"
    if (verb.ge.5) then
        do i=1,n
        do j=1,n
          write(*,*) i,j,Q1_r(i,j)
        enddo
        enddo
    endif
    call check_unit_mtx_r(Q1_r)
    ! Finally check A = Q R
    ! frst zero lower triangle of R contained in A_r (GPU data)
    do i=2,n
    do j=1,i
      if (i.ne.j) then
       A_r(i,j) = 0.0e0
      endif
    enddo
    enddo
    ! Finally check A = Q R, what is A1_r = Q_r * A_r - A2_r = 0
    write(*,'(3X,A)') "--- FINAL CHECK of GPU CULA_SGEQRF results through A-QR=0: ---"
    stat = cula_sgemm('n','n',n,n,n,+1.0e0,Q_r,n,A_r,n,-1.0e0,A2_r,n);call check_cula_status(stat)
    stat = CULA_SHUTDOWN(); call check_cula_status(stat)
    call mtx_check(A1_r)

  endif
#else
 print *,'sorry, CULA library is not activated (through -D USE_CULA)'
#endif
end subroutine sgeqrf_test

subroutine sgesv_test
! test of solving linear equations;  A.x=b
#if defined USE_CULA
#if defined VAR_IFORT
  real(kind=size_of_real), external :: rand
#endif
  real(kind=size_of_real) :: diffs = 0.0e0
  write(*,'(2x,a)') 'cula_sgesv GPU vs lapack CPU test'
  allocate (A_r(n,n),B_r(n,n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in allocations of A_r,B_r !"
  else
    if (verb.ge.1) then
      write(*,*) "A_r, B_r matrixes allocated in CPU memory"
    endif
  endif
  allocate(TAU_r(n),TAU1_r(n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in TAU_r, TAU1_r array allocation"
  else
    if (verb.ge.1) then
      print *,"allocated real*4 arrays TAU_r(n), TAU1_r(n) for solving linear equations"
    endif
  endif
  allocate(ipiv(n),ipiv1(n),STAT=AllocateStatus)
  if (AllocateStatus.ne.0) then
    stop "error in ipiv, ipiv1 integer arrays allocation"
  else
    if (verb.ge.1) then
      print *,"allocated integer arrays ipiv(n), ipiv1(n) for solving linear equations"
    endif
  endif

  if (verb.ge.1) then
     allocate(A1_r(n,n),A2_r(n,n),b1_r(n),b2_r(n),STAT=AllocateStatus)
     if (AllocateStatus.ne.0) then
        stop "error in A1_r(n,n),A2_r(n,n),b1_r(n),b2_r(n) allocations"
     else
        print *,"A1_r(n,n),A2_r(n,n),b1_r(n),b2_r(n) additional allocations OK." 
     endif
  endif

! fill allocated objects A_r, B_r, TAU1_r, TAU_r with random data for later checking
  !write(*,"(2X,A)") "Matrix A[i,j]=-i+j*j and vector b[i]=i*i for solving linear equations Ax=b:"
  !write(*,"(2X,A)") "Matrix A[i,j]=i+j and vector b[i]=i for solving linear equations Ax=b:"
  write(*,"(2X,A)") "Matrix A[i,j]=rand(0) and vector b[i]=rand(0) for solving linear equations Ax=b:"
  do i=1,n
  do j=1,n
     A_r(i,j) = rand(0);  ! <= ONLY RANDOM WORKS !!!
     !A_r(i,j) = float(i)+float(j);  
     !A_r(i,j) = float(i)-(float(j)*float(j)); 
     !A_r(i,j) = float(i)-(float(j)*float(j))-(float(i)*float(j)); 
     B_r(i,j)=A_r(i,j)
     if (verb.ge.3) then
        write(*,"(7X,A,i2,A,i2,a,f12.7)")  "A[",i,",",j,"]=",A_r(i,j)
     endif
     if (verb.ge.1) then
       A1_r(i,j)=A_r(i,j);A2_r(i,j)=A_r(i,j)
     endif
  enddo
     TAU_r(i) = rand(0) ! <= ONLY RANDOM WORKS!!!
     !TAU_r(i) = float(i)*float(i)-1
     !TAU_r(i) = 2*float(i)*float(i)-1
     TAU1_r(i) = TAU_r(i)
     if (verb.ge.3) then
        write(*,"(7X,A,I2,A,f12.7)") "b[",i,"]=",TAU_r(i)
     endif
     if (verb.ge.1) then
       b1_r(i) = TAU_r(i); b2_r(i) = TAU_r(i)
     endif
  enddo
  stat = CULA_INITIALIZE(); call check_cula_status(stat)

  ! ... start the time measurement
  call system_clock( c1, cr, cm )
  stat = cula_sgesv(n, 1 , A_r, n , ipiv, TAU_r, n, info )
  call system_clock( count=c2 )
  call check_cula_status(stat)
  call print_time("GPU CULA_SGESV")
  stat = CULA_SHUTDOWN(); call check_cula_status(stat)

  ! now CPU lapack sgesv with time measurement...
  call system_clock( c1, cr, cm )
  call sgesv(n, 1 , B_r, n , ipiv1, TAU1_r, n, stat )
  call system_clock( count=c2 )
  IF (stat.ne.0) then
    print *,'CPU LAPACK SGESV returning error argument (check manual):',stat 
    stop "error in CPU LAPACK SGESV ..."
  else
    if (verb.ge.1) then
      print *,"LAPACK SGESV finished successfully ..."
    endif
  endif
  
  call print_time("CPU LAPACK SGESV")
  ! check  results - roots - must be identical 
  if (verb.ge.2) then
    !print *," roots and their differences  "
    write(*,"(7X,A)")  "root #     GPU           CPU           diff"
  endif
  do i=1,n
    diffs = diffs + abs(TAU_r(i)-TAU1_r(i))  
    if (verb.ge.3) then
      print *,i,TAU_r(i),TAU1_r(i),TAU_r(i)-TAU1_r(i)
    endif
  enddo
  write(*,"(1X,A,D11.6)") & 
  '*** FINAL CONTROL OF RESULTS - roots GPU CULA_SGESV vs CPU SGESV averaged absolute diffs are ',diffs/float(n)
  ! -------------------------------------------------------- !
  ! ---- HARD CHECK OF ROOTS against relation A.x-b=0   ----
  ! -------------------------------------------------------- !
  if (verb.ge.1) then
    !now perform CPU lib routine:  y = beta*y+alpha*A*x
    call sgemv('n',n,n,+1.0e0,A1_r,n,TAU1_r,1,-1.0e0,b1_r,1) ! CPU - vector b1_r should be zero
    write(*,"(2X,A,D11.6)") 'CPU SGESV Ax-b=0 vector control output, norm: ',snrm2(n,b1_r,1)

    stat = CULA_INITIALIZE(); call check_cula_status(stat)
    stat=cula_sgemv('n',n,n,+1.0e0,A2_r,n,TAU_r,1,-1.0e0,b2_r,1) ! GPU - vector b2_r should be zero
    stat = CULA_SHUTDOWN(); call check_cula_status(stat)
    write(*,"(2X,A,D11.6)") 'GPU SGESV Ax-b=0 vector control output, norm: ',snrm2(n,b2_r,1)
  endif

#else
 print *,'sorry, CULA library (at least free version) is not active'
#endif
end subroutine sgesv_test

#if defined USE_CULA
SUBROUTINE check_cula_status(STATUS)
! TOOD: make it more general wrt cula&lapack routines exit codes
         INTEGER, intent(in) :: STATUS
         INTEGER INFO
         external :: CULA_GETERRORINFO
         INTEGER  :: CULA_GETERRORINFO
         INFO = CULA_GETERRORINFO()
         IF (STATUS .NE. 0) THEN
            write(*,'(2x,a,$)') "CULA system status: "
            IF (STATUS .EQ. 7) THEN
!*              culaArgumentError
               WRITE(*,*) 'Invalid value for parameter ', INFO
            ELSE IF (STATUS .EQ. 8) THEN
!*              culaDataError
               WRITE(*,*) 'Data error (', INFO ,')'
            ELSE IF (STATUS .EQ. 9) THEN
!*              culaBlasError
               WRITE(*,*) 'Blas error (', INFO ,')'
            ELSE IF (STATUS .EQ. 10) THEN
!*              culaRuntimeError
               WRITE(*,*) 'Runtime error (', INFO ,')'
            ELSE
!*              others
               call CULA_GETSTATUSSTRING(STATUS)
            ENDIF
            STOP 1
         END IF
END subroutine check_cula_status
#endif

subroutine mtx_check(Mtx)
! simple check of the real matrix
real*4, intent(in) :: Mtx(n,n)
real*4 ::  diag,offdiag
diag=0.0e0;offdiag=diag
do i=1,n
do j=1,n
  if (i.eq.j) then
    diag = diag + Mtx(i,j)
  else
    offdiag = offdiag + Mtx(i,j)
  endif
enddo
enddo
write(*,'(1X,A,D6.2,A,D6.2)') & 
 'aver. sum of diag. elem.:',diag/float(n),' aver.sum of offdiag elem.:',offdiag/(float(n*n)-float(n))
write(*,'(1X,4(A,D12.4,1X))') '1,1:',Mtx(1,1),'1,n:',Mtx(1,n),'n,1:',Mtx(n,1),'n,n:',Mtx(n,n)
end subroutine mtx_check

subroutine check_unit_mtx_r(Mtx)
real*4, intent(in) :: Mtx(n,n)
real*4 ::  diag, offdiag
diag = 0.0e0; offdiag = 0.0e0
do i=1,n
do j=1,n
  if (i.eq.j) then
     diag = diag + Mtx(i,j)
  else
     offdiag = offdiag + Mtx(i,j)
  endif
enddo
enddo
write(*,"(1X,A,F8.6,A,F8.6)") "aver.sum of diag:",diag/dfloat(n)," aver.sum of offdiag:",offdiag/(dfloat(n*n)-dfloat(n))
end subroutine check_unit_mtx_r


subroutine mtx_check_d(Mtx)
real*8, intent(in) :: Mtx(n,n)
real*8 ::  diag
diag=0.0d0
do i=1,n
  diag = diag + Mtx(i,i)
enddo
write(*,'(1X,A,D8.2,2X,4(A,D8.2,1X))')  & 
 'aver.sum of diag el.:',diag/dfloat(n),'1,1:',Mtx(1,1),'1,n:',Mtx(1,n),'n,1:',Mtx(n,1),'n,n:',Mtx(n,n)
end subroutine mtx_check_d

end module gpu_stuff
