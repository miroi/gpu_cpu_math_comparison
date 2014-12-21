module fortran_timing
! http://software.intel.com/sites/products/documentation/hpc/compilerpro/en-us/fortran/lin/compiler_f/lref_for/source_files/rfsystmc.htm
! If the type is INTEGER(2), count_rate is 1000. If the type is INTEGER(4), count_rate is 10000. If the type is INTEGER(8), count_rate is 1000000.
! but gfortran 4.8.0 ([trunk revision 185723]) has COUNT_RATE=1000000000 for INTEGER(8)...
implicit none
integer(kind=8) :: c1,c2,cr,cm ! accuracy in nanoseconds
!integer(kind=4) :: c1,c2,cr,cm ! accuracy in miliseconds
! default control values - these differ from compiler to compiler !
integer, parameter :: count_rate_kind4=1000, count_rate_kind8=1000000
contains

subroutine print_time_info
call system_clock( count_rate=cr, count_max=cm )
write(*,'(1x,a,i8,a,i20)') "--- Fortran system_clock timing routine, COUNT_RATE=",cr," COUNT_MAX=",cm
!if (cr.ne.count_rate_kind8) then
!  write(*,'(2x,a,i7)') 'wrong system value of system COUNT_RATE parameter, must be (kind=8) ',count_rate_kind8
  !write(*,'(2x,a,i7)') 'wrong system value of system COUNT_RATE parameter, must be (kind=4) ',count_rate_kind4
!endif
end subroutine print_time_info

subroutine  print_time(message)
! print spent time in routine between c1 and c2 
! http://gcc.gnu.org/onlinedocs/gfortran/SYSTEM_005fCLOCK.html
! print spent time in routine ...
character(*), intent(in) :: message
write(*,'(2x,a,a,a,f13.8,a)') 'Routine <',message,'> spent runtime:', (dfloat(c2-c1) / dfloat(cr)), ' seconds.'
end subroutine  print_time

end module fortran_timing
