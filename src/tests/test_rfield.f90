!****h* XR/Test_rfield
! NAME
!   Test rfield
!
! DESCRIPTION
!   Basic test.
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   20.09.2007
!******

program test_rfield
!$ use omp_lib
   use base
   use rfield

  implicit none

  integer,  parameter :: n_fields = 1
  integer,  parameter :: fres     = 500
  real(FD), parameter :: width    = 10.0

  real(fd)            :: p1       = 0.5
  real(fd)            :: p2       = 1.0

  character(len=100)  :: tempStr

  type(rf_randomField)  :: rf
  real :: sTime, cTime, eTime

  !$ call omp_set_num_threads(4)

  call dfftw_init_threads()
  call dfftw_plan_with_nthreads(4)

  if(command_argument_count() == 2) then
     call get_command_argument(1, tempStr)
     read(tempStr,*) p1
     call get_command_argument(2, tempStr)
     read(tempStr,*) p2
  end if

  call cpu_time(sTime)

  call rf_init(rf, fres, width)

  call cpu_time(sTime)
  write(*,'("Generating spectrum:")')
  call rf_generateSpectrum(rf, RF_SPECTRUM_FBM, p1, p2)
  call cpu_time(cTime)
  write(*,'("  Time: ", F6.3)') cTime - sTime

  call cpu_time(sTime)
  write(*,'("Generating field:")')
  call rf_generateField(rf,1)
  call cpu_time(cTime)
  write(*,'("  Time: ", F6.3)') cTime - sTime

  call dfftw_cleanup_threads()

  call rf_computeStatistics(rf)

  call rf_writeFits(rf, "!testField.fits")

  open(1, file="rf.dat")
  write(1,*) real(rf%field,fd)
  close(1)

end program test_rfield

