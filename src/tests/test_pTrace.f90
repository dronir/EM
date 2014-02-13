!****h* XR/Test_ptrace
! NAME
!   Test trace
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

program test_ptrace
  use base
  use geometry
  use particle
  use hemisphere
  use medium
  use trace

  implicit none

  integer,  parameter :: nParts   = 1
  integer,  parameter :: gResH    = 2
  integer,  parameter :: gResV    = 1
  integer,  parameter :: hres_t   = 500

  real(FD), parameter :: wsize    = 6.0
  real(FD), parameter :: wheight  = 6.0

  integer, parameter  :: nSamples = 10

  type(gth_hemisphere) :: H
  type(med_medium)     :: M

  real :: sTime, cTime, eTime

  character (len=FNAME_LENGTH) :: fname_h = "pTrace_h.eps"

  call med_medium_init(M, wsize, wHeight, nParts, "medium")
  call gth_hemisphere_init(H, hres_t)

  M%parts(1)%P  = (/1.0_fd, 0.0_fd, 1.0_fd/)
  M%parts(1)%r  = 0.5_fd
  M%parts(1)%rr = M%parts(1)%r ** 2

  call med_gridAssign(M,gResH, gResV)
  call med_gridFit(M)

  CALL CPU_TIME(sTime)
  call sample_medium  (M, H)
  CALL CPU_TIME(cTime)
  WRITE(*,'(F6.3)') cTime - sTime

  CALL gth_hemisphere_plot_data(H, fname_h)
  call gth_hemisphere_save          ( H,        "pTest.nc" )

contains
  subroutine sample_medium(m, h)
    type(med_medium)     :: m
    type(gth_hemisphere) :: h

    type(intersection_geometry) :: isect

    type(ray) :: r
    integer :: i, j, k
    real(FD) :: t, a
    logical  :: iFound

!!$    r % rayID = r % rayID + 1
!!$    r % P     = (/1e-6_fd, 1e-6_fd, 1e-6_fd/)
!!$    r % D     = (/-2.0_fd,  0.0_fd,  1.0_fd/)
!!$
!!$    call vec_normalize(r%D)
!!$
!!$    iFound = trc_traceNearest(M%grid,r,iSect)
!!$    write(*,*) iFound 
!!$    stop

    a = 1.0_fd / real(nSamples, fd)

     do i = 1, h % res_theta
        do j = 1, h % rows(i) % res_phi
           do k = 1, nSamples
 
              r % rayID = r % rayID + 1
              r % P     = (/1e-6_fd, 1e-6_fd, 0.1_fd/)
              r % D     = gth_cellRandomSampleCar(h, i, j)

              if(r%D(3) < 1e-2) then
                 r%D(3) = 1e-2
              end if

              if(trc_traceNearest(M%grid,r,iSect)) then
!!$              if(.NOT.trc_traceOcclusion(c,r)) then
                 t = abs(iSect%N(3) * a)
              else
                 t = 0.0
              end if

              !$omp atomic
              h % rows(i) % cells(j) = h % rows(i) % cells(j) + t

           end do
        end do
     end do

  end subroutine sample_medium

end program test_ptrace

