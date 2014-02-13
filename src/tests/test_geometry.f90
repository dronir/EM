!****h* XR/Test_geometry
! NAME
!   Test geometry
!
! DESCRIPTION
!   Basic geometry tests.
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   19.11.2007
!******

program test_xr
  use base
  use geometry
  use hemisphere
  use quartersphere

  implicit none

  integer, parameter      :: hres_t   = 10
  integer, parameter      :: nSamples = 500000

  type(gth_hemisphere)    :: hs
  type(gth_quartersphere) :: qs
  type(ray)               :: r
  real                    :: sTime, cTime, eTime
  real(FD)                :: v
  integer                 :: i

  character(len=FNAME_LENGTH)       :: fname_hs = "test_hemisphere_hs.eps"
  character(len=FNAME_LENGTH)       :: fname_qs = "test_hemisphere_qs.eps"

  call gth_hemisphere_init(hs, hres_t)

  v = 1.0_fd / (real(nSamples, KIND=FD) * hs%dA_mean)

  CALL CPU_TIME(sTime)
  DO i=1, nSamples
     !r % D = vec_cart_random_hemispherical_uniform()
     !call gth_addDataCar(hs, r%D, v)
  END DO
  CALL CPU_TIME(cTime)
  WRITE(*,'(F6.3)') cTime - sTime

  CALL gth_distributeValue(hs, 2.0_fd)

  !CALL gth_hemisphere_normalize(hs)
  WRITE(*,*) gth_hemisphere_integrate(hs)

  CALL gth_hemisphere_plot_data(hs, fname_hs)

end program test_xr

