!****h* XR/Test_trace
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

program test_trace
!$  use omp_lib

  use base
  use geometry
  use particle
  use hemisphere
  use quartersphere
  use container_grid3d
  use trace

  implicit none

  
  integer,  parameter :: n_fields = 1
  integer,  parameter :: n_parts  = 6000
  integer,  parameter :: fres     = 512
  integer,  parameter :: gres_xy  = 200
  integer,  parameter :: gres_z   = 20
  integer,  parameter :: hres_t   = 180
  integer,  parameter :: nDrops   = 100

  real(FD), parameter :: wsize    = 40.0
  real(FD), parameter :: wheight  = 6.0

  integer, parameter  :: nPntSmpA = 2

  integer, parameter  :: nSamples = 1
  integer, parameter  :: nPosSamples = nPntSmpA**2

  type(intersection_geometry) :: isect
  type(gth_hemisphere) :: hemis
  type(gth_quartersphere) :: quart
  type(cnt_grid3d)     :: world
  type(ray)            :: r

  integer :: i,j, k, l, tCeDir(2)
  real(FD) :: tSpDir(2), tCaDir(3), hError

  real :: sTime, cTime, eTime

  real(FD), DIMENSION(:,:), ALLOCATABLE :: pSamples

  character (len=FNAME_LENGTH) :: fname_h1 = "hemisphere_1.eps"
  character (len=FNAME_LENGTH) :: fname_h2 = "hemisphere_2.eps"
  character (len=FNAME_LENGTH) :: fname_h3 = "quartersphere.eps"

!$  call omp_set_num_threads(2)

  call cnt_grid3d_init(world, wsize, wHeight, gres_xy, gres_z, n_parts)
  call cnt_grid3d_set_distribution_invGamma(world, 14.05_fd, 1.1_fd, 0.075_fd, 0.250_fd)

  call gth_hemisphere_init(hemis, hres_t)
  call gth_quartersphere_init(quart, hres_t)

  call cnt_grid3d_pack_spheres(world, nDrops)
  
  CALL CPU_TIME(sTime)

  ALLOCATE(pSamples(2, nPosSamples))
  pSamples = utl_griddedSamples2D(nPntSmpA)

  !$omp parallel do
  do l = 1, nPosSamples
     call sample_medium  (world, hemis, pSamples(:,l))
     call sample_medium_q(world, quart, pSamples(:,l))
  end do
  !$omp end parallel do

  !$ write(*,*) omp_get_thread_num()

  CALL CPU_TIME(cTime)
  WRITE(*,'(F6.3)') cTime - sTime

  call cnt_grid3d_export_rib(world)
  CALL gth_hemisphere_plot_data(hemis, fname_h2)
  CALL gth_hemisphere_plot_grid(hemis)

  CALL gth_quartersphere_plot_data(quart, fname_h3)

contains
  subroutine sample_medium(c, h, pSample)
    type(cnt_grid3D)     :: c
    type(gth_hemisphere) :: h
    real(FD), DIMENSION(2) :: pSample

    type(ray) :: r

    integer :: i, j, k

    real(FD) :: t, a

    a = 1.0_fd / real(nSamples) / nPosSamples

     CALL RANDOM_NUMBER(r%P)
     r % P = (r%P - 0.5) * c % width
     r % P(1:2) = pSample
     r % P(3) = cnt_grid3D_getHeight_xy(c, r%P(1:2)) + 1e-4

     do i = 1, h % res_theta
        do j = 1, h % rows(i) % res_phi
           do k = 1, nSamples
              r % rayID = r % rayID + 1
              r % D     = gth_cellRandomSampleCar(h, i, j)
             !if(.NOT.trc_traceNearest(c,r,iSect)) then
!!$              if(.NOT.trc_traceOcclusion(c,r)) then
!!$                 t = a
!!$              else
!!$                 t = 0.0
!!$              end if
              !t         = trc_trace(c, r) * a

              t = trc_traceOpticalDepth(c,r)

              if(t/=0.0) t =  EXP(-t)

              !$omp atomic
              h % rows(i) % cells(j) = h % rows(i) % cells(j) + t

           end do
        end do
     end do

  end subroutine sample_medium

  subroutine sample_medium_q(c, h, pSample)
    type(cnt_grid3D)     :: c
    type(gth_quartersphere) :: h
    real(FD), DIMENSION(2) :: pSample

    type(ray) :: r

    integer :: i, j, k

    real(FD) :: t, a

    a = 1.0_fd / real(nSamples) / nPosSamples

     CALL RANDOM_NUMBER(r%P)
     r % P = (r%P - 0.5) * c % width
     r % P(1:2) = pSample
     r % P(3) = cnt_grid3D_getHeight_xy(c, r%P(1:2)) + 1e-4

     do i = 1, h % res_theta
        do j = 1, h % rows(i) % res_phi
           do k = 1, nSamples
              r % rayID = r % rayID + 1
              r % D     = gth_qrts_cellRandomSampleCar(h, i, j)
              !if(.NOT.trc_traceNearest(c,r,iSect)) then
!!$              if(.NOT.trc_traceOcclusion(c,r)) then
!!$                 t = a
!!$              else
!!$                 t = 0.0
!!$              end if

              t = trc_traceOpticalDepth(c,r)

              !t         = trc_trace(c, r) * a

              !$omp atomic
              h % rows(i) % cells(j) % I(1) = h % rows(i) % cells(j) % I(1) + t
              
              !$omp atomic
              h % rows(i) % cells(j) % weight = h % rows(i) % cells(j) % weight + 1.0_fd
              
           end do
        end do
     end do

  end subroutine sample_medium_q

end program test_trace

