!****h* XR/Test_xr
! NAME
!   Test xr
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

program test_xr
!$  use omp_lib

  use base
  use elements
  use material
  use geometry
  use particle
  use hemisphere
  use quartersphere
  use container_grid3d
  use trace
  use mtprng

  use random

  implicit none

  integer,  parameter :: n_fields = 1
  integer,  parameter :: n_parts  = 5
  integer,  parameter :: fres     = 512
  integer,  parameter :: gres_xy  = 200
  integer,  parameter :: gres_z   = 20
  integer,  parameter :: hres_t   = 180
  integer,  parameter :: nDrops   = 100

  real(FD), parameter :: wsize    = 40.0
  real(FD), parameter :: wheight  = 6.0

  integer, parameter  :: nSamples = 1

  type(gth_hemisphere)        :: hemis
  type(cnt_grid3d)            :: world
  type(ray)                   :: r
  type(spc_spectrum)          :: spc

  type(elm_element)           :: e(1)
  real(FD)                    :: eFrac(1) = (/1.0_fd/)
  type(mat_material), TARGET  :: mat

  real :: sTime, cTime, eTime

  integer :: i

  character (len=FNAME_LENGTH) :: fname_h = "hemisphere_xr.eps"

  !character (len=FNAME_LENGTH) :: fname_FE = "../data/FE.dat"
  character (len=2)            :: fname_FE = "fe"
  character (len=FNAME_LENGTH) :: fname_SS = "../data/Solsim.dat"
  character (len=256)          :: matName  = "Iron"

  real(FD) :: TT
  !$ call omp_set_num_threads(2)

  call cnt_grid3d_init(world, wsize, wHeight, gres_xy, gres_z, n_parts)
  call gth_hemisphere_init(hemis, hres_t)

  call elm_element_init(e(1), fname_FE)
  call mat_material_init(mat, matName, e, eFrac)

  call spc_spectrum_read_discrete(spc, fname_SS)

  call cnt_grid3d_set_distribution_constant(world, 0.5_fd)
  world % parts(1) % P = 0.0_fd
  world % parts(2) % P = (/-1.1, 0.0, 0.0/)
  world % parts(2) % P = (/0.0, 0.0, 1.1/)
  world % parts(5) % P = (/0.0, 0.0, 4.2/)
  world % parts(4) % P = (/0.0, 0.0, 3.2/)
  world % parts(3) % P = (/0.0, 0.0, 2.2/)

  !call cnt_grid3d_pack_spheres(world, nDrops)
  call cnt_grid3D_optimize_grid(world)

  mat % muTable = 1.0_fd

  DO i=1, world%n_parts
     world%parts(i)%mat => mat
  END DO

  CALL CPU_TIME(sTime)
  call sample_medium(world, hemis)
  CALL CPU_TIME(cTime)

  WRITE(*,'(F6.3)') cTime - sTime

  CALL gth_hemisphere_plot_data(hemis, fname_h)


contains
  subroutine sample_medium(c, h)
    type(cnt_grid3D)     :: c
    type(gth_hemisphere) :: h
    
    type(ray) :: r
   
    integer  :: i, j, k
    real(FD) :: t, a

    a = 1.0_fd / real(nSamples)

    r%energy = 5100.0_fd

    r % P = (/0.0, 0.0, 0.1/)
    r % D = (/0.0, 0.0, 1.0/)
    r % rayID = 100

    call trc_traceRayToOpticalDepth(c, r, 1.5_fd)

    r % P = 0.0_fd
    r % P(1) = -0.25_fd

    do i = 1, h % res_theta
        do j = 1, h % rows(i) % res_phi
           do k = 1, nSamples
              r % rayID = r % rayID + 1
              r % D     = gth_cellRandomSampleCar(h, i, j)

              t = trc_traceOpticalDepth(c,r)

              !write(*,*) t

              !$omp atomic
              h % rows(i) % cells(j) = h % rows(i) % cells(j) + exp(-t)
           end do
        end do
     end do

  end subroutine sample_medium

end program test_xr

