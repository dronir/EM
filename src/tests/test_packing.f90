!****h* XR/Test_packing
! NAME
!   Test packing.
!
! DESCRIPTION
!   A test for the packing routines.
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   25.08.2007
!******

program test_packing

!  use hdf5

!  use netcdf

  use base

  use particle
  use rfield
  use container_grid3d
  use igdev
  use hemisphere

  implicit none

  integer,  parameter :: n_fields = 1
  integer,  parameter :: n_parts  = 45000
  integer,  parameter :: fres     = 512
  integer,  parameter :: gres_xy  = 10
  integer,  parameter :: gres_z   = 20
  integer,  parameter :: nDrops   = 700

  integer,  parameter :: hmRes    = 512

  real(FD), parameter :: wsize    = 1.2
  real(FD), parameter :: wheight  = 0.5

  type(cnt_grid3d)   :: world
  type(random_field) :: f

  real, dimension(fres, fres) :: rf

  real(FD), dimension(10000) :: rtest

  real(FD) :: a(hmRes,hmRes), b(hmRes,hmRes,3)

  real(FD) :: hMean, hStd, hMin, hMax

  integer :: i,j

  real :: sTime, cTime, eTime

  CHARACTER (LEN=FNAME_LENGTH) :: hfName = "!pHeight.fits"

  real(FD), DIMENSION(:,:), ALLOCATABLE :: pSamples

  CALL CPU_TIME(sTime)

  ! Generate the medium
  !
  call cnt_grid3d_init(world, wsize, wHeight, gres_xy, gres_z, n_parts)
  !call cnt_grid3d_set_distribution_logNormal(world, 0.035_fd, 0.01_fd)
  call cnt_grid3D_set_distribution_invGamma(world, 5.5_fd, 0.035_fd, 0.001_fd, 0.075_fd)
  !call cnt_grid3d_set_distribution_uniform(world, 0.025_fd, 0.10_fd)
  call cnt_grid3d_pack_spheres(world, nDrops)

  write(*,*) SUM(world % parts(:) % r) / world % n_parts

  ! Compute height map of the medium
  !
  b(:,:,1) = cnt_grid3D_computeHeightMap(world, hmRes)

  ! Compute medium statistics
  !
  call cnt_grid3D_computeStats(world, hmRes, hMean, hStd, hMin, hMax)
  write(*,*) hMean, hStd, hMin, hMax

  write(*,*) "Porosity:", cnt_grid3d_compute_porosity(world)
  
  ! Copmute random field
  !
  call dfftw_init_threads()
  call dfftw_plan_with_nthreads(2)
  call initialize_field(f, fres, wsize)
  call generate_spectrum_m(f,0.30_fd)
  call generate_field(f,0)
  call dfftw_cleanup_threads()

  rf       = real(f%field) * 0.2_fd
  b(:,:,3) = rf

  CALL CPU_TIME(cTime)
  write(*,*) "time ", cTime - sTime

  ! Mask the medium with the random field
  !
  !call cnt_grid3D_mask_height(world, rf)


  ! Compute medium statistics
  !
  call cnt_grid3D_computeStats(world, hmRes, hMean, hStd, hMin, hMax)
  write(*,*) hMean, hStd, hMin, hMax


  ! Compute height map of the new medium
  !
  b(:,:,2) = cnt_grid3D_computeHeightMap(world, hmRes)

  CALL CPU_TIME(cTime)
  write(*,*) "time ", cTime - sTime

  call cnt_grid3d_export_rib(world)
  call utl_writeFits3D(b, hfName)

  !write(*,'(9F6.2)') pSamples

end program test_packing

