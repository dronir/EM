!****h* XR/Testi
! NAME
!   Testi
!
! DESCRIPTION
!   Basic test.
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   25.08.2007
!******

program testi

!  use hdf5

!  use netcdf

  use base

  use particle
  use rfield
!  use container_grid2d
  use container_grid3d
  use igdev
  use hemisphere

  implicit none

  REAL TIMEARRAY (2)
  REAL ELAPSE
  
  integer,  parameter :: n_fields = 1
  integer,  parameter :: n_parts  = 400
  integer,  parameter :: fres     = 512
  integer,  parameter :: gres_xy  = 100
  integer,  parameter :: gres_z   = 20
  integer,  parameter :: nDrops   = 500

  real(FD), parameter :: wsize    = 20.0
  real(FD), parameter :: wheight  = 6.0

!  type(cnt_grid2d)   :: grid2
  type(cnt_grid3d)   :: grid3
  type(random_field) :: f

  type(gth_hemisphere) :: hp

  real, dimension(fres, fres) :: rf

  real(FD), dimension(10000) :: rtest

  integer :: i,j, tCeDir(2)
  real(FD) :: tSpDir(2), tCaDir(3)

  real :: sTime, cTime, eTime

  CHARACTER (LEN=FNAME_LENGTH) :: hfname = "hemisphere.eps"

  CALL CPU_TIME(sTime)

  call cnt_grid3d_init(grid3, wsize, wHeight, gres_xy, gres_z, n_parts)
  call cnt_grid3d_pack_spheres(grid3, nDrops)

   !call cnt_grid2d_init(grid2, wsize, gres_xy, n_parts)
  !call cnt_grid2d_pack_spheres(grid2)

  write(*,*) "Porosity:", cnt_grid3d_compute_porosity(grid3)
  
  call dfftw_init_threads()
  call dfftw_plan_with_nthreads(2)

  call initialize_field(f, fres, wsize)
  call generate_spectrum(f,0.70_fd)
  !call generate_spectrum_m(f,0.60_fd)

  do i = 1, n_fields
     call generate_field(f,0)
  end do

  call dfftw_cleanup_threads()

  rf = real(f%field) * 0.2

  CALL CPU_TIME(cTime)
  write(*,*) "time ", sTime - cTime

  !call cnt_grid2d_mask_height(grid, rf)

  CALL CPU_TIME(cTime)
  write(*,*) "time ", sTime - cTime


  open(1, FILE='rfdata_f.txt', status='REPLACE')

  do i = 0, f % res
     write(1,*) real(f % field(i,:))
  end do

  close(1)

  call dst_generate_invGamma(rtest, 15.5_fd, 1.042_fd, 0.075_fd, 0.25_fd)
  open(1,FILE='gammq.txt', status='REPLACE')

  do i = 1, size(rtest)
     write(1,*) rtest(i)
  end do

  close(1)

  call cnt_grid3d_export_rib(grid3)
  !call cnt_grid2d_export_rib(grid2)

  call gth_hemisphere_init(hp, 10)

  tCeDir = gth_sphDirToCell(hp, (/ 0.0_fd * PI , 1.0_fd * PI /))
  write(*,*) tCeDir

  tCeDir = gth_carDirToCell(hp, (/ 0.0_fd, 0.7_fd, 0.7_fd /))
  write(*,*) tCeDir

  CALL gth_hemisphere_plot_data(hp, hfname)

end program testi

