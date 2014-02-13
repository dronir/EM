!****h* XR/lib/simulation/simulation
! NAME
!   simulation
!
! DESCRIPTION
!
! NOTES
!
! AUTHOR
!   Hannu Parviainen
!
! USES
!   Base
!   Particle
!
!
! CREATION DATE
!   19.12.2007
!******
module simulation
implicit none

  !****d* container_grid3D/pMin
  ! NAME
  !   pMin
  !
  ! DESCRIPTION
  !   Minimum number of particles allocated for a vector.
  !
  !******
  integer, parameter :: pMin = 5

  !****s* container_grid3D/g_cell
  ! NAME
  !   g_cell
  !
  ! DESCRIPTION
  !
  !******
  type :: cnt_gridCell
     integer :: n_parts, lsize
     integer, dimension(:), pointer :: plist
  end type cnt_gridCell

  !****s* container_grid3D/cnt_grid3D
  ! NAME
  !   cnt_grid3D
  !
  ! DESCRIPTION
  !
  !******
  type :: cnt_grid3D
     real(fd) :: rw, wr, rh, hr
     integer  :: resHor, resVer

     !****d* container_grid3D/maxHeight_g
     ! NAME
     !   maxHeight_g
     !
     ! DESCRIPTION
     !   Maximum global gridded height of the container.  
     !
     !******
     integer :: maxHeight_g 

     !****d* container_grid3D/maxHeight_l
     ! NAME
     !   maxHeight_l
     !
     ! DESCRIPTION
     !   Maximum gridded height of the container per (x,y)-cell.  
     !
     !******
     integer, dimension(:,:), pointer :: maxHeight_l

     !****d* container_grid3D/m
     ! NAME
     !   m
     !
     ! DESCRIPTION
     !   Pointer to the medium.
     !
     !******
     type(med_medium), pointer :: m

     !****d* container_grid3D/cells
     ! NAME
     !   cells
     !
     ! DESCRIPTION
     !   Three-dimensional array of cnt_gridCell structures holding the actual
     !   grid-structure.
     !
     !******
     type(cnt_gridCell), dimension(:,:,:), pointer :: cells
  end type cnt_grid3D

end module simulation
