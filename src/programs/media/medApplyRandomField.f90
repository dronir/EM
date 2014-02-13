!****h* XR/programs/medApplyRandomField
! NAME
!   vScatter
!
! DESCRIPTION
!   
!
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   15.11.2007
!******

program medApplyRandomField
!  use iso_binding_c
  use base
  use util
  use geometry
   use particle
  use hemisphere
  use container_grid3d
  use random
  use medium
  use trace
  use sampler
  use rfield
  use brdf

  !$ use omp_lib

  implicit none 

  character (len=FNAME_LENGTH) :: parFilename        = ""
  character (len=FNAME_LENGTH) :: mediumFileName     = "medium.nc"

  character (len = 100)        :: tmpStr             = ""

  integer                      :: gridResHorizontal  = 200
  integer                      :: gridResVertical    = 10

  type(med_medium)                   :: M
  type(med_mediumFile)               :: Mf 
  type(rf_randomField)               :: RF
  type(utl_timer)                    :: timer

  real(fd) :: H, std 
  real     :: sTime, cTime, eTime
  integer  :: i


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! INITIALIZE SIMULATION PARAMETERS
  !!

  if(command_argument_count() < 3) then
     write(*,'("Usage: medApplyRandomField mediumFile H std")')
     stop
  else 
     call get_command_argument(1, mediumFileName)

     call get_command_argument(2, tmpStr)
     read(tmpStr,*) H

     call get_command_argument(3, tmpStr)
     read(tmpStr,*) std
  end if
  
  !! Select the media for the simulation from a file.
  !!
  call med_mediumFileOpen(Mf, mediumFilename)

  !! Generate and save a test random field.
  !! 
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! RUN SIMULATION
  !!

  do i = 1, 1 !Mf%nMedia
     call med_mediumFileRead(M, Mf, i)

     call med_gridAssign(M, gridResHorizontal, gridResVertical)
     call med_gridFit(M)

     call rf_init             (RF, 500, M%width)
     call rf_generateSpectrum (RF, RF_SPECTRUM_FBM, H, std)
     !call rf_generateSpectrum (RF, RF_SPECTRUM_GAUSSIAN, H, std)
     call rf_generateField    (RF, 0)

     write(*,*) M%nParts, M%width

     call med_maskHeight(m, real(rf%field, fs))

     write(*,*) M%nParts

     call med_updateStatistics(M)

     call med_exportRib(M, trim(mediumFilename)//".rib", .true.)

  end do

  write(*,'("Time taken: " ,(F8.3))') cTime - sTime

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! WRITE FILES AND CLEAN UP
  !!

 
  call med_mediumFileClose(Mf)

end program medApplyRandomField

