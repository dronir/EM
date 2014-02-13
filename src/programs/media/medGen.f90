!!--- BEGIN GPL --- 
!!
!! medgen -  a program to generate particulate media packings.
!! Copyright (C) 2009 Hannu Parviainen
!!
!! This program is part of UHOEM.
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! Contributor(s): Hannu Parviainen
!!
!!--- END GPL ---


!> A program to generate particulate media packings.
!! 
!! \use     medgen parameter_file
!!
!! \author  Hannu Parviainen
!! \date    4.12.2007
!!
!! Parameter file options:
!!
!!  General properties
!!   mediumFilename
!!   mediumBasename
!!   nParticles
!!   nMedia
!!
!!  Accelerating grid
!!   resHorizontal
!!   resVertical
!!
!!  Medium dimensions
!!   mediumWidth
!!   mediumHeight
!!
!!  Size distribution
!!   sDistribution
!!   sDistParams
!!
!!  Dropping
!!   nnDrops
!!   nDrops
!!   dropSeed
!!
program medgen
  !$ use omp_lib
     use base
     use medium

  implicit none

  !!--- INPUT VARIABLES ---
  !!
  character (len=FNAME_LENGTH) :: parFileName

  character (len=FNAME_LENGTH) :: mediumFileName  = "medium.nc" !< Name of the medium file to be written.
  character (len=256)          :: mediumBaseName  = "medium"    !< Base name for the media.

  real(fd)                     :: mediumWidth   = 1.0_fd        !< Width of the media.
  real(fd)                     :: mediumHeight  = 0.05_fd       !< Expected maximum height of the media.

  integer                      :: nParticles    = 25000         !< Number of particles.
  integer                      :: resHorizontal = 400           !< Horizontal resolution of the accelerating grid.
  integer                      :: resVertical   = 10            !< Vertical resolution of the accelerating grid.
  integer                      :: nMedia        = 1             !< Number of media to create.

  character (len=100)          :: sDistribution                 !< Name of the particle size distribution.
  character (len=100)          :: sDistParams   = "0.02"        !< Particle size distribution parameters.

  integer                      :: nnDrops       = 1             !< Size of the drop table.
  character (len=200)          :: nDrops        = "50"          !< Number of drops.

  logical                      :: exportRib     = .false.       !< Should we generate a rib of the media.

  !!--- SIMULATION VARIABLES ---
  !!
  real(fd)                            :: sDistParArr(4)         !< Size distribution parameters.
  integer, dimension(:), allocatable  :: nd                     !< Number of drops. Read from the nDrops array.
  character (len=256)                 :: mediumName             !< Name of each generated medium

  type(med_medium)                    :: M                      !< The medium.
  type(med_mediumFile)                :: Mf                     !< The medium file.

  type(utl_timer)                     :: simCpuTimer
  real                                :: sTime, cTime, eTime
 
  integer :: dropSeed = 1
  integer :: distSeed = 1
 
  namelist /params/ nParticles, nDrops,  resHorizontal, resVertical, mediumWidth, mediumHeight, &
       & mediumFileName, mediumBaseName, sDistribution, sDistParams, nMedia, &
       & nnDrops, dropSeed, exportRib

  !!
  !!--- INITIALIZE MEDGEN ---
  !!
  !!
  if(command_argument_count() == 0) then
     write(*,'("No parameter file given, using ""medium.in"" ")')
     parFileName = "medium.in"
  else if(command_argument_count() == 1) then
     call get_command_argument(1, parFileName)
     write(*,'(("Using parameter file "),(A))') parFileName
  else
     write(*,'("Usage: medGen paramFile")')
     stop
  end if
  !!
  !!- Read simulation parameters from the input file
  !!
  open(1,file=parFileName,status="old", form="formatted")
  read(1,NML=params)
  close(1)
  !!
  !!- Trim input data
  !!
  mediumFileName = trim(mediumFileName)
  mediumBaseName = trim(mediumBaseName)
  sDistribution  = trim(sDistribution)
  sDistParams    = trim(sDistParams)
  !!
  !!- Allocate and read drop table
  !!
  allocate(nd(nnDrops))
  read(nDrops,*) nd
  !!
  !!- Write medium name
  !!
  write(mediumName,'((A),("_"),(I3.3),("_"),(I3.3),("_"),(I3.3))') &
       & trim(adjustl(mediumBaseName)), nd(1), distSeed, dropSeed
  !!
  !!- Initialise medium
  !!
  call med_medium_init(M, mediumWidth, mediumHeight, nParticles, mediumName)
  !!
  !!
  !!--- SET SIZE DISTRIBUTION ---
  !!
  !!
  call rnd_init(1, 4, .false.)
  !!
  select case(sDistribution)
     case('constant')
        write(*,'("Using constant particle size distribution.")')
        read(sDistParams,*)  sDistParArr(1)
        call med_setDistribution_constant(M,  sDistParArr(1))
     case('uniform')
        write(*,'("Using uniform particle size distribution.")')
        read(sDistParams,*)  sDistParArr(1:2)
        call med_setDistribution_uniform(M, sDistParArr(1), sDistParArr(2), distSeed)
     case default
        write(*,'("Error: Unsupported particle size distribution. Exiting.")')
        stop
  end select
  !!
  !!- initialise timer
  !!
  call utl_timer_init(simCpuTimer, 5.0_fd, nParticles * nMedia)
  !!
  !!--- START PACKING ---
  !!
  !call cpu_time(sTime) 
  do dropSeed = 1, nMedia
     call med_gridAssign(M, resHorizontal, resVertical)
     call setName()     

     call med_pack(M, nd(1), dropSeed, nd, simCpuTimer)

     call med_updateStatistics(M)
     call med_gridFit(M)

     if(exportRib) call med_exportRib(M, cut=.true.)

     call med_medium_save(M, mediumFileName)
     call med_gridRemove(M)
  end do
  !call cpu_time(cTime)
  !!
  !!--- PACKING FINISHED
  !!
  !print *, "Total CPU time: " // utl_secs_to_hms(cTime - sTime)

contains
  subroutine setName()
     write(mediumName,'((A),("_"),(I3.3),("_"),(I3.3),("_"),(I3.3))') &
          & trim(adjustl(mediumBaseName)), nd(1), distSeed, dropSeed
     call med_setName(M, mediumName)
  end subroutine setName

end program medgen

