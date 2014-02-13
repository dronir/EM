!****h* XR/programs/medRender
! NAME
!   medRender
!
! DESCRIPTION
!   
!
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   28.01.2008
!******

program medRender
!  use iso_binding_c
  use base
  use util
  use geometry
  use particle
!  use container_grid3d
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
  character (len=100)          :: sDistribution      = "constant"
  character (len=100)          :: brdfType           = "shadowing"

  character (len=100)          :: cBuf               = ""

  character (len=300)          :: thetaInTable       = ""

  real(fd)                     :: thetaIn            = 45.0_fd
  real(fd)                     :: rho                = 0.50_fd
  real(fd)                     :: rhoAllowedError    = 1e3_fd

  integer                      :: nOrders            = 1
  integer, dimension(:), allocatable :: nSamples  
  integer                      :: sampleSeed         = 0

  integer                      :: resolution         = 200
  integer                      :: gridResHorizontal  = 200
  integer                      :: gridResVertical    = 10

  integer                      :: thetaInTableN      = 0
  logical                      :: useThetaInTable    = .false.

  type(med_medium)             :: M
  type(med_mediumFile)         :: Mf
  type(utl_timer)              :: timer

  real(fd), dimension(:), allocatable :: thetaInTableV

  real    :: sTime, cTime, eTime
  integer :: i

  namelist /params/ &
       & resolution, gridResHorizontal, gridResVertical, &
       & mediumFilename, &
       & sDistribution, rho, rhoAllowedError, thetaIn, sampleSeed, &
       & thetaInTable, useThetaInTable, thetaInTableN, brdfType


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! INITIALIZE SIMULATION PARAMETERS
  !!

  if(command_argument_count() < 5) then
     write(*,'("Usage: medRender mediumFile resolution nOrders nSamples[]")')
     stop
  else
     call get_command_argument(1, mediumFileName)
     call get_command_argument(2, cBuf)
     read(cBuf,*) resolution
     call get_command_argument(3, cBuf)
     read(cBuf,*) nOrders
     
     allocate(nSamples(nOrders))

     do i =1, nOrders
        call get_command_argument(3+i, cBuf)
        read(cBuf,*) nSamples(i)
     end do
  end if

  if(parFilename /= "") then
     open(1,file=parFilename,status="old", form="formatted")
     read(1,NML=params)
     close(1)
  end if

  if(useThetaInTable) then
     allocate(thetaInTableV(thetaInTableN))
     read(thetaInTable,*) thetaInTableV
  end if

  call med_mediumFileOpen(Mf, mediumFilename)
  call med_mediumFileSelectMedia(Mf, rho, rhoAllowedError, sDistribution)

  !call utl_timer_init(timer, 5.0_fd, Mf%nSelectedMedia * H%nCells) !* nSamplesPerMedium)
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! RUN SIMULATION
  !!

  call utl_message("BRDF type used: " // brdfType)

  do i = 1, Mf%nSelectedMedia
     call med_mediumFileRead(M, Mf)

     call med_gridAssign(M, gridResHorizontal, gridResVertical)
     call med_gridFit(M)

     call med_updateStatistics(M)

     call utl_message("Beginning trace.")
     call cpu_time(sTime) 
     call renderSurface("!testOut.fits", resolution, nOrders, nSamples, 55.0_fd, M)
     call cpu_time(cTime)

  end do

  write(*,'("Time taken: " ,(F6.3))') cTime - sTime

  call med_mediumFileClose(Mf)

contains

  subroutine renderSurface(fitsFileName, resolution, nOrders, nSamples, thetaI, M)

    character(len=*) :: fitsFileName
    integer          :: resolution, nOrders
    integer          :: nSamples(nOrders)
    real(fd)         :: thetaI, thtIn
    type(med_medium) :: M
    real(fd)         :: Il(nOrders)

    type(ray) :: rC
    type(ray) :: rS
    type(intersection_geometry) :: iSect

    real(fd), dimension(2,resolution**2) :: samples
    real(fd), dimension(3)               :: pSurface

    real(fd), dimension(resolution,resolution,nOrders+1) :: image

    integer :: x, y, o

    logical :: pFound, pLit
 
    Il = 0.0_fd

    image = 0.0_fd

    thtIn   = thetaI * (PI/180.0_fd)

    call smpl_griddedSamples2D(samples, resolution**2, 0)

    samples = (samples * M%width - M%hWidth) * 0.99_fd

    !$omp parallel default(none) &
    !$omp shared(M, samples, resolution, thtIn, image, nOrders, nSamples) &
    !$omp private(y,x,Il,o,rS,rC,iSect,pFound)

    call ray_init(rC, RAY_TYPE_CAMERA)
    call ray_init(rS, RAY_TYPE_SHADOW)

    rS%D = [sin(thtIn), 0.0_fd, cos(thtIn)]

    !$omp do 
    do y = 1,resolution
       do x = 1, resolution

          Il = 0.0_fd
          o  = 1

          !!write(*,*) omp_get_thread_num(), x, y

          rC%P(1:2) = samples(:, x + (y-1) * resolution)
          rC%P(3)   = M%grid%height - TRACE_EPS
          rC%D      = [0.0_fd, 0.0_fd, -1.0_fd]
          
          rC%rayID  = rC%rayID + RAY_ID_INCR

          pFound    = trc_traceNearest(M%grid, rC, iSect)

          if(pFound) then

             !write(*,*) "aa", iSect%P1
             call trc_gatherRadiance(M%grid, rC%D, rS%D, iSect%P1 + TRACE_EPS * iSect%N, iSect%N, 1.0_fd, &
                  & nSamples, nOrders, o, Il, brdf_Lambert, brdf_Lambert_f)

             image(x,y,1:nOrders) = Il
             image(x,y,nOrders+1) = sum(image(x,y,1:nOrders))

          end if
 
       end do
       if(mod(y,20) == 0) write(*,*) y, " of ", resolution
    end do
    !$omp end do
    !$omp end parallel

    print*,"Done"
    call utl_writeFits3D(image, fitsFileName)

  end subroutine renderSurface

end program medRender

