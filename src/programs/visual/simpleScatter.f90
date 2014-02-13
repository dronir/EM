!****h* XR/programs/simpleScatter
! NAME
!   simpleScatter
!
! DESCRIPTION
!   
!
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   20.02.2008
!******

module scatterStats
  use base 

  type :: scStat
     integer  :: raysMissed = 0

     real(fd) :: timeMediaLoad = 0.0_fd
     real(fd) :: timeFieldGen  = 0.0_fd
     real(fd) :: timeFieldApp  = 0.0_fd
     real(fd) :: timeSample    = 0.0_fd

  end type scStat
end module scatterStats

program simpleScatter
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

  use scatterStats

  !$ use omp_lib

  implicit none 

  character (len=FNAME_LENGTH) :: parFilename        = ""
  character (len=FNAME_LENGTH) :: mediumFileName     = "medium.nc"
  character (len=FNAME_LENGTH) :: outFileName        = "scatterData.sct"
  character (len=100)          :: brdfType           = "LommelSeeliger"

  character (len=10)           :: tmpStr

  integer                      :: gridResHorizontal  = 200
  integer                      :: gridResVertical    = 10

  integer                      :: sampleSeed         = 0


  integer                      :: nOrders            = 1
  character(len=100)           :: nSamplesPerOrder   = "1"
  integer, allocatable         :: nSamplesPerOrderTable(:)


  !! Surface roughness variables
  !!
  !!

  integer                      :: rf_nFieldsPerMed  = 1
  character(len=100)           :: rf_spectrumType   = "Gaussian"

  real(fd)                     :: thetaIMax         = 89.0_fd
  real(fd)                     :: thetaIMin         = 0.0_fd
  integer                      :: thetaIN           = 90

  real(fd)                     :: thetaEMax         = 10.0_fd
  real(fd)                     :: thetaEMin         = 0.0_fd
  integer                      :: thetaEN           = 11

  real(fd)                     :: roughPMax         = 1.0_fd
  real(fd)                     :: roughPMin         = 1.0_fd
  integer                      :: roughPN           = 1

  real(fd)                     :: srfStdMax         = 0.0_fd
  real(fd)                     :: srfStdMin         = 0.0_fd
  integer                      :: srfStdN           = 1


  real(fd)                     :: medDensity        = 0.2

  type(med_medium)             :: M0, M
  type(med_mediumFile)         :: Mf
  type(rf_randomField)         :: RF
  type(utl_timer)              :: timer

  type(scStat)                 :: stats

  real(fd), dimension(:,:,:,:), allocatable :: scatterData

  real(fd) :: dStd, cStd, dRough, cRough

  real    :: sTime, cTime, eTime
  integer :: i, iField, iStd, iRough

  integer :: fileID, dataID, dimID(4)

  namelist /params/ &
       & gridResHorizontal, gridResVertical, &
       & rf_spectrumType, rf_nFieldsPerMed, &
       & mediumFilename,  brdfType, &
       & sampleSeed, nOrders, nSamplesPerOrder,&
       & thetaIMax, thetaIMin, thetaIN, &
       & thetaEMax, thetaEMin, thetaEN, &
       & roughPMax, roughPMin, roughPN, &
       & srfStdMax, srfStdMin, srfStdN, &
       & medDensity, outFileName

  !$ call omp_set_num_threads(1)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! INITIALIZE SIMULATION PARAMETERS
  !!

  if(command_argument_count() == 0) then
     write(*,'("No parameter file given, using default values.")')
  else if(command_argument_count() == 1) then
     call get_command_argument(1, parFilename)
     write(*,'(("Using parameter file "),(A))') parFilename
  else
     write(*,'("Usage: vScatter paramFile")')
     stop
  end if

  !! Read the namelist from the parameter file.
  !!
  if(parFilename /= "") then
     open(1,file=parFilename,status="old", form="formatted")
     read(1,NML=params)
     close(1)
  end if
  
 
  !! Read the number of samples for each scattering order into a table.
  !!
  allocate(nSamplesPerOrderTable(nOrders))
  read(nSamplesPerOrder,*) nSamplesPerOrderTable

  !! Allocate and initialize 
  !!
  allocate(scatterData(thetaIN, thetaEN, srfStdN, roughPN))

  scatterData = 0.0_fd


  !! Select the media for the simulation from a file.
  !!
  call med_mediumFileOpen(Mf, mediumFilename)
  call med_mediumFileSelectMedia(Mf)

 
  call rf_init(RF, 1000, 10.0_fd)

  call utl_timer_init(timer, 5.0_fd, Mf%nSelectedMedia * thetaEN * srfStdN * roughPN * rf_nFieldsPerMed)
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! RUN SIMULATION
  !!

  call openDataFile(outFileName, fileID, dataID, dimID)

  call utl_message("BRDF type used: " // brdfType)

  call utl_message("Beginning trace.")

  call med_gridAssign(M, gridResHorizontal, gridResVertical)

  dStd   = (srfStdMax - srfStdMin) / real(srfStdN-1, fd) 
  dRough = (roughPMax - roughPMin) / real(roughPN-1, fd) 

  do iRough = 1, roughPN
     cRough = roughPMin + (iRough-1) * dRough

     do iStd = 1, srfStdN
        cStd = srfStdMin + (iStd-1) * dStd

        select case(rf_spectrumType)
        case('Gaussian')
           call rf_generateSpectrum (RF, RF_SPECTRUM_GAUSSIAN, cRough, cStd)
        case('fBm')
           call rf_generateSpectrum (RF, RF_SPECTRUM_FBM, cRough, cStd)
        end select

        do i = 1, Mf%nSelectedMedia

           call med_mediumFileRead(M0, Mf, Mf%varSelection(i))

           do iField = 1, rf_nFieldsPerMed

              call rf_generateField(RF, iField + rf_nFieldsPerMed*i)

              call med_maskHeightToNew(M0, M, real(rf%field(1:1000,1:1000), fs))

              select case(brdfType)
              case("shadowing")
                 call sampleAngles(brdf_shadowing, brdf_shadowing_f, iRough, iStd)
              case("Lambert")
                 call sampleAngles(brdf_Lambert, brdf_Lambert_f, iRough, iStd)
              case("LommelSeeliger")
                 call sampleAngles(brdf_LommelSeeliger, brdf_LommelSeeliger_f, iRough, iStd)
              case default
                 call utl_fatal_error("Unsupported BRDF type.")
              end select

           end do
        end do

        scatterData = (1.0_fd / (real(rf_nFieldsPerMed, fd) * real(MF%nSelectedMedia,fd))) * scatterData

        call saveData(scatterData, iStd, iRough)

     end do
  end do


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! WRITE FILES AND CLEAN UP
  !!

  !write(*,'(11F6.2)') scatterData(:,:,1,1)
  do iRough = 1, roughPN
     write(tmpStr,'(I1)'), iRough
     call utl_writeFits3D(scatterData(:,:,:,iRough), "!scatterData"//trim(tmpStr)//".fits")
  end do

  call med_mediumFileClose(Mf)

  call closeDataFile(fileID)

contains

  subroutine sampleAngles(brdf, brdf_f, iRough, iStd)

    real(fd), external  :: brdf, brdf_f

    type(ray) :: rC
    type(intersection_geometry) :: iSect

    real(fd), dimension(2,nSamplesPerOrderTable(1)) :: samples
    real(fd), dimension(3)          :: pSurface, tempD
    real(fd), dimension(nOrders)    :: IllTemp
    real(fd), dimension(3,thetaIN)  :: D
    real(fd)                        :: dz, dTI, dTE
    integer  :: i, j, k
    logical  :: pFound, pLit

    integer  :: iE, iI, iRough, iStd, iSample


    !! Precompute the angles of incidence.
    !!
    dTI = (thetaIMax - thetaIMin) / real(thetaIN-1,fd)
    dTE = (thetaEMax - thetaEMin) / real(thetaEN-1,fd)

    do iI = 1, thetaIN

       D(1,iI) = sin((thetaIMin + (iI-1) * dTI) * d_r)
       D(2,iI) = 0.0_fd
       D(3,iI) = cos((thetaIMin + (iI-1) * dTI) * d_r)

    end do

    call smpl_griddedSamples2D(samples, nSamplesPerOrderTable(1), sampleSeed)

    samples = (samples * M%width - M%hWidth) * 0.5_fd

    !$omp parallel default(none)                  &
    !$omp shared(M, D, samples, nOrders, nSamplesPerOrderTable, thetaEN, thetaEMin, thetaIN) &
    !$omp shared(iStd, iRough, dTI, dTE, brdf, brdf_f, timer, scatterData, stats) &
    !$omp private(rC, dz, pSurface, iE, iI, iSample, pFound, pLit, iSect, illTemp, tempD)
 
    dz = M%grid%height - M%hMean - TRACE_EPS

    pSurface(3) = M%hMean

    call ray_init(rC, RAY_TYPE_CAMERA)
  
    !$omp do schedule(dynamic)
    do iE = 1, thetaEN
       do iSample= 1, nSamplesPerOrderTable(1)
          pFound = .false.

          !! Select the sample point of the incident camera ray from
          !! the mean medium surface.
          !!
          pSurface(1:2) = samples(:,iSample)

          do while (.not. pFound)

             !! Find the intersection of the camera ray and the top of 
             !! the bounding box of the periodic medium.
             !!

             rC % D(1)  = sin((thetaEMin + (iE-1) * dTE) * d_r)
             rC % D(2)  = 0.0_fd
             rC % D(3)  = cos((thetaEMin + (iE-1) * dTE) * d_r)

             if(rC%D(3) < 1e-2_fd) then
                rC%D(3) = 1e-2
                call vec_normalize(rC%D)
             end if

             !if (iSample == 1) print *, rC%D

             rC % P(1)     = pSurface(1) + dz * (rC%D(1) / rC%D(3))
             rC % P(2)     = pSurface(2) + dz * (rC%D(2) / rC%D(3))
             rC % P(3)     = M%grid%height - TRACE_EPS

             rC % P(1)     = modulo(rC%P(1)+M%hWidth, M%width) - M%hWidth
             rC % P(2)     = modulo(rC%P(2)+M%hWidth, M%width) - M%hWidth

             rC % rayID    = rC % rayID + RAY_ID_INCR

             !! Negate the ray.
             !!
             rC % D        = - rC % D


             !! Find the true intersection point of the camera ray and medium.
             !!
             pFound        = trc_traceNearest(M%grid, rC, iSect)

             !! If an intersection is found, check if the point is shadowed.
             !!
             if(pFound) then

                do iI = 1, thetaIN
                   illTemp = 0.0_fd

                   tempD = D(:,iI)

                   call trc_gatherRadiance(M%grid, rC%D, tempD, iSect%P1 + TRACE_EPS * iSect%N, &
                        & iSect%N, 1.0_fd, &
                        & nSamplesPerOrderTable, nOrders, 1, illTemp, brdf, brdf_f)

                   scatterData(iI,iE,iStd,iRough) = scatterData(iI,iE,iStd,iRough) + illTemp(1)

                end do

             else

                call dst_generate_uniform(pSurface(1:2), 0.0_fd, 1.0_fd)
                pSurface(1:2) = (pSurface(1:2) * M%width - M%hWidth) * 0.5_fd

                stats%raysMissed = stats%raysMissed + 1

             end if

          end do

       end do
       call utl_timerIncrease(timer)

    end do
    !$omp end do
    !$omp end parallel

  end subroutine sampleAngles


  subroutine rayToBBTop(r, M, pSurface, dz)
    type(ray)        :: r
    type(med_medium) :: M
    real(fd)         :: dz, pSurface(3)

    r % P(1)     = pSurface(1) + dz * (r%D(1) / r%D(3))
    r % P(2)     = pSurface(2) + dz * (r%D(2) / r%D(3))
    r % P(3)     = M%grid%height - TRACE_EPS

    r % P(1)     = modulo(r%P(1)+M%hWidth, M%width) - M%hWidth
    r % P(2)     = modulo(r%P(2)+M%hWidth, M%width) - M%hWidth

    r % rayID    = r % rayID + RAY_ID_INCR

    r % D        = - r % D

  end subroutine rayToBBTop

 subroutine openDataFile(fName, fileID, dataID, dimID)

    character (len = *),          intent(IN)  :: fName

    real(fd)         :: eTable(thetaEN), iTable(thetaIN), sTable(srfStdN), rTable(roughPN)

    integer          :: fileID, dataID, dimID(4), eId, eDimId, iId, iDimId, i
    integer          :: rId, rDimId, sId, sDimId
    real(fd)         :: t

    t = (thetaIMax - thetaIMin) / real(thetaIN-1, fd)
    do i = 1, thetaIN
       iTable(i) = thetaIMin + (i-1)*t
    end do

    t = (thetaEMax - thetaEMin) / real(thetaEN-1, fd)
    do i = 1, thetaEN
       eTable(i) = thetaEMin + (i-1)*t
    end do

    t = (srfStdMax - srfStdMin) / real(srfStdN-1, fd)
    do i = 1, srfStdN
       sTable(i) = srfStdMin + (i-1)*t
    end do

    t = (roughPMax - roughPMin) / real(roughPN-1, fd)
    do i = 1, roughPN
       rTable(i) = roughPMin + (i-1)*t
    end do


    call nfCheck( nf90_create(fName, NF90_CLOBBER, fileID) )

    call nfCheck( nf90_def_dim(fileID, "thetaI", thetaIN, dimID(1)) )
    call nfCheck( nf90_def_dim(fileID, "thetaE", thetaEN, dimID(2)) )
    call nfCheck( nf90_def_dim(fileID, "std",    srfStdN, dimID(3)) )
    call nfCheck( nf90_def_dim(fileID, "roughP", roughPN, dimID(4)) )

    call nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "Author",  "Hannu Parviainen" ) )
    call nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "Program", "simpleScatter"    ) )
    call nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "Version",  "1.0"             ) )

    call nfCheck( nf90_def_var(fileID, "scatterData", NF90_DOUBLE, dimID, dataID  ) )
    call nfCheck( nf90_def_var(fileID, "thetaI",      NF90_DOUBLE, dimID(1), iId  ) )
    call nfCheck( nf90_def_var(fileID, "thetaE",      NF90_DOUBLE, dimID(2), eId  ) )
    call nfCheck( nf90_def_var(fileID, "std",         NF90_DOUBLE, dimID(3), sId  ) )
    call nfCheck( nf90_def_var(fileID, "roughP",      NF90_DOUBLE, dimID(4), rId  ) )

    call nfCheck( nf90_put_att(fileID, iId, "units",  "degrees"                   ) )
    call nfCheck( nf90_put_att(fileID, eId, "units",  "degrees"                   ) )

    call nfCheck( nf90_put_att(fileID, dataID, "brdfType",  brdfType              ) )
    call nfCheck( nf90_put_att(fileID, dataID, "roughnessType",  rf_spectrumType  ) )
    call nfCheck( nf90_put_att(fileID, dataID, "volumeDensity",  [medDensity]     ) )

    call nfCheck( nf90_enddef (fileID) )

    call nfCheck( nf90_put_var(fileID, iID, iTable) )
    call nfCheck( nf90_put_var(fileID, eID, eTable) )
    call nfCheck( nf90_put_var(fileID, sID, sTable) )
    call nfCheck( nf90_put_var(fileID, rID, rTable) )

  end subroutine openDataFile


  subroutine saveData(data, sIdx, rIdx)
    real(fd), dimension(:,:,:,:) :: data
    integer sIdx, rIdx
    write(*,*) "Writing data", sIdx, rIdx
    call nfCheck( nf90_put_var(fileID, dataID, data(:,:,sIdx,rIdx), start=[1,1,sIdx,rIdx], count=[thetaIN, thetaEN,1,1]) )
  end subroutine saveData


  subroutine closeDataFile(fileID)
    integer          :: fileID
    call nfCheck( nf90_close(fileID) )
  end subroutine closeDataFile

end program simpleScatter

