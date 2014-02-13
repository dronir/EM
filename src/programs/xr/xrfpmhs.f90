!****h* XR/programs/xrfpmhs
! NAME
!   xrfpmhs
!
! DESCRIPTION
!   Computes soft X-ray fluorescence of particulate media.
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   26.07.2008
!******

program xrfpmhs
  use base
  use material
  use geometry
  use particle
  use hemisphere
  use random
  use medium
  use rfield
  use trace
  use sampler


  !$ use omp_lib

  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! INPUT PARAMETERS
  !!

  !! FILES
  !!
  character (len=FNAME_LENGTH) :: parFilename        = ""
  character (len=FNAME_LENGTH) :: outFilename        = "xrfpn.nc"
  character (len=FNAME_LENGTH) :: spcFilename        = ""
  character (len=FNAME_LENGTH) :: mediumFileName     = ""


  !! GRID ACCELERATOR
  !!
  integer                      :: gridResHorizontal  = 200
  integer                      :: gridResVertical    = 10

  !!
  !!
  integer                      :: nThetaIn           = 1
  character (len=10000)        :: thetaIn            = "45.0"
  real(fd), allocatable        :: thetaInTable(:)


  !! material PARAMETERS
  !!
  character (len=256)          :: matName            = "Iron"

  real(fd)                     :: numberDensity      = 1.0

  integer                      :: nElements          = 2
  character (len=10000)        :: elements           = "Fe Ca"
  character (len=10000)        :: elemFracs          = "0.5 0.5"
  character(2), allocatable    :: elementTable(:)
  real(fd), allocatable        :: elemFracTable(:)


  !! SIMULATION
  !!
  character (len=100)          :: method             = "MonteCarlo"
  integer                      :: nOrders            = 1
  character(len=100)           :: nSamplesPerOrder   = "1"
  integer, allocatable         :: nSamplesPerOrderTable(:)
  integer                      :: sampleSeed         = 0

  integer                      :: nThreads          = 4


  !! SPECTRUM
  !!
  character (len=100)          :: sourceType         = "Spectrum"


  !! HEMISPHERE
  !!
  integer                      :: resThetaE          = 25
  integer(il)                  :: nSamples           = 50000


  !! RANDOM FIELDS
  !!
  logical                      :: rf_applyFields    = .false.
  integer                      :: rf_nFieldsPerMed  = 1
  character(len=100)           :: rf_spectrumType   = "Gaussian"
  real(fd)                     :: rf_P              = 0.5_fd
  real(fd)                     :: rf_std            = 1.0_fd

  integer                      :: mediumMapRes        = 128
  integer                      :: mediumDensitymapRes = 300


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! SIMULATION VARIABLES
  !!


  type(gth_hemisphere)         :: H
  type(spc_spectrum)           :: spc
  type(med_mediumFile)         :: Mf
  type(rf_randomField)         :: RF

  type(med_medium)             :: med
  type(mat_material), target   :: mat

  real, dimension(:,:,:), allocatable :: mediumHeightMap
  real(fd), dimension(:,:,:), allocatable :: mediumDensityStructure

  real    :: sTime, cTime, eTime
  integer(il) :: i, iField
  
  namelist /params/ &
       & mediumFileName, outFilename, spcFilename, &
       & nThetaIn, thetaIn, &
       & nElements, elements, elemFracs, numberDensity, &
       & resThetaE, &
       & gridResHorizontal,  gridResVertical, &
       & nSamples, &
       & matName, &
       & nOrders,   &
       & method,    &
       & sourceType, &
       & rf_applyFields, rf_nFieldsPerMed, rf_spectrumType, rf_P, rf_std, &
       & nThreads, &
       & mediumMapRes, mediumDensitymapRes



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! READ THE INPUT FILE
  !!

  !$ call omp_set_num_threads(nThreads)

  if(command_argument_count() == 0) then
     write(*,'("No parameter file given, using default values.")')
  else if(command_argument_count() == 1) then
     call get_command_argument(1, parFilename)
     write(*,'(("Using parameter file "),(A))') parFilename
  else
     write(*,'("Usage: xrfpm parameterFile")')
     stop
  end if

 !! Read the namelist from the parameter file.
  !!
  if(parFilename /= "") then
     open(1,file=parFilename,status="old", form="formatted")
     read(1,NML=params)
     close(1)
  end if

  !! Read the incidence angles.
  !!
  allocate(thetaInTable(nThetaIn))
  read(thetaIn,*) thetaInTable
  thetaInTable = thetaInTable * CNV_DEGTORAD
  
  !! Read the elements used by the material.
  !!
  allocate(elementTable(nElements))
  allocate(elemFracTable(nElements))

  read(elements,*) elementTable
  read(elemFracs,*) elemFracTable


  !! Read the number of samples for each scattering order into a table.
  !!
  !allocate(nSamplesPerOrderTable(nOrders))
  !read(nSamplesPerOrder,*) nSamplesPerOrderTable

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! INITIALIZE VARIABLES
  !!

  call dst_init(0)
  call mat_material_init(mat, matName, elementTable, numberDensity, elemFracTable, 1.5e4_fd)

  call spc_spectrum_read(spc, spcFilename, minval(mat%eEnergy), 1.45e4_fd)
  spc%I = 1e4_fd*spc%I

  call gth_hemisphere_init(H, resThetaE, nThetaIn, mat%nLines, GTH_TYPE_HS)
  write(H%name,'("Hemisphere")')


  !! Select the media for the simulation from a file.
  !!
  call med_mediumFileOpen(Mf, mediumFilename)
  call med_mediumFileSelectMedia(Mf)

  if(rf_applyFields) then

     call rf_init(RF, 500, 5.0_fd)

     select case(rf_spectrumType)
     case('Gaussian')
        call rf_generateSpectrum (RF, RF_SPECTRUM_GAUSSIAN, rf_P, rf_std)
     case('fBm')
        call rf_generateSpectrum (RF, RF_SPECTRUM_FBM, rf_P, rf_std)
     end select
  end if


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! RUN SIMULATION
  !!

  allocate(mediumHeightMap(MF%nSelectedMedia, mediumMapRes, mediumMapRes))
  allocate(mediumDensityStructure(MF%nSelectedMedia, mediumDensitymapRes, 2))

  do i = 1, MF%nSelectedMedia
     do iField = 1, rf_nFieldsPerMed
        call med_mediumFileRead(med, Mf, Mf%varSelection(i))

  print *, med%pDistribution
  print *, med%pDistParms

        call med_gridAssign(med, gridResHorizontal, gridResVertical)
        call med_gridFit(med)

        if(rf_applyFields) then
           call rf_generateField(RF, iField)
           call med_maskHeight(med, real(rf%field, fs))
        end if
 
        print *,"A"
        mediumHeightMap(i,:,:) = cnt_grid3D_cmpHeightMap(med%grid, mediumMapRes)
        call med_computePorosityStructure(med, 5000, mediumDensityStructure(i,:,:))
        print *,"B"


        call cpu_time(sTime) 
        select case(method)
        case('MonteCarlo')
           call simulation_method_montecarlo(med, mat, thetaInTAble, H, nSamples)
        case('FirstOrder')
           call simulation_method_firstorder(med, mat, thetaInTAble, H, nSamples)
        case('Analytic')
           call simulation_method_analytic(mat, spc, H, thetaInTAble)
        case default
           call utl_fatal_error("Error: Unknown simulation method: " // method)
        end select
        call cpu_time(cTime)
     end do
  end do

  write(*,'("Time taken: " ,(F6.3))') cTime - sTime

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! WRITE FILE
  !!

  call saveHemisphere(H,  mat, outFilename)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! UTILITY FUNCTIONS AND SUBROUTINES
  !!
  !!

  !! First order analytic
  !!

  subroutine simulation_method_analytic(m, s, hs, theta_i)
    type(mat_material)       :: m
    type(spc_spectrum)       :: s
    type(gth_hemisphere)     :: hs

    real(fd), dimension(:)   :: theta_i

    real(fd)                 :: e, spcNorm

    real(fd)    :: muInTot, muInExt, muInLine(m%nLines), muEmLine(m%nLines), fYield(m%nLines)
    real(fd)    :: A(m%nLines), B(m%nLines), fOut(m%nLines), xEn

    integer :: iTht, iThtE, iLine, iSpec, iEn

 
    do iSpec = 1, s%nPoints

       e = s%E(iSpec)

       spcNorm = s%I(iSpec) / sum(s%I)

       muInTot  = mat_evalMuAbs(m, e, iEn, xEn)
       muInLine = ((1.0_fd - xEn) * m%muFluorLine(:,iEn) + xEn * m%muFluorLine(:,iEn+1))

       muInExt  = mat_evalMuExt(m, e) - muInTot

       do iLine = 1, m%nLines
          if(m%lEnergy(iLine) > 0.0_fd) then
             muEmLine(iLine) = mat_evalMuExt(m, m%lEnergy(iLine))
             fYield(iLine)   = m%fYield((iLine-1)/2 + 1)
          end if
       end do

       A = spcNorm * muInLine / muInTot * fYield
       B = muEmLine / muInTot

       do iTht = 1, size(theta_i)
          do iThtE = 1, hs%resTheta

             fOut = INV_FOUR_PI * A * cos(hs%mTheta(iThtE)) / &
                  & ( cos(hs%mTheta(iThtE)) + B * cos(theta_i(iTht)))

             !$omp critical
             do iLine = 1, m%nLines
                call gth_hs_addToRow(hs, iThtE, fOut(iLine), iTht, iLine)
             end do
             !$omp end critical

          end do
       end do
    end do

  end subroutine simulation_method_analytic


  subroutine simulation_method_firstorder(med, mat, theta_i, hs, nSamples)
    type(med_medium)       :: med
    type(mat_material)     :: mat
    real(FD), dimension(:) :: theta_i
    type(gth_hemisphere)   :: hs

    integer(il)            :: nSamples

    integer, parameter :: RAY_OK   = 0
    integer, parameter :: RAY_EXIT = 1
    integer, parameter :: RAY_LOW  = 2
    integer, parameter :: SRTABLE  = 20000


    ! Simulation (logical) variables
    !
    integer(il)                     :: iRay, iRand, iLine, iTht
    integer                         :: rayStat, order

    real(fd), dimension(2,nSamples) :: samples
    real(fd), dimension(3)          :: pSurface
    real(fd)                        :: dz, xEn, pTemp(3)
    integer                         :: iEn, iThtE, iPhiE

    logical                         :: isInside


    type(intersection_geometry)     :: iSect
    type(utl_timer)                 :: timer

    ! Physical variables
    !
    real(fd)    :: muIonization, muExtinction, l(1), tau
    real(fd), dimension(mat%nLines) :: lineExtinction, lineFluorescenceYield, lineRatios, lineFluorescence
    type(ray)   :: r_s, r, rI

    ! Random number generator variables
    !
    real(fd)    :: rTable(SRTABLE)
    integer     :: rInfo, rSeed(1), rState(640)

    !$call omp_set_num_threads(nThreads)

    call smpl_griddedSamples2D(samples, nSamples, sampleSeed)

    !! Compute the per-line coefficients (total extinction coefficient and fluorescence yield)
    !!
    do iLine = 1, mat%nLines
       if(mat%lEnergy(iLine) > 0.0_fd) then
          lineExtinction(iLine)        = mat_evalMuExt(mat, mat%lEnergy(iLine))
          lineFluorescenceYield(iLine) = mat%fYield((iLine-1)/2 + 1)
       end if
    end do
 
    call utl_timer_init(timer, 1.0_fd, size(theta_i))
    do iTht = 1, size(theta_i)

       pSurface(3)  = med%hMean
       dz           = med%grid%height - med%hMean - TRACE_EPS

       r_s % P      = med%grid%height - TRACE_EPS
       r_s % D      = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]

       iRand = 1
       rSeed = 1
       
       !$omp parallel default(none)&
       !$omp shared(r_s, iTht, nSamples, med, mat, spc, hs, theta_i, samples, lineExtinction, lineFluorescenceYield) &
       !$omp private(iRay, iLine, rayStat, dz, r, tau, muExtinction, muIonization, iEn, xEn, l, rTable) &
       !$omp private(rSeed, rState, rInfo, lineRatios, lineFluorescence, pTemp, iThtE, iPhiE) &
       !$omp firstprivate(iRand)

       !$ rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)
       
 
       !$omp do schedule(dynamic) 
       do iRay = 1, nSamples
          
          if(iRand == 1) call dranduniform(SRTABLE, 0.0_fd, 1.0_fd, rState, rTable, rInfo)
          iRand = mod(iRand, SRTABLE) + 1

          r % I      = 1.0_fd / real(nSamples, fd)
          r % energy = spc_getSample(spc, rTable(iRand))

          r % status = RAY_ACTIVE
          r % D      = r_s%D
          r % P(3)   = r_s%P(3)
          r % P(1:2) = samples(:, iRay) + dz * (r%D(1:2) / r%D(3))
          r % P(1:2) = modulo(r%P(1:2)+med%hWidth, med%width) - med%hWidth

          r % rayID  = r % rayID + RAY_ID_INCR

          rayStat    = RAY_OK
          
          muIonization  = mat_evalMuAbs(mat, r%energy, iEn, xEn)
          muExtinction  = mat_evalMuExt(mat, r%energy) - muIonization

          call drandexponential(1, 1.0_fd, rState, l, rInfo)
          call trc_traceRayToOpticalDepth(med%grid, r, l(1), mat)

          r % rayID  = r % rayID + RAY_ID_INCR

          if(r%status == RAY_EXIT_D) rayStat = RAY_EXIT

          if(rayStat == RAY_OK) then

             lineRatios    = ((1.0_fd - xEn) * mat%muFluorLine(:,iEn) + xEn * mat%muFluorLine(:,iEn+1)) / muIonization
             pTemp         = r%P

             do iThtE = 1, hs%resTheta
                do iPhiE = 1, hs%resPhi(iThtE)

                   r % status = RAY_ACTIVE
                   r%P = pTemp
                   r%D = gth_cellRandomSampleCar(hs, iThtE, iPhiE)

                   !! r%D    = vec_cart_random_hemispherical_uniform()
                   
                   !! We need to enable (and fix) this if we want to use multiple materials in the medium.
                   !!
                   !! tau    = trc_traceOpticalDepth(med%grid, mat, r)

                   tau        = trc_tracePhysicalDepth(med%grid, r)

                   lineFluorescence = INV_FOUR_PI * r%I * lineFluorescenceYield * lineRatios * exp(-lineExtinction * tau)

!!$                   print *,r_s%P
!!$                   print *,pTemp
!!$                   print *,r%D

                   !$omp critical 
                   do iLine = 1, mat%nLines
                      call gth_addDataCar(hs, r%D, lineFluorescence(iLine), iTht, iLine)
                   end do
                   !$omp end critical
                end do
             end do

          end if
       end do
       !$omp end do
       !$omp end parallel

       call utl_timerIncrease(timer)
    end do
  end subroutine simulation_method_firstorder


  subroutine simulation_method_montecarlo(med, mat, theta_i, hs, nSamples)
    type(med_medium)       :: med
    type(mat_material)     :: mat
    real(FD), dimension(:) :: theta_i
    type(gth_hemisphere)   :: hs

    integer(il)            :: nSamples

    integer, parameter :: RAY_OK   = 0
    integer, parameter :: RAY_EXIT = 1
    integer, parameter :: RAY_LOW  = 2
    integer, parameter :: SRTABLE  = 20000


    ! Simulation (logical) variables
    !
    integer(il)                     :: iRay, iRand, iLine, iTht
    integer                         :: rayStat, order

    real(fd), dimension(2,nSamples) :: samples
    real(fd), dimension(3)          :: pSurface
    real(fd)                        :: dz

    logical                         :: isInside


    type(intersection_geometry)     :: iSect
    type(utl_timer)                 :: timer

    ! Physical variables
    !
    real(fd)    :: muFluorLine, muExt, fYield, l(1), tau
    type(ray)   :: r_s, r, rI

    ! Random number generator variables
    !
    real(fd)    :: rTable(SRTABLE)
    integer     :: rInfo, rSeed(1), rState(640)

    !$ call omp_set_num_threads(nThreads)

    call smpl_griddedSamples2D(samples, nSamples, sampleSeed)

    call utl_timer_init(timer, 1.0_fd, size(theta_i))
    do iTht = 1, size(theta_i)

       pSurface(3)  = med%hMean
       dz           = med%grid%height - med%hMean - TRACE_EPS

       r_s % P      = med%grid%height - TRACE_EPS
       r_s % D      = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]

       iRand = 1
       rSeed = 1
       
       !$omp parallel default(none)&
       !$omp shared(r_s, iTht, nSamples, med, mat, spc, hs, theta_i, nOrders, pSurface, samples) &
       !$omp private(iRay, iLine, rayStat, dz, r, tau, muExt, muFluorLine, l, fYield, rTable, rSeed, rState, rInfo, order) &
       !$omp firstprivate(iRand)

       !$rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)
       
       !$omp do schedule(dynamic) 
       do iRay = 1, nSamples
          
          if(iRand == 1) call dranduniform(SRTABLE, 0.0_fd, 1.0_fd, rState, rTable, rInfo)
          iRand = mod(iRand, SRTABLE) + 1


          !! INITIALIZE RAY
          !!

          pSurface(1:2) = samples(:,iRay)

          r          = r_s
          r % I      = 1.0_fd / real(nSamples, fd)
          r % energy = spc_getSample(spc, rTable(iRand))

          r % P(1)   = pSurface(1) + dz * (r%D(1) / r%D(3))
          r % P(2)   = pSurface(2) + dz * (r%D(2) / r%D(3))

          r % P(1)   = modulo(r%P(1)+med%hWidth, med%width) - med%hWidth
          r % P(2)   = modulo(r%P(2)+med%hWidth, med%width) - med%hWidth

          r % rayID  = r % rayID + RAY_ID_INCR

          rayStat    = RAY_OK
          order      = 1                

          do while(order <= nOrders .and. rayStat == RAY_OK)

             muFluorLine  = mat_evalMuAbs(mat, r%energy)
             muExt        = mat_evalMuExt(mat, r%energy) - muFluorLine

             if(muFluorLine > 1e-22_fd) then

                call drandexponential(1, 1.0_fd, rState, l, rInfo)
                call trc_traceRayToOpticalDepth(med%grid, r, l(1), mat)
                r % rayID  = r % rayID + RAY_ID_INCR
          

                if(r%status == RAY_EXIT_D) then
                   rayStat = RAY_EXIT
                   exit
                end if

                !! Attenuate intinsity due to the extinction methods not used in the simulation
                !!
                r%I  = r%I * exp(- l(1) * muExt)


                !! Fluorescence
                !!
                if(rayStat == RAY_OK) then

                   call mat_evalFluorescence(mat, r%energy, r%energy, fYield, iLine)

                   !! Peeling
                   !!
                   muExt  = mat_evalMuExt(mat, r%energy)
                   r%D    = vec_cart_random_hemispherical_uniform()

                   tau    = trc_traceOpticalDepth(med%grid, mat, r)
                   !tau    = trc_tracePhysicalDepth(med%grid, r)
                   r % rayID  = r % rayID + RAY_ID_INCR
          

                   !$omp critical 
                   call gth_addDataCar(hs, r%D, INV_FOUR_PI * r%I * fYield * exp(-tau), iTht, iLine)
                   !$omp end critical

                   !! Fluorescence
                   !!
                   r%D  = vec_cart_random_spherical_uniform()
                   r%I  = r%I * fYield 

                end if

             !! In the case that the absorbtion coefficient would be
             !! very small, consider it zero.
             !!
             else
                if(r%D(3) <= 0.0_fd) then
                   rayStat = RAY_LOW
                else
                   rayStat = RAY_EXIT
                end if
             end if

             order = order + 1

          end do

!!$          if(rayStat == RAY_EXIT) then
!!$             !$omp critical 
!!$             call gth_addDataCar(hs, r%D, r%I, iTht, iLine)
!!$             !$omp end critical
!!$          end if

       end do
       !$omp end do
       !$omp end parallel

       call utl_timerIncrease(timer)
    end do
  end subroutine simulation_method_montecarlo


  subroutine saveHemisphere(h, m, fName)
    type(gth_hemisphere), INTENT(IN)  :: h
    type(mat_material), intent(IN)    :: m
    character (len = *)               :: fName

    type(gth_hsFile)   :: hf
    integer            :: i

    integer            :: dimMatMuEnergy, dimFlrLine, dimSpecEnergy, dimSpecICDF
    integer            :: idMatMuIon, idMatMuFluor, idMatMuFluorCdf,idMatMuRay, idMatMuExt, idMatMuEnergy
    integer            :: idSpecEnergy, idSpecIntensity, idSpecCdf, idSpecIcdf

    integer            :: dimMedNum, dimMedMapRes, dimMedDenRes, idMedMap, idMedDen

    integer, parameter :: RTP = NF90_DOUBLE
    integer, parameter :: RTS = NF90_FLOAT


    call gth_hsFileOpen(hf, fName, "w")
    call gth_hsWriteHeader(hf, h, "Hannu Parviainen", "xrfpmhs", "1.0")

    !! X-RAY SPECIFIC STUFF
    !!
    !!

    !! Dimensions
    !!
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Fluorescence_line",     m%nLines,           dimFlrLine     ) )
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Material_energy",       m%nE,               dimMatMuEnergy ) )

    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Spectrum_energy",       spc%nPoints,        dimSpecEnergy  ) )
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Spectrum_ICDF_energy",  spc%nPoints*5,      dimSpecICDF    ) )

    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Number_of_media",       MF%nSelectedMedia,  dimMedNum      ) )
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Medium_map_resolution", mediumMapRes,       dimMedMapRes   ) )
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Medium_densitymap_resolution", mediumDensitymapRes, dimMedDenRes   ) )


    !! Variables
    !!
!!$    Call gth_nfCheck( nf90_def_var(hf%fileID, "Photoionization_coefficient", RTP, [dimFlrLine,dimMatMuEnergy], idMatMuIon) )
!!$    Call gth_nfCheck( nf90_def_var(hf%fileID, "Fluorescence_line_coefficient", RTP, [dimFlrLine,dimMatMuEnergy], idMatMuFluor) )
!!$    Call gth_nfCheck( nf90_def_var(hf%fileID, "Fluorescence_line_cdf",         RTP, [dimFlrLine,dimMatMuEnergy], idMatMuFluorCdf))
!!$    !Call gth_nfCheck( nf90_def_var(hf%fileID, "Photoionization_cdf",         RTP, [dimFlrLine,dimMatMuEnergy], idMatMuIonCdf))
!!$    Call gth_nfCheck( nf90_def_var(hf%fileID, "Rayleigh_coefficient",        RTP, [dimMatMuEnergy],            idMatMuRay) )
!!$    Call gth_nfCheck( nf90_def_var(hf%fileID, "Extinction_coefficient",      RTP, [dimMatMuEnergy],            idMatMuExt) )
!!$    Call gth_nfCheck( nf90_def_var(hf%fileID, "Material_energy",             RTP, [dimMatMuEnergy],            idMatMuEnergy) )
!!$
!!$    Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_energy",             RTP, [dimSpecEnergy], idSpecEnergy) )
!!$    Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_intensity",          RTP, [dimSpecEnergy], idSpecIntensity) )
!!$    Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_cdf",                RTP, [dimSpecEnergy], idSpecCdf) )
!!$    Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_inverse_cdf",        RTP, [dimSpecICDF],   idSpecIcdf) )

    Call gth_nfCheck( nf90_def_var(hf%fileID, "Photoionization_coefficient",   RTP, [dimMatMuEnergy],            idMatMuIon) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Fluorescence_line_coefficient", RTP, [dimFlrLine,dimMatMuEnergy], idMatMuFluor) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Fluorescence_line_cdf",         RTP, [dimFlrLine,dimMatMuEnergy], idMatMuFluorCdf))
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Rayleigh_coefficient",          RTP, [dimMatMuEnergy],            idMatMuRay) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Extinction_coefficient",        RTP, [dimMatMuEnergy],            idMatMuExt) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Material_energy",               RTP, [dimMatMuEnergy],            idMatMuEnergy) )

    Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_energy",               RTP, [dimSpecEnergy], idSpecEnergy) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_intensity",            RTP, [dimSpecEnergy], idSpecIntensity) )
    if(sourceType == 'Spectrum') then
       Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_cdf",               RTP, [dimSpecEnergy], idSpecCdf) )
       Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_inverse_cdf",       RTP, [dimSpecICDF],   idSpecIcdf) )
    end if

    Call gth_nfCheck( nf90_def_var(hf%fileID, "Medium_heightmap",            RTS, [dimMedNum,dimMedMapRes,dimMedMapRes], idMedMap) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Medium_densitymap",           RTS, [dimMedNum,dimMedDenRes,2],            idMedDen) )

    !! Global attributes
    !!
    call gth_hsAddAttS(hf, NF90_GLOBAL, "Elements",           elements)
    call gth_hsAddAttS(hf, NF90_GLOBAL, "Spectrum_type",      sourceType)
    call gth_hsAddAttS(hf, NF90_GLOBAL, "Simulation_method",  method)
    call gth_hsAddAttF(hf, NF90_GLOBAL, "Theta_in",           thetaInTable)

    call gth_hsAddAttF(hf, NF90_GLOBAL, "Fluorescence_line_energy", m%lEnergy)
    call gth_hsAddAttF(hf, NF90_GLOBAL, "Absorbtion_edge_energy",   m%eEnergy)

    call gth_hsSaveData(h, hf)


    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuIon,       m%muPhotoIon) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuFluor,     m%muFluorLine) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuFluorCdf,  m%muFluorLineCDF) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuRay,       m%muRayleigh) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuExt,       m%muExtTotal) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuEnergy,    m%E) )

    Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecEnergy,     spc%E) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecIntensity,  spc%I) )
    if(sourceType == 'Spectrum') then
       Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecCdf,     spc%cdf) )
       Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecIcdf,    spc%icdf) )
    end if


!!$    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuIon,       mat%muPhotoIon) )
!!$    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuFluor,     mat%muFluorLine) )
!!$    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuFluorCdf,  mat%muFluorLineCDF) )
!!$!    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuIon,       m%muFluorLine) )
!!$!    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuIonCdf,    m%muFluorLineCDF) )
!!$    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuRay,       m%muRayleigh) )
!!$    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuExt,       m%muExtTotal) )
!!$    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuEnergy,    m%E) )
!!$
!!$    Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecEnergy,     spc%E) )
!!$    Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecIntensity,  spc%I) )
!!$    Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecCdf,        spc%cdf) )
!!$    Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecIcdf,       spc%icdf) )

    Call gth_nfCheck( nf90_put_var(hf%fileID, idMedMap,         mediumHeightMap) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMedDen,         mediumDensityStructure) )


    call gth_hsFileClose(hf)

  end subroutine saveHemisphere

end program xrfpmhs

