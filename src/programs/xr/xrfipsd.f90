!****h* XR/programs/xrfipsd
! NAME
!   xrfipsd
!
! DESCRIPTION
!   Computes fluorescence for an infinite slab.
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   15.11.2007
!******

program xrfipsd
  use base
  use material
  use geometry
  use particle
  use random
  use hemisphere
  use bfdf

!$use omp_lib

  implicit none

  logical                      :: printDeviation     = .false.

  character (len=FNAME_LENGTH) :: parFilename        = ""
  character (len=FNAME_LENGTH) :: outFilename        = "xrfipsd.nc"

  character (len=100)          :: sourceType         = "Spectrum"
  character (len=FNAME_LENGTH) :: spcFilename        = ""
  real(fd)                     :: spcLineEnergy      = 9.0e3

  integer                      :: nThetaIn           = 1
  character (len=10000)        :: thetaIn            = "45.0"
  real(fd), allocatable        :: thetaInTable(:)

  real(fd)                     :: thetaE             = 0.0

  real(fd)                     :: obsHalfAngle       = 5.0
  real(fd)                     :: obsSolidAngle      = 1.0

  character (len=100)          :: method             = "MonteCarlo"

  character (len=256)          :: matName            = "Default Iron/Calcium"
  integer                      :: nElements          = 2
  character (len=10000)        :: elements           = "Fe Ca"
  character (len=10000)        :: elemFracs          = "0.5 0.5"
  character(2), allocatable    :: elementTable(:)
  real(fd), allocatable        :: elemFracTable(:)
  real(fd)                     :: molarDensity       = 1.0_fd

  real(fd)                     :: energyMin          = 0.0_fd
  real(fd)                     :: energyMax          = 1.0e4_fd

  integer(il)                  :: nSampleMultiplier  = 1
  integer(il)                  :: nSamples           = 5000
  integer                      :: nOrders            = 1

  real(fd), allocatable        :: fResults(:,:), fResultsAn(:,:), FresultsMC(:,:)

  real(fd)                     :: testSpec(500)

  type(spc_spectrum)           :: spc
  type(mat_material), target   :: mat

  character (len=100)          :: devFormat

  real    :: sTime, cTime, eTime
  integer(il) :: i
  
  namelist /params/ &
       & printDeviation, outFilename, spcFilename, &
       & nThetaIn, thetaIn, thetaE, obsHalfAngle, &
       & nElements, elements, elemFracs, molarDensity,&
       & nSamples, nSampleMultiplier, &
       & matName,   &
       & nOrders,   &
       & method,    &
       & sourceType, spcLineEnergy, &
       & energyMin, energyMax


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! READ THE INPUT FILE
  !!

  if(command_argument_count() == 0) then
     write(*,'("No parameter file given, using default values.")')
  else if(command_argument_count() == 1) then
     call get_command_argument(1, parFilename)
     write(*,'(("Using parameter file "),(A))') parFilename
  else
     write(*,'("Usage: xr_InfSlab paramFile")')
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
  thetaInTable  = thetaInTable * CNV_DEGTORAD
  
  thetaE        = thetaE * CNV_DEGTORAD

  obsHalfAngle  = obsHalfAngle * CNV_DEGTORAD
  obsSolidAngle = (1.0_fd - cos(obsHalfAngle)) * TWO_PI

  !! Read the elements used by the material.
  !!
  allocate(elementTable(nElements))
  allocate(elemFracTable(nElements))

  read(elements,*) elementTable
  read(elemFracs,*) elemFracTable


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! INITIALIZE
  !!

  call dst_init(0)
  call mat_material_init(mat, matName, elementTable, molarDensity, elemFracTable, energyMax)

  select case(sourceType)
     case('Spectrum')
        call spc_spectrum_read(spc, spcFilename, mat%edgeEnergyMin, energyMax)
     case('Line')
        call spc_spectrum_initline(spc, spcLineEnergy, 1.0_fd)
     case default
        call utl_fatal_error("Error: Unknown spectrum type: " // sourceType)
     end select


  allocate(fResults(nThetaIn, mat%nLines))
  allocate(fResultsAn(nThetaIn, mat%nLines))

  fResults   = 0.0_fd
  fResultsAn = 0.0_fd

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! RUN SIMULATION
  !!

  call cpu_time(sTime) 
  select case(method)
     case('MonteCarlo')
        call infslab_method_montecarlo(mat, thetaInTAble, thetaE, nSamples)
     case('MonteCarlo_nopeel')
        call infslab_method_montecarlo_bruteforce(mat, thetaInTAble, thetaE, nSamples)
     case('FirstOrder')
        call infslab_method_firstorder(mat, spc, thetaInTAble, thetaE, nSamples)
   case('Analytic')
        call infslab_method_analytic(mat, spc, thetaInTAble, thetaE, fResults)
     case default
        call utl_fatal_error("Error: Unknown simulation method: " // method)
     end select
  call cpu_time(cTime)

  write(*,'("Time taken: " ,(F10.3))') cTime - sTime

  if(printDeviation) then

     write (devFormat,'("((A2),(2F10.1),(",(I2.2),"F8.2))")') nThetaIn

     call infslab_method_analytic(mat, spc, thetaInTable, thetaE, fResultsAn)

     print *
     print '("Simulation:")'
     do i = 1, mat%nLines
        print devFormat, elementTable((i+1)/2), mat%eEnergy((i+1)/2), mat%lEnergy(i), 1e6*fResults(:,i)
     end do

     print *
     print '("Analytical: ")'
     do i = 1, mat%nLines
        print devFormat, elementTable((i+1)/2), mat%eEnergy((i+1)/2), mat%lEnergy(i), 1e6*fResultsAn(:,i)
     end do

     print *
     print '("Error: ")'
     do i = 1, mat%nLines
        print devFormat, elementTable((i+1)/2), mat%eEnergy((i+1)/2), mat%lEnergy(i), &
             & 1e2 * (fResults(:,i) - fResultsAn(:,i)) / fResultsAn(:,i)
     end do

  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! WRITE FILE
  !!

  call testSpectrum(spc, testSpec, 100000)

  call saveData(fResults, mat, spc, outFilename)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! UTILITY FUNCTIONS AND SUBROUTINES
  !!
  !!


  !****f* xrfipsd/infslab_method_analytic
  ! NAME
  !   infslab_method_analytic
  !
  ! DESCRIPTION
  !   Computes the fluorescent signal from an infinite plane using first-order analytic approximation.
  !
  ! INPUTS
  !   mat      : mat_material   ... Material
  !   spc      : spc_scpectrum  ... Spectrum
  !   theta_i  : real(:)        ... Angle(s) of incidence
  !   theta_e  : real           ... Angle of emergence
  !
  ! OUTPUT
  !   resutls  : real(:,:)      ... Fluorescent signal for each line and incidence angle
  !
  ! SOURCE
  subroutine infslab_method_analytic(mat, spc, theta_i, theta_e, results)
    type(mat_material)       :: mat
    type(spc_spectrum)       :: spc
    real(fd), dimension(:)   :: theta_i
    real(fd)                 :: theta_e
    real(fd), dimension(:,:) :: results

    real(fd)    :: muLineTotal(mat%nLines), fYield(mat%nLines)
    integer     :: iTht, iLine

    do iLine = 1, mat%nLines
       if(mat%lEnergy(iLine) > 0.0_fd) then
          call mat_evalMu(mat, mat%lEnergy(iLine), muTotal = muLineTotal(iLine))
          fYield(iLine) = mat%fYield((iLine-1)/2 + 1)
       end if
    end do

    do iTht = 1, size(theta_i)
       results(iTht, :) = bfdf_Parviainen(mat, spc, theta_i(iTht), theta_e, muLineTotal, fYield)
    end do

  end subroutine infslab_method_analytic
  !******


  !****f* xrfipsd/infslab_method_firstorder
  ! NAME
  !   infslab_method_firstorder
  !
  ! DESCRIPTION
  !   Computes the fluorescent signal from an infinite plane using first-order Monte-Carlo approximation.
  !
  ! INPUTS
  !   mat      : mat_material ... Material
  !   theta_i  : real(:)      ... Angle(s) of incidence
  !   theta_e  : real         ... Angle of emergence
  !   nSamples : integer      ... Number of samples
  !
  ! TODO
  !   The sampling over the spectrum is using a brute force method, we sample over all of the spectrum points. 
  !   This should be changed to weighted Monte Carlo sampling.
  !
  ! SOURCE
  subroutine infslab_method_firstorder(mat, spc, theta_i, theta_e, nSamples)
    type(mat_material)     :: mat
    type(spc_spectrum)     :: spc
    real(FD), dimension(:) :: theta_i
    real(fd)               :: theta_e
    integer(il)            :: nSamples

    ! Simulation (logical) variables
    !
    integer(il) :: iSpec, iRay, iLine, iTht
    real(fd)    ::  x

    type(utl_timer) :: timer

    ! Physical variables
    !
    real(fd)    :: muFluor, muTotal, fluorescence_ratio, l(1)
    real(fd)    :: muLineExt(mat%nLines), fluorescence_intensity(mat%nLines), fYield(mat%nLines), line_ratio(mat%nLines)
    type(ray)   :: r_ref, r

    ! Random number generator variables
    !
    integer     :: rInfo, rSeed(1), rState(640)

    muLineExt = 0.0_fd

    do iLine = 1, mat%nLines
       if(mat%lEnergy(iLine) > 0.0_fd) then
          call mat_evalMu(mat, mat%lEnergy(iLine), muTotal = muLineExt(iLine))
          fYield(iLine)      = mat%fYield((iLine-1)/2 + 1)
       end if
    end do

    call utl_timer_init(timer, 10.0_fd, size(theta_i))

    do iTht = 1, size(theta_i)

       r_ref % P = 0.0_fd
       r_ref % D = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]
       r_ref % I = 1.0_fd / real(nSamples, fd)

       rSeed = 1

       !$omp parallel default(none) &
       !$omp shared       ( mat, spc, r_ref, nSamples, theta_i, theta_e, fResults, iTht ) &
       !$omp private      ( muFluor, muTotal,   iRay, iSpec                             ) &
       !$omp private      ( fluorescence_ratio, line_ratio, fluorescence_intensity      ) &
       !$omp private      ( r, l, x, rSeed, rState, rInfo                               ) &
       !$omp firstprivate ( muLineExt, fYield                                           )

       !$rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)

       !$omp do schedule(dynamic) 
       do iSpec = 1, spc%nPoints
          do iRay = 1, nSamples

             !call dranduniform     (1, 0.0_fd, 1.0_fd, rState, x, rInfo)
             call drandexponential (1, 1.0_fd,         rState, l, rInfo)

             r          = r_ref
             r % energy = spc%E(iSpec) !spc_getSample(spc, x)
             r % I      = spc%I(iSpec) / real(nSamples, fd)

             call mat_evalMu(mat, r%energy, muFluorLineTotal = muFluor, muTotal = muTotal)

             l    = l / muTotal
             r%P  = r%P + r%D * l(1)

             fluorescence_ratio     = muFluor / muTotal
             line_ratio             = mat_evalFluorLineFractions(mat, r%energy)

             fluorescence_intensity = INV_FOUR_PI * r%I * fYield * fluorescence_ratio * line_ratio * exp(-muLineExt * abs(r%P(3)) / cos(theta_e))

             !$omp critical
             fResults(iTht, :) = fResults(iTht, :) + fluorescence_intensity(:)
             !$omp end critical
          end do
       end do
       !$omp end do
       !$omp end parallel

       call utl_timerIncrease(timer)
    end do

    fResults= fResults / sum(spc%I)

  end subroutine infslab_method_firstorder
  !******


  subroutine infslab_method_montecarlo_o(m, theta_i, nSamples)
    type(mat_material)     :: m
    real(FD), dimension(:) :: theta_i
    integer(il)            :: nSamples

    integer, parameter :: RAY_OK   = 0
    integer, parameter :: RAY_EXIT = 1
    integer, parameter :: RAY_LOW  = 2
    integer, parameter :: SRTABLE  = 20000

    ! Simulation (logical) variables
    !
    integer(il) :: iRay, iRand, iLine, iTht
    integer     :: rayStat, order

    type(utl_timer) :: timer

    ! Physical variables
    !
    real(fd)    :: muAbs, muExt, fYield, l(1)
    type(ray)   :: r_ref, r

    ! Random number generator variables
    !
    real(fd)    :: rTable(SRTABLE)
    integer     :: rInfo, rSeed(1), rState(640)

    call utl_timer_init(timer, 1.0_fd, size(theta_i))
    do iTht = 1, size(theta_i)

       r_ref % P = -1.0e-8_fd
       r_ref % D = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]

       iRand     = 1
       rSeed     = 1
       
       !$omp parallel &
       !$omp shared(r_ref, nSamples, m, theta_i, fResults, nOrders) &
       !$omp private(iRay, rayStat, r, muAbs, muExt, l, fYield, rTable, rSeed, rState, rInfo, order) &
       !$omp firstprivate(iRand)

       !$rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)

       !$omp do schedule(static) 
       do iRay = 1, nSamples
   
          if(iRand == 1) call dranduniform(SRTABLE, 0.0_fd, 1.0_fd, rState, rTable, rInfo)
          iRand = mod(iRand, SRTABLE) + 1

          r          = r_ref
          r % I      = 1.0_fd / real(nSamples, fd)
          r % energy = spc_getSample(spc, rTable(iRand))

          rayStat    = RAY_OK
          order      = 1

          do while(rayStat == RAY_OK)

             muAbs  = mat_evalMuAbs(m, r%energy)
             muExt  = mat_evalMuExt(m, r%energy) - muAbs

             if(muAbs > 1e-22_fd) then

                call drandexponential(1, 1.0_fd, rState, l, rInfo)
                l = l / muAbs

                r%P  = r%P + r%D * l(1)
                r%I  = r%I * exp(- l(1) * muExt)

                !! Test whether we are still inside the medium.
                !!
                if(r%P(3) > 0.0_fd) then
                   rayStat = RAY_EXIT
                else if(r%I < 1.0e-10_fd) then
                   rayStat = RAY_LOW
                end if

                !! Fluorescence
                !!
                if(rayStat == RAY_OK) then

                   call mat_evalFluorescence(m, r%energy, r%energy, fYield, iLine)

                   if (order <= nOrders) then

                      !! Peeling
                      !!
                      muExt  = mat_evalMuExt(m, r%energy)

                      !$omp atomic
                      fResults(iTht, iLine) = fResults(iTht, iLine) &
                           & + INV_FOUR_PI * r%I * fYield * exp(-muExt * abs(r%P(3)))

                      !! Fluorescence
                      !!
                      r%D  = vec_cart_random_spherical_uniform()
                      r%I  = r%I * fYield 

                   else
                      rayStat = RAY_EXIT
                   end if

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

       end do
       !$omp end do
       !$omp end parallel

       call utl_timerIncrease(timer)

    end do
  end subroutine infslab_method_montecarlo_o


  !****f* xrfipsd/infslab_method_montecarlo
  ! NAME
  !   infslab_method_montecarlo
  !
  ! DESCRIPTION
  !   Computes the fluorescent signal from an infinite plane using Monte Carlo ray tracing with peeling.
  !
  ! INPUTS
  !   mat      : mat_material ... Material
  !   theta_i  : real(:)      ... Angle(s) of incidence
  !   theta_e  : real         ... Angle of emergence
  !   nSamples : integer      ... Number of samples
  !
  ! TODO
  !   The sampling over the spectrum is using a brute force method, we sample over all of the spectrum points. 
  !   This should be changed to weighted Monte Carlo sampling.
  !
  ! SOURCE
  subroutine infslab_method_montecarlo(mat, theta_i, theta_e, nSamples)
    type(mat_material)     :: mat
    real(FD), dimension(:) :: theta_i
    real(fd)               :: theta_e
    integer(il)            :: nSamples

    integer,  parameter    :: SRTABLE       = 20000
    real(fd), parameter    :: MIN_MU        = 1e-20_fd
    real(fd), parameter    :: MIN_INTENSITY = 1e-23_fd

    ! Simulation (logical) variables
    !
    integer(il) :: iSampleset, iSpec, iRay, iRand, iLine, iTht
    integer     :: rayStat, order
    real(fd)    :: dObs(3), cos_half_angle

    type(utl_timer) :: timer

    ! Physical variables
    !
    real(fd)    :: muFluorLine, muTotal, fYield, l(1)
    type(ray)   :: r_ref, r

    ! Random number generator variables
    !
    real(fd)    :: rTable(SRTABLE)
    integer     :: rInfo, rSeed(1), rState(640)

    dObs = [sin(theta_e), 0.0_fd, cos(theta_e)]

    cos_half_angle = cos(obsHalfAngle)

    call utl_timer_init(timer, 5.0_fd, size(theta_i) * nSampleMultiplier * spc%nPoints) !nSampleMultiplier)
    do iTht = 1, size(theta_i)

       r_ref % P      = 0.0_fd
       r_ref % D      = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]

       iRand = 1
       rSeed = 1
       
       !$omp parallel default(none) &
       !$omp shared(fResults, iTht, r_ref, nSamples, mat, spc, theta_i, nOrders, nSampleMultiplier, timer, obsSolidAngle) &
       !$omp private(iRay,iSpec,iSampleSet,iLine,rayStat, r, rTable, muFluorLine, muTotal, l, fYield, rSeed, rState, rInfo, order) &
       !$omp firstprivate(iRand, dObs, cos_half_angle)

       !$rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)

       do iSampleSet = 1, nSampleMultiplier
          !$omp do schedule(guided, 100)
          do iSpec = 1, spc%nPoints
             do iRay = 1, nSamples

                if(iRand == 1) call dranduniform(SRTABLE, 0.0_fd, 1.0_fd, rState, rTable, rInfo)
                iRand      = mod(iRand, SRTABLE) + 1

                order      = 0
                r          = r_ref
                r % I      = spc%I(iSpec) / real(nSamples*nSampleMultiplier, fd) ! 1.0_fd / real(nSamples*nSampleMultiplier, fd)
                r % energy = spc%E(iSpec) !spc_getSample(spc, rTable(iRand))
                r % status = RAY_ACTIVE

                do while(r%status == RAY_ACTIVE .AND. order <= nOrders)

                   if(r%energy < 1e-2) then
                      print *, r%energy, order, iLine, fYield
                      stop
                   end if

                   call mat_evalMu(mat, r%energy, muFluorLineTotal = muFluorLine, muTotal = muTotal)

                   if(muFluorLine > MIN_MU) then

                      call drandexponential(1, 1.0_fd, rState, l, rInfo)
                      l = l / muTotal

                      r%P  = r%P + r%D * l(1)  

                      if(r%P(3) < 0.0) then    
                         call mat_evalFluorescence(mat, r%energy, r%energy, fYield, iLine)
                         r%D  = vec_cart_random_spherical_uniform()
                         r%I  = r%I * fYield * (muFluorLine/muTotal)

                         call mat_evalMu(mat, r%energy, muTotal = muTotal)
                         !$omp atomic
                         fResults(iTht, iLine) = fResults(iTht, iLine) + INV_FOUR_PI * r%I * exp(-muTotal * abs(r%P(3)/r%D(3)))
                      else
                         r%status = RAY_EXIT_U
                      end if

                   else
                      if(r%D(3) <= 0.0) then
                         r%status = RAY_EXIT_D
                      else
                         r%I      = r%I * exp(- abs(r%P(3)/r%D(3)) * muTotal)
                         r%status = RAY_EXIT_U
                      end if
                   end if

                   order = order + 1
                end do

                if(r%status == RAY_EXIT_U) then
                   if(dot_product(r%D, dObs) > cos_half_angle) then
                      !$omp atomic
                      fResults(iTht, iLine) = fResults(iTht, iLine) + r%I / obsSolidAngle
                   end if
                end if
             end do
             call utl_timerIncrease(timer)
          end do
          !$omp end do
       end do
       !$omp end parallel
    end do

  end subroutine infslab_method_montecarlo
  !******


 subroutine infslab_method_montecarlo_bruteforce(mat, theta_i, theta_e, nSamples)
    type(mat_material)     :: mat
    real(FD), dimension(:) :: theta_i
    real(fd)               :: theta_e
    integer(il)            :: nSamples

    integer,  parameter    :: SRTABLE       = 20000
    real(fd), parameter    :: MIN_MU        = 1e-20_fd
    real(fd), parameter    :: MIN_INTENSITY = 1e-23_fd

    ! Simulation (logical) variables
    !
    integer(il) :: iSampleset, iSpec, iRay, iRand, iLine, iTht
    integer     :: rayStat, order
    real(fd)    :: dObs(3), cos_half_angle

    type(utl_timer) :: timer

    ! Physical variables
    !
    real(fd)    :: muFluorLine, muTotal, fYield, l(1)
    type(ray)   :: r_ref, r

    ! Random number generator variables
    !
    real(fd)    :: rTable(SRTABLE)
    integer     :: rInfo, rSeed(1), rState(640)

    dObs = [sin(theta_e), 0.0_fd, cos(theta_e)]

    cos_half_angle = cos(obsHalfAngle)

    call utl_timer_init(timer, 5.0_fd, size(theta_i) * nSampleMultiplier * spc%nPoints) !nSampleMultiplier)
    do iTht = 1, size(theta_i)

       r_ref % P      = 0.0_fd
       r_ref % D      = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]

       iRand = 1
       rSeed = 1

       !$omp parallel default(none) &
       !$omp shared(fResults, iTht, r_ref, nSamples, mat, spc, theta_i, nOrders, nSampleMultiplier, timer, obsSolidAngle) &
       !$omp private(iRay,iSpec,iSampleSet,iLine,rayStat, r, rTable, muFluorLine, muTotal, l, fYield, rSeed, rState, rInfo, order) &
       !$omp firstprivate(iRand, dObs, cos_half_angle)

       !$ rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)

       do iSampleSet = 1, nSampleMultiplier
          !$omp do schedule(guided, 100)
          do iSpec = 1, spc%nPoints
             do iRay = 1, nSamples

                if(iRand == 1) call dranduniform(SRTABLE, 0.0_fd, 1.0_fd, rState, rTable, rInfo)
                iRand      = mod(iRand, SRTABLE) + 1

                order      = 0
                r          = r_ref
                r % I      = spc%I(iSpec) / real(nSamples*nSampleMultiplier, fd) ! 1.0_fd / real(nSamples*nSampleMultiplier, fd)
                r % energy = spc%E(iSpec) !spc_getSample(spc, rTable(iRand))
                r % status = RAY_ACTIVE

                do while(r%status == RAY_ACTIVE .AND. order <= nOrders)

                   if(r%energy < 1e-2) then
                      print *, r%energy, order, iLine, fYield
                      stop
                   end if

                   call mat_evalMu(mat, r%energy, muFluorLineTotal = muFluorLine, muTotal = muTotal)

                   if(muFluorLine > MIN_MU) then

                      call drandexponential(1, 1.0_fd, rState, l, rInfo)
                      l = l / muTotal

                      r%P  = r%P + r%D * l(1)  

                      if(r%P(3) < 0.0) then    
                         call mat_evalFluorescence(mat, r%energy, r%energy, fYield, iLine)
                         r%D  = vec_cart_random_spherical_uniform()
                         r%I  = r%I * fYield * (muFluorLine/muTotal)
                      else
                         r%status = RAY_EXIT_U
                      end if

                   else
                      if(r%D(3) <= 0.0) then
                         r%status = RAY_EXIT_D
                      else
                         r%I      = r%I * exp(- abs(r%P(3)/r%D(3)) * muTotal)
                         r%status = RAY_EXIT_U
                      end if
                   end if

                   order = order + 1
                end do

                if(r%status == RAY_EXIT_U) then
                   if(dot_product(r%D, dObs) > cos_half_angle) then
                      !$omp atomic
                      fResults(iTht, iLine) = fResults(iTht, iLine) + r%I / obsSolidAngle
                   end if
                end if
             end do
             call utl_timerIncrease(timer)
          end do
          !$omp end do
       end do
       !$omp end parallel
    end do

  end subroutine infslab_method_montecarlo_bruteforce

 
  subroutine saveData(fData, m, s, fName)
    real(fd), dimension(:,:)          :: fData
    type(mat_material), intent(IN)    :: m
    type(spc_spectrum), intent(IN)    :: s
    character (len = *)               :: fName

    integer          :: i

    integer          :: fileID
    integer          :: dimL, dimE, dimThetaIn, muAbsID, muCdfID, muExtID
    integer          :: dimSpcE,  dimSpcE2, fluorID, fluorAnID, spcID, spcCdfID, spcIcdfID

    Call gth_nfCheck( nf90_create(fName, NF90_CLOBBER, fileID) )

    Call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "Author",  "Hannu Parviainen" ) )
    Call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "Program", "xrfipsd" ) )
    Call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "Version", "0.9" ) )

    Call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "Elements", elements ) )
    Call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "LineEnergy", mat%lEnergy ) )
    Call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "EdgeEnergy", mat%eEnergy ) )
 
   if(sourceType == 'Spectrum') then
       Call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "SpectrumEnergy", s%E ) )
       Call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "TestSpectrum", testSpec ) )
    end if
       

    Call gth_nfCheck( nf90_def_dim(fileID, "Line", m%nLines, dimL) )
    Call gth_nfCheck( nf90_def_dim(fileID, "Energy",   m%nE,  dimE) )
    Call gth_nfCheck( nf90_def_dim(fileID, "ThetaIn",  nThetaIn,  dimThetaIn) )

    Call gth_nfCheck( nf90_def_dim(fileID, "SpectrumEnergy",   s%nPoints,  dimSpcE) )
    Call gth_nfCheck( nf90_def_dim(fileID, "SpectrumICDFEnergy",   5*s%nPoints,  dimSpcE2) )

    Call gth_nfCheck( nf90_def_var(fileID, "fluorescenceData", NF90_FLOAT, [dimThetaIn,dimL], fluorID) )
    Call gth_nfCheck( nf90_def_var(fileID, "fluorescenceData_analytic", NF90_FLOAT, [dimThetaIn,dimL], fluorAnID) )


    Call gth_nfCheck( nf90_def_var(fileID, "muAbs", NF90_FLOAT, [dimL,dimE], muAbsID) )
    Call gth_nfCheck( nf90_def_var(fileID, "muAbsCdf", NF90_FLOAT, [dimL,dimE], muCdfID) )
    Call gth_nfCheck( nf90_def_var(fileID, "muExt", NF90_FLOAT, [dimE], muExtID) )

    Call gth_nfCheck( nf90_def_var(fileID, "Spectrum", NF90_FLOAT, [dimSpcE], spcID) )

    if(sourceType == 'Spectrum') then
      Call gth_nfCheck( nf90_def_var(fileID, "SpectrumCdf", NF90_FLOAT, [dimSpcE], spcCdfID) )
      Call gth_nfCheck( nf90_def_var(fileID, "SpectrumCdfInv", NF90_FLOAT, [dimSpcE2], spcIcdfID) )
    end if


    call gth_nfCheck( nf90_enddef(fileID) )

    Call gth_nfCheck( nf90_put_var(fileID, fluorID, fData) )
    Call gth_nfCheck( nf90_put_var(fileID, fluorAnID, fResultsAn) )


    Call gth_nfCheck( nf90_put_var(fileID, muAbsID, m%muFluorLine) )
    Call gth_nfCheck( nf90_put_var(fileID, muCdfID, m%muFluorLineCDF) )
    Call gth_nfCheck( nf90_put_var(fileID, muExtID, m%muExtTotal) )
    Call gth_nfCheck( nf90_put_var(fileID, spcID,   s%I) )

    if(sourceType == 'Spectrum') then
       Call gth_nfCheck( nf90_put_var(fileID, spcCdfID,   s%cdf) )
       Call gth_nfCheck( nf90_put_var(fileID, spcIcdfID,   s%icdf) )
    end if

    Call gth_nfCheck( nf90_close(fileID) )

  end subroutine saveData

  subroutine testSpectrum(s, d, nSamples)
    type(spc_spectrum), intent(IN)    :: s
    real(fd), dimension(:)            :: d
    integer(il)                       :: nSamples

    integer(il)  :: res, r, i
    real(fd)     :: minE, maxE, dE, dEInv, E, x(1)
 
    d = 0.0_fd

    res   = size(d)
    minE  = s%E(1)
    maxE  = s%E(size(s%E))
    dE    = (maxE - minE) / real(res,fd)

    dEInv = 1._fd / dE

    do i = 1, nSamples
       call dst_generate_uniform_n(x)
       E = spc_getSample(s, x(1))
       r = int(floor( (E-minE) * dEInv ), il)+1

       d(r) = d(r) + 1.0_fd
    end do

    d = d / real(nSamples, fd)

  end subroutine testSpectrum

end program xrfipsd

