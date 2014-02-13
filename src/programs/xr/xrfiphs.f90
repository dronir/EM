!****h* XR/programs/xrfiphs
! NAME
!   xrfiphs
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

program xrfiphs
  use base
  use material
  use geometry
  use particle
  use hemisphere
  use random
  use bfdf
  !use random

  !$ use omp_lib

  implicit none

  character (len=FNAME_LENGTH) :: parFilename        = ""
  character (len=FNAME_LENGTH) :: hsFilename         = "xrfiphs.nc"

  character (len=100)          :: srcType            = "Spectrum"
  character (len=FNAME_LENGTH) :: srcSpectrumfile    = ""
  real(fd)                     :: srcLineEnergy      = 9e3

  integer                      :: nThetaIn           = 1
  character (len=10000)        :: thetaIn            = "45.0"
  real(fd), allocatable        :: thetaInTable(:)

  character (len=100)          :: method             = "MonteCarlo"

  character (len=256)          :: matName            = "Iron"
  integer                      :: nElements          = 2
  character (len=10000)        :: elements           = "Fe Ca"
  character (len=10000)        :: elemFracs          = "0.5 0.5"
  character(2), allocatable    :: elementTable(:)
  real(fd), allocatable        :: elemFracTable(:)
  real(fd)                     :: molarDensity       = 1.0_fd

  real(fd)                     :: energyMin          = 0.0_fd
  real(fd)                     :: energyMax          = 1.0e4_fd

  integer                      :: resThetaE          = 25
  integer(il)                  :: nSamples           = 50000
  integer                      :: nSampleMultiplier  = 1
  integer                      :: nOrders            = 1

  type(gth_hemisphere)         :: H
  type(gth_hemisphere)         :: hs_analytic

  type(spc_spectrum)           :: spc

  type(mat_material), target   :: mat

  real    :: sTime, cTime, eTime
  integer(il) :: i
  
  integer :: pid

  namelist /params/ &
       & hsFilename,&
       & srcType, srcSpectrumfile, srcLineEnergy, &
       & nThetaIn, thetaIn, &
       & nElements, elements, elemFracs, molarDensity,&
       & resThetaE, &
       & nSamples, nSampleMultiplier, &
       & matName, &
       & nOrders,   &
       & method,    &
       & energyMin, energyMax

  !pid = getpid()
  !print *,pid
  !call kill(pid, 9)

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
  thetaInTable = thetaInTable * CNV_DEGTORAD
  
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

  select case(srcType)
     case('Spectrum')
        call spc_spectrum_read(spc, srcSpectrumfile, mat%edgeEnergyMin, energyMax)
     case('Line')
        call spc_spectrum_initline(spc, srcLineEnergy, 1.0_fd)
     case default
        call utl_fatal_error("Error: Unknown spectrum type: " // srcType)
     end select

  call gth_hemisphere_init(H, resThetaE, nThetaIn, mat%nLines, GTH_TYPE_HS)
  write(H%name,'("Hemisphere")')

  call gth_hemisphere_init(hs_analytic, resThetaE, nThetaIn, mat%nLines, GTH_TYPE_HS)
  write(hs_analytic%name,'("Hemisphere_analytic")')


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! RUN SIMULATION
  !!

  call cpu_time(sTime) 
  select case(method)
     case('MonteCarlo_no_peeling')
        call infslab_method_montecarlo_pure(thetaInTAble, H, mat, nSamples)
     case('MonteCarlo')
        call infslab_method_montecarlo(thetaInTAble, H, mat, nSamples)
     case('FirstOrder')
        call infslab_method_firstorder(mat, H, thetaInTAble, nSamples)
     case('Analytic')
        call infslab_method_analytic(mat, spc, H, thetaInTAble)
     case default
        call utl_fatal_error("Error: Unknown simulation method: " // method)
     end select
  call cpu_time(cTime)

  write(*,'("Time taken: " ,(F8.3))') cTime - sTime


  call cpu_time(sTime) 
  call infslab_method_analytic(mat, spc, hs_analytic, thetaInTable)
  call cpu_time(cTime)
  write(*,'("Time taken: " ,(F8.3))') cTime - sTime


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! WRITE FILE
  !!

  call saveHemisphere(H,  mat, hsFilename)

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! UTILITY FUNCTIONS AND SUBROUTINES
  !!
  !!

  !! First order analytic
  !!


 subroutine infslab_method_analytic(mat, spc, hs, theta_i)
    type(mat_material)       :: mat
    type(spc_spectrum)       :: spc
    type(gth_hemisphere)     :: hs

    real(fd), dimension(:)   :: theta_i

    real(fd)                 :: energyIn, intensityIn
    real(fd), dimension(mat%nLines) :: fluorescence

    real(fd)    :: muInTotal, muEmTotal, muInFluor
    real(fd)    :: A, B, fOut, xEn, c_thtI

    integer :: iTht, iThtE, iLine, iSpec, iEn

 
    do iTht = 1, size(theta_i)
       do iThtE = 1, hs%resTheta
          fluorescence = bfdf_Parviainen(mat, spc, theta_i(iTht), hs%mTheta(iThtE))
          
          do iLine = 1, mat%nLines
             call gth_hs_addToRow(hs, iThtE, fluorescence(iLine), iTht, iLine)
          end do
       end do
    end do

  end subroutine infslab_method_analytic



  subroutine infslab_method_analytic_o(mat, s, hs, theta_i)
    type(mat_material)       :: mat
    type(spc_spectrum)       :: s
    type(gth_hemisphere)     :: hs

    real(fd), dimension(:)   :: theta_i

    real(fd)                 :: e, spcNorm

    real(fd)    :: muInTot, muInExt, muInLine(mat%nLines), muEmLine(mat%nLines), fYield(mat%nLines)
    real(fd)    :: A(mat%nLines), B(mat%nLines), fOut(mat%nLines), xEn

    integer :: iTht, iThtE, iLine, iSpec, iEn

 
    do iSpec = 1, s%nPoints

       e = s%E(iSpec)

       spcNorm = s%I(iSpec) / sum(s%I)

       muInTot  = mat_evalMuAbs(mat, e, iEn, xEn)
       muInLine = ((1.0_fd - xEn) * mat%muFluorLine(:,iEn) + xEn * mat%muFluorLine(:,iEn+1))

       muInExt  = mat_evalMuExt(mat, e) - muInTot

       do iLine = 1, mat%nLines
          if(mat%lEnergy(iLine) > 0.0_fd) then
             muEmLine(iLine) = mat_evalMuExt(mat, mat%lEnergy(iLine))
             fYield(iLine)   = mat%fYield((iLine-1)/2 + 1)
          end if
       end do

       A = spcNorm * muInLine / muInTot * fYield
       B = muEmLine / muInTot

       do iTht = 1, size(theta_i)
          do iThtE = 1, hs%resTheta

             fOut = INV_FOUR_PI * A * cos(hs%mTheta(iThtE)) / &
                  & ( cos(hs%mTheta(iThtE)) + B * cos(theta_i(iTht)))

             !$omp critical
             do iLine = 1, mat%nLines
                call gth_hs_addToRow(hs, iThtE, fOut(iLine), iTht, iLine)
             end do
             !$omp end critical

          end do
       end do
    end do

  end subroutine infslab_method_analytic_o


  !! First order Monte Carlo approximation
  !!
  subroutine infslab_method_firstorder(mat, hs, theta_i, nSamples)
    type(mat_material)     :: mat
    type(gth_hemisphere)   :: hs

    real(FD), dimension(:) :: theta_i
    integer(il)            :: nSamples

    ! Simulation (logical) variables
    !
    integer(il) :: iSpc, iRay, iRand, iLine, iTht, iThtE
    real(fd)    :: fluorescence_ratio, spcNorm, x

    type(utl_timer) :: timer

    ! Physical variables
    !
    real(fd), dimension(mat%nLines) :: muLineExt, fluorescence_intensity, fYield, line_ratio
    real(fd)                        :: muFluor, muTotal,  l(1)
    type(ray)   :: r_ref, r

    ! Random number generator variables
    !
    integer     :: rInfo, rSeed(1), rState(640)

    spcNorm   = 1.0_fd / sum(spc%I)
    muLineExt = 0.0_fd

    do iLine = 1, mat%nLines
       if(mat%lEnergy(iLine) > 0.0_fd) then
          muLineExt(iLine) = mat_evalMuExt(mat, mat%lEnergy(iLine))
          fYield(iLine)    = mat%fYield((iLine-1)/2 + 1)
       end if
    end do

    call utl_timer_init(timer, 5.0_fd, size(theta_i) * nSamples * spc%nPoints)

    do iTht = 1, size(theta_i)

       r_ref % P = 0.0_fd
       r_ref % D = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]
       r_ref % I = 1.0_fd / real(nSamples, fd)

       iRand   = 1
       rSeed   = 1

       !$omp parallel default(none)&
       !$omp shared(r_ref, nSamples, hs, mat, spc, theta_i, iTht, spcNorm, timer) &
       !$omp private(iSpc, iThtE, iRay, iLine, r, fluorescence_ratio, fluorescence_intensity, line_ratio, muFluor, muTotal, l, rSeed, rState, rInfo) &
       !$omp firstprivate(iRand, muLineExt, fYield)

       !$ rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)

       do iSpc = 1, spc%nPoints
          !$omp do schedule(dynamic) 
          do iRay = 1, nSamples

             !call dranduniform     (1, 0.0_fd, 1.0_fd, rState, x, rInfo)
             call drandexponential(1, 1.0_fd, rState, l, rInfo)

             r          = r_ref
             r % energy = spc%E(iSpc)
             r % I      = spc%I(iSpc)/real(nSamples,fd)*spcNorm

             call mat_evalMu(mat, r%energy, muFluorLineTotal = muFluor, muTotal = muTotal)

             l    = l / muTotal
             r%P  = r%P + r%D * l(1)

             fluorescence_ratio     = muFluor / muTotal
             line_ratio             = mat_evalFluorLineFractions(mat, r%energy)

             do iThtE = 1, hs%resTheta
                fluorescence_intensity = &
                     & INV_FOUR_PI * r%I * fYield * fluorescence_ratio * line_ratio &
                     & * exp(-muLineExt * abs(r%P(3)) / cos(hs%mTheta(iThtE)))

                !$omp critical
                do iLine = 1, mat%nLines
                   call gth_hs_addToRow(hs, iThtE, fluorescence_intensity(iLine), iTht, iLine)
                end do
                !$omp end critical
             end do

             call utl_timerIncrease(timer)

          end do
          !$omp end do
       end do
       !$omp end parallel

    end do
  end subroutine infslab_method_firstorder

  !! First order Monte Carlo approximation
  !!
  subroutine infslab_method_firstorder_b(mat, hs, theta_i, nSamples)
    type(mat_material)     :: mat
    type(gth_hemisphere)   :: hs

    real(FD), dimension(:) :: theta_i
    integer(il)            :: nSamples

    integer, parameter :: SRTABLE  = 20000

    ! Simulation (logical) variables
    !
    integer(il) :: iRay, iRand, iLine, iTht, iThtE, iEn
    real(fd)    :: xEn

    type(utl_timer) :: timer

    ! Physical variables
    !
    real(fd)    :: muFluorLine, muExt, muLineExt(mat%nLines), fOut(mat%nLines), fYield(mat%nLines), fInt(mat%nLines), l(1)
    type(ray)   :: r_ref, r

    ! Random number generator variables
    !
    real(fd)    :: rTable(SRTABLE)
    integer     :: rInfo, rSeed(1), rState(640)

    muLineExt = 0.0_fd

    do iLine = 1, mat%nLines
       if(mat%lEnergy(iLine) > 0.0_fd) then
          muLineExt(iLine) = mat_evalMuExt(mat, mat%lEnergy(iLine))
          fYield(iLine)    = mat%fYield((iLine-1)/2 + 1)
       end if
    end do

    call utl_timer_init(timer, 5.0_fd, size(theta_i) * nSamples)

    do iTht = 1, size(theta_i)

       r_ref % P = 0.0_fd
       r_ref % D = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]
       r_ref % I = 1.0_fd / real(nSamples, fd)

       iRand   = 1
       rSeed   = 1

       !$omp parallel &
       !$omp shared(r_ref, nSamples, mat, theta_i) &
       !$omp private(iRay, iLine, r, muFluorLine, muExt, fOut, fInt, l, rTable, rSeed, rState, rInfo, iEn, xEn, iThtE) &
       !$omp firstprivate(iRand, muLineExt, fYield)

       !$ rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)

       !$omp do schedule(dynamic) 
       do iRay = 1, nSamples

          if(iRand == 1) call dranduniform(SRTABLE, 0.0_fd, 1.0_fd, rState, rTable, rInfo)
          iRand = mod(iRand, SRTABLE) + 1

          r          = r_ref
          r % energy = spc_getSample(spc, rTable(iRand))

          muFluorLine  = mat_evalMuAbs(mat, r%energy, iEn, xEn)
          muExt        = mat_evalMuExt(mat, r%energy)

          call drandexponential(1, 1.0_fd, rState, l, rInfo)
          l = l / muFluorLine

          r%P  = r%P + r%D * l(1)

          fInt = ((1.0_fd - xEn) * mat%muFluorLine(:,iEn) + xEn * mat%muFluorLine(:,iEn+1)) / muFluorLine

          do iThtE = 1, hs%resTheta
             fOut = INV_FOUR_PI * r%I * fYield * fInt * exp(-muLineExt * abs(r%P(3) / cos(hs%mTheta(iThtE))))

             do iLine = 1, mat%nLines
                call gth_hs_addToRow(hs, iThtE, fOut(iLine), iTht, iLine)
             end do
          end do

          call utl_timerIncrease(timer)

       end do
       !$omp end do
       !$omp end parallel

    end do
  end subroutine infslab_method_firstorder_b

  subroutine infslab_method_montecarlo(theta_i, hs, mat, nSamples)
    real(FD), dimension(:) :: theta_i
    type(gth_hemisphere)   :: hs
    type(mat_material)     :: mat
    integer(il)            :: nSamples

    integer, parameter :: RAY_OK   = 0
    integer, parameter :: RAY_EXIT = 1
    integer, parameter :: RAY_LOW  = 2
    integer, parameter :: SRTABLE  = 20000


    ! Simulation (logical) variables
    !
    integer(il) :: iSpec, iRay, iRand, iLine, iTht
    integer     :: rayStat, order

    type(utl_timer) :: timer

    ! Physical variables
    !
    real(fd)    :: muFluor, muTotal, muExt, fYield, l(1)
    type(ray)   :: r_ref, r

    ! Random number generator variables
    !
    real(fd)    :: rTable(SRTABLE)
    integer     :: rInfo, rSeed(1), rState(640)

    call utl_timer_init(timer, 10.0_fd, size(theta_i) * nSamples*spc%nPoints)
    do iTht = 1, size(theta_i)

       r_ref % P      = 0.0_fd
       r_ref % D      = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]

       iRand = 1
       rSeed = 1

       !$omp parallel &
       !$omp shared(r_ref, nSamples, mat, hs, theta_i, nOrders) &
       !$omp private(iRay, iSpec, rayStat, r, muFluor, muTotal, muExt, l, fYield, rSeed, rState, rInfo, order) &
       !$omp firstprivate(iRand)

       !$ rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)
       
       !$omp do schedule(dynamic)
       do iSpec = 1, spc%nPoints
          do iRay = 1, nSamples

             !if(iRand == 1) call dranduniform(SRTABLE, 0.0_fd, 1.0_fd, rState, rTable, rInfo)
             !iRand = mod(iRand, SRTABLE) + 1

             r          = r_ref
             r % I      = spc%I(iSpec) / real(nSamples, fd)
             r % energy = spc%E(iSpec)

             !r % I      = 1.0_fd / real(nSamples, fd)
             !r % energy = spc_getSample(spc, rTable(iRand))

             rayStat    = RAY_OK
             order      = 1

             do while(rayStat == RAY_OK)

                call mat_evalMu(mat, r%energy, muFluorLineTotal = muFluor, muTotal = muTotal)

                if(muFluor > 1e-22_fd) then

                   call drandexponential(1, 1.0_fd, rState, l, rInfo)

                   l    = l / muTotal
                   r%P  = r%P + r%D * l(1)

                   !! Test whether we are still inside the medium.
                   !!
                   if(r%P(3) > 0.0_fd) then
                      rayStat = RAY_EXIT
                   else if(r%I < 1.0e-23_fd) then
                      rayStat = RAY_LOW
                   end if

                   !! Fluorescence
                   !!
                   if(rayStat == RAY_OK) then

                      call mat_evalFluorescence(mat, r%energy, r%energy, fYield, iLine)
                      r%I  = r%I * fYield  * (muFluor/muTotal)

                      if (order <= nOrders) then

                         !! Peeling
                         !!
                         call mat_evalMu(mat, r%energy, muTotal = muExt)

                         r%D    = vec_cart_random_hemispherical_uniform()

                         !$omp critical 
                         call gth_addDataCar(hs, r%D, INV_FOUR_PI * r%I * exp(-muExt * abs(r%P(3)/r%D(3))), iTht, iLine)
                         !$omp end critical

                         !! Fluorescence
                         !!
                         r%D  = vec_cart_random_spherical_uniform()

                      else
                         rayStat = RAY_LOW
                      end if

                   end if

                else
                   if(r%D(3) <= 0.0_fd) then
                      rayStat = RAY_LOW
                   else
                      rayStat = RAY_EXIT
                   end if
                end if

                order = order + 1

             end do

             !if(rayStat == RAY_EXIT) then
             !   !$omp critical 
             !   call gth_addDataCar(hs, r%D, r%I, iTht, iLine)
             !   !$omp end critical
             !end if

             call utl_timerIncrease(timer)
          end do
       end do
       !$omp end do
       !$omp end parallel

    end do
  end subroutine infslab_method_montecarlo

  subroutine infslab_method_montecarlo_pure(theta_i, hs, mat, nSamples)
    real(FD), dimension(:) :: theta_i
    type(gth_hemisphere)   :: hs
    type(mat_material)     :: mat
    integer(il)            :: nSamples

    integer,  parameter :: SRTABLE       = 20000
    real(fd), parameter :: MIN_MU        = 1e-20_fd
    real(fd), parameter :: MIN_INTENSITY = 1e-23_fd

    ! Simulation (logical) variables
    !
    integer(il) :: iSampleset, iRay, iRand, iLine, iTht
    integer     :: rayStat, order

    type(utl_timer) :: timer

    ! Physical variables
    !
    real(fd)    :: muFluorLine, muTotal, fYield, l(1)
    type(ray)   :: r_ref, r

    ! Random number generator variables
    !
    real(fd)    :: rTable(SRTABLE)
    integer     :: rInfo, rSeed(1), rState(640)

    call utl_timer_init(timer, 5.0_fd, size(theta_i) * nSamples * nSampleMultiplier)
    do iTht = 1, size(theta_i)

       r_ref % P      = 0.0_fd
       r_ref % D      = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]

       iRand = 1
       rSeed = 1

       !$omp parallel &
       !$omp shared(r_ref, nSamples, mat, hs, theta_i, nOrders, nSampleMultiplier) &
       !$omp private(iRay,iSampleSet, rayStat, r, muFluorLine, muTotal, l, fYield, rSeed, rState, rInfo, order) &
       !$omp firstprivate(iRand)

       !$ rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)

       do iSampleSet = 1, nSampleMultiplier
          !$omp do schedule(guided, 100)
          do iRay = 1, nSamples

             if(iRand == 1) call dranduniform(SRTABLE, 0.0_fd, 1.0_fd, rState, rTable, rInfo)
             iRand      = mod(iRand, SRTABLE) + 1

             order      = 0
             r          = r_ref
             r % I      = 1.0_fd / real(nSamples*nSampleMultiplier, fd)
             r % energy = spc_getSample(spc, rTable(iRand))
             r % status = RAY_ACTIVE

             do while(r%status == RAY_ACTIVE .AND. order <= nOrders)

                muFluorLine  = mat_evalMuAbs(mat, r%energy)
                muTotal      = mat_evalMuExt(mat, r%energy)

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
                !$omp critical
                call gth_hs_addToRowCar(hs, r%D, r%I, iTht, iLine)
                !$omp end critical
             end if

             call utl_timerIncrease(timer)
          end do
          !$omp end do
       end do
       !$omp end parallel
    end do

    call gth_hs_normalizeByTheta(hs)

  end subroutine infslab_method_montecarlo_pure



 subroutine saveHemisphere(h, mat, fName)
    type(gth_hemisphere), INTENT(IN)  :: h
    type(mat_material), intent(IN)    :: mat
    character (len = *)               :: fName

    type(gth_hsFile)   :: hf
    integer            :: i

    integer            :: dimMatMuEnergy, dimFlrLine, dimSpecEnergy, dimSpecICDF
    integer            :: idMatMuIon, idMatMuFluor, idMatMuFluorCdf, idMatMuRay, idMatMuExt, idMatMuEnergy
    integer            :: idSpecEnergy, idSpecIntensity, idSpecCdf, idSpecIcdf
    integer            :: idAnalytic

    integer, parameter :: RTP = NF90_DOUBLE
    integer, parameter :: RTS = NF90_FLOAT


    call gth_hsFileOpen(hf, fName, "w")
    call gth_hsWriteHeader(hf, h, "Hannu Parviainen", "xrfiphs", "1.0")

    Call gth_nfCheck( nf90_def_var(hf%fileID, "Hemisphere_analytic", NF90_DOUBLE, hf%dimHs, idAnalytic) )
    Call gth_nfCheck( nf90_put_att(hf%fileID,  idAnalytic, "long_name",  "Analytic approximation" ) )

    !! Dimensions
    !!
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Fluorescence_line",     mat%nLines,           dimFlrLine     ) )
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Material_energy",       mat%nE,               dimMatMuEnergy ) )

    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Spectrum_energy",       spc%nPoints,          dimSpecEnergy  ) )
    if(srcType == 'Spectrum') then
       Call gth_nfCheck( nf90_def_dim(hf%fileID, "Spectrum_ICDF_energy",  spc%nPoints*5,      dimSpecICDF    ) )
    end if

    !! Variables
    !!
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Photoionization_coefficient",   RTP, [dimMatMuEnergy],            idMatMuIon) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Fluorescence_line_coefficient", RTP, [dimFlrLine,dimMatMuEnergy], idMatMuFluor) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Fluorescence_line_cdf",         RTP, [dimFlrLine,dimMatMuEnergy], idMatMuFluorCdf))
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Rayleigh_coefficient",          RTP, [dimMatMuEnergy],            idMatMuRay) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Extinction_coefficient",        RTP, [dimMatMuEnergy],            idMatMuExt) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Material_energy",               RTP, [dimMatMuEnergy],            idMatMuEnergy) )

    Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_energy",               RTP, [dimSpecEnergy], idSpecEnergy) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_intensity",            RTP, [dimSpecEnergy], idSpecIntensity) )
    if(srcType == 'Spectrum') then
       Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_cdf",               RTP, [dimSpecEnergy], idSpecCdf) )
       Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_inverse_cdf",       RTP, [dimSpecICDF],   idSpecIcdf) )
    end if

    !! Global attributes
    !!
    call gth_hsAddAttS(hf, NF90_GLOBAL, "Elements",           elements)
    call gth_hsAddAttS(hf, NF90_GLOBAL, "Spectrum_type",      srcType)
    call gth_hsAddAttS(hf, NF90_GLOBAL, "Simulation_method",  method)
    call gth_hsAddAttF(hf, NF90_GLOBAL, "Theta_in",           thetaInTable)

    call gth_hsAddAttF(hf, NF90_GLOBAL, "Fluorescence_line_energy", mat%lEnergy)
    call gth_hsAddAttF(hf, NF90_GLOBAL, "Absorbtion_edge_energy",   mat%eEnergy)

    call gth_hsSaveData(h, hf)
    Call gth_nfCheck( nf90_put_var(hf%fileID, idAnalytic, hs_analytic%data) )

    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuIon,       mat%muPhotoIon) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuFluor,     mat%muFluorLine) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuFluorCdf,  mat%muFluorLineCDF) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuRay,       mat%muRayleigh) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuExt,       mat%muExtTotal) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuEnergy,    mat%E) )

    Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecEnergy,     spc%E) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecIntensity,  spc%I) )
    if(srcType == 'Spectrum') then
       Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecCdf,     spc%cdf) )
       Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecIcdf,    spc%icdf) )
    end if

    call gth_hsFileClose(hf)

  end subroutine saveHemisphere


end program xrfiphs

