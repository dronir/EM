!!--- BEGIN GPL --- 
!!
!! xrfpm -  a program to compute soft X-ray fluorescence of particulate media.
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

!> Computes light scatterig from particulate media.
!!
!! \author Hannu Parviainen
!! \version 1.0
!! \date 15.11.2007

program vScatter
  use iso_c_binding
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
  character (len=FNAME_LENGTH) :: hsCarPlotFilename  = "vScatter_car.eps"
  character (len=FNAME_LENGTH) :: hsSphPlotFilename  = "vScatter_sph.eps"
  character (len=FNAME_LENGTH) :: hsPlotFilename     = ""
  character (len=FNAME_LENGTH) :: hsFilename         = "hemisphere.nc"
  character (len=FNAME_LENGTH) :: mediumFileName     = "medium.nc"
  character (len=100)          :: sDistribution      = "constant"
  character (len=100)          :: brdfType           = "shadowing"
  character (len=100)          :: brdfPhaseFunction  = "constant"

  integer                      :: gridResHorizontal  = 200
  integer                      :: gridResVertical    = 10

  real(fd)                     :: rho                = 0.50_fd
  real(fd)                     :: rhoAllowedError    = 0.05_fd

  integer                      :: resThetaE          = 45
  integer                      :: sampleSeed         = 0


  integer                      :: nOrders            = 1
  character(len=100)           :: nSamplesPerOrder   = "1"
  integer, allocatable         :: nSamplesPerOrderTable(:)

  integer                      :: nThetaIn           = 1
  character (len=10000)        :: thetaIn            = "45.0"
  real(fd), allocatable        :: thetaInTable(:)

  real(fd)                     :: mu                 = 1.0
  real(fd)                     :: w                  = 0.1

  integer                      :: nPfParams          = 1
  character (len=10000)        :: pfParams           = "0.0"
  real(fd)                     :: pfParamsTable(10)

  logical                      :: rf_applyFields    = .false.
  integer                      :: rf_nFieldsPerMed  = 1
  character(len=100)           :: rf_spectrumType   = "Gaussian"
  real(fd)                     :: rf_P              = 0.5_fd
  real(fd)                     :: rf_std            = 1.0_fd

  integer                      :: mediumMapRes        = 256
  integer                      :: mediumDensitymapRes = 300

  integer                      :: nThreads            = 1

  type(gth_hemisphere)  :: H
  type(med_medium)      :: M
  type(med_mediumFile)  :: Mf
  type(rf_randomField)  :: rndField
  type(utl_timer)       :: timer

  real, dimension(:,:,:), allocatable :: mediumHeightMap
  real(fd), dimension(:,:,:), allocatable :: mediumDensityStructure

  real(fd), external, pointer    :: Pf

  real    :: sTime, cTime, eTime
  integer :: i, iField

  namelist /params/ &
       & resThetaE, gridResHorizontal, gridResVertical, &
       & hsCarPlotFilename, hsSphPlotFilename, hsFilename, mediumFilename, &
       & sDistribution, rho, rhoAllowedError,  sampleSeed, &
       & nThetaIn, thetaIn, &
       & nOrders, nSamplesPerOrder, &
       & mu, w, nPfParams, pfParams, &
       & brdfType, brdfPhaseFunction, &
       & rf_applyFields, rf_nFieldsPerMed, rf_spectrumType, rf_P, rf_std, &
       & nThreads, &
       & mediumMapRes, mediumDensitymapRes


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

  !$ call omp_set_num_threads(nThreads)

  !!- Initialize random number generator
  !!
  call rnd_init(1, nThreads, .true.)

  !! Read the angles of incidence into a table.
  !!
  allocate(thetaInTable(nThetaIn))
  read(thetaIn,*) thetaInTable
 
  thetaInTable = CNV_DEGTORAD * thetaInTable

  !! Read the number of samples for each scattering order into a table.
  !!
  allocate(nSamplesPerOrderTable(nOrders))
  read(nSamplesPerOrder,*) nSamplesPerOrderTable

  !! Read the phase function parameters into a table.
  !!
  !allocate(pfParamsTable(nPfParams))
  read(PfParams,*) pfParamsTable(1:nPfParams)


  !! Initialize hemisphere
  !!

  !if(brdfType == "MonteCarlo") then
  !   call gth_hemisphere_init(H, resThetaE, nThetaIn, 1, GTH_TYPE_QS)
  !else
     call gth_hemisphere_init(H, resThetaE, nThetaIn, nOrders, GTH_TYPE_QS)
  !end if

  write(H%name,'("Hemisphere")')

  !! Select the media for the simulation from a file.
  !!
  call med_mediumFileOpen(Mf, mediumFilename)
  !!call med_mediumFileSelectMedia(Mf, rho, rhoAllowedError, sDistribution)
  call med_mediumFileSelectMedia(Mf)

  if(rf_applyFields) then

     call rf_init(rndField, 1000, 5.0_fd)

     select case(rf_spectrumType)        
     case('Gaussian')
        call rf_generateSpectrum (rndField, RF_SPECTRUM_GAUSSIAN, rf_P, rf_std)
     case('fBm')
        call rf_generateSpectrum (rndField, RF_SPECTRUM_FBM, rf_P, rf_std)
     end select
  end if

  call utl_timer_init(timer, 5.0_fd, Mf%nSelectedMedia * rf_nFieldsPerMed * H%nCells)
     
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! RUN SIMULATION
  !!

  allocate(mediumHeightMap(MF%nSelectedMedia, mediumMapRes, mediumMapRes))
  allocate(mediumDensityStructure(MF%nSelectedMedia, mediumDensitymapRes, 2))

  call utl_message("BRDF type used: " // brdfType)

  do i = 1, MF%nSelectedMedia
     do iField = 1, rf_nFieldsPerMed
        call med_mediumFileRead(M, Mf, Mf%varSelection(i))

        call med_gridAssign(M, gridResHorizontal, gridResVertical)
        call med_gridFit(M)
        call med_updateStatistics(M)


        if(rf_applyFields) then

           !!- Smooth the medium surface: 
           !!  Remove a constant percentage of the medium height from the top. 
           !!
           if(rf_spectrumType == 'constant') then
              rndField%field = rf_P * M%height

           !!- Apply macroscale surface roughness.
           !!
           else
              call rf_generateField(rndField, iField)
              call med_maskHeight(m, real(rndField%field, fs))
           end if
        end if

        mediumHeightMap(i,:,:) = cnt_grid3D_cmpHeightMap(M%grid, mediumMapRes)
        call med_computePorosityStructure(M, 5000, mediumDensityStructure(i,:,:))

        !call med_gridFit(M)
        !call med_updateStatistics(M)

        !call utl_message("Beginning trace.")
 
        call cpu_time(sTime) 
        select case(brdfType)
        case("shadowing")
           call sampleHemisphere(brdf_shadowing, nSamplesPerOrderTable(1), phase_function_constant, pfParamsTable)
        case("Lambert")
           call sampleHemisphere(brdf_Lambert, nSamplesPerOrderTable(1), phase_function_constant, pfParamsTable)
        case("LommelSeeliger")
           select case(brdfPhaseFunction)
           case("constant")
              call sampleHemisphere(brdf_LommelSeeliger, nSamplesPerOrderTable(1), phase_function_constant, pfParamsTable)
           case("HG1")
              call sampleHemisphere(brdf_LommelSeeliger, nSamplesPerOrderTable(1), phase_function_HG1, pfParamsTable)
           case("HG2")
              call sampleHemisphere(brdf_LommelSeeliger, nSamplesPerOrderTable(1), phase_function_HG2, pfParamsTable)
           end select
       case("RadTransfer")
           call sampleHemisphere(brdf_RadTrans, nSamplesPerOrderTable(1), phase_function_constant, pfParamsTable)
        case("MonteCarlo")
           call sampleHemisphere_MC(nSamplesPerOrderTable(1))
           !call simulation_method_montecarlo(M, thetaInTable, H, mu, w, nSamplesPerOrderTable(1))
        case default
           call utl_fatal_error("Unsupported BRDF type.")
        end select

        call cpu_time(cTime)
     end do
  end do

  H % data = H % data / (real(MF%nSelectedMedia,fd)*(real(rf_nFieldsPerMed,fd)))

  write(*,'("Time taken: " ,(F8.3))') cTime - sTime

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! WRITE FILES AND CLEAN UP
  !!
  call saveHemisphere(H, hsFilename)
  call med_mediumFileClose(Mf)


contains
  !! TODO: Include single scattering phase function to computations!
  !!
  subroutine sampleHemisphere(f, nSamples, Pf, Pp)

    real(fd), external :: f
    integer            :: nSamples
    real(fd), external :: Pf
    real(fd)           :: Pp(10)

    type(ray) :: rC
    type(intersection_geometry) :: iSect

    real(fd), dimension(2,nSamples) :: samples
    real(fd), dimension(3)          :: pSurface(3)!, pElement(3)
    real(fd), dimension(3,nThetaIn) :: D
    real(fd)                        :: dz, thtIn
    integer  :: i, j, k, iTheta
    logical  :: pFound, pLit

    real     :: tstRnd

    D(1,:) = sin(thetaInTable)
    D(2,:) = 0.0_fd
    D(3,:) = cos(thetaInTable)

    call smpl_griddedSamples2D(samples, nSamples)

    samples = (samples * M%width - M%hWidth) * 0.5_fd

    ! $omp parallel default(none)                 &
    ! $omp shared(M, H, D, samples, f, w, Pf, Pp, nOrders, nSamplesPerOrderTable, nThetaIn, thetaInTable, timer) &
    ! $omp private(rC, dz, thtIn, pSurface, i, j, k, iTheta, pFound, pLit, iSect, tstRnd)
 
    dz = M%grid%height - M%hMean - TRACE_EPS

    pSurface(3) = M%hMean

    call ray_init(rC, RAY_TYPE_CAMERA)

    ! $omp do schedule(dynamic) 
    do j = 1, H % resTheta
       do k = 1, H % resPhi(j)

          ! $omp do schedule(dynamic) 
          do i= 1, nSamplesPerOrderTable(1)
             
             !! Select the sample point of the incident camera ray from
             !! the mean medium surface.
             !!
             pSurface(1:2) = samples(:,i)

             pFound = .false.
             do while (.not. pFound)
                
                !! Find the intersection of the camera ray and the top of 
                !! the bounding box of the periodic medium.
                !!
                rC%D        = gth_cellRandomSampleCar(H, j, k)

                if(h%type == GTH_TYPE_QS) then
                   call RANDOM_NUMBER(tstRnd)
                   if(tstRnd > 0.5) rC%D(2) = -rC%D(2)
                end if

                if(rC%D(3) < 1e-2_fd) then
                   rC%D(3) = 1e-2
                   call vec_normalize(rC%D)
                end if

                rC % P(1)     = pSurface(1) + dz * (rC%D(1) / rC%D(3))
                rC % P(2)     = pSurface(2) + dz * (rC%D(2) / rC%D(3))
                rC % P(3)     = M%grid%height - TRACE_EPS

                rC % P(1)     = modulo(rC%P(1)+M%hWidth, M%width) - M%hWidth
                rC % P(2)     = modulo(rC%P(2)+M%hWidth, M%width) - M%hWidth

                rC % rayID    = rC%rayID + RAY_ID_INCR

                !! Negate the ray.
                !!
                rC % D     = - rC%D

                !! Find the true intersection point of the camera ray and medium.
                !!
                pFound        = trc_traceNearest(M%grid, rC, iSect)

                !! If an intersection is found, check if the point is shadowed.
                !!
                if(pFound) then
                   do iTheta = 1, nThetaIn
                      call trc_gatherRadiance(M%grid, rC%D, D(:,iTheta), iSect%P1 + TRACE_EPS * iSect%N, &
                           & iSect%N, 1.0_fd / real(nSamplesPerOrderTable(1), fd), &
                           & nSamplesPerOrderTable, nOrders, 1, H % data(iTheta, h%cIdx(j)+k-1,:), w, f, Pf, Pp)
                   end do
                else
                   call rnd_generate_uniform(0.0_fd, 1.0_fd, pSurface(1:2))
                   pSurface(1:2) = (pSurface(1:2) * M%width - M%hWidth) * 0.5_fd

                end if

             end do

          end do
          ! $omp end do
          call utl_timerIncrease(timer)  

       end do
    end do
    ! $omp end do
    ! $omp end parallel

  end subroutine sampleHemisphere

  subroutine sampleHemisphere_MC(nSamples)
    integer   :: nSamples

    integer, parameter :: RAY_OK   = 0
    integer, parameter :: RAY_EXIT = 1
    integer, parameter :: RAY_LOW  = 2

    type(ray) :: rC, rPeel
    type(intersection_geometry) :: iSect

    real(fd), dimension(2,nSamples) :: samples
    real(fd), dimension(3)          :: pSurface(3)
    real(fd), dimension(3,nThetaIn) :: D
    real(fd)                        :: dz, thtIn
    integer  :: i, j, k, iTheta, rayStat, order
    logical  :: pFound, pLit

    real(fd) :: l(1), optDepth
    real     :: tstRnd

    integer  :: rInfo, rSeed(1), rState(640)

    D(1,:) = -sin(thetaInTable)
    D(2,:) = 0.0_fd
    D(3,:) = cos(thetaInTable)

    call smpl_griddedSamples2D(samples, nSamples)

    samples = (samples * M%width - M%hWidth) * 0.5_fd

    !$omp parallel default(none)                  &
    !$omp shared(M, H, D, nSamples, samples, nOrders, nSamplesPerOrderTable, nThetaIn, thetaInTable, w, mu, timer) &
    !$omp private(rC, dz, thtIn, pSurface, i, j, k, iTheta, tstRnd, optDepth, l, rSeed, rState, rInfo, order, rPeel, rayStat)
 
    dz = M%grid%height - M%hMean - TRACE_EPS
    pSurface(3) = M%hMean
    call ray_init(rC, RAY_TYPE_CAMERA)
    !$ rSeed(1) = omp_get_thread_num()+1
    !call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)

    !$omp do schedule(dynamic) 
    do j = 1, H % resTheta
       do k = 1, H % resPhi(j)

          ! $omp do schedule(dynamic) 
          do i= 1, nSamplesPerOrderTable(1)
             
             !! Select the sample point of the incident camera ray from
             !! the mean medium surface.
             !!
             pSurface(1:2) = samples(:,i)

                
             !! Find the intersection of the camera ray and the top of 
             !! the bounding box of the periodic medium.
             !!
             
             rC % D        = gth_cellRandomSampleCar(H, j, k)

             if(h%type == GTH_TYPE_QS) then
                call RANDOM_NUMBER(tstRnd)
                if(tstRnd > 0.5) rC%D(2) = -rC%D(2)
             end if

             if(rC%D(3) < 1e-2_fd) then
                rC%D(3) = 1e-2
                call vec_normalize(rC%D)
             end if

             rC % status   = RAY_ACTIVE
             rC % I        = 1.0 / real(nSamplesPerOrderTable(1),fd)
             rC % P(1)     = pSurface(1) + dz * (rC%D(1) / rC%D(3))
             rC % P(2)     = pSurface(2) + dz * (rC%D(2) / rC%D(3))
             rC % P(3)     = M%grid%height - TRACE_EPS

             rC % P(1)     = modulo(rC%P(1)+M%hWidth, M%width) - M%hWidth
             rC % P(2)     = modulo(rC%P(2)+M%hWidth, M%width) - M%hWidth

             rC % rayID    = rC % rayID + RAY_ID_INCR

             !! Negate the ray.
             !!
             rC % D        = - rC % D

             rayStat       = RAY_OK
             order         = 1                
             do while(order <= nOrders .and. rayStat == RAY_OK)

                call rnd_generate_exponential(1.0_fd, l)
                !!call drandexponential(1, 1.0_fd, rState, l, rInfo)
                call trc_traceRayToOpticalDepth( M%grid, rC, l(1), matMu=mu )
                rC % rayID  = rC % rayID + RAY_ID_INCR

                if(rC%status == RAY_EXIT_U) then
                   rayStat = RAY_EXIT
                   exit
                else if(rC%status == RAY_EXIT_D) then
                   rayStat = RAY_LOW
                   exit
                end if

                if(rayStat == RAY_OK) then

                   do iTheta = 1, nThetaIn
                      rPeel        = rC
                      rPeel%status = RAY_ACTIVE
                      rPeel%D      = D(:,iTheta)
                      rC % rayID   = rPeel % rayID + RAY_ID_INCR

                      optDepth     = trc_tracePhysicalDepth(M%grid, rPeel) * mu

                      !$omp atomic
                      H%data(iTheta, H%cIdx(j)+k-1, order) = H%data(iTheta, H%cIdx(j)+k-1, order) &
                           & + INV_FOUR_PI * rC%I * w * exp(-optDepth)
                      
                   end do

                   rC%D  = vec_cart_random_spherical_uniform()
                   rC%I  = rC%I * w

                end if
                order = order + 1
             end do

          end do
          ! $omp end do
          call utl_timerIncrease(timer)  

       end do
    end do
    !$omp end do
    !$omp end parallel

  end subroutine sampleHemisphere_MC


  subroutine sampleDiscrete(nSurfaceSamples, iTheta, iPhi, s)
    integer                           :: nSurfaceSamples
    real(fd), dimension(:)            :: iTheta, iPhi, s
    real(fd), dimension(size(iTheta)) :: sWeight

    type(ray) :: rC
    type(ray) :: rS
    type(intersection_geometry) :: iSect

    real(fd), dimension(2,nSurfaceSamples) :: samples
    real(fd), dimension(3)                 :: pSurface
    real(fd)                               :: dz, thtIn

    integer  :: nAngles, i, k
    logical  :: pFound, pLit


    nAngles = size(iTheta)

    call smpl_griddedSamples2D(samples, nSurfaceSamples)

    samples = (samples * M%width - M%hWidth) * 0.5_fd
    thtIn   = thetaInTable(1)
    rS%D    = (/sin(thtIn), 0.0_fd, cos(thtIn)/)
    dz      = M%grid%height - M%hMean - TRACE_EPS

    pSurface(3) = M%hMean

    call ray_init(rC, RAY_TYPE_CAMERA)
    call ray_init(rS, RAY_TYPE_SHADOW)

    do k = 1, size(iTheta)
       do i= 1, nSurfaceSamples

          pSurface(1:2) = samples(:,i)

          rC % D        = (/sin(iTheta(k)) * cos(iPhi(k)), sin(iTheta(k)) * sin(iPhi(k)), cos(iTheta(k))/)

          if(rC%D(3) < 1e-2_fd) then
             rC%D(3) = 1e-2
             call vec_normalize(rC%D)
          end if

          call rayToBBTop(rC, M, pSurface, dz)

          pFound        = trc_traceNearest(M%grid, rC, iSect)

          if(pFound) then
             rS % P        = iSect % P1 + TRACE_EPS * iSect % N
             rS % rayID    = rS % rayID + RAY_ID_INCR

             pLit          = .not. trc_traceOcclusion(M%grid, rS)

             ! $omp atomic
             sWeight(k) = sWeight(k) + 1.0_fd

             if(pLit) then
                ! $omp atomic
                s(k)  = s(k)  + 1.0_fd
             end if
          end if
       end do
    end do

    where(sWeight > 0) 
       s = s / sWeight
    end where

  end subroutine sampleDiscrete

  subroutine simulation_method_montecarlo(med, theta_i, hs, mu, w, nSamples)
    type(med_medium)       :: med
    real(FD), dimension(:) :: theta_i
    type(gth_hemisphere)   :: hs
    real(fd)               :: mu, w
    integer(il)            :: nSamples

    integer, parameter :: RAY_OK   = 0
    integer, parameter :: RAY_EXIT = 1
    integer, parameter :: RAY_LOW  = 2

    ! Simulation (logical) variables
    !
    integer(il)                     :: iRay, iTht, iThtE, iPhiE
    integer                         :: rayStat, order

    real(fd), dimension(2,nSamples) :: samples
    real(fd), dimension(3)          :: pSurface
    real(fd)                        :: dz

    logical                         :: isInside

    type(intersection_geometry)     :: iSect
    type(utl_timer)                 :: timer

    ! Physical variables
    !
    real(fd)    :: l(1), optDepth
    type(ray)   :: r_s, r, r_peel
    real(fd)    :: pTemp(3)

    ! Random number generator variables
    !
    integer     :: rInfo, rSeed(1), rState(640)

    !$call omp_set_num_threads(nThreads)
    call smpl_griddedSamples2D(samples, nSamples)

    call utl_timer_init(timer, 1.0_fd, size(theta_i))
    do iTht = 1, size(theta_i)

       pSurface(3)  = med%hMean
       dz           = med%grid%height - med%hMean - TRACE_EPS
       
       r_s % status = RAY_ACTIVE
       r_s % P      = med%grid%height - TRACE_EPS
       r_s % D      = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]
       r_s % I      = 1.0_fd / real(nSamples, fd)
       r_s % energy = 1.0_fd
       
       !$omp parallel default(none) &
       !$omp shared(r_s, iTht, nSamples, med, hs, mu, w, theta_i, nOrders, pSurface, samples) &
       !$omp private(iRay, iThtE, iPhiE, rayStat, dz, r, r_peel, optDepth, l, rSeed, rState, rInfo, order)
       
       !$rSeed(1) = omp_get_thread_num()+1
       call drandinitialize(3, 0, rSeed, 1, rState, 640, rInfo)
       
       !$omp do schedule(dynamic) 
       do iRay = 1, nSamples
          
          pSurface(1:2) = samples(:,iRay)

          r          = r_s

          r % P(1)   = pSurface(1) + dz * (r%D(1) / r%D(3))
          r % P(2)   = pSurface(2) + dz * (r%D(2) / r%D(3))

          r % P(1)   = modulo(r%P(1)+med%hWidth, med%width) - med%hWidth
          r % P(2)   = modulo(r%P(2)+med%hWidth, med%width) - med%hWidth

          r % rayID  = r % rayID + RAY_ID_INCR

          rayStat    = RAY_OK
          order      = 1                

          do while(order <= nOrders .and. rayStat == RAY_OK)
             
             call drandexponential(1, 1.0_fd, rState, l, rInfo)
             call trc_traceRayToOpticalDepth( med%grid, r, l(1), matMu=mu )
             r % rayID  = r % rayID + RAY_ID_INCR

             if(r%status == RAY_EXIT_U) then
                rayStat = RAY_EXIT
                exit
             else if(r%status == RAY_EXIT_D) then
                rayStat = RAY_LOW
                exit
             end if

             if(rayStat == RAY_OK) then

                do iThtE = 1, hs%resTheta
                   do iPhiE = 1, hs%resPhi(iThtE)

                      r_peel        = r
                      r_peel%status = RAY_ACTIVE
                      r_peel%D      = gth_cellRandomSampleCar(hs, iThtE, iPhiE)
                      r % rayID     = r_peel % rayID + RAY_ID_INCR

                      optDepth = trc_tracePhysicalDepth(med%grid, r_peel) * mu

                      ! $omp atomic
                      !hs%data(iTht, hs%cIdx(iThtE)+iPhiE, order) = hs%data(iTht, h%cIdx(iThtE)+iPhiE, order) &
                      !     & + 0.5 * INV_FOUR_PI * r%I * w * exp(-optDepth)

                      !$omp critical 
                      call gth_addDataCar(hs, r_peel%D, 0.5 * INV_FOUR_PI * r%I * w * exp(-optDepth), iTht, order)
                      !call gth_addDataCar(hs, r_peel%D, r_peel%D(1), iTht, order)
                      !$omp end critical
                   end do
                end do

                r%D  = vec_cart_random_spherical_uniform()
                r%I  = r%I * w

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

  subroutine saveHemisphere(h, fName)
    type(gth_hemisphere), INTENT(IN)  :: h
    character (len = *)               :: fName

    type(gth_hsFile) :: hf
    integer          :: i

    integer            :: dimMedNum, dimMedMapRes, dimMedDenRes, idMedMap, idMedDen

    call gth_hsFileOpen(hf, fName, "w")
    call gth_hsWriteHeader(hf, h, "Hannu Parviainen", "vScatter", "0.95")
 
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Number_of_media",       MF%nSelectedMedia,  dimMedNum      ) )
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Medium_map_resolution", mediumMapRes,       dimMedMapRes   ) )
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Medium_densitymap_resolution", mediumDensitymapRes, dimMedDenRes   ) )

    Call gth_nfCheck( nf90_def_var(hf%fileID, "Medium_heightmap",  NF90_FLOAT, [dimMedNum,dimMedMapRes,dimMedMapRes], idMedMap) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Medium_densitymap", NF90_FLOAT, [dimMedNum,dimMedDenRes,2],            idMedDen) )

    call gth_hsAddAttS(hf, NF90_GLOBAL, "Brdf_type",                brdfType)
    call gth_hsAddAttF(hf, NF90_GLOBAL, "Theta_in",                 thetaInTable)
    call gth_hsAddAttS(hf, NF90_GLOBAL, "medium_filename",          mediumFileName)
    call gth_hsAddAttS(hf, NF90_GLOBAL, "medium_size_distribution", sDistribution)

    call gth_hsAddAttF(hf, NF90_GLOBAL, "Single_scattering_albedo", [w])
    call gth_hsAddAttF(hf, NF90_GLOBAL, "extinction_coefficient", [mu])


    call gth_hsAddAttF(hf, NF90_GLOBAL, "nMedia", [real(MF%nSelectedMedia,fd)] )

    if(rf_applyfields) then
       call gth_hsAddAttS(hf, NF90_GLOBAL, "srfType",    rf_spectrumType  )
       call gth_hsAddAttF(hf, NF90_GLOBAL, "nFields",    [real(rf_nFieldsPerMed,fd)] )
       call gth_hsAddAttF(hf, NF90_GLOBAL, "srfStd",     [rf_P]            )
       call gth_hsAddAttF(hf, NF90_GLOBAL, "srfRoughP",  [rf_std]          )
    end if

    call gth_hsSaveData(h, hf)

    Call gth_nfCheck( nf90_put_var(hf%fileID, idMedMap,         mediumHeightMap) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMedDen,         mediumDensityStructure) )

    call gth_hsFileClose(hf)

  end subroutine saveHemisphere

end program vScatter

