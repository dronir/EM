  subroutine simulation_method_montecarlo(mat, med, spc, theta_i, hs, nSamples)
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

    type(mat_material)     :: mat
    type(med_medium)       :: med
    type(spc_spectrum)     :: spc
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
    real(fd)    :: muFluorLine, muExt, fYield, l, tau
    type(ray)   :: r_s, r, rI

    ! Random number generator variables
    !
    real(fd)    :: rTable(SRTABLE)
    integer     :: rInfo, rSeed(1), rState(640)
    integer     :: sampleSeed

    sampleSeed = 0

    !$ call omp_set_num_threads(nThreads)

    call smpl_griddedSamples2D(samples, nSamples)

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

       !$omp do schedule(dynamic) 
       do iRay = 1, nSamples
          
          if(iRand == 1) call rnd_generate_uniform_n(rTable)
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

                l = rnd_exponential(1.0_fd)
                call trc_traceRayToOpticalDepth(med%grid, r, l, mat)
                r % rayID  = r % rayID + RAY_ID_INCR
          

                if(r%status == RAY_EXIT_D) then
                   rayStat = RAY_EXIT
                   exit
                end if

                !! Attenuate intinsity due to the extinction methods not used in the simulation
                !!
                r%I  = r%I * exp(- l * muExt)


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
