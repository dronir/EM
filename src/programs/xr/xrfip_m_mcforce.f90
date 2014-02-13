!!--- BEGIN GPL --- 
!!
!! This file is part of UHOEM.
!!
!! UHOEM is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! UHOEM is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with UHOEM.  If not, see <http://www.gnu.org/licenses/>.
!!
!! Copyright (C) 2009 Hannu Parviainen
!!
!! Contributor(s): Hannu Parviainen
!!
!!--- END GPL ---

!****f* xrfipsd/infslab_method_montecarlo_force
! NAME
!   infslab_method_montecarlo_force
!
! DESCRIPTION
!   Computes the fluorescent signal from an infinite plane using Monte Carlo ray tracing without peeling.
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
subroutine infslab_method_montecarlo_force(mat, spc, nSamples, theta_i, theta_e, outType, hs, tb)
  use base
  use material
  use hemisphere
  use geometry

  type(mat_material)         :: mat
  type(spc_spectrum)         :: spc
  integer(il)                :: nSamples
  real(fd), dimension(:)     :: theta_i
  real(fd), dimension(:)     :: theta_e
  integer                    :: outType
  type(gth_hemisphere)       :: hs
  real(fd), dimension(:,:,:) :: tb


  integer,  parameter    :: SRTABLE       = 20000
  real(fd), parameter    :: MIN_MU        = 1e-20_fd
  real(fd), parameter    :: MIN_INTENSITY = 1e-23_fd

  ! Simulation (logical) variables
  !
  integer(il) :: iSampleset, iSpec, iRay, iRand, iLine, iTht, iThtE
  integer     :: rayStat, order
  real(fd)    :: dObs(size(theta_e),3), cos_half_angle
  real(fd)    :: zObs(size(theta_e))

  type(utl_timer) :: timer

  ! Physical variables
  !
  real(fd)    :: muFluorLine, muTotal, fYield, l(1)
  type(ray)   :: r_ref, r

  ! Random number generator variables
  !
  real(fd)    :: rTable(SRTABLE)

  zObs = cos(theta_e)

  do iThtE = 1, size(theta_e)
     dObs(iThtE, :) = [sin(theta_e(iThtE)), 0.0_fd, cos(theta_e(iThtE))]
  end do

  cos_half_angle = cos(obsHalfAngle)

  call utl_timer_init(timer, 5.0_fd, size(theta_i) * nSampleMultiplier * spc%nPoints) !nSampleMultiplier)
  do iTht = 1, size(theta_i)

     r_ref % P      = 0.0_fd
     r_ref % D      = [ sin(theta_i(iTht)), 0.0_fd, -cos(theta_i(iTht)) ]

     iRand = 1
     rSeed = 1

     !$omp parallel default(none) &
     !$omp shared(tb, iTht, r_ref, nSamples, mat, spc, theta_i, nOrders, nSampleMultiplier, timer, obsSolidAngle) &
     !$omp private(iRay,iSpec,iSampleSet,iLine,iThtE, rayStat, r, rTable, muFluorLine, muTotal, l, fYield, order) &
     !$omp firstprivate(iRand, dObs, cos_half_angle)


     do iSampleSet = 1, nSampleMultiplier
        !$omp do schedule(guided, 100)
        do iSpec = 1, spc%nPoints
           do iRay = 1, nSamples

              if(iRand == 1) call rnd_generate_uniform_n(rTable) 
              iRand      = mod(iRand, SRTABLE) + 1

              order      = 0
              r          = r_ref
              r % I      = spc%I(iSpec) / real(nSamples*nSampleMultiplier, fd) ! 1.0_fd / real(nSamples*nSampleMultiplier, fd)
              r % energy = spc%E(iSpec) !spc_getSample(spc, rTable(iRand))
              r % status = RAY_ACTIVE

              do while(r%status == RAY_ACTIVE .AND. order <= nOrders)
                 call mat_evalMu(mat, r%energy, muFluorLineTotal = muFluorLine, muTotal = muTotal)

                 if(r%energy > mat%edgeEnergyMin .AND. muFluorLine > MIN_MU) then

                    l = rnd_exponential(1.0_fd)
                    l = l / muTotal

                    r%P  = r%P + r%D * l(1)  

                    if(r%P(3) < 0.0) then    
                       if(order < nOrders) then
                          call mat_evalFluorescence(mat, r%energy, r%energy, fYield, iLine)
                          r%D  = vec_cart_random_spherical_uniform()
                          r%I  = r%I * fYield * (muFluorLine/muTotal)
                       else
                          r%status = RAY_EXIT_D
                       end if
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

!!$                !$omp critical
!!$                call gth_hs_addToRowCar(hs, r%D, r%I, iTht, iLine)
!!$                !$omp end critical


                 do iThtE = 1, size(theta_e)
                    if(dot_product(r%D, dObs(iThtE,:)) > cos_half_angle) then
                       !$omp atomic
                       tb(iTht, iThtE, iLine) = tb(iTht, iThtE, iLine) + r%I / obsSolidAngle
                       exit
                    end if
                 end do
              end if

           end do
           call utl_timerIncrease(timer)
        end do
        !$omp end do
     end do
     !$omp end parallel
  end do
end subroutine infslab_method_montecarlo_force
