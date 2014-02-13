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
subroutine infslab_method_montecarlo(mat, spc, nPosSamples, nSrcSamples, nOrders, obsHalfAngle, theta_i, theta_e, outType, hs, tb, nThreads)
  use base
  use material
  use hemisphere
  use geometry
!$ use omp_lib

  implicit none

  type(mat_material)         :: mat
  type(spc_spectrum)         :: spc
  integer(il)                :: nPosSamples
  integer(il)                :: nSrcSamples
  integer(il)                :: nOrders
  real(fd)                   :: obsHalfAngle
  real(fd), dimension(:)     :: theta_i
  real(fd), dimension(:)     :: theta_e
  integer                    :: outType
  type(gth_hemisphere)       :: hs
  real(fd), dimension(:,:,:) :: tb
  integer, optional          :: nThreads

  !!- Constants
  !!
  integer,  parameter    :: SRTABLE       = 20000
  real(fd), parameter    :: MIN_MU        = 1e-20_fd
  real(fd), parameter    :: MIN_INTENSITY = 1e-23_fd

  !!- Simulation 
  !!
  integer(il) :: iSampleset, iSpec, iRay, iRand, iLine, iThtI, iThtE
  integer     :: rayStat, order
  real(fd)    :: dObs(size(theta_e),3), cos_half_angle
  real(fd)    :: zObs(size(theta_e))
  real(fd)    :: obsSolidAngle !!FIXME: Compute me!!

  type(utl_timer) :: timer

  !!- Physics
  !!
  real(fd)    :: muFluorLine, muTotal, fYield, l(1), fluorTemp(size(theta_e))
  type(ray)   :: r_ref, r

  !!- Random number generator 
  !!
  real(fd)    :: rTable(SRTABLE)
!  integer     :: rInfo, rSeed(1), rState(640)

  !!--- INITIALIZATION PHASE ---
  !!
  !!- Precompute emergence angle cosines
  !!
  if(outType == OUTPUT_TYPE_HEMISPHERE) then
     zObs = cos(hs%mTheta)
  else
     zObs = cos(theta_e)
  end if
  !!
  !!- Initialize incidence angle vectors
  !!
  do iThtE = 1, size(theta_e)
     dObs(iThtE, :) = [sin(theta_e(iThtE)), 0.0_fd, cos(theta_e(iThtE))]
  end do
  !!
  cos_half_angle = cos(obsHalfAngle)
  !!
  !!- Force number of threads if using OpenMP
  !!
  !$ if(present(nThreads)) then
  !$   call omp_set_num_threads(nThreads)
  !$ end if
  !!
  !!- Initialize timer
  !!
  call utl_timer_init(timer, 5.0_fd, size(theta_i) * nPosSamples * nSrcSamples)


  !!--- SIMULATION PHASE ---
  !!
  do iThtI = 1, size(theta_i)

     !$omp parallel default(none) &
     !$omp shared(mat, spc, hs, tb, outType, nPosSamples, nOrders, nSrcSamples, timer,  zObs, dObs) &
     !$omp shared(iThtI, theta_i, obsSolidAngle, cos_half_angle) &
     !$omp private(iRay, iSpec, iSampleSet, iLine, iThtE, r_ref, fluorTemp) &
     !$omp private(rayStat, r, rTable, muFluorLine, muTotal, l, fYield, iRand)

     r_ref % P      = 0.0_fd
     r_ref % D      = [ sin(theta_i(iThtI)), 0.0_fd, -cos(theta_i(iThtI)) ]
     r_ref % order  = 0_il
     r_ref % status = RAY_ACTIVE

     iRand = 1

     !$omp do schedule(guided, 100)
     do iSpec = 1, nSrcSamples
        do iRay = 1, nPosSamples

           if(iRand == 1) call rnd_generate_uniform_n(rTable)
           iRand      = mod(iRand, SRTABLE) + 1

           !!--- INITIALIZE RAY ---
           !!
           r          = r_ref
           r % I      = 1.0_fd / real(nPosSamples*nSrcSamples, fd)
           r % energy = spc_getSample(spc, rTable(iRand))

           !!--- RAY TRAVEL PHASE ---
           !!
           !! First, evaluate the total extinction coefficient and the fluorescence line coefficient for the
           !! current ray energy. If the ray energy is above the minimum edge energy, and if the fluorescence
           !! line coefficient is above the minimum allowed value, compute random length for the ray to travel
           !! from exponenetial distribution scaled by the total extinction coefficient.
           !!
           do while(r%status == RAY_ACTIVE .AND. r%order <= nOrders)

              call mat_evalMu(mat, r%energy, muFluorLineTotal = muFluorLine, muTotal = muTotal)

              if(r%energy > mat%edgeEnergyMin .AND. muFluorLine > MIN_MU) then

                 !call drandexponential(1, 1.0_fd, rState, l, rInfo)
                 l = rnd_exponential(1.0_fd)
                 l = l / muTotal

                 !!--- FLUORESCENCE PHASE ---
                 !!
                 !! If the ray will be inside the medium after the travel, and we still want to compute more
                 !! fluorescence orders, evaluate fluorescence for the current ray energy. This will change
                 !! the ray energy to correspond the chosen fluorecence line, and return the line index and
                 !! and fluorescence yield.
                 !!
                 !! Next, we update the position of the ray and the intensity of the ray, and compute a random
                 !! fluorscence direction.
                 !!
                 if(r%order < nOrders .AND. r%P(3) < -r%D(3)*l(1)) then  
                    call mat_evalFluorescence(mat, r%energy, r%energy, fYield, iLine)

                    r%P  = r%P + r%D * l(1)  
                    r%D  = vec_cart_random_spherical_uniform()
                    r%I  = r%I * fYield * (muFluorLine/muTotal)

                    call mat_evalMu(mat, r%energy, muTotal = muTotal)

                    !!--- PEELING PHASE ---
                    !!
                    !! We compute the fluorescence signal from the position of the ray to the gathering elements.
                    !!
                    !!- Hemisphere branch
                    !!
                    if(outType == OUTPUT_TYPE_HEMISPHERE) then
                       do iThtE = 1, hs%resTheta
                          !$omp critical
                          call gth_hs_addToRow(hs, iThtE, &
                               & INV_FOUR_PI * r%I * exp(-muTotal * abs(r%P(3)) / zObs(iThtE)), iThtI, iLine)
                          !$omp end critical

                          !!--- OPTIONAL FORM ---
                          !!
                          !! For some strange reason, this seems to be faster for multithreaded simulation...
                          !!
                          !! !$omp critical
                          !! call gth_hs_addToRow(hs, iThtE, &
                          !!     & INV_FOUR_PI * r%I * exp(-muTotal * abs(r%P(3)) / cos(hs%mTheta(iThtE))),&
                          !!     & iThtI, iLine)
                          !! !$omp end critical

                       end do
                       !!
                       !!- Tabulated branch
                       !!
                    else
                       !fluorTemp = INV_FOUR_PI * r%I * exp(-muTotal * abs(r%P(3)/zObs))
                       !$omp critical
                       tb(iThtI, :, iLine) = tb(iThtI, :, iLine) + INV_FOUR_PI * r%I * exp(-muTotal * abs(r%P(3)/zObs))
                       !$omp end critical
                    end if
                    !!
                    !!--- END PEELING PHASE ---
                    !!
                    !! If order > nOrders or the ray escapes from the medium, kill the ray.
                    !!
                 else
                    r%status = RAY_DEAD
                 end if
                 !!
                 !! If ray energy < minimum edge energy or fluorescence line coefficient < minimum allowed, kill the ray.
                 !!
              else
                 r%status = RAY_DEAD
              end if

              r%order = r%order + 1
           end do
           !!
           !!--- END RAY TRAVEL PHASE ---
           !!
        end do
        call utl_timerIncrease(timer)
     end do
     !$omp end do
     !$omp end parallel
  end do

end subroutine infslab_method_montecarlo
!******
