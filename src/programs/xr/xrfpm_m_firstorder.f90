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

!> Computes the fluorescent signal using a first-order Monte-Carlo approximation.
!!
!! \param[in]  mat          Material
!! \param[in]  med          Medium
!! \param[in]  spc          Spectrum
!! \param[in]  nPosSamples  Number of position samples
!! \param[in]  nSrcSamples  Number of spectrum samples
!! \param[in]  inTheta      Angles of incidence
!! \param[in]  emTheta      Emergent theta angles
!! \param[in]  emPhi        Emergent phi angles
!! \param[in]  outType      Type of output, can be hemisphere or tabulated
!! \param[in]  verbose      Should we be verbose?
!!
!! \param[out] hs           Output hemisphere
!! \param[out] tb           Output table
!!
!! \param[in]  nThreads     Number of threads, automatic if < 0
!!
subroutine simulation_method_firstorder(mat, med, spc, nPosSamples, nSrcSamples, inTheta, emTheta, emPhi, outType, verbose, timer, hs, tb, nThreads)
  use base
  use material
  use geometry
  use particle
  use hemisphere
  use random
  use medium
  use trace
  use sampler
  !$ use omp_lib

  implicit none

  type(mat_material)       :: mat
  type(med_medium)         :: med
  type(spc_spectrum)       :: spc
  integer(il)              :: nPosSamples
  integer(il)              :: nSrcSamples
  real(fd), dimension(:)   :: inTheta
  real(fd), dimension(:)   :: emTheta, emPhi
  integer                  :: outType
  logical                  :: verbose
  type(utl_timer)          :: timer

  type(gth_hemisphere),       optional :: hs
  real(fd), dimension(:,:,:), optional :: tb

  integer, optional        :: nThreads

  !!- Constants
  !!
  integer,  parameter    :: RAY_OK   = 0
  integer,  parameter    :: RAY_EXIT = 1
  integer,  parameter    :: RAY_LOW  = 2
  integer,  parameter    :: SRTABLE  = 20000
  real(fd), parameter    :: MIN_MU        = 1e-20_fd
  real(fd), parameter    :: MIN_INTENSITY = 1e-23_fd


  !!- Simulation
  !!
  integer(il)                     :: iSpc, iRay, iRand, iLine, iThtI
  integer                         :: order
  real(fd)                        :: fluorescence_ratio, spcNorm, x

  real(fd), allocatable           :: xrf_partial(:,:)
  real(fd)                        :: partial_limit
  integer(il)                     :: partial_idx
  integer                         :: i, j, k

  real(fd), dimension(2,nPosSamples) :: samples
  real(fd)                        :: dz, xEn, pTemp(3)
  integer                         :: iEn, iThtE, iPhiE, iEmAngle

  logical                         :: isInside !< Is the ray inside a particle (medium?)
  type(intersection_geometry)     :: iSect

  integer                         :: nEmAngles
  real(fd), allocatable           :: em_vectors(:,:)
  integer, allocatable            :: em_vector_idx(:,:)

  real(fd)                        :: dObs(size(emTheta),3)
  real(fd)                        :: zObs(size(emTheta))

  !!- Physics
  !!
  real(fd), dimension(mat%nLines) :: muLineExt, line_ratio, lFluorescenceYield
  real(fd), dimension(mat%nLines) :: fluorescence_intensity
  real(fd)                        :: l, tau, muFluor, muTotal

  type(ray)   :: r, rI, r_ref

  !!- Random number generator
  !!
  real(fd)    :: rTable(SRTABLE)
  integer     :: rInfo, rSeed(1), rState(640)

  integer     :: sampleSeed

  sampleSeed = 0

  !!--- INITIALIZATION PHASE ---
  !!
  !!- Precompute emergence vectors
  !!
  if(outType == OUTPUT_TYPE_TABULATED) then
     nEmAngles = size(emTheta)
     allocate(em_vectors(nEmAngles, 3))
     do i = 1, nEmAngles
        em_vectors(i, :) = [sin(emTheta(i)) * cos(emPhi(i)), &
             &              sin(emTheta(i)) * sin(emPhi(i)), &
             &              cos(emTheta(i))]
     end do
  else
     nEmAngles = hs%nCells
     allocate(em_vectors    (nEmAngles, 3) )
     allocate(em_vector_idx (nEmAngles, 2) )

     k = 1
     do i = 1, hs%resTheta
        do j = 1, hs%resPhi(i)             
           em_vectors    (k, :) = gth_cellCenterCar(hs, i, j)
           em_vector_idx (k, :) = [i,j]
           k = k + 1
        end do
     end do
  end if
  !!
  !!- Precompute incidence vectors
  !!
  do iThtE = 1, size(emTheta)
     dObs(iThtE, :) = [sin(emTheta(iThtE)), 0.0_fd, cos(emTheta(iThtE))]
  end do

  !!
  !!- Precompute the per-line coefficients (total extinction coefficient and fluorescence yield)
  !!
  do iLine = 1, mat%nLines
     if(mat%lEnergy(iLine) > 0.0_fd) then
        muLineExt(iLine)          = mat_evalMuExt(mat, mat%lEnergy(iLine))
        lFluorescenceYield(iLine) = mat%fYield((iLine-1)/2 + 1)
     end if
  end do
  !!
  !!- Force number of threads if using OpenMP
  !! 
  !$ if(present(nThreads) .and. nThreads > 0) then
  !$   call omp_set_num_threads(nThreads)
  !$ end if
  !!
  !!- Generate mean-surface sample point locations
  !!
  call smpl_griddedSamples2D(samples, nPosSamples)

  partial_limit = 1.0_fd/real(nPosSamples * nSrcSamples,fd) * 1.0e0_fd

  !!--- SIMULATION PHASE ---
  !!
  do iThtI = 1, size(inTheta)
     !$omp parallel default(none) &
     !$omp shared (med, mat, spc, hs, timer, outType) &
     !$omp shared (nPosSamples, nSrcSamples, nEmAngles) &
     !$omp shared (iThtI, inTheta, emTheta, emPhi, em_vectors) &
     !$omp shared (samples, muLineExt, lFluorescenceYield, r_ref, partial_limit) &
     !$omp private(iRay, iLine, iRand, iSpc, iEmAngle, dz, r, tau, iEn, xEn, l, fluorescence_ratio, fluorescence_intensity) &
     !$omp private(rTable, rSeed, rState, rInfo, pTemp, iThtE, iPhiE, muFluor, muTotal, line_ratio, xrf_partial, partial_idx, i, j, k) 

     !!---   ---
     !!
     !! 
     !!
     !$omp single
     r_ref % P(3)   = med%grid%height - TRACE_EPS
     r_ref % D      = [ sin(inTheta(iThtI)), 0.0_fd, -cos(inTheta(iThtI)) ]
     r_ref % I      = 1.0_fd/real(nPosSamples * nSrcSamples,fd)
     r_ref % order  = 0_il
     r_ref % status = RAY_ACTIVE
     r_ref % rayID  = 0
     !$omp end single

     dz             = med%grid%height - med%hMean - TRACE_EPS

     !!- Initialize the random number generator.
     !!
     iRand = 1

     !!- Initialize the partial sum table if not yet initialized
     !!
     if(.not. allocated(xrf_partial)) then
        if(outType == OUTPUT_TYPE_TABULATED) then
           allocate(xrf_partial(nEmAngles, mat%nLines))
        else
           allocate(xrf_partial(hs%nCells, mat%nLines))
        end if
     end if
     !!
     xrf_partial = 0.0_fd
     partial_idx = 1
     !!
     !!

     !!--- LOOP OVER SPECTRUM SAMPLES  ---
     !!
     !! 
     !$omp do schedule(dynamic) 
     do iSpc = 1, nSrcSamples

        !!---  LOOP OVER MEDIUM SAMPLES ---
        !!
        do iRay = 1, nPosSamples

           !!--- INCIDENT RAY GENERATION ---
           !!
           !! Generate random numbers if random number table is empty
           !!
           if(iRand == 1) call rnd_generate_uniform_n(rTable)
           !!
           !!--- INITIALIZE RAY ---
           !!
           call ray_init(r, RAY_TYPE_ENERGY)
           r     % D      = r_ref % D
           r     % I      = r_ref % I
           r     % energy = spc_getSample(spc, rTable(iRand))
           iRand          = mod(iRand, SRTABLE) + 1
           !!
           !! Compute the incident ray starting (x,y) position at the top of the bounding box of the medium.
           !!
           r     % P(1:2) = samples(:, iRay) + dz * (r%D(1:2) / r%D(3))
           r     % P(1:2) = modulo(r%P(1:2)+med%hWidth, med%width) - med%hWidth
           r     % P(3)   = r_ref % P(3)
           !!
           !!--- TRACE RAY TO MEDIA ---
           !!
           !! Evaluate the fluorescence and total extinction coefficients, generate a random optical depth, 
           !! and trace the ray to the optical depth generated.
           !!
           call mat_evalMu(mat, r%energy, muFluorLineTotal = muFluor, muTotal = muTotal)

           l = rnd_exponential(1.0_fd) 
           call trc_traceRayToOpticalDepth(med%grid, r, l, mat)

           if(r%status == RAY_ACTIVE) then

              !!--- fLUORESCENCE EVALUATION ---
              !!
              !! Compute the ratio of fluorescence over total extinction, and the fluorescence for each
              !! fluorescent line.
              !!
              fluorescence_ratio = muFluor / muTotal
              line_ratio         = mat_evalFluorLineFractions(mat, r%energy)
              pTemp              = r%P

              !!--- LOOP OVER EMISSION ANGLES ---
              !!
              !!
              do iEmAngle = 1, nEmAngles
                 
                 r     % rayid = ray_get_newid(r%raytype)
                 r     % status = RAY_ACTIVE
                 r     % P      = pTemp
                 r     % D      = em_vectors(iEmAngle, :)

                 !! We need to enable (and fix) this if we want to use multiple materials in the medium.
                 !!
                 !! tau    = trc_traceOpticalDepth(med%grid, mat, r)
                 
                 tau       = trc_tracePhysicalDepth(med%grid, r)

                 fluorescence_intensity = &
                      & INV_FOUR_PI * r%I * lFluorescenceYield * fluorescence_ratio * line_ratio &
                      & * exp(-muLineExt * tau)

                 xrf_partial(iEmAngle, :) = xrf_partial(iEmAngle, :) + fluorescence_intensity(:)

              end do

              !!--- UPDATE PARTIAL SUMS ---
              !!
              if(mod(partial_idx, 1000) == 0) then
                 if(count(xrf_partial > partial_limit) > 5) then

                    !!--- TABULATED BRANCH ---
                    !!
                    if(outType == OUTPUT_TYPE_TABULATED) then
                       !$omp critical
                       forall(j = 1 : nEmAngles, k = 1 : mat%nLines, xrf_partial(j,k) > partial_limit )
                          tb(iThtI, j, k)  = tb(iThtI,j,k) + xrf_partial(j,k)
                          xrf_partial(j,k) = 0.0_fd
                       end forall
                       !$omp end critical

                       !!--- HEMISPHERE BRANCH ---
                       !!
                    else
                       !$omp critical
                       forall(j = 1 : nEmAngles, k = 1 : mat%nLines, xrf_partial(j,k) > partial_limit )
                          hs%data(iThtI, j, k) = hs%data(iThtI, j, k) + xrf_partial(j,k)
                          xrf_partial(j,k)      = 0.0_fd
                       end forall
                       !$omp end critical
                    end if
                 end if
              end if
              partial_idx = partial_idx + 1
           end if

           ! $omp critical
           call utl_timerIncrease(timer)
           ! $omp end critical

        end do
     end do
     !$omp end do

     !!--- FINAL PARTIAL SUM UPDATE ---
     !!
     if(outType == OUTPUT_TYPE_TABULATED) then
        !$omp critical
        tb(iThtI, :, :) = tb(iThtI, :, :) + xrf_partial(:, :)
        !$omp end critical
     else
        !$omp critical
        hs%data(iThtI, :, :) = hs%data(iThtI, :, :) + xrf_partial(:, :)
        !$omp end critical
     end if
     !$omp end parallel

  end do
end subroutine simulation_method_firstorder
