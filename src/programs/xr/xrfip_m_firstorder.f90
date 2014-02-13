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
!! \param[in]  spc          Spectrum
!! \param[in]  nPoSamples   Number of samples
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
subroutine infslab_method_firstorder(mat, spc, nPosSamples, nSrcSamples, inTheta, emTheta, emPhi, outType, verbose, hs, tb, nThreads)
  use base
  use material
  use geometry
  use hemisphere
  use random
  !$ use omp_lib

  implicit none

  type(mat_material)       :: mat
  type(spc_spectrum)       :: spc
  integer(il)              :: nPosSamples
  integer(il)              :: nSrcSamples
  real(fd), dimension(:)   :: inTheta
  real(fd), dimension(:)   :: emTheta
  real(fd), dimension(:)   :: emPhi
  integer                  :: outType
  logical                  :: verbose

  type(gth_hemisphere),       optional :: hs
  real(fd), dimension(:,:,:), optional :: tb

  integer, optional        :: nThreads

  !!- Simulation
  !!
  integer(il)            :: iSpc, iRay, iRand, iLine, iThtI, iThtE, iEmAngle
  real(fd)               :: fluorescence_ratio, x

  real(fd), allocatable  :: xrf_partial(:,:)
  real(fd)               :: partial_limit
  integer(il)            :: partial_idx
  integer                :: i, j, k

  integer                :: nEmAngles
  real, allocatable      :: em_angles(:,:)
  real, allocatable      :: em_vectors(:,:)
  integer, allocatable   :: em_vector_idx(:,:)

  type(utl_timer)        :: timer

  !!- Physics
  !!
  real(fd), dimension(mat%nLines) :: muLineExt, fluorescence_intensity, fYield, line_ratio
  real(fd)                        :: muFluor, muTotal,  l(1)
  type(ray)                       :: r_ref, r

  !!- Random number generator
  !!
  integer     :: rInfo, rSeed(1), rState(640)

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
  !!- Precompute the per-line coefficients (total extinction coefficient and fluorescence yield)
  !!
  muLineExt = 0.0_fd
  do iLine = 1, mat%nLines
     if(mat%lEnergy(iLine) > 0.0_fd) then
        muLineExt(iLine) = mat_evalMuExt(mat, mat%lEnergy(iLine))
        fYield(iLine)    = mat%fYield((iLine-1)/2 + 1)
     end if
  end do
  !!
  !!- Force number of threads if using OpenMP
  !!
  !$ if(present(nThreads) .and. nThreads > 0) then
  !$   call omp_set_num_threads(nThreads)
  !$ end if
  !!
  !!- Initialize timer
  !!
  if(verbose) then
     call utl_timer_init(timer, 5.0_fd, size(inTheta) * nPosSamples * nSrcSamples)
  end if

  partial_limit = 1.0_fd/real(nPosSamples * nSrcSamples,fd) * 1.0e0_fd

  !!--- SIMULATION PHASE ---
  !!
  do iThtI = 1, size(inTheta)
     !$omp parallel default(none) &
     !$omp shared       ( hs, tb, mat, spc, timer, verbose, nPosSamples, nSrcSamples ) &
     !$omp shared       ( r_ref, inTheta, nEmAngles, em_vectors, em_vector_idx, iThtI, partial_limit  ) &
     !$omp private      ( iSpc, iThtE, iRay, iLine, iEmAngle, r, fluorescence_ratio, x, xrf_partial, partial_idx, i, j, k) &
     !$omp private      ( fluorescence_intensity, line_ratio, muFluor, muTotal, l, rSeed, rState, rInfo ) &
     !$omp firstprivate ( iRand, muLineExt, fYield, outType)

     
     !!---   ---
     !!
     !! 
     !!
     !$omp single
     r_ref % P = 0.0_fd
     r_ref % D = [ sin(inTheta(iThtI)), 0.0_fd, -cos(inTheta(iThtI)) ]
     r_ref % I = 1.0_fd / real(nPosSamples, fd)
     !$omp end single

     !!- Initialize the random number generator.
     !!
     iRand   = 1
     
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

     !!--- LOOP OVER SPECTRUM SAMPLES  ---
     !!
     !! 
     !$omp do schedule(dynamic) 
     do iSpc = 1, nSrcSamples

        !!--- LOOP OVER SURFACE SAMPLES  ---
        !!
        !! 
        ! $omp do schedule(dynamic) 
        do iRay = 1, nPosSamples

           !!--- INCIDENT RAY GENERATION ---
           !!
           !!--- INITIALIZE RAY ---
           !!
           r          = r_ref
           r % I      = 1.0_fd/real(nPosSamples * nSrcSamples,fd) 
           r % energy = spc_getSample(spc, rnd_uniform_n())

           !!--- TRACE RAY TO MEDIA ---
           !!
           !! Evaluate the fluorescence and total extinction coefficients, generate a random optical depth, 
           !! and trace the ray to the optical depth generated.
           !!
           call mat_evalMu(mat, r%energy, muFluorLineTotal = muFluor, muTotal = muTotal)

           l    = rnd_exponential(1.0_fd)
           l    = l / muTotal
           r%P  = r%P + r%D * l(1)

           !!--- FLUORESCENCE EVALUATION ---
           !!
           !! Compute the ratio of fluorescence over total extinction, and the fluorescence for each
           !! fluorescent line.
           !!
           fluorescence_ratio     = muFluor / muTotal
           line_ratio             = mat_evalFluorLineFractions(mat, r%energy)


           !!--- LOOP OVER EMISSION ANGLES ---
           !!
           !!
           do iEmAngle = 1, nEmAngles
              fluorescence_intensity = &
                   & INV_FOUR_PI * r%I * fYield * fluorescence_ratio * line_ratio &
                   & * exp(-muLineExt * abs(r%P(3)) / em_vectors(iEmAngle, 3))

              xrf_partial(iEmAngle, :)  = xrf_partial(iEmAngle, :)  + fluorescence_intensity(:)

           end do

           if (verbose) then
              call utl_timerIncrease(timer)
           end if

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
        end do
        ! $omp end do
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

end subroutine infslab_method_firstorder
!******
