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


!> Computes the fluorescent signal from a particulate medium plane using first-order analytic approximation.
!!
!! \param[in]  mat          Material
!! \param[in]  med          Medium
!! \param[in]  spc          Spectrum
!! \param[in]  nPosSamples  Number of position samples
!! \param[in]  inTheta      Angle(s) of incidence
!! \param[in]  emTheta      Angle(s) of emergence
!! \param[in]  outType      Output type
!!
!! \param[out] hs           Fluorescent signal for each line and incidence angle
!! \param[out] tb           Fluorescent signal for each line and incidence angle
!!
subroutine simulation_method_analytic(mat, med, spc, nPosSamples, inTheta, emTheta, emPhi, outType, timer, hs, tb)
  use base
  use material
  use geometry
  use hemisphere
  use particle
  use bfdf
  use medium
  use trace_xr
  use sampler

  implicit none

  type(mat_material)       :: mat
  type(med_medium)         :: med
  type(spc_spectrum)       :: spc
  integer(il)              :: nPosSamples
  real(fd), dimension(:)   :: inTheta
  real(fd), dimension(:)   :: emTheta
  real(fd), dimension(:)   :: emPhi
  integer                  :: outType
  type(utl_timer)          :: timer

  type(gth_hemisphere)               :: hs
  real(fd), dimension(:,:,:), target :: tb

  real(fd), dimension(mat%nLines)      :: fluorescence, muLineTotal, fYield

  real(fd)                             :: posSamples(2, nPosSamples)

  type(ray)                            :: rShd
  type(intersection_geometry)          :: iSect
  logical                              :: isFound, pShadowed

  real(fd)                             :: in_vectors(size(inTheta),3)

  integer                              :: nEmAngles
  real(fd), allocatable                :: em_vectors(:,:)
  integer, allocatable                 :: em_vector_idx(:,:)

  real(fd)                             :: inAngle, emAngle

  integer     :: iInAngle, iEmAngle, iLine, iRand, iPosSample, i, j, k
  integer     :: sampleSeed

  real(fd), dimension(:,:,:), pointer  :: xrfr !< Pointer to the result array
 
  !!--- INITIALIZATION ---
  !!
  sampleSeed = 0
  !!
  !!- Precompute incidence vectors
  !!
  do iInAngle = 1, size(inTheta)
     in_vectors(iInAngle, :) = [sin(inTheta(iInAngle)), 0.0_fd, cos(inTheta(iInAngle))]
  end do
  !!
  !!- Precompute emergence vectors
  !!
  if(outType == OUTPUT_TYPE_TABULATED) then
     xrfr      => tb
     nEmAngles = size(emTheta)
     allocate(em_vectors(nEmAngles, 3))
     do i = 1, nEmAngles
        em_vectors(i, :) = [sin(emTheta(i)) * cos(emPhi(i)), &
             &              sin(emTheta(i)) * sin(emPhi(i)), &
             &              cos(emTheta(i))]
     end do
  else
     xrfr      => hs%data
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
  !!- Precompute fluorescence yield and total mu
  !!
  do iLine = 1, mat%nLines
     if(mat%lEnergy(iLine) > 0.0_fd) then
        call mat_evalMu(mat, mat%lEnergy(iLine), muTotal = muLineTotal(iLine))
        fYield(iLine) = mat%fYield((iLine-1)/2 + 1)
     end if
  end do
  !!
  !!- Generate mean-surface sample point locations
  !!
  call smpl_griddedSamples2D(posSamples, nPosSamples)
  posSamples = (posSamples * med%width - med%hWidth) * 0.5_fd
  !!
  !!--- END INITIALIZATION ---

  !!--- SIMULATION ---
  !!
  do iEmAngle = 1, nEmAngles
     !$omp parallel do default(none) &
     !$omp shared(iEmAngle, timer, med, em_vectors, in_vectors, nPosSamples, posSamples, inTheta, xrfr, mat, spc) &
     !$omp private(iPosSample, isFound, iInAngle, isect, emAngle, inAngle, rShd, pShadowed)
     do iPosSample = 1, nPosSamples

        !! Find the observer-medium intersection point
        !!
        isFound = trc_findSurfaceIntersection(med%grid, [posSamples(:,iPosSample), med%hMean], -em_vectors(iEmAngle,:), isect)

        if(isFound) then
           emAngle   = acos(dot_product(em_vectors(iEmAngle, :), isect%N))

           do iInAngle = 1, size(inTheta)

              call ray_init(rShd, RAY_TYPE_SHADOW)
              rShd % P  = isect%P1
              rShd % D  = in_vectors(iInAngle, :)

              pShadowed = trc_traceOcclusion(med%grid, rShd)

              if(.not. pShadowed) then
                 InAngle   = acos(dot_product(in_vectors(iInAngle, :), isect%N))

                 xrfr(iInAngle, iEmAngle, :) = xrfr(iInAngle, iEmAngle, :)                 &
                      &                      + bfdf_Parviainen(mat, spc, inAngle, EmAngle)
               end if
           end do
        end if

        call utl_timerIncrease(timer)

     end do
     !$omp end parallel do
  end do
  !!
  !!--- END SIMULATION ---

  xrfr = xrfr / nPosSamples


!!$  !!--- SIMULATION ---
!!$  !!
!!$  do iInAngle = 1, size(inTheta)
!!$     do iPosSample = 1, nPosSamples
!!$
!!$        isFound = trc_findSurfaceIntersection(med%grid, [posSamples(iPosSample,:), med%hMean], in_vectors(iInAngle,:), isect)
!!$
!!$        if(isFound) then
!!$
!!$           InAngle   = acos(dot_product(-in_vectors(iInAngle, :), isect%N))
!!$
!!$           do iEmAngle = 1, nEmAngles
!!$
!!$              call ray_init(rEm, RAY_TYPE_SHADOW)
!!$              rEm % P = isect%P1
!!$              rEm % D = em_vectors(iEmAngle, :)
!!$
!!$              pVisible  = trc_traceOcclusion(med%grid, rEm)
!!$
!!$              if(pVisible) then
!!$                 emAngle   = acos(dot_product( em_vectors(iEmAngle, :), isect%N))
!!$
!!$                 xrfr(iInAngle, iEmAngle, :) = xrfr(iInAngle, iEmAngle, :)                 &
!!$                      &                      + bfdf_Parviainen(mat, spc, inAngle, EmAngle)
!!$               end if
!!$
!!$           end do
!!$        end if
!!$     end do
!!$  end do
  !!
  !!--- END SIMULATION ---

end subroutine simulation_method_analytic



!> Computes the fluorescent signal from a particulate medium plane using first-order analytic approximation.
!!
!! \param[in]  mat       Material
!! \param[in]  spc       Spectrum
!! \param[in]  inTheta   Angle(s) of incidence
!! \param[in]  emTheta   Angle(s) of emergence
!! \param[in]  outType   Output type
!!
!! \param[out] hs        Fluorescent signal for each line and incidence angle
!! \param[out] tb        Fluorescent signal for each line and incidence angle
!!
subroutine simulation_method_analytic_mc(mat, spc, nSpcSamples, inTheta, emTheta, outType, hs, tb)
  use base
  use material
  use geometry
  use hemisphere
  use particle
  use bfdf

  implicit none

  type(mat_material)       :: mat
  type(spc_spectrum)       :: spc
  real(fd), dimension(:)   :: inTheta
  real(fd), dimension(:)   :: emTheta
  integer                  :: outType
  integer                  :: nSpcSamples

  type(gth_hemisphere),       optional :: hs
  real(fd), dimension(:,:,:), optional :: tb

  real(fd), dimension(mat%nLines)      :: fluorescence, muLineTotal, fYield

  integer     :: iInAngle, iEmAngle, iLine

  ! Precompute fluorescence yield and total mu
  !
  do iLine = 1, mat%nLines
     if(mat%lEnergy(iLine) > 0.0_fd) then
        call mat_evalMu(mat, mat%lEnergy(iLine), muTotal = muLineTotal(iLine))
        fYield(iLine) = mat%fYield((iLine-1)/2 + 1)
     end if
  end do
  !!
  !!- Hemisphere branch
  !!
  if(outType == OUTPUT_TYPE_HEMISPHERE) then
     if(.not. present(hs)) then
        call utl_fatal_error("Fatal error in (xrfip.infslab_method_analytic): outputType set to 'Hemisphere', but no hemisphere given.")
     end if

     do iInAngle = 1, size(inTheta)
        do iEmAngle = 1, hs%resTheta
           fluorescence = bfdf_Parviainen_mc(mat, spc, inTheta(iInAngle), hs%mTheta(iEmAngle), nSpcSamples)

           do iLine = 1, mat%nLines
              call gth_hs_addToRow(hs, iEmAngle, fluorescence(iLine), iInAngle, iLine)
           end do
        end do
     end do
  end if
  !!
  !!- Tabulated branch
  !!
  if(outType == OUTPUT_TYPE_TABULATED) then
     if(.not. present(tb)) then
        call utl_fatal_error("Fatal error in (xrfip.infslab_method_analytic): outputType set to 'Tabulated', but no output table given.")
     end if

     do iInAngle = 1, size(inTheta)
        do iEmAngle = 1, size(emTheta)
           tb(iInAngle, iEmAngle, :) = bfdf_Parviainen_mc(mat, spc, inTheta(iInAngle), emTheta(iEmAngle), nSpcSamples, muLineTotal, fYield)
        end do
     end do
  end if

end subroutine simulation_method_analytic_mc
