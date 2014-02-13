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


!>   Computes the fluorescent signal from an infinite plane using first-order analytic approximation.
!!
!!   \param[in]  mat       Material
!!   \param[in]  spc       Spectrum
!!   \param[in]  theta_i   Angle(s) of incidence
!!   \param[in]  theta_e   Angle(s) of emergence
!!   \param[in]  outType   Output type
!!
!!   \param[out] hs        Fluorescent signal for each line and incidence angle
!!   \param[out] tb        Fluorescent signal for each line and incidence angle
!!
subroutine infslab_method_analytic(mat, spc, theta_i, theta_e, outType, hs, tb)
  use base
  use material
  use geometry
  use hemisphere
  use bfdf

  implicit none

  type(mat_material)       :: mat
  type(spc_spectrum)       :: spc
  real(fd), dimension(:)   :: theta_i
  real(fd), dimension(:)   :: theta_e
  integer                  :: outType

  type(gth_hemisphere),       optional :: hs
  real(fd), dimension(:,:,:), optional :: tb

  real(fd), dimension(mat%nLines)      :: fluorescence, muLineTotal, fYield

  integer     :: iThtI, iThtE, iLine

  ! Precompute fluorescence yield and total mu
  !
  do iLine = 1, mat%nLines
     if(mat%lEnergy(iLine) > 0.0_fd) then
        call mat_evalMu(mat, mat%lEnergy(iLine), muTotal = muLineTotal(iLine))
        fYield(iLine) = mat%fYield((iLine-1)/2 + 1)
     end if
  end do

  ! Hemisphere branch
  !
  if(outType == OUTPUT_TYPE_HEMISPHERE) then
     if(.not. present(hs)) then
        call utl_fatal_error("Fatal error in (xrfip.infslab_method_analytic): outputType set to 'Hemisphere', but no hemisphere given.")
     end if

     do iThtI = 1, size(theta_i)
        do iThtE = 1, hs%resTheta
           fluorescence = bfdf_Parviainen(mat, spc, theta_i(iThtI), hs%mTheta(iThtE))

           do iLine = 1, mat%nLines
              call gth_hs_addToRow(hs, iThtE, fluorescence(iLine), iThtI, iLine)
           end do
        end do
     end do
  end if

  ! Tabulated branch
  !
  if(outType == OUTPUT_TYPE_TABULATED) then
     if(.not. present(tb)) then
        call utl_fatal_error("Fatal error in (xrfip.infslab_method_analytic): outputType set to 'Tabulated', but no output table given.")
     end if

     do iThtI = 1, size(theta_i)
        do iThtE = 1, size(theta_e)
           tb(iThtI, iThtE, :) = bfdf_Parviainen(mat, spc, theta_i(iThtI), theta_e(iThtE), muLineTotal, fYield)
        end do
     end do
  end if

end subroutine infslab_method_analytic


!>   Computes the fluorescent signal from an infinite plane using first-order analytic approximation.
!!
!!   \param[in]  mat          Material
!!   \param[in]  spc          Spectrum
!!   \param[in]  nSpcSamplse  Number of spectrum samples
!!   \param[in]  theta_i      Angle(s) of incidence
!!   \param[in]  theta_e      Angle(s) of emergence
!!   \param[in]  outType      Output type
!!
!!   \param[out] hs           Fluorescent signal for each line and incidence angle
!!   \param[out] tb           Fluorescent signal for each line and incidence angle
!!
subroutine infslab_method_analytic_mc(mat, spc, nSpcSamples, theta_i, theta_e, outType, hs, tb)
  use base
  use material
  use geometry
  use hemisphere
  use bfdf

  implicit none

  type(mat_material)       :: mat
  type(spc_spectrum)       :: spc
  real(fd), dimension(:)   :: theta_i
  real(fd), dimension(:)   :: theta_e
  integer                  :: outType
  integer                  :: nSpcSamples

  type(gth_hemisphere),       optional :: hs
  real(fd), dimension(:,:,:), optional :: tb

  real(fd), dimension(mat%nLines)      :: fluorescence, muLineTotal, fYield

  integer     :: iThtI, iThtE, iLine

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

     do iThtI = 1, size(theta_i)
        do iThtE = 1, hs%resTheta
           fluorescence = bfdf_Parviainen_mc(mat, spc, theta_i(iThtI), hs%mTheta(iThtE), nSpcSamples)

           do iLine = 1, mat%nLines
              call gth_hs_addToRow(hs, iThtE, fluorescence(iLine), iThtI, iLine)
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

     do iThtI = 1, size(theta_i)
        do iThtE = 1, size(theta_e)
           tb(iThtI, iThtE, :) = bfdf_Parviainen_mc(mat, spc, theta_i(iThtI), theta_e(iThtE), nSpcSamples, muLineTotal, fYield)
        end do
     end do
  end if

end subroutine infslab_method_analytic_mc
