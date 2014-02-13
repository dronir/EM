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

!> Defines numerical types, mathematical and physical constants useful in computations, conversion
!! constants, and physical conversion functions.
!!
module base

  use, intrinsic :: ISO_C_BINDING

  implicit none

  integer,  parameter :: FS = C_FLOAT                     !< Single precision float
  integer,  parameter :: FD = C_DOUBLE                    !< Double precision float

  integer,  parameter :: IS = C_SHORT                     !< Short integer
  integer,  parameter :: IL = C_INT                       !< Long  integer

  ! Mathematical constants
  !
  real(FD), parameter :: PI         = 3.1415926535        !< \f$ \pi \f$
  real(FD), parameter :: INV_PI     = 1.0_fd / PI         !< \f$ \pi^{-1} \f$
  real(FD), parameter :: TWO_PI     = 2.0_fd * PI         !< \f$ 2\pi \f$
  real(FD), parameter :: INV_TWO_PI = 1.0_fd / TWO_PI     !< \f$ \frac{1}{2\pi} \f$
  real(FD), parameter :: FOUR_PI    = 4.0_fd * PI         !< \f$ 4\pi \f$
  real(FD), parameter :: INV_FOUR_PI= 1.0_fd / FOUR_PI    !< \f$ \frac{1}{4\pi} \f$
  real(FD), parameter :: HALF_PI    = 0.5_fd * PI         !< \f$ \frac{\pi}{2} \f$
  real(FD), parameter :: INV_HALF_PI= 1.0_fd / HALF_PI    !< \f$ \frac{2}{\pi} \f$

  !!- Trace constants.
  !!
  real(FD), parameter :: TRACE_EPS  = 1e-7_fd             !< Ray tracing epsilon
  real(FD), parameter :: TRACE_HUGE = huge(PI)            !< Ray tracing infinity


  !!- Output type 
  !!
  integer, parameter  :: OUTPUT_TYPE_HEMISPHERE = 0       !< Use hemisphere
  integer, parameter  :: OUTPUT_TYPE_TABULATED  = 1       !< Use tabulated emergence angle values

  !
  !
  integer, parameter  :: FNAME_LENGTH    = 200            !< Maximum length of a filename



  ! Physical constants.
  !
  real(FD), parameter :: PHS_C    = 2.99792458e8_fd        !< Speed of light [m / s]
  real(FD), parameter :: PHS_H_J  = 6.62606896e-34_fd      !< Planck's constant [J s]
  real(FD), parameter :: PHS_HL_J = PHS_H_J / TWO_PI       !< Planck constant over 2 pi [J s]
  real(FD), parameter :: PHS_R_E  = 2.8179402894e-15_fd    !< Classical electron radius [m]

  real(FD), parameter :: PHS_M_E  = 9.109389e-31_fd        !< Electron rest mass [kg]

  real(FD), parameter :: PHS_NA = 6.02214179e23_fd         !< Avogardo constant [1 / mol]
  real(FD), parameter :: PHS_MU = 1.660538782e-27_fd       !< Atomic mass unit [kg]

  real(FD), parameter :: micron = 1.0e-6_fd                !< Micron [m]


  ! Physical conversion factors.
  !
  real(FD), parameter :: eV_keV    = 1e-3_fd               !< [eV]  to [keV]
  real(FD), parameter :: keV_eV    = 1e3_fd                !< [keV] to [eV]

  real(FD), parameter :: eV_J      = 1.602176462e-19_fd    !< [eV] to [J]
  real(FD), parameter :: J_eV      = 1.0_fd / eV_J         !< [J]  to [eV]

  real(FD), parameter :: gcm3_kgm3 = 1e3_fd                !< [g/cm^3] to [kg/m^3]
  real(FD), parameter :: kgm3_gcm3 = 1e-3_fd               !< [kg/m^3] to [g/cm^3]

  real(FD), parameter :: barn_to_cm2  = 1e-24_fd           !< [barn] to [cm^2]
  real(FD), parameter :: barn_to_m2   = 1e-28_fd           !< [barn] to [m^2]

  real(FD), parameter :: CNV_DEGTORAD    = PI / 180.0_fd   !< Degrees to radians
  real(FD), parameter :: CNV_RADTODEG    = 180.0_fd / PI   !< Radians to degrees

  real(FD), parameter :: d_r             = PI / 180.0_fd   !< Degrees to radians
  real(FD), parameter :: r_d             = 180.0_fd / PI   !< Radians to degrees


  ! http://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html

  real(FD), parameter :: ELM_NA_CA = 1.550e-3_fd / (PHS_MU * 20.0_fd)
  real(FD), parameter :: ELM_NA_FE = 7.874e-3_fd / (PHS_MU * 26.0_fd)
  

contains
  
  !> Energy [eV] to wavelength [m] conversion
  !!
  pure real(FD) function eToLambda(e)
    real(FD), intent(IN) :: e
    eToLambda = PHS_H_J * PHS_C / (e * eV_J)
  end function eToLambda
 
  !> Wavelength [m] to energy [eV] conversion
  !!
  pure real(FD) function lambdaToE_eV(l)
    real(FD), intent(IN) :: l
    lambdaToE_eV = J_eV * PHS_H_J * PHS_C / l
  end function lambdaToE_eV

  !> Wavelength [m] to energy [J] conversion
  !!
  pure real(FD) function lambdaToE_J(l)
    real(FD), intent(IN) :: l
    lambdaToE_J = PHS_H_J * PHS_C / l
  end function lambdaToE_J

end module base
