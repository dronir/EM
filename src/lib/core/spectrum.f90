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
!! Copyright (C) 2007 - 2009 Hannu Parviainen
!!
!! Contributor(s): Hannu Parviainen
!!
!!--- END GPL ---

!>  Spectrum module defining the spc_spectrum type. 
!!
!! \date 02.11.2007
!!
module spectrum
  use base
  use util

  !>   Type defining a "continuous" spectrum.
  !!
  type :: spc_spectrum

     !>          Numnber of data points.
     integer  :: nPoints

     !>          Minimum energy level of the spectrum.
     real(FD) :: lMin

     !>          The width of energy bucket for continuous spectra.
     real(FD) :: dl

     !>          Weight factor.
     real(FD) :: weight

     !>          Type of the spectrum. This can have values:
     !!
     !!             SPC_TYPE_LINE
     !!             SPC_TYPE_SPECTRUM
     !!
     integer  :: type

     !>                                 Energy level table for discrete spectra.
     real(FD), dimension(:), pointer :: E   => null()

     !>                                 Intensity table.
     real(FD), dimension(:), pointer :: I   => null()

     !>                                 Spectrum cdf.
     real(FD), dimension(:), pointer :: CDF   => null()

     !>                                 Inverse of the spectrum cdf.
     real(FD), dimension(:), pointer :: ICDF   => null()
     real(FD), dimension(:), pointer :: ICDFX  => null()


  end type spc_spectrum

  integer, parameter :: SPC_TYPE_LINE       = 0
  integer, parameter :: SPC_TYPE_SPECTRUM   = 1

contains
  
  subroutine spc_spectrum_init(s, nPoints, E)
    type(spc_spectrum) :: s
    integer  :: nPoints, i
    real(FD), dimension(:), optional :: E


    if(present(E) .and. size(E) /= nPoints) then
       write(*,*) "Fatal error: E /= nPoints."
       stop
    end if

    allocate(s % E     (    nPoints ))
    allocate(s % I     (    nPoints ))
    allocate(s % CDF   (    nPoints ))
    allocate(s % ICDF  ( 10*nPoints ))
    allocate(s % ICDFX ( 10*nPoints ))


    s % nPoints = nPoints
    s % weight  = 0.0_fd
    s % I       = 0.0_fd
    s % E       = 0.0_fd
    s % type    = SPC_TYPE_SPECTRUM

    if(present(E)) s%E = E

    s % ICDFX = [(i*(1.0 / real(size(s%ICDFX)-1)) , i = 0, size(s%ICDFX)-1)]

  end subroutine spc_spectrum_init

  subroutine spc_spectrum_initline(s, E, I)
    type(spc_spectrum) :: s
    real(FD)           :: E
    real(FD), optional :: I

    call spc_spectrum_FREE(s)

    allocate(s % E(1))
    allocate(s % I(1))

    s % nPoints = 1
    s % E       = E
    s % type    = SPC_TYPE_LINE

    if(present(I)) then
       s % I       = I
    else
       s % I       = 1.0_fd
    end if

  end subroutine spc_spectrum_initline


  subroutine spc_spectrum_read(s, dFileName, eMin, eMax)
    type(spc_spectrum)          :: s    
    character(LEN=FNAME_LENGTH) :: dFileName
    real(FD)                    :: eMin, eMax

    real(FD), dimension(:,:), allocatable :: temp

    integer :: nPoints, iMin, iMax

    iMin = 1
    iMax = 1

    call spc_spectrum_free(s)

    open(1, FILE=dFileName, STATUS="OLD")
    read(1,*) nPoints
    
    allocate(temp(2,nPoints))
    read(1,*) temp
    close(1)

    !! Convert from keV to eV
    !!
    temp(1,:) = 1e3_fd * temp(1,:)

    do
       if (iMin > nPoints) call utl_fatal_error("Initializing spectrum - Spectrum energy minimum higher than maximum spectrum energy.")
       if ((temp(1,iMin) > eMin)) exit
       iMin = iMin + 1
    end do

    iMax = iMin
    do
       if (iMax > nPoints) then
          call utl_message("Initializing spectrum - Maximum energy greater than maximum spectrum energy, using the latter!")
          exit
       end if
       if ((temp(1,iMax+1) > eMax)) exit
       iMax = iMax + 1
    end do
 
    call spc_spectrum_init(s, iMax-iMin+1)

    s % type = SPC_TYPE_SPECTRUM
    s%E      = temp(1,iMin:iMax)
    s%I      = temp(2,iMin:iMax)

    s%I      = s%I / sum(s%I)
    s%lMin   = minval(s%E)

    !!- Compute the cumulative distribution function.
    !!
    s%CDF(1) = 0.0_fd
    do i=2, s%nPoints
       s%CDF(i) = s%CDF(i-1) + 0.5_fd * (s%I(i-1) + s%I(i)) * (s%E(i) - s%E(i-1))
    end do
    !!
    !!- Normalize spectrum cdf
    !!
    s%CDF = s%CDF / s%CDF(s%nPoints)
    !!
    !!- Invert spectrum cdf
    !!
    call utl_invert_cdf(s%E, s%cdf, s%icdf)


    deallocate(temp)

  end subroutine spc_spectrum_read


  real(FD) function spc_getSample(s, r)
    type(spc_spectrum)          :: s    
    real(fd)                    :: r

    if(s%type == SPC_TYPE_LINE) then
       spc_getSample = s%E(1)
    else
       !spc_getSample =  utl_lerp_lin_array(s%ICDF,r)
       spc_getSample =  utl_lerp_free_array(r, s%ICDFX, s%ICDF)
    end if

  end function spc_getSample

  real(FD) function spc_getLineSample(s)
    type(spc_spectrum)          :: s    

    spc_getLineSample = s%E(1)
  end function spc_getLineSample


  subroutine spc_spectrum_free(s)
    type(spc_spectrum) :: s

    if(associated(s%I)) deallocate(s%I)
    if(associated(s%E)) deallocate(s%E)

  end subroutine spc_spectrum_free

end module spectrum
