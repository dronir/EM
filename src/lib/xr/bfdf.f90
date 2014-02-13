!****h* XR/lib/xr/bfdf.f90
! NAME
!   bfdf
!
! DESCRIPTION
!   Analytic bidirectional fluorescence distribution functions
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   19.08.2008
!******

module bfdf
  use base
  use material
  use spectrum

  implicit none

contains
  
  function bfdf_Parviainen(mat, spc, theta_i, theta_e, muLine, fluorYield) result(fluorescence)
    type(mat_material) :: mat
    type(spc_spectrum) :: spc
    real(fd)           :: theta_i
    real(fd)           :: theta_e
    real(fd)           :: fluorescence(mat%nLines)

    real(fd), optional :: muLine(:), fluorYield(:)

    real(fd)           :: muLineTotal(mat%nLines), fYield(mat%nLines), line_fraction(mat%nLines)
    real(fd)           :: muInTotal, muInPhoto, muInFluor
    real(fd)           :: fluor_output, cos_theta_i, cos_theta_e

    integer            :: iLine, iSpec

    fluorescence = 0.0_fd

    if(present(muLine) .and. present(fluorYield)) then
       muLineTotal = muLine
       fYield      = fluorYield
    else
       do iLine = 1, mat%nLines
          if(mat%lEnergy(iLine) > 0.0_fd) then
             call mat_evalMu(mat, mat%lEnergy(iLine), muTotal = muLineTotal(iLine))
             fYield(iLine) = mat%fYield((iLine-1)/2 + 1)
          end if
       end do
    end if

    do iSpec = 1, spc%nPoints
       
       call mat_evalMu(mat, spc%E(iSpec), muPhotoIon       = muInPhoto, &
            &                             muFluorLineTotal = muInFluor, &
            &                             muTotal          = muInTotal  ) 

       fluor_output     = spc%I(iSpec) * INV_FOUR_PI * muInFluor / muInTotal
       line_fraction    = mat_evalFluorLineFractions(mat, spc%E(iSpec)) * fYield(:)

       fluorescence(:) = fluorescence(:)                                &
            & + fluor_output * line_fraction(:)                         &
            & * cos(theta_e) / (cos(theta_e)  +  muLineTotal(:)/muInTotal * cos(theta_i))
    end do

    fluorescence = fluorescence / sum(spc%I)

  end function bfdf_Parviainen


  function bfdf_Parviainen_mc(mat, spc, theta_i, theta_e, nSpcSamples, muLine, fluorYield) result(fluorescence)
    type(mat_material) :: mat
    type(spc_spectrum) :: spc
    real(fd)           :: theta_i
    real(fd)           :: theta_e
    real(fd)           :: fluorescence(mat%nLines)
    integer            :: nSpcSamples

    real(fd), optional :: muLine(:), fluorYield(:)

    real(fd)           :: I, energy, e
    real(fd)           :: muLineTotal(mat%nLines), fYield(mat%nLines), line_fraction(mat%nLines)
    real(fd)           :: muInTotal, muInPhoto, muInFluor
    real(fd)           :: fluor_output, cos_theta_i, cos_theta_e

    integer            :: iLine, iSpec

    fluorescence = 0.0_fd

    if(present(muLine) .and. present(fluorYield)) then
       muLineTotal = muLine
       fYield      = fluorYield
    else
       do iLine = 1, mat%nLines
          if(mat%lEnergy(iLine) > 0.0_fd) then
             call mat_evalMu(mat, mat%lEnergy(iLine), muTotal = muLineTotal(iLine))
             fYield(iLine) = mat%fYield((iLine-1)/2 + 1)
          end if
       end do
    end if

    do iSpec = 1, nSpcSamples

       call random_number(e)

       energy = spc_getSample(spc, e)
       I      = 1.0_fd / real(nSpcSamples, fd) 

       !  utl_lerp_lin_array(spc%I,(energy-spc%E(1)) / (spc%E(spc%nPoints) - spc%E(1)) )
       
       call mat_evalMu(mat, energy, muPhotoIon       = muInPhoto, &
            &                       muFluorLineTotal = muInFluor, &
            &                       muTotal          = muInTotal  ) 

       fluor_output     = I * INV_FOUR_PI * muInFluor / muInTotal
       line_fraction    = mat_evalFluorLineFractions(mat, energy) * fYield(:)

       fluorescence(:) = fluorescence(:)                                &
            & + fluor_output * line_fraction(:)                         &
            & * cos(theta_e) / (cos(theta_e)  +  muLineTotal(:)/muInTotal * cos(theta_i))
    end do

  end function bfdf_Parviainen_mc

!!$  subroutine bfdf_Brunetti(mat, s, theta_i, theta_e, results)
!!$    type(mat_material)       :: mat
!!$    type(spc_spectrum)       :: s
!!$    real(fd), dimension(:)   :: theta_i
!!$    real(fd)                 :: theta_e
!!$    real(fd), dimension(:,:) :: results
!!$
!!$    real(fd)                 :: energyIn, intensityIn
!!$
!!$    real(fd)    :: muInTotal, muEmTotal, muInPhoto, muInFluor, muLineExt(mat%nLines), fYield(mat%nLines)
!!$    real(fd)    :: A, B, fOut, xEn, cos_theta_i, cos_theta_e, r(mat%nLines)
!!$
!!$    integer :: iTht, iThtE, iLine, iSpec, iEn
!!$
!!$    do iLine = 1, mat%nLines
!!$       if(mat%lEnergy(iLine) > 0.0_fd) then
!!$          call mat_evalMu(mat, mat%lEnergy(iLine), muTotal = muLineExt(iLine))
!!$          fYield(iLine) = mat%fYield((iLine-1)/2 + 1)
!!$          r             = (JumpFactor(mat%Z((iLine-1)/2 + 1), K_SHELL) - 1.0) / JumpFactor(mat%Z((iLine-1)/2 + 1), K_SHELL)
!!$       end if
!!$    end do
!!$
!!$    cos_theta_e = cos(theta_e)
!!$
!!$    do iTht = 1, size(theta_i)
!!$       cos_theta_i = cos(theta_i(iTht))
!!$
!!$       do iSpec = 1, s%nPoints
!!$          energyIn    = s%E(iSpec)
!!$          IntensityIn = s%I(iSpec)
!!$
!!$          call mat_evalMu(mat, energyIn, muPhotoIon = muInPhoto, muFluorLineTotal = muInFluor, muTotal = muInTotal, iEn=iEn, xEn=xEn)
!!$
!!$          do iLine = 1, mat%nLines
!!$             if(mat%lEnergy(iLine) > 0.0) then
!!$                
!!$                !A = IntensityIn * obsSolidAngle * INV_FOUR_PI * muInFluor / ( muInTotal * cos_theta_i)
!!$                !B = ((1.0_fd - xEn) * mat%muFluorLine(iLine,iEn) + xEn * mat%muFluorLine(iLine,iEn+1)) * fYield(iLine)
!!$
!!$                A = intensityIn * fYield(iLine) * obsSolidAngle * INV_FOUR_PI * muInPhoto / muInTotal
!!$                B = ((1.0_fd - xEn) * mat%muFluorLine(iLine,iEn) + xEn * mat%muFluorLine(iLine,iEn+1)) / muInFluor
!!$
!!$                !results(iTht, iLine) = results(iTht, iLine) + A * B / (( muInTotal/c_thtI) + (muLineExt(iLine)/cos(theta_e)))
!!$                !results(iTht, iLine) = results(iTht, iLine) + A * B / ( muLineExt(iLine) / muInTotal * cos_theta_i + cos_theta_e)
!!$                results(iTht, iLine) = results(iTht, iLine) + A * B * cos_theta_e / ( cos_theta_e + muLineExt(iLine) / muInTotal * cos_theta_i)
!!$
!!$             end if
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    results = results / sum(s%I)
!!$
!!$  end subroutine bfdf_Brunetti
!!$
!!$
!!$ subroutine bfdf_Carpenter(mat, s, theta_i, theta_e, results)
!!$    type(mat_material)       :: mat
!!$    type(spc_spectrum)       :: s
!!$    real(fd), dimension(:)   :: theta_i
!!$    real(fd)                 :: theta_e
!!$    real(fd), dimension(:,:) :: results
!!$
!!$    real(fd)                 :: energyIn, intensityIn
!!$
!!$    real(fd)    :: muInTotal, muPhotoIon, muEmTotal, muInFluor, muLineExt(mat%nLines), fYield(mat%nLines)
!!$    real(fd)    :: A, muFluorLine, fOut, xEn, cos_theta_i, cos_theta_e
!!$
!!$    integer :: iTht, iThtE, iLine, iSpec, iEn
!!$
!!$    do iLine = 1, mat%nLines
!!$       if(mat%lEnergy(iLine) > 0.0_fd) then
!!$          call mat_evalMu(mat, mat%lEnergy(iLine), muTotal = muLineExt(iLine))
!!$          fYield(iLine) = mat%fYield((iLine-1)/2 + 1)
!!$       end if
!!$    end do
!!$
!!$    A = INV_FOUR_PI 
!!$ 
!!$    cos_theta_e = cos(theta_e)
!!$
!!$    do iTht = 1, size(theta_i)
!!$       cos_theta_i = cos(theta_i(iTht))
!!$
!!$       do iSpec = 1, s%nPoints
!!$          energyIn    = s%E(iSpec)
!!$          IntensityIn = s%I(iSpec)
!!$
!!$          do iLine = 1, mat%nLines
!!$             if(mat%lEnergy(iLine) > 0.0) then
!!$                
!!$                call mat_evalMu(mat, energyIn, muPhotoIon = muPhotoIon, muFluorLineTotal = muInFluor, muTotal = muInTotal, iEn=iEn, xEn=xEn)
!!$
!!$                !A = IntensityIn * obsSolidAngle * INV_FOUR_PI * muInFluor / muInTotal
!!$                A = IntensityIn * obsSolidAngle * INV_FOUR_PI * muInFluor / muPhotoIon
!!$
!!$                muFluorLine = ((1.0_fd - xEn) * mat%muFluorLine(iLine,iEn) + xEn * mat%muFluorLine(iLine,iEn+1)) * fYield(iLine)
!!$                
!!$                results(iTht, iLine) = results(iTht, iLine) + A * muFluorLine / ( muInTotal + muLineExt(iLine) * cos_theta_i / cos_theta_e )
!!$             end if
!!$          end do
!!$
!!$       end do
!!$    end do
!!$
!!$    results = results / sum(s%I)
!!$
!!$  end subroutine bfdf_Carpenter


end module bfdf
  
