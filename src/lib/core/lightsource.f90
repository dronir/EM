!****h* XR/Lightsource
! NAME
!   Lightsource
!
! DESCRIPTION
!   Lightsource module defining the lsc_sun type. 
!
! NOTES
!
!
! AUTHOR
!   Hannu Parviainen
!
! USES
!   Base
!   Util
!   Spectrum
!
! CREATION DATE
!   08.11.2007
!******
MODULE LIGHTSOURCE
  USE BASE
  USE UTIL
  USE SPECTRUM

  !****s* Lightsource/lsc_sun
  ! NAME
  !   spc_spectrum
  !
  ! DESCRIPTION
  !   Type defining a "continuous" spectrum.
  !
  !******
  TYPE :: LSC_SUN
     !****d* lsc_sun/spectrum
     ! NAME
     !   spectrum
     !
     ! DESCRIPTION
     !   The spectrum of the sun.
     !
     !******
     type(spc_spectrum) :: spectrum

     !****d* lsc_sun/icdf
     ! NAME
     !   icdf
     !
     ! DESCRIPTION
     !   The inverse cumulative distribution function.
     !
     !******
     real(FD), DIMENSION(:), POINTER :: icdf

  END TYPE LSC_SUN

CONTAINS
  
  SUBROUTINE LSC_SPECTRUM_READ_CONTINUOUS(s, dFileName)
    type(LSC_SPECTRUM)          :: s    
    CHARACTER(LEN=FNAME_LENGTH) :: dFileName
    
    real(FD), DIMENSION(:,:), ALLOCATABLE :: temp

    INTEGER :: nPoints

    CALL LSC_SPECTRUM_FREE(s)

    OPEN(1, FILE=dFileName, STATUS="OLD")
    READ(1,*) nPoints
    
    ALLOCATE(temp(4, nPoints))
    READ(1,*) temp
    CLOSE(1)

    CALL LSC_SPECTRUM_INIT_CONTINUOUS(s, nPoints, temp(1,1), temp(1,2) - temp(1,1))

    s%I = temp(2,:)

    DEALLOCATE(temp)
  
  END SUBROUTINE LSC_SPECTRUM_READ_CONTINUOUS

  SUBROUTINE LSC_SPECTRUM_READ_DISCRETE(s, dFileName)
    type(LSC_SPECTRUM)          :: s    
    CHARACTER(LEN=FNAME_LENGTH) :: dFileName
    
    real(FD), DIMENSION(:,:), ALLOCATABLE :: temp

    INTEGER :: nPoints

    CALL LSC_SPECTRUM_FREE(s)

    OPEN(1, FILE=dFileName, STATUS="OLD")
    READ(1,*) nPoints

    CALL LSC_SPECTRUM_INIT_DISCRETE(s, nPoints)
    
    ALLOCATE(temp(4, nPoints))
    READ(1,*) temp
    CLOSE(1)

    s%E = temp(1,:)
    s%I = temp(2,:)

    s%lMin = MINVAL(s%E)

    DEALLOCATE(temp)
  
  END SUBROUTINE LSC_SPECTRUM_READ_DISCRETE


  !****d* lsc_spectrum/ICDF
  ! NAME
  !   ICDF
  !
  ! DESCRIPTION
  !   Inverse cumulative distribution function of the spectrum.
  !
  !******
  !     real(FD), DIMENSION(:), POINTER :: ICDF => NULL()
  FUNCTION LSC_SPECTRUM_COMPUTEICDF(s, nPoints) RESULT(ICDF)
    type(LSC_SPECTRUM) :: s
    INTEGER            :: nPoints, i

    real(FD), DIMENSION(SIZE(s%I)) :: CDF
    real(FD), DIMENSION(nPoints)   :: ICDF

    ! Compute the cumulative distribution function.
    !
    CDF(1) = s%I(1)
    DO i=2, s%nPoints
       CDF(i) = CDF(i-1) + s%I(i)
    END DO

    ! Normalization
    !
    CDF = CDF / CDF(s%nPoints)

    CALL utl_invertCDF(cdf, icdf, s%E(1), s%E(s%nPoints))

  END FUNCTION LSC_SPECTRUM_COMPUTEICDF

  SUBROUTINE LSC_SPECTRUM_FREE(s)
    type(LSC_SPECTRUM) :: s

    IF(ASSOCIATED(s%I)) DEALLOCATE(s%I)
    IF(ASSOCIATED(s%E)) DEALLOCATE(s%E)

  END SUBROUTINE LSC_SPECTRUM_FREE

END MODULE LIGHTSOURCE
