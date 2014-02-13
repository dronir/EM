!****h* XR/hemisphere
! NAME
!   hemisphere
!
! DESCRIPTION
!
!           () -> []
!           () -> [][]
!     ^     () -> [][][]      
!     |     () -> [][][][][]  
!  theta_e  () -> [][][][][][][][]
!           () -> [][][][][][][][][][][]
!
!           0      phi_e ->   2 PI
!
!           [] = 
!
! NOTES
!
!
! AUTHOR
!   Hannu Parviainen
!
! USES
!   Base
!
!
! CREATION DATE
!   19.09.2007
!******
MODULE hemisphere
  USE base
  USE netcdf
  USE spectrum

  IMPLICIT NONE

  TYPE :: gth_cellTable
     real(FD) :: dPhi, dPhiInv, dA, mTheta
     INTEGER  :: res_phi

     real(FD), DIMENSION(:,:), POINTER :: cells
     real(FD), DIMENSION(:), POINTER :: weight
  END TYPE gth_cellTable

  TYPE :: gth_hemisphere
     INTEGER  :: res_theta
     real(FD) :: dTheta, dThetaInv

     real(FD) :: dA_mean
     INTEGER  :: nCells, nLevels

     character(len=36) :: name
     type(gth_cellTable), DIMENSION(:), POINTER :: rows
  END TYPE gth_hemisphere

  type :: gth_hsFile
     integer                 :: fileID  
     integer, dimension(100) :: dataID
     integer                 :: nDatasets
     logical                 :: isDef
  end type gth_hsFile

CONTAINS

  SUBROUTINE gth_hemisphere_init(h, rT, nDataLevels)
    type(gth_hemisphere) :: h
    INTEGER              :: rT
    integer, optional    :: nDataLevels

    real(FD)             :: dPhi, dA0
    integer              :: rPhi, i, j

    h % res_theta = rT
    h % dTheta    = 0.5_fd * PI / real(rT)
    h % dThetaInv = 1.0_fd / h % dTheta

    h % nCells    = nDataLevels
    h % nLevels   = nDataLevels

    ALLOCATE(h%rows(rT))

    ! Compute the solid angle of the [0 .. dTheta] row (single cell).
    ! Also, initialize the first row.

    dA0 = TWO_PI * (1.0_fd - COS(h%dTheta))

    h % rows(1) % res_phi = 1
    h % rows(1) % dPhi    = TWO_PI
    h % rows(1) % dPhiInv = INV_TWO_PI
    h % rows(1) % dA      = dA0

    ALLOCATE(h % rows(1) % cells(1,nDataLevels))
    ALLOCATE(h % rows(1) % weight(1))

    DO i = 2, rT
       dPhi = dA0 / (COS((i-1) * h%dTheta) - COS(i * h%dTheta))
       rPhi = NINT(TWO_PI / dPhi)
       dPhi = TWO_PI / real(rPhi)

       h % rows(i) % res_phi = rPhi
       h % rows(i) % dPhi    = dPhi
       h % rows(i) % dPhiInv = 1.0_fd / dPhi
       h % rows(i) % dA      = dPhi * (COS((i-1) * h%dTheta) - COS(i * h%dTheta))

       h % nCells = h % nCells + rPhi * nDataLevels

       ALLOCATE(h % rows(i) % cells(rPhi, nDataLevels))
       ALLOCATE(h % rows(i) % weight(rPhi))
    END DO

    DO i = 1, rT
       h % rows(i) % cells     = 0.0_fd
       h % rows(i) % weight    = 0.0_fd
       h % rows(i) % mTheta    = h % dTheta * (real(i, kind=fd) - 0.5_fd)
    END DO

    h % dA_mean = TWO_PI / real(h%nCells, KIND=FD)

  END SUBROUTINE gth_hemisphere_init

  !****f* hemisphere/gth_sphDirToCell
  ! NAME
  !   gth_sphDirToCell
  !
  ! DESCRIPTION
  !   
  !
  ! INPUTS
  !   h : hemisphere ... 
  !   p : real(2)    ... 
  !
  ! SOURCE

  PURE FUNCTION gth_sphDirToCell(h, P) RESULT(C)
    type(gth_hemisphere), INTENT(IN) :: h
    real(FD), DIMENSION(2), INTENT(IN)  :: P
    INTEGER,  DIMENSION(2) :: C
    
    C(1) = FLOOR( P(1) * h % dThetaInv            ) + 1
    C(2) = FLOOR( P(2) * h % rows(C(1)) % dPhiInv ) + 1

  END FUNCTION gth_sphDirToCell
  !******

  !****f*  hemisphere/gth_carDirToCell
  ! NAME
  !   gth_carDirToCell
  !
  ! DESCRIPTION
  !   
  !
  ! INPUTS
  !   h : hemisphere ... 
  !   p : real(3)    ... 
  !
  ! SOURCE
  PURE FUNCTION gth_carDirToCell(h, P) RESULT(C)
    type(gth_hemisphere), INTENT(IN) :: h
    real(FD), DIMENSION(3), INTENT(IN)  :: P
    INTEGER,  DIMENSION(2) :: C

    real(FD) :: theta, phi, r

    r = SQRT(SUM(P**2))

    theta = ACOS(P(3) / r);
    phi   = ATAN2(P(2) / r, P(1) / r);
    
    IF(phi < 0) phi = TWO_PI + phi;
    
    C(1) = FLOOR( theta * h % dThetaInv            ) + 1
    C(2) = FLOOR( phi   * h % rows(C(1)) % dPhiInv ) + 1

  END FUNCTION gth_carDirToCell
  !******

  PURE FUNCTION gth_cellCenterCar(h, t, p) RESULT(d)
    type(gth_hemisphere), INTENT(IN) :: h
    INTEGER, INTENT(IN) :: t, p

    REAL, DIMENSION(3)   :: d
    REAL :: theta, phi

    theta = (real(t) - 0.5_fd) * h % dTheta
    phi   = (real(p) - 0.5_fd) * h % rows(t) % dPhi
    
    d(1) = SIN(theta) * COS(phi);
    d(2) = SIN(theta) * SIN(phi);
    d(3) = COS(theta);

  END FUNCTION gth_cellCenterCar

  FUNCTION gth_cellRandomSampleCar(h, t, p) RESULT(d)
    type(gth_hemisphere), INTENT(IN) :: h
    INTEGER, INTENT(IN) :: t, p

    real(FD), DIMENSION(3)   :: d
    real(FD) :: theta, phi, r(2)

    CALL RANDOM_NUMBER(r)

    theta = (real(t,fd) - r(1)) * h % dTheta
    phi   = (real(p,fd) - r(2)) * h % rows(t) % dPhi
    
    d(1) = SIN(theta) * COS(phi);
    d(2) = SIN(theta) * SIN(phi);
    d(3) = COS(theta);

  END FUNCTION gth_cellRandomSampleCar

  subroutine gth_hs_addToRow(h,r,v)
    type(gth_hemisphere),   INTENT(INOUT)  :: h
    INTEGER                                :: r
    real(FD)                               :: v

    integer :: i

    do i=1, h%rows(r) % res_phi
       h % rows(r) % cells(i,1) = h % rows(r) % cells(i,1) + v
    end do

  end subroutine gth_hs_addToRow

  SUBROUTINE gth_addDataCar(h, d, v)
    type(gth_hemisphere),   INTENT(INOUT)  :: h
    real(FD), DIMENSION(3), INTENT(IN)     :: d
    real(FD),               INTENT(IN)     :: v

    INTEGER,  DIMENSION(2) :: C

    C = gth_carDirToCell(h, d)

    h % rows(C(1)) % cells(C(2), 1) = h % rows(C(1)) % cells(C(2), 1) + v 

  END SUBROUTINE gth_addDataCar

  SUBROUTINE gth_distributeValue(h, v)
    type(gth_hemisphere),   INTENT(INOUT)  :: h
    real(FD),               INTENT(IN)     :: v

    INTEGER  :: i, j
    real(FD) :: vDist

    DO i=1, h%res_theta

       vDist = v * INV_TWO_PI !* (h%dA_mean / h%rows(i)%dA)

       DO j=1, h%rows(i)%res_phi
          h%rows(i)%cells(j,1) = h%rows(i)%cells(j,1) + vDist
       END DO
    END DO

  END SUBROUTINE gth_distributeValue

  real(FD) FUNCTION gth_hemisphere_integrate(h, level) RESULT(sum)
    type(gth_hemisphere),   INTENT(IN)  :: h
    integer, optional                   :: level

    INTEGER :: i, j, lvl

    lvl = 1
    if(present(level)) lvl = level

    sum = 0.0_fd

    DO i=1, h%res_theta
       DO j=1, h%rows(i)%res_phi
          sum = sum + h%rows(i)%dA * h%rows(i)%cells(j, lvl)
       END DO
    END DO

  END FUNCTION gth_hemisphere_integrate


!!$  SUBROUTINE gth_hemisphere_normalize(h)
!!$    type(gth_hemisphere), INTENT(INOUT)  :: h
!!$
!!$    INTEGER  :: i, j
!!$    real(FD) :: nFact
!!$
!!$    DO i=1, h%res_theta
!!$       DO j = 1, h % rows(i) % res_phi
!!$          h%rows(i)%cells(j,:) = h%rows(i)%cells(j,:) * (h%dA_mean / h%rows(i)%dA)
!!$       END DO
!!$    END DO
!!$
!!$  END SUBROUTINE gth_hemisphere_normalize


  subroutine gth_hsFileOpen(hf, fName, mode)
    type(gth_hsFile), intent(INOUT) :: hf
    character(len=*), intent(IN)    :: fName
    character(len=1), intent(IN)    :: mode

    if(mode == 'r') then
       call gth_nfCheck( nf90_open(fName, NF90_NOWRITE, hf%fileID) )
    else
       call gth_nfCheck( nf90_create(fName, NF90_CLOBBER, hf%fileID) )
    end if
    
    hf%dataID    = 0
    hf%nDatasets = 0
    hf%isDef     = .true.
    
  end subroutine gth_hsFileOpen

  subroutine gth_hsFileClose(hf)
    type(gth_hsFile), intent(INOUT) :: hf

    call gth_nfCheck( nf90_close(hf%fileID) )
    hf%fileID = 0
  end subroutine gth_hsFileClose

  subroutine gth_hsWriteHeader(hf, h, auth, prog, vers)
    type(gth_hsFile),     intent(INOUT) :: hf
    type(gth_hemisphere), INTENT(IN), dimension(:) :: h

    integer                             :: nHemis

    character(len=*) :: auth, prog, vers

    integer,  dimension(h(1)%res_theta)  :: nPhi
    real(fd), dimension(h(1)%res_theta)  :: dPhi, dA, mTheta

    integer  :: dimID(1)
    integer  :: i, j, k
    
    k = 1
    do i=1, h(1)%res_theta

       nPhi(i)   = h(1)%rows(i)%res_phi
       dPhi(i)   = h(1)%rows(i)%dPhi
       dA(i)     = h(1)%rows(i)%dA
       mTheta(i) = h(1)%rows(i)%mTheta

    end do

    call gth_nfCheck( nf90_def_dim(hf%fileID, "nData", h(1)%nCells, dimID(1)) )

    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "Author",  auth ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "Program", prog ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "Version", vers ) )

    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "nDatasets",  (/size(h)/)    ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "nCells",  (/h(1)%nCells/)    ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "nTheta",  (/h(1)%res_theta/) ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "nLevels", (/h(1)%nLevels/)   ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "nPhi",    nPhi            ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "dTheta",  (/h(1)%dTheta/)    ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "dPhi",    dPhi            ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "dA_mean", (/h(1)%dA_mean/)   ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "dA",      dA              ) )
    call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "mTheta",  mTheta          ) )

    do i=1,size(h)
       write(*,*) h(i)%name
       call gth_nfCheck( nf90_def_var(hf%fileID, h(i)%name, NF90_DOUBLE, dimID, hf%dataID(i)) )
    end do

    hf%nDatasets = size(h)

  end subroutine gth_hsWriteHeader

  subroutine gth_hsAddAttS(hf, varIdx, attName, attData)
    type(gth_hsFile),     intent(INOUT) :: hf
    integer                             :: varIdx
    character(len=*)                    :: attName
    character(len=*)                    :: attData

    integer                             :: vID
    
    if(.not. hf%isDef) call gth_nfCheck( nf90_redef(hf%fileID) )

    vID = varIdx
    if(vID /= NF90_GLOBAL) vID = hf%dataID(varIdx)

    call gth_nfCheck( nf90_put_att(hf%fileID, vID, attName, attData ) )

  end subroutine gth_hsAddAttS

  subroutine gth_hsAddAttF(hf, varIdx, attName, attData)
    type(gth_hsFile),     intent(INOUT) :: hf
    integer                             :: varIdx
    character(len=*)                    :: attName
    real(fd), dimension(:)              :: attData
    
    integer                             :: vID
    
    if(.not. hf%isDef) call gth_nfCheck( nf90_redef(hf%fileID) )

    vID = varIdx
    if(vID /= NF90_GLOBAL) vID = hf%dataID(varIdx)

    call gth_nfCheck( nf90_put_att(hf%fileID, vID, attName, attData ) )

  end subroutine gth_hsAddAttF


  subroutine gth_hsSaveData(h, hf, norm)
    type(gth_hemisphere), INTENT(IN), dimension(:)  :: h
    type(gth_hsFile), intent(INOUT) :: hf
  
    real(fd), optional :: norm

    real(fd), dimension(h(1)%nCells)     :: dataOut
    integer  :: i, j, k, setIdx

    if(hf%isDef) call gth_nfCheck( nf90_enddef(hf%fileID) )
    hf%isDef = .false.

    do setIdx = 1, hf%nDataSets 
       k = 1
       do i = 1, h(setIdx)%res_theta
          do j = 1, h(setIdx)%rows(i)%res_phi
             dataOut(k:k+h(setIdx)%nLevels-1) = h(setIdx)%rows(i)%cells(j,:)
             k = k + h(setIdx)%nLevels
          end do
       end do

       if(present(norm)) dataOut = dataOut * norm

       call gth_nfCheck( nf90_put_var(hf%fileID, hf%dataID(setIdx), dataOut) )
    end do
    
  end subroutine gth_hsSaveData

  SUBROUTINE gth_hemisphere_save(h, fName)
    type(gth_hemisphere), INTENT(IN)  :: h
    character (len = *)               :: fName

    real(fd) :: dataOut(h%nCells)
    integer,  dimension(h%res_theta)  :: nPhi
    real(fd), dimension(h%res_theta)  :: dPhi, dA, mTheta

    integer  :: fileID, dataID, dimID(1)
    integer  :: i, j, k
    
    k = 1
    do i=1, h%res_theta

       nPhi(i)   = h%rows(i)%res_phi
       dPhi(i)   = h%rows(i)%dPhi
       dA(i)     = h%rows(i)%dA
       mTheta(i) = h%rows(i)%mTheta

       do j=1, h%rows(i)%res_phi
          dataOut(k) = h%rows(i)%cells(j,1)
          k = k+1
       end do
    end do

    call check( nf90_create (fName, NF90_CLOBBER,  fileID) )
    call check( nf90_def_dim(fileID, "Cells", h%nCells, dimID(1)) )
    call check( nf90_def_var(fileID, "Hemisphere", NF90_DOUBLE, dimID, dataID) )

    call check( nf90_put_att(fileID, dataID, "nCells",  (/h%nCells/)    ) )
    call check( nf90_put_att(fileID, dataID, "nTheta",  (/h%res_theta/) ) )
    call check( nf90_put_att(fileID, dataID, "nPhi",    nPhi            ) )
    call check( nf90_put_att(fileID, dataID, "dTheta",  (/h%dTheta/)    ) )
    call check( nf90_put_att(fileID, dataID, "dPhi",    dPhi            ) )
    call check( nf90_put_att(fileID, dataID, "dA_mean", (/h%dA_mean/)   ) )
    call check( nf90_put_att(fileID, dataID, "dA",      dA              ) )
    call check( nf90_put_att(fileID, dataID, "mTheta",  mTheta          ) )

    call check( nf90_enddef (fileID) )

    call check( nf90_put_var(fileID, dataID, dataOut) )

    call check( nf90_close  (fileID) )

  contains
      subroutine check(status)
        integer, intent ( in) :: status

        if(status /= nf90_noerr) then
           print *, trim(nf90_strerror(status))
           stop "Stopped"
        end if
      end subroutine check
  
  END SUBROUTINE gth_hemisphere_save
  

  SUBROUTINE gth_hemisphere_plot_data(h, fileName, level)
    type(gth_hemisphere) :: h
    CHARACTER(LEN=*) :: fileName
    integer, optional :: level
    INTEGER      :: i, j, lvl
    REAL         :: margin, scale, cWidth, cHeight
 
    margin  =    10.0
    scale   =  1000.0

    lvl = 1
    if(present(level)) lvl = level

    OPEN(1, FILE=fileName, status="REPLACE")
    
    WRITE(1,'(A)') "%!PS-Adobe-2.0 EPSF-1.2"
    WRITE(1,'(A)') "%%Title: hemisphere"
    WRITE(1,'(A)') "%%Creator: XR, Hannu Parviainen."
    WRITE(1,'(A)') "%%Oriantation: Portrait"
    WRITE(1,'(A)') "%%Pages: 1"
    WRITE(1,'(A, 2(I5))') "%%BoundingBox: 0 0", int(scale + 2.0 * margin), int(.25 * scale + 2.0 * margin)
    WRITE(1,'(A)') "%%EndComments\n"
    WRITE(1,'(A)') ""

    cHeight = h%dTheta / TWO_PI
    DO i = 1, h%res_theta
       cWidth = h%rows(i)%dPhi / TWO_PI
       DO j = 1, h % rows(i) % res_phi
          WRITE(1,*) h % rows(i) % cells(j,lvl), " setgray"
          WRITE(1,*) &
               & (             j-1) * cWidth  * scale + margin, &
               & (h% res_theta - i) * cHeight * scale + margin, &
               & cWidth  * scale,                               &
               & cHeight * scale,                               &
               & " rectfill"
       END DO
    END DO
 
    CLOSE(1)

  END SUBROUTINE gth_hemisphere_plot_data

  SUBROUTINE gth_hs_plot(h, fileName)
    type(gth_hemisphere) :: h
    CHARACTER(LEN=*) :: fileName
    INTEGER      :: i, j, k
    REAL         :: margin, scale, cWidth, cHeight
    REAL         :: width, height, lvlHeight

    margin      =    10.0
    scale       =  1000.0

    width       = scale + 2.0 * margin
    lvlHeight   = .25 * scale + 2.0 * margin
    height      = h%nLevels * lvlHeight


    OPEN(1, FILE=fileName, status="REPLACE")
    
    WRITE(1,'(A)') "%!PS-Adobe-2.0 EPSF-1.2"
    WRITE(1,'(A)') "%%Title: hemisphere"
    WRITE(1,'(A)') "%%Creator: XR, Hannu Parviainen."
    WRITE(1,'(A)') "%%Oriantation: Portrait"
    WRITE(1,'(A)') "%%Pages: 1"
    WRITE(1,'(A, 2(I5))') "%%BoundingBox: 0 0", int(width), int(height)
    WRITE(1,'(A)') "%%EndComments\n"
    WRITE(1,'(A)') ""

    do k = 1, h%nLevels
       cHeight = h%dTheta * INV_TWO_PI
       DO i = 1, h%res_theta
          cWidth = h%rows(i)%dPhi * INV_TWO_PI
          DO j = 1, h % rows(i) % res_phi
             WRITE(1,*) h % rows(i) % cells(j,k), " setgray"
             WRITE(1,*) &
                  &     (             j-1) * cWidth  * scale + margin, &
                  & ((h%nLevels - k) * lvlHeight) + (h% res_theta - i) * cHeight * scale + margin, &
                  & cWidth  * scale,                                   &
                  & cHeight * scale,                                   &
                  & " rectfill"
          END DO
       END DO
    end do

    CLOSE(1)

  END SUBROUTINE gth_hs_plot

  SUBROUTINE gth_hemisphere_plotSpherical(h, fileName, level)
    type(gth_hemisphere) :: h
    CHARACTER(LEN=*) :: fileName
    integer, optional :: level
    INTEGER      :: i, j, lvl
    REAL         :: margin, scale, cWidth, cHeight
 
    real         :: cX, cY

    lvl = 1
    if(present(level)) lvl = level

    margin  =    10.0
    scale   =   250.0

    cX = 0.5 * (scale*4.0 + margin)
    cY = 0.5 * (scale*4.0 + margin)

    OPEN(1, FILE=fileName, status="REPLACE")
    
    WRITE(1,'(A)') "%!PS-Adobe-2.0 EPSF-1.2"
    WRITE(1,'(A)') "%%Title: Hemisphere"
    WRITE(1,'(A)') "%%Creator: XR, Hannu Parviainen."
    WRITE(1,'(A)') "%%Oriantation: Portrait"
    WRITE(1,'(A)') "%%Pages: 1"
    WRITE(1,'(A, 2(I5))') "%%BoundingBox: 0 0", int(scale*4.0 + 2.0 * margin), int(scale*4.0 + 2.0 * margin)
    WRITE(1,'(A)') "%%EndComments\n"
    WRITE(1,'(A)') ""

    cHeight = h%dTheta / TWO_PI
    DO i = 1, h%res_theta
       cWidth = h%rows(i)%dPhi / TWO_PI
       DO j = 1, h % rows(i) % res_phi
          write(1,'("newpath")')
          write(1,*) cX + scale * (i-1) * h%dTheta * cos((j-1) * h%rows(i)%dPhi), cY + (i-1) * &
               & scale * h%dTheta * sin((j-1) * h%rows(i)%dPhi), "moveto"
          write(1,*) cX + scale * (i-1) * h%dTheta * cos((j  ) * h%rows(i)%dPhi), cY + (i-1) * &
               & scale * h%dTheta * sin((j  ) * h%rows(i)%dPhi), "lineto"
          write(1,*) cX + scale *  i    * h%dTheta * cos((j  ) * h%rows(i)%dPhi), cY +  i    * &
               & scale * h%dTheta * sin((j  ) * h%rows(i)%dPhi), "lineto"
          write(1,*) cX + scale *  i    * h%dTheta * cos((j-1) * h%rows(i)%dPhi), cY +  i    * &
               & scale * h%dTheta * sin((j-1) * h%rows(i)%dPhi), "lineto"
          write(1,*) "closepath"

          write(1,*)  h % rows(i) % cells(j, lvl), " setgray"

          !write(1,'("stroke")')
          write(1,'("fill")')

       END DO
    END DO
 
    write(1,'("showpage")')

    CLOSE(1)

  END SUBROUTINE gth_hemisphere_plotSpherical

  subroutine gth_nfCheck(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine gth_nfCheck

END MODULE hemisphere
