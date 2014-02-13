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
Module hemisphere
  Use base
  Use netcdf
  Use spectrum

  Implicit None

  Type :: gth_hemisphere
     Character(36) :: name

     Integer (IS)  :: resTheta
     Real    (FD)  :: dTheta, dThetaInv

     Real    (FD)  :: dAMean
     Integer (IL)  :: nCells, nThetaI, nLevels

     Real    (FD), Dimension(:),     Pointer :: dPhi, dPhiInv, dA, mTheta
     Integer (IL), Dimension(:),     Pointer :: resPhi, cIdx

     Real    (FD), Dimension(:,:,:), Pointer :: data
     Real    (FD), Dimension(:,:,:), Pointer :: weight

  End Type gth_hemisphere

  Type :: gth_hsFile
     Integer :: fileID  
     Integer :: dataID
     Logical :: isDef
  End Type gth_hsFile

Contains

  Subroutine gth_hemisphere_init(h, resTheta, nThetaI_, nDataLevels_)
    Type(gth_hemisphere) :: h
    Integer              :: resTheta
    Integer, Optional    :: nThetaI_, nDataLevels_

    integer              :: nThetaI = 1, nDatalevels = 1

    Real(FD)             :: dPhi, dA0
    Integer              :: rPhi, i, j

    if( present ( nThetaI_     ) ) nThetaI = nThetaI_ 
    if( present ( nDataLevels_ ) ) nDataLevels = nDataLevels_

    h % resTheta  = resTheta
    h % dTheta    = HALF_PI / Real(resTheta, fd)
    h % dThetaInv = 1.0_fd / h % dTheta

    h % nCells    = 1
    h % nThetaI   = nThetaI
    h % nLevels   = nDataLevels

    Allocate( h % dPhi    ( resTheta ) )
    Allocate( h % dPhiInv ( resTheta ) )
    Allocate( h % dA      ( resTheta ) )
    Allocate( h % mTheta  ( resTheta ) )
    Allocate( h % resPhi  ( resTheta ) )
    Allocate( h % cIdx    ( resTheta ) )

    dA0 = TWO_PI * (1.0_fd - Cos(h%dTheta))

    h % resPhi  (1) = 1
    h % cIdx    (1) = 1
    h % dPhi    (1) = TWO_PI
    h % dPhiInv (1) = INV_TWO_PI
    h % dA      (1) = dA0
    h % mTheta  (1) = 0.5_fd * h % dTheta


    Do i = 2, resTheta
       dPhi = dA0 / (Cos((i-1) * h%dTheta) - Cos(i * h%dTheta))
       rPhi = Nint(TWO_PI / dPhi)
       dPhi = TWO_PI / Real(rPhi)

       h % resPhi  (i) = rPhi
       h % dPhi    (i) = dPhi
       h % dPhiInv (i) = 1.0_fd / dPhi
       h % dA      (i) = dPhi * (Cos((i-1) * h%dTheta) - Cos(i * h%dTheta))
       h % mTheta  (i) = h % dTheta * (Real(i, fd) - 0.5_fd)

       h % cIdx    (i) = h % cIdx(i-1) + h%resPhi(i-1)

       h % nCells = h % nCells + rPhi   
    end do

    h % dAMean = TWO_PI / Real(h%nCells, fd)

    allocate( h % data   ( nThetaI, h%nCells, nDataLevels ) )
    allocate( h % weight ( nThetaI, h%nCells, nDataLevels ) )

    h % data   = 0.0_fd
    h % weight = 0.0_fd

  End Subroutine gth_hemisphere_init


  !****f* hemisphere/gth_sphDirToCell
  ! NAME
  !   gth_sphDirToCell
  !
  ! DESCRIPTION
  !   Returns the cell index corresponding to a direction
  !   vector given in spherical coordinates.
  !
  ! INPUTS
  !   h : hemisphere ... 
  !   p : real(2)    ... 
  !
  ! SOURCE

  Pure Function gth_sphDirToCell(h, D) Result(C)
    Type(gth_hemisphere),   Intent(IN) :: h
    Real(FD), Dimension(2), Intent(IN) :: D
    Integer                            :: C
    
    integer :: t, p

    t = Floor( D(1) * h % dThetaInv  ) + 1
    p = Floor( D(2) * h % dPhiInv(t) ) + 1

    C = h % cIdx(t) + p

  End Function gth_sphDirToCell
  !******


  !****f*  hemisphere/gth_carDirToCell
  ! NAME
  !   gth_carDirToCell
  !
  ! DESCRIPTION
  !   Returns the cell index corresponding to a direction
  !   vector given in cartesian coordinates.
  !
  ! INPUTS
  !   h : hemisphere ... 
  !   p : real(3)    ... 
  !
  ! SOURCE
  Pure Function gth_carDirToCell(h, D) Result(C)
    Type(gth_hemisphere),   Intent(IN) :: h
    Real(FD), Dimension(3), Intent(IN) :: D
    Integer                            :: C

    Real(FD) :: theta, phi, r
    integer :: t, p

    r     = Sqrt  ( Sum(D**2)          )
    theta = Acos  ( D(3) / r           )
    phi   = Atan2 ( D(2) / r, D(1) / r )
    
    If( phi < 0.0 ) phi = TWO_PI + phi
    
    t = Floor( theta * h % dThetaInv     ) + 1
    p = Floor( phi   * h % dPhiInv(t) ) + 1

    C = h % cIdx(t) + p

  End Function gth_carDirToCell
  !******

  !****f*  hemisphere/gth_cellCenterCar
  ! NAME
  !   gth_cellCenterCar
  !
  ! DESCRIPTION
  !   Returns a cartesian direction vector towards the
  !   center of the cell(theta, phi).
  !
  ! INPUTS
  !   h : hemisphere ... 
  !   p : real(3)    ... 
  !
  ! SOURCE
  Pure Function gth_cellCenterCar(h, t, p) Result(d)
    Type(gth_hemisphere), Intent(IN) :: h
    Integer,              Intent(IN) :: t, p

    Real(fd)                         :: d(3)
    Real(fd)                         :: theta, phi

    theta = (Real(t,fd) - 0.5_fd) * h % dTheta
    phi   = (Real(p,fd) - 0.5_fd) * h % dPhi(t)
    
    d(1)  = Sin(theta) * Cos(phi)
    d(2)  = Sin(theta) * Sin(phi)
    d(3)  = Cos(theta)

  End Function gth_cellCenterCar
  !******

  !****f*  hemisphere/gth_cellCenterCar
  ! NAME
  !   gth_cellCenterCar
  !
  ! DESCRIPTION
  !   Returns a random cartesian direction
  !   sample for a given cell(theta, phi).
  !
  ! INPUTS
  !   h : hemisphere ... 
  !   p : real(3)    ... 
  !
  ! SOURCE
  Function gth_cellRandomSampleCar(h, t, p) Result(d)
    Type(gth_hemisphere), Intent(IN) :: h
    Integer,              Intent(IN) :: t, p

    Real(FD)                         :: d(3)
    Real(FD)                         :: theta, phi
    
    real(FD)                         :: r(2)

    Call Random_number(r)

    theta = (Real(t,fd) - r(1)) * h % dTheta
    phi   = (Real(p,fd) - r(2)) * h % dPhi(t)
    
    d(1) = Sin(theta) * Cos(phi);
    d(2) = Sin(theta) * Sin(phi);
    d(3) = Cos(theta);

  End Function gth_cellRandomSampleCar
  !******


  Subroutine gth_hs_addToRow(h,t,v)
    Type(gth_hemisphere),   Intent(INOUT)  :: h
    Integer                                :: t
    Real(FD)                               :: v

    h % data(1, h%cIdx(t) : h%cIdx(t) + h%resPhi(t), 1) = h % data(1, h%cIdx(t) : h%cIdx(t) + h%resPhi(t), 1) + v
 
  End Subroutine gth_hs_addToRow


  Subroutine gth_addDataCar(h, d, v, tIn, l)
    Type(gth_hemisphere),   Intent(INOUT)  :: h
    Real(FD), Dimension(3), Intent(IN)     :: d
    Real(FD),               Intent(IN)     :: v
    Integer, Optional,      Intent(IN)     :: tIn, l

    Integer                :: c
    Integer                :: thtIn = 1, lvl = 1

    If(Present(l))   lvl   = l
    If(Present(tIn)) thtIn = tIn

    c = gth_carDirToCell(h, d)

    h % data(thtIn, c, lvl) = h % data(thtIn, c, lvl) + v

  End Subroutine gth_addDataCar


  Subroutine gth_hsFileOpen(hf, fName, mode)
    Type(gth_hsFile), Intent(INOUT) :: hf
    Character(len=*), Intent(IN)    :: fName
    Character(len=1), Intent(IN)    :: mode

    If(mode == 'r') Then
       Call gth_nfCheck( nf90_open(fName, NF90_NOWRITE, hf%fileID) )
    Else
       Call gth_nfCheck( nf90_create(fName, NF90_CLOBBER, hf%fileID) )
    End If
    
    hf%dataID    = 0
    hf%isDef     = .True.
    
  End Subroutine gth_hsFileOpen

  Subroutine gth_hsFileClose(hf)
    Type(gth_hsFile), Intent(INOUT) :: hf

    Call gth_nfCheck( nf90_close(hf%fileID) )
    hf%fileID = 0
  End Subroutine gth_hsFileClose

  Subroutine gth_hsWriteHeader(hf, h, auth, prog, vers)
    Type(gth_hsFile),     Intent(INOUT) :: hf
    Type(gth_hemisphere), Intent(IN)    :: h

    Character(len=*) :: auth, prog, vers

    Integer  :: dimID(3)
  
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "thetaI", h%nThetaI, dimID(1)) )
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "data",   h%nCells,  dimID(2)) )
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "level",  h%nLevels, dimID(3)) )

    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "Author",  auth ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "Program", prog ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "Version", vers ) )

    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "nThetaI", [ h%nThetaI  ]  ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "nThetaE", [ h%resTheta ]  ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "nCells",  [ h%nCells   ]  ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "nLevels", [ h%nLevels  ]  ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "dTheta",  [ h%dTheta   ]  ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "dAMean",  [ h%dAMean   ]  ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "nPhi",    h%resPhi        ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "cIdx",    h%cIdx          ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "dPhi",    h%dPhi          ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "dA",      h%dA            ) )
    Call gth_nfCheck( nf90_put_att(hf%fileID, NF90_GLOBAL, "mTheta",  h%mTheta        ) )

    Call gth_nfCheck( nf90_def_var(hf%fileID, "Hemisphere", NF90_DOUBLE, dimID, hf%dataID) )

    Call gth_nfCheck( nf90_put_att(hf%fileID,  hf%dataID, "long_name",  "Hemisphere data" ) )

  End Subroutine gth_hsWriteHeader

  Subroutine gth_hsAddAttS(hf, varIdx, attName, attData)
    Type(gth_hsFile),     Intent(INOUT) :: hf
    Integer                             :: varIdx
    Character(len=*)                    :: attName
    Character(len=*)                    :: attData

    Integer                             :: vID
    
    If(.Not. hf%isDef) Call gth_nfCheck( nf90_redef(hf%fileID) )

    vID = varIdx
    If(vID /= NF90_GLOBAL) vID = hf%dataID

    Call gth_nfCheck( nf90_put_att(hf%fileID, vID, attName, attData ) )

  End Subroutine gth_hsAddAttS

  Subroutine gth_hsAddAttF(hf, varIdx, attName, attData)
    Type(gth_hsFile),     Intent(INOUT) :: hf
    Integer                             :: varIdx
    Character(len=*)                    :: attName
    Real(fd), Dimension(:)              :: attData
    
    Integer                             :: vID
    
    If(.Not. hf%isDef) Call gth_nfCheck( nf90_redef(hf%fileID) )

    vID = varIdx
    If(vID /= NF90_GLOBAL) vID = hf%dataID

    Call gth_nfCheck( nf90_put_att(hf%fileID, vID, attName, attData ) )

  End Subroutine gth_hsAddAttF


  Subroutine gth_hsSaveData(h, hf, norm)
    Type(gth_hemisphere), Intent(IN)    :: h
    Type(gth_hsFile),     Intent(INOUT) :: hf
  
    Real(fd), Optional :: norm

    !If(Present(norm)) dataOut = dataOut * norm
    If(hf%isDef) Call gth_nfCheck( nf90_enddef(hf%fileID) )

    Call gth_nfCheck( nf90_put_var(hf%fileID, hf%dataID, h%data) )
    
  End Subroutine gth_hsSaveData

  Subroutine gth_hemisphere_plot_data(h, fileName, thetaI, level)
    Type(gth_hemisphere) :: h
    Character(LEN=*) :: fileName
    Integer, Optional :: thetaI, level
    Integer      :: i, j, thtI, lvl
    Real         :: margin, scale, cWidth, cHeight
 
    margin  =    10.0
    scale   =  1000.0

    thtI = 1
    If(Present(thetaI)) thtI = thetaI

    lvl = 1
    If(Present(level)) lvl = level

    Open(1, FILE=fileName, status="REPLACE")
    
    Write(1,'(A)') "%!PS-Adobe-2.0 EPSF-1.2"
    Write(1,'(A)') "%%Title: hemisphere"
    Write(1,'(A)') "%%Creator: XR, Hannu Parviainen."
    Write(1,'(A)') "%%Oriantation: Portrait"
    Write(1,'(A)') "%%Pages: 1"
    Write(1,'(A, 2(I5))') "%%BoundingBox: 0 0", Int(scale + 2.0 * margin), Int(.25 * scale + 2.0 * margin)
    Write(1,'(A)') "%%EndComments\n"
    Write(1,'(A)') ""

    cHeight = h%dTheta / TWO_PI
    Do i = 1, h%resTheta
       cWidth = h%dPhi(i) / TWO_PI
       Do j = 1, h%resPhi(i)
          Write(1,*) h % data(thtI, h%cIdx(i) + j, lvl), " setgray"
          Write(1,*) &
               & (             j-1) * cWidth  * scale + margin, &
               & (h% resTheta - i)  * cHeight * scale + margin, &
               & cWidth  * scale,                               &
               & cHeight * scale,                               &
               & " rectfill"
       End Do
    End Do
 
    Close(1)

  End Subroutine gth_hemisphere_plot_data

  Subroutine gth_hs_plot(h, fileName)
    Type(gth_hemisphere) :: h
    Character(LEN=*) :: fileName
    Integer      :: i, j, k
    Real         :: margin, scale, cWidth, cHeight
    Real         :: width, height, lvlHeight

    margin      =    10.0
    scale       =  1000.0

    width       = scale + 2.0 * margin
    lvlHeight   = .25 * scale + 2.0 * margin
    height      = h%nLevels * lvlHeight


    Open(1, FILE=fileName, status="REPLACE")
    
    Write(1,'(A)') "%!PS-Adobe-2.0 EPSF-1.2"
    Write(1,'(A)') "%%Title: hemisphere"
    Write(1,'(A)') "%%Creator: XR, Hannu Parviainen."
    Write(1,'(A)') "%%Oriantation: Portrait"
    Write(1,'(A)') "%%Pages: 1"
    Write(1,'(A, 2(I5))') "%%BoundingBox: 0 0", Int(width), Int(height)
    Write(1,'(A)') "%%EndComments\n"
    Write(1,'(A)') ""

!!$    Do k = 1, h%nLevels
!!$       cHeight = h%dTheta * INV_TWO_PI
!!$       Do i = 1, h%res_theta
!!$          cWidth = h%rows(i)%dPhi * INV_TWO_PI
!!$          Do j = 1, h % rows(i) % res_phi
!!$             Write(1,*) h % rows(i) % cells(j,k), " setgray"
!!$             Write(1,*) &
!!$                  &     (             j-1) * cWidth  * scale + margin, &
!!$                  & ((h%nLevels - k) * lvlHeight) + (h% res_theta - i) * cHeight * scale + margin, &
!!$                  & cWidth  * scale,                                   &
!!$                  & cHeight * scale,                                   &
!!$                  & " rectfill"
!!$          End Do
!!$       End Do
!!$    End Do

    Close(1)

  End Subroutine gth_hs_plot

  Subroutine gth_hemisphere_plotSpherical(h, fileName, level)
    Type(gth_hemisphere) :: h
    Character(LEN=*) :: fileName
    Integer, Optional :: level
    Integer      :: i, j, lvl
    Real         :: margin, scale, cWidth, cHeight
 
    Real         :: cX, cY

    lvl = 1
    If(Present(level)) lvl = level

    margin  =    10.0
    scale   =   250.0

    cX = 0.5 * (scale*4.0 + margin)
    cY = 0.5 * (scale*4.0 + margin)

    Open(1, FILE=fileName, status="REPLACE")
    
    Write(1,'(A)') "%!PS-Adobe-2.0 EPSF-1.2"
    Write(1,'(A)') "%%Title: Hemisphere"
    Write(1,'(A)') "%%Creator: XR, Hannu Parviainen."
    Write(1,'(A)') "%%Oriantation: Portrait"
    Write(1,'(A)') "%%Pages: 1"
    Write(1,'(A, 2(I5))') "%%BoundingBox: 0 0", Int(scale*4.0 + 2.0 * margin), Int(scale*4.0 + 2.0 * margin)
    Write(1,'(A)') "%%EndComments\n"
    Write(1,'(A)') ""

!!$    cHeight = h%dTheta / TWO_PI
!!$    Do i = 1, h%res_theta
!!$       cWidth = h%rows(i)%dPhi / TWO_PI
!!$       Do j = 1, h % rows(i) % res_phi
!!$          Write(1,'("newpath")')
!!$          Write(1,*) cX + scale * (i-1) * h%dTheta * Cos((j-1) * h%rows(i)%dPhi), cY + (i-1) * &
!!$               & scale * h%dTheta * Sin((j-1) * h%rows(i)%dPhi), "moveto"
!!$          Write(1,*) cX + scale * (i-1) * h%dTheta * Cos((j  ) * h%rows(i)%dPhi), cY + (i-1) * &
!!$               & scale * h%dTheta * Sin((j  ) * h%rows(i)%dPhi), "lineto"
!!$          Write(1,*) cX + scale *  i    * h%dTheta * Cos((j  ) * h%rows(i)%dPhi), cY +  i    * &
!!$               & scale * h%dTheta * Sin((j  ) * h%rows(i)%dPhi), "lineto"
!!$          Write(1,*) cX + scale *  i    * h%dTheta * Cos((j-1) * h%rows(i)%dPhi), cY +  i    * &
!!$               & scale * h%dTheta * Sin((j-1) * h%rows(i)%dPhi), "lineto"
!!$          Write(1,*) "closepath"
!!$
!!$          Write(1,*)  h % rows(i) % cells(j, lvl), " setgray"
!!$
!!$          !write(1,'("stroke")')
!!$          Write(1,'("fill")')
!!$
!!$       End Do
!!$    End Do
 
    Write(1,'("showpage")')

    Close(1)

  End Subroutine gth_hemisphere_plotSpherical

  Subroutine gth_nfCheck(status)
    Integer, Intent ( in) :: status

    If(status /= nf90_noerr) Then
       Print *, Trim(nf90_strerror(status))
       Stop "Stopped"
    End If
  End Subroutine gth_nfCheck

End Module hemisphere
