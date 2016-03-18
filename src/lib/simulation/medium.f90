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

!> 
!!
!! \author  Hannu Parviainen
!! \date    26.11.2007
!!
module medium
!$ use omp_lib
   use netcdf
   use base
   use particle
   use random
   use container_grid3D

  implicit none

  private :: rgen

  !****s* Medium/med_medium
  ! NAME
  !   med_medium
  !
  ! DESCRIPTION
  !
  !******
  type :: med_medium
     real(FD) :: width, hwidth, height
     integer  :: nParts

     real(fd) :: hMin, hMax, hMean, hStd
     real(fd) :: rho

     character(len=256) :: name
     character(len=64)  :: pDistribution
     character(len=128) :: pDistParms

     integer  :: nDrops

     integer  :: dropSeed, sDistSeed

     !****d* Medium/parts
     ! NAME
     !   parts
     !
     ! DESCRIPTION
     !   Particle array.
     !
     !******
     type(prt_sphere), dimension(:),  pointer :: parts => null()

     !****d* Medium/grid
     ! NAME
     !   grid
     !
     ! DESCRIPTION
     !
     !******
     type(cnt_grid3D), pointer :: grid => null()
  end type med_medium

  type :: med_mediumFile
     integer :: fileID
     integer :: selectionIdx
     integer :: nMedia
     integer :: nSelectedMedia = 0
     integer, dimension(:), pointer :: varSelection => NULL()
  end type med_mediumFile

  type(rnd_gen), save :: rgen

contains
  
  !****f* Medium/med_medium_init
  ! NAME
  !   med_medium_initt
  !
  ! DESCRIPTION
  !   Initializes the medium.
  !
  ! INPUTS
  !   c : cnt_grid3D ... Container
  !   w : real       ... Width of the container
  !   r : real       ... Resolution of the grid
  !   n : integer    ... Number of particles
  !
  ! SOURCE
  subroutine med_medium_init(m, w, h, n, name, verbose)
    type(med_medium),  intent(INOUT) :: m
    integer,           intent(IN)    :: n
    real(FD),          intent(IN)    :: w, h
    character(len=*),  intent(in)    :: name
    logical, optional, intent(in)    :: verbose

    real(FD), dimension(n) :: temp_r
    integer :: i, j, k

    logical v

    v = .false.
    if(present(verbose)) v = verbose
    
    if(v) write(*,'("Initializing medium.")')

    m % name = name

    m % width   = w
    m % height  = h
    m % nParts  = n

    m % hWidth  = w * 0.5
 
    if(v) write(*,'(" -- Allocating space for ",(I7)," particles")') m%nParts

    allocate(m%parts(n))

    if(v) write(*,'("Medium initialization complete")')

  end subroutine med_medium_init
  !******


  subroutine med_medium_destroy(m)
    type(med_medium), intent(INOUT) :: m

    if(associated(m%grid)) call med_gridRemove(m)

    if(associated(m%parts)) then
       deallocate(m%parts)
       nullify(m%parts)
    end if
      
  end subroutine med_medium_destroy

  subroutine med_setName(m, name)
    type(med_medium), intent(INOUT) :: m
    character(len=*)   :: name

    write(m%name, '(A)') trim(adjustl(name))

  end subroutine med_setName

  subroutine med_pack(m, nDrops, seed, dropTable, timer)
    type(med_medium), intent(INOUT) :: m
    integer :: nDrops
    integer, optional :: seed
    integer, dimension(:), optional :: dropTable
    type(utl_timer), optional       :: timer


    if(present(timer)) then
       call cnt_grid3D_packParticles(M%grid, nDrops, dropTable, timer)
    else
       call cnt_grid3D_packParticles(M%grid, nDrops, dropTable)
    end if

  end subroutine med_pack

  logical function med_is_point_inside_medium(m, p) result(isInside)
    type(med_medium), intent(IN) :: m
    real(fd), intent(IN)         :: p(3)

    if((abs(p(1)) > m%hWidth) .or. (abs(p(2)) > m%hWidth)) then
       isInside = .false.
       return
    else if(abs(p(3)) > m%height) then
       isInside = .false.
       return
    else 
       isInside = .true.
    end if
       
  end function med_is_point_inside_medium


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! GRID
  !!
  !!

  subroutine med_gridAssign(m, r_xy, r_z, verbose)
    type(med_medium), intent(INOUT) :: m
    integer :: r_xy, r_z
    logical, optional :: verbose

    logical :: v

    v = .false.

    if(present(verbose)) v = verbose

    if(associated(m%grid)) call med_gridRemove(m)

    allocate(m%grid)
    call cnt_grid3D_init(m%grid, m%width, m%height, r_xy, r_z, m%parts, v)

  end subroutine med_gridAssign

  subroutine med_gridRemove(m)
    type(med_medium), intent(INOUT) :: m

    call cnt_grid3D_destroy(m%grid)
    deallocate(m%grid)
    nullify(m%grid)

  end subroutine med_gridRemove


  subroutine med_gridFit(m)
    type(med_medium), intent(INOUT) :: m

    real(fd) :: maxX, minX, maxY, minY, maxZ, minZ
    integer  :: x, y, z, nv, nh
    real(fd) :: dx, dy, dz
    real(fd) :: xlim(2), ylim(2), zLim(2)
    integer(il) :: k, i

    type(prt_sphere), dimension(:),  pointer, save :: parts

!!$    nh = 10
!!$    nv = 3
!!$
!!$    allocate(parts(m%nParts))
!!$
!!$    maxX = MAXVAL(m % parts(:) % P(1) + m % parts(:) % r)
!!$    minX = MINVAL(m % parts(:) % P(1) - m % parts(:) % r)
!!$    maxY = MAXVAL(m % parts(:) % P(2) + m % parts(:) % r)
!!$    minY = MINVAL(m % parts(:) % P(2) - m % parts(:) % r)
!!$    maxZ = MAXVAL(m % parts(:) % P(3) + m % parts(:) % r)
!!$    minZ = MINVAL(m % parts(:) % P(3) - m % parts(:) % r)
!!$
!!$    dx = (maxX - minX) / real(nh,fd)
!!$    dy = (maxY - minY) / real(nh,fd)
!!$    dz = (maxZ - minZ) / real(nv,fd)
!!$
!!$    print *, "Co Start"
!!$
!!$    !print *, m%parts(:)%P(1) < 2.5_fd
!!$
!!$    k = 1
!!$    do z = 1, nv
!!$       zLim = [minZ + dz*(z-1), minZ + dz*z]
!!$       do y = 1, nh
!!$          yLim = [minY + dy*(y-1), minY + dy*y]
!!$          do x = 1, nh
!!$             xlim = [minX + dx*(x-1), minX + dx*x]
!!$             do i = 1, m%nParts
!!$                if( &
!!$                     & m%parts(i)%P(1) > xLim(1) .and. m%parts(i)%P(1) < xLim(2) .and. &
!!$                     & m%parts(i)%P(2) > yLim(1) .and. m%parts(i)%P(2) < yLim(2) .and. &
!!$                     & m%parts(i)%P(3) > zLim(1) .and. m%parts(i)%P(3) < zLim(2) ) then
!!$                   
!!$                   parts(k) = m%parts(i)
!!$                   k = k + 1
!!$                end if
!!$             end do
!!$          end do
!!$       end do
!!$    end do
!!$
!!$    print *,"CO Done\n"
!!$
!!$    deallocate(m%parts)
!!$    m%parts => parts
!!$    m%grid%parts => parts

    call cnt_grid3D_optimizeGrid(m%grid)

    m % height = m%grid%height

  end subroutine med_gridFit


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! STATISTICS
  !!
  !!

  subroutine med_updateStatistics(m)
    type(med_medium), intent(INOUT) :: m

    call cnt_grid3D_cmpStats(m%grid, 300, m%hMean, m%hStd, m%hMin, m%hMax)
 
    m % rho    = 1.0 - cnt_grid3D_cmpDensity(m%grid)

  end subroutine med_updateStatistics

  subroutine med_computePorosityStructure(m, nSamples, rho)
    type(med_medium), intent(IN) :: m
    integer                      :: nSamples
    real(fd), dimension(:,:)     :: rho

    real(fd)                     :: dz, hdz
    integer                      :: i

    dz  = m%grid%height / real(size(rho(:,1)), fd)
    hdz = 0.5_fd * dz

    do i=1, size(rho(:,1))
       rho(i,1) = i * dz - hdz
    end do

    call cnt_grid3D_cmpDensityStructure(m%grid, nSamples, rho(:,2))

    rho(:,2) = 1.0_fd - rho(:,2)

  end subroutine med_computePorosityStructure

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! DISTRIBUTIONS
  !!
  !!

  subroutine med_setDistribution_constant(m, v)
    type(med_medium) :: m
    real(FD) :: v

    m % parts(:) % r  = v
    m % parts(:) % rr = m % parts % r ** 2

    write(m % pDistribution, '("constant")')
    write(m % pDistParms,    '(E15.7)') v

  end subroutine med_setDistribution_constant

  subroutine med_setDistribution_uniform(m, vMin, vMax, seed)
    type(med_medium) :: m
    real(FD) :: vMin, vMax, rTemp(m%nParts)
    integer, optional :: seed

    call rnd_generate_uniform(vMin, vMax, rTemp)

    m % parts(:) % r  = rTemp
    m % parts(:) % rr = m % parts % r ** 2

    write(m % pDistribution, '("uniform")')
    write(m % pDistParms,    '(2(E15.7))') vMin, vMax

  end subroutine med_setDistribution_uniform

  subroutine med_setDistribution_invGamma(m, a, b, vMin, vMax)
    type(med_medium) :: m
    real(FD) :: a, b, vMin, vMax
    real(FD) :: rTemp(m%nParts)

    call utl_fatal_error("Inverse gamma distributin needs fixing, duh.")

!!$    call rnd_generate_invGamma(rGen, a, b, vMin, vMax, rTemp)
!!$
!!$    m % parts(:) % r  = rTemp
!!$    m % parts(:) % rr = m % parts % r ** 2
!!$for
!!$    write(m % pDistribution, '("invGamma")')

  end subroutine med_setDistribution_invGamma

  subroutine med_setDistribution_logNormal(m, mean, std, seed)
    type(med_medium) :: m
    real(FD) :: mean, std
    real(FD) :: rTemp(m%nParts)
    integer, optional :: seed

    call rnd_generate_logNormal(mean, std, rTemp)

    m % parts(:) % r  = rTemp
    m % parts(:) % rr = m % parts % r ** 2

    write(m % pDistribution, '("logNormal")')
    write(m % pDistParms,    '(2(E15.7))') mean, std

  end subroutine med_setDistribution_logNormal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! FILE IO
  !!
  !!

  subroutine med_medium_save(m, fName)
    type(med_medium), INTENT(IN) :: m
    character(len=*), intent(in) :: fName

    real(fd), dimension(m%nParts, 4) :: dataOut

    real(fd), dimension(300,2)       :: zRho

    integer  :: fileID, dataID, rhoID, dimID(2), rhoDimID(2)
    integer  :: i, j, k, status

    character(len=300) :: rhoName, dataName

    dataOut(:,1) = m%parts(:)%p(1)
    dataOut(:,2) = m%parts(:)%p(2)
    dataOut(:,3) = m%parts(:)%p(3)
    dataOut(:,4) = m%parts(:)%r

    call med_computePorosityStructure(m, 5000, zRho)

!    write(rhoName,'((A),("_rho"))') trim(adjustl(m%name))
!    write(dataName,'((A),("_particleData"))') trim(adjustl(m%name))
    write(rhoName,'("rho")') 
    write(dataName,'("particleData")')


    status = nf90_create (fName, NF90_NOCLOBBER,  fileID)
 
    if(status /= NF90_NOERR) THEN
       if(status == NF90_EEXIST) THEN
          call nfCheck ( nf90_open (fName, NF90_WRITE, fileID) )
          call nfCheck ( nf90_inq_dimid(fileID, "pIdx", dimID(1)) )
          call nfCheck ( nf90_inq_dimid(fileID, "xyzr", dimID(2)) )
          call nfCheck ( nf90_inq_dimid(fileID, "hrho", rhoDimID(1)) )
          call nfCheck ( nf90_inq_dimid(fileID,   "hr", rhoDimID(2)) )
          call nfCheck ( nf90_redef(fileID) )
       else 
          call nfCheck(status)
       end if
    else
       call nfCheck( nf90_def_dim(fileID, "pIdx",     m%nParts, dimID(1)) )
       call nfCheck( nf90_def_dim(fileID, "xyzr",            4, dimID(2)) )
       call nfCheck( nf90_def_dim(fileID, "hrho", size(zRho,1), rhoDimID(1)) )
       call nfCheck( nf90_def_dim(fileID,   "hr", size(zRho,2), rhoDimID(2)) )
    end if

    status = nf90_inq_varid(fileID, dataName, dataID)

    if(status == NF90_NOERR) then
       write(*,'("Warning: given medium name exists, skipping.")')
       call nfCheck( nf90_close  (fileID) )
       return
    end if

!    call nfCheck( nf90_def_var(fileID, dataName, NF90_DOUBLE, dimID, dataID) )
!    call nfCheck( nf90_def_var(fileID,  rhoName, NF90_DOUBLE, dimID, rhoID) )

    call nfCheck( nf90_def_var(fileID, dataName, NF90_FLOAT,    dimID, dataID) )
    call nfCheck( nf90_def_var(fileID,  rhoName, NF90_FLOAT, rhoDimID, rhoID) )

    call nfCheck( nf90_put_att(fileID, dataID, "width",     (/m%width/)  ) )
    call nfCheck( nf90_put_att(fileID, dataID, "height",    (/m%height/) ) )

    call nfCheck( nf90_put_att(fileID, dataID, "distribution", m%pDistribution ) )
    call nfCheck( nf90_put_att(fileID, dataID, "distParams",   adjustl(m%pDistParms) ) )
    call nfCheck( nf90_put_att(fileID, dataID, "distSeed",     (/m%sDistSeed/) ) )
    call nfCheck( nf90_put_att(fileID, dataID, "nDrops",       (/m%nDrops/) ) )
    call nfCheck( nf90_put_att(fileID, dataID, "dropSeed",     (/m%dropSeed/) ) )
    call nfCheck( nf90_put_att(fileID, dataID, "porosityMean", (/m%rho/)    ) )
    call nfCheck( nf90_put_att(fileID, dataID, "hMean",        (/m%hMean/)  ) )
    call nfCheck( nf90_put_att(fileID, dataID, "hMin",         (/m%hMin/)   ) )
    call nfCheck( nf90_put_att(fileID, dataID, "hMax",         (/m%hMax/)   ) )
    call nfCheck( nf90_put_att(fileID, dataID, "hStd",         (/m%hStd/)   ) )
    
    call nfCheck( nf90_enddef (fileID) )
    call nfCheck( nf90_put_var(fileID, dataID, dataOut) )
    call nfCheck( nf90_put_var(fileID, rhoID, zRho ) )
    call nfCheck( nf90_close  (fileID) )

  end subroutine med_medium_save


  subroutine med_mediumFileOpen(mf, fName)
    type(med_mediumFile), intent(INOUT) :: mf
    character(len=*),     intent(IN)    :: fName

    call nfCheck( nf90_open(fName, NF90_NOWRITE, mf%fileID) )
    call nfCheck( nf90_inquire(mf%fileID, nVariables = mf%nMedia) )

  end subroutine  med_mediumFileOpen


  subroutine med_mediumFileClose(mf)
    type(med_mediumFile), intent(INOUT) :: mf

    call nfCheck( nf90_close(mf%fileID) )
    mf%fileID = 0

  end subroutine  med_mediumFileClose

  subroutine med_mediumFileSelectMedia(mf, rho, rhoEps, distType)
    type(med_mediumFile), intent(INOUT) :: mf
    real(fd), optional,   intent(IN)    :: rho
    real(fd), optional,   intent(IN)    :: rhoEps
    character(len=*), optional          :: distType

    character(len=NF90_MAX_NAME) :: varName
    character(len=64)            :: medSizeDistribution

    real(fd) :: medRho(1), rhoE
    integer  :: varID, nVarsTotal, nMedia, varType, ii

    integer  :: mediaTable(50)

    logical  :: mediumOk

    if(associated(mf%varSelection)) deallocate(mf%varSelection)

    mf%selectionIdx = 1

    nMedia   = 0
    mediumOk = .true.
    
    if(present(rhoEps)) then
       rhoE = rhoEps
    else
       rhoE = 0.05_fd
    end if

    call nfCheck ( nf90_inquire(mf%fileID, nVariables = nVarsTotal ) )

    do varID = 1, nVarsTotal
       mediumOk = .true.

       call nfCheck ( nf90_inquire_variable(mf%fileID, varID, varName, varType) )

       if( index(varName, "particleData") < 1 ) mediumOk = .false.
  
       if(mediumOk .and. present(rho)) then
          call nfCheck ( nf90_get_att(mf%fileID, varID, "porosityMean", medRho) )            
          if(abs(medRho(1) - rho) > rhoE) mediumOk = .false.
       end if
  
       if(mediumOk .and. present(distType)) then
          call nfCheck ( nf90_get_att(mf%fileID, varID, "distribution", medSizeDistribution) ) 
          if(distType /= medSizeDistribution) mediumOk = .false.
       end if

       if(mediumOk) then
          nMedia = nMedia + 1
          mediaTable(nMedia) = varID
       end if

    end do

    if(nMedia > 0) then
       allocate(mf%varSelection(nMedia))
       mf%varSelection = mediaTable(1:nMedia)
       mf%nSelectedMedia = nMedia
    else
       write(*,*) "ERROR: no suitable media found from the file."
       stop
    end if
    
  end subroutine  med_mediumFileSelectMedia


  subroutine med_mediumFileRead(m, mf, idx)
    type(med_medium),     intent(INOUT) :: m
    type(med_mediumFile), intent(INOUT) :: mf
    integer, optional                   :: idx

    real(fd), dimension(:,:), allocatable :: dataIn

    real(fd), dimension(1) :: height, width
    integer  :: fileID, dataID, dimID(2), dims(2)
    integer  :: mIdx, i, j, k, status

    character(len=300) :: dataName

    call med_medium_destroy(m)
 
    if(present(idx)) then
       mIdx   = idx
       dataID = mIdx
   else
       if(mf%nSelectedMedia == 0) then
          write(*,*) "ERROR: Trying to read the next selected medium with nothing selected."
          stop
       end if

       mIdx = mf%selectionIdx
       dataID = mf%varSelection(mIdx)

       mf%selectionIdx = mf%selectionIdx + 1

    end if

    if(mIdx > mf%nMedia) then
       write(*,*) "ERROR: Trying to read a nonexistent variable."
       stop
    end if
 
    call nfCheck ( nf90_get_att(mf%fileID, dataID, "width",  width) ) 
    call nfCheck ( nf90_get_att(mf%fileID, dataID, "hMax", height) ) 
 
    call nfCheck ( nf90_inquire_variable(mf%fileID, dataID, name=dataName, dimids=dimID) )
    call nfCheck ( nf90_inquire_dimension(mf%fileID, dimID(1), len = dims(1)) )
    call nfCheck ( nf90_inquire_dimension(mf%fileID, dimID(2), len = dims(2)) )

    allocate(dataIn(dims(1),dims(2)))

    call nfCheck ( nf90_get_var(mf%fileID, dataID, dataIn) )

    call med_medium_init(m, width(1), height(1), dims(1), dataName)

    m%parts(:)%p(1) = dataIn(:,1)
    m%parts(:)%p(2) = dataIn(:,2)
    m%parts(:)%p(3) = dataIn(:,3)
    m%parts(:)%r    = dataIn(:,4)

    m%parts(:)%rr   = m%parts(:)%r**2

    deallocate(dataIn)

  end subroutine med_mediumFileRead

  !****f* medium/med_maskHeight
  ! NAME
  !   med_maskHeight
  !
  ! DESCRIPTION
  !   Masks particles with a height map.
  !
  ! INPUTS
  !   m : med_medium ... Container
  !
  !******
  subroutine med_maskHeight(m, hmap)
    type(med_medium)     :: m
    real, dimension(:,:) :: hmap

    real(fd)             :: rw, mean, std
    integer              :: p(2)
    integer              :: i, j
    integer              :: nNew


    type(prt_sphere), dimension(:), pointer :: pTemp
    logical, dimension(m%nParts)            :: pMask
    
    write(*,*) "Medium width", m%width
    write(*,*) "Medium max height", m%height
    hmap = hmap - maxval(hmap)
    mean = sum(hmap) / size(hmap)
    std = sqrt(sum((hmap - mean)**2) / size(hmap))
    
    pMask = .true.
    rw = real(size(hmap,1))/m%width
    hmap = hmap + m%height
    write(*,*) "Height map max   ", maxval(hmap)
    write(*,*) "Height map min   ", minval(hmap)
    write(*,*) "Height map ampl. ", maxval(hmap)-minval(hmap)
    
    write(*,*) "Height map mean  ", mean
    write(*,*) "Height map std   ", std
    write(*,*) "Height map size  ", size(hmap)
    

    do i=1,m%nParts
       p = floor((m % parts(i) % P(1:2) + m % hwidth) * rw) + 1
       if( m % parts(i) % P(3) + m % parts(i) % r > hmap(p(1),p(2)))  pMask(i) = .false.
    end do
    
    nNew = count(pMask)

    allocate(pTemp(nNew))

    j = 1
    do i=1,m%nParts
       if(pMask(i)) then
          pTemp(j) = m%parts(i)
          j = j + 1
       end if
    end do

    deallocate(m%parts)

    m%nParts      =  nNew
    m%grid%nParts =  nNew
    m%parts       => pTemp
    m%grid%parts  => pTemp

    call med_gridFit(m)

  end subroutine med_maskHeight


  !****f* medium/med_maskHeight
  ! NAME
  !   med_maskHeight
  !
  ! DESCRIPTION
  !   Masks particles with a height map.
  !
  ! INPUTS
  !   m : med_medium ... Container
  !
  !******
  subroutine med_clipHeight(m, v, type)
    type(med_medium)           :: m
    real(fd)                   :: v
    character(len=*), optional :: type

    real(fd)             :: clipHeight

    integer              :: i, j
    integer              :: nNew

    type(prt_sphere), dimension(:), pointer :: pTemp
    logical, dimension(m%nParts)            :: pMask
    
    pMask = .true.

    if(present(type)) then
       if (type == 'absolute') then
          clipHeight = v
       else
          clipHeight = v * m%Height
       end if
    else
       clipHeight = v * m%Height
    end if

    do i=1,m%nParts
       if( m % parts(i) % P(3) + m % parts(i) % r > clipHeight)  pMask(i) = .false.
    end do
    
    nNew = count(pMask)

    allocate(pTemp(nNew))

    j = 1
    do i=1,m%nParts
       if(pMask(i)) then
          pTemp(j) = m%parts(i)
          j = j + 1
       end if
    end do

    deallocate(m%parts)

    m%nParts      =  nNew
    m%grid%nParts =  nNew
    m%parts       => pTemp
    m%grid%parts  => pTemp

    call med_gridFit(m)

  end subroutine med_clipHeight


  !****f* medium/med_maskHeight
  ! NAME
  !   med_maskHeight
  !
  ! DESCRIPTION
  !   Masks particles with a height map.
  !
  ! INPUTS
  !   m : med_medium ... Container
  !
  !******
  subroutine med_maskHeightToNew(m1, m2, hmap)
    type(med_medium)     :: m1, m2
    real, dimension(:,:) :: hmap

    real(fd)             :: rw
    integer              :: p(2)
    integer              :: i, j
    integer              :: nNew

    type(prt_sphere), dimension(:), pointer :: pTemp
    type(cnt_grid3D),               pointer :: gTemp
    logical, dimension(m1%nParts)           :: pMask

    if(associated(m2%parts)) then
       deallocate(m2%parts)
       nullify(m2%parts)
    end if

    gTemp => m2%grid

    m2 = m1

    m2 % grid => gTemp

    pMask = .true.
    rw = real(size(hmap,1))/m1%width
    hmap = hmap + m1%height - maxval(hmap)

    !print*, "minHeight", minval(hmap)

    do i=1,m1%nParts
       p = floor((m1 % parts(i) % P(1:2) + m1 % hwidth) * rw) + 1
       if( m1 % parts(i) % P(3) + m1 % parts(i) % r > hmap(p(1),p(2)))  pMask(i) = .false.
    end do
    
    nNew = count(pMask)

    allocate(pTemp(nNew))

    j = 1
    do i=1,m1%nParts
       if(pMask(i)) then
          pTemp(j) = m1%parts(i)
          j = j + 1
       end if
    end do

    m2%nParts      =  nNew
    m2%grid%nParts =  nNew
    m2%parts       => pTemp
    m2%grid%parts  => pTemp

    call med_gridFit(m2)

  end subroutine med_maskHeightToNew



  !****f* medium/med_exportRib
  ! NAME
  !   med_exportRib
  !
  ! DESCRIPTION
  !   
  !
  ! INPUTS
  !   m : med_medium ... Medium
  !
  !******
  subroutine med_exportRib(m, fNameIn, cut)
    type(med_medium)  :: m
    character(len=*), optional  :: fNameIn
    logical, optional :: cut
    integer           :: i
    character(len=256) :: fName

    logical, dimension(m%nParts) :: export

    export = .true.

    if(present(fNameIn)) then
       open(1, file=fNameIn, status='REPLACE')
    else
       write(fName,'((A),(".rib"))'), trim(m%name)
       open(1, file=fName, status='REPLACE')
    end if

    if(present(cut)) then
       do i=1, m % nParts
          if(sum(m%parts(i)%P(1:2)**2) > m%hWidth**2) export(i) = .false. 
       end do
    end if

    write(1,*) 'Points "P" [' 

    do i=1, m % nParts
       if(export(i)) write(1,*) m%parts(i) % P
    end do

    write(1,*) '] "width" ['

    do i=1, m % nParts
       if(export(i)) write(1,*) 2.0 * m%parts(i)%r
    end do

    write(1,*) '] "uniform string type" "sphere"'

    close(1)

  end subroutine med_exportRib

  subroutine nfCheck(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine nfCheck

end module medium
