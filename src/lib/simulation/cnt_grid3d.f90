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

!> Three-dimensional regular grid container for spherical particles.
!!
!! \notes
!!   Care must be taken to set the height and vertical resolution of
!!   the container, since this has a notable effect on the packing
!!   speed. Container height of the expected packing height is fine,
!!   and height/resolution should be greater than the diameter of
!!   the largest particles.
!!
!! \author  Hannu Parviainen
!! \date    04.09.2007
!!
module container_grid3D
!$ use omp_lib
   use base
   use particle
   use random
   use netcdf

  implicit none

  private

  public  :: cnt_grid3D, cnt_gridCell

  public  :: cnt_grid3D_init,                cnt_grid3D_destroy, &
       &     cnt_grid3D_optimizeGrid,        cnt_grid3D_packParticles, &
       &     cnt_grid3D_getHeight_xy,        cnt_grid3D_export_rib, &
       &     cnt_grid3D_cmpDensity,          cnt_grid3D_cmpDensityMC, &
       &     cnt_grid3D_cmpDensityStructure, cnt_grid3D_cmpHeightMap, &
       &     cnt_grid3D_cmpStats,            cnt_grid3D_quantizePosition, &
       &     cnt_grid3D_isPointInsideParticle

  real(fd), parameter :: CNT_BORDER = 1e-3_fd

  !>   Minimum number of particles allocated for a vector.
  !!
  integer, parameter :: pMin = 5

  !>
  !!
  type :: cnt_gridCell
     integer(IL)  :: nParts, lSize
     integer, dimension(:), pointer :: pList
  end type cnt_gridCell

  !> Regular grid container structure.
  !!
  type :: cnt_grid3D
     real(fd) :: width, hWidth, height, rw, wr, rh, hr
     integer  :: res_xy, res_z, nParts

     !>   Maximum global gridded height of the container.  
     !!
     integer :: maxHeight_g 

     !>   Maximum gridded height of the container per (x,y)-cell.  
     !!
     integer, dimension(:,:), pointer :: maxHeight_l

     !>   Particle array.
     !!
     type(prt_sphere), dimension(:),     pointer :: parts

     !>   Three-dimensional array of cnt_gridCell structures holding the actual
     !!   grid-structure.
     !!
     type(cnt_gridCell), dimension(:,:,:), pointer :: grid
  end type cnt_grid3D

contains
  
  !> Initializes a cnt_grid3D container.
  !!
  !! \param[inout] c  Container
  !! \param[in]    w  Width of the container
  !! \param[in]    r  Resolution of the grid
  !! \param[in]    n  Number of particles
  !!
  !! SOURCE
  subroutine cnt_grid3D_init(c, w, h, r_xy, r_z, parts, verbose)
    type(cnt_grid3D), intent(INOUT) :: c
    integer,          intent(IN)    :: r_xy, r_z
    real(FD),         intent(IN)    :: w, h
    logical, optional, intent(in)   :: verbose

    type(prt_sphere), dimension(:), pointer :: parts

    integer :: i, j, k

    logical :: v

    v = .false.

    if(present(verbose)) v = verbose

    if(v) write(*,'(A)') "Initializing grid3D"

    c % width   = w
    c % height  = h
    c % res_xy  = r_xy
    c % res_z   = r_z
    c % nParts = size(parts)

    c % hwidth  = w * 0.5
    c % rw      = real(r_xy) / w
    c % wr      = w / real(r_xy)
    c % rh      = real(r_z) / h
    c % hr      = h / real(r_z)

    c%parts => parts

    allocate(c%grid(r_xy,r_xy,r_z))
    allocate(c%maxHeight_l(r_xy,r_xy))

    c % maxHeight_g    = 1
    c % maxHeight_l    = 1

    c % grid % nParts  = 0

    do k=1,r_z
       do j=1,r_xy
          do i=1,r_xy
             allocate(c % grid(i,j,k) % plist(pMin))
          end do
       end do
    end do

    c % grid % lsize = pMin

    if(v) write(*,'(A)') "... Complete"

  end subroutine cnt_grid3D_init
  !******

  subroutine cnt_grid3D_destroy(c)
    type(cnt_grid3D), intent(INOUT) :: c
    integer :: i, j, k

    do k=1,c%res_z
       do j=1,c%res_xy
          do i=1,c%res_xy
             deallocate(c % grid(i,j,k) % plist)
          end do
       end do
    end do

    deallocate(c%maxHeight_l)
    deallocate(c%grid)
    
    nullify(c%maxHeight_l)
    nullify(c%grid)

  end subroutine cnt_grid3D_destroy



  function cnt_grid3D_isPointInsideParticle(c, p, pInd) result(isInside)
    type(cnt_grid3D),       intent(IN)  :: c
    real(FD), dimension(3), intent(IN)  :: p 
    integer(IL), optional,  intent(OUT) :: pInd

    logical :: isInside
    integer :: i, pG(3)

    pG = 0

    isInside = .false.

    ! $omp critical
    !pG = cnt_grid3D_quantizePosition(c, p)

    !write(*,*) "T:", omp_get_thread_num(), p, pG, c%rh,  FLOOR((p(3) ) * c%rh) + 1

    pG(1)   = floor((p(1) + c%hwidth) * c%rw) + 1
    pG(2)   = floor((p(2) + c%hwidth) * c%rw) + 1
    pG(3)   = floor((p(3) ) * c%rh) + 1

    !write(*,*) "T:", omp_get_thread_num(), p, pG, c%rh,  FLOOR((p(3) ) * c%rh) + 1

    ! $omp end critical

    if(pG(3) < c%res_z .and. pG(3) > 1) then
 
       do i=1, c%grid(pG(1), pG(2), pG(3))%nParts
          if(sum((p-c%parts(c%grid(pG(1), pG(2), pG(3))%plist(i))%p)**2) < c%parts(c%grid(pG(1), pG(2), pG(3))%plist(i))%rr) then
             isInside = .true.
             if(present(pInd)) pInd     = i
             exit
          end if
       end do

    end if

  end function cnt_grid3D_isPointInsideParticle

  !>
  !!
  !!
  !!   \param[inout] c  Container
  !!   \param        si Sphere index
  !!   \param        la 
  !!   \param        lb 
  !!
  !! SOURCE
  subroutine cnt_grid3D_addParticleToCells(c, si, la, lb)
    type(cnt_grid3D) :: c    
    integer :: si, x, y, z, i, j, k
    integer, dimension(:), pointer :: tmp
    integer, dimension(3) :: la, lb



    do k = la(3), lb(3)
		z = modulo(k-1, c%res_xy) + 1
       do j = la(2), lb(2)
		   y = modulo(j-1, c%res_xy) + 1
          do i = la(1), lb(1)
			  x = modulo(i-1, c%res_xy) + 1
             !! Manual reallocation of the data array if the array is full.
             !!
             if(c % grid(x,y,z) % nParts == c % grid(x,y,z) % lsize) then
                tmp => c % grid(x,y,z) % plist
                allocate(c % grid(x,y,z) % plist(c % grid(x,y,z) % lsize + pmin))
                c % grid(x,y,z) % plist(:c % grid(x,y,z) % nParts) = tmp
                deallocate(tmp)
                c % grid(x,y,z) % lsize = c % grid(x,y,z) % lsize + pmin
             end if
             !! Add the sphere index to the array.
             !!

             c % grid(x,y,z) % plist(c % grid(x,y,z) % nParts + 1) = si
             c % grid(x,y,z) % nParts = c % grid(x,y,z) % nParts + 1

             if(z > c % maxHeight_l(x,y)) c % maxHeight_l(x,y) = z
             tmp => null()
          end do
       end do

       if(z > c % maxHeight_g) c % maxHeight_g = z

    end do
  end subroutine cnt_grid3D_addParticleToCells
  !******


  !****f* container_grid3D/cnt_grid3D_optimize_grid
  ! NAME
  !   cnt_grid3D_optimize_grid
  !
  ! DESCRIPTION
  !   
  !
  ! INPUTS
  !   c          : cnt_grid3D  ... Container
  !
  ! SOURCE
  subroutine cnt_grid3D_optimizeGrid(c, verbose)
    type(cnt_grid3D) :: c
    logical, optional :: verbose

    real(FD) :: minX, maxX, minY, maxY, pp(3)
    integer :: i, x, y, z, la(3), lb(3), perBin

    logical v

    v = .false.

    if(present(verbose)) v = verbose

    if(v) write(*,*) "Recalculating world bounds and optimizing grid."

    maxX = abs(maxval(c % parts(:) % P(1) + c % parts(:) % r))
    minX = abs(minval(c % parts(:) % P(1) - c % parts(:) % r))
    maxY = abs(maxval(c % parts(:) % P(2) + c % parts(:) % r))
    minY = abs(minval(c % parts(:) % P(2) - c % parts(:) % r))

    c % width  = 2.0_fd * (maxval((/maxX, minX, maxY, minY/)) + CNT_BORDER)
    c % hWidth = c % width * 0.5
    c % rw     = real(c%res_xy) / c%width
    c % wr     = c%width / real(c%res_xy)

    c % height = maxval(c % parts(:) % P(3) + c % parts(:) % r) + CNT_BORDER
    c % rh     = real(c%res_z) / c%height
    c % hr     = c%height / real(c%res_z)

    c % grid % nParts = 0
    c % grid % lsize   = pMin

    do z = 1, c % res_z
       do y = 1, c % res_xy
          do x = 1, c % res_xy
             deallocate(c % grid(x,y,z) % plist)
             allocate(c % grid(x,y,z) % plist(pMin))
             c % grid(x,y,z) % plist = 0
          end do
       end do
    end do

    do i = 1, c % nParts
       pp = c%parts(i)%P

       la(1:2) = floor((pp(1:2) - c%parts(i)%r + c%hwidth) * c%rw) + 1
       lb(1:2) = floor((pp(1:2) + c%parts(i)%r + c%hwidth) * c%rw) + 1

!       where(la < 1)        la = 1
!       where(lb > c%res_xy) lb = c%res_xy

       la(3) = floor((pp(3) - c%parts(i) % r) * c % rh) + 1
       lb(3) = floor((pp(3) + c%parts(i) % r) * c % rh) + 1

       if(la(3) < 1) la(3) = 1
       if(la(3) > c%res_z) la(3) = c%res_z

       if(lb(3) > c%res_z) lb(3) = c%res_z

       call cnt_grid3D_addParticleToCells(c, i, la, lb)
    end do

    if(v) write(*,*) "Done."
  end subroutine cnt_grid3D_optimizeGrid
  !******


  function cnt_grid3D_cmpDropHeight(c, cPa, cPaIdx, p, lim2D_a, lim2D_b, iDrop) result(height)
    type(cnt_grid3D), intent(IN) :: c
    type(prt_sphere), intent(IN) :: cPa
    integer, intent(in) :: cPaIdx
    real(FD), dimension(:,:), intent(IN) :: p
    integer, dimension(:,:), intent(IN) :: lim2D_a, lim2D_b
    integer, intent(IN) :: iDrop

    type(prt_sphere), pointer :: oPa

    real(FD) :: height

    integer :: gx, gy, gz, iPart, nDown, np, k, x,y

    integer, dimension(:), pointer :: pList 

    integer, dimension(300) :: ppList

    real(FD) :: l  !! Squared distance of the sphere centers in (x,y)-plane.
    real(FD) :: d  !! Squared sum of the sphere radii.
    real(FD) :: h  !! Height of the new sphere center.
	
	integer :: maxh, th

    logical  :: iFound
 
    height = 0.0_fd
    iFound = .false.
    nDown  = 0

    ! Starting from the top, loop through the bins intersecting
    ! the xy-projection of the sphere searching for the resting
    ! height.
    !
	maxh = 0
    do y = lim2D_a(2, iDrop), lim2D_b(2, iDrop)
		gy = modulo(y-1, c%res_xy) + 1
       	do x = lim2D_a(1, iDrop), lim2D_b(1, iDrop)
       		gx = modulo(x-1, c%res_xy) + 1
			th = c % maxHeight_l(gx, gy)
			if (th > maxh) maxh = th
		end do
	end do

    do gz = maxh,1,-1
          do y = lim2D_a(2, iDrop), lim2D_b(2, iDrop)
			 gy = modulo(y-1, c%res_xy) + 1
             do x = lim2D_a(1, iDrop), lim2D_b(1, iDrop)
             	gx = modulo(x-1, c%res_xy) + 1
                pList => c % grid(gx,gy,gz) % plist
                np = c % grid(gx,gy,gz) % nParts

                k = 0
                do iPart = 1, np
                   if(cPaIdx > pList(iPart)) then
                      k = k + 1
                      ppList(k) = pList(iPart)
                   end if
                end do


                ! Loop through the particles in each bin.
                ! 

                do iPart = 1, k

                   oPa => c % parts( ppList(iPart) )

                   l = sum( ( p(1:2,iDrop) - oPa%P(1:2) )**2 )
                   d = cPa%rr + oPa%rr + 2.0 * cPa%r * oPa%r

                   if(l < d) then
                      iFound = .true.
                      
                      h = sqrt(d - l) + oPa % P(3)
                      if(height < h) then
                         height = h
                      end if
                   end if

                end do

             end do
          end do

       ! How many cells down do we check after the first intersection?
       !
       if(iFound) nDown  = nDown + 1
       if(nDown > 1) exit

    end do

    !! If no intersection is found, set height to particle radius.
    !! 
    if(.not. iFound) height = cPa%r

  end function cnt_grid3D_cmpDropHeight

  !****f* container_grid3D/cnt_grid3D_dropParticle
  ! NAME
  !   cnt_grid3D_dropParticle
  !
  ! DESCRIPTION
  !   
  !
  ! INPUTS
  !   c          : cnt_grid3D  ... Container
  !   sI         : integer     ... Sphere index
  !   nDrops     : integer     ... Number of drops.
  !
  ! RESULT
  !    cnt_grid3D_drop_sphere (integer)
  !
  ! SOURCE
  !
  integer function cnt_grid3D_dropParticle(c, sI, nDrops)
    type(cnt_grid3D), intent(INOUT) :: c
    integer, intent(IN)             :: sI, nDrops

    integer,  dimension(2, nDrops)  :: lim2D_a, lim2D_b
    integer,  dimension(3)          :: lim3D_a, lim3D_b
    integer                         :: iDrop, iMinHeight(1)

    real(FD), dimension(nDrops)     :: height
    real(FD), dimension(2,nDrops)   :: p
 
    type(prt_sphere), pointer       :: cPa, oPa


    cPa => c%parts(sI)

    !$omp parallel
    !$omp single

    call rnd_generate_uniform(-c%hWidth, c%hWidth, p(1,:))
    call rnd_generate_uniform(-c%hWidth, c%hWidth, p(2,:))
      
    lim2D_a = floor((p - cPa%r + c%hwidth) * c%rw) + 1
    lim2D_b = floor((p + cPa%r + c%hwidth) * c%rw) + 1


    !where(lim2D_a < 1)        lim2D_a = 1
    !where(lim2D_b > c%res_xy) lim2D_b = c%res_xy

    height = 0.0

    !$omp end single copyprivate(height, p, lim2d_a, lim2d_b)

    !$omp do private(iDrop) schedule(static)
    do iDrop = 1, nDrops
       height(iDrop) = cnt_grid3D_cmpDropHeight(c, cPa, sI, p, lim2D_a, lim2D_b, iDrop)
    end do
    !$omp end do

    !$omp workshare
    iMinHeight = minloc(height)
    !$omp end workshare


    !$omp sections
    cPa%P(1:2)   = p(:,iMinHeight(1))
    cPa%P(3)     = height(iMinHeight(1))
    !$omp section
    lim3D_a(1:2) = lim2D_a(:,iMinHeight(1))
    lim3D_b(1:2) = lim2D_b(:,iMinHeight(1))
    !$omp end sections
 
    !$omp single
    if(c%height < cPa%P(3)) then
       c%height = cPa%P(3) + cPa % r
    end if
    !$omp end single


    !$omp  sections
    lim3D_a(3) = floor((cPa%P(3) - cPa % r) * c % rh) + 1
    !$omp section
    lim3D_b(3) = floor((cPa%P(3) + cPa % r) * c % rh) + 1
    !$omp end sections

    !$omp single
    if(lim3D_a(3) < 1) lim3D_a(3) = 1
    if(lim3D_a(3) > c%res_z) lim3D_a(3) = c%res_z
    if(lim3D_b(3) > c%res_z) lim3D_b(3) = c%res_z

    call cnt_grid3D_addParticleToCells(c, sI, lim3D_a, lim3D_b)

    cnt_grid3D_dropParticle = 1

    !$omp end single
    !$omp end parallel
  end function cnt_grid3D_dropParticle
  !******


  !>   Subroutine to pack spherical particles inside the cnt_grid_2D
  !!   container. Calls cnt_grid3D_drop_sphere for each particle in
  !!   the container, and computes the statistics.
  !! 
  !!   \param[inout] c Container
  !!
  !! SOURCE
  !!
  subroutine cnt_grid3D_packParticles(c, nDrops, dropTable, timer_tot)
    type(cnt_grid3D)          :: c
    integer                   :: nDrops
    integer, dimension(:)     :: dropTable
    type(utl_timer), optional :: timer_tot

    integer         :: i, success
    type(utl_timer) :: timer

    write(*,'(A)') 'Packing media:'
    
    if(present(timer_tot)) then
       timer = timer_tot
       write(*,*) " Total CPU time"
    else
       call utl_timer_init(timer, 5.0_fd, c%nParts)
       write(*,*) " CPU time per medium"
    end if

    write(*,*) "  Completed          Elapsed                   Estimated"

    do i = 1, c%nParts
       success = cnt_grid3D_dropParticle(c, i, dropTable(1 + mod(i,size(dropTable))))
       call utl_timerIncrease(timer)
    end do

    call cnt_grid3D_optimizeGrid(c)

    if(present(timer_tot)) timer_tot = timer

  end subroutine cnt_grid3D_packParticles
  !******

  function cnt_grid3D_getHeight_xy(c, P) result(h)
    type(cnt_grid3D), intent(IN) :: c
    real(FD), dimension(2), intent(IN) :: p
    real(FD) :: h

    type(prt_sphere), pointer :: cPa
    integer, dimension(2) :: gp
    integer  :: gx, gy, gz, iPart
    real(FD) :: xy_sqr, ht

    gp = floor((P + c%hWidth) * c % rw) + 1

    gx = gp(1)
    gy = gp(2)

    ht = 0.0_fd
    h  = 0.0_fd

!    DO gz = c % maxHeight_l(gx, gy), 1, -1
    do gz = c % maxHeight_g, 1, -1
       do iPart = 1, c % grid(gx,gy,gz) % nParts

          cPa => c % parts( c % grid(gx,gy,gz) % plist(iPart) )
          xy_sqr = sum((cPa%P(1:2) - P)**2)

          if(xy_sqr < cPa%rr) then
             ht = sqrt(cPa%rr - xy_sqr) + cPa % P(3)
             if(ht > h) h = ht
          end if

       end do
       if(h>0.0_fd) exit
    end do

  end function cnt_grid3D_getHeight_xy


  !****f* container_grid3D/cnt_grid3D_export_rib
  ! NAME
  !   cnt_grid3D_export_rib
  !
  ! DESCRIPTION
  !   
  !
  ! INPUTS
  !   c : cnt_grid3D ... Container
  !
  !******
  subroutine cnt_grid3D_export_rib(c)
    type(cnt_grid3D) :: c
    integer          :: i
    
    !Open(1, FILE="particles.rib", status='REPLACE')
    !Write(1,*) 'Points "P" [', c%parts(:) % P, '] "width" [', 2.0 * c%parts%r, '] "uniform string type" "sphere"'
    !Close(1)

    open(1, FILE="particles.rib", status='REPLACE')
    open(2, FILE="spheres.dat", status='REPLACE')

    write(1,*) 'Points "P" [' 

    do i=1, c % nParts
       write(1,*) c%parts(i) % P
    end do

    write(1,*) '] "width" ['

    do i=1, c % nParts
       write(1,*) 2.0 * c%parts(i)%r
    end do

    write(1,*) '] "uniform string type" "sphere"'

    do i=1, c % nParts
       write(2,*) c%parts(i)%P(1), c%parts(i)%P(2), c%parts(i)%P(3), c%parts(i)%r
    end do

    close(1)
    close(2)

  end subroutine cnt_grid3D_export_rib


  !****f* container_grid3D/cnt_grid3D_cmp_porosity
  ! NAME
  !   cnt_grid3D_cmp_porosity
  !
  ! DESCRIPTION
  !   Computes the porosity of the particulate medium inside
  !   the container. Simply, the volume of the container - the volume of the
  !   spheres inside the container.
  !
  ! INPUTS
  !   c : cnt_grid3D ... Container
  !
  ! RESULT
  !    cnt_grid3D_cmp_porosity : real
  !
  !******
  real(fd) function cnt_grid3D_cmpDensity(c) result(rho)
    type(cnt_grid3D), intent(in) :: c

    real :: r, volParts, volCont

    r = 4.0/3.0 * 3.1416
    volParts = sum(r * c%parts%r**3)
    volCont  = c%width**2 * c%height
    rho = volParts / volCont

  end function cnt_grid3D_cmpDensity

  !****f* container_grid3D/cnt_grid3D_cmp_porosity_mc
  ! NAME
  !   cnt_grid3D_cmp_porosity_mc
  !
  ! DESCRIPTION
  !   Computes the porosity of the particulate medium inside
  !   the container using Monte-Carlo sampling. 
  !
  ! INPUTS
  !   c : cnt_grid3D ... Container
  !
  ! RESULT
  !    cnt_grid3D_cmp_porosity : real
  !
  !******
  real(fd) function cnt_grid3D_cmpDensityMC(c, nSamples) result(rho)
    type(cnt_grid3D), intent(in) :: c
    integer :: nSamples, i, x, y, z
    integer(il) :: pId

    type(rnd_gen) :: rgen
    real(fd) :: sample(3)

    !call rnd_init(rgen, 0)

    rho = 0.0_fd

    do z = 1, c % res_z
       do y = 1, c % res_xy
          do x = 1, c % res_xy

             do i=1, nSamples
                call rnd_generate_uniform(0.0_fd, 1.0_fd, sample)
                sample(1) = (real(x, fd) - sample(1)) * c%wr - c%hWidth
                sample(2) = (real(y, fd) - sample(2)) * c%wr - c%hWidth
                sample(3) = (real(z, fd) - sample(3)) * c%hr

                if(cnt_grid3D_isPointInsideParticle(c, sample, pId)) rho = rho + 1.0_fd
             end do

          end do
       end do
    end do

    rho = rho / real(nSamples * c%res_z * c%res_xy * c%res_xy, fd)

  end function cnt_grid3D_cmpDensityMC

  !****f* container_grid3D/cnt_grid3D_cmpPorosityStructure
  ! NAME
  !   cnt_grid3D_cmpPorosityStructure
  !
  ! DESCRIPTION
  !   Computes the porosity of the particulate medium inside
  !   the container using Monte-Carlo sampling. 
  !
  ! INPUTS
  !   c : cnt_grid3D ... Container
  !
  ! RESULT
  !    cnt_grid3D_cmp_porosity : real
  !
  !******
  subroutine cnt_grid3D_cmpDensityStructure(c, nSamples, rho)
    type(cnt_grid3D), intent(in) :: c
    integer :: nSamples, iSample, iHeight
    real(fd), dimension(:) :: rho

    real(fd)    :: sample(3), dz, hdz
    integer(il) :: pId

    dz  = c%height / real(size(rho), fd)
    hdz = 0.5_fd * dz

    rho = 0.0_fd
    do iHeight=1, size(rho)
       sample(3)   = iHeight * dz - hdz
       do iSample=1, nSamples
          call rnd_generate_uniform(0.0_fd, 1.0_fd, sample(1:2))
          sample(1:2) = sample(1:2) * c % width - c % hWidth

          if(cnt_grid3D_isPointInsideParticle(c, sample, pId)) rho(iHeight) = rho(iHeight) + 1.0_fd
       end do

       rho(iHeight) = 1.0_fd - rho(iHeight) / real(nSamples, fd)
    end do

  end subroutine cnt_grid3D_cmpDensityStructure

  !****f* container_grid3D/cnt_grid3D_quantizePosition
  ! NAME
  !   cnt_grid3D_quantizePosition
  !
  ! DESCRIPTION
  !   Returns the cell 
  !
  ! INPUTS
  !   c     : cnt_grid3D ... Container
  !   pos_r : real       ... Position of the point
  !
  !
  ! RESULT
  !   pos_g : real       ... Quantized position
  !
  ! SOURCE
  !PURE 
  function cnt_grid3D_quantizePosition(c, pos_r) result(pos_g)
    type(cnt_grid3D), intent(IN) :: c
    real(FD), dimension(3), intent(IN) :: pos_r
    integer, dimension(3) :: pos_g

    !pos_g(1:2) = FLOOR((pos_r(1:2) + c%hwidth) * c%rw) + 1
    pos_g(1)   = floor((pos_r(1) + c%hwidth) * c%rw) + 1
    pos_g(2)   = floor((pos_r(2) + c%hwidth) * c%rw) + 1
    pos_g(3)   = floor((pos_r(3) ) * c%rh) + 1
	pos_g(1) = max(min(pos_g(1), c%res_xy), 1)
	pos_g(2) = max(min(pos_g(2), c%res_xy), 1)
	pos_g(3) = max(min(pos_g(3), c%res_z), 1)

  end function cnt_grid3D_quantizePosition
  !******

  function cnt_grid3D_cmpHeightMap(c, res) result(hMap)
    type(cnt_grid3D) :: c
    integer :: res

    real(FD), dimension(res,res) :: hMap

    integer  :: i, j
    real(FD) :: dl, hdl

    dl  = c % width / res
    hdl = 0.5_fd * dl

    do i=1,res
       do j=1,res
          hMap(j,i) = cnt_grid3D_getHeight_xy(c, (/ j*dl - hdl - c%hwidth, i*dl - hdl - c%hwidth /))
       end do
    end do

  end function cnt_grid3D_cmpHeightMap

  subroutine cnt_grid3D_cmpStats(c, res, mean, std, hMin, hMax)
    type(cnt_grid3D), intent(IN) :: c
    integer, intent(IN) :: res

    real(FD) ::  mean
    real(FD), optional :: std, hMin, hMax
    real(FD), dimension(res,res) :: hMap

    hMap = cnt_grid3D_cmpHeightMap(c,res)

    mean = sum(hMap) / res**2

    if(present(std))  std  = sqrt(sum(hMap**2) / res**2 - mean**2)
    if(present(hMin)) hMin = minval(hMap)
    if(present(hMax)) hMax = maxval(hMap)

  end subroutine cnt_grid3D_cmpStats

end module container_grid3D
