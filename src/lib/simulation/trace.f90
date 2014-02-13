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

!****h* XR/trace
! NAME
!   trace
!
! DESCRIPTION
!
! NOTES
!
!
! AUTHOR
!   Hannu Parviainen
!
! USES
!   Base
!   Particle
!
!
! CREATION DATE
!   04.09.2007
!******
module trace
  use base
  use geometry
  use particle
  use random
  use container_grid3D
  use material

  !$ USE omp_lib

  implicit none

  type :: trc_gridTrace
     type(cnt_grid3D), pointer :: c
     type(ray),        pointer :: r 

     real(FD) :: dydx, dzdx, dxdy, dzdy, dxdz, dydz
     real(FD) :: tMaxX, tMaxY, tMaxZ
     real(FD) :: tDeltaX, tDeltaY, tDeltaZ
     real(FD) :: lx, ly, lz
     real(FD) :: l(3), step(3)
     
     integer  :: pGrid(3)
     integer  :: gMax(3)

  end type trc_gridTrace

  type :: trc_depthListElem
     type(trc_depthListElem), pointer :: prev => null()
     type(trc_depthListElem), pointer :: next => null()
     
     type(prt_sphere), pointer    :: part     => null()
     real(FD)                     :: t, l, mu 
  end type trc_depthListElem

contains
  
  subroutine trc_gridTrace_init(g, c, r)
    type(trc_gridTrace)      :: g
    type(cnt_grid3D), target :: c
    type(ray),        target :: r

    g%c => c
    g%r => r

    g%gMax(1:2) = c % res_xy
    g%gMax(3)   = c % res_z
    
    g%pGrid = cnt_grid3D_quantizePosition(c, r%P)

    where(r%D >= 0.0_fd)
       g % step = 1.0_fd
    elsewhere
       g % step = -1.0_fd
    end where

    if(g % step(1) > 0) then
       g % lx = g%pGrid(1) * c%wr - r%P(1) - c%hWidth
    else
       g % lx = r%P(1) + c%hWidth - (g%pGrid(1)-1) * c%wr 
    end if

    if(g % step(2) > 0) then
       g % ly = g%pGrid(2) * c%wr - r%P(2) - c%hWidth
    else
       g % ly = r%P(2) + c%hWidth - (g%pGrid(2)-1) * c%wr 
    end if

    if(g % step(3) > 0) then
       g % lz = g%pGrid(3) * c%hr - r%P(3)
    else
       g % lz = r%P(3) - (g%pGrid(3)-1) * c%hr
    end if

    g % tMaxX = sqrt((r%D(2)/r%D(1) * g%lx)**2 + (r%D(3)/r%D(1) * g%lx)**2 + g%lx**2)
    g % tMaxY = sqrt((r%D(1)/r%D(2) * g%ly)**2 + (r%D(3)/r%D(2) * g%ly)**2 + g%ly**2)
    g % tMaxZ = sqrt((r%D(1)/r%D(3) * g%lz)**2 + (r%D(2)/r%D(3) * g%lz)**2 + g%lz**2)

    g % tDeltaX = sqrt((r%D(2)/r%D(1) * c%wr)**2 + (r%D(3)/r%D(1) * c%wr)**2 + c%wr**2)
    g % tDeltaY = sqrt((r%D(1)/r%D(2) * c%wr)**2 + (r%D(3)/r%D(2) * c%wr)**2 + c%wr**2)
    g % tDeltaZ = sqrt((r%D(1)/r%D(3) * c%hr)**2 + (r%D(2)/r%D(3) * c%hr)**2 + c%hr**2)

    if(g%tMaxX /= g%tMaxX) g%tMaxX = huge(g%tMaxX)
    if(g%tMaxY /= g%tMaxY) g%tMaxY = huge(g%tMaxY)
    if(g%tMaxZ /= g%tMaxZ) g%tMaxZ = huge(g%tMaxZ)

    if(g%tDeltaX /= g%tDeltaX) g%tDeltaX = huge(g%tDeltaX)
    if(g%tDeltaY /= g%tDeltaY) g%tDeltaY = huge(g%tDeltaY)
    if(g%tDeltaZ /= g%tDeltaZ) g%tDeltaZ = huge(g%tDeltaZ)

  end subroutine trc_gridTrace_init

  logical function trc_gridTrace_step(g, x, y, z, outAxis) result(inside)
    type(trc_gridTrace) :: g
    integer :: x, y, z
    integer, optional :: outAxis

    inside = .true.

    if(g% tMaxX < g% tMaxY) then
       if(g% tMaxX < g% tMaxZ) then

          X     = X + g% step(1)

          if(x > g%gMax(1) .or. x < 1) then
             if(present(outAxis)) outAxis = 1
             inside = .false.
             return
          end if
          g % tMaxX = g% tMaxX + g% tDeltaX;
       else

          Z     = Z + g% step(3)

          if(z > g%gMax(3) .or. z < 1) then
             if(present(outAxis)) outAxis = 3
             inside = .false.
             return
          end if
          g % tMaxZ = g% tMaxZ + g% tDeltaZ;
       end if
    else
       if(g% tMaxY < g% tMaxZ) then

          Y     = Y + g% step(2)

          if(y > g%gMax(2) .or. y < 1) then
             if(present(outAxis)) outAxis = 2
             inside = .false.
             return
          end if
          g % tMaxY = g% tMaxY + g% tDeltaY;
       else 

          Z     = Z + g% step(3)

          if(z > g%gMax(3) .or. z < 1) then
             if(present(outAxis)) outAxis = 3
             inside = .false.
             return
          end if
          g % tMaxZ = g% tMaxZ + g% tDeltaZ;
       end if
    end if

  end function  trc_gridTrace_step

  logical function trc_gridTrace_stepPeriodic(g, x, y, z) result(isInside)
    type(trc_gridTrace) :: g
    integer :: x, y, z
    integer :: outAxis

    integer :: pGrid(3)

    real(fd) :: ds

  
    isInside = trc_gridTrace_step(g,x,y,z, outAxis)

    if(.not. isInside) then

       if(outAxis == 3) then
          return
       else if(outAxis == 1) then

          if(g%r%D(1) > 0) then
             ds = (g%c%hWidth - g%r%P(1))

             g%r%P(3) =  g%r%P(3) + ds * (g%r%D(3) / abs(g%r%D(1)))
             g%r%P(2) =  g%r%P(2) + ds * (g%r%D(2) / abs(g%r%D(1)))
             g%r%P(1) = - g%c%hwidth + TRACE_EPS
          else

             ds = (g%c%hWidth + g%r%P(1))

             g%r%P(3) =  g%r%P(3) + ds * (g%r%D(3) / abs(g%r%D(1)))
             g%r%P(2) =  g%r%P(2) + ds * (g%r%D(2) / abs(g%r%D(1)))
             g%r%P(1) =   g%c%hWidth - TRACE_EPS
 
          end if
       else if(outAxis == 2) then
         
          if(g%r%D(2) > 0) then

             ds = (g%c%hWidth - g%r%P(2))

             g%r%P(3) =  g%r%P(3) + ds * (g%r%D(3) / abs(g%r%D(2)))
             g%r%P(1) =  g%r%P(1) + ds * (g%r%D(1) / abs(g%r%D(2)))
             g%r%P(2) = - g%c%hwidth + TRACE_EPS
          else
             ds = (g%c%hWidth + g%r%P(2))
             g%r%P(3) =  g%r%P(3) + ds * (g%r%D(3) / abs(g%r%D(2)))
             g%r%P(1) =  g%r%P(1) + ds * (g%r%D(1) / abs(g%r%D(2)))
             g%r%P(2) =   g%c%hWidth - TRACE_EPS
          end if

       end if

       if(g%r%P(3) < 0.0_fd .or. g%r%P(3) > g%c%height) then

          isInside = .false.
          return
       end if

       if(.not. cnt_grid3D_isPointInsideParticle(g%c, g%r%p)) then

          pGrid = cnt_grid3D_quantizePosition(g%c, g%r%P)
          x = pGrid(1)
          y = pGrid(2)
          z = pGrid(3)

          g%r%rayID = g%r%rayID + RAY_ID_INCR

          call trc_gridTrace_init(g, g%c, g%r)
          isInside = .true.
       end if

    end if


  end function  trc_gridTrace_stepPeriodic


  subroutine trc_addDepthListElem(root, s, t, l, mu)
    type(trc_depthListElem), pointer :: root
    type(prt_sphere),        pointer :: s
    real(FD)                         :: t      ! The distance to the intersection point.
    real(FD)                         :: l      ! The distance traveled inside the particle.
    real(FD)                         :: mu     ! The linear attenuation coefficient for ray with energy e.

    type(trc_depthListElem), pointer :: curr
    type(trc_depthListElem), pointer :: new_e

    allocate(new_e)
    new_e%part => s
    new_e%t    =  t
    new_e%l    =  l
    new_e%mu   =  mu

    if(.not. associated(root)) then
       root      => new_e
       root%prev => null()
       root%next => null()
       return
    end if

    curr => root

    if(t < curr%t) then
       do
          if(associated(curr%prev)) then
             curr => curr%prev
             if(t > curr%t) then
                new_e%prev      => curr
                new_e%next      => curr%next
                curr%next%prev  => new_e
                curr%next       => new_e
                exit
             end if
          else
             new_e%next => curr
             new_e%prev => null()
             curr%prev  => new_e
             exit
          end if
       end do
    else
       do
          if(associated(curr%next)) then
             curr => curr%next
             if(t < curr%t) then
                new_e%prev      => curr%prev
                new_e%next      => curr
                curr%prev%next  => new_e
                curr%prev       => new_e
                exit
             end if
          else
             new_e%prev => curr
             new_e%next => null()
             curr%next  => new_e
             exit
          end if
       end do
    end if


  end subroutine trc_addDepthListElem

  subroutine trc_setFirstElem(e)
    type(trc_depthListElem), pointer :: e

    if(associated(e)) then 
       do
          if(associated(e%prev)) then
             e => e%prev
          else
             exit
          end if
       end do
    end if

  end subroutine trc_setFirstElem

  subroutine trc_deleteElemTree(e)
    type(trc_depthListElem), pointer :: e

    call trc_setFirstElem(e)

    if(associated(e)) then 
       do
          if(associated(e%next)) then
             e => e%next
             deallocate(e % prev)
          else
             exit
          end if
       end do

       deallocate(e)
       nullify(e)

    end if

  end subroutine trc_deleteElemTree

!!$
!!$  SUBROUTINE trc_nextElem(e)
!!$    type(trc_depthListElem), POINTER :: e
!!$
!!$    IF(ASSOCIATED(e%next)) THEN
!!$       e => e%next
!!$    ELSE
!!$       EXIT
!!$    END IF
!!$ 
!!$  END SUBROUTINE trc_nextElem

  real(FD) function trc_traceDumb(c, r) result(I)
    type(cnt_grid3d) :: c
    type(ray)        :: r

    type(intersection_geometry) :: isect

    integer :: k, x, y, z
    logical :: isFound

    isFound = .false.
    r % tEnd = TRACE_HUGE
    I = 0.0_fd

    do k = 1, c % nParts
       isFound = ( prt_sph_trace(c % parts(k), r, isect) .or. isFound)
    end do

    if(isFound) then
       call vec_normalize(isect%N)
       !I = 1.0 - SQRT(SUM((isect % P1 - r % P)**2)) / SQRT(2.0 * c%width**2) 
       I = dot_product(isect % N, -r%D)
    end if

  end function trc_traceDumb

  logical function trc_traceNearest(c, r, iSect) result(isFound)
    type(cnt_grid3d) :: c
    type(ray)        :: r
    type(intersection_geometry) :: iSect
    type(trc_gridTrace) :: gTrace

    integer(IL) :: k

    integer  :: pGrid(3)
    integer  :: x, y, z

    pGrid = cnt_grid3D_quantizePosition(c, r%P)

    x = pGrid(1)
    y = pGrid(2)
    z = pGrid(3)

    isFound = .false.
    r % tEnd = TRACE_HUGE

    do k = 1, c % grid(x, y, z) % nParts
       if(c%parts(c % grid(x, y, z) % plist(k)) % rayID /= r % rayID) then
          isFound = (prt_sph_trace(c%parts(c % grid(x, y, z) % plist(k)), r, isect) .or. isFound)
       end if
    end do

    !! No intersection found in the starting cell, starting the grid-tracing
    !! routine.
    if(.not. isFound) then
       call trc_gridTrace_init(gTrace, c, r)
       do
          if(.not. trc_gridTrace_stepPeriodic(gTrace, x, y, z)) exit

          do k = 1, c % grid(x, y, z) % nParts
             if(c%parts(c % grid(x, y, z) % plist(k)) % rayID /= r % rayID) then
                isFound = (prt_sph_trace(c%parts(c % grid(x, y, z) % plist(k)), r, isect) .or. isFound)
             end if
          end do

          if(isFound) then
             exit
          end if

       end do
    end if

  end function trc_traceNearest

  logical function trc_traceOcclusion(c, r) result(isFound)
    type(cnt_grid3d) :: c
    type(ray)        :: r
    type(trc_gridTrace) :: gTrace

    integer :: k

    integer  :: pGrid(3)
    integer  :: x, y, z
   
    isFound = .false.

    pGrid = cnt_grid3D_quantizePosition(c, r%P)

    x = pGrid(1)
    y = pGrid(2)
    z = pGrid(3)

    isFound = .false.
    r % tEnd = TRACE_HUGE

    do k = 1, c % grid(x, y, z) % nParts !kMax
       if(c%parts(c % grid(x, y, z) % plist(k)) % rayID /= r % rayID) then
          isFound = prt_sph_trace_shadow(c%parts(c % grid(x, y, z) % plist(k)), r)
          if(isFound) return
       end if
    end do

    call trc_gridTrace_init(gTrace, c, r)

    do
       if(.not. trc_gridTrace_stepPeriodic(gTrace, x, y, z)) exit 

       do k = 1, c % grid(x, y, z) % nParts
          if(c%parts(c % grid(x, y, z) % plist(k)) % rayID /= r % rayID) then
             isFound = prt_sph_trace_shadow(c%parts(c % grid(x, y, z) % plist(k)), r)
             if(isFound) return
          end if
       end do

    end do

  end function trc_traceOcclusion


  real(FD) function trc_traceOpticalDepth(c, m, r) result(oDepth)
    type(cnt_grid3d)    :: c
    type(mat_material)  :: m
    type(ray)           :: r
    type(trc_gridTrace) :: gTrace

    type(intersection_geometry) :: isect

    type(prt_sphere), pointer :: cPart

    integer(IL) :: k

    real(FD) :: muTotal

    integer  :: pGrid(3)
    integer  :: x, y, z

    logical  :: isFound

    pGrid = cnt_grid3D_quantizePosition(c, r%P)

    x = pGrid(1)
    y = pGrid(2)
    z = pGrid(3)

    oDepth  = 0.0_fd
    isFound = .false.
    r % tEnd = TRACE_HUGE

    r % rayID  = r % rayID + RAY_ID_INCR

    !! An optimization applicable if we have only a single material, must be adjusted
    !! if multiple materials (heterogeneous media) are needed.
    !!
    call mat_evalMu(m, r%energy, muTotal = muTotal)

    do k = 1, c % grid(x, y, z) % nParts
       cPart =>  c%parts(c % grid(x, y, z) % plist(k))
       if(cPart % rayID /= r % rayID) then
          isFound = prt_sph_trace(cPart, r, iSect)
          r % tEnd = TRACE_HUGE
          if(isFound) then
             if(iSect%nIsects == 2) then
                oDepth = oDepth + sqrt(sum((iSect%P1 - iSect%P2)**2)) * muTotal
             else
                oDepth = oDepth + sqrt(sum((r%P - iSect%P1)**2)) * muTotal
             end if
          end if
       end if
    end do

    call trc_gridTrace_init(gTrace, c, r)
    do
       if(.not. trc_gridTrace_step(gTrace, x, y, z)) exit

       do k = 1, c % grid(x, y, z) % nParts
          cPart =>  c%parts(c % grid(x, y, z) % plist(k))
          if(cPart % rayID /= r % rayID) then
             isFound = prt_sph_trace(cPart, r, iSect)
             r % tEnd = TRACE_HUGE
             if(isFound) then
                oDepth = oDepth + sqrt(sum((iSect%P1 - iSect%P2)**2)) * muTotal
             end if
          end if
       end do

    end do

  end function trc_traceOpticalDepth

  real(FD) function trc_tracePhysicalDepth(c, r) result(depth)
    type(cnt_grid3d)    :: c
    type(ray)           :: r
    type(trc_gridTrace) :: gTrace

    type(intersection_geometry) :: isect

    type(prt_sphere), pointer :: cPart

    integer(IL) :: k

    integer  :: pGrid(3)
    integer  :: x, y, z

    logical  :: isFound

    pGrid = cnt_grid3D_quantizePosition(c, r%P)

    x = pGrid(1)
    y = pGrid(2)
    z = pGrid(3)

    depth    = 0.0_fd
    isFound  = .false.
    r % tEnd = TRACE_HUGE
    r % rayID  = r % rayID + RAY_ID_INCR
    
    do k = 1, c % grid(x, y, z) % nParts

       cPart =>  c%parts(c % grid(x, y, z) % plist(k))
       if(cPart % rayID /= r % rayID) then
          isFound = prt_sph_trace(cPart, r, iSect)
          r % tEnd = TRACE_HUGE
          if(isFound) then
             if(iSect%nIsects == 2) then
                depth = depth + sqrt(sum((iSect%P1 - iSect%P2)**2))
             else
                depth = depth + sqrt(sum((r%P - iSect%P1)**2))
             end if
          else
          end if
       end if
    end do

    call trc_gridTrace_init(gTrace, c, r)
    do
       if(.not. trc_gridTrace_step(gTrace, x, y, z)) exit

       do k = 1, c % grid(x, y, z) % nParts
          cPart =>  c%parts(c % grid(x, y, z) % plist(k))
          if(cPart % rayID /= r % rayID) then
             isFound = prt_sph_trace(cPart, r, iSect)
             r % tEnd = TRACE_HUGE
             if(isFound) then
                depth = depth + sqrt(sum((iSect%P1 - iSect%P2)**2)) 
             end if
          end if
       end do
    end do


  end function trc_tracePhysicalDepth

  subroutine trc_traceRayToOpticalDepth(c, r, gamma, mat, matMu)
    type(cnt_grid3d) :: c
    type(ray)        :: r
    real(FD)         :: gamma

    type(mat_material), optional :: mat
    real(fd),           optional :: matMu

    type(trc_gridTrace)                  :: gTrace
    type(intersection_geometry)          :: isect
    type(prt_sphere),            pointer :: cPart
    type(trc_depthListElem),     pointer :: root

    real(FD)    :: oDepth, oDepth_t, s, t, muTotal

    integer(IL) :: k
    integer     :: pGrid(3)
    integer     :: x, y, z

    logical     :: isFound

    cPart => null()
    root  => null()

    pGrid = cnt_grid3D_quantizePosition(c, r%P)
 
    x = pGrid(1)
    y = pGrid(2)
    z = pGrid(3)

    if (z > c%res_z) then
       call utl_fatal_error("Ray out of accelerating structure in trc_traceRayToOpticalDepth.")
    end if

    oDepth  = 0.0_fd
    isFound = .false.
    r % tEnd = TRACE_HUGE
    r % rayID  = r % rayID + RAY_ID_INCR


    if(present(mat)) then
       call mat_evalMu(mat, r%energy, muTotal = muTotal)
    else if(present(matMu)) then
       muTotal = matMu
    else
       call utl_fatal_error("Error: no mu or mat given to trc_traceRayToOpticalDepth")
    end if

    !! We start by computing the total optical depth oDepth for the starting cell.
    !! Particles intersected by the ray are added to the linked list and arranged
    !! by the intersection distance. If the total optical depth for the cell exceeds
    !! the given optical depth, we continue directly to search the exact distance
    !! traveled by the ray, else we continue to trace trough the cells until the 
    !! given optical depth has been reached.
    !!

    do k = 1, c % grid(x, y, z) % nParts

       cPart =>  c%parts(c % grid(x, y, z) % plist(k))

       if(cPart % rayID /= r % rayID) then
          isFound = prt_sph_trace(cPart, r, iSect)
          r % tEnd = TRACE_HUGE
          if(isFound) then
             if(iSect%nIsects == 2) then
                s = sqrt(sum((iSect%P1 - iSect%P2)**2))
             else
                s = sqrt(sum((r%P - iSect%P1)**2))
             end if

             !!mu     = mat_evalMuAbs(mat, r%energy)
             oDepth = oDepth + s * muTotal
             call trc_addDepthListElem(root, cPart, sum((r%P - iSect%P1)**2), s, muTotal)
            
          end if
       end if
    end do


    !! If optical depth is smaller than the given value, continue trace.
    !!
    if(gamma > oDepth) then
       call trc_gridTrace_init(gTrace, c, r)
       do
          if(.not. trc_gridTrace_stepPeriodic(gTrace, x, y, z)) exit

          do k = 1, c % grid(x, y, z) % nParts
             cPart =>  c%parts(c % grid(x, y, z) % plist(k))
              if(cPart % rayID /= r % rayID) then
                isFound = prt_sph_trace(cPart, r, iSect)
                r % tEnd = TRACE_HUGE
                if(isFound) then
                   s      = sqrt(sum((iSect%P1 - iSect%P2)**2))
                   !!mu     = mat_evalMuAbs(mat, r%energy)
                   oDepth = oDepth + s * muTotal
                   call trc_addDepthListElem(root, cPart, sum((r%P - iSect%P1)**2), s, muTotal)
                end if
             end if
          end do

          if(gamma < oDepth ) exit
       end do

       !! The ray goes out of bounds without reaching the given optical depth.
       !! We mark it as exited, and pass the exact search.
       !!
       if(gamma > oDepth) then
          if(r%D(3) < 0.0_fd) then
             r % status = RAY_EXIT_D
          else
             r % status = RAY_EXIT_U
          end if
       end if

    end if

    !if(oDepth < 1e-6) then
    !   print *, r%D
    !end if

    !! Next, we find the value of t = gamma / mu by traversing the
    !! found intersections starting from the nearest.
    !!
    oDepth = 0.0_fd
    if(associated(root) .and. r%status == RAY_ACTIVE) then 

       call trc_setFirstElem(root)
       do
          !! Compute the maximum optical depth for the ray intersecting
          !! the particle.
          !!
          !!
          !!    |#######|
          !!    |---l---|
          !!    |#######|
          !!
          !! gamma = l * mu
          !!
          oDepth_t = root%l * root%mu

          !! When we find the particle p where the maximum optical depth
          !! is larger than the wanted optical depth, we compute the 
          !! distance to travel inside the particle corresponding the
          !! given optical depth.
          !!
          !!   |#######|    |#####|    |#######|
          !!   |#*-----|----|-----|----|--->###|
          !!   |#######|    |#####|    |-l-|###|
          !!   
          !!     |----------- t -----------|
          !!
          !!          l = g_max,p / mu
          !!
          !! The value for t is now the sum of l and the closest intersection distance (p_i1)
          !! with p.
          !!
          !!          t = dist_p_i1 + l
          !!
          if(oDepth + oDepth_t < gamma) then
             oDepth = oDepth + oDepth_t
          else
             t = sqrt(root%t) + (gamma - oDepth) / root%mu
             exit
          end if

          if(associated(root%next)) then
             root => root%next
          else
             exit
          end if

       end do

       r%P = r%P + t*r%D
       
       if(r%P(1) >  c%hWidth) r%P(1)  =  modulo(r%P(1), c%hWidth) - c%hWidth
       if(r%P(1) < -c%hWidth) r%P(1)  = -modulo(r%P(1), c%hWidth) + c%hWidth
       if(r%P(2) >  c%hWidth) r%P(2)  =  modulo(r%P(2), c%hWidth) - c%hWidth
       if(r%P(2) < -c%hWidth) r%P(2)  = -modulo(r%P(2), c%hWidth) + c%hWidth
       if(r%P(3) <      0.0) r%status = RAY_EXIT_D
       if(r%P(3) > c%height) r%status = RAY_EXIT_U

    end if


    if(.not. associated(root)) then
       r % status = RAY_EXIT_D
    end if

    call trc_deleteElemTree(root)

  end subroutine trc_traceRayToOpticalDepth

  recursive subroutine trc_gatherRadiance(G, E, L, P, N, W, nSamples, nOrders, o, Il, albedo, brdf, Pf, Pp)
    type(cnt_grid3d), intent(inout)        :: G
    real(fd), dimension(3), intent(IN)     :: E, L, P, N
    real(fd), dimension(nOrders)           :: Il
    real(fd), intent(IN)                   :: W
    integer,  intent(IN)                   :: nSamples(nOrders), nOrders, o
    real(fd)                               :: albedo
    real(fd), external                     :: brdf

    real(fd), external                     :: Pf
    real(fd)                               :: Pp(10)

  
    type(intersection_geometry) :: iSect
    type(ray) :: r, rS
    integer   :: i, oo

    real(fd)  :: ww, dw

    logical   :: pFound, pLit

    call ray_init(rS, RAY_TYPE_SHADOW)

    rS%D      = L
    rS%P      = P

    pLit      = .not. trc_traceOcclusion(G, rS)

    if(pLit) then
       !$omp atomic
       Il(o) = Il(o) + W * brdf(L, E, N, albedo, Pf, Pp)
    end if

    if((o < nOrders) .and. (W > 1e-9)) then
       oo  = o + 1

       dw  = TWO_PI / real(nSamples(oo), fd)

       call ray_init(r,  RAY_TYPE_CAMERA)

       do i=1, nSamples(oo)

          do
             r%D = vec_cart_random_spherical_uniform()
             ww  = brdf(r%D, E, N, albedo, Pf, Pp)
             if(ww > 0.0_fd) exit
          end do

          ww        = dw * W * ww
          r%P       = P
  
          pFound    = trc_traceNearest(G, r, iSect)
  
          if(pFound) then 
             call trc_gatherRadiance(G, r%D, L, iSect%P1 + TRACE_EPS * iSect%N, iSect%N, ww, &
                  & nSamples, nOrders, oo, Il, albedo, brdf, Pf, Pp)
          end if

          !$omp atomic
          r%rayID   = r%rayID + RAY_ID_INCR

       end do
    end if

  end subroutine trc_gatherRadiance

  !> Finds the true medium surface intersection point for a ray coming from direction d to
  !! the point (x,y,z).
  !!
  !! Given point p and direction d, traces the intersection point of a ray coming from the
  !! upper boundary of the medium.
  !!
  !! \param[in]   c      The acceleration grid.     
  !! \param[in]   p      The intersection point (x,y) of the ray in the mean medium plane
  !! \param[in]   d      The direction vector of the incident ray
  !! \param[out]  isect  Intersection geoemtry
  !!
  logical function trc_findSurfaceIntersection(c, p, d, isect) result(isFound)
    type(cnt_grid3d), intent(in) :: c
    real(fd),         intent(in) :: p(3)
    real(fd),         intent(in) :: d(3)

    type(intersection_geometry), intent(out) :: iSect

    type(ray) :: r

    !! Compute the incident ray starting (x,y) position at the top of the bounding box of the medium.
    !!
    call ray_init(r, RAY_TYPE_ENERGY)
    r % D      = d
    r % P(3)   = c%height - TRACE_EPS
    r % P(1:2) = p(1:2) + (r%P(3) - p(3)) * (r%D(1:2) / r%D(3))
    r % P(1:2) = modulo(r%P(1:2)+c%hWidth, c%width) - c%hWidth

    isFound = trc_traceNearest(c, r, isect)

  end function trc_findSurfaceIntersection

end module trace
