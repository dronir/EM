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

!****h* XR/trace_xr
! NAME
!   trace_xr
!
! DESCRIPTION
!
! NOTES
!
!
! AUTHOR
!   Olli Wilkman
!
! USES
!   Base
!   Particle
!
!
! CREATION DATE
!   25.03.2014
!******
module trace_xr
  use base
  use geometry
  use particle
  use random
  use container_grid3D
  use material
  use trace

  !$ USE omp_lib

  implicit none

contains
 


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



end module trace_xr
