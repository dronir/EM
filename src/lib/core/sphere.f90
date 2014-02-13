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

!>   Provides type definitions of particles and the routines to work with them.
!!
!!
!! HISTORY
!!   16.11.2007
!!    - Modified the sphere type to include a pointer to a material.
!!
!! \date 04.09.2007
!!
MODULE particle
  USE base
  USE geometry
!  USE material
!$  USE omp_lib


  IMPLICIT NONE
  
CONTAINS
  
  !>   Initializes a sphere.
  !!
  !!   \param[in]  s  Sphere
  !!   \param[in]  p  Sphere centre
  !!   \param[in]  r  Sphere radius
  !!
  SUBROUTINE prt_init_sphere(s, p, r)
    type(prt_sphere) :: s
    real(FD)     :: p(3), r

    s%P  = p
    s%r  = r
    s%rr = r**2

    s%rayID = 0_il

  END SUBROUTINE prt_init_sphere

  !>   Ray-sphere intersectin test.
  !!
  !! \param[in]   s  Sphere.
  !! \param[in]   r  Ray.
  !! \param[in]   i  Intersection geometry.
  !!
  LOGICAL FUNCTION prt_sph_trace(s, r, i)
    type(intersection_geometry), INTENT(inout) :: i
    type(prt_sphere), TARGET,    INTENT(inout) :: s
    type(ray),                   INTENT(inout) :: r

    real(FD) :: A, B, C, d, t

       s%rayID = r%rayID

       A = 1.0_fd !! |ray.D^2|
       B = 2.0_fd * SUM(r%D * (r%P - s%P))
       C = SUM((r%P - s%P)**2) - s%rr

       d = B*B - 4.0_fd * C

       IF (d < 0.0) THEN
          prt_sph_trace = .FALSE.
          i%nIsects     = 0
       ELSE
          d  = SQRT(d)
          t = (-B -d) * 0.5_fd

          i%nIsects = 2

          IF(t < 0.0_fd) THEN
             t = (-B + d) * 0.5_fd
             i%nIsects = 1
          END IF

          IF(t < 0.0_fd .OR. t > r%tEnd) THEN
             prt_sph_trace = .FALSE.
             i%nIsects     = 0
          ELSE
             r % tEnd = t
             i % P1   = r%P + r%D * t                  ! t1
             i % P2   = r%P + r%D * (-B + d) * 0.5_fd  ! t2
             i % N    = (i%P1 - s%P) / s%r

             i % part => s

             prt_sph_trace = .TRUE.
          END IF
       END IF

  END FUNCTION prt_sph_trace
  !******
  
  !****f* particle/prt_sph_trace_shadow
  ! NAME
  !   sph_trace_shadow
  !
  ! DESCRIPTION
  !   Ray-sphere intersectin test.
  !
  ! INPUTS
  !   s : sphere   ... Sphere.
  !   r : ray      ... Ray.
  !   i : intersection_geometry ... Intersection geometry.
  !
  ! SOURCE
  LOGICAL FUNCTION prt_sph_trace_shadow(s, r)
    type(prt_sphere),            INTENT(inout) :: s
    type(ray),                   INTENT(in) :: r

    real(FD) :: A, B, C, d

    s%rayID = r%rayID
    
    A  = 1.0_fd 
    B  = 2.0_fd * SUM(r%D * (r%P - s%P))
    C  = SUM((r%P - s%P)**2) - s%rr

    d  = B*B - 4.0_fd * C

    IF (d < 0.0_fd) THEN
       prt_sph_trace_shadow = .FALSE.
    ELSE
       IF(B > SQRT(d)) THEN
          prt_sph_trace_shadow = .FALSE.
       ELSE
          prt_sph_trace_shadow = .TRUE.
       END IF
    END IF

  END FUNCTION prt_sph_trace_shadow
  !******

  real(FD) FUNCTION prt_sph_evalMu(s, energy)
    type(prt_sphere), INTENT(in) :: s
    real(FD), INTENT(IN)         :: energy
    
    prt_sph_evalMu = 0.0_fd !mat_evalMuAbs(s%mat, energy)

    !write(*,*) prt_sph_evalMu

  END FUNCTION prt_sph_evalMu

  !****f* particle/prt_sphere_sphere_intersection
  ! NAME
  !   prt_sphere_sphere_intersection
  !
  ! DESCRIPTION
  !   Tests whether two spheres intersect.
  !
  ! INPUTS
  !   s1 : sphere   ... First sphere.
  !   s2 : sphere   ... Second sphere.
  !
  ! SOURCE
  FUNCTION prt_sphere_sphere_intersection(s1, s2) RESULT(isect)
    type(prt_sphere) :: s1, s2
    LOGICAL      :: isect

    IF(SUM((s1%P - s2%P)**2) > (s1%rr + s2%rr)) THEN
       isect = .FALSE.
    ELSE
       isect = .TRUE.
    END IF

  END FUNCTION prt_sphere_sphere_intersection
  !******

END MODULE particle
