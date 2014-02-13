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
!! \author Hannu Parviainen gannu@astro.helsinki.fi
!!
!!
!! \history
!!   16.11.2007
!!    - Modified the sphere type to include a pointer to a material.
!!
!! \date 04.09.2007
!!
module geometry
  use base

  !$ use omp_lib

  implicit none
  
  integer(is), parameter :: RAY_ACTIVE = 0
  integer(is), parameter :: RAY_EXIT_U = 1
  integer(is), parameter :: RAY_EXIT_D = 2
  integer(is), parameter :: RAY_EXIT_S = 3
  integer(is), parameter :: RAY_DEAD   = 4


  integer(is), parameter :: RAY_TYPE_ENERGY = 1
  integer(is), parameter :: RAY_TYPE_CAMERA = 2
  integer(is), parameter :: RAY_TYPE_SHADOW = 3
  integer(is), parameter :: RAY_TYPE_OTHER1 = 4 

  integer(is), parameter :: RAY_ID_INCR     = 4

  integer(il)            :: gmt_rayid(4) = [1,2,3,4]

  type :: ray
     real(fd), dimension(3) :: P = [ 0,  0,  0 ]
     real(fd), dimension(3) :: D = [ 0,  0, -1 ]

     real(fd) :: tStart = 1e-10_fd
     real(fd) :: tEnd   = 1e18_fd

     real(fd) :: I      = 1.0_fd

     real(fd) :: energy = 1.0_fd
     real(fd) :: lambda = PHS_H_J * PHS_C / (1.0_fd * eV_J)

     integer(il)  :: order   = 0_il
     integer(il)  :: rayid   = 0_il
     integer(is)  :: status  = RAY_ACTIVE
     integer(is)  :: rayType

  end type ray


  !> Spherical particle with location P, radius r and squared radius rr.
  !!
  type :: prt_sphere
     real(fd), dimension(3)      :: P
     real(fd)                    :: r  = 1.0_fd
     real(fd)                    :: rr = 1.0_fd
     real(fd)                    :: k  = 1.0_fd
     integer(il)                 :: rayID = 0_il 
  end type prt_sphere


  type :: intersection_geometry
     type(prt_sphere), pointer :: part
     real(fd), dimension(3)    :: N, P1, P2
     real(fd)                  :: l
     integer                   :: nIsects
  end type intersection_geometry

contains

  subroutine ray_init(r, rtype)
    type    ( ray ), intent(INOUT) :: r
    integer (  is ), intent(IN)    :: rType

    r%raytype = rtype
    r%rayID   = ray_get_newid(rtype)
    r%status  = RAY_ACTIVE
    r%order   = 0_il

  end subroutine ray_init

  integer(il) function ray_get_newid(rtype) result(rid)
    integer(is), intent(in) :: rtype

    !$omp critical
    rid              = gmt_rayid(rtype)
    gmt_rayid(rtype) = gmt_rayid(rtype) + RAY_ID_INCR
    !$omp end critical
  end function ray_get_newid


  pure real(fd) function vec_length(v)
    real(fd), dimension(:), intent(IN) :: v

    vec_length = sqrt(sum(v**2))

  end function vec_length

  pure real(fd) function vec_length_sqr(v)
    real(fd), dimension(:), intent(IN) :: v

    vec_length_sqr = sum(v**2)

  end function vec_length_sqr

  pure subroutine vec_normalize(v)
    real(fd), dimension(:), intent(INOUT) :: v

    v = v / sqrt(sum(v**2))
  end subroutine vec_normalize

  function vec_cart_random_hemispherical_uniform() result(v)
    real(fd), dimension(3) :: v
    real(fd) :: phi, r

    call random_number(phi)
    call random_number(v(3))

    phi  = phi * TWO_PI
    r    = sqrt(1.0_fd - v(3)**2)

    v(1) = r * cos(phi)
    v(2) = r * sin(phi)

  end function vec_cart_random_hemispherical_uniform

  function vec_cart_random_spherical_uniform() result(v)
    real(fd), dimension(3) :: v
    real(fd) :: phi, r

    call random_number(phi)
    call random_number(v(3))

    phi  = phi * TWO_PI
    r    = sqrt(1.0_fd - v(3)**2)
    v(1) = r * cos(phi)
    v(2) = r * sin(phi)

    call random_number(phi)

    if(phi < 0.5) v(3) = -v(3)

  end function vec_cart_random_spherical_uniform
 
end module geometry
