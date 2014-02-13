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
!! Copyright (C) 2009 Hannu Parviainen
!!
!! Contributor(s): Hannu Parviainen
!!
!!--- END GPL ---

!> Module containing routines to compute random variates following different distributions.
!!
!! \date 04.09.2007
!!
module random
!$ use omp_lib
   use base
   use util
   use stdtypes
   use mtprng
   use igdev
   implicit none

   private :: rnd_state

  !> A front-end for the ACML random number generator. Uses Mersenne Twister as a base generator.
  !!
  type rnd_gen
     type(mtprng_state) :: state_in
     integer            :: state(640)
     integer            :: seed
     integer            :: info
     logical            :: initialized
  end type rnd_gen

  logical                                          :: rnd_initialized = .false.
  type(rnd_gen), dimension(:), allocatable, save   :: rnd_state

contains

  !> Initialises the random number generator.
  !!
  subroutine rnd_init(seed, n_threads, verbose, gen)
    integer, intent(in)     :: seed
    integer, intent(in)     :: n_threads
    logical, intent(in)     :: verbose
    type(rnd_gen), optional :: gen

    integer :: s, i

    i = 1
    s = seed

    if(.not. rnd_initialized) then
       allocate(rnd_state(n_threads))
       rnd_initialized = .true.
    end if

#ifdef WITH_ACML
    if(present(gen)) then
       gen % seed = seed
       call drandinitialize(3, 0, [gen%seed], 1, gen%state, 640, gen%info)
    else
       !$ do i = 1, n_threads
         if (verbose) call utl_message("Initializing random number generator.")
         rnd_state(i) % seed = s
         call drandinitialize(3, 0, [s], 1, rnd_state(i)%state, 640, rnd_state(i)%info)
         !$ s = s + 1
       !$ end do
    end if
#else
    if(present(gen)) then
       gen % seed = seed
       call mtprng_init(seed, gen%state_in)
    else
       !$ do i = 1, n_threads
         if (verbose) call utl_message("Initializing random number generator.")
         rnd_state(i) % seed = s
         call mtprng_init(s, rnd_state(i)%state_in)
         !$ s = s + 1
       !$ end do
    end if
#endif
  end subroutine rnd_init


  !> Completely useless subroutine to assign constant values to a real array.
  !!
  !! \param[in]  a  Value
  !! \param[out] r  Output array
  !!
  subroutine rnd_generate_constant(a, r)
    real(fd), dimension(:), intent(OUT) :: r
    real(fd), intent(IN) :: a

    r = a
  end subroutine rnd_generate_constant


  !> Fills the given array with random values drawn from the uniform distribution [0 .. 1].
  !!
  !! \param[in]  gen  The random number generator
  !! \param[out] r    Output array
  !!
  subroutine rnd_generate_uniform_n(r)
    real(fd), dimension(:), intent(OUT) :: r
    integer                             :: i

    integer                             :: thread

    thread = 1
    !$ thread = omp_get_thread_num() + 1

#ifdef WITH_ACML
    call dranduniform(size(r), 0.0_fd, 1.0_fd, rnd_state(thread)%state, r, rnd_state(thread)%info)
#else
    do i = 1, size(r)
       r(i) = mtprng_rand_real1(rnd_state(thread)%state_in)
    end do
#endif
  end subroutine rnd_generate_uniform_n


  !> Draws a random value from the uniform distribution [0 .. 1].
  !!
  !! \param[in]  gen  The random number generator
  !!
  real(fd) function rnd_uniform_n() result(r)
    real(fd)                            :: rt(1)
    integer                             :: thread

    thread = 1
    !$ thread = omp_get_thread_num() + 1

#ifdef WITH_ACML
    call dranduniform(1, 0.0_fd, 1.0_fd, rnd_state(thread)%state, rt, rnd_state(thread)%info)
    r = rt(1)
#else
    r = mtprng_rand_real1(rnd_state(thread)%state_in)
#endif

  end function rnd_uniform_n


  !> Fills the given array with random values drawn from the uniform distribution [vmin .. vmax].
  !!
  !! \param[in]  gen  The random number generator
  !! \param[in]  vmin Minimum value
  !! \param[in]  vmax Maximum value
  !! \param[out] r    Output array
  !!
  subroutine rnd_generate_uniform(vmin, vmax, r)
    real(fd), intent(IN)                :: vmin
    real(fd), intent(IN)                :: vmax
    real(fd), dimension(:), intent(OUT) :: r

    integer                             :: i
    integer                             :: thread

    thread = 1
    !$ thread = omp_get_thread_num() + 1

#ifdef WITH_ACML
    call dranduniform(size(r), vmin, vmax, rnd_state(thread)%state, r, rnd_state(thread)%info)
#else
    do i = 1, size(r)
       r(i) = mtprng_rand_real1(rnd_state(thread)%state_in)
    end do
    r = r * (vMax - vMin) + vMin
#endif
  end subroutine rnd_generate_uniform


  !> Draws a random value from the uniform distribution [vmin .. vmax].
  !!
  !! \param[in]  gen  The random number generator
  !! \param[in]  vmin Minimum value
  !! \param[in]  vmax Maximum value
  !!
  real(fd) function rnd_uniform(vmin, vmax) result(r)
    real(fd), intent(IN)                :: vmin
    real(fd), intent(IN)                :: vmax
    real(fd)                            :: rt(1)

    integer                             :: thread

    thread = 1
    !$ thread = omp_get_thread_num() + 1

#ifdef WITH_ACML
    call dranduniform(1, vmin, vmax, rnd_state(thread)%state, rt, rnd_state(thread)%info)
    r = rt(1)
#else
    r = mtprng_rand_real1(rnd_state(thread)%state_in) * (vMax - vMin) + vMin
#endif    

  end function rnd_uniform


  !> Fills the given array with random values drawn from the exponential distribution.
  !!
  !! \param[in]  gen    The random number generator
  !! \param[in]  lambda Lambda
  !! \param[out] r      Output array
  !!
  subroutine rnd_generate_exponential(lambda, r)
    real(fd), intent(IN)                :: lambda
    real(fd), dimension(:), intent(OUT) :: r

    integer                             :: i
    integer                             :: thread

    thread = 1
    !$ thread = omp_get_thread_num() + 1

#ifdef WITH_ACML
    call drandexponential(size(r), lambda, rnd_state(thread)%state, r, rnd_state(thread)%info)
#else
    do i = 1, size(r)
       r(i) = -log(mtprng_rand_real1(rnd_state(thread)%state_in)) / lambda
    end do
#endif

  end subroutine rnd_generate_exponential


  !> Draws a random value from the exponential distribution.
  !!
  !! \param[in]  gen    The random number generator
  !! \param[in]  lambda Lambda
  !!
  real(fd) function rnd_exponential(lambda) result(r)
    real(fd), intent(IN)                :: lambda
    real(fd)                            :: rt(1)

    integer                             :: thread

    thread = 1
    !$ thread = omp_get_thread_num() + 1

#ifdef WITH_ACML
    call drandexponential(1, lambda, rnd_state(thread)%state, rt, rnd_state(thread)%info)
    r = rt(1)
#else
    r = -log(mtprng_rand_real1(rnd_state(thread)%state_in)) / lambda
#endif    

  end function rnd_exponential


  !> Fills the given array with random values drawn from the lognormal distribution.
  !!
  !! \param[in]  gen    The random number generator
  !! \param[in]  mean   Mean of the underlying Gaussian distribution
  !! \param[in]  std    Standard deviation of the underlying Gaussian distribution.
  !! \param[out] r      Output array
  !!
  subroutine rnd_generate_logNormal(mean, std, r)
    real(fd), intent(IN)                :: mean
    real(fd), intent(IN)                :: std
    real(fd), dimension(:), intent(OUT) :: r


    real(fd), dimension(size(r)) :: v1, v2, rsq
    integer  :: i, n, ng, nn, m
    logical, dimension(size(r)) :: mask

    integer                             :: thread

    thread = 1
    !$ thread = omp_get_thread_num() + 1

#ifdef WITH_ACML
    call drandlognormal(size(r), mean, std, rnd_state(thread)%state, r, rnd_state(thread)%info)
#else
    ng = 1
    n  = size(r)

    do
       if(ng > n) exit

       do i = ng, n
          v1(i) = mtprng_rand_real1(rnd_state(thread)%state_in)
          v2(i) = mtprng_rand_real1(rnd_state(thread)%state_in)
       end do
       
       v1(ng:n)  = 2.0_fd * v1(ng:n) - 1.0_fd
       v2(ng:n)  = 2.0_fd * v2(ng:n) - 1.0_fd
 
       rsq(ng:n) = v1(ng:n)**2 + v2(ng:n)**2

       mask(ng:n) = (rsq(ng:n) > 0.0 .and. rsq(ng:n) < 1.0)

       call utl_array_copy(pack(v1(ng:n), mask(ng:n)), v1(ng:), nn, m)
       v2(ng:ng+nn-1) = pack(v2(ng:n), mask(ng:n))
       rsq(ng:ng+nn-1) = pack(rsq(ng:n), mask(ng:n))

       ng = ng+nn
    end do

    rsq = sqrt(-2.0_fd * log(rsq) / rsq)
    r = v1 * rsq

    r = exp(r*std + mean)
#endif
  end subroutine rnd_generate_logNormal



!!$  !****f* distributions/dst_generate_invGamma
!!$  ! NAME
!!$  !   dst_generate_invGamma
!!$  !
!!$  ! DESCRIPTION
!!$  !   Fills the given array with random values from the clipped inverse Gamma distribution.
!!$  !
!!$  ! INPUTS
!!$  !   a    : real(fd)    ... a parameter of the inverse gamma distribution.  
!!$  !   b    : real(fd)    ... b parameter of the inverse gamma distribution.
!!$  !   xMin : real(fd)    ... Minimum value.
!!$  !   xMax : real(fd)    ... Maximum value.
!!$  !
!!$  ! OUTPUTS
!!$  !   r    : real(fd)[]  ... Output array
!!$  !
!!$  ! SOURCE
!!$  subroutine dst_generate_invGamma(r, a, b, xMin, xMax)
!!$    real(fd), dimension(:), intent(OUT) :: r
!!$    real(fd), intent(IN) :: a, b, xMin, xMax
!!$    real(fd), dimension(5000) :: cdf, icdf
!!$
!!$    real(fd) :: dx
!!$    integer  :: i
!!$
!!$    integer                             :: thread
!!$
!!$    thread = 1
!!$    !$ thread = omp_get_thread_num()
!!$
!!$
!!$    dx = (xMax - xMin) / real(size(cdf) - 1)
!!$
!!$    do i = 1, size(cdf)
!!$       cdf(i) = gammq_s(a, b / (xMin + dx * (i-1)))
!!$    end do
!!$
!!$    call utl_invertCDF(cdf, icdf, xMin, xMax)
!!$        
!!$    do i = 1, size(r)
!!$       r(i) = mtprng_rand_real1(rnd_state(thread)%state_in)
!!$    end do
!!$
!!$    forall(i = 1 : size(r))
!!$       r(i) = utl_lerp_lin_array(icdf, r(i))
!!$    end forall
!!$
!!$  end subroutine dst_generate_invGamma
!!$  !******

end module random
