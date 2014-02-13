!****h* XR/lib/core/brdf
! NAME
!   brdf
!
! DESCRIPTION
!   
!
! AUTHOR
!   Hannu Parviainen
!
! CREATION DATE
!   15.01.2007
!******module brdf
module brdf
  use base

contains
  
  real(fd) function brdf_shadowing(L, I, N, w, Pf, Pp) result(c)
    real(fd), dimension(3)       :: L, I, N
    real(fd)                     :: w
    real(fd), optional, external :: Pf    
    real(fd), optional           :: Pp(:)

    c = 1.0_fd
    
  end function brdf_shadowing
  

  real(fd) function brdf_Lambert(L, I, N, w, Pf, Pp) result(c)
    real(fd), dimension(3)       :: L, I, N
    real(fd)                     :: w
    real(fd), optional, external :: Pf    
    real(fd), optional           :: Pp(:)

    c = w * INV_PI * dot_product(L, N)
 
    if(c < 0.0_fd) c = 0.0_fd

  end function brdf_Lambert


  real(fd) function brdf_LommelSeeliger(L, I, N, w, Pf, Pp) result(c)
    real(fd), dimension(3)       :: L, I, N
    real(fd)                     :: w

    real(fd), optional, external :: Pf       !! Phase function
    real(fd), optional           :: Pp(10)   !! Phase function parameters

    real(fd) :: mu, mu0, a, Pv

    mu  = dot_product(-I, N)
    mu0 = dot_product( L, N)

    !! Compute the value of the phase function, if present.
    !!
    Pv = 1.0
    if(present(Pf) .AND. present(Pp)) then
       a  = dot_product(-I,L) !! a = cos(phase)
       Pv = Pf(a, Pp)
    end if

    !! Compute the Lommel-Seeliger brdf
    !!
    if(mu > 0.0_fd .and. mu0 > 0.0_fd) then
       c = w * INV_FOUR_PI * mu0 / (mu + mu0) * Pv
    else
       c = 0.0_fd
    end if

  end function brdf_LommelSeeliger


  real(fd) function brdf_RadTrans(L, I, N, w, Pf, Pp) result(c)
    real(fd), dimension(3)       :: L, I, N
    real(fd)                     :: w 
    real(fd), optional, external :: Pf    
    real(fd), optional           :: Pp(:)

    real(fd) :: mu, mu0, H_mu, H_mu0, a, Pv

    mu  = dot_product(-I, N)
    mu0 = dot_product( L, N)

    Pv = 1.0
    if(present(Pf) .AND. present(Pp)) then
       a  = abs(dot_product(-I,L))
       Pv = Pf(a, Pp)
    end if

    if(mu > 0.0_fd .and. mu0 > 0.0_fd) then

       H_mu0 = H(mu0, w)
       H_mu  = H(mu,  w)

       c = w * INV_FOUR_PI * mu0 / (mu + mu0) * H_mu0 * H_mu * Pv
    else
       c = 0.0_fd
    end if
    
    contains
      real(fd) function H(x, w)
        real(fd) :: x, w
        real(fd) :: g, r0
    
        g  = sqrt(1.0_fd - w)
        r0 = (1.0-g)/(1.0+g)
        H = 1.0 / (1.0 - w*x*(r0 + 0.5*(1.0 - 2.0*r0*x) * log((1.0+x)/x)))

      end function H

  end function brdf_RadTrans

  pure real(fd) function phase_function_constant(a, par) result(P)
    real(fd), intent(in) :: a, par(:)
    P = 1.0_fd
  end function phase_function_constant

  pure real(fd) function phase_function_HG1(a, par) result(P)
    real(fd), intent(in) :: a, par(1)
    P = (1.0 - par(1)**2) / (1.0 + par(1)**2 + 2.0*par(1)*a)**1.5
  end function phase_function_HG1

  pure real(fd) function phase_function_HG2(a, par) result(P)
    real(fd), intent(in) :: a, par(3)
    real(fd) :: w, g1, g2
    w  = par(1)
    g1 = par(2)
    g2 = par(3)

!    P = w  * (1.0_fd - g1**2) / (1.0_fd + g1**2 + 2.0_fd*g1*cos(PI - a))**1.5_fd + &
!         & (1.0_fd - w) * (1.0_fd - g2**2) / (1.0_fd + g2**2 + 2.0_fd*g2*cos(PI - a))**1.5_fd

    P = w  * (1.0_fd - g1**2) / (1.0_fd + g1**2 + 2.0_fd*g1*a)**1.5_fd + &
         & (1.0_fd - w) * (1.0_fd - g2**2) / (1.0_fd + g2**2 + 2.0_fd*g2*a)**1.5_fd


  end function phase_function_HG2

end module brdf
