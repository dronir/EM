MODULE igdev
USE BASE
IMPLICIT NONE

CONTAINS

  real(FD) FUNCTION gammln_s(xx)
    real(FD), INTENT(IN) :: xx
    real(FD) :: tmp, x

    real(FD)  :: stp = 2.5066282746310005_fd
    real(FD), DIMENSION(6) :: coef = &
         & (/ 76.18009172947146_fd,      -86.50532032941677_fd, &
         &    24.01409824083091_fd,      -1.231739572450155_fd, &
         &      .1208650973866179e-2_fd, -.5395239384953e-5_fd /)

    real(FD), DIMENSION(6) :: y 

    y = (/1.0, 2.0, 3.0, 4.0, 5.0, 6.0/)

    x   = xx
    y   = y+x
    tmp = x + 5.5_fd
    tmp = (x+0.5_fd) * LOG(tmp) - tmp
    gammln_s = tmp + LOG(stp * ( 1.000000000190015_fd + SUM(coef / y)) / x)

  END FUNCTION gammln_s

  FUNCTION gammln_v(xx)
    real(FD), DIMENSION(:), INTENT(IN) :: xx
    real(FD), DIMENSION(SIZE(xx))      :: gammln_v
    real(FD), DIMENSION(SIZE(xx))      :: ser, tmp, x, y

    real(FD)  :: stp = 2.5066282746310005_fd
    real(FD), DIMENSION(6) :: coef = &
         & (/ 76.18009172947146_fd,      -86.50532032941677_fd, &
         &    24.01409824083091_fd,      -1.231739572450155_fd, &
         &      .1208650973866179e-2_fd, -.5395239384953e-5_fd /)

    INTEGER :: i

    x   = xx
    y   = xx
    tmp = x + 5.5_fd
    tmp = (x+0.5_fd) * LOG(tmp) - tmp

    ser =  1.000000000190015_fd

    DO i = 1, 6
       y   = y + 1.0_fd 
       ser = ser+coef(i)/y
    END DO

    gammln_v = tmp + LOG(stp * ser / x)

  END FUNCTION gammln_v

  real(FD) FUNCTION gcf_s(a,x,gln)
    real(FD), INTENT(in) :: a, x
    real(FD), OPTIONAL, INTENT(out) :: gln

    INTEGER, PARAMETER :: ITMAX = 100
    real(FD), PARAMETER :: EPS=EPSILON(x), FPMIN=TINY(x)/EPS

    INTEGER  :: i
    real(FD) :: an,b,c,d,del,h

    IF(x == 0.0_fd) THEN
       gcf_s = 1.0_fd
       RETURN
    END IF

    b = x + 1.0_fd - a
    c = 1.0_fd / FPMIN
    d = 1.0_fd / b
    h = d
    
    DO i = 1, ITMAX
       an  = -i * (i - a)
       b   = b + 2.0_fd
       d   = an * d + b
       IF(ABS(d) < FPMIN) d = FPMIN
       c   = b + an / c
       IF(ABS(c) < FPMIN) c = FPMIN
       d   = 1.0_fd / d
       del = d * c
       h   = h * del
       
       IF(ABS(del-1.0_fd) <= EPS) EXIT
    END DO
    IF(PRESENT(gln)) THEN
       gln = gammln_s(a)
       gcf_s = EXP(-x+a*LOG(x)-gln)*h
    ELSE
       gcf_s = EXP(-x+a*LOG(x)-gammln_s(a))*h
    END IF

  END FUNCTION gcf_s

  FUNCTION gcf_v(a,x,gln)
    real(FD), DIMENSION(:), INTENT(in) :: a, x
    real(FD), DIMENSION(:), OPTIONAL, INTENT(out) :: gln

    real(FD), DIMENSION(SIZE(a)) :: gcf_v

    INTEGER, PARAMETER :: ITMAX = 100
    real(FD), PARAMETER :: EPS=EPSILON(x), FPMIN=TINY(x)/EPS

    INTEGER  :: i
    real(FD), DIMENSION(SIZE(a)) :: an,b,c,d,del,h
    LOGICAL, DIMENSION(SIZE(a)) :: converged, zero

    zero = (x == 0.0_fd)

    WHERE(zero)
       gcf_v = 1.0_fd
    ELSEWHERE
       b = x + 1.0_fd - a
       c = 1.0_fd / FPMIN
       d = 1.0_fd / b
       h = d
    END WHERE

    converged = zero
    
    DO i = 1, ITMAX
       WHERE(.NOT. converged)
          an  = -i * (i - a)
          b   = b + 2.0_fd
          d   = an * d + b
          d   = MERGE(FPMIN, d, ABS(d) < FPMIN)
          c   = b + an / c
          c   = MERGE(FPMIN, c, ABS(c) < FPMIN)
          d   = 1.0_fd / d
          del = d * c
          h   = h * del

          converged = (ABS(del-1.0_fd) <= EPS)
       END WHERE
       IF(ALL(converged)) EXIT
    END DO

    IF(PRESENT(gln)) THEN
       gln = gammln_v(a)
       WHERE(.NOT. zero) gcf_v = EXP(-x+a*LOG(x)-gln)*h
    ELSE
       WHERE(.NOT. zero) gcf_v = EXP(-x+a*LOG(x)-gammln_v(a))*h
    END IF

  END FUNCTION gcf_v

  real(FD)  FUNCTION gser_s(a,x,gln)
    real(FD), INTENT(IN) :: a,x
    real(FD), OPTIONAL, INTENT(OUT) :: gln

    INTEGER, PARAMETER  :: ITMAX = 100
    real(FD), PARAMETER :: EPS = epsilon(x)

    INTEGER :: n
    real(FD) :: ap,del,summ

    IF(x == 0.0_fd) THEN
       gser_s = 0.0_fd
       RETURN
    END IF

    ap  = a
    summ = 1.0_fd / a
    del = summ
    DO n=1,ITMAX
       ap   = ap + 1.0_fd
       del  = del*x/ap
       summ = summ+del
       IF(ABS(del) < ABS(summ)*EPS) EXIT
    END DO

    IF(present(gln)) THEN
       gln = gammln_s(a)
       gser_s = summ*EXP(-x+a*LOG(x)-gln)
    ELSE
       gser_s = summ*EXP(-x+a*LOG(x)-gammln_s(a))
    end if

  END FUNCTION GSER_S


  FUNCTION gammq_s(a,x)
    real(FD), INTENT(IN) :: a, x
    real(FD) :: gammq_s

    IF(x < a+1.0_fd)THEN
       gammq_s = 1.0_fd - gser_s(a,x)
    ELSE
       gammq_s = gcf_s(a,x)
    ENDIF

  END FUNCTION GAMMQ_S


END MODULE igdev
