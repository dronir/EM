module rfield

  use base
  use util
  use mtprng

  implicit none
  include 'fftw3.f'

  type :: rf_randomField
     integer      :: res
     real(FD)     :: width

     integer*8    :: fftw_plan

     complex(FD), dimension(:,:), pointer :: spectrum
     complex(FD), dimension(:,:), pointer :: field
  end type rf_randomField
  
  integer, parameter :: RF_SPECTRUM_GAUSSIAN = 0
  integer, parameter :: RF_SPECTRUM_FBM      = 1

contains

  elemental real function power_fbm(k, H)
    real(FD), intent(IN) :: k, H
    power_fbm = k ** ( -(H + 1.0) )
  end function power_fbm

  elemental real function power_fbm_sqr(k_sqr, H)
    real(FD), intent(IN) :: k_sqr, H
    power_fbm_sqr = k_sqr ** ( -(H + 1.0) * 0.5 )
  end function power_fbm_sqr

  elemental real function power_Gaussian_sqr(k_sqr, l)
    real(FD), intent(IN) :: k_sqr, l
    power_Gaussian_sqr = l/(2.0*sqrt(pi)) * exp(-k_sqr * l**2 * 0.25)
  end function power_Gaussian_sqr


  subroutine rf_init(f, r, w)
    type(rf_randomField) :: f
    integer :: r
    real(FD) :: w

    f % res = r
    f % width = w

    allocate(f % spectrum (0:r-1, 0:r-1))
    allocate(f % field    (0:r-1, 0:r-1))

    f % spectrum = (0.,0.)
    f % field    = (0.,0.)

    call dfftw_plan_dft_2d(f%fftw_plan, f%res, f%res, f%field, f%field, FFTW_BACKWARD, FFTW_ESTIMATE)

  end subroutine rf_init


  subroutine rf_generateSpectrum(f, spectrumType, p1, p2)
    type(rf_randomField) :: f
    integer              :: spectrumType
    real(fd)             :: p1, p2

    integer :: i, j, prty, hres
    real(FD) :: power, dk, std_factor
    real(FD), dimension(:,:), allocatable :: rad
    complex(FD), dimension(:,:), allocatable :: tSpectrum

    allocate(rad(0:f%res-1, 0:f%res-1))

    dk = TWO_PI / f % width

    hres  = f%res / 2
    power = 0.0_fd
    prty  = 0
    rad   = 0.0_fd

    if(mod(f%res,2) == 0) prty = 1

    select case(spectrumType)
    case(RF_SPECTRUM_GAUSSIAN)

       !$omp parallel workshare
       forall( i = 0 : hres-1, j = 0 : hres-1, i /= 0 .and. j /= 0 )
          rad(i,j) = power_Gaussian_sqr((i*dk)**2 + (j*dk)**2, p1)
       end forall
       forall( i = 1 : hres-1, j = hres+1 : f%res-1 )
          rad(i,j) = power_Gaussian_sqr((i*dk)**2 + ((f%res - j)*dk)**2, p1)
       end forall
       !$omp end parallel workshare
       
    case(RF_SPECTRUM_FBM)
       
       !$omp parallel workshare
       forall( i = 0 : hres-1, j = 0 : hres-1, i /= 0 .and. j /= 0 )
          rad(i,j) = power_fbm_sqr((i*dk)**2 + (j*dk)**2, p1)
       end forall
       forall( i = 1 : hres-1, j = hres+1 : f%res-1 )
          rad(i,j) = power_fbm_sqr((i*dk)**2 + ((f%res - j)*dk)**2, p1) 
       end forall
       !$omp end parallel workshare
       
    case default
       write(*,'("Error: Unsupported power spectrum type. Exiting.")')
       stop

    end select

    !$omp parallel workshare
    f % spectrum = cmplx(rad,rad)

    f % spectrum(hres+1 : ,      0          ) = conjg(f%spectrum( hres-1 : 0 : -1,          0     ))
    f % spectrum(     0   , hres+1 :        ) = conjg(f%spectrum(          0     , hres-1 : 1 :-1 ))

    f % spectrum(hres+1 : , hres+1 :        ) = conjg(f%spectrum( hres-1 : 1 : -1, hres-1  : 1    : -1 ))
    f % spectrum(hres+1 : ,      1 : hres-1 ) = conjg(f%spectrum( hres-1 : 1 : -1, f%res-1 : hres : -1 ))

    power      = 4.0_fd * sum(rad**2)
    std_factor = p2     / sqrt(power)

    f % spectrum = f % spectrum * std_factor
    !$omp end parallel workshare

    deallocate(rad)

  end subroutine rf_generateSpectrum

  subroutine rf_generateField(f, seed)
    type(rf_randomField) :: f
    integer            :: seed

    type(mtprng_state) :: mtp_state
    real(FD) :: phase, rad, rnd(2)
    integer   :: i, i0, j, j0

    f%field = f%spectrum
    
    call mtprng_init(seed, mtp_state) 

    do j = 0, f % res / 2 - 1
         do i = 0, f % res / 2 - 1

          rnd(1) = mtprng_rand_real1(mtp_state)
          rnd(2) = mtprng_rand_real1(mtp_state)

          phase = 2.0 * pi * RND(1)
          rad   = sqrt(-2.0 * log(RND(2)))

          f%field(i,j) = cmplx(real(f%field(i,j)) * rad * cos(phase), aimag(f%field(i,j)) * rad * sin(phase))
          
          if(i==0) then
             i0 = 0
          else
             i0 = f % res - i
          end if

          if(j==0) then
             j0 = 0
          else
             j0 = f % res - j
          end if

          f%field(i0,j0) = conjg(f%field(i,j))

       end do
    end do

    do j = f % res / 2 + 1, f % res - 1
       do i = 1, f % res / 2 - 1

          rnd(1) = mtprng_rand_real1(mtp_state)
          rnd(2) = mtprng_rand_real1(mtp_state)

          phase = 2.0 * pi * RND(1)
          rad   = sqrt(-2.0 * log(RND(2)))

          f%field(i,j) = cmplx(real(f%field(i,j)) * rad * cos(phase), aimag(f%field(i,j)) * rad * sin(phase))

          i0 = f % res - i
          j0 = f % res - j

          f % field(i0,j0) = conjg(f%field(i,j))
       end do
    end do

    call dfftw_execute(f%fftw_plan)

  end subroutine rf_generateField

  !!
  !! For the full definition of different quantities computed here, see 
  !!   Shepard, M., K., Icarus 141, 156-171 (1999)
  !!   Shepard, M., K., Journal of Geophysical Research, vol. 106, no. E12, pages 32,777-32,795, 2001
  !!
  subroutine rf_computeStatistics(f)
    type(rf_randomField) :: f

    real(fd) :: mean, var, std
 
    real(fd) :: xi_0 = 0.0_fd  ! Unit RMS height
    real(fd) :: mu_0 = 0.0_fd  ! Unit RMS deviation
    real(fd) :: s_0  = 0.0_fd  ! Unit RMS slope

    real(fd) :: rx(100), ry(100)
    integer  :: x(100), y(100)
    integer  :: i
 
    real(fd) :: a, b, c

    !$omp parallel workshare
    mean = sum(real(f%field,fd)) / real(f%res**2, fd)
    var  = sum((real(f%field,fd) - mean)**2) /  (real(f%res**2, fd) - 1.0_fd)
    std  = sqrt(var)
    !$omp end parallel workshare

    call random_number(rx)
    call random_number(ry)

    x = floor(rx * real(f%res- 10,fd)) + 5
    y = floor(ry * real(f%res- 10,fd)) + 5

    ! note to self, this is actually mu_0...
    do i = 1, 100
       a = real(f%field(  x(i), y(i)), fd)
       b = real(f%field(x(i)+1, y(i)), fd)
       c = 0.5_fd * (a + b)
       xi_0 = xi_0 +  0.5_fd * ((a - c)**2 + (b - c)**2)
    end do

    xi_0 = sqrt(xi_0 / 100.0_fd)

    print '(3f8.2)', mean, sqrt(var), xi_0

  end subroutine rf_computeStatistics

  subroutine rf_writeFits(f, fName)
    type(rf_randomField) :: f
    character(len=*)   :: fName

    integer :: status,unit, blocksize,bitpix,naxis,naxes(2), nelements, group,fpixel
    logical :: simple, extend

    real(FD), dimension(f%res,f%res) :: data

#ifdef WITH_CFITSIO
    status = 0

    simple=.true.
    bitpix=-32
    naxis=2
    naxes(1)=f%res
    naxes(2)=f%res
    extend=.true.

    blocksize = 1
    group=1
    fpixel=1
    nelements=naxes(1)*naxes(2)

    data = real(f % field, FD)

    call ftgiou(unit,status)
    call ftinit(unit,fName,blocksize,status)

    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    call ftpprd(unit,group,fpixel,nelements,data,status)

    call ftpkyd(unit,'WIDTH',f%width,3,'Field width',status)

    call ftclos(unit, status)
    call ftfiou(unit, status)
#else
     call utl_message("Warning: Cannot write fits-file, not compiled with fits-support.")
#endif
  end subroutine rf_writeFits

end module rfield
