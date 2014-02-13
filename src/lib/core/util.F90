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

!> Module containing miscellanceous utility functions.
!!
!!
!! \date 16.11.2007
!!
module util
  use base
  
  implicit none

  !> Cpu timer
  !! 
  type :: utl_timer
     real(fd)           :: dc, dcc, state, msgPeriod, tStart
     character(len=100) :: name
  end type utl_timer

  !> The intent level of text output.
  !!
  integer               :: utl_output_indent = 0

contains

  !> Converts time given in seconds into a string "h hours m min s s."
  !! 
  character(len=22) function utl_secs_to_hms(s) result(t)
    real(fd) :: s
    real(fd) :: tHours, tMins, tSecs

    tHours = floor(s / 3600.0)
    tMins  = floor((s - tHours * 3600.0) / 60.0)
    tSecs  = s - tHours * 3600.0 - tMins * 60.0

    write(t, '((I4)," h ",(I2)," min ",(F6.3)," s")') int(tHours), int(tMins), tSecs

  end function utl_secs_to_hms


  !!--- TIMER ROUTINES ---
  !!
  !!

  !> Initialises cpu timer.
  !!
  subroutine utl_timer_init(c, msgPeriod, nMaxSteps)
    type(utl_timer) :: c
    integer         :: nMaxSteps
    real(fd)        :: msgPeriod, tStart

    c%dc        = 100.0_fd / real(nMaxSteps-1, fd)
    c%dcc       = msgPeriod 
    c%msgPeriod = msgPeriod
    c%state     = 0.0_fd

    call cpu_time(c%tStart)

  end subroutine utl_timer_init

  !> Increases cpu timer
  !!
  subroutine utl_timerIncrease(c)
    type(utl_timer) :: c

    real(fd) :: tCurr, tElap

    !$omp critical
    c%state = c%state + c%dc

    if(c%state >= c%dcc) then
       call cpu_time(tCurr)
       tElap = tCurr - c%tStart

       write(*,'("   " F6.1 "    " A22 "    " A22)') c%state, utl_secs_to_hms(tElap), utl_secs_to_hms(tElap * (100.0 / c%state))

       c%dcc = c%dcc + c%msgPeriod
    end if
    !$omp end critical

  end subroutine utl_timerIncrease


  !!--- MESSAGE ROUTINES ---
  !!
  !!
  subroutine utl_increase_indent(i)
    integer, optional :: i
    integer           :: ii

    ii = 1
    if(present(i)) ii = i
    !$omp atomic
    utl_output_indent = utl_output_indent + ii
  end subroutine utl_increase_indent


  subroutine utl_decrease_indent(i)
    integer, optional :: i
    integer           :: ii

    ii = 1
    if(present(i)) ii = i
    !$omp atomic
    utl_output_indent = utl_output_indent - ii
    if(utl_output_indent < 0) utl_output_indent = 0
  end subroutine utl_decrease_indent


  subroutine utl_message(msg)
    character(len=*)  :: msg
    character(len = utl_output_indent) :: indent

    indent = ''

    write(*,'(A)') indent // "Message: " // msg
    
  end subroutine utl_message


  subroutine utl_warning(msg)
    character(len=*) :: msg

    write(*,'("Warning: ",(A))') msg
  end subroutine utl_warning


  subroutine utl_fatal_error(msg)
    character(len=*) :: msg

    write(*,'("Fatal error: ",(A))') msg
    stop

  end subroutine utl_fatal_error


  !!--- ARRAY ROUTINES ---
  !!
  !!

  !****f* Util/utl_bisectSearch
  ! NAME
  !   utl_bisectSearch
  !
  ! DESCRIPTION
  !  
  !   Basic bisection search. 
  !   "Numerical Recipes", 3rd edition, 2007, p. 115
  !   
  ! SOURCE
  !
  pure integer function utl_bisectSearch(x, d) result(il)
    real(FD), dimension(:), intent(IN) :: d
    real(fd), intent(IN)               :: x

    integer  :: iu, im
    logical  :: ascending
    
    ascending = (d(2) >= d(1))

    il = 1
    iu = size(d)

    do while (iu - il > 1)
       im = (iu + il) / 2
       if((x > d(im)) .eqv. ascending) then
          il = im
       else  
          iu = im
       end if
    end do

  end function utl_bisectSearch


  !****f* Util/utl_lerp_free_array
  ! NAME
  !   utl_lerp_free_array
  !
  ! DESCRIPTION
  !   Evaluates a tabulated function v(x).
  !   
  ! SOURCE
  !
  pure real(FD) function utl_lerp_free_array(x, xArr, yArr) result(f)
    real(FD), dimension(:), intent(IN) :: xArr, yArr
    real(FD), intent(IN) :: x
    
    integer  :: i

    i = utl_bisectSearch(x, xArr)
    f = yArr(i) + ((x - xArr(i)) / (xArr(i+1) - xArr(i))) * (yArr(i+1) - yArr(i))
 
  end function utl_lerp_free_array


  !****f* Util/utl_lerp_lin_array
  ! NAME
  !   utl_lerp_lin_array
  !
  ! DESCRIPTION
  !   Evaluates a tabulated function v(x), where x is in range [0,1].
  !   
  ! SOURCE
  !
  pure real(FD) function utl_lerp_lin_array(data, x)
    real(FD), dimension(:), intent(IN) :: data
    real(FD), intent(IN) :: x
    
    real(FD) :: dtx
    integer  :: dti
  
    if(x <  0.0_fd) then
       utl_lerp_lin_array = data(1)
       return
    end if
    
    if(x >= 1.0_fd) then
       utl_lerp_lin_array = data(size(data))
       return
    end if

    dtx = x * real(size(data)-1, FD) + 1.0_fd
    dti = floor(dtx)

    if(dti < size(data)) then
       utl_lerp_lin_array = data(dti) * (1.0_fd - dtx + dti) + data(dti+1) * (dtx - dti)
    else
       utl_lerp_lin_array = data(dti)
    end if

  end function utl_lerp_lin_array
  !******

  !****f* Util/utl_lerp_log_array
  ! NAME
  !   utl_lerp_log_array
  !
  ! DESCRIPTION
  !   Evaluates a tabulated function v(x), where the x-axis is
  !   logarithmic, and also may contain additional data points.
  !
  ! SOURCE
  !
  real(FD) function utl_lerp_log_array(xTable, vTable, x)
    real(FD), intent(IN), dimension(:) :: xTable, vTable
    real(FD), intent(IN)               :: x

    integer       :: nX, i
    real(FD)      :: dX, a

    nX = size(xTable)

    if(x <= xTable(1)) then
       utl_lerp_log_array = vTable(1)
       return
    end if

    if(x >= xTable(nX)) then
       utl_lerp_log_array = vTable(nX)
       return
    end if

    dX = (log10(xTable(nX)) - log10(xTable(1))) / real(nX - 1, KIND=FD)
    i  = aint((log10(x) - log10(xTable(1))) / dX) + 1

    if(x < xTable(i)) then
       do while(x < xTable(i))
          i = i - 1
       end do
    else if(x > xTable(i+1)) then
       do while(x > xTable(i+1))
          i = i + 1
       end do
    end if

    a = (x - xTable(i)) / (xTable(i+1) - xTable(i))
   
    utl_lerp_log_array = (1.0_fd - a) * vTable(i) + a * vTable(i+1) 

  end function utl_lerp_log_array
  !******

  !****f* Util/utl_array_copy
  ! NAME
  !   utl_array_copy
  !
  ! DESCRIPTION
  !   
  !   
  !
  ! SOURCE
  !
  subroutine utl_array_copy(src, dest, n_copied, n_not_copied)
    real(FD), dimension(:), intent(IN)  :: src
    real(FD), dimension(:), intent(OUT) :: dest
    integer, intent(OUT)                :: n_copied, n_not_copied

    n_copied         = min(size(src), size(dest))
    n_not_copied     = size(src) - n_copied
    
    dest(1:n_copied) = src(1:n_copied)
  end subroutine utl_array_copy
 !******

  !****f* Util/utl_invert_cdf
  ! NAME
  !   utl_invert_cdf
  !
  ! DESCRIPTION
  !   
  !   
  !
  ! SOURCE
  !
  subroutine utl_invertCDF(cdf, icdf, xMin, xMax)
    real(FD), dimension(:), intent(IN)  ::  cdf
    real(FD), dimension(:), intent(OUT) :: icdf

    real(FD), intent(IN) :: xMin, xMax

    real(FD), dimension(size(cdf)) :: tcdf

    real(FD) :: y, dx, dy, idy, yMin, yMax

    integer  :: i, ix_a, ix_b

    dx = (xMax - xMin) / real(size(cdf) - 1)
    dy = 1.0_fd / real(size(icdf) - 1)

    yMin = minval(cdf)
    yMax = maxval(cdf)

    tcdf = (cdf - yMin) * 1.0_fd / (yMax - yMin)

    icdf(1) = xMin
    icdf(size(icdf)) = xMax

    ix_a = 1
    do i = 2, size(icdf) - 1
       y = yMin + dy * (i-1)
       ix_b = size(tcdf)

       do while (ix_b - ix_a > 1)
          if(tcdf((ix_a + ix_b) / 2) > y) then
             ix_b = (ix_a + ix_b) / 2
          else
             ix_a = (ix_a + ix_b) / 2
          end if
       end do

       ! Solve x from y = a(1-x) + bx
       !
       icdf(i) = xMin + (ix_a - 1)*dx + (y - tcdf(ix_a)) / (tcdf(ix_b) - tcdf(ix_a))

    end do

  end subroutine utl_invertCDF
   

  !****f* Util/utl_invert_cdf
  ! NAME
  !   utl_invert_cdf
  !
  ! DESCRIPTION
  !   
  !   
  !
  ! SOURCE
  !
  subroutine utl_invert_cdf(xArr, yArr, iArr)
    real(FD), dimension(:), intent(IN)  :: xArr, yArr
    real(FD), dimension(:), intent(OUT) :: iArr

    real(FD) :: y, dy, yMin, yMax
    integer  :: i

    yMin = yArr(1)
    yMax = yArr(size(yArr))

    dy   = (yMax - yMin) / real(size(iArr) - 1)

    do i = 1, size(iArr)
       y       = yMin + dy * (i-1)
       iArr(i) = utl_lerp_free_array(y, yArr, xArr)
    end do

   end subroutine utl_invert_cdf



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! FITS HANDLING ROUTINES
  !!
  !!

   subroutine utl_writeFits2D(data, fName)
     real(FD), dimension(:,:)    :: data
     character(len=*)            :: fName

     integer :: status,unit, blocksize,bitpix,naxis,naxes(2), nelements, group,fpixel
     logical :: simple, extend

#ifdef WITH_CFITSIO
     status = 0

     simple=.true.
     bitpix=-32
     naxis=2
     naxes(1)=size(data,1)
     naxes(2)=size(data,2)
     extend=.true.

     blocksize = 1
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)

     call ftgiou(unit,status)
     call ftinit(unit,fName,blocksize,status)
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
     call ftpprd(unit,group,fpixel,nelements,data,status)

     call ftclos(unit, status)
     call ftfiou(unit, status)
#else
     call utl_message("Warning: Cannot write fits-file, not compiled with fits-support.")
#endif
   end subroutine utl_writeFits2D

   subroutine utl_writeFits3D(data, fName)
     real(FD), dimension(:,:,:)  :: data
     character(len=*)            :: fName

     character(len=80) :: e

     integer :: status,unit, blocksize,bitpix,naxis,naxes(3), nelements, group,fpixel
     logical :: simple, extend

#ifdef WITH_CFITSIO
     status = 0

     simple=.true.
     bitpix=-32
     naxis=3
     naxes(1)=size(data,1)
     naxes(2)=size(data,2)
     naxes(3)=size(data,3)
     extend=.true.

     blocksize = 1
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)*naxes(3)

     call ftgiou(unit,status)
     call ftinit(unit,fName,blocksize,status)
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
     call ftpprd(unit,group,fpixel,nelements,data,status)

     call ftclos(unit, status)
     call ftfiou(unit, status)
#else
     call utl_message("Warning: Cannot write fits-file, not compiled with fits-support.")
#endif
   end subroutine utl_writeFits3D

 end module util
