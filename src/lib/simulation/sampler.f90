!****h* XR/lib/simulation/sampler
! NAME
!   sampler
!
! DESCRIPTION
!
! NOTES
!
! AUTHOR
!   Hannu Parviainen
!
! USES
!   Distributions
!
!
! CREATION DATE
!   21.12.2007
!******
module sampler
  use random 
  implicit none

contains

   subroutine smpl_griddedSamples2D(samples, nSamples)
     real(fd), dimension(:,:)      :: samples
     integer                       :: nSamples
     
     integer  :: res, i, j
     real(FD) :: dl

     if(size(samples,1) /= 2) then
        call utl_fatal_error("Size of the sample array dim(1) is wrong, should be 2.")
     end if

     res = floor(sqrt(real(nSamples, fd)))

     if(size(samples,2) < res**2) then
        call utl_fatal_error("Not enough space in the given sample array.")
     end if

     dl = 1.0_fd / res

     call rnd_generate_uniform_n(samples(1,:))
     call rnd_generate_uniform_n(samples(2,:))

     forall(i = 1:res, j = 1:res)
        samples(1, i + res*(j-1)) =  (samples(1, i + res*(j-1)) + i - 1) * dl
        samples(2, i + res*(j-1)) =  (samples(2, i + res*(j-1)) + j - 1) * dl
     end forall

   end subroutine smpl_griddedSamples2D

end module sampler
