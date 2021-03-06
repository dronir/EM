!!--- BEGIN GPL --- 
!!
!! xrfpm -  a program to compute soft X-ray fluorescence of particulate media.
!! Copyright (C) 2009 Hannu Parviainen
!!
!! This program is part of UHOEM.
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! Contributor(s): Hannu Parviainen
!!
!!--- END GPL ---

!> Computes light scatterig from particulate media.
!!
!! \author Hannu Parviainen
!! \version 1.0
!! \date 15.11.2007

program fit_2hg
  use iso_c_binding
  use netcdf
  use base
  use util
  use geometry
  use particle
  use hemisphere
  use container_grid3d
  use random
  use medium
  use trace
  use sampler
  use rfield
  use brdf

  !$ use omp_lib

  implicit none 

  character (len=FNAME_LENGTH) :: parFilename        = ""
  character (len=FNAME_LENGTH) :: outFilename        = "dscatter.scangles"
  character (len=FNAME_LENGTH) :: mediumFileName     = "medium.nc"
  character (len=100)          :: sDistribution      = "constant"
  character (len=100)          :: brdfType           = "shadowing"
  character (len=100)          :: brdfPhaseFunction  = "constant"

  integer                      :: gridResHorizontal  = 200
  integer                      :: gridResVertical    = 10

  real(fd)                     :: rho                = 0.50_fd
  real(fd)                     :: rhoAllowedError    = 0.05_fd

  integer                      :: resThetaE          = 45
  integer                      :: sampleSeed         = 0


  integer                      :: nOrders            = 1
  character(len=100)           :: nSamplesPerOrder   = "1"
  integer, allocatable         :: nSamplesPerOrderTable(:)

  integer                      :: nThetaIn           = 1
  character (len=10000)        :: thetaIn            = "45.0"
  real(fd), allocatable        :: thetaInTable(:)

  integer                      :: nPfParams          = 1
  character (len=10000)        :: pfParams           = "0.0"
  real(fd)                     :: pfParamsTable(10)

  logical                      :: full_output       = .true.

  logical                      :: rf_applyFields    = .false.
  integer                      :: rf_nFieldsPerMed  = 1
  character(len=100)           :: rf_spectrumType   = "Gaussian"
  real(fd)                     :: rf_P              = 0.5_fd
  real(fd)                     :: rf_std            = 1.0_fd

  integer                      :: mediumMapRes        = 256
  integer                      :: mediumDensitymapRes = 300

  integer                      :: nThreads            = 1

  type(med_medium)      :: M
  type(med_mediumFile)  :: Mf
  type(rf_randomField)  :: rndField
  type(utl_timer)       :: timer

  real    :: sTime, cTime, eTime
  integer :: i, iField

  character (len=FNAME_LENGTH) :: photFilename = "photometry.nc"

  real(fd), dimension(:), allocatable :: iTheta, eTheta, ePhi, intensity
  real(fd), dimension(:,:), allocatable :: mu, mut, mu0, mu0t, cosa, cosat
  integer, dimension(:), allocatable :: nIlluminated, nIlluminatedt
  integer :: n_datapoints, nSurfaceSamples, ipta, iptb
  
  namelist /params/ &
       & resThetaE, gridResHorizontal, gridResVertical, &
       & outFilename, mediumFilename, photFilename, &
       & sDistribution, rho, rhoAllowedError, sampleSeed, &
       & nThetaIn, thetaIn, &
       & nOrders, nSamplesPerOrder, &
       & nPfParams, pfParams, &
       & brdfType, brdfPhaseFunction, &
       & rf_applyFields, rf_nFieldsPerMed, rf_spectrumType, rf_P, rf_std, &
       & nThreads, &
       & mediumMapRes, mediumDensitymapRes, full_output


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! INITIALIZE SIMULATION PARAMETERS
  !!

  if(command_argument_count() == 0) then
     write(*,'("No parameter file given, using default values.")')
  else if(command_argument_count() == 1) then
     call get_command_argument(1, parFilename)
     write(*,'(("Using parameter file "),(A))') parFilename
  else
     write(*,'("Usage: fit_2hg paramFile")')
     stop
  end if

  !! Read the namelist from the parameter file.
  !!
  if(parFilename /= "") then
     open(1,file=parFilename,status="old", form="formatted")
     read(1,NML=params)
     close(1)
  end if

  !!- Initialize random number generator
  !!
  call rnd_init(sampleSeed, nThreads, .true.)


  !! Read the number of samples for each scattering order into a table.
  !!
  allocate(nSamplesPerOrderTable(nOrders))
  read(nSamplesPerOrder,*) nSamplesPerOrderTable
  nSurfaceSamples = nSamplesPerOrderTable(1)

  !! Select the media for the simulation from a file.
  !!
  call med_mediumFileOpen(Mf, mediumFilename)
  call med_mediumFileSelectMedia(Mf)

  if(rf_applyFields) then
     call rf_init(rndField, 1000, 5.0_fd)
     select case(rf_spectrumType)        
     case('Gaussian')
        call rf_generateSpectrum (rndField, RF_SPECTRUM_GAUSSIAN, rf_P, rf_std)
     case('fBm')
        call rf_generateSpectrum (rndField, RF_SPECTRUM_FBM, rf_P, rf_std)
     end select
  else
     rf_nFieldsPerMed = 1
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! READ PHOTOMETIC DATA
  !!
  call load_datapoints(intensity, iTheta, eTheta, ePhi, n_datapoints, photFilename)
  allocate(nIlluminated(n_datapoints), nIlluminatedt(n_datapoints))

  allocate(mut(nSurfaceSamples, n_datapoints))
  allocate(mu0t(nSurfaceSamples, n_datapoints))
  allocate(cosat(nSurfaceSamples, n_datapoints))

  allocate(mu(nSurfaceSamples * rf_nFieldsPerMed, n_datapoints))
  allocate(mu0(nSurfaceSamples * rf_nFieldsPerMed, n_datapoints))
  allocate(cosa(nSurfaceSamples * rf_nFieldsPerMed, n_datapoints))

  nIlluminated = 0
  mu = -1.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! RUN SIMULATION
  !!

  call utl_timer_init(timer, 5.0_fd, Mf%nSelectedMedia * rf_nFieldsPerMed * n_datapoints)

  do i = 1, MF%nSelectedMedia
     do iField = 1, rf_nFieldsPerMed
        call med_mediumFileRead(M, Mf, Mf%varSelection(i))
        call med_gridAssign(M, gridResHorizontal, gridResVertical)
        call med_gridFit(M)
        call med_updateStatistics(M)

        if(rf_applyFields) then
           !!- Smooth the medium surface
           !!
           if(rf_spectrumType == 'constant') then
              rndField%field = rf_P * M%height

           !!- Apply macroscale surface roughness.
           !!
           else
              call rf_generateField(rndField, iField)
              call med_maskHeight(m, real(rndField%field, fs))
           end if
        end if

        nIlluminatedt = 0
        mut = -1.
        mu0t = 0.
        cosat = 0.

        call cpu_time(sTime) 
        call sampleDiscrete(nSurfaceSamples, iTheta, eTheta, ePhi, nIlluminatedt, mut, mu0t, cosat)
        call cpu_time(cTime)

        ipta = (iField-1)*nSurfaceSamples + 1
        iptb = (iField)*nSurfaceSamples + 1
        mu(ipta:iptb,:) = mut
        mu0(ipta:iptb,:) = mu0t
        cosa(ipta:iptb,:) = cosat
        nIlluminated = nIlluminated + nIlluminatedt
     end do
  end do
  write(*,'("Time taken: " ,(F8.3))') cTime - sTime

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! WRITE FILES AND CLEAN UP
  !!
  call save(outFilename, nSurfaceSamples * rf_nFieldsPerMed, n_datapoints, nIlluminated, mu, mu0, cosa)
  call med_mediumFileClose(Mf)
  

contains
  subroutine sampleDiscrete(nSurfaceSamples, iTheta, eTheta, ePhi, nIlluminated, mu, mu0, cosa)
    integer, intent(in)               :: nSurfaceSamples
    real(fd), dimension(:), intent(in):: iTheta, eTheta, ePhi

    integer, dimension(size(iTheta)), intent(out) :: nIlluminated
    real(fd), dimension(nSurfaceSamples,size(iTheta)), intent(out) :: mu, mu0, cosa

    type(ray) :: rC
    type(ray) :: rS
    type(intersection_geometry) :: iSect

    real(fd), dimension(2,nSurfaceSamples) :: samples
    real(fd), dimension(3)                 :: pSurface
    real(fd)                               :: dz, tmu0, L(3)

    integer  :: nDatapoints, idp, iss
    logical  :: pFound, pLit

    nDatapoints = size(iTheta)

    call smpl_griddedSamples2D(samples, nSurfaceSamples)

    samples = (samples * M%width - M%hWidth) * 0.5_fd
    dz      = M%grid%height - M%hMean - TRACE_EPS

    do idp = 1, nDatapoints
       !!--- Initialize the lightsource direction ---
       !!
       !!             |  A
       !!             | /
       !!             |/
       !!       --------------
       !! 
       L    = [sin(iTheta(idp)), 0.0_fd, cos(iTheta(idp))]

       !!--- Initialize the observer ray ---
       !!
       !!             |  /
       !!             | /
       !!             |v
       !!       --------------
       !! 
       call ray_init(rC, RAY_TYPE_CAMERA)
       rC%D = [-cos(ePhi(idp))*sin(eTheta(idp)), -sin(ePhi(idp))*sin(eTheta(idp)), -cos(eTheta(idp))]
       rC%rayID = rC%rayID + 100

       if(rC%D(3) > -1e-2_fd) then
          rC%D(3) = -1e-2
          call vec_normalize(rC%D)
       end if

       do iss = 1, nSurfaceSamples
          pSurface(1:2) = samples(:,iss)
          pFound = .false.
          do while (.not. pFound)
             call rayToBBTop(rC, M, pSurface, dz)
             pFound = trc_traceNearest(M%grid, rC, iSect)
             if(pFound) then
                tmu0 = dot_product(L, iSect%N)

                if(tmu0 > 0.) then
                   call ray_init(rS, RAY_TYPE_SHADOW)
                   rS%D      = L
                   rS%P      = iSect%P1 + TRACE_EPS * iSect%N
                   pLit      = .not. trc_traceOcclusion(M%grid, rS)

                   if(pLit) then
                      nIlluminated(idp) = nIlluminated(idp) + 1
                      mu0(nIlluminated(idp),idp)  = tmu0
                      mu(nIlluminated(idp),idp)   = dot_product(-rC%D, iSect%N)
                      cosa(nIlluminated(idp),idp) = dot_product(-rC%D, rS%D)
                   end if
                end if
             else
                call rnd_generate_uniform(0.0_fd, 1.0_fd, pSurface(1:2))
                pSurface(1:2) = (pSurface(1:2) * M%width - M%hWidth) * 0.5_fd   
             end if
          end do
       end do
       call utl_timerIncrease(timer)
    end do
  end subroutine sampleDiscrete

  subroutine rayToBBTop(r, M, pSurface, dz)
    type(ray)        :: r
    type(med_medium) :: M
    real(fd)         :: dz, pSurface(3)
    
    r % P(1)     = pSurface(1) + dz * (r%D(1) / r%D(3))
    r % P(2)     = pSurface(2) + dz * (r%D(2) / r%D(3))
    r % P(3)     = M%grid%height - TRACE_EPS
    r % P(1)     = modulo(r%P(1)+M%hWidth, M%width) - M%hWidth
    r % P(2)     = modulo(r%P(2)+M%hWidth, M%width) - M%hWidth
    r % rayID    = r % rayID + RAY_ID_INCR
  end subroutine rayToBBTop

  subroutine load_datapoints(intensity, iTheta, eTheta, ePhi, n_datapoints, photFilename)
    real(fd), dimension(:), allocatable, intent(out) :: intensity, iTheta, eTheta, ePhi
    integer :: n_datapoints

    character (len=FNAME_LENGTH) :: photFilename

    integer :: fileID, dimID, iThetaID, eThetaID, ePhiID, intensityID, ncstatus

    ncstatus = nf90_open(photFilename, NF90_NOWRITE, fileID)
    ncstatus = nf90_inq_dimid(fileid, 'n_datapoints', dimID)
    ncstatus = nf90_inquire_dimension(fileid, dimid, len=n_datapoints)

    allocate(intensity(n_datapoints), iTheta(n_datapoints), eTheta(n_datapoints), ePhi(n_datapoints))
    iTheta = 0.
    iTheta = 0.

    ncstatus = nf90_inq_varid(fileid, 'intensity', intensityID)
    ncstatus = nf90_get_var(fileid, intensityID, intensity)
    ncstatus = nf90_inq_varid(fileid, 'i_theta', iThetaID)
    ncstatus = nf90_get_var(fileid, iThetaID, iTheta)
    ncstatus = nf90_inq_varid(fileid, 'e_theta', eThetaID)
    ncstatus = nf90_get_var(fileid, eThetaID, eTheta)
    ncstatus = nf90_inq_varid(fileid, 'e_phi', ePhiID)
    ncstatus = nf90_get_var(fileid, ePhiID, ePhi)
    ncstatus = nf90_close(fileID)

    if (ncstatus /= 0) then 
       print *, 'Fatal error: could not read input file'
       call exit()
    end if

  end subroutine load_datapoints

  subroutine save(fname, samples, datapoints, nIll, mu, mu0, cosa)
    character(len = *), intent(in) :: fName
    integer :: samples, datapoints
    integer, dimension(:) :: nIll
    real(fd), dimension(:,:) :: mu, mu0, cosa

    integer            :: fileID, dimSamples, dimDatapoints, muID, mu0ID, cosaID, nIllID
    integer            :: status

    status = nf90_create(fName, NF90_CLOBBER, fileID)
    status = nf90_def_dim(fileID, "n_samples", samples, dimSamples)
    status = nf90_def_dim(fileID, "n_datapoints", datapoints, dimDatapoints)
    status = nf90_put_att(fileID, NF90_GLOBAL, "n_samples", [samples])
    status = nf90_put_att(fileID, NF90_GLOBAL, "n_datapoints", [datapoints])
    status = nf90_def_var(fileID, 'n_illuminated', NF90_INT, [dimDatapoints], nIllID)
    if (full_output) then
        status = nf90_def_var(fileID, 'mu', NF90_FLOAT, [dimSamples, dimDatapoints], muID)
        status = nf90_def_var(fileID, 'mu0', NF90_FLOAT, [dimSamples, dimDatapoints], mu0ID)
        status = nf90_def_var(fileID, 'cosa', NF90_FLOAT, [dimSamples, dimDatapoints], cosaID)
    end if

    status = nf90_enddef(fileID)
    status = nf90_put_var(fileID, nIllID, nIll) 
    if (full_output) then
        status = nf90_put_var(fileID, muID, mu) 
        status = nf90_put_var(fileID, mu0ID, mu0) 
        status = nf90_put_var(fileID, cosaID, cosa) 
    end if
    status = nf90_close(fileID)
  end subroutine save

end program fit_2hg

