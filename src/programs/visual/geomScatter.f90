!!--- BEGIN GPL ---
!!
!! geomScatter - a program to compute visual scatter in given geometries
!! Copyright (C) 2013 Olli Wilkman
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
!! Contributor(s): Olli Wilkman
!!
!!--- END GPL ---

!> Computes light scatterig from particulate media.
!!
!! \author Olli Wilkman
!! \version 0.1
!! \date 12.2.2013

program geomScatter
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

    implicit none

    character (len=FNAME_LENGTH) :: parFilename        = ""
    character (len=FNAME_LENGTH) :: outFilename        = "dscatter.scangles"
    character (len=FNAME_LENGTH) :: mediumFileName     = "medium.nc"
    character (len=100)          :: sDistribution      = "constant"

    integer                      :: gridResHorizontal  = 200
    integer                      :: gridResVertical    = 10

    real(fd)                     :: rho                = 0.50_fd
    real(fd)                     :: rhoAllowedError    = 0.05_fd

    integer                      :: resThetaE          = 45
    integer                      :: sampleSeed         = 0

    real(fd)                     :: mu                 = 1.0
    real(fd)                     :: w                  = 0.1

    integer                      :: nOrders            = 1
    character(len=100)           :: nSamplesPerOrder   = "1"
    integer, allocatable         :: nSamplesPerOrderTable(:)

    character (len=100)          :: brdfType           = "shadowing"
    character (len=100)          :: brdfPhaseFunction  = "constant"
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

    real(fd), external, pointer    :: Pf

    character (len=FNAME_LENGTH) :: photFilename = "photometry.nc"

    real(fd), dimension(:), allocatable :: iTheta, eTheta, ePhi, intensity
    real(fd), dimension(:,:), allocatable :: totalIntensity, totalIntensity_t
    integer :: n_datapoints, nFirstOrderSamples

    namelist /params/ &
        & resThetaE, gridResHorizontal, gridResVertical, &
        & outFilename, mediumFilename, photFilename, &
        & sDistribution, rho, rhoAllowedError, sampleSeed, &
        & nOrders, nSamplesPerOrder, &
        & mu, w, nPfParams, pfParams, &
        & brdfType, brdfPhaseFunction, &
        & rf_applyFields, rf_nFieldsPerMed, rf_spectrumType, rf_P, rf_std, &
        & nThreads, &
        & mediumMapRes, mediumDensitymapRes

    !! INITIALIZE SIMULATION PARAMETERS
    !!

    if(command_argument_count() == 0) then
        write(*,'("No parameter file given, using default values.")')
    else if(command_argument_count() == 1) then
        call get_command_argument(1, parFilename)
        write(*,'(("Using parameter file "),(A))') parFilename
    else
        write(*,'("Usage: geomScatter paramFile")')
        stop
    end if

    !! Read the namelist from the parameter file.
    if(parFilename /= "") then
        open(1,file=parFilename,status="old", form="formatted")
        read(1,NML=params)
        close(1)
    end if

    !!- Initialize random number generator
    call rnd_init(sampleSeed, nThreads, .true.)

    !! Read the phase function parameters into a table.
    read(PfParams,*) pfParamsTable(1:nPfParams)

    !! Read the number of samples for each scattering order into a table.
    allocate(nSamplesPerOrderTable(nOrders))
    read(nSamplesPerOrder,*) nSamplesPerOrderTable
    nFirstOrderSamples = nSamplesPerOrderTable(1)

    !! Select the media for the simulation from a file.
    write(*,'(("Using medium file "),(A))') mediumFilename
    call med_mediumFileOpen(Mf, mediumFilename)
    call med_mediumFileSelectMedia(Mf)

    if(rf_applyFields) then
        call rf_init(rndField, 1000, 5.0_fd)
        call utl_message("Roughness spectrum type used: " // rf_spectrumType)
        select case(rf_spectrumType)
        case('Gaussian')
            call rf_generateSpectrum (rndField, RF_SPECTRUM_GAUSSIAN, rf_P, rf_std)
        case('fBm')
            call rf_generateSpectrum (rndField, RF_SPECTRUM_FBM, rf_P, rf_std)
        end select
    else
        rf_nFieldsPerMed = 1
    end if

    !! READ PHOTOMETIC DATA
    call load_datapoints(intensity, iTheta, eTheta, ePhi, n_datapoints, photFilename)
    allocate(totalIntensity(n_datapoints,nOrders), totalIntensity_t(n_datapoints,nOrders))
    totalIntensity = 0.0

    !! RUN SIMULATION
    call utl_message("BRDF type used: " // brdfType)
    call utl_timer_init(timer, 5.0_fd, Mf%nSelectedMedia * rf_nFieldsPerMed * n_datapoints)

    do i = 1, Mf%nSelectedMedia
        call med_mediumFileRead(M, Mf, Mf%varSelection(i))
        call med_gridAssign(M, gridResHorizontal, gridResVertical)
        call med_gridFit(M)
        call med_updateStatistics(M)
        do iField = 1, rf_nFieldsPerMed
            if(rf_applyFields) then
                !!- Smooth the medium surface
                !!
                if(rf_spectrumType == 'constant') then
                    rndField%field = rf_P * M%height
                !!- Apply macroscale surface roughness.
                else
                    call rf_generateField(rndField, iField)
                    call med_maskHeight(m, real(rndField%field, fs))
               end if
            end if
            call cpu_time(sTime)
            ! START SIMULATING BASED ON BRDF TYPE
            totalIntensity_t = 0.0
            select case(brdfType)
            case("shadowing")
                call sampleGeometries(brdf_shadowing, nFirstOrderSamples, phase_function_constant, pfParamsTable, &
                                     iTheta, eTheta, ePhi, totalIntensity_t)
            case("Lambert")
                call sampleGeometries(brdf_Lambert, nFirstOrderSamples, phase_function_constant, pfParamsTable, &
                                     iTheta, eTheta, ePhi, totalIntensity_t)
            case("LommelSeeliger")
                select case(brdfPhaseFunction)
                case("constant")
                    call utl_message("Using constant phase function")
                    call sampleGeometries(brdf_LommelSeeliger, nFirstOrderSamples, phase_function_constant, &
                                          pfParamsTable, iTheta, eTheta, ePhi, totalIntensity_t)
                case("HG1")
                    call utl_message("Using HG phase function")
                    call sampleGeometries(brdf_LommelSeeliger, nFirstOrderSamples, phase_function_HG1, pfParamsTable, &
                                         iTheta, eTheta, ePhi, totalIntensity_t)
                case("HG2")
                    call utl_message("Using 2HG phase function")
                    call sampleGeometries(brdf_LommelSeeliger, nFirstOrderSamples, phase_function_HG2, pfParamsTable, &
                                         iTheta, eTheta, ePhi, totalIntensity_t)
                case default
                    call utl_fatal_error("Unsupported phase function type.")
                end select
            case("RadTransfer")
                call sampleGeometries(brdf_RadTrans, nFirstOrderSamples, phase_function_constant, pfParamsTable, &
                                     iTheta, eTheta, ePhi, totalIntensity_t)
            case default
                call utl_fatal_error("Unsupported BRDF type.")
            end select
            call cpu_time(cTime)
            totalIntensity = totalIntensity + totalIntensity_t
        end do
    end do
    write(*,'("Time taken: " ,(F8.3))') cTime - sTime

    ! Normalize by number of fields
    totalIntensity = totalIntensity / (real(MF%nSelectedMedia,fd)*(real(rf_nFieldsPerMed,fd)))

    !! WRITE FILES AND CLEAN UP
    call save(outFilename, nFirstOrderSamples * rf_nFieldsPerMed, n_datapoints, sum(totalIntensity,2))
    call med_mediumFileClose(Mf)
contains


subroutine sampleGeometries(f, nFirstOrderSamples, Pf, Pp, iTheta, eTheta, ePhi, totalIntensity_t)
    real(fd), external :: f
    real(fd), external :: Pf
    real(fd) :: Pp(10)
    integer, intent(in) :: nFirstOrderSamples
    real(fd), dimension(:), intent(in):: iTheta, eTheta, ePhi
    real(fd), dimension(size(iTheta),nOrders), intent(out) :: totalIntensity_t
    type(ray) :: rC
    type(ray) :: rS
    type(intersection_geometry) :: iSect
    real(fd), dimension(2,nFirstOrderSamples) :: samples
    real(fd), dimension(3) :: pSurface
    real(fd) :: dz, L(3)
    integer  :: nDatapoints, idp, iss
    logical  :: pFound, pLit

    nDatapoints = size(iTheta)

    call smpl_griddedSamples2D(samples, nFirstOrderSamples)

    samples = (samples * M%width - M%hWidth) * 0.5_fd
    dz      = M%grid%height - M%hMean - TRACE_EPS

    do idp = 1, nDatapoints
        !! Initialize the lightsource direction
        L = [sin(iTheta(idp)), 0.0_fd, cos(iTheta(idp))]

        !! Initialize the camera ray (pointing at surface)
        call ray_init(rC, RAY_TYPE_CAMERA)
        rC%D = [-cos(ePhi(idp))*sin(eTheta(idp)), -sin(ePhi(idp))*sin(eTheta(idp)), -cos(eTheta(idp))]

        ! Prevent the camera direction from being too close to horizontal.
        if(rC%D(3) > -1e-2_fd) then
            rC%D(3) = -1e-2
            call vec_normalize(rC%D)
        end if

        ! MAIN RAY LOOP
        do iss = 1, nFirstOrderSamples
            pSurface(1:2) = samples(:,iss)
            pFound = .false.
            do while (.not. pFound)
                ! Compute the ray starting point
                rC%P(1) = pSurface(1) - dz * (rC%D(1) / rC%D(3))
                rC%P(2) = pSurface(2) - dz * (rC%D(2) / rC%D(3))
                rC%P(3) = M%grid%height - TRACE_EPS
                rC%P(1) = modulo(rC%P(1) + M%hWidth, M%width) - M%hWidth
                rC%P(2) = modulo(rC%P(2) + M%hWidth, M%width) - M%hWidth

                rC%rayID = rC%rayID + RAY_ID_INCR

                ! Find the true intersection point of the camera ray and medium.
                ! This call needs the ray to point at the surface
                pFound = trc_traceNearest(M%grid, rC, iSect)

                ! If an intersection is found, compute the radiance,
                ! else randomize a new surface location and try again.
                if(pFound) then
                    call trc_gatherRadiance(M%grid, rC%D, L, iSect%P1+TRACE_EPS*iSect%N, &
                           & iSect%N, 1.0_fd/real(nFirstOrderSamples, fd), &
                           & nSamplesPerOrderTable, nOrders, 1, totalIntensity_t(idp,:), w, f, Pf, Pp)
                else
                    call rnd_generate_uniform(0.0_fd, 1.0_fd, pSurface(1:2))
                    pSurface(1:2) = (pSurface(1:2) * M%width - M%hWidth) * 0.5_fd
                end if
            end do
        end do
        call utl_timerIncrease(timer)
    end do
end subroutine sampleGeometries


! This subroutine loads observational geometries from a NetCDF file.
subroutine load_datapoints(intensity, iTheta, eTheta, ePhi, n_datapoints, photFilename)
    real(fd), dimension(:), allocatable, intent(out) :: intensity, iTheta, eTheta, ePhi
    integer :: n_datapoints
    character (len=FNAME_LENGTH) :: photFilename
    integer :: fileID, dimID, iThetaID, eThetaID, ePhiID, intensityID, ncstatus

    ! Open NetCDF file and find out the number of data points
    ncstatus = nf90_open(photFilename, NF90_NOWRITE, fileID)
    ncstatus = nf90_inq_dimid(fileid, 'n_datapoints', dimID)
    ncstatus = nf90_inquire_dimension(fileid, dimid, len=n_datapoints)

    ! Allocate tables based on number of data points
    allocate(intensity(n_datapoints), iTheta(n_datapoints), eTheta(n_datapoints), ePhi(n_datapoints))
    iTheta = 0.
    eTheta = 0.
    ePhi = 0.

    ! Load numbers from file to allocated tables
    ncstatus = nf90_inq_varid(fileid, 'intensity', intensityID)
    ncstatus = nf90_get_var(fileid, intensityID, intensity)
    ncstatus = nf90_inq_varid(fileid, 'i_theta', iThetaID)
    ncstatus = nf90_get_var(fileid, iThetaID, iTheta)
    ncstatus = nf90_inq_varid(fileid, 'e_theta', eThetaID)
    ncstatus = nf90_get_var(fileid, eThetaID, eTheta)
    ncstatus = nf90_inq_varid(fileid, 'e_phi', ePhiID)
    ncstatus = nf90_get_var(fileid, ePhiID, ePhi)
    ncstatus = nf90_close(fileID)

    ! Check for errors
    if (ncstatus /= 0) then
        print *, 'Fatal error: could not read input file'
        call exit()
    end if
end subroutine load_datapoints


! This subroutine saves the results of the computation.
subroutine save(fname, samples, datapoints, nIll)
    character(len = *), intent(in) :: fName
    integer :: samples, datapoints
    real(fd), dimension(:) :: nIll
    integer            :: fileID, dimSamples, dimDatapoints, nIllID
    integer            :: status

    status = nf90_create(fName, NF90_CLOBBER, fileID)
    status = nf90_def_dim(fileID, "n_rays", samples, dimSamples)
    status = nf90_def_dim(fileID, "n_datapoints", datapoints, dimDatapoints)
    status = nf90_put_att(fileID, NF90_GLOBAL, "n_rays", [samples])
    status = nf90_put_att(fileID, NF90_GLOBAL, "n_datapoints", [datapoints])
    status = nf90_def_var(fileID, 'intensity', NF90_FLOAT, [dimDatapoints], nIllID)
    status = nf90_enddef(fileID)
    status = nf90_put_var(fileID, nIllID, nIll)
    status = nf90_close(fileID)
end subroutine save

end program geomScatter

