!!--- BEGIN GPL --- 
!!
!! hemiScatter
!! Copyright (C) 2014 Olli Wilkman
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
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.    If not, see <http://www.gnu.org/licenses/>.
!!
!! Contributor(s): Olli Wilkman
!!
!!--- END GPL ---

!> Computes light scattering from particulate media.
!!
!! \author Olli Wilkman
!! \version 0.1
!! \date 20.5.2014

program vScatter
    use iso_c_binding
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

    character (len=FNAME_LENGTH) :: parFilename          = ""
    character (len=FNAME_LENGTH) :: hsCarPlotFilename    = "vScatter_car.eps"
    character (len=FNAME_LENGTH) :: hsSphPlotFilename    = "vScatter_sph.eps"
    character (len=FNAME_LENGTH) :: hsPlotFilename       = ""
    character (len=FNAME_LENGTH) :: hsFilename           = "hemisphere.nc"
    character (len=FNAME_LENGTH) :: mediumFileName       = "medium.nc"
    character (len=100)          :: brdfType             = "shadowing"
    character (len=100)          :: brdfPhaseFunction    = "constant"

    integer                      :: gridResHorizontal    = 200
    integer                      :: gridResVertical      = 10

    integer                      :: resTheta            = 45
    integer                      :: sampleSeed           = 0


    integer                      :: nOrders              = 1
    character(len=100)           :: nSamplesPerOrder     = "1"
    integer, allocatable         :: nSamplesPerOrderTable(:)

    real(fd)                     :: mu                   = 1.0
    real(fd)                     :: w                    = 0.1

    integer                      :: nPfParams            = 1
    character (len=10000)        :: pfParams             = "0.0"
    real(fd)                     :: pfParamsTable(10)

    logical                      :: rf_applyFields       = .false.
    integer                      :: rf_nFieldsPerMed     = 1
    character(len=100)           :: rf_spectrumType      = "Gaussian"
    real(fd)                     :: rf_P                 = 0.5_fd
    real(fd)                     :: rf_std               = 1.0_fd

    integer                      :: mediumMapRes         = 256
    integer                      :: mediumDensitymapRes  = 300
    logical                      :: saveMediumMap        = .true.

    integer                      :: nThreads             = 1

    type(gth_hemisphere)    :: H
    type(med_medium)        :: M
    type(med_mediumFile)    :: Mf
    type(rf_randomField)    :: rndField
    type(utl_timer)         :: timer

    real, dimension(:,:,:), allocatable :: mediumHeightMap
    real(fd), dimension(:,:,:), allocatable :: mediumDensityStructure

    real(fd), external, pointer :: Pf


    real(fd) :: dTheta

    real    :: sTime, cTime, eTime
    integer :: i, iField

    namelist /params/ &
             & resTheta, gridResHorizontal, gridResVertical, &
             & hsCarPlotFilename, hsSphPlotFilename, hsFilename, mediumFilename, &
             & sampleSeed, &
             & nOrders, nSamplesPerOrder, &
             & mu, w, nPfParams, pfParams, &
             & brdfType, brdfPhaseFunction, &
             & rf_applyFields, rf_nFieldsPerMed, rf_spectrumType, rf_P, rf_std, &
             & nThreads, &
             & mediumMapRes, mediumDensitymapRes, saveMediumMap


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
        write(*,'("Usage: vScatter paramFile")')
        stop
    end if

    !! Read the namelist from the parameter file.
    !!
    if(parFilename /= "") then
        open(1,file=parFilename,status="old", form="formatted")
        read(1,NML=params)
        close(1)
    end if

    !$ call omp_set_num_threads(nThreads)

    !!- Initialize random number generator
    !!
    call rnd_init(sampleSeed, nThreads, .true.)

    !! Generate the angles of incidence into a table.
    !!
    dTheta = HALF_PI/real(resTheta, fd)

    !! Read the number of samples for each scattering order into a table.
    !!
    allocate(nSamplesPerOrderTable(nOrders))
    read(nSamplesPerOrder,*) nSamplesPerOrderTable

    !! Read the phase function parameters into a table.
    !!
    !allocate(pfParamsTable(nPfParams))
    read(PfParams,*) pfParamsTable(1:nPfParams)


    !! Initialize hemisphere
    !!
    call gth_hemisphere_init(H, resTheta, resTheta, nOrders, GTH_TYPE_QS)
    write(H%name,'("Hemisphere")')

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
    end if

    call utl_timer_init(timer, 5.0_fd, Mf%nSelectedMedia * rf_nFieldsPerMed * H%nCells)
         
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! RUN SIMULATION
    !!

    allocate(mediumHeightMap(MF%nSelectedMedia, mediumMapRes, mediumMapRes))
    allocate(mediumDensityStructure(MF%nSelectedMedia, mediumDensitymapRes, 2))

    call utl_message("BRDF type used: " // brdfType)

    do i = 1, MF%nSelectedMedia
        do iField = 1, rf_nFieldsPerMed
            call med_mediumFileRead(M, Mf, Mf%varSelection(i))
            call med_gridAssign(M, gridResHorizontal, gridResVertical)
            call med_gridFit(M)
            call med_updateStatistics(M)

            if(rf_applyFields) then

                !!- Smooth the medium surface: 
                !!    Remove a constant percentage of the medium height from the top. 
                !!
                if(rf_spectrumType == 'constant') then
                    rndField%field = rf_P * M%height
                else
                    call rf_generateField(rndField, iField)
                    call med_maskHeight(m, real(rndField%field, fs))
                end if
            end if

            mediumHeightMap(i,:,:) = cnt_grid3D_cmpHeightMap(M%grid, mediumMapRes)
            call med_computePorosityStructure(M, 5000, mediumDensityStructure(i,:,:))

            !call med_gridFit(M)
            !call med_updateStatistics(M)

            !call utl_message("Beginning trace.")
 
            call cpu_time(sTime)
            
            select case(brdfType)
            case("shadowing")
                call sampleHemisphere(brdf_shadowing, nSamplesPerOrderTable(1), phase_function_constant, pfParamsTable)
            case("Lambert")
                call sampleHemisphere(brdf_Lambert, nSamplesPerOrderTable(1), phase_function_constant, pfParamsTable)
            case("LommelSeeliger")
                select case(brdfPhaseFunction)
                case("constant")
                    call sampleHemisphere(brdf_LommelSeeliger, nSamplesPerOrderTable(1), phase_function_constant, pfParamsTable)
                case("HG1")
                    call sampleHemisphere(brdf_LommelSeeliger, nSamplesPerOrderTable(1), phase_function_HG1, pfParamsTable)
                case("HG2")
                    call sampleHemisphere(brdf_LommelSeeliger, nSamplesPerOrderTable(1), phase_function_HG2, pfParamsTable)
                case default
                    call utl_fatal_error("Unsupported BRDF type.")
                end select
            case("RadTransfer")
                call sampleHemisphere(brdf_RadTrans, nSamplesPerOrderTable(1), phase_function_constant, pfParamsTable)
            case default
                call utl_fatal_error("Unsupported BRDF type.")
            end select

            call cpu_time(cTime)
         end do
    end do

    H % data = H % data / (real(MF%nSelectedMedia,fd)*(real(rf_nFieldsPerMed,fd)))

    write(*,'("Time taken: " ,(F8.3))') cTime - sTime

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!
    !! WRITE FILES AND CLEAN UP
    !!
    call saveHemisphere(H, hsFilename)
    call med_mediumFileClose(Mf)


contains
    !! TODO: Include single scattering phase function to computations!
    !!
    subroutine sampleHemisphere(f, nSamples, Pf, Pp)
        real(fd), external :: f
        integer            :: nSamples
        real(fd), external :: Pf
        real(fd)           :: Pp(10)

        type(ray) :: rC
        type(intersection_geometry) :: iSect

        real(fd), dimension(2,nSamples) :: samples
        real(fd), dimension(3)          :: pSurface(3)
        real(fd), dimension(3) :: D
        real(fd)                        :: dz, thetaIn, thetaInOffset
        integer :: i, j, k, iTheta
        logical :: pFound, pLit

        real :: tstRnd

        call smpl_griddedSamples2D(samples, nSamples)

        samples = (samples * M%width - M%hWidth) * 0.5_fd

        ! $omp parallel default(none)                                 &
        ! $omp shared(M, H, D, samples, f, w, Pf, Pp, nOrders, nSamplesPerOrderTable, resTheta, thetaInTable, timer) &
        ! $omp private(rC, dz, thetaIn, pSurface, i, j, k, iTheta, pFound, pLit, iSect, tstRnd)
 
        dz = M%grid%height - M%hMean - TRACE_EPS

        pSurface(3) = M%hMean
        
        D = 0.0_fd

        call ray_init(rC, RAY_TYPE_CAMERA)

        ! $omp do schedule(dynamic) 
        do j = 1, H % resTheta
            do k = 1, H % resPhi(j)
                ! $omp do schedule(dynamic) 
                do i= 1, nSamplesPerOrderTable(1)
                    !! Select the sample point of the incident camera ray from
                    !! the mean medium surface.
                    !!
                    pSurface(1:2) = samples(:,i)
                    call RANDOM_NUMBER(tstRnd)
                    thetaInOffset = dTheta*tstRnd

                    pFound = .false.
                    do while (.not. pFound)            
                        !! Find the intersection of the camera ray and the top of 
                        !! the bounding box of the periodic medium.
                        !!
                        rC%D = gth_cellRandomSampleCar(H, j, k)

                        if(h%type == GTH_TYPE_QS) then
                            call RANDOM_NUMBER(tstRnd)
                            if(tstRnd > 0.5) rC%D(2) = -rC%D(2)
                        end if

                        if(rC%D(3) < 1e-2_fd) then
                             rC%D(3) = 1e-2
                             call vec_normalize(rC%D)
                        end if

                        rC % P(1)  = pSurface(1) + dz * (rC%D(1) / rC%D(3))
                        rC % P(2)  = pSurface(2) + dz * (rC%D(2) / rC%D(3))
                        rC % P(3)  = M%grid%height - TRACE_EPS
                                    
                        rC % P(1)  = modulo(rC%P(1)+M%hWidth, M%width) - M%hWidth
                        rC % P(2)  = modulo(rC%P(2)+M%hWidth, M%width) - M%hWidth
                                    
                        rC % rayID = rC%rayID + RAY_ID_INCR

                        !! Negate the ray.
                        !!
                        rC % D = - rC%D

                        !! Find the true intersection point of the camera ray and medium.
                        !!
                        pFound = trc_traceNearest(M%grid, rC, iSect)

                        !! If an intersection is found, check if the point is shadowed.
                        !!
                        if(pFound) then
                            do iTheta = 1, resTheta
                                thetaIn = (iTheta-1)*dTheta + thetaInOffset
                                D(1) = sin(thetaIn)
                                D(3) = cos(thetaIn)
                                call trc_gatherRadiance(M%grid, rC%D, D, iSect%P1 + TRACE_EPS * iSect%N, &
                                    & iSect%N, 1.0_fd / real(nSamplesPerOrderTable(1), fd), &
                                    & nSamplesPerOrderTable, nOrders, 1, H % data(iTheta, h%cIdx(j)+k-1,:), w, f, Pf, Pp)
                            end do
                        else
                            call rnd_generate_uniform(0.0_fd, 1.0_fd, pSurface(1:2))
                            pSurface(1:2) = (pSurface(1:2) * M%width - M%hWidth) * 0.5_fd
                        end if
                    end do
                end do
                ! $omp end do
                call utl_timerIncrease(timer)    
            end do
        end do
        ! $omp end do
        ! $omp end parallel

    end subroutine sampleHemisphere




    subroutine saveHemisphere(h, fName)
        type(gth_hemisphere), INTENT(IN) :: h
        character (len = *) :: fName

        type(gth_hsFile) :: hf
        integer :: i

        integer :: dimMedNum, dimMedMapRes, dimMedDenRes, idMedMap, idMedDen

        call gth_hsFileOpen(hf, fName, "w")
        call gth_hsWriteHeader(hf, h, "Olli Wilkman", "hemiScatter", "0.1")
 
        Call gth_nfCheck( nf90_def_dim(hf%fileID, "Number_of_media", MF%nSelectedMedia, dimMedNum))
        Call gth_nfCheck( nf90_def_dim(hf%fileID, "Medium_map_resolution", mediumMapRes, dimMedMapRes))
        Call gth_nfCheck( nf90_def_dim(hf%fileID, "Medium_densitymap_resolution", mediumDensitymapRes, dimMedDenRes))

        if(saveMediumMap) then
            Call gth_nfCheck( nf90_def_var(hf%fileID, "Medium_heightmap", NF90_FLOAT, [dimMedNum,dimMedMapRes,dimMedMapRes], idMedMap) )
            Call gth_nfCheck( nf90_def_var(hf%fileID, "Medium_densitymap", NF90_FLOAT, [dimMedNum,dimMedDenRes,2], idMedDen) )
        end if

        call gth_hsAddAttS(hf, NF90_GLOBAL, "Brdf_type", brdfType)
        call gth_hsAddAttS(hf, NF90_GLOBAL, "medium_filename", mediumFileName)
        call gth_hsAddAttF(hf, NF90_GLOBAL, "Single_scattering_albedo", [w])
        call gth_hsAddAttF(hf, NF90_GLOBAL, "extinction_coefficient", [mu])
        call gth_hsAddAttF(hf, NF90_GLOBAL, "nMedia", [real(MF%nSelectedMedia,fd)])

        if(rf_applyfields) then
            call gth_hsAddAttS(hf, NF90_GLOBAL, "srfType", rf_spectrumType)
            call gth_hsAddAttF(hf, NF90_GLOBAL, "nFields", [real(rf_nFieldsPerMed,fd)])
            call gth_hsAddAttF(hf, NF90_GLOBAL, "srfStd", [rf_P])
            call gth_hsAddAttF(hf, NF90_GLOBAL, "srfRoughP", [rf_std])
        end if

        call gth_hsSaveData(h, hf)

        if(saveMediumMap) then
            Call gth_nfCheck( nf90_put_var(hf%fileID, idMedMap, mediumHeightMap))
            Call gth_nfCheck( nf90_put_var(hf%fileID, idMedDen, mediumDensityStructure))
        end if

        call gth_hsFileClose(hf)
    end subroutine saveHemisphere

end program vScatter

