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

!> Computes soft X-ray fluorescence of particulate media.
!!
!! \author Hannu Parviainen
!! \version 1.0
!! \date 20.01.2009
!!
!! \see simulation_method_firstorder
!!
program xrfpm
  use base
  use material
  use geometry
  use particle
  use hemisphere
  use random
  use medium
  use rfield
  use trace
  use sampler


  !$ use omp_lib

  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! VERSION CONSTANTS
  !!
  !!
  integer, parameter :: MAJOR_VERSION = 0
  integer, parameter :: MINOR_VERSION = 9
  integer, parameter :: FILE_VERSION  = 100

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! INPUT PARAMETERS
  !!
  !!- Files
  !!
  character (len=FNAME_LENGTH) :: parameterFile      = ""           !< The input parameter filename
  character (len=FNAME_LENGTH) :: resultFile         = "xrfpn.nc"   !< The output results filename
  character (len=FNAME_LENGTH) :: srcSpectrumFile    = ""           !< The input x-ray source spectrum filename
  character (len=FNAME_LENGTH) :: mediumFile         = ""           !< The medium filename

  !!- General
  !!
  logical                      :: verbose            = .True.
  integer                      :: nThreads           = 4

  !!- X-ray source
  !!
  character (len=100)          :: srcType            = "Spectrum"
  real(fd)                     :: srcLineEnergy      = 9e3
  real(fd)                     :: energyMin          = 0.0_fd
  real(fd)                     :: energyMax          = 1.0e4_fd


  !!- Incidence angles
  !!
  integer                      :: nThetaIn           = 1
  character (len=10000)        :: thetaIn            = "45.0"
  real(fd), allocatable        :: thetaInTable(:)

  !!- Simulation
  !!
  character (len=100)          :: method             = "mc_force"
  integer(il)                  :: nPosSamples        = 50000
  integer(il)                  :: nSrcSamples        = 1000
  integer                      :: nOrders            = 1
  integer                      :: sampleSeed         = 0
  real(fd)                     :: obsHalfAngle       = 5.0
  real(fd)                     :: obsSolidAngle      = 1.0

  !!- Output
  !!
  character(len=100)           :: outputType         = "hemisphere"
  integer                      :: hemisphereThetaRes = 25
  integer                      :: nEmergenceAngles   = 1
  character (len=10000)        :: thetaEm            = "45.0"
  real(fd), allocatable        :: thetaEmTable(:)
  character (len=10000)        :: phiEm              = "0.0"
  real(fd), allocatable        :: phiEmTable(:)

  !!- Material
  !!
  character (len=256)          :: matName            = "Default Fe/Ca"
  integer                      :: nElements          = 2
  character (len=10000)        :: elements           = "Fe Ca"
  character (len=10000)        :: elemFracs          = "0.5 0.5"
  character(2), allocatable    :: elementTable(:)
  real(fd), allocatable        :: elemFracTable(:)
  real(fd)                     :: molarDensity       = 1.0_fd

  !!- Grid accelerator
  !!
  integer                      :: gridResHorizontal  = 200
  integer                      :: gridResVertical    = 10

  !!- Random fields
  !!
  logical                      :: rf_applyFields    = .false.
  integer                      :: rf_nFieldsPerMed  = 1
  character(len=100)           :: rf_spectrumType   = "Gaussian"
  real(fd)                     :: rf_P              = 0.5_fd
  real(fd)                     :: rf_std            = 1.0_fd

  !!- Medium auxiliary data
  !!
  integer                      :: mediumMapRes        = 128
  integer                      :: mediumDensitymapRes = 300


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! SIMULATION VARIABLES
  !!
  !!
  !!- Constants
  !!
  !!
  !!- Internal variables
  !!
  type(gth_hemisphere)         :: fhs_simulated
  type(gth_hemisphere)         :: fhs_analytic
  real(fd), allocatable        :: ftb_simulated(:,:,:)
  real(fd), allocatable        :: ftb_analytic(:,:,:)
  !!
  type(spc_spectrum)           :: spc
  type(mat_material), target   :: mat
  type(med_medium)             :: med
  !!
  type(med_mediumFile)         :: medFile
  type(rf_randomField)         :: rndField
  !!
  real                         :: sTime, cTime, eTime, tHours, tMins, tSecs
  integer(il)                  :: i, iField, outType
  
  type(utl_timer)              :: simCpuTimer


  real,     dimension(:,:,:), allocatable :: mediumHeightMap
  real(fd), dimension(:,:,:), allocatable :: mediumDensityStructure


  namelist /params/ &
       & verbose,                                                                                 &
       & outputType,         resultFile,                                                          &
       & srcType,            srcSpectrumfile,  srcLineEnergy, nSrcSamples,                        &
       & nThetaIn,           thetaIn,                                                             &
       & method,             nOrders,          nPosSamples,        nThreads,                      &
       & nEmergenceAngles,   thetaEm,          phiEm,           hemisphereThetaRes,               &
       & matName,            nElements,        elements,        elemFracs,          molarDensity, &
       & energyMin,          energyMax,        obsHalfAngle,                                      &
       & mediumFile,                                                                              &
       & gridResHorizontal,  gridResVertical,                                                     &
       & rf_applyFields,     rf_nFieldsPerMed, rf_spectrumType, rf_P, rf_std,                     &
       & mediumMapRes,       mediumDensitymapRes


  !!--- INTERFACES ---
  !!
  interface 
     subroutine simulation_method_analytic(mat, med, spc, nPosSamples, inTheta, emTheta, emPhi, outType, timer, hs, tb)
       use material
       use geometry
       use hemisphere
       use medium

       implicit none

       type(mat_material)       :: mat
       type(med_medium)         :: med
       type(spc_spectrum)       :: spc
       integer(il)              :: nPosSamples
       real(fd), dimension(:)   :: inTheta
       real(fd), dimension(:)   :: emTheta
       real(fd), dimension(:)   :: emPhi
       integer                  :: outType
       type(utl_timer)          :: timer

       type(gth_hemisphere),       optional :: hs
       real(fd), dimension(:,:,:), optional :: tb
     end subroutine simulation_method_analytic
  end interface

  interface
     subroutine simulation_method_analytic_mc(mat, spc, nSpcSamples, inTheta, emTheta, emPhi, outType, hs, tb)
       use material
       use geometry
       use hemisphere
       use medium

       implicit none

       type(mat_material)       :: mat
       type(med_medium)         :: med
       type(spc_spectrum)       :: spc
       integer(il)              :: nSpcSamples
       real(fd), dimension(:)   :: inTheta
       real(fd), dimension(:)   :: emTheta
       real(fd), dimension(:)   :: emPhi
       integer                  :: outType

       type(gth_hemisphere),       optional :: hs
       real(fd), dimension(:,:,:), optional :: tb
     end subroutine simulation_method_analytic_mc
  end interface

  interface 
     subroutine simulation_method_firstorder(mat, med, spc, nPosSamples, nSrcSamples, inTheta, emTheta, emPhi, outType, verbose, timer, hs, tb, nThreads)
       use medium
       use material
       use geometry
       use hemisphere

       type(mat_material)       :: mat
       type(med_medium)         :: med
       type(spc_spectrum)       :: spc
       integer(il)              :: nPosSamples
       integer(il)              :: nSrcSamples
       real(fd), dimension(:)   :: inTheta
       real(fd), dimension(:)   :: emTheta
       real(fd), dimension(:)   :: emPhi
       integer                  :: outType
       logical                  :: verbose
       type(utl_timer)          :: timer

       type(gth_hemisphere),       optional :: hs
       real(fd), dimension(:,:,:), optional :: tb
       
       integer, optional        :: nThreads

     end subroutine simulation_method_firstorder
  end interface

  interface
     subroutine simulation_method_montecarlo(mat, med, spc, theta_i, hs, nPosSamples)
       use medium
       use material
       use geometry
       use hemisphere

       type(mat_material)     :: mat
       type(med_medium)       :: med
       type(spc_spectrum)     :: spc
       real(FD), dimension(:) :: theta_i
       type(gth_hemisphere)   :: hs
     end subroutine simulation_method_montecarlo
  end interface

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! READ THE INPUT FILE
  !!

  !$ call omp_set_num_threads(nThreads)

  if(command_argument_count() == 0) then
     write(*,'("No parameter file given, using default values.")')
  else if(command_argument_count() == 1) then
     call get_command_argument(1, parameterFile)
     write(*,'(("Using parameter file "),(A))') parameterFile
  else
     write(*,'("Usage: xrfpm parameterFile")')
     stop
  end if

  !! Read the namelist from the parameter file.
  !!
  if(parameterFile /= "") then
     open(1,file=parameterFile,status="old", form="formatted")
     read(1,NML=params)
     close(1)
  end if

  !! Read the incidence angles.
  !!
  allocate(thetaInTable(nThetaIn))
  read(thetaIn,*) thetaInTable
  thetaInTable = thetaInTable * CNV_DEGTORAD
  
  !! Read the emergence angles.
  !!
  allocate(thetaEmTable (nEmergenceAngles) )
  allocate(phiEmTable   (nEmergenceAngles) )
  read(thetaEm,*) thetaEmTable
  read(phiEm,  *) phiEmTable
  thetaEmTable  = thetaEmTable * CNV_DEGTORAD
  phiEmTable    = phiEmTable   * CNV_DEGTORAD

  !! Read the elements used by the material.
  !!
  allocate(elementTable(nElements))
  allocate(elemFracTable(nElements))

  read(elements,*) elementTable
  read(elemFracs,*) elemFracTable

  !! Read the output type
  !!
  if(outputType == 'hemisphere') then
     outType = OUTPUT_TYPE_HEMISPHERE
  else if(outputType == 'tabulated') then
     outType = OUTPUT_TYPE_TABULATED
  else
     call utl_fatal_error("Fatal Error: unknown output type " // outputType)
  end if

  !! Convert degrees to radians
  !!
  obsHalfAngle  = obsHalfAngle * CNV_DEGTORAD
  obsSolidAngle = (1.0_fd - cos(obsHalfAngle)) * TWO_PI

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! INITIALIZE VARIABLES
  !!

  !!- Initialize random number generator
  !!
  call rnd_init(1, nThreads, .true.)

  !!- Initialize material
  !!
  call mat_material_init(mat, matName, elementTable, molarDensity, elemFracTable, energyMax)

  !!- Initialize spectrum
  !!
  select case(srcType)
  case('spectrum')
     call spc_spectrum_read(spc, srcSpectrumfile, mat%edgeEnergyMin, energyMax)
  case('line')
     call spc_spectrum_initline(spc, srcLineEnergy, 1.0_fd)
  case default
     call utl_fatal_error("Error: Unknown spectrum type: " // srcType)
  end select

  !!- Initialize hemisphere
  !!
  if(outType == OUTPUT_TYPE_HEMISPHERE) then
     call gth_hemisphere_init(fhs_simulated, hemisphereThetaRes, nThetaIn, mat%nLines, GTH_TYPE_HS)
     call gth_hemisphere_init(fhs_analytic,  hemisphereThetaRes, nThetaIn, mat%nLines, GTH_TYPE_HS)
     write(fhs_simulated%name, '("res__fld_hemisphere_simulated")')
     write(fhs_analytic%name,  '("res__fld_hemisphere_analytical")' )
  else if(outType == OUTPUT_TYPE_TABULATED) then
     allocate(ftb_simulated(nThetaIn, nEmergenceAngles, mat%nLines))
     allocate(ftb_analytic(nThetaIn, nEmergenceAngles, mat%nLines))
     ftb_simulated = 0.0
     ftb_analytic  = 0.0
  end if


  !!- Initialize medium 
  !!
  call med_mediumFileOpen(medFile, mediumFile)
  call med_mediumFileSelectMedia(medFile)

  !!- Initialize random field generator
  !!
  if(rf_applyFields) then

     call rf_init(rndField, 500, 5.0_fd)

     select case(rf_spectrumType)
     case('Gaussian')
        call rf_generateSpectrum (rndField, RF_SPECTRUM_GAUSSIAN, rf_P, rf_std)
     case('fBm')
        call rf_generateSpectrum (rndField, RF_SPECTRUM_FBM, rf_P, rf_std)
     end select
  end if

  !!- Initialize extra medium data to be saved
  !!
  allocate(mediumHeightMap(medFile%nSelectedMedia, mediumMapRes, mediumMapRes))
  allocate(mediumDensityStructure(medFile%nSelectedMedia, mediumDensitymapRes, 2))


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! RUN SIMULATION
  !!

  !!- Initialize timer
  !!
  if(verbose) then
     call utl_timer_init(simCpuTimer, 5.0_fd, nThetaIn * nPosSamples * nSrcSamples * medFile%nSelectedMedia * rf_nFieldsPerMed)
  end if

  call cpu_time(sTime) 
  do i = 1, medFile%nSelectedMedia
     do iField = 1, rf_nFieldsPerMed
        call med_mediumFileRead(med, medFile, medFile%varSelection(i))

        call med_gridAssign(med, gridResHorizontal, gridResVertical)
        call med_gridFit(med)
        call med_updateStatistics(med)

        if(rf_applyFields) then

           !!- Smooth the medium surface: 
           !!  Remove a constant percentage of the medium height from the top. 
           !!
           if(rf_spectrumType == 'constant') then
              call med_clipHeight(med, rf_P, 'relative')

           !!- Apply macroscale surface roughness.
           !!
           else
              call rf_generateField(rndField, iField)
              call med_maskHeight(med, real(rndField%field, fs))
           end if
        end if
 
        mediumHeightMap(i,:,:) = cnt_grid3D_cmpHeightMap(med%grid, mediumMapRes)
        call med_computePorosityStructure(med, 5000, mediumDensityStructure(i,:,:))

        select case(method)
        case('mc_force')
           call simulation_method_montecarlo(mat, med, spc, thetaInTAble, fhs_simulated, nPosSamples)
        case('mc_peel')
           call simulation_method_montecarlo(mat, med, spc, thetaInTAble, fhs_simulated, nPosSamples)
        case('first_order')
           call simulation_method_firstorder(mat, med, spc, nPosSamples, nSrcSamples, thetaInTable, thetaEmTable, phiEmTable, outType, verbose, simCpuTimer, fhs_simulated, ftb_simulated, nThreads)
        case('analytic')
           call simulation_method_analytic(mat, med, spc, nPosSamples, thetaInTAble, thetaEmTable, phiEmTable, outType, simCpuTimer, fhs_simulated, ftb_simulated)
        case default
           call utl_fatal_error("Error: Unknown simulation method: " // method)
        end select
     end do
  end do
  call cpu_time(cTime)

  !! Divide the results by the number of media and number of fields
  !!
  if(outType == OUTPUT_TYPE_HEMISPHERE) then
     fhs_simulated%data = fhs_simulated%data / real(medFile%nSelectedMedia * rf_nFieldsPerMed,fd)
  else if(outType == OUTPUT_TYPE_TABULATED) then
     ftb_simulated      = ftb_simulated / real(medFile%nSelectedMedia * rf_nFieldsPerMed,fd)
  end if

  tHours = floor((cTime - sTime) / 3600.0)
  tMins  = floor(((cTime - sTime) - tHours * 3600.0) / 60.0)
  tSecs  = (cTime - sTime) - tHours * 3600.0 - tMins * 60.0

  write(*,'("Cpu time: " ,(I4)," h ",(I2)," min ",(F6.3)," s")') int(tHours), int(tMins), tSecs

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! WRITE OUTPUT
  !!

  call save_results(outType, mat, spc, resultFile, fhs_simulated, ftb_simulated)

contains

  !****f* xrfip/save_results
  ! NAME
  !   save_results
  !
  ! DESCRIPTION
  !   
  !
  ! INPUTS
  !   
  !   
  !   
  !   
  !
  ! TODO
  !   
  !   
  !
  ! SOURCE
  subroutine save_results(outType, mat, spc, filename, hs_sim, fDataSim)
    integer, intent(IN)                         :: outType
    type(gth_hemisphere), intent(IN), optional  :: hs_sim
    real(fd), dimension(:,:,:), optional        :: fDataSim
    type(mat_material), intent(IN)              :: mat
    type(spc_spectrum), intent(IN)              :: spc
    character (len = *)                         :: filename

    integer            :: i

    integer            :: fileID
    integer            :: fluorID, fluorAnID

    integer            :: dimMatMuEnergy, dimFlrLine, dimSpecEnergy, dimSpecICDF, dimThetaIn, dimEmAngle, dimMedNum, dimMedMapRes, dimMedDenRes, dimMedDenVars
    integer            :: idMatMuIon, idMatMuFluor, idMatMuFluorCdf, idMatMuRay, idMatMuExt, idMatMuEnergy
    integer            :: idSpecEnergy, idSpecIntensity, idSpecCdf, idSpecIcdf, idMedMap, idMedDen
    integer            :: idAnalytic

    integer, parameter :: RTP = NF90_DOUBLE
    integer, parameter :: RTS = NF90_FLOAT

    if(outType == OUTPUT_TYPE_HEMISPHERE .and. (.not. present(hs_sim))) then
       call utl_fatal_error("Incorrect data given in saveTabulated")
    end if

   if(outType == OUTPUT_TYPE_TABULATED .and. (.not. present(fDataSim))) then
       call utl_fatal_error("Incorrect data given in saveTabulated")
    end if


    call gth_nfCheck( nf90_create(filename, NF90_CLOBBER, fileID) )

    !!--- CORE ---
    !!
    !!- Core attributes
    !!
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "author",           "Hannu Parviainen" ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "program",          "xrfip"            ) ) !! Program name
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "version",          "1.0"              ) ) !! Program version
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "file_version",     FILE_VERSION       ) ) !! File format version
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "simulation_type",  method             ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "result_type",      OutputType         ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "medium_type",      "infinite_plane"   ) )


    !!--- SPECTRUM ---
    !!
    !!- Spectrum attributes
    !!
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "spc__n_spectrum_samples",  spc%nPoints     ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "spc__spectrum_type",       srcType         ) )
    !!
    !!- Spectrum dimensions
    !!
    call gth_nfCheck( nf90_def_dim(fileID,    "spc__spectrum_energy",              spc%nPoints,   dimSpecEnergy  ) )
    if(srcType == 'spectrum') then
       call gth_nfCheck( nf90_def_dim(fileID, "spc__spectrum_inverse_cdf_energy",  size(spc%icdf), dimSpecICDF    ) )
    end if
    !!
    !!- Spectrum variables
    !!
    call gth_nfCheck( nf90_def_var(fileID, "spc__spectrum_energy",               RTP, [dimSpecEnergy], idSpecEnergy) )
    call gth_nfCheck( nf90_def_var(fileID, "spc__spectrum_intensity",            RTP, [dimSpecEnergy], idSpecIntensity) )
    if(srcType == 'spectrum') then
       call gth_nfCheck( nf90_def_var(fileID, "spc__spectrum_cdf",               RTP, [dimSpecEnergy], idSpecCdf) )
       call gth_nfCheck( nf90_def_var(fileID, "spc__spectrum_inverse_cdf",       RTP, [dimSpecICDF],   idSpecIcdf) )
    end if


    !!--- material ---
    !!
    !!- Material attributes
    !!
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "mat__n_elements",                 nElements    ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "mat__n_fluorescence_lines",       mat%nLines   ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "mat__n_material_energy_samples",  mat%nE       ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "mat__element_names",              elements     ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "mat__element_fractions",          elemFracs    ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "mat__molar_density",              molarDensity ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "mat__fluorescence_line_energy",   mat%lEnergy  ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "mat__absorbtion_edge_energy",     mat%eEnergy  ) )
    !!
    !!- Dimensions
    !!
    call gth_nfCheck( nf90_def_dim(fileID,              "mat__fluorescence_line",          mat%nLines, dimFlrLine     ) )
    call gth_nfCheck( nf90_def_dim(fileID,              "mat__material_energy",            mat%nE,     dimMatMuEnergy ) )
    !!
    !!- Variables
    !!
    call gth_nfCheck( nf90_def_var(fileID, "mat__photoionization_coefficient",   RTP, [dimMatMuEnergy],            idMatMuIon) )
    call gth_nfCheck( nf90_def_var(fileID, "mat__fluorescence_line_coefficient", RTP, [dimFlrLine,dimMatMuEnergy], idMatMuFluor) )
    call gth_nfCheck( nf90_def_var(fileID, "mat__fluorescence_line_cdf",         RTP, [dimFlrLine,dimMatMuEnergy], idMatMuFluorCdf))
    call gth_nfCheck( nf90_def_var(fileID, "mat__rayleigh_coefficient",          RTP, [dimMatMuEnergy],            idMatMuRay) )
    call gth_nfCheck( nf90_def_var(fileID, "mat__extinction_coefficient",        RTP, [dimMatMuEnergy],            idMatMuExt) )
    call gth_nfCheck( nf90_def_var(fileID, "mat__material_energy",               RTP, [dimMatMuEnergy],            idMatMuEnergy) )

    
    !!--- PARTICULATE MEDIUM ---
    !!
    !!- Attributes
    !!
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "med__medium_filename",             mediumFile             ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "med__n_media",                     medFile%nSelectedMedia ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "med__n_random_fields_per_medium",  rf_nFieldsPerMed       ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "med__rf_spectrum_type",            rf_spectrumType        ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "med__rf_parameter",                rf_P                   ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "med__rf_standard_deviation",       rf_std                 ) )
    !!
    !!- Dimensions
    !!
    Call gth_nfCheck( nf90_def_dim(fileID, "med__media",                        medFile%nSelectedMedia,  dimMedNum      ) )
    Call gth_nfCheck( nf90_def_dim(fileID, "med__medium_map_resolution",        mediumMapRes,            dimMedMapRes   ) )
    Call gth_nfCheck( nf90_def_dim(fileID, "med__medium_densitymap_resolution", mediumDensitymapRes,     dimMedDenRes   ) )
    Call gth_nfCheck( nf90_def_dim(fileID, "med__medium_densitymap_variables",  2,                       dimMedDenVars  ) )
    !!
    !!- Variables
    !!
    Call gth_nfCheck( nf90_def_var(fileID, "med__medium_heightmap",           RTS, [dimMedNum,dimMedMapRes,dimMedMapRes],  idMedMap) )
    Call gth_nfCheck( nf90_def_var(fileID, "med__medium_densitymap",          RTS, [dimMedNum,dimMedDenRes,dimMedDenVars], idMedDen) )


    !!--- RESULTS ---
    !!
    !!- Result attributes
    !!
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "res__n_incidence_angles",       nThetaIn         ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "res__n_emergence_angles",       nEmergenceAngles ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "res__n_samples",                nPosSamples      ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "res__n_orders",                 nOrders          ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "res__incident_theta",           thetaInTable     ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "res__emergent_theta",           ThetaEmTable     ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "res__emergent_phi",             PhiEmTable       ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "res__bin_solid_angle",          obsSolidAngle    ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "res__bin_half_angle",           obsHalfAngle     ) )
    !!
    !!- Result dimensions
    !!
    call gth_nfCheck( nf90_def_dim(fileID, "res__incident_theta",        nThetaIn,             dimThetaIn     ) )
    call gth_nfCheck( nf90_def_dim(fileID, "res__emergence_angle",       nEmergenceAngles,     dimEmAngle     ) )
    !!
    !!- Result Variables
    !!
    if(outType == OUTPUT_TYPE_HEMISPHERE) then
       fluorID   = gth_writeNetcdfHeader(hs_sim, fileID)
    else
       call gth_nfCheck( nf90_def_var(fileID, "res__fld_table_simulated",   RTP, [dimThetaIn,dimEmAngle,dimFlrLine], fluorID   ) )
    end if

    call gth_nfCheck( nf90_enddef(fileID) )

    !!--- WRITE VARIABLES ---
    !!
    if(outType == OUTPUT_TYPE_HEMISPHERE) then
       call gth_nfCheck( nf90_put_var(fileID, fluorID,   hs_sim%data ) )
    else
       call gth_nfCheck( nf90_put_var(fileID, fluorID,   fDataSim    ) )
    end if

    call gth_nfCheck( nf90_put_var(fileID, idMatMuIon,       mat%muPhotoIon) )
    call gth_nfCheck( nf90_put_var(fileID, idMatMuFluor,     mat%muFluorLine) )
    call gth_nfCheck( nf90_put_var(fileID, idMatMuFluorCdf,  mat%muFluorLineCDF) )
    call gth_nfCheck( nf90_put_var(fileID, idMatMuRay,       mat%muRayleigh) )
    call gth_nfCheck( nf90_put_var(fileID, idMatMuExt,       mat%muExtTotal) )
    call gth_nfCheck( nf90_put_var(fileID, idMatMuEnergy,    mat%E) )

    Call gth_nfCheck( nf90_put_var(fileID, idMedMap,      mediumHeightMap) )
    Call gth_nfCheck( nf90_put_var(fileID, idMedDen,      mediumDensityStructure) )

    call gth_nfCheck( nf90_put_var(fileID, idSpecEnergy,     spc%E) )
    call gth_nfCheck( nf90_put_var(fileID, idSpecIntensity,  spc%I) )
    if(srcType == 'Spectrum') then
       call gth_nfCheck( nf90_put_var(fileID, idSpecCdf,     spc%cdf) )
       call gth_nfCheck( nf90_put_var(fileID, idSpecIcdf,    spc%icdf) )
    end if

    call gth_nfCheck( nf90_close(fileID) )

  end subroutine save_results
  !******


end program xrfpm

