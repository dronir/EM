!!--- BEGIN GPL --- 
!!
!! xrfip -  a program to compute soft X-ray fluorescence of an infinite slab.
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

!> A program to compute soft X-ray fluorescence of an infinite slab.
!!
!! \author  Hannu Parviainen
!! \version              1.0
!! \date           7.01.2009
!!
program xrfip
  use base
  use material
  use geometry
  use particle
  use hemisphere
  use random
  use bfdf

  !$ use omp_lib

  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! INPUT PARAMETERS
  !!

  !! FILES
  !!
  character (len=FNAME_LENGTH) :: parameterFile      = ""             !< The parameter file filename
  character (len=FNAME_LENGTH) :: resultFile         = "xrfiphs.nc"   !< The result file filename
  character (len=FNAME_LENGTH) :: srcSpectrumfile    = ""             !< The Spectrum file filename


  !! GENERAL
  !!
  logical                      :: verbose            = .True.         !< Print simulation runtime information.
  integer                      :: nThreads           = 1              !< Number of OpenMP threads


  !! X-RAY SOURCE 
  !!
  character (len=100)          :: srcType            = "spectrum"     !< Type of the source. Can be either "spectrum" or "line"
  real(fd)                     :: srcLineEnergy      = 9e3            !< Energy of the line source
  real(fd)                     :: energyMin          = 0.0_fd         !< Minimum energy
  real(fd)                     :: energyMax          = 1.0e4_fd       !< Maximum energy


  !! INCIDENCE ANGLES
  !!
  integer                      :: nThetaIn           = 1              !< Number of incidence angles
  character (len=10000)        :: thetaIn            = "45.0"         !< Incidence angles in degrees as a string
  real(fd), allocatable        :: thetaInTable(:)                     !< Incidence angles in radians


  !! SIMULATION
  !!
  character (len=100)          :: method             = "mc_force"     !< Simulation method
  integer(il)                  :: nPosSamples        = 50000          !< Number of position samples
  integer(il)                  :: nSrcSamples        = 1000           !< Number of spectrum samples
  integer                      :: nOrders            = 1              !< Number of fluorescence orders
  integer                      :: sampleSeed         = 0              !< Sample seed

  !! OUTPUT
  !!
  character(len=100)           :: outputType         = "hemisphere"   !< Type of the output file. Can be either "hemisphere" or "tabulated" 
  integer                      :: hemisphereThetaRes = 25             !< Theta resolution of the hemisphere

  integer                      :: nEmergenceAngles   = 1              !< Number of emergence angles
  character (len=10000)        :: thetaEm            = "45.0"         !< Emergent theta angles in degrees as a string
  real(fd), allocatable        :: thetaEmTable(:)                     !< Emergent theta angles in radians
  character (len=10000)        :: phiEm              = "0.0"          !< Emergent phi angles in degrees as a string
  real(fd), allocatable        :: phiEmTable(:)                       !< Emergent phi angles in radians

  real(fd)                     :: obsHalfAngle       = 5.0            !< Half angle of the gathering element for tabulated Monte-Carlo simulations
  real(fd)                     :: obsSolidAngle      = 1.0            !< Solid angle of the gathering element for tabulated Monte-Carlo simulations


  !! material
  !!
  character (len=256)          :: matName            = "Default Fe/Ca" !< Material name
  integer                      :: nElements          = 2               !< Number of elements in the material
  character (len=10000)        :: elements           = "Fe Ca"         !< List of elements as a string
  character (len=10000)        :: elemFracs          = "0.5 0.5"       !< Elemental abundancy fractions as a string
  character(2), allocatable    :: elementTable(:)                      !<
  real(fd), allocatable        :: elemFracTable(:)                     !< Elemental abundancy fractions
  real(fd)                     :: molarDensity       = 1.0_fd          !< Total molar density of the material


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! SIMULATION VARIABLES
  !!

  !- Constants
  !


  !- Internal variables
  !
  type(gth_hemisphere)         :: fhs_simulated
  type(gth_hemisphere)         :: fhs_analytic
  real(fd), allocatable        :: ftb_simulated(:,:,:)
  real(fd), allocatable        :: ftb_analytic(:,:,:)

  type(spc_spectrum)           :: spc
  type(mat_material), target   :: mat

  real                         :: sTime, cTime, eTime, tHours, tMins, tSecs
  integer(il)                  :: outType, i

  namelist /params/ &
       & verbose,                                                                                 &
       & outputType,         resultFile,                                                          &
       & srcType,            srcSpectrumfile,  srcLineEnergy,   nSrcSamples,                      &
       & nThetaIn,           thetaIn,                                                             &
       & method,             nOrders,          nPosSamples,        nThreads,                         &
       & nEmergenceAngles,   thetaEm,          phiEm,           hemisphereThetaRes,               &
       & matName,            nElements,        elements,        elemFracs,          molarDensity, &
       & energyMin,          energyMax,        obsHalfAngle


  !!--- INTERFACES ---
  !!
  interface 
     subroutine infslab_method_analytic(mat, spc, theta_i, theta_e, outType, hs, tb)
       use material
       use geometry
       use hemisphere

       type(mat_material)       :: mat
       type(spc_spectrum)       :: spc
       real(fd), dimension(:)   :: theta_i
       real(fd), dimension(:)   :: theta_e
       integer                  :: outType

       type(gth_hemisphere),       optional :: hs
       real(fd), dimension(:,:,:), optional :: tb

     end subroutine infslab_method_analytic
  end interface

  interface
     subroutine infslab_method_analytic_mc(mat, spc, nSpcSamples, theta_i, theta_e, outType, hs, tb)   
       use material
       use geometry
       use hemisphere

       type(mat_material)       :: mat
       type(spc_spectrum)       :: spc
       real(fd), dimension(:)   :: theta_i
       real(fd), dimension(:)   :: theta_e
       integer                  :: outType
       integer                  :: nSpcSamples

       type(gth_hemisphere),       optional :: hs
       real(fd), dimension(:,:,:), optional :: tb
     end subroutine infslab_method_analytic_mc
  end interface

  interface 
     subroutine infslab_method_firstorder(mat, spc, nPosSamples, nSrcSamples, inTheta, emTheta, emPhi, outType, verbose, hs, tb)
       use material
       use geometry
       use hemisphere

       type(mat_material)       :: mat
       type(spc_spectrum)       :: spc
       integer(il)              :: nPosSamples
       integer(il)              :: nSrcSamples
       real(fd), dimension(:)   :: inTheta
       real(fd), dimension(:)   :: emTheta
       real(fd), dimension(:)   :: emPhi
       integer                  :: outType
       logical                  :: verbose

       type(gth_hemisphere),       optional :: hs
       real(fd), dimension(:,:,:), optional :: tb
     end subroutine infslab_method_firstorder
  end interface

  interface
     subroutine infslab_method_montecarlo(mat, spc, nPosSamples, nSrcamples, nOrders, obsHalfAngle, theta_i, theta_e, outType, hs, tb)
       use material
       use geometry
       use hemisphere

       type(mat_material)         :: mat
       type(spc_spectrum)         :: spc
       integer(il)                :: nPosSamples
       integer(il)                :: nSrcSamples
       integer(il)                :: nOrders
       real(fd)                   :: obsHalfAngle
       real(fd), dimension(:)     :: theta_i
       real(fd), dimension(:)     :: theta_e
       integer                    :: outType
       type(gth_hemisphere)       :: hs
       real(fd), dimension(:,:,:) :: tb
     end subroutine infslab_method_montecarlo
  end interface

  interface
     subroutine infslab_method_montecarlo_force(mat, spc, nPosSamples, theta_i, theta_e, outType, hs, tb)
       use material
       use geometry
       use hemisphere

       type(mat_material)         :: mat
       type(spc_spectrum)         :: spc
       integer(il)                :: nPosSamples
       real(fd), dimension(:)     :: theta_i
       real(fd), dimension(:)     :: theta_e
       integer                    :: outType
       type(gth_hemisphere)       :: hs
       real(fd), dimension(:,:,:) :: tb
     end subroutine infslab_method_montecarlo_force
  end interface


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! READ THE INPUT FILE
  !!

  if(command_argument_count() == 0) then
     write(*,'("No parameter file given, using default values.")')
  else if(command_argument_count() == 1) then
     call get_command_argument(1, parameterFile)
     write(*,'(("Using parameter file "),(A))') parameterFile
  else
     write(*,'("Usage: xr_InfSlab paramFile")')
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! RUN SIMULATION
  !!

  call cpu_time(sTime) 
  select case(method)
  case('mc_force')
     call infslab_method_montecarlo_force(mat, spc, nPosSamples, thetaInTable, ThetaEmTable, outType, fhs_simulated, ftb_simulated)
  case('mc_peel')
     call infslab_method_montecarlo(mat, spc, nPosSamples, nSrcSamples, nOrders, obsHalfAngle, thetaInTable, ThetaEmTable, outType, fhs_simulated, ftb_simulated)
  case('first_order')
     call infslab_method_firstorder(mat, spc, nPosSamples, nSrcSamples, thetaInTable, ThetaEmTable, phiEmTable, outType, verbose, fhs_simulated, ftb_simulated)
  case('analytic')
     call infslab_method_analytic_mc(mat, spc, nPosSamples, thetaInTable, ThetaEmTable, outType, fhs_simulated, ftb_simulated)
  case default
     call utl_fatal_error("Error: Unknown simulation method: " // method)
  end select
  call cpu_time(cTime)

  tHours = floor((cTime - sTime) / 3600.0)
  tMins  = floor(((cTime - sTime) - tHours * 3600.0) / 60.0)
  tSecs  = (cTime - sTime) - tHours * 3600.0 - tMins * 60.0

  write(*,'("Cpu time: " ,(I4)," h ",(I2)," min ",(F6.3)," s")') int(tHours), int(tMins), tSecs

  call infslab_method_analytic(mat, spc, thetaInTable, ThetaEmTable, outType, fhs_analytic, ftb_analytic)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! WRITE OUTPUT
  !!

  call save_results(outType, mat, spc, resultFile, fhs_simulated, fhs_analytic, ftb_simulated, ftb_analytic)


contains


  !>
  !!
  !!
  subroutine save_results(outType, mat, spc, filename, hs_sim, hs_an, fDataSim, fDataAn)
    integer, intent(IN)                         :: outType
    type(gth_hemisphere), intent(IN), optional  :: hs_sim, hs_an
    real(fd), dimension(:,:,:), optional        :: fDataSim, fDataAn
    type(mat_material), intent(IN)              :: mat
    type(spc_spectrum), intent(IN)              :: spc
    character (len = *)                         :: filename

    integer            :: i

    integer            :: fileID
    integer            :: fluorID, fluorAnID

    integer            :: dimMatMuEnergy, dimFlrLine, dimSpecEnergy, dimSpecICDF, dimThetaIn, dimEmAngle
    integer            :: idMatMuIon, idMatMuFluor, idMatMuFluorCdf, idMatMuRay, idMatMuExt, idMatMuEnergy
    integer            :: idSpecEnergy, idSpecIntensity, idSpecCdf, idSpecIcdf
    integer            :: idAnalytic

    integer, parameter :: RTP = NF90_DOUBLE
    integer, parameter :: RTS = NF90_FLOAT

    if(outType == OUTPUT_TYPE_HEMISPHERE .and. (.not. present(hs_sim) .or. .not. present(hs_an))) then
       call utl_fatal_error("Incorrect data given in saveTabulated")
    end if

    if(outType == OUTPUT_TYPE_TABULATED .and. (.not. present(fDataSim) .or. .not. present(fDataAn))) then
       call utl_fatal_error("Incorrect data given in saveTabulated")
    end if


    call gth_nfCheck( nf90_create(filename, NF90_CLOBBER, fileID) )

    !!--- CORE ---
    !!
    !!- Core attributes
    !!
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "author",           "Hannu Parviainen" ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "program",          "xrfip"            ) )
    call gth_nfCheck( nf90_put_att(fileID, NF90_GLOBAL, "version",          "1.0"              ) )
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
       fluorAnID = gth_writeNetcdfHeader(hs_an, fileID)
    else
       call gth_nfCheck( nf90_def_var(fileID, "res__fld_table_simulated",   RTP, [dimThetaIn,dimEmAngle,dimFlrLine], fluorID   ) )
       call gth_nfCheck( nf90_def_var(fileID, "res__fld_table_analytical",  RTP, [dimThetaIn,dimEmAngle,dimFlrLine], fluorAnID ) )
    end if

    call gth_nfCheck( nf90_enddef(fileID) )


    !!--- WRITE VARIABLES ---
    !!
    if(outType == OUTPUT_TYPE_HEMISPHERE) then
       call gth_nfCheck( nf90_put_var(fileID, fluorID,   hs_sim%data ) )
       call gth_nfCheck( nf90_put_var(fileID, fluorAnID, hs_an%data  ) )
    else
       call gth_nfCheck( nf90_put_var(fileID, fluorID,   fDataSim    ) )
       call gth_nfCheck( nf90_put_var(fileID, fluorAnID, fDataAn     ) )
    end if

    call gth_nfCheck( nf90_put_var(fileID, idMatMuIon,       mat%muPhotoIon) )
    call gth_nfCheck( nf90_put_var(fileID, idMatMuFluor,     mat%muFluorLine) )
    call gth_nfCheck( nf90_put_var(fileID, idMatMuFluorCdf,  mat%muFluorLineCDF) )
    call gth_nfCheck( nf90_put_var(fileID, idMatMuRay,       mat%muRayleigh) )
    call gth_nfCheck( nf90_put_var(fileID, idMatMuExt,       mat%muExtTotal) )
    call gth_nfCheck( nf90_put_var(fileID, idMatMuEnergy,    mat%E) )

    call gth_nfCheck( nf90_put_var(fileID, idSpecEnergy,     spc%E) )
    call gth_nfCheck( nf90_put_var(fileID, idSpecIntensity,  spc%I) )
    if(srcType == 'Spectrum') then
       call gth_nfCheck( nf90_put_var(fileID, idSpecCdf,     spc%cdf) )
       call gth_nfCheck( nf90_put_var(fileID, idSpecIcdf,    spc%icdf) )
    end if


    call gth_nfCheck( nf90_close(fileID) )

  end subroutine save_results
  !******

end program xrfip
