module xrFileIO
implicit none

use hemisphere
use material
use spectrum
use netcdf

contains
  subroutine save_xr(filename, author, creator, version, dataHemisphere, anDataHemisphere, dataGrid, material, spectrum)
    character(len=*)               :: filename
    character(len=*)               :: author
    character(len=*)               :: creator
    character(len=*)               :: version
    type(mat_material)             :: material
    type(spc_spectrum)             :: spectrum

    type(gth_hemisphere), optional :: dataHemisphere
    type(gth_hemisphere), optional :: anDataHemisphere
    real(fd),             optional :: dataGrid(:,:)

    type(gth_hsFile)   :: hf
    integer            :: i

    integer            :: dimMatMuEnergy, dimFlrLine, dimSpecEnergy, dimSpecICDF
    integer            :: idMatMuIon, idMatMuFluor, idMatMuFluorCdf, idMatMuRay, idMatMuExt, idMatMuEnergy
    integer            :: idSpecEnergy, idSpecIntensity, idSpecCdf, idSpecIcdf
    integer            :: idAnalytic

    integer, parameter :: RTP = NF90_DOUBLE
    integer, parameter :: RTS = NF90_FLOAT

    call gth_hsFileOpen(hf, fName, "w")

    if(present(dataHemisphere)) then
       call gth_hsWriteHeader(hf, dataHemisphere, author, creator, version)
    end if

    if(present(anDataHemisphere)) then
       Call gth_nfCheck( nf90_def_var(hf%fileID, "Hemisphere_analytic", NF90_DOUBLE, hf%dimHs, idAnalytic) )
       Call gth_nfCheck( nf90_put_att(hf%fileID,  idAnalytic, "long_name",  "Analytic approximation" ) )
    end if

    !! DEFINE DIMENSIONS
    !!
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Fluorescence_line",     material%nLines,           dimFlrLine     ) )
    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Material_energy",       material%nE,               dimMatMuEnergy ) )
    

    Call gth_nfCheck( nf90_def_dim(hf%fileID, "Spectrum_energy",       spectrum%nPoints,          dimSpecEnergy  ) )
    if(spectrum%nPoints > 1) then
       Call gth_nfCheck( nf90_def_dim(hf%fileID, "Spectrum_ICDF_energy",  spectrum%nPoints*5,      dimSpecICDF    ) )
    end if


    !! DEFINE VARIABLES
    !!

    Call gth_nfCheck( nf90_def_var(hf%fileID, "Photoionization_coefficient",   RTP, [dimMatMuEnergy],            idMatMuIon) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Fluorescence_line_coefficient", RTP, [dimFlrLine,dimMatMuEnergy], idMatMuFluor) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Fluorescence_line_cdf",         RTP, [dimFlrLine,dimMatMuEnergy], idMatMuFluorCdf))
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Rayleigh_coefficient",          RTP, [dimMatMuEnergy],            idMatMuRay) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Extinction_coefficient",        RTP, [dimMatMuEnergy],            idMatMuExt) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Material_energy",               RTP, [dimMatMuEnergy],            idMatMuEnergy) )

    Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_energy",            RTP, [dimSpecEnergy], idSpecEnergy) )
    Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_intensity",         RTP, [dimSpecEnergy], idSpecIntensity) )
    if(spectrum%nPoints > 1) then
       Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_cdf",            RTP, [dimSpecEnergy], idSpecCdf) )
       Call gth_nfCheck( nf90_def_var(hf%fileID, "Spectrum_inverse_cdf",    RTP, [dimSpecICDF],   idSpecIcdf) )
    end if
 

    !! WRITE GLOBAL ATTRIBUTES
    !!
    call gth_hsAddAttS(hf, NF90_GLOBAL, "Elements",           elements)
    call gth_hsAddAttS(hf, NF90_GLOBAL, "Spectrum_type",      sourceType)
    call gth_hsAddAttS(hf, NF90_GLOBAL, "Simulation_method",  method)
    call gth_hsAddAttF(hf, NF90_GLOBAL, "Theta_in",           thetaInTable)

    call gth_hsAddAttF(hf, NF90_GLOBAL, "Fluorescence_line_energy", material%lEnergy)
    call gth_hsAddAttF(hf, NF90_GLOBAL, "Absorbtion_edge_energy",   material%eEnergy)

    call gth_hsSaveData(h, hf)
    Call gth_nfCheck( nf90_put_var(hf%fileID, idAnalytic, hs_analytic%data) )

    !! WRITE VARIABLES
    !!
    
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuIon,       material%muPhotoIon) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuFluor,     material%muFluorLine) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuFluorCdf,  material%muFluorLineCDF) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuRay,       material%muRayleigh) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuExt,       material%muExtTotal) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idMatMuEnergy,    material%E) )
    
    Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecEnergy,     spectrum%E) )
    Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecIntensity,  spectrum%I) )
    if(spectrum%nPoints > 1) then
       Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecCdf,     spectrum%cdf) )
       Call gth_nfCheck( nf90_put_var(hf%fileID, idSpecIcdf,    spectrum%icdf) )
    end if
 
    call gth_hsFileClose(hf)

  end subroutine saveHemisphere

end module xrFileIO
