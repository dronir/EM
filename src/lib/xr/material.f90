!--- BEGIN GPL --- 
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

!> Everything needed to compute soft X-ray fluorescence from a material
!! of arbitrary elementary composition. Contains material type and functions
!! to compute total and per-line photoionisation coefficients, as well as
!! everything else useful.
!!
!! \date  12.11.2007
module material
  use base
  use random
  use xraylib

  implicit none

  integer,  parameter :: NAMELENGTH = 256
  real(fd), parameter :: EPS_EDGE   = 1.e-2_fd


  !> Material
  !!
  type :: mat_material

     !!> Material name.
     !!
     character(LEN=NAMELENGTH) :: name

     !> Number of elements.
     !!
     integer     :: nElems
 
     !> Number of fluorescence lines.
     !!
     integer     :: nLines

     !!---------- ELEMENT PROPERTIES ------------------------------------------------------------
     !!
 
     !> Atomic number table.
     !!
     !! Table containing the atomic numbers of the elements used by the material.
     !!
     integer,  dimension(:), pointer    :: Z

     !> Number density fraction table.
     !! 
     !! Table containing the number density fractions of the elements used by the material.
     !! Normalised to sum up to 1  when the material is initialised.
     !!
     real(fd), dimension(:), pointer    :: fract


     !> Number density [1/cm^3] for each elements in the material.
     !!
     real(fd), dimension(:), pointer    :: elNumberDensity


     !> Atomic weight for each element in the material.
     !!
     !!
     real(fd), dimension(:), pointer    :: elAtomicWeight


     !> Mass density [g/cm^3] for each elements in the material.
     !!
     !! The mass density for each element is computed from the number density and atomic weights as
     !!
     !! 
     real(fd), dimension(:), pointer    :: elMassDensity


     !!---------- MATERIAL VARIABLES ------------------------------------------------------------
     !!

     !> Mass density [g/cm^3] for the material.
     !!
     real(fd) :: massDensity


     !> Number density [1/cm^3] for the material.
     !!
     real(fd) :: numberDensity


     !> Line energy [eV] table.
     !!
     !! Fluorescent line energy [eV] for all of the used fluorescence lines from all of the elements.
     !! \note Has two values per element, since currently only Ka and Kb lines are included.
     !! 
     real(fd), dimension(:),   pointer :: lEnergy


     !> Edge energy [eV] table.
     !!
     !! Fluorescence edge energy [eV] from all of the elements used in the material.
     !! \note Has one value per element, since currently only the K-shell ionization is considered.
     !!
     real(fd), dimension(:),   pointer :: eEnergy


     !> Energy [eV] table.
     !!
     real(fs), dimension(:), pointer :: E


     !> Material photoionization coefficient \f$\mu_{phot}\f$ [1/cm] as a function of energy.
     !!
     real(fd), dimension(:), pointer :: muPhotoIon


     !> Total fluorescence line photoionization coefficient \f$\mu_{abs,tot}\f$ [1/cm] as a function of energy.
     !!
     !! \note Should be very close to  \f$\mu_{phot}\f$ for most e.
     !!
     real(fd), dimension(:),   pointer :: muFluorLineTotal


     !> Fluorescence line photoionization coefficient \muPhoto \f$\mu_{ion,e,l}\f$ [1/cm] for each line as a function of energy.
     !!
     !! Fluorescence line photoionization coefficient \f$\mu_{ion,e,l}\f$ for the K shell is the product of 
     !! the K shell photoionisation cross-section, K shell fluorescence yield and line radiation rate, and 
     !! element number density, that is
     !! \f[ \mu_{ion,e,l} = \sigma_{ion,e,K} \; f_{e,K} \; r_{e,l} \; \rho_{N,e} \f]
     !!
     !! \note Should sum up to \f$\mu_{abs,tot}\f$ and \f$\mu_{phot}\f$ - higher level line absorbtion for each e. 
     !!
     real(fd), dimension(:,:), pointer :: muFluorLine


     !> Fluorescence line photoionization cdf for each line as a function of energy.
     !!
     !! \note Should sum up to 1 for each e.
     !!
     real(fd), dimension(:,:), pointer :: muFluorLineCdf


     !> Rayleigh scattering coefficient \f$\mu_{Rl}\f$ [1/cm] as a function of energy.
     !!
     real(fd), dimension(:),   pointer :: muRayleigh


     !> Total linear extinction coefficient \f$\mu_{tot}\f$ [1/cm] as a function of energy.
     !!
     real(fd), dimension(:),   pointer :: muExtTotal


     !> The length of the mean free path of a photon [cm] as a function of energy 
     !!
     real(fd), dimension(:),   pointer :: meanFreePath


     !> Fluorescence yield per shell.
     !!
     real(fd), dimension(:),   pointer :: fYield

 
     !> Minimum energy [eV] that can be handled by the material (= the energy of the lowest energy fluorescence line).
     !!
     real(fd) :: eMin


     !> Maximum energy [eV].
     !!
     real(fd) :: eMax


     !> Energy resolution [eV].
     !!
     real(fd) :: dE

     !> Inverse of the energy resolution.
     !!
     real(fd) :: dEInv

     !> size of the energy table.
     !!
     integer(il) :: nE


     !> Miminum edge energy [eV].
     !!
     !! Energies smaller than edgeEnergyMin cannot cause fluorescence in the material.
     !!
     real(fd) :: edgeEnergyMin

  end type mat_material

contains

  !> Initialises the material. Needs the material name, names of the elements, total molar density of the material,
  !! number density fractions of the elements, and the maximum energy to tabulate the different precomputed
  !! variables.
  !!
  !! Computes the total linear extinction coefficient  \f$ \mu_{tot} \f$, the total linear 
  !! absorption (ionization) coefficient \f$ \mu_{abs} \f$, and per-line absorption probabilities.
  !! 
  !! Linear absorbtion coefficient for an element e and a line l is computed as
  !!    \f[ \mu_l = N_e * \sigma_{e,l} = (\rho_e * N_a) / m_e * \sigma_{e,l} \quad [1/cm], \f]
  !! where \f$ N_e \f$ is the number density, \f$ \sigma_{e,l} \f$ is the ionization cross-section, 
  !! \f$ \rho_e \f$ is the mass density, \f$ N_a \f$ is the Avogardo constant [mol^-1], and \f$ m_e \f$ is the atomic 
  !! mass of the element [au], if \f$ \sigma_{e,l} \f$ is given in [barns/atom], or as
  !!    \f[ \mu_l = \rho_e * \sigma_{e,l} \quad [1/cm], \f]
  !! if \f$ \sigma_{e,l} \f$ is given in [cm^2/g] and \f$ \rho_e \f$ in [g/cm^3]. In this code we use
  !! the latter, since the xraylib-package has functions to obtain the cross-sections in [cm^2/g].
  !!
  !! \param[out] mat          Material
  !! \param[in]  matName      Material name
  !! \param[in]  eNames       Element names
  !! \param[in]  numberFracs  Element number/mole density fractions
  !! \param[in]  molarDensity Total molar density [mol/cm^3] of the material
  !! \param[in]  maxEnergy    Maximum energy [eV]
  !!
  subroutine mat_material_init(mat, matName, eNames, molarDensity, numberFracs, maxEnergy)

    type(mat_material)      :: mat
    character(LEN=*)        :: matName
    character(len=2)        :: eNames(:)
    real(fd)                :: numberFracs(:)
    real(fd)                :: molarDensity
    real(fd)                :: maxEnergy

    character(len=2)        :: eNameTbl(size(eNames))
    integer                 :: Z(size(eNames))

    real(fs)                :: energy

    integer(is) :: i, j

    eNameTbl = eNames

    !! Force element fraction sum to unity
    !!
    numberFracs = numberFracs / sum(numberFracs)

    mat % name   = matName
    mat % nElems = size(Z)
    mat % nLines = mat % nElems * 2
 
    do i = 1, mat%nElems
       call mat_elemInit(eNames(i), Z(i))
    end do

    allocate( mat % eEnergy (mat%nElems) )
    allocate( mat % lEnergy (mat%nLines) )
    allocate( mat % Z       (mat%nElems) )

    mat % Z      = Z

    !! Compute the edge and line energies for the elements.
    !!
    mat%edgeEnergyMin = 1e10_fd
    mat%eMin          = 1e10_fd
    do i = 1, mat%nElems
       mat%lEnergy(i*2-1) = LineEnergy(Z(i),KA_LINE) * keV_eV
       mat%lEnergy(i*2)   = LineEnergy(Z(i),KB_LINE) * keV_eV

       if(sum(mat%lEnergy(i*2-1:i*2)) > 0.0) then
          mat%eEnergy(i)     = EdgeEnergy(Z(i),K_SHELL) * keV_eV
       end if
          
       if(mat%eEnergy(i)     < mat%edgeEnergyMin .and. mat%eEnergy(i    ) > EPS_EDGE) mat%edgeEnergyMin = mat%eEnergy(i)
       if(mat%lEnergy(i*2-1) < mat%eMin          .and. mat%lEnergy(i*2-1) > EPS_EDGE) mat%eMin          = mat%lEnergy(i*2-1)
       if(mat%lEnergy(i*2)   < mat%eMin          .and. mat%lEnergy(i*2  ) > EPS_EDGE) mat%eMin          = mat%lEnergy(i*2)
    end do
    

    !! We optimize the size of the precomputed tables by setting the
    !! minimum energy to the minimum line energy. We don't need to
    !! consider energies smaller than this when computing fluorescence.
    !! If scattering is included to the simulation, this should be
    !! changed.
    !!

    if(mat%eMin > mat%edgeEnergyMin) mat%eMin = mat%edgeEnergyMin

    mat % eMin = floor(mat%eMin)
    mat % eMax  = maxEnergy
   
    mat % dE    = 10.0_fd
    mat % dEInv = 1.0_fd / mat % dE
    
    mat % nE    = ceiling((mat % eMax - mat % eMin) / mat % dE) + 1
    
    allocate(mat%muPhotoIon(mat%nE), mat%muFluorLine(mat%nLines, mat%nE), mat%muFluorLineCdf(mat%nLines, mat%nE))
    allocate(mat%muExtTotal(mat%nE), mat%muFluorLineTotal(mat%nE),  mat%muRayleigh(mat%nE))
    allocate(mat%meanFreePath(mat%nE))
    allocate(mat%elNumberDensity(mat%nElems))
    allocate(mat%fYield(mat%nElems))
    
    allocate(mat%E(mat%nE))

    !! Precompute the number density N_e [1/cm^3] and fluorescence yields for each element.
    !!
    !! PHS_NA [1/mol]
    !!
    do i = 1, mat % nElems
       mat % elNumberDensity = numberFracs(i) * molarDensity * PHS_NA  !! [1/cm^3] = [1] [mol/cm^3] [1/mol]
       mat % fYield(i)       = FluorYield(Z(i), K_SHELL)
    end do

    !! Compute the linear absorption coefficients for each line as a function of energy, and 
    !! the total extinction coefficient as a function of energy.
    !!

    mat % muPhotoIon        = 0.0_fd
    mat % muExtTotal        = 0.0_fd
    mat % muFluorLineTotal  = 0.0_fd
    mat % muRayleigh        = 0.0_fd

    do j = 1, mat % nE
       mat % E(j) = (mat % eMin + (j-1) * mat % dE) 
       energy     = real(mat % E(j) * eV_keV, fs)
      do i = 1, mat % nElems

          if(mat%lEnergy(i*2-1) > 0.0) then
             mat % muFluorLine(i*2-1,j) = CSb_FluorLine_Kissel(Z(i), KA_LINE, energy ) * barn_to_cm2 &
                  & * mat % elNumberDensity(i) / mat%fYield(i)
          end if

          if(mat%lEnergy(i*2) > 0.0) then
             mat % muFluorLine(i*2  ,j) = CSb_FluorLine_Kissel(Z(i), KB_LINE, energy ) * barn_to_cm2  &
                  & * mat % elNumberDensity(i)  / mat%fYield(i)
          end if

          mat % muPhotoIon(j) = mat % muPhotoIon(j) + CSb_Photo ( Z(i), energy ) * barn_to_cm2 * mat%elNumberDensity(i)
          mat % muRayleigh(j) = mat % muRayleigh(j) + CSb_Rayl  ( Z(i), energy ) * barn_to_cm2 * mat%elNumberDensity(i)
          mat % muExtTotal(j) = mat % muExtTotal(j) + CSb_Total ( Z(i), energy ) * barn_to_cm2 * mat%elNumberDensity(i)

       end do  

       mat % meanFreePath(j) = (1.0 / mat%muExtTotal(i))

    end do

    !! Compute the total absorption coefficient as a function of energy, and the line absorption cdf.
    !!
    do j = 1, mat % nE
       mat % muFluorLineTotal(j) = sum(mat % muFluorLine(:,j))
       
       if(mat % muFluorLineTotal(j) > 0.0_fd) then
          mat % muFluorLineCdf(1,j) = mat % muFluorLine(1,j) / mat % muFluorLineTotal(j)
          do i = 2, mat % nLines
             mat % muFluorLineCdf(i,j) = mat % muFluorLineCdf(i-1,j) + mat % muFluorLine(i,j) / mat % muFluorLineTotal(j)
          end do
          
       else
          mat % muFluorLineCdf(:,j) = 0.0_fd
       end if

    end do

  end subroutine mat_material_init

  subroutine mat_evalMu(m, energy, muFluorLineTotal, muPhotoIon, muTotal, iEn, xEn)
    type(mat_material) :: m
    real(fd)      :: energy
    real(fd), optional, intent(out) :: muFluorLineTotal, muPhotoIon, muTotal
    real(fd), optional :: xEn
    integer, optional  :: iEn

    integer(il)   :: i
    real(fd)      :: x

    if(energy >= m%eMin .and. energy <= m%eMax) then

       x = (energy - m % eMin) * m % dEInv
       i = floor(x)
       x = x - i
       i = i + 1

       if( present(muFluorLineTotal) ) muFluorLineTotal = (1.0_fd - x) * m%muFluorLineTotal(i) +  x * m%muFluorLineTotal(i+1)
       if( present(muPhotoIon)       ) muPhotoion       = (1.0_fd - x) * m%muPhotoion(i)       +  x * m%muPhotoIon(i+1)
       if( present(muTotal)          ) muTotal          = (1.0_fd - x) * m%muExtTotal(i)       +  x * m%muExtTotal(i+1)

       if(present(iEn)) iEn = i
       if(present(xEn)) xEn = x

    else
       print *, energy, m%eMax, m%emin
       call utl_fatal_error("Fatal error evaluating mat_evalMu: energy index out of bounds.")
    end if

  end subroutine  mat_evalMu

!!$  subroutine mat_evalFluorLineFractions(m, energy, lineFractions)
!!$    type(mat_material)     :: m
!!$    real(fd)               :: energy
!!$    real(fd), dimension(:) :: lineFractions
!!$
!!$    integer(il)   :: i
!!$    real(fd)      :: x
!!$
!!$    if(size(lineFractions) == m%nLines) then
!!$       if(energy >= m%eMin .AND. energy <= m%eMax) then
!!$
!!$          x = (energy - m % eMin) * m % dEInv
!!$          i = floor(x)
!!$          x = x - i
!!$          i = i + 1
!!$
!!$          lineFractions  = ((1.0_fd - x) * m%muFluorLine(:,i) + x * m%muFluorLine(:,i+1))
!!$          lineFractions  = lineFractions / sum(lineFractions)
!!$
!!$       else
!!$          call utl_fatal_error("Fatal error evaluating fluorescence line fractions: energy index out of bounds.")
!!$       end if
!!$    else
!!$       call utl_fatal_error("Fatal error evaluating fluorescence line fractions: Wrong array size.")
!!$    end if
!!$  END subroutine  mat_evalFluorLineFractions

  function mat_evalFluorLineFractions(m, energy) result(lineFractions)
    type(mat_material) :: m
    real(fd)           :: energy
    real(fd)           :: lineFractions(m%nLines)

    integer(il)   :: i
    real(fd)      :: x

    if(energy >= m%eMin .and. energy <= m%eMax) then

       x = (energy - m % eMin) * m % dEInv
       i = floor(x)
       x = x - i
       i = i + 1

       lineFractions  = ((1.0_fd - x) * m%muFluorLine(:,i) + x * m%muFluorLine(:,i+1))
       lineFractions  = lineFractions / sum(lineFractions)

    else
       print *, energy, m%eMax, m%emin
       call utl_fatal_error("Fatal error evaluating fluorescence line fractions: energy index out of bounds.")
    end if

  end function mat_evalFluorLineFractions


  real(fd) function mat_evalMuAbs(m, energy, iEn, xEn)
    type(mat_material) :: m
    real(fd)           :: energy
    integer, optional  :: iEn
    real(fd), optional :: xEn

    integer(il)   :: i
    real(fd)      :: x

    if(energy >= m%eMin .and. energy <= m%eMax) then

       x = (energy - m % eMin) * m % dEInv
       i = floor(x)
       x = x - i
       i = i + 1

       mat_evalMuAbs = (1.0_fd - x) * m%muFluorLineTotal(i) +  x * m%muFluorLineTotal(i+1)

       if(present(iEn)) iEn = i
       if(present(xEn)) xEn = x

    else
       print *, energy, m%eMin, m%eMax
       call utl_fatal_error("Fatal error evaluating mat_evalMuAbs: energy index out of bounds.")
    end if

  end function mat_evalMuAbs


  real(fd) function mat_evalMuExt(m, energy)
    type(mat_material) :: m
    real(fd)      :: energy

    integer(il)   :: i
    real(fd)      :: x

    if(energy >= m%eMin .and. energy <= m%eMax) then

       x = (energy - m % eMin) * m % dEInv
       i = floor(x)
       x = x - i
       i = i + 1

       mat_evalMuExt = (1.0_fd - x) * m%muExtTotal(i) +  x * m%muExtTotal(i+1)

    else
       call utl_fatal_error("Fatal error evaluating mat_evalMuExt: energy index out of bounds.")
    end if

  end function mat_evalMuExt

  subroutine mat_evalFluorescence(m, energy, fEnergy, fYield, lineIdx)
    type(mat_material), intent(IN) :: m
    real(fd),           intent(IN) :: energy
    
    real(fd), intent(OUT)          :: fEnergy, fYield
    integer, intent(OUT)           :: lineIdx

    real(fd) :: x(1)
    integer  :: i

    call rnd_generate_uniform_n(x)   

    i = 1
    do while (m%muFluorLineCDF(i, floor((energy-m%eMin)*m%dEInv)+1 ) < x(1) .and. i <= m%nLines) 
       i = i+1
    end do

    if(i <= m%nLines ) then

       fEnergy = m % lEnergy(i)
       fYield  = m % fYield((i+1)/2)
       lineIdx = i

    else
       call utl_fatal_error("Fatal error evaluating fluorescence: The ray is below the lowest edge energy.")
    end if
       

  end subroutine mat_evalFluorescence

  subroutine mat_elemInit(eName, Z)
    character(len=2) :: eName

    integer  :: Z
    logical  :: eFound = .false.

    character(len=2), dimension(54), parameter :: eNameTable = [ &
         & 'H ', 'He', & 
         & 'Li', 'Be',  'B ',  'C ',  'N ',  'O ',  'F ', 'Ne', & 
         & 'Na', 'Mg', 'Al', 'Si',  'P ',  'S ', 'Cl', 'Ar', & 
         &  'K ', 'Ca', 'Sc', 'Ti',  'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', & 
         & 'Rb', 'Sr',  'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',  'I ', 'Xe'  &
         & ]

    Z = 1
    do Z = 1, size(eNameTable)
       if (eNameTable(Z) == eName) then
          eFound = .true.
          exit
       end if
    end do

    if(.not. eFound) then
       call utl_fatal_error("Fatal Error: Element " // eName // " was not found.")
    end if

  end subroutine mat_elemInit

end module material


