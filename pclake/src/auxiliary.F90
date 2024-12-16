#include "fabm_driver.h"
!-----------------------------------------------------------------------
!
!BOP
!
! !INTERFACE:
   module pclake_auxiliary
! !USES:
   use fabm_types
   use fabm_expressions
   use pclake_utility, ONLY: uFunTmAbio,uFunTmBio

   implicit none

!  default: all is private.
      private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_auxiliary
!  diagnostic variables for local output
!  state dependencies identifiers
!  SW: Sediment to Water
!  dependencies to pclake_abiotic_water state variables
   type (type_state_variable_id)            :: id_SWNH4,id_SWNO3,id_SWPO4,id_SWPAIM,id_SWO2,id_SWSiO2
   type (type_state_variable_id)            :: id_SWDIM,id_SWDPOM,id_SWNPOM,id_SWPPOM,id_SWSiPOM
!  dependencies to phytoplankton_water state variables
   type (type_state_variable_id)            :: id_SWDBlue,id_SWNBlue,id_SWPBlue
   type (type_state_variable_id)            :: id_SWDDiat,id_SWNDiat,id_SWPDiat  !,id_SWSiDiat
   type (type_state_variable_id)            :: id_SWDGren,id_SWNGren,id_SWPGren
!  WS: Water to Sediment
!  dependencies to pclake_abiotic_sediment state variables
   type (type_bottom_state_variable_id)     :: id_WSPO4,id_WSPAIM,id_WSNH4,id_WSNO3
   type (type_bottom_state_variable_id)     :: id_WSDIM,id_WSDPOM,id_WSNPOM,id_WSPPOM,id_WSSiPOM
   type (type_bottom_state_variable_id)     :: id_WSDHum,id_WSNHum,id_WSPHum
!  dependencies to phytoplankton_sediment state variables
   type (type_bottom_state_variable_id)     :: id_WSDBlue,id_WSNBlue,id_WSPBlue
   type (type_bottom_state_variable_id)     :: id_WSDDiat,id_WSNDiat,id_WSPDiat  !,id_WSSiDiat
   type (type_bottom_state_variable_id)     :: id_WSDGren,id_WSNGren,id_WSPGren
!  dependencies to zooplankton
   type (type_state_variable_id)            :: id_DTranZoo,id_NTranZoo,id_PTranZoo
   type (type_bottom_state_variable_id)     :: id_TurbFish
!  dependencies to macrophytes state variables
   type (type_horizontal_dependency_id)     :: id_DragVeg
!   type (type_bottom_state_variable_id)     :: id_DragVeg
!  environmental dependencies
   type (type_global_dependency_id)         :: id_Day
   type (type_dependency_id)                :: id_uTm ,id_dz
   type (type_horizontal_dependency_id)     :: id_sDepthW,id_shear
   type (type_horizontal_diagnostic_variable_id)       :: id_tDBurIM,id_shearstress
   type (type_horizontal_diagnostic_variable_id)       :: id_aFunDimSusp,id_aFunTauSet,id_tDResusDead
#ifdef _DEVELOPMENT_
!  diagnostic variables for resuspension fluxes
   type (type_horizontal_diagnostic_variable_id)  :: id_tAuxDIMW,id_tDAuxPOMW,id_tNAuxPOMW
   type (type_horizontal_diagnostic_variable_id)  :: id_tPAuxPOMW,id_tSiAuxPOMW,id_tAuxPAIMW
   type (type_horizontal_diagnostic_variable_id)  :: id_tNAuxNH4W,id_tNAuxNO3W,id_tPAuxPO4W
   type (type_horizontal_diagnostic_variable_id)  :: id_tDAxuDiatW,id_tNAuxDiatW,id_tPAuxDiatW
   type (type_horizontal_diagnostic_variable_id)  :: id_tDAuxGrenW,id_tNAuxGrenW,id_tPAuxGrenW
   type (type_horizontal_diagnostic_variable_id)  :: id_tDAuxBlueW,id_tNAuxBlueW,id_tPAuxBlueW
   type (type_horizontal_diagnostic_variable_id)  :: id_tAuxDIMS,id_tDAuxPOMS,id_tNAuxPOMS
   type (type_horizontal_diagnostic_variable_id)  :: id_tPAuxPOMS,id_tSiAuxPOMS,id_tPAuxPO4S
   type (type_horizontal_diagnostic_variable_id)  :: id_tAuxPAIMS,id_tNAuxNH4S,id_tNAuxNO3S
   type (type_horizontal_diagnostic_variable_id)  :: id_tDAuxHumS,id_tPAuxHumS,id_tNAuxHumS
   type (type_horizontal_diagnostic_variable_id)  :: id_tDAxuDiatS,id_tNAuxDiatS,id_tPAuxDiatS
   type (type_horizontal_diagnostic_variable_id)  :: id_tDAuxGrenS,id_tNAuxGrenS,id_tPAuxGrenS
   type (type_horizontal_diagnostic_variable_id)  :: id_tDAuxBlueS,id_tNAuxBlueS,id_tPAuxBlueS
#endif
!  diagnostic dependencies,due to burial process
   type ( type_horizontal_dependency_id)             :: id_tDAbioHumS
   type ( type_horizontal_dependency_id)             :: id_tDAbioPOMS,id_tDPrimPOMS,id_tDWebPOMS,id_tDBedPOMS
!  Model parameters
!  logical variables for choosing linking dependencies
!  diagnostic dependencies, due to resuspension
   logical      :: fish_module,macrophytes_module
   real(rk)                   :: cDepthS
!  sediment properties parameters
   real(rk)                   :: fLutum,fLutumRef,bPorS
!  resuspension process parameters
   real(rk)                   :: kVegResus,kTurbFish
   real(rk)                   :: cSuspRef,cSuspMin,cSuspMax,cSuspSlope
   real(rk)                   :: hDepthSusp,cFetchRef,cFetch
   real(rk)                   :: cResusPhytExp,kResusPhytMax  !,cSiDDiat
!  sedimentation parameters
   real(rk)                   :: cThetaSed,cVSedIM,cVSedPOM
   real(rk)                   :: cVSedDiat,cVSedGren,cVSedBlue
!  Burial process parameters
   real(rk)                   :: cRhoIM,cRhoOM,fDOrgSoil
   real(rk)                   :: cPO4Ground,cNH4Ground,cNO3Ground
!  parameter for fish temperature function
   real(rk)                   :: cSigTmFish,cTmOptFish
!  nutrient ratios parameter
   real(rk)   :: cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin
   real(rk)   :: cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax
!  parameters relating to new resuspension method
   real(rk)   :: crt_shear,ref_shear,alpha,eta,cVSedMain
!  variable for selecting different resuspension methods
   integer    :: resusp_meth
!  variables for atmospheric depostions
   real(rk)                   :: tDDepoIM,tDDepoPOM,tNDepoPOM,tPDepoPO4
   REAL(rk)                   :: tPDepoPOM,tNDepoNH4,tNDepoNO3

   contains

!  Model procedures
   procedure ::  initialize
   procedure ::  do_bottom
   PROCEDURE ::  do_surface
   end type type_pclake_auxiliary

!  private data members (API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
   real(rk),parameter :: Pi=3.14159265358979_rk
   real(rk), parameter :: mmPerm=1000.0_rk
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_pclake_auxiliary), intent(inout), target :: self
   integer,                     intent(in)            :: configunit

!  Store parameter values in derived type
!  NB: all rates must be provided in values per day in .yaml input,
!  and are converted here to values per second.
   call self%get_parameter(self%cDepthS,      'cDepthS',      'm',                 'sediment depth',                                         default=0.1_rk)
   call self%get_parameter(self%cThetaSed,    'cThetaSed',    '1/e^oC',            'temperature parameter for sedimentation',                 default=1.01_rk)
   call self%get_parameter(self%kVegResus,    'kVegResus',    'm2/gDW',            'relative resuspension reduction per gram macrophytes',    default=0.01_rk)
   call self%get_parameter(self%kTurbFish,    'kTurbFish',    'g/gfish/d',         'relative resuspension by adult fish browsing',           default=1.0_rk)
!  don't convert here,due to PCLake empirical function. kTurbFish
   call self%get_parameter(self%cSuspRef,     'cSuspRef',     '[-]',               'reference suspended matter function [-]',                default=0.0_rk)
   call self%get_parameter(self%cSuspMin,     'cSuspMin',     '[-]',               'minimum value of logistic function',                     default=6.1_rk)
   call self%get_parameter(self%cSuspMax,     'cSuspMax',     '[-]',               'maximum value of logistic function',                     default=25.2_rk)
   call self%get_parameter(self%cSuspSlope,   'cSuspSlope',   '[-]',               'slope of logistic function',                             default=2.1_rk)
   call self%get_parameter(self%hDepthSusp,   'hDepthSusp',   '[-]',               'half-sat. value of depth in logistic function',          default=2.0_rk)
   call self%get_parameter(self%cFetchRef,    'cFetchRef',    'm',                 'reference fetch',                                        default=1000.0_rk)
   call self%get_parameter(self%cFetch,       'cFetch',       'm',                 'length of water body in prevailing wind direction',            default=1000.0_rk)
   call self%get_parameter(self%fLutum ,      'fLutum',       '[-]',               'lutum content of inorganic matter',                         default=0.1_rk)
   call self%get_parameter(self%fLutumRef,    'fLutumRef',    '[-]',               'reference lutum fraction',                               default=0.2_rk)
   call self%get_parameter(self%bPorS,        'bPorS',        'm3water/m3sediment','sediment porosity',                                      default=0.847947_rk)
   call self%get_parameter(self%cResusPhytExp,'cResusPhytExp','(gDW m-2 d-1)-1',   'exponent parameter for phytoplankton resuspension',      default=-0.379_rk)
   call self%get_parameter(self%kResusPhytMax,'kResusPhytMax','d-1',               'maximum phytoplankton resuspension rate',                   default=0.25_rk)
   call self%get_parameter(self%cVSedIM,      'cVSedIM',      'm d-1',             'maximum sedimentation velocity of inert organic matter',    default=-1.0_rk)
   call self%get_parameter(self%cVSedPOM,     'cVSedPOM',     'm d-1',             'maximum sedimentation velocity of POM',                  default=0.25_rk)
   call self%get_parameter(self%cVSedDiat,    'cVSedDiat',    'm d-1',             'sedimentation velocity of diatoms',                         default=0.5_rk)
   call self%get_parameter(self%cVSedGren,    'cVSedGren',    'm d-1',             'sedimentation velocity of greens',                       default=0.2_rk)
   call self%get_parameter(self%cVSedBlue ,   'cVSedBlue',    'm d-1',             'sedimentation velocity of blue-greens',                     default=0.06_rk)
   call self%get_parameter(self%cRhoIM,       'cRhoIM',       'g m-3',             'density of sediment IM',                                 default=2500000.0_rk)
   call self%get_parameter(self%cRhoOM,       'cRhoOM',       'g m-3',             'density of sediment organic matter',                           default=1400000.0_rk)
   call self%get_parameter(self%fDOrgSoil,    'fDOrgSoil',    '[-]',               'fraction soil organic matter',                           default=0.1_rk)
   call self%get_parameter(self%cPO4Ground,   'cPO4Ground',   'mgP l-1',           'PO4 concentration in groundwater',                                default=0.1_rk)
   call self%get_parameter(self%cNH4Ground,   'cNH4Ground',   'mgN l-1',           'NH4 concentration in groundwater',                                default=1.0_rk)
   call self%get_parameter(self%cNO3Ground,   'cNO3Ground',   'mgN l-1',           'NO3 concentration in groundwater',                                default=0.1_rk)
   call self%get_parameter(self%cTmOptFish,   'cTmOptFish',   'degree C',          'optimal temperature for fish',                            default=25.0_rk)
   call self%get_parameter(self%cSigTmFish,   'cSigTmFish',   'degree C',          'temperature constant for fish (sigma)',                   default=10.0_rk)
   call self%get_parameter(self%cNDDiatMin,   'cNDDiatMin',   'mgN/mgDW',          'minimum N/DW ratio for diatoms',                            default=0.01_rk)
   call self%get_parameter(self%cPDDiatMin,   'cPDDiatMin',   'mgP/mgDW',          'minimum P/DW ratio for diatoms',                            default=0.0005_rk)
   call self%get_parameter(self%cNDGrenMin,   'cNDGrenMin',   'mgN/mgDW',          'minimum N/DW ratio for greens',                             default=0.02_rk)
   call self%get_parameter(self%cPDGrenMin,   'cPDGrenMin',   'mgP/mgDW',          'minimum P/DW ratio for greens',                             default=0.0015_rk)
   call self%get_parameter(self%cNDBlueMin,   'cNDBlueMin',   'mgN/mgDW',          'minimum N/DW ratio for blue-greens',                         default=0.03_rk)
   call self%get_parameter(self%cPDBlueMin,   'cPDBlueMin',   'mgP/mgDW',          'minimum P/DW ratio for blue-greens',                         default=0.0025_rk)
   call self%get_parameter(self%cNDBlueMax,   'cNDBlueMax',   'mgN/mgDW',          'maximum N/DW ratio for blue-greens',                         default=0.15_rk)
   call self%get_parameter(self%cNDDiatMax,   'cNDDiatMax',   'mgN/mgDW',          'maximum N/DW ratio for diatoms',                            default=0.05_rk)
   call self%get_parameter(self%cNDGrenMax,   'cNDGrenMax',   'mgN/mgDW',          'maximum N/DW ratio for greens',                             default=0.1_rk)
   call self%get_parameter(self%cPDBlueMax,   'cPDBlueMax',   'mgP/mgDW',          'maximum P/DW ratio for blue-greens',                        default=0.025_rk)
   call self%get_parameter(self%cPDDiatMax,   'cPDDiatMax',   'mgP/mgDW',          'maximum P/DW ratio for diatoms',                            default=0.005_rk)
   call self%get_parameter(self%cPDGrenMax,   'cPDGrenMax',   'mgP/mgDW',          'maximum P/DW ratio for greens',                             default=0.015_rk)
!  resuspension related to shear stress
   call self%get_parameter(self%crt_shear,    'crt_shear',    'N m-2',             'critical shear stress',                                  default=0.005_rk)
   call self%get_parameter(self%ref_shear,    'ref_shear',    'N m-2',             'reference shear stress',                                 default=1.0_rk)
   call self%get_parameter(self%alpha,        'alpha',        'g m-2 d-1',         'gross rate of sediment erosion (9000-9900g m-2/d)',      default=9000.0_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%eta,          'eta',          '[-]',               'shear stress correction factor',                         default=1.0_rk)
   call self%get_parameter(self%cVSedMain,    'cVSedMain',    'm d-1',             'depth averaged settling velocity (between 0.5-1.5m/d)',  default=0.5_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%resusp_meth,  'resusp_meth',  '[-]',               '1=original PCLake resuspension function',                default=2)
   call self%get_parameter(self%tDDepoIM,     'tDDepoIM',     'g m-2 d-1',         'inorganic matter deposition',                            default=0.0_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%tDDepoPOM,    'tDDepoPOM',    'g m-2 d-1',         'organic matter deposition, dry weight',                   default=0.0_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%tNDepoPOM,    'tNDepoPOM',    'g m-2 d-1',         'organic matter deposition, nitrogen',                     default=0.0_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%tPDepoPOM,    'tPDepoPOM',    'g m-2 d-1',         'organic matter deposition, phosphorus',                   default=0.0_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%tPDepoPO4,    'tPDepoPO4',    'g m-2 d-1',         'phosphate deposition',                                   default=0.0_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%tNDepoNH4,    'tNDepoNH4',    'g m-2 d-1',         'ammonium deposition',                                    default=0.0_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%tNDepoNO3,    'tNDepoNO3',    'g m-2 d-1',         'nitrate deposition',                                     default=0.0_rk, scale_factor=1.0_rk/secs_pr_day)
!  Register dependencies to abiotic water module
   call self%register_state_dependency(self%id_SWNH4,   'ammonium_pool_in_water',         'g m-3',  'ammonium pool in water')
   call self%register_state_dependency(self%id_SWNO3,   'nitrate_pool_in_water',          'g m-3',  'nitrate pool in water')
   call self%register_state_dependency(self%id_SWPO4,   'phosphate_pool_in_water',        'g m-3',  'phosphate pool in water')
   call self%register_state_dependency(self%id_SWPAIM,  'adsorbed_phosphorus_in_water'  , 'g m-3',  'adsorbed phosphorus in water')
   call self%register_state_dependency(self%id_SWO2,    'oxygen_pool_in_water',           'g m-3',  'oxygen pool in water')
   call self%register_state_dependency(self%id_SWDIM,   'inorg_pool_in_water',            'g m-3',  'inorganic matter pool in water')
   call self%register_state_dependency(self%id_SWDPOM,  'POM_DW_in_water',                'g m-3',  'POM DW in water')
   call self%register_state_dependency(self%id_SWNPOM,  'POM_N_in_water',                 'g m-3',  'POM N in water')
   call self%register_state_dependency(self%id_SWPPOM,  'POM_P_in_water',                 'g m-3',  'POM P in water')
   call self%register_state_dependency(self%id_SWSiPOM, 'POM_Si_in_water',                'g m-3',  'POM Si in water')
   call self%register_state_dependency(self%id_SWSiO2,  'SiO2_pool_water',                'g m-3',  'SiO2 pool water')
!  Register dependencies to abiotic sediment module
   call self%register_state_dependency(self%id_WSNH4,   'ammonium_pool_in_sediment',      'g m-2', 'ammonium pool in sediment')
   call self%register_state_dependency(self%id_WSNO3,   'nitrate_pool_in_sediment',       'g m-2', 'nitrate pool in sediment')
   call self%register_state_dependency(self%id_WSPO4,   'phosphate_pool_in_sediment',     'g m-2', 'phosphate pool in sediment')
   call self%register_state_dependency(self%id_WSPAIM,  'adsorbed_phosphorus_in_sediment','g m-2', 'adsorbed phosphorus in sediment')
   call self%register_state_dependency(self%id_WSDIM,   'inorg_pool_in_sediment',         'g m-2', 'inorganic matter pool in sediment')
   call self%register_state_dependency(self%id_WSDPOM,  'POM_DW_in_sediment',             'g m-2', 'POM DW in sediment')
   call self%register_state_dependency(self%id_WSNPOM,  'POM_N_in_sediment',              'g m-2', 'POM N in sediment')
   call self%register_state_dependency(self%id_WSPPOM,  'POM_P_in_sediment',              'g m-2', 'POM P in sediment')
   call self%register_state_dependency(self%id_WSSiPOM, 'POM_Si_in_sediment',             'g m-2', 'POM Si in sediment')
   call self%register_state_dependency(self%id_WSDHum,  'humus_DW_in_sediment',           'g m-2', 'humus DW in sediment')
   call self%register_state_dependency(self%id_WSNHum,  'humus_N_in_sediment',            'g m-2', 'humus N in sediment')
   call self%register_state_dependency(self%id_WSPHum,  'humus_P_in_sediment',            'g m-2', 'humus P in sediment')
!  Register dependencies to phytoplankton in water column
   call self%register_state_dependency(self%id_SWDDiat, 'diatom_DW_in_water',             'g m-3', 'diatom DW in water')
   call self%register_state_dependency(self%id_SWDGren, 'green_DW_in_water',              'g m-3', 'green DW in water')
   call self%register_state_dependency(self%id_SWDBlue, 'blue_DW_in_water',               'g m-3', 'blue-green DW in water')
   call self%register_state_dependency(self%id_SWNDiat, 'diatom_N_in_water',              'g m-3', 'diatom N in water')
   call self%register_state_dependency(self%id_SWNGren, 'green_N_in_water',               'g m-3', 'green N in water')
   call self%register_state_dependency(self%id_SWNBlue, 'blue_N_in_water',                'g m-3', 'blue-green N in water')
   call self%register_state_dependency(self%id_SWPDiat, 'diatom_P_in_water',              'g m-3', 'diatom P in water')
   call self%register_state_dependency(self%id_SWPGren, 'green_P_in_water',               'g m-3', 'green P in water')
   call self%register_state_dependency(self%id_SWPBlue, 'blue_P_in_water',                'g m-3', 'blue-green P in water')
!  Register dependencies to phytoplankton in sediment
   call self%register_state_dependency(self%id_WSDDiat, 'diatom_DW_in_sediment',          'g m-2', 'diatom DW in sediment')
   call self%register_state_dependency(self%id_WSDGren, 'green_DW_in_sediment',           'g m-2', 'green DW in sediment')
   call self%register_state_dependency(self%id_WSDBlue, 'blue_DW_in_sediment',            'g m-2', 'blue-green DW in sediment')
   call self%register_state_dependency(self%id_WSNDiat, 'diatom_N_in_sediment',           'g m-2', 'diatom N in sediment')
   call self%register_state_dependency(self%id_WSNGren, 'green_N_in_sediment',            'g m-2', 'green N in sediment')
   call self%register_state_dependency(self%id_WSNBlue, 'blue_N_in_sediment',             'g m-2', 'blue-green N in sediment')
   call self%register_state_dependency(self%id_WSPDiat, 'diatom_P_in_sediment',           'g m-2', 'diatom P in sediment')
   call self%register_state_dependency(self%id_WSPGren, 'green_P_in_sediment',            'g m-2', 'green P in sediment')
   call self%register_state_dependency(self%id_WSPBlue, 'blue_P_in_sediment',             'g m-2', 'blue-green P in sediment')
!  register macrophytes and fish for resuspension dependency
   call self%register_dependency(self%id_DragVeg, 'macrophytes_DW',                  'g m-2', 'macrophytes DW')
!   call self%register_state_dependency(self%id_DragVeg, 'macrophytes_DW',                  'g m-2', 'macrophytes DW')
   call self%register_state_dependency(self%id_TurbFish,'adult_fish_DW',                  'g m-2', 'adult fish DW')
!  register zooplankton for transport purpose
   call self%register_state_dependency(self%id_DTranZoo,'zooplankton_DW',                 'g m-3', 'zooplankton DW')
   call self%register_state_dependency(self%id_PTranZoo,'zooplankton_P',                  'g m-3', 'zooplankton P')
   call self%register_state_dependency(self%id_NTranZoo,'zooplankton_N',                  'g m-3', 'zooplankton N')
!  register environmental dependencies
   call self%register_dependency(self%id_uTm,    standard_variables%temperature)
   call self%register_dependency(self%id_dz,standard_variables%cell_thickness)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)
   call self%register_dependency(self%id_Day,    standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_shear,  standard_variables%bottom_stress)
!  Register dependencies on external diagnostic variables
   call self%register_dependency(self%id_tDAbioPOMS, 'POM_abiotic_update',      'g m-2 s-1', 'POM abiotic update')
   call self%register_dependency(self%id_tDAbioHumS, 'humus_abiotic_update',    'g m-2 s-1', 'humus abiotic update')
   call self%register_dependency(self%id_tDPrimPOMS, 'POM_from_algae',     '[-]',       'POM from algae')
   call self%register_dependency(self%id_tDWebPOMS,  'POM_from_zoobenthos','[-]',       'POM from zoobenthos')
   call self%register_dependency(self%id_tDBedPOMS,  'POM_from_macrophytes','[-]',       'POM from macrophytes')
!  register diagnostic variables
   call self%register_diagnostic_variable(self%id_tDBurIM,     'tDBurIM',     'g m-2 s-1','tDBurIM',                output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_shearstress, 'shearstress', 'N m-2',    'shearstress',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aFunDimSusp, 'aFunDimSusp', '[-]',      'aFunDimSusp',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDResusDead, 'tDResusDead', '[-]',      'tDResusDead',            output=output_time_step_integrated)
   call self%register_diagnostic_variable(self%id_aFunTauSet,  'aFunTauSet',  '[-]',      'basic settling rate',    output=output_time_step_integrated)
#ifdef _DEVELOPMENT_
!  register diagnostic variables for resuspension fluxes
   call self%register_diagnostic_variable(self%id_tAuxDIMW,    'tAuxDIMW',    'g m-2 s-1','auxiliary_DIMW_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAuxPOMW,   'tDAuxPOMW',   'g m-2 s-1','auxiliary_DPOMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxPOMW,   'tNAuxPOMW',   'g m-2 s-1','auxiliary_NPOMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxPOMW,   'tPAuxPOMW',   'g m-2 s-1','auxiliary_PPOMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tSiAuxPOMW,  'tSiAuxPOMW',  'g m-2 s-1','auxiliary_SiPOMW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tAuxPAIMW,   'tAuxPAIMW',   'g m-2 s-1','auxiliary_PAIMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxNH4W,   'tNAuxNH4W',   'g m-2 s-1','auxiliary_NH4W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxNO3W,   'tNAuxNO3W',   'g m-2 s-1','auxiliary_NO3W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxPO4W,   'tPAuxPO4W',   'g m-2 s-1','auxiliary_PO4W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAxuDiatW,  'tDAxuDiatW',  'g m-2 s-1','auxiliary_DDiatW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxDiatW,  'tNAuxDiatW',  'g m-2 s-1','auxiliary_NDiatW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxDiatW,  'tPAuxDiatW',  'g m-2 s-1','auxiliary_PDiatW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAuxGrenW,  'tDAuxGrenW',  'g m-2 s-1','auxiliary_DGrenW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxGrenW,  'tNAuxGrenW',  'g m-2 s-1','auxiliary_NGrenW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxGrenW,  'tPAuxGrenW',  'g m-2 s-1','auxiliary_PGrenW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAuxBlueW,  'tDAuxBlueW',  'g m-2 s-1','auxiliary_DBlueW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxBlueW,  'tNAuxBlueW',  'g m-2 s-1','auxiliary_NBlueW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxBlueW,  'tPAuxBlueW',  'g m-2 s-1','auxiliary_PBlueW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tAuxDIMS,    'tAuxDIMS',    'g m-2 s-1','auxiliary_DIMS_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAuxPOMS,   'tDAuxPOMS',   'g m-2 s-1','auxiliary_DPOMS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxPOMS,   'tNAuxPOMS',   'g m-2 s-1','auxiliary_NPOMS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxPOMS,   'tPAuxPOMS',   'g m-2 s-1','auxiliary_PPOMS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tSiAuxPOMS,  'tSiAuxPOMS',  'g m-2 s-1','auxiliary_SiPOMS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxPO4S,   'tPAuxPO4S',   'g m-2 s-1','auxiliary_PO4S_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tAuxPAIMS,   'tAuxPAIMS',   'g m-2 s-1','auxiliary_PAIMS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxNH4S,   'tNAuxNH4S',   'g m-2 s-1','auxiliary_NH4S_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxNO3S,   'tNAuxNO3S',   'g m-2 s-1','auxiliary_NO3S_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAuxHumS,   'tDAuxHumS',   'g m-2 s-1','auxiliary_DHumS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxHumS,   'tPAuxHumS',   'g m-2 s-1','auxiliary_PHumS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxHumS,   'tNAuxHumS',   'g m-2 s-1','auxiliary_NHumS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAxuDiatS,  'tDAxuDiatS',  'g m-2 s-1','auxiliary_DDiatS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxDiatS,  'tNAuxDiatS',  'g m-2 s-1','auxiliary_NDiatS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxDiatS,  'tPAuxDiatS',  'g m-2 s-1','auxiliary_PDiatS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAuxGrenS,  'tDAuxGrenS',  'g m-2 s-1','auxiliary_DGrenS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxGrenS,  'tNAuxGrenS',  'g m-2 s-1','auxiliary_NGrenS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxGrenS,  'tPAuxGrenS',  'g m-2 s-1','auxiliary_PGrenS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAuxBlueS,  'tDAuxBlueS',  'g m-2 s-1','auxiliary_DBlueS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAuxBlueS,  'tNAuxBlueS',  'g m-2 s-1','auxiliary_NBlueS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAuxBlueS,  'tPAuxBlueS',  'g m-2 s-1','auxiliary_PBlueS_change', output=output_instantaneous)
#endif

   return

   end subroutine initialize
!
!EOC
!-----------------------------------------------------------------------
!BOP
! !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_pclake_auxiliary), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
! !LOCAL VARIABLES:
!  carriers for state dependencies in different modules
!  in abiotic water column module
   real(rk)                   :: sDIMW,sDPOMW,sNPOMW,sPPOMW,sSiPOMW
   real(rk)                   :: sPO4W,sPAIMW,sNH4W,sNO3W
!  in phytoplankton water column module
   real(rk)                   :: sDDiatW,sDGrenW,sDBlueW
   real(rk)                   :: sNDiatW,sNGrenW,sNBlueW
   real(rk)                   :: sPDiatW,sPGrenW,sPBlueW
!  in abiotic water column module
   real(rk)                   :: sDIMS,sDPOMS,sNPOMS,sPPOMS,sSiPOMS
   real(rk)                   :: sPO4S,sPAIMS,sNH4S,sNO3S
   real(rk)                   :: sDHumS,sNHumS,sPHumS
!  in phytoplankton sediment module
   real(rk)                   :: sDDiatS,sDGrenS,sDBlueS
   real(rk)                   :: sNDiatS,sNGrenS,sNBlueS
   real(rk)                   :: sPDiatS,sPGrenS,sPBlueS
!  in macrophytes module
   real(rk)                   :: sDVeg
!  in fish module
   real(rk)                   :: sDFiAd
!  carriers for environmental dependencies
!  depth for empirical resuspension function, should be bottom_depth
   real(rk)                   :: uTm,sDepthW ,dz
!  carriers for diagnostic dependencies
   real(rk)                   :: tDAbioHumS
   real(rk)                   :: tDAbioPOMS,tDPrimPOMS,tDWebPOMS,tDBedPOMS
!  variables for nutrient rations for POM in sediment
   real(rk)                   :: rPDPOMS,rNDPOMS,rSiDPOMS
   real(rk)                   :: rPDHumS,rNDHumS
!  variables for nutrient rations for phytoplankton in water
   real(rk)                   :: rNDDiatW,rNDGrenW,rNDBlueW
   real(rk)                   :: rPDDiatW,rPDGrenW,rPDBlueW
!  variables for nutrient rationsfor phytoplankton in sediment
   real(rk)                   :: rNDDiatS,rNDGrenS,rNDBlueS
   real(rk)                   :: rPDDiatS,rPDGrenS,rPDBlueS
!  temperature related variables
   real(rk)                   :: uFunTmSed,uFunTmFish
!  variables related to resuspension (in the order of appearance)
   real(rk)                   :: aFunVegResus,tDTurbFish
   real(rk)                   :: aFunDimSusp,tDResusTauDead,tDResusBareDead
   real(rk)                   :: tDResusDead,tDResusIM,tDResusPOM,tPResusPOM
   real(rk)                   :: tNResusPOM,tSiResusPOM,tPResusPO4,tPResusAIM
   real(rk)                   :: tNResusNO3,tNResusNH4
!  variables for phytoplankton resuspension
   real(rk)                   :: akResusPhytRef
   real(rk)                   :: tDResusDiat,tDResusGren,tDResusBlue
   real(rk)                   :: tNResusDiat,tNResusGren,tNResusBlue
   real(rk)                   :: tPResusDiat,tPResusGren,tPResusBlue  !,tSiResusDiat
!  variables related to sedimentation (in the order of appearance)
   real(rk)                   :: aFunTauSet
   real(rk)                   :: uCorVSedIM,tDSetIM,tPSetAIM
   real(rk)                   :: uCorVSedPOM,tDSetPOM,tPSetPOM,tNSetPOM,tSiSetPOM
!  variables for phytoplankton sedimentation
   real(rk)                   :: uCorVSedDiat,uCorVSedGren,uCorVSedBlue
   real(rk)                   :: tDSetDiat,tDSetGren,tDSetBlue
   real(rk)                   :: tNSetDiat,tNSetGren,tNSetBlue
   real(rk)                   :: tPSetDiat,tPSetGren,tPSetBlue !,tSiSetDiat
!  Variables related to burial process (in the order of appearance)
   real(rk)                   :: tDIMS,tDPOMS,vDeltaS
   real(rk)                   :: tDBurIM,tDBurPOM,tPBurPOM,tPBurAIM,tPBurPO4
   real(rk)                   :: tNBurPOM,tNBurNH4,tNBurNO3,tSiBurPOM
!  Humus variables
   real(rk)                   :: tDHumS,tDBurHum,tDBurOM,tNBurHum,tPBurHum
!  variables of new resuspension method
   real(rk)                   :: shear
!
!EOP
!-----------------------------------------------------------------------
!BOC
!------------------------------------------------------------------------
!  Spatial loop
   _FABM_HORIZONTAL_LOOP_BEGIN_
!-----------------------------------------------------------------------
!  Retrieve dependencies  value
!-----------------------------------------------------------------------
!  Retrieve state dependencies value
!  from abiotic water column
   _GET_(self%id_SWNH4,sNH4W)
   _GET_(self%id_SWNO3,sNO3W)
   _GET_(self%id_SWPO4,sPO4W)
   _GET_(self%id_SWPAIM,sPAIMW)
   _GET_(self%id_SWDIM,sDIMW)
   _GET_(self%id_SWDPOM,sDPOMW)
   _GET_(self%id_SWNPOM,sNPOMW)
   _GET_(self%id_SWPPOM,sPPOMW)
   _GET_(self%id_SWSiPOM,sSiPOMW)
!  from phytoplankton in water column
   _GET_(self%id_SWDDiat,sDDiatW)
   _GET_(self%id_SWDGren,sDGrenW)
   _GET_(self%id_SWDBlue,sDBlueW)
   _GET_(self%id_SWNDiat,sNDiatW)
   _GET_(self%id_SWNGren,sNGrenW)
   _GET_(self%id_SWNBlue,sNBlueW)
   _GET_(self%id_SWPDiat,sPDiatW)
   _GET_(self%id_SWPGren,sPGrenW)
   _GET_(self%id_SWPBlue,sPBlueW)
!  from abiotic sediment
   _GET_HORIZONTAL_(self%id_WSNH4,sNH4S)
   _GET_HORIZONTAL_(self%id_WSNO3,sNO3S)
   _GET_HORIZONTAL_(self%id_WSPO4,sPO4S)
   _GET_HORIZONTAL_(self%id_WSPAIM,sPAIMS)
   _GET_HORIZONTAL_(self%id_WSDIM,sDIMS)
   _GET_HORIZONTAL_(self%id_WSDPOM,sDPOMS)
   _GET_HORIZONTAL_(self%id_WSNPOM,sNPOMS)
   _GET_HORIZONTAL_(self%id_WSPPOM,sPPOMS)
   _GET_HORIZONTAL_(self%id_WSSiPOM,sSiPOMS)
   _GET_HORIZONTAL_(self%id_WSDHum,sDHumS)
   _GET_HORIZONTAL_(self%id_WSNHum,sNHumS)
   _GET_HORIZONTAL_(self%id_WSPHum,sPHumS)
!  from phytoplankton in sediment
   _GET_HORIZONTAL_(self%id_WSDDiat,sDDiatS)
   _GET_HORIZONTAL_(self%id_WSDGren,sDGrenS)
   _GET_HORIZONTAL_(self%id_WSDBlue,sDBlueS)
   _GET_HORIZONTAL_(self%id_WSNDiat,sNDiatS)
   _GET_HORIZONTAL_(self%id_WSNGren,sNGrenS)
   _GET_HORIZONTAL_(self%id_WSNBlue,sNBlueS)
   _GET_HORIZONTAL_(self%id_WSPDiat,sPDiatS)
   _GET_HORIZONTAL_(self%id_WSPGren,sPGrenS)
   _GET_HORIZONTAL_(self%id_WSPBlue,sPBlueS)
!  vegatation influence on macrophytes
   _GET_HORIZONTAL_(self%id_DragVeg,sDVeg)
!  fish predation influence on resuspension
   _GET_HORIZONTAL_(self%id_TurbFish,sDFiAd)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_(self%id_dz,dz)
   _GET_HORIZONTAL_(self%id_shear,shear)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
!  fish biomass converted to g/m^2
!   sDFiAd=sDFiAd*sDepthW
!  retrieve diagnostic dependencies
   _GET_HORIZONTAL_(self%id_tDAbioPOMS,tDAbioPOMS)
   _GET_HORIZONTAL_(self%id_tDPrimPOMS,tDPrimPOMS)
   _GET_HORIZONTAL_(self%id_tDWebPOMS,tDWebPOMS)
   _GET_HORIZONTAL_(self%id_tDBedPOMS,tDBedPOMS)
   _GET_HORIZONTAL_(self%id_tDAbioHumS,tDAbioHumS)
!-----------------------------------------------------------------------
!  Current nutrients ratios (check the current state)
!-----------------------------------------------------------------------
   rPDPOMS=sPPOMS/(sDPOMS+NearZero)
   rNDPOMS=sNPOMS/(sDPOMS+NearZero)
   rSiDPOMS=sSiPOMS/(sDPOMS+NearZero)
   rPDHumS=sPHumS/(sDHumS+NearZero)
   rNDHumS=sNHumS/(sDHumS+NearZero)
!  external source states
!  for phytoplankton in water
   rPDDiatW = sPDiatW /(sDDiatW+NearZero)
   rNDDiatW = sNDiatW /(sDDiatW+NearZero)
   rPDGrenW = sPGrenW /(sDGrenW+NearZero)
   rNDGrenW = sNGrenW /(sDGrenW+NearZero)
   rPDBlueW = sPBlueW /(sDBlueW+NearZero)
   rNDBlueW = sNBlueW /(sDBlueW+NearZero)
!  for phytoplankton in sediment
   rPDDiatS = sPDiatS /(sDDiatS+NearZero)
   rNDDiatS = sNDiatS /(sDDiatS+NearZero)
   rPDGrenS = sPGrenS /(sDGrenS+NearZero)
   rNDGrenS = sNGrenS /(sDGrenS+NearZero)
   rPDBlueS = sPBlueS /(sDBlueS+NearZero)
   rNDBlueS = sNBlueS /(sDBlueS+NearZero)
!-----------------------------------------------------------------------
!  Temperature functions for sediment abiotic process
!-----------------------------------------------------------------------
!  temperature_correction_of_sedimentation
   uFunTmSed= uFunTmAbio(uTm,self%cThetaSed)
   uFunTmFish= uFunTmBio(uTm,self%cSigTmFish,self%cTmOptFish)
!-----------------------------------------------------------------------
!  Process related to other modules
!-----------------------------------------------------------------------
!  macrophytes_dependence_of_resuspension
   aFunVegResus=max(0.0_rk,1.0_rk-self%kVegResus*sDVeg)
!-----------------------------------------------------------------------
!  Resuspension and sedimentation 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  Calculate resuspension rate
!-----------------------------------------------------------------------
!  bioturbation_by_fish
      tDTurbFish=self%kTurbFish*uFunTmFish*sDFiAd
!  calculate resuspension rate, two methods
   select case(self%resusp_meth)
      case(1) ! (original PCLake method)
!     Empirical_resuspension_function_(logistic_fit_to_data), in day-1
          if (uTm >= 0.1_rk) then
              aFunDimSusp=self%cSuspRef*((self%cSuspMin+self%cSuspMax/(1.0_rk+exp(self%cSuspSlope*&
              & (sDepthW-self%hDepthSusp))))*((((self%cFetch +NearZero)/ self%cFetchRef) )** (0.5_rk)))
          else
              aFunDimSusp=0.0_rk
          endif
      case(2) ! (resuspension based on shear stress from physical driver model)
         if(shear <= self%crt_shear)then
            aFunDimSusp=0.0_rk
         else
            aFunDimSusp=(self%alpha/self%cVSedMain)*((shear-self%crt_shear)/self%ref_shear)**self%eta
         endif
   end select
   tDResusTauDead=min(aFunDimSusp, ((aFunDimSusp +NearZero )**(0.5_rk))) &
   &*((self%fLutum/ self%fLutumRef )** (0.5_rk))*self%bPorS
!  resuspension_due_to_shear_stress_AND_fish, in day-1
   tDResusBareDead=tDResusTauDead+tDTurbFish
!  resuspension,_corrected_for_macrophytes_effect, in secs-1
!   tDResusDead=tDResusBareDead*aFunVegResus
   tDResusDead=tDResusBareDead*aFunVegResus/secs_pr_day
!------------------------------------------------------------------------------------------------------------
!  Resuspension rates of individual components in the sediment
!------------------------------------------------------------------------------------------------------------
!  The inorganic matter resuspension
   tDResusIM=self%fLutum*sDIMS/(self%fLutum*sDIMS+sDPOMS)*tDResusDead
!  POMrital_resuspension_DW
   tDResusPOM=sDPOMS/(self%fLutum*sDIMS+sDPOMS)*tDResusDead
!  POMrital_resuspension_P
   tPResusPOM=rPDPOMS*tDResusPOM
!  POM_resuspension_N
   tNResusPOM=rNDPOMS*tDResusPOM
!  POM_resuspension_SI
   tSiResusPOM=rSiDPOMS*tDResusPOM
!  resuspension_nutrient_P
   tPResusPO4=sPO4S/sDPOMS*tDResusPOM
!  resuspension_adsorbed_PAIM
   tPResusAIM=sPAIMS/sDIMS*tDResusIM
!  resuspension_nutrient_NO3
   tNResusNO3=sNO3S/sDPOMS*tDResusPOM
!  resuspension_nutrient_NH4
   tNResusNH4=sNH4S/sDPOMS*tDResusPOM
!  convert the seconds rate to daily rate, for the equation purpose
   tDResusDead=tDResusDead*secs_pr_day
!  phytoplankton_resuspension_rate_constant, in day
   akResusPhytRef = self%kResusPhytMax * (1.0_rk - exp(self%cResusPhytExp * tDResusDead))
!  convert to secs
   akResusPhytRef=akResusPhytRef/secs_pr_day
!  Algae group resuspension
!  resuspension of Diatoms,DW
   tDResusDiat=akResusPhytRef*sDDiatS
!  resuspension of Green algae,DW
   tDResusGren=akResusPhytRef*sDGrenS
!  resuspension of Blue-green algae,DW
   tDResusBlue=akResusPhytRef*sDBlueS
!  resuspension of Diatoms,N
   tNResusDiat = rNDDiatS * tDResusDiat
!  resuspension of Green algae,N
   tNResusGren = rNDGrenS * tDResusGren
!  resuspension of Blue-green algae,N
   tNResusBlue = rNDBlueS * tDResusBlue
!  resuspension of Diatom,P
   tPResusDiat = rPDDiatS * tDResusDiat
!  resuspension of Green algae,P
   tPResusGren = rPDGrenS * tDResusGren
!  resuspension of Blue-green algae,P
   tPResusBlue = rPDBlueS * tDResusBlue
!-----------------------------------------------------------------------
!  The net sedimentation calculation, based on resuspension
!-----------------------------------------------------------------------
!  correction_factor_for_settling_rate_(<=_1),net settling rate, in day
   aFunTauSet=min(1.0_rk,1.0_rk/((aFunDimSusp +NearZero )**(0.5_rk)))
!-----------------------------------------------------------------------
!  Sedimentation rates of individual components in the water column based on the net settling rate
!-----------------------------------------------------------------------
!  sedimentation_velocity_of_IM, in day
   uCorVSedIM=aFunTauSet*((self%fLutumRef/self%fLutum)**(0.5))*uFunTmSed*self%cVSedIM
!  convert to seconds
   uCorVSedIM=uCorVSedIM/secs_pr_day
!  sedimentation_IM
   tDSetIM=uCorVSedIM*sDIMW
!  sedimentation_PAIM
   tPSetAIM=sPAIMW/(sDIMW +NearZero)*tDSetIM
!  sedimentation_velocity_of_POM, in day
   uCorVSedPOM=self%cVSedPOM*aFunTauSet*uFunTmSed
!  convert to seconds
   uCorVSedPOM=uCorVSedPOM/secs_pr_day
!  sedimentation_flux_of_POM
   tDSetPOM=uCorVSedPOM*sDPOMW
!  sedimentation_POM_P
   tPSetPOM=uCorVSedPOM*sPPOMW
!  sedimentation_POM_N
   tNSetPOM=uCorVSedPOM*sNPOMW
!  sedimentation_POM_Si
   tSiSetPOM=uCorVSedPOM*sSiPOMW
!  corrected_sedimentation_velocity_of_Algae, in day
   uCorVSedDiat = self%cVSedDiat * aFunTauSet * uFunTmSed
!  convert to seconds
   uCorVSedDiat=uCorVSedDiat/secs_pr_day
!  sedimentation_flux_of_Diatom
   tDSetDiat = uCorVSedDiat * sDDiatW
!  corrected_sedimentation_velocity_of_Algae,in day
   uCorVSedGren = self%cVSedGren * aFunTauSet * uFunTmSed
!  convert to seconds
   uCorVSedGren=uCorVSedGren/secs_pr_day
!  sedimentation_flux_of_Algae
   tDSetGren = uCorVSedGren * sDGrenW
!  corrected_sedimentation_velocity_of_Algae, in day
   uCorVSedBlue = self%cVSedBlue * aFunTauSet * uFunTmSed
!  convert to seconds
   uCorVSedBlue=uCorVSedBlue/secs_pr_day
!  sedimentation_flux_of_Algae
   tDSetBlue = uCorVSedBlue * sDBlueW
!  sedimentation
   tNSetDiat = rNDDiatW * tDSetDiat
!  sedimentation
   tNSetGren = rNDGrenW * tDSetGren
!  sedimentation
   tNSetBlue = rNDBlueW * tDSetBlue
!  sedimentation
   tPSetDiat = rPDDiatW * tDSetDiat
!  sedimentation
   tPSetGren = rPDGrenW * tDSetGren
!  sedimentation
   tPSetBlue = rPDBlueW * tDSetBlue
!  Diatoms_sedimentation
!   tSiSetDiat = self%cSiDDiat * tDSetDiat
!-----------------------------------------------------------------------
!  Burial of sediment components
!-----------------------------------------------------------------------
!  increase_in_inorganic_matter_in_sediment
   tDIMS=  tDSetIM - tDResusIM
!  increase_in_sediment_humus_in_lake
   tDHumS = tDAbioHumS  ! uDErosOM+
!  increase_in_sediment_POM_in_lake
   tDPOMS= tDSetPOM - tDResusPOM+tDAbioPOMS+ tDPrimPOMS + tDWebPOMS + tDBedPOMS
!  turnover_depth_in_lake
   vDeltaS = (tDIMS / self%cRhoIM + (tDHumS + tDPOMS) / self%cRhoOM)/(1.0_rk - self%bPorS)
!  burial_flux_of_DW_in_inorganic_matter_in_lake
   if (vDeltaS >= 0.0_rk) then
      tDBurIM = ((tDHumS + tDPOMS) +(self%cRhoOM / self%cRhoIM) * tDIMS) / ((sDHumS + sDPOMS) /sDIMS &
      & + self%cRhoOM / self%cRhoIM)
   else
      tDBurIM = ( (tDHumS + tDPOMS) +(self%cRhoOM / self%cRhoIM) * tDIMS) / (self%fDOrgSoil &
      &/(1.0_rk - self%fDOrgSoil) + self%cRhoOM / self%cRhoIM)
   endif
!  burial_flux_of_DW_in_organic_matter_in_lake
   if (vDeltaS >= 0.0) then
      tDBurOM = (sDHumS + sDPOMS) / sDIMS * tDBurIM
   else
      tDBurOM = self%fDOrgSoil /(1.0 - self%fDOrgSoil) * tDBurIM
   endif
!  burial_flux_of_DW_in_POM_in_lake
   if (vDeltaS >= 0.0) then
      tDBurPOM = sDPOMS /(sDHumS + sDPOMS) * tDBurOM
   else
      tDBurPOM = 0.0_rk
   endif
!  burial_flux_of_P_in_POM_in_lake
   if (vDeltaS >= 0.0_rk) then
      tPBurPOM = rPDPOMS * tDBurPOM
   else
      tPBurPOM = 0.0_rk
   endif
!  burial_flux_of_P_adsorbed_onto_inorganic_matter_in_lake
   if (vDeltaS >= 0.0_rk) then
      tPBurAIM = sPAIMS / sDIMS * tDBurIM
   else
      tPBurAIM = 0.0_rk
   endif
!  burial_flux_of_dissolved_P_in_lake
   if (vDeltaS >= 0.0_rk) then
      tPBurPO4 = sPO4S *(vDeltaS / self%cDepthS)
   else
      tPBurPO4 = self%cPO4Ground *(self%bPorS * vDeltaS)
   endif
!  burial_flux_of_N_in_POMritus_in_lake
   if (vDeltaS >= 0.0_rk) then
      tNBurPOM =rNDPOMS * tDBurPOM
   else
      tNBurPOM = 0.0_rk
   endif
!  burial_flux_of_dissolved_NH4_in_lake
   if (vDeltaS >= 0.0_rk) then
      tNBurNH4 = sNH4S *(vDeltaS /self%cDepthS)
   else
      tNBurNH4 = self%cNH4Ground *(self%bPorS * vDeltaS)
   endif
!  burial_flux_of_dissolved_NO3_in_lake
   if (vDeltaS >= 0.0_rk) then
      tNBurNO3 = sNO3S *(vDeltaS / self%cDepthS)
   else
      tNBurNO3 = self%cNO3Ground *(self%bPorS * vDeltaS)
   endif

!  burial_flux_of_Si_in_POM_in_lake
   if (vDeltaS >= 0.0_rk) then
      tSiBurPOM = rSiDPOMS * tDBurPOM
   else
      tSiBurPOM = 0.0_rk
   endif
! Humus burial fluxes
!  burial_flux_of_DW_in_humus_in_lake
   if (vDeltaS >= 0.0) then
     tDBurHum = tDBurOM - tDBurPOM
   else
     tDBurHum = tDBurOM
   endif
!  burial_flux_of_P_in_humus_in_lake
   if (vDeltaS >= 0.0) then
      tPBurHum = rPDHumS * tDBurHum
   else
      tPBurHum = 0.001_rk * tDBurHum   ! cPDSoilOM=0.001
   endif
!  burial_flux_of_N_in_humus_in_lake
   if (vDeltaS >= 0.0) then
      tNBurHum = rNDHumS * tDBurHum
   else
      tNBurHum = 0.01_rk * tDBurHum   !cNDSoilOM =0.01
   endif
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
!  update inorganic and organic matters in water column
   _SET_BOTTOM_EXCHANGE_(self%id_SWDIM,  tDResusIM-tDSetIM)
   _SET_BOTTOM_EXCHANGE_(self%id_SWDPOM, tDResusPOM-tDSetPOM)
   _SET_BOTTOM_EXCHANGE_(self%id_SWNPOM, tNResusPOM-tNSetPOM)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPPOM, tPResusPOM-tPSetPOM)
   _SET_BOTTOM_EXCHANGE_(self%id_SWSiPOM,tSiResusPOM-tSiSetPOM)
!  update dissolved nutrients in water column
   _SET_BOTTOM_EXCHANGE_(self%id_SWNH4,  tNResusNH4)
   _SET_BOTTOM_EXCHANGE_(self%id_SWNO3,  tNResusNO3)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPO4,  tPResusPO4)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPAIM, tPResusAIM-tPSetAIM)
!  update phytoplankton in water column
   _SET_BOTTOM_EXCHANGE_(self%id_SWDDiat,tDResusDiat-tDSetDiat)
   _SET_BOTTOM_EXCHANGE_(self%id_SWNDiat,tNResusDiat-tNSetDiat)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPDiat,tPResusDiat-tPSetDiat)
!   _SET_BOTTOM_EXCHANGE_(self%id_SWSiDiat,tSiResusDiat-tSiSetDiat)
   _SET_BOTTOM_EXCHANGE_(self%id_SWDGren,tDResusGren-tDSetGren)
   _SET_BOTTOM_EXCHANGE_(self%id_SWNGren,tNResusGren-tNSetGren)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPGren,tPResusGren-tPSetGren)
   _SET_BOTTOM_EXCHANGE_(self%id_SWDBlue,tDResusBlue-tDSetBlue)
   _SET_BOTTOM_EXCHANGE_(self%id_SWNBlue,tNResusBlue-tNSetBlue)
   _SET_BOTTOM_EXCHANGE_(self%id_SWPBlue,tPResusBlue-tPSetBlue)
!  update abiotic variables in sediment
   _SET_ODE_BEN_(self%id_WSDIM,  tDSetIM-tDResusIM-tDBurIM)
   _SET_ODE_BEN_(self%id_WSDPOM, tDSetPOM-tDResusPOM-tDBurPOM)
   _SET_ODE_BEN_(self%id_WSPPOM, tPSetPOM-tPResusPOM-tPBurPOM)
   _SET_ODE_BEN_(self%id_WSNPOM, tNSetPOM-tNResusPOM-tNBurPOM)
   _SET_ODE_BEN_(self%id_WSSiPOM,tSiSetPOM-tSiResusPOM-tSiBurPOM)
   _SET_ODE_BEN_(self%id_WSPO4,  -tPResusPO4-tPBurPO4)
   _SET_ODE_BEN_(self%id_WSPAIM, tPSetAIM-tPResusAIM-tPBurAIM)
   _SET_ODE_BEN_(self%id_WSNH4,  -tNResusNH4-tNBurNH4)
   _SET_ODE_BEN_(self%id_WSNO3,  -tNResusNO3-tNBurNO3)
   _SET_ODE_BEN_(self%id_WSDHum, tDBurHum)
   _SET_ODE_BEN_(self%id_WSPHum, tPBurHum)
   _SET_ODE_BEN_(self%id_WSNHum, tNBurHum)
!  update settled phytoplankton
   _SET_ODE_BEN_(self%id_WSDDiat, tDSetDiat-tDResusDiat)
   _SET_ODE_BEN_(self%id_WSNDiat, tNSetDiat-tNResusDiat)
   _SET_ODE_BEN_(self%id_WSPDiat, tPSetDiat-tPResusDiat)
   _SET_ODE_BEN_(self%id_WSDGren, tDSetGren-tDResusGren)
   _SET_ODE_BEN_(self%id_WSNGren, tNSetGren-tNResusGren)
   _SET_ODE_BEN_(self%id_WSPGren, tPSetGren-tPResusGren)
   _SET_ODE_BEN_(self%id_WSDBlue, tDSetBlue-tDResusBlue)
   _SET_ODE_BEN_(self%id_WSNBlue, tNSetBlue-tNResusBlue)
   _SET_ODE_BEN_(self%id_WSPBlue, tPSetBlue-tPResusBlue)
!  update diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBurIM,     tDBurIM)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_shearstress, shear)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aFunDimSusp, aFunDimSusp)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDResusDead, tDResusDead)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aFunTauSet,  aFunTauSet)
#ifdef _DEVELOPMENT_
!  output diagonostic variable for resuspension fluxes
!  fluxes for abiotic water state variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tAuxDIMW,   (tDResusIM-tDSetIM)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAuxPOMW,  (tDResusPOM-tDSetPOM)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxPOMW,  (tNResusPOM-tNSetPOM)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxPOMW,  (tPResusPOM-tPSetPOM)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tSiAuxPOMW, (tSiResusPOM-tSiSetPOM)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tAuxPAIMW,  (tPResusAIM-tPSetAIM)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxNH4W,  tNResusNH4/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxNO3W,  tNResusNO3/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxPO4W,  tPResusPO4/dz*secs_pr_day)
!  fluxes for phytoplankton water state variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAxuDiatW, (tDResusDiat-tDSetDiat)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxDiatW, (tNResusDiat-tNSetDiat)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxDiatW, (tPResusDiat-tPSetDiat)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAuxGrenW, (tDResusGren-tDSetGren)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxGrenW, (tNResusGren-tNSetGren)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxGrenW, (tPResusGren-tPSetGren)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAuxBlueW, (tDResusBlue-tDSetBlue)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxBlueW, (tNResusBlue-tNSetBlue)/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxBlueW, (tPResusBlue-tPSetBlue)/dz*secs_pr_day)
!  fluxes for abiotic sediment state variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tAuxDIMS,   (tDSetIM-tDResusIM-tDBurIM) *secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAuxPOMS,  (tDSetPOM-tDResusPOM-tDBurPOM)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxPOMS,  (tPSetPOM-tPResusPOM-tPBurPOM)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxPOMS,  (tNSetPOM-tNResusPOM-tNBurPOM)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tSiAuxPOMS, (tSiSetPOM-tSiResusPOM-tSiBurPOM)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxPO4S,  (-tPResusPO4-tPBurPO4)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tAuxPAIMS,  (tPSetAIM-tPResusAIM-tPBurAIM)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxNH4S,  (-tNResusNH4-tNBurNH4)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxNO3S,  (-tNResusNO3-tNBurNO3)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAuxHumS,  tDBurHum*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxHumS,  tPBurHum*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxHumS,  tNBurHum*secs_pr_day)
!  fluxes for phytoplankton sediment state variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAxuDiatS, (tDSetDiat-tDResusDiat)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxDiatS, (tNSetDiat-tNResusDiat)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxDiatS, (tPSetDiat-tPResusDiat)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAuxGrenS, (tDSetGren-tDResusGren)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxGrenS, (tNSetGren-tNResusGren)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxGrenS, (tPSetGren-tPResusGren)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAuxBlueS, (tDSetBlue-tDResusBlue)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAuxBlueS, (tNSetBlue-tNResusBlue)*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAuxBlueS, (tPSetBlue-tPResusBlue)*secs_pr_day)
#endif


   _FABM_HORIZONTAL_LOOP_END_
! Spatial loop end
   end subroutine do_bottom
!EOC
!
!  IROUTINE: this subroutine deals with the atmospheric depositions
!  including POM nitrogen, ammonium, nitrate, phosphate,
!  POM phosphorus
 subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_pclake_auxiliary),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!  local variables

!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_

   _SET_SURFACE_EXCHANGE_(self%id_SWDIM,  self%tDDepoIM)
   _SET_SURFACE_EXCHANGE_(self%id_SWDPOM, self%tDDepoPOM)
   _SET_SURFACE_EXCHANGE_(self%id_SWNPOM, self%tNDepoPOM)
   _SET_SURFACE_EXCHANGE_(self%id_SWPPOM, self%tPDepoPOM)
   _SET_SURFACE_EXCHANGE_(self%id_SWNH4,  self%tNDepoNH4)
   _SET_SURFACE_EXCHANGE_(self%id_SWNO3,  self%tNDepoNO3)
   _SET_SURFACE_EXCHANGE_(self%id_SWPO4,  self%tPDepoPO4)
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC

!------------------------------------------------------------------------------
   end module pclake_auxiliary
!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
