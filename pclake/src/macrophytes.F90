#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module pclake_macrophytes
! !USES:
   use fabm_types
   use fabm_expressions
   use fabm_standard_variables
   use pclake_utility, ONLY:uFunTmVeg, DayLength

   implicit none

!  default: all is private.
   private

! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model), public :: type_pclake_macrophytes
!  local state variable identifiers
!  id_sDVeg, macrophytes in dry-weight, gDW/m**2
!  id_sPVeg, macrophytes in nitrogen element, gN/m**2
!  id_sPVeg, macrophytes in phosphorus element, gP/m**2
   type (type_bottom_state_variable_id)            :: id_sDVeg,id_sNVeg,id_sPVeg
!  diagnostic variables for dependencies (not for output)
   type (type_horizontal_diagnostic_variable_id)       :: id_aDSubVeg,id_aCovVeg
   type (type_horizontal_diagnostic_variable_id)       :: id_tDBedPOMS,id_afCovSurfVeg
!  diagonostic variables for light attenuation coefficient for macrophytes
   type (type_horizontal_diagnostic_variable_id)       :: id_aDayInitVeg
!  diagnostic variables for modular fluxes
   type (type_horizontal_diagnostic_variable_id)       :: id_tO2BedW,id_tDBedVeg
   type (type_horizontal_diagnostic_variable_id)       :: id_tNBedVeg,id_tPBedVeg
#ifdef _DEVELOPMENT_
   type (type_horizontal_diagnostic_variable_id)       :: id_wNBedNH4W,id_wNBedNO3W,id_wPBedPO4W
   type (type_horizontal_diagnostic_variable_id)       :: id_wDBedPOMW,id_wNBedPOMW
   type (type_horizontal_diagnostic_variable_id)       :: id_wPBedPOMW,id_tNBedNH4S,id_tNBedNO3S
   type (type_horizontal_diagnostic_variable_id)       :: id_tPBedPO4S,id_tPBedPOMS,id_tNBedPOMS
   type (type_horizontal_diagnostic_variable_id)       :: id_tDBedPOMSflux
#endif
!  state dependencys identifiers
   type (type_state_variable_id)                :: id_NH4poolW,id_NO3poolW,id_PO4poolW,id_O2poolW
   type (type_state_variable_id)                :: id_DPOMpoolW,id_NPOMpoolW,id_PPOMpoolW
   type (type_bottom_state_variable_id)         :: id_NH4poolS,id_NO3poolS,id_PO4poolS
   type (type_bottom_state_variable_id)         :: id_DPOMpoolS,id_NPOMpoolS,id_PPOMpoolS
   type (type_state_variable_id)                :: id_DDOMpoolW,id_NDOMpoolW,id_PDOMpoolW
   type (type_bottom_state_variable_id)         :: id_DDOMpoolS,id_NDOMpoolS,id_PDOMpoolS
!  environmental dependencies
   type (type_global_dependency_id)         :: id_Day
   type (type_dependency_id)                :: id_uTm,id_extc,id_dz
   
   type (type_horizontal_dependency_id)     :: id_I_0
   type (type_horizontal_dependency_id)     :: id_sDepthW,id_lat
!  output light variables
   type (type_horizontal_diagnostic_variable_id)       :: id_par,id_meanpar
   type (type_horizontal_diagnostic_variable_id)       :: id_aLPAR1Veg,id_aLPAR2Veg
   type (type_horizontal_diagnostic_variable_id)       :: id_aLLimVeg,id_aNutLimVeg
   type (type_horizontal_diagnostic_variable_id)       :: id_macroextinction
!  diagnostic dependencies
   type (type_horizontal_dependency_id)     :: id_afOxySed
!  Model parameters
!  Primary production parameters
   real(rk)    :: cDVegIn,kMigrVeg,cMuMaxVeg,cDCarrVeg,kDRespVeg
   real(rk)    :: hLRefVeg
!  nutrient ratio parameters
   real(rk)    :: cNDVegMin,cNDVegMax,cPDVegMin,cPDVegMax
!  macrophytes shoots and roots allocation parameters
   real(rk)    :: cDayWinVeg,cLengAllo
   real(rk)    :: fRootVegWin,fRootVegSum,cTmInitVeg
   real(rk)    :: fEmergVeg,fFloatVeg,cDLayerVeg,cCovSpVeg
   real(rk)    :: kMortVegSum,cLengMort,fWinVeg
!  temperature function parameters
   real(rk)    :: cQ10ProdVeg,cQ10RespVeg
!  parameters for nitrogen and phosphorus equations
   real(rk)    :: cPDVeg0,cNDVeg0
   real(rk)    :: fSedUptVegMax,fSedUptVegCoef,fSedUptVegExp,cAffNUptVeg,cVNUptMaxVeg
   integer     :: UseEmpUpt
   real(rk)    :: cVPUptMaxVeg,cAffPUptVeg,fDissMortVeg
!  paremters for sediment properties(pore water concentration)
   real(rk)   :: cDepthS,bPorS,cCPerDW,hO2BOD
   real(rk)   :: cExtSpVeg,fTOMWMortVeg
!  plant height
   real(rk)   :: cHeightVeg,fPAR,fRefl
!  minimun state variable vaules
   real(rk)   :: cDVegMin, cNVegMin,cPVegMin
!  fraction of dissolved organic matter from macrophytes
   real(rk)   :: fVegDOMW,fVegDOMS
!  special local parameter aDSubVeg, borrowed from do_bottom
!  subroutine to get_light_extinction subroutine
   real(rk)   :: aDSubVeg

   contains

!  Model procedure
   procedure :: initialize
   procedure :: do_bottom
   procedure :: get_light_extinction
   end type type_pclake_macrophytes

!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
   real(rk),parameter :: Pi=3.14159265358979_rk
!  ratio of mol.weights, = 32/12 [gO2/gC],
   real(rk),parameter :: molO2molC = 2.6667_rk
!  mol_O2_formed_per_mol_NO3-_ammonified
   real(rk),parameter ::O2PerNO3 = 1.5_rk
!  ratio of mol.weights,32/14 [gO2/gN],
   real(rk),parameter :: molO2molN = 2.2857_rk
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
   class (type_pclake_macrophytes), intent(inout), target :: self
   integer,                          intent(in)            :: configunit
!EOP
!-----------------------------------------------------------------------
!BOC
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%fPAR,       'fPAR',               '[-]',                  'the fraction photosynthetically active radiation',                                                        default=0.48_rk)
   call self%get_parameter(self%fRefl,      'fRefl',              '[-]',                  'the fraction reflected at the surface',                                                                   default=0.1_rk)
   call self%get_parameter(self%cDVegIn,       'cDVegIn',        'gD/m2',               'external macrophytes density',                                                                            default=1.0_rk)
   call self%get_parameter(self%kMigrVeg,      'kMigrVeg',       'd-1',                 'macrophytes migration rate',                                                                              default=0.001_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cNDVegMin,     'cNDVegMin',      'mgN/mgD',             'minimum N/DW ratio macrophytes',                                                                          default=0.01_rk)
   call self%get_parameter(self%cNDVegMax,     'cNDVegMax',      'mgN/mgD',             'maximum N/DW ratio macrophytes',                                                                          default=0.035_rk)
   call self%get_parameter(self%cPDVegMin,     'cPDVegMin',      'mgP/mg',              'minimum P/DW ratio macrophytes',                                                                          default=0.0008_rk)
   call self%get_parameter(self%cPDVegMax,     'cPDVegMax',      'mgP/mgD',             'maximum P/DW ratio macrophytes',                                                                          default=0.0035_rk)
   call self%get_parameter(self%cMuMaxVeg,     'cMuMaxVeg',      'd-1',                 'maximum growth rate of macrophytes at 20 degree C',                                                       default=0.2_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDCarrVeg,     'cDCarrVeg',      'gDW m-2',             'maximum macrophytes standing crop',                                                                       default=400.0_rk)
   call self%get_parameter(self%kDRespVeg,     'kDRespVeg',      'd-1',                 'dark respiration rate of macrophytes',                                                                    default=0.02_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDayWinVeg,    'cDayWinVeg',     'day of the year',     'day of the year for the end of growing season',                                                                                   default=259.0_rk)
   call self%get_parameter(self%cLengAllo,     'cLengAllo',      'days',                'duration of allocation and reallocation phase',                                                           default=15.0_rk)
   call self%get_parameter(self%fRootVegWin,   'fRootVegWin',    'groot/gveg',          'root fraction outside growing season',                                                                    default=0.6_rk)
   call self%get_parameter(self%fRootVegSum,   'fRootVegSum',    'groot/gveg',          'root fraction outside growing season',                                                                    default=0.1_rk)
   call self%get_parameter(self%cTmInitVeg,    'cTmInitVeg',     'degree C',            'temperature for onset of initial growth',                                                                          default=9.0_rk)
   call self%get_parameter(self%fEmergVeg,     'fEmergVeg',      'gfloating/gshoot',    'emergent fraction of shoot',                                                                              default=0.0_rk)
   call self%get_parameter(self%fFloatVeg,     'fFloatVeg',      'gfloating/gshoot',    'floating fraction of shoot',                                                                              default=0.0_rk)
   call self%get_parameter(self%cDLayerVeg,    'cDLayerVeg',     'gD/m2',               'biomass of a single layer floating leaves',                                                               default=0.0_rk)
   call self%get_parameter(self%cCovSpVeg,     'cCovSpVeg',      'l/gDW/m2',            'specific cover',                                                                                          default=0.5_rk)
   call self%get_parameter(self%hLRefVeg,      'hLRefVeg',       'W/m2PAR',             'half-saturtion for influence of light on macrophytes',                                                                          default=17.0_rk)
   call self%get_parameter(self%cQ10ProdVeg,   'cQ10ProdVeg',    '[-]',                 'temperature quotient of production',                                                                      default=1.2_rk)
   call self%get_parameter(self%cQ10RespVeg,   'cQ10RespVeg',    '[-]',                 'temperature quotient of respiration',                                                                     default=2.0_rk)
   call self%get_parameter(self%kMortVegSum,   'kMortVegSum',    'day-1',               'macrophytes mortality rate in Spring and Summer (low)',                                                   default=0.005_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cLengMort,     'cLengMort',      'days',                'length of shoot mort. period',                                                                            default=15.0_rk)
   call self%get_parameter(self%fWinVeg,       'fWinVeg',        '[-]',                 'fraction surviving in winter',                                                                            default=0.3_rk)
   call self%get_parameter(self%cPDVeg0,       'cPDVeg0',        'gP/gD',               'initial P fraction in veg',                                                                               default=0.002_rk)
   call self%get_parameter(self%cNDVeg0,       'cNDVeg0',        'gN/gD',               'initial N fraction in veg',                                                                               default=0.02_rk)
   call self%get_parameter(self%fSedUptVegMax, 'fSedUptVegMax',  '[-]',                 'maximum sediment fraction of nutrient uptake',                                                            default=0.998_rk)
   call self%get_parameter(self%fSedUptVegCoef,'fSedUptVegCoef', '[-]',                 'sigma coef. for sediment fraction of nutrient uptake',                                                    default=2.66_rk)
   call self%get_parameter(self%fSedUptVegExp, 'fSedUptVegExp',  '[-]',                 'exponent in regression for sediment fraction of nutrient uptake',                                         default=-0.83_rk)
   call self%get_parameter(self%cAffNUptVeg,   'cAffNUptVeg',    'l/mgDW/d',            'initial N uptake affinity macrophytes',                                                                   default=0.2_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVNUptMaxVeg,  'cVNUptMaxVeg',   'mgN/mgDW/d',          'maximum N uptake capacity of macrophytes',                                                                default=0.1_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDepthS,       'cDepthS',        '[m]',                 'sediment depth',                                                                                          default=0.1_rk)
   call self%get_parameter(self%bPorS,         'bPorS',          '[m3water/m3sediment]','sediment porosity',                                                                                       default=0.847947_rk)
   call self%get_parameter(self%UseEmpUpt,     'UseEmpUpt',      '[-]',                 'false=do not use this empirical relation',                                                                default=0)
   call self%get_parameter(self%cAffPUptVeg,   'cAffPUptVeg',    'l/mgDW/d',            'initial P uptake affinity macrophytes',                                                                   default=0.2_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVPUptMaxVeg,  'cVPUptMaxVeg',   'mgP/mgDW/d',          'maximum P uptake capacity of macrophytes',                                                                default=0.01_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDissMortVeg,  'fDissMortVeg',   '[-]',                 'fraction dissolved nutrients from dead plants',                                                           default=0.25_rk)
   call self%get_parameter(self%cCPerDW,       'cCPerDW',        'gC/gDW',              'C content of organic matter',                                                                             default=0.4_rk)
   call self%get_parameter(self%fTOMWMortVeg,  'fTOMWMortVeg',   '[-]',                 'fraction of shoot mortality becoming water organic matter',                                               default=0.1_rk)
   call self%get_parameter(self%hO2BOD,        'hO2BOD',         'mgO2/l',              'half-saturation constant for oxygen conc. on BOD',                                                        default=1.0_rk)
   call self%get_parameter(self%cHeightVeg,    'cHeightVeg',     'm',                   'macrophytes height',                                                                                      default=1.0_rk)
   call self%get_parameter(self%cExtSpVeg,     'cExtSpVeg',      'm2/gDW',              'specific extinction of macrophytes',                                                                                     default=0.01_rk)
!  user defined minumum value for state variables
   call self%get_parameter(self%cDVegMin,      'cDVegMin',        'gDW/m2',             'minimum dry weight of macrophytes in system',                                                             default=0.00001_rk)
   call self%get_parameter(self%fVegDOMS,  'fVegDOMS',      '[-]',                'dissolved organic fraction from benthic macrophytes',                                                           default=0.5_rk)
   call self%get_parameter(self%fVegDOMW,  'fVegDOMW',      '[-]',                'dissolved organic fraction from pelagic macrophytes',                                                           default=0.5_rk)
!  Register local state variable
   call self%register_state_variable(self%id_sDVeg,'sDVeg','gDW m-2','macrophytes DW',    &
                                    initial_value=1.0_rk,minimum=self%cDVegMin)
   call self%register_state_variable(self%id_sNVeg,'sNVeg','gN m-2','macrophytes N',     &
                                    initial_value=0.02_rk,minimum=self%cDVegMin * self%cNDVegMin)
   call self%register_state_variable(self%id_sPVeg,'sPVeg','gP m-2','macrophytes P',     &
                                    initial_value=0.002_rk,minimum=self%cDVegMin * self%cPDVegMin)
!  register diagnostic variables
   call self%register_diagnostic_variable(self%id_aDSubVeg,       'aDSubVeg',       'gDW m-2',    'aDSubVeg',                         output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aCovVeg,        'aCovVeg',        '%',          'aCovVeg',                          output=output_instantaneous,domain=domain_bottom)
   call self%register_diagnostic_variable(self%id_afCovSurfVeg,   'afCovSurfVeg',   '[-]',        'afCovSurfVeg',                     output=output_none)
   call self%register_diagnostic_variable(self%id_aDayInitVeg,    'aDayInitVeg',    'd',          'aDayInitVeg',                      output=output_none)
   call self%register_diagnostic_variable(self%id_aNutLimVeg,     'aNutLimVeg',     '[-]',        'nutrient limitation factor',       output=output_instantaneous,domain=domain_bottom)
   call self%register_diagnostic_variable(self%id_macroextinction,'macroextinction','[-]',        'macroextinction',                  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aLLimVeg,       'aLLimVeg',       '[-]',        'light limitation factor',          output=output_instantaneous,domain=domain_bottom)
   call self%register_diagnostic_variable(self%id_aLPAR1Veg,      'aLPAR1Veg',      'W m-2',      'light at top of the macrophytes',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aLPAR2Veg,      'aLPAR2Veg',      'W m-2',      'light at bottom of the macrophytes',output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDBedPOMS,      'tDBedPOMS',      'g m-2 s-1',  'tDBedPOMS',                        output=output_none)
!  register diagnostic variables for modular fluxes
   call self%register_diagnostic_variable(self%id_tDBedVeg,       'tDBedVeg',      'g m-2 s-1',  'macrophytes_DVeg_exchange',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBedVeg,       'tNBedVeg',      'g m-2 s-1',  'macrophytes_NVeg_exchange',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBedVeg,       'tPBedVeg',      'g m-2 s-1',  'macrophytes_PVeg_exchange',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tO2BedW,        'tO2BedW',       'g m-2 s-1',  'macrophytes_tO2BedW_exchange', output=output_instantaneous)
#ifdef _DEVELOPMENT_
   call self%register_diagnostic_variable(self%id_wNBedNH4W,      'wNBedNH4W',     'g m-2 s-1',  'macrophytes_NH4W_exchange',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNBedNO3W,      'wNBedNO3W',     'g m-2 s-1',  'macrophytes_NO3W_exchange',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPBedPO4W,      'wPBedPO4W',     'g m-2 s-1',  'macrophytes_PO4W_exchange',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDBedPOMW,      'wDBedPOMW',     'g m-2 s-1',  'macrophytes_DPOMW_exchange',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNBedPOMW,      'wNBedPOMW',     'g m-2 s-1',  'macrophytes_NPOMW_exchange',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPBedPOMW,      'wPBedPOMW',     'g m-2 s-1',  'macrophytes_PPOMW_exchange',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBedNH4S,      'tNBedNH4S',     'g m-2 s-1',  'macrophytes_NH4S_exchange',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBedNO3S,      'tNBedNO3S',     'g m-2 s-1',  'macrophytes_NO3S_exchange',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBedPO4S,      'tPBedPO4S',     'g m-2 s-1',  'macrophytes_PO4S_exchange',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDBedPOMSflux,  'tDBedPOMSflux', 'g m-2 s-1',  'macrophytes_DPOMS_exchange',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBedPOMS,      'tNBedPOMS',     'g m-2 s-1',  'macrophytes_NPOMS_exchange',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBedPOMS,      'tPBedPOMS',     'g m-2 s-1',  'macrophytes_PPOMS_exchange',   output=output_instantaneous)
#endif
!  register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNVeg)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPVeg)
!  register state variables dependencies
   call self%register_state_dependency(self%id_NH4poolW,     'ammonium_pool_water',              'g m-3', 'ammonium pool for nutrient uptake')
   call self%register_state_dependency(self%id_NO3poolW,     'nitrate_pool_water',               'g m-3', 'nitrate pool for nutrient uptake')
   call self%register_state_dependency(self%id_PO4poolW,     'phosphate_pool_water',             'g m-3', 'phosphate pool for nutrient uptake')
   call self%register_state_dependency(self%id_O2poolW,      'oxygen_pool_water',                'g m-3', 'oxygen pool in water')
   call self%register_state_dependency(self%id_DPOMpoolW,    'POM_DW_pool_water',                'g m-3', 'POM DW pool in water')
   call self%register_state_dependency(self%id_NPOMpoolW,    'POM_N_pool_water',                 'g m-3', 'POM N pool in water')
   call self%register_state_dependency(self%id_PPOMpoolW,    'POM_P_pool_water',                 'g m-3', 'POM P pool in water')
   call self%register_state_dependency(self%id_NH4poolS,     'ammonium_pool_sediment',           'g m-2', 'ammonium pool for nutrient uptake')
   call self%register_state_dependency(self%id_NO3poolS,     'nitrate_pool_sediment',            'g m-2', 'nitrate pool for nutrient uptake')
   call self%register_state_dependency(self%id_PO4poolS,     'phosphate_pool_sediment',          'g m-2', 'phosphate pool for nutrient uptake')
   call self%register_state_dependency(self%id_DPOMpoolS,    'POM_DW_pool_sediment',             'g m-2', 'POM DW pool sediment')
   call self%register_state_dependency(self%id_NPOMpoolS,    'POM_N_pool_sediment',              'g m-2', 'POM N pool sediment')
   call self%register_state_dependency(self%id_PPOMpoolS,    'POM_P_pool_sediment',              'g m-2', 'POM P pool sediment')
   call self%register_state_dependency(self%id_DDOMpoolW,    'DOM_DW_pool_water',                'g m-3', 'DOM DW in water')
   call self%register_state_dependency(self%id_NDOMpoolW,    'DOM_N_pool_water',                 'g m-3', 'DOM N in water')
   call self%register_state_dependency(self%id_PDOMpoolW,    'DOM_P_pool_water',                 'g m-3', 'DOM P in water')
   call self%register_state_dependency(self%id_DDOMpoolS,    'DOM_DW_pool_sediment',             'g m-2', 'DOM DW in sediment')
   call self%register_state_dependency(self%id_NDOMpoolS,    'DOM_N_pool_sediment',              'g m-2', 'DOM N in sediment')
   call self%register_state_dependency(self%id_PDOMpoolS,    'DOM_P_pool_sediment',              'g m-2', 'DOM P in sediment')
!------------------------------------------------------------------------------------------------------------
!  register diagnostic dependencies
!------------------------------------------------------------------------------------------------------------
!  register dependencies on external diagnostic variables
   call self%register_dependency(self%id_afOxySed,        'oxic_layer_fraction',          '[-]',  'oxic layer fraction')
!------------------------------------------------------------------------------------------------------------
!  register environmental dependencies
!
!   call self%register_dependency(self%id_par,    standard_variables%downwelling_photosynthetic_radiative_flux)
!   call self%register_dependency(self%id_meanpar,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk,missing_value=0.0_rk))
!
!------------------------------------------------------------------------------------------------------------
   call self%register_dependency(self%id_uTm,    standard_variables%temperature)
   call self%register_dependency(self%id_extc,   standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
   call self%register_dependency(self%id_Day,    standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)
   call self%register_dependency(self%id_lat,    standard_variables%latitude)
   call self%register_dependency(self%id_dz,     standard_variables%cell_thickness)
   call self%register_dependency(self%id_I_0,    standard_variables%surface_downwelling_shortwave_flux)
!
   !
   call self%register_diagnostic_variable(self%id_par,    'par','W m-2','downwelling_photosynthetic_radiative_flux')
   call self%register_diagnostic_variable(self%id_meanpar, 'meanpar','W m-2','surface_downwelling_photosynthetic_radiative_flux')

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
   class (type_pclake_macrophytes), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
! !LOCAL VARIABLES:
!  Carriers for environment dependencies
   real(rk)     :: uTm,extc,dz,swr
   real(rk)     :: day,sDepthW,lat
!  state variable value carriers
   real(rk)     :: sDVeg,sNVeg,sPVeg
!  external state variable carriers
   real(rk)     :: sNH4W,sNO3W,sPO4W,sNH4S,sNO3S,sPO4S
   real(rk)     :: sO2W
   real(rk)     :: aCorO2BOD,afOxySed
!  variables with o- prefix
   real(rk)     :: oNH4S,oNO3S,oNDissS,oNDissW,oPO4S
!  varaibles for re-allocation of the roots and shoots
   real(rk)     :: bfRootVeg,bfShootVeg,aDRootVeg,aDShootVeg
   real(rk)     :: aDEmergVeg,aDFloatVeg,bfSubVeg,aDSubVeg
   real(rk)     :: bkMortVeg
   real(rk),save :: aDayInitVeg = -1.0_rk
!  variables for coverage of macrophytes (light function and fish feeding)
   real(rk)     :: afCovEmergVeg,aCovVeg
!  temperature function variables
   real(rk)    ::  uFunTmProdVeg,uFunTmRespVeg
!  light function variables
   real(rk)    :: afCovSurfVeg,aLLimShootVeg
   real(rk)    :: uhLVeg,aMuTmLVeg,ufDay,par_bott
   real(rk)    :: aLPAR1Veg,aLPAR2Veg
!  Light function variables
   real(rk)   :: par,fPAR,meanpar,fRefl     !do 过程中要用到的局部变量都需在这里定义
!  nutrient ratio variables
   real(rk)    :: rPDVeg,rNDVeg
!  variables for dry-weight change of macrophytes
   real(rk)    :: tDBedVeg,tDMigrVeg,tDProdVeg,tDRespVeg
   real(rk)    :: aMuVeg,aNutLimVeg,aPLimVeg,aNLimVeg
   real(rk)    :: tDEnvProdVeg,akDIncrVeg,tDEnvVeg
   real(rk)    :: tDMortVeg,tDEnvMortVeg
!  for macrophytes light attenuation
   real(rk)    :: wDBedVeg
!  variables for nitrogen change of macrophytes
   real(rk)    :: tNBedVeg,tNMigrVeg,tNUptVeg
   real(rk)    :: tNUptVegW,tNUptVegS,afNUptVegS,ahNUptVeg
   real(rk)    :: aVNUptMaxCrVeg
   real(rk)    :: aVNUptVegW,aVNUptVegS
   real(rk)    :: tNExcrVeg,tNMortVeg
!  variables for phosphorus change of macrophytes
   real(rk)   :: tPBedVeg,tPMigrVeg,tPUptVeg
   real(rk)   :: tPUptVegW,tPUptVegS,aVPUptMaxCrVeg,ahPUptVeg
   real(rk)   :: aVPUptVegW,aVPUptVegS,afPUptVegS
   real(rk)   :: tPExcrVeg,tPMortVeg
!  variables for external links
!  NH4W & NO3W
   real(rk)   :: wNBedNH4W,tNUptNH4VegW,afNH4UptVegW,tNMortVegNH4W
   real(rk)   :: tNMortVegNH4,wNBedNO3W,tNUptNO3VegW,tNExcrVegS,tNExcrVegW
!  NH4S & NO3S
   real(rk)   :: tNBedNH4S,tNUptNH4VegS,afNH4UptVegS,tNUptNO3VegS
   real(rk)   :: tNMortVegNH4S,tNBedNO3S
!  PO4W&PO4S
   real(rk)   :: wPBedPO4W,tPExcrVegS,tPExcrVegW,tPBedPO4S
   real(rk)   :: tPMortVegPO4W,tPMortVegPO4,tPMortVegPO4S
!  O2W
   real(rk)   :: tO2BedW,tO2ProdVegW,tO2ProdVeg,tO2ProdVegS
   real(rk)   :: tO2RespVegS,tO2RespVegW,tO2UptNO3VegW
!  Orgnics in water,DW,Nitrogen,phosphorus
   real(rk)   :: wDBedTOMW,tDMortVegW
   real(rk)   :: wNBedTOMW,tNMortVegTOMW,tNMortVegTOM
   real(rk)   :: wPBedTOMW,tPMortVegTOMW,tPMortVegTOM
!  organic matter in sediment,DW,Nitrogen,phosphorus
   real(rk)   :: tDBedPOMS,tDMortVegS
   real(rk)   :: tNBedPOMS,tNMortVegS,tNMortVegTOMS
   real(rk)   :: tPBedPOMS,tPMortVegS,tPMortVegTOMS
!  Dissolved organic matter in sediment
   real(rk)   :: tDBedDOMS,tNBedDOMS,tPBedDOMS
   real(rk)   :: wDBedDOMW,wNBedDOMW,wPBedDOMW
   real(rk)   :: wDBedPOMW,wNBedPOMW,wPBedPOMW
   _GET_HORIZONTAL_(self%id_I_0,swr)                   !self%id_I_0，获取self下的标识符id_I_0的值，并赋值给swr （swr表面短波辐射）， _GET_ function可以获取环境变量(environmental dependencies variables)和初始场状态变量（state_variables），这些由yaml文件或物理模型传递给FABM的有具体数值的变量的值。
   if (swr.lt.0._rk) swr=0._rk
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
   _FABM_HORIZONTAL_LOOP_BEGIN_
!-----------------------------------------------------------------------
!  Retrieve current (local) state variable values.
!-----------------------------------------------------------------------
   _GET_HORIZONTAL_(self%id_sDVeg,sDVeg)
   _GET_HORIZONTAL_(self%id_sNVeg,sNVeg)
   _GET_HORIZONTAL_(self%id_sPVeg,sPVeg)
!-----------------------------------------------------------------------
!  Retrieve dependencies
!-----------------------------------------------------------------------
!  Retrieve state dependencies 
   _GET_(self%id_NH4poolW,sNH4W)
   _GET_(self%id_NO3poolW,sNO3W)
   _GET_(self%id_PO4poolW,sPO4W)
   _GET_(self%id_O2poolW,sO2W)
   _GET_HORIZONTAL_(self%id_NH4poolS,sNH4S)
   _GET_HORIZONTAL_(self%id_NO3poolS,sNO3S)
   _GET_HORIZONTAL_(self%id_PO4poolS,sPO4S)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_(self%id_dz,dz)
!  current fabm0d, benthic retrieve meanpar in the center of water column
!  benthic par retrieve par at the bottom of the water column
!  this is up to December 3rd, 2014
   par=swr*self%fPAR
   meanpar=swr*self%fPAR*self%(1.0-fRefl)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_par,par)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_meanpar,meanpar)
   _GET_GLOBAL_(self%id_Day,Day)
   _GET_(self%id_extc,extc)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
   _GET_HORIZONTAL_(self%id_lat,lat)
!  retrieve diagnostic dependency
   _GET_HORIZONTAL_(self%id_afOxySed,afOxySed)
!  auxiliaries for nutrient variables
   oNDissW = sNO3W + sNH4W
   oNDissS = (sNO3S + sNH4S)/self%cDepthS / self%bPorS
   oNH4S=sNH4S/self%cDepthS / self%bPorS
   oNO3S=sNO3S/self%cDepthS / self%bPorS
   oPO4S=sPO4S/self%cDepthS / self%bPorS
!-----------------------------------------------------------------------
!  Current local nutrient ratios in phytoplankton (check the current state)
!-----------------------------------------------------------------------
!  P/D_ratio_of_macrophytes
   rPDVeg=sPVeg/(sDVeg+NearZero)
!  N/D_ratio_of_Diatom
   rNDVeg = sNVeg /(sDVeg+NearZero)
!-----------------------------------------------------------------------
!  temperature functions for macrophytes
!-----------------------------------------------------------------------
!  temperature_function_of_macrophytes_production
   uFunTmProdVeg =  uFunTmVeg(uTm,self%cQ10ProdVeg )
!  temperature_function_of_macrophytes_respiration
   uFunTmRespVeg = uFunTmVeg(uTm,self%cQ10RespVeg)
!-----------------------------------------------------------------------
!  the germination, allocation and reallocation processes
!-----------------------------------------------------------------------
!  Initial_growth_only_once_a_year
!  reset every New Year (careful only NH)
   if (Day .lt. 1._rk) then
      aDayInitVeg = -1._rk
   end if
!  check for vegitation start
   if (uTm >= self%cTmInitVeg .and. aDayInitVeg .eq. -1) then
      aDayInitVeg = Day
   end if
!  cases: northern hemisphere and southern hemisphere

!  setting_root_fraction, shorten the code length
   bfRootVeg = self%fRootVegWin
   if (Day .lt. aDayInitVeg + self%cLengAllo) then
      bfRootVeg = 0.5*(self%fRootVegWin + self%fRootVegSum) + 0.5*(self%fRootVegWin - self%fRootVegSum) * &
      &cos(Pi/self%cLengAllo * (Day - aDayInitVeg))
   else if (Day .lt. self%cDayWinVeg) then
      bfRootVeg = self%fRootVegSum
   else if (Day .lt. self%cDayWinVeg + self%cLengAllo) then
      bfRootVeg = 0.5*(self%fRootVegWin + self%fRootVegSum) - 0.5*(self%fRootVegWin - self%fRootVegSum) * &
      &cos(Pi/self%cLengAllo * (Day - self%cDayWinVeg))
   endif
!  mortality_constant
   if (Day < self%cDayWinVeg) then
   bkMortVeg = self%kMortVegSum
   else if (Day < self%cDayWinVeg + self%cLengMort) then
   bkMortVeg = - log(self%fWinVeg) / self%cLengMort/secs_pr_day
   else
   bkMortVeg = self%kMortVegSum
   endif
!-----------------------------------------------------------------------
!  fractions of roots and shoots
!-----------------------------------------------------------------------
!  shoot_fraction
   bfShootVeg = 1.0_rk - bfRootVeg
!  root_biomass
   aDRootVeg = bfRootVeg * sDVeg
!  shoot_biomass
   aDShootVeg = bfShootVeg * sDVeg
!  emergent_biomass
   aDEmergVeg = self%fEmergVeg * aDShootVeg
!  floating_biomass
   aDFloatVeg = self%fFloatVeg * aDShootVeg
!  submerged_fraction_of_shoot
   bfSubVeg = 1.0_rk - self%fFloatVeg - self%fEmergVeg
!  put in control,telling users
!  (self%fFloatVeg + self%fEmergVeg) should not be more than 1.0
!  submerged_biomass
   aDSubVeg = bfSubVeg * aDShootVeg
!-----------------------------------------------------------------------
!  macrophytes migration
!-----------------------------------------------------------------------
!  migration_flux_dry_weight
   tDMigrVeg = self%kMigrVeg * (self%cDVegIn - sDVeg)
!  net_migration_flux_P
   tPMigrVeg = self%kMigrVeg * (self%cPDVeg0*self%cDVegIn - sPVeg)
!  net_migration_flux_N
   tNMigrVeg = self%kMigrVeg * (self%cNDVeg0*self%cDVegIn - sNVeg)
!-----------------------------------------------------------------------
!  water coverage by macrophytes
!-----------------------------------------------------------------------
!  fraction_of_water_SURFACE_covered_by_plant_species
   afCovSurfVeg = min(1.0_rk, max(aDFloatVeg / (self%cDLayerVeg + NearZero), aDEmergVeg / (&
   &self%fEmergVeg * self%cDCarrVeg + NearZero) ))
!  fraction_emergent_coverage
   afCovEmergVeg = min(1.0_rk, 0.01_rk * self%cCovSpVeg * aDEmergVeg)
!  percent_cover
   aCovVeg = min(100.0_rk, self%cCovSpVeg * aDShootVeg)
!=======================================================================
!  Assimilation section
!=======================================================================
!-----------------------------------------------------------------------
!  Light function of macrophytes  , original PCLake method
!-----------------------------------------------------------------------
!   feh
!   introduced plant height to subtitute the aDpeht1Veg,if fDepth1Veg=0.5, then plant height
!   is sDepthW-sDepthW*fDepth1Veg
!   Since the bottom light can be retrieved, then light on top of the plant can be calculated.
!  half-sat._light_for_macrophytes_production_at_current_temp.
   uhLVeg = self%hLRefVeg * uFunTmProdVeg
   par_bott= meanpar
   aLPAR2Veg = par_bott*exp(- extc *dz/2.0_rk)
   aLPAR1Veg = aLPAR2Veg / exp(- extc * self%cHeightVeg)
!  feh: July 1st, 2016
!  replace sDepthW with plant height, currently the best solution
!   aLLimShootVeg = self%fEmergVeg + self%fFloatVeg * (1.0 - afCovEmergVeg) + bfSubVeg * (1.0 &
!   &- afCovSurfVeg) * 1.0 / (extc * sDepthW) * log( (1.0 + aLPAR1Veg / uhLVeg) /&
!   & (1.0 + aLPAR2Veg / uhLVeg))
   aLLimShootVeg = self%fEmergVeg + self%fFloatVeg * (1.0 - afCovEmergVeg) &
                 + bfSubVeg* (1.0 - afCovSurfVeg) * 1.0 / (extc * self%cHeightVeg) &
                 * log( (1.0 + aLPAR1Veg / uhLVeg) / (1.0 + aLPAR2Veg / uhLVeg))
!=======================================================================
#if 0
   ufDay = 0.5_rk - 0.2_rk * cos(2.0_rk*Pi*Day / 365.0_rk)
#else
   ufDay = DayLength(lat,int(Day))
!   write(*,*) 'lat, Day, ufDay: ',lat,int(Day),ufDay
#endif
!  max._growth_rate_at_current_temp._AND_light
   aMuTmLVeg =ufDay * bfShootVeg * aLLimShootVeg * uFunTmProdVeg * self%cMuMaxVeg
!-----------------------------------------------------------------------
!  Nutrient limitation functions
!-----------------------------------------------------------------------
!  Droop_function_(P)_for_macrophytes
   aPLimVeg = max(0.0_rk, (1.0_rk - self%cPDVegMin / rPDVeg) * self%cPDVegMax / (self%cPDVegMax - self%cPDVegMin) )
!  Droop_function_(N)_for_macrophytes
   aNLimVeg = max(0.0_rk, (1.0_rk - self%cNDVegMin / rNDVeg) * self%cNDVegMax / (self%cNDVegMax - self%cNDVegMin) )
!  nutrient_limitation_function_of_macrophytes
   aNutLimVeg = min( aPLimVeg,aNLimVeg)
!  actual_growth_rate_of_macrophytes
   aMuVeg = aMuTmLVeg * aNutLimVeg
!-----------------------------------------------------------------------
!  macrophytes growth rate correction
!-----------------------------------------------------------------------
!  intrinsic_net_increase_rate_of_macrophytes
   akDIncrVeg = aMuTmLVeg - self%kDRespVeg * uFunTmRespVeg -bkMortVeg
!  logistic_correction_of_macrophytes
   tDEnvVeg = max(0.0_rk, sDVeg**2.0_rk*akDIncrVeg / (self%cDCarrVeg+NearZero))
!  logistic_correction_of_production
   tDEnvProdVeg = aMuVeg / self%cMuMaxVeg * tDEnvVeg
!  macrophytes_production
   tDProdVeg = max(0.0_rk, aMuVeg * sDVeg - tDEnvProdVeg)
!-----------------------------------------------------------------------
!  macrophytes nutrient uptake ---Nitrogen
!-----------------------------------------------------------------------
!  fraction_of_N_uptake_from_sediment
   if (self%UseEmpUpt==0) then
       afNUptVegS = 0.0_rk
   elseif (bfRootVeg <= NearZero) then
      afNUptVegS = 0.0_rk
   else if (self%fFloatVeg + bfSubVeg <= NearZero) then
      afNUptVegS = 1.0_rk
   else
      afNUptVegS = self%fSedUptVegMax / (1.0_rk + self%fSedUptVegCoef * ((((oNDissS+NearZero) / (oN&
      &DissW+NearZero)) )** self%fSedUptVegExp))
   endif
!  fraction_of_P_uptake_from_sediment
   if (self%UseEmpUpt==0) then
      afPUptVegS = 0.0_rk
   elseif (bfRootVeg <= NearZero) then
      afPUptVegS = 0.0_rk
   else if (self%fFloatVeg + bfSubVeg <= NearZero) then
      afPUptVegS = 1.0_rk
   else
      afPUptVegS = self%fSedUptVegMax / (1.0_rk + self%fSedUptVegCoef * ((((oPO4S+NearZero) / (sPO4&
      &W+NearZero)) )** self%fSedUptVegExp))
   endif
!  maximum_P_uptake_rate_of_macrophytes,_corrected_for_P/D_ratio
   aVPUptMaxCrVeg = max( 0.0_rk, self%cVPUptMaxVeg * uFunTmProdVeg * (self%cPDVegMax-rPDVeg) / (&
   &self%cPDVegMax-self%cPDVegMin))
!    P_uptake_RATE_by_subm_AND_floating_parts
   if (self%UseEmpUpt==0) then
      aVPUptVegW = sPO4W * aVPUptMaxCrVeg / (aVPUptMaxCrVeg / self%cAffPUptVeg + sPO4W)
   else
      aVPUptVegW = 0.0_rk
   endif
!    P_uptake_rate_by_roots
   if  (self%UseEmpUpt==0) then
      aVPUptVegS = oPO4S * aVPUptMaxCrVeg / (aVPUptMaxCrVeg / self%cAffPUptVeg + oPO4S)
   else
      aVPUptVegS = 0.0_rk
   endif
!  P_uptake_from_water
   if (self%UseEmpUpt==0) then
      tPUptVegW = aVPUptVegW * (aDSubVeg + aDFloatVeg)
   else
      tPUptVegW = (1.0 - afPUptVegS) * aVPUptMaxCrVeg * sPO4W / (aVPUptMaxCrVeg / self%cAff&
      &PUptVeg + sPO4W) * sDVeg
   endif
! P_uptake_from_pore_water_(by_root_fraction)
   if (self%UseEmpUpt==0) then
      tPUptVegS = aVPUptVegS * aDRootVeg
   else
      tPUptVegS = afPUptVegS * aVPUptMaxCrVeg * oPO4S / (aVPUptMaxCrVeg / self%cAffPUptVeg &
      &+ oPO4S) * sDVeg
   endif

!    total_P_uptake_macrophytes
      tPUptVeg = tPUptVegW + tPUptVegS

!  maximum_N_uptake_rate_of_macrophytes,_corrected_for_N/D_ratio
   aVNUptMaxCrVeg = max( 0.0_rk, self%cVNUptMaxVeg * uFunTmProdVeg * (self%cNDVegMax - rNDVeg) /&
   & (self%cNDVegMax - self%cNDVegMin))
!  half-sat._'constant'_for_N_uptake
   ahNUptVeg = aVNUptMaxCrVeg / self%cAffNUptVeg
!  N_uptake_RATE_by_subm_AND_floating_parts
   if (self%UseEmpUpt==0) then
      aVNUptVegW = oNDissW * aVNUptMaxCrVeg / (ahNUptVeg + oNDissW)
   else
      aVNUptVegW = 0.0_rk
   endif
!  N_uptake_from_water_(by_shoots)
   if (self%UseEmpUpt==0) then
      tNUptVegW = aVNUptVegW * (aDSubVeg + aDFloatVeg)
   else
      tNUptVegW = (1.0_rk - afNUptVegS) * aVNUptMaxCrVeg * oNDissW / (aVNUptMaxCrVeg / self%cA&
      &ffNUptVeg + oNDissW) * sDVeg
   endif
!  N_uptake_RATE_of_roots
   if (self%UseEmpUpt==0) then
      aVNUptVegS = oNDissS * aVNUptMaxCrVeg / (ahNUptVeg + oNDissS)
   else
      aVNUptVegS = 0.0
   endif
!  N_uptake_from_pore_water_(by_roots)
   if (self%UseEmpUpt==0)  then
      tNUptVegS = aVNUptVegS * aDRootVeg
   else
      tNUptVegS = afNUptVegS * aVNUptMaxCrVeg * oNDissS / (aVNUptMaxCrVeg / self%cAffNUptVeg&
      & + oNDissS) * sDVeg
   endif
!  total Nitrogen uptake of water column and sediment
   tNUptVeg = tNUptVegW+ tNUptVegS

!-----------------------------------------------------------------------
!  macrophytes nutrient uptake ---phosphorus
!-----------------------------------------------------------------------
!  half-sat._'constant'_for_P_uptake
   ahPUptVeg = aVPUptMaxCrVeg / self%cAffPUptVeg
!  P_uptake_RATE_by_subm._AND_floating_parts
   aVPUptVegW = sPO4W * aVPUptMaxCrVeg / (ahPUptVeg + sPO4W)
!  P_uptake_rate_by_roots
   aVPUptVegS = oPO4S * aVPUptMaxCrVeg / (ahPUptVeg + oPO4S)
!=======================================================================
!  Dissimilation section
!=======================================================================
!-----------------------------------------------------------------------
!  macrophytes respiration
!-----------------------------------------------------------------------
!  dark_respiration_of_macrophytes
   tDRespVeg = self%kDRespVeg * uFunTmRespVeg * sDVeg
!-----------------------------------------------------------------------
!  macrophytes excretion, nitrogen
!-----------------------------------------------------------------------
#ifdef _V509_
   tNExcrVeg = rNDVeg / (self%cNDVegMin + rNDVeg) * rNDVeg * tDRespVeg
#endif
!  N_excretion_by_macrophytes
   tNExcrVeg = (2.0_rk * rNDVeg) / (self%cNDVegMax + rNDVeg) * rNDVeg * tDRespVeg
!-----------------------------------------------------------------------
!  macrophytes excretion, phosphorus
!-----------------------------------------------------------------------
!  P_excretion_by_macrophytes
   tPExcrVeg = rPDVeg / (self%cPDVegMin + rPDVeg) * rPDVeg * tDRespVeg
!-----------------------------------------------------------------------
!  macrophytes mortality, dry-weight
!-----------------------------------------------------------------------
!  logistic_correction_of_mortality
   tDEnvMortVeg = tDEnvVeg - tDEnvProdVeg
!  total_mortality_flux_DW_macrophytes
   tDMortVeg = bkMortVeg * sDVeg + tDEnvMortVeg
!-----------------------------------------------------------------------
!  macrophytes Mortality, nitrogen
!-----------------------------------------------------------------------
!  N_mortality_flux_of_macrophytes
   tNMortVeg = rNDVeg * tDMortVeg
!-----------------------------------------------------------------------
!  macrophytes Mortality, phosphorus
!-----------------------------------------------------------------------
!  P_mortality_flux_of_macrophytes
   tPMortVeg = rPDVeg * tDMortVeg
!-----------------------------------------------------------------------
!  derivative_of_macrophytes
!-----------------------------------------------------------------------
!  derivative_of_macrophytes_biomass
   tDBedVeg = tDMigrVeg + tDProdVeg- tDMortVeg - tDRespVeg
!  total_macrophytes_N_flux_in_bed_module
   tNBedVeg = tNMigrVeg + tNUptVeg - tNExcrVeg - tNMortVeg
!  total_macrophytes_P_flux_in_bed_module
   tPBedVeg = tPMigrVeg + tPUptVeg - tPExcrVeg - tPMortVeg
!=======================================================================
!  Macrophytes influence on other modules
!=======================================================================
!-----------------------------------------------------------------------
!  Update NH4 in water
!-----------------------------------------------------------------------
!  fraction_ammonium_uptake_from_water_column_(adapted from_WASP_model,_EPA)
   afNH4UptVegW = sNH4W * sNO3W / ((ahNUptVeg + sNH4W) * (ahNUptVeg + sNO3W + NearZero))&
   & + sNH4W * ahNUptVeg / ((sNH4W + sNO3W + NearZero) * (ahNUptVeg + sNO3W + NearZero))
!  NH4_uptake_of_macrophytes_from_water
   tNUptNH4VegW = afNH4UptVegW * tNUptVegW
!  N_excretion_by_macrophytes_to_sediment
   tNExcrVegS = bfRootVeg * tNExcrVeg
!  N_excretion_by_macrophytes_to_water
   tNExcrVegW = tNExcrVeg - tNExcrVegS
!  mortality_flux_of_macrophytes_becoming_dissolved_N
   tNMortVegNH4 = self%fDissMortVeg * tNMortVeg
!  mortality_flux_of_macrophytes_becoming_dissolved_N_in_sediment
   tNMortVegNH4S = bfRootVeg * tNMortVegNH4
!  mortality_flux_of_macrophytes_becoming_dissolved_N_in_water
   tNMortVegNH4W = tNMortVegNH4 - tNMortVegNH4S
!  total_N_flux_from_macrophytes_module_to_NH4_in_water
   wNBedNH4W = - tNUptNH4VegW + tNExcrVegW + tNMortVegNH4W
!-----------------------------------------------------------------------
!  Update NH4 in sediment
!-----------------------------------------------------------------------
!  fraction_ammonium_uptake_from_pore_water_(from_WASP_model,_EPA)
   afNH4UptVegS = oNH4S * oNO3S / ((ahNUptVeg + oNH4S +NearZero) * (ahNUptVeg + oNO&
   &3S +NearZero)) + oNH4S * ahNUptVeg / ((oNH4S + oNO3S+NearZero) * (ahNUptVeg + oN&
   &O3S+NearZero))
!  NH4_uptake_of_macrophytes_from_sediment
   tNUptNH4VegS = afNH4UptVegS * tNUptVegS
!  total_N_flux_from_macrophytes_module_to_NH4_in_pore_water
   tNBedNH4S = - tNUptNH4VegS + tNExcrVegS + tNMortVegNH4S
!-----------------------------------------------------------------------
!  Update NO3 in water
!-----------------------------------------------------------------------
!  NO3_uptake_of_macrophytes_from_water
   tNUptNO3VegW = tNUptVegW - tNUptNH4VegW
!  total_N_flux_from_macrophytes_module_to_NO3_in_water
   wNBedNO3W = - tNUptNO3VegW
!-----------------------------------------------------------------------
!  Update NO3 in sediment
!-----------------------------------------------------------------------
! NO3_uptake_of_macrophytes_from_sediment
   tNUptNO3VegS = tNUptVegS - tNUptNH4VegS
! total_N_flux_from_macrophytes_module_to_NO3_in_pore_water
   tNBedNO3S = - tNUptNO3VegS
!-----------------------------------------------------------------------
!  Update PO4 in water
!-----------------------------------------------------------------------
!  P_excretion_by_macrophytes_in_sediment
   tPExcrVegS = bfRootVeg * tPExcrVeg
!  P_excretion_by_macrophytes_in_water
   tPExcrVegW = tPExcrVeg - tPExcrVegS
!  mortality_flux_of_macrophytes_becoming_dissolved_P
   tPMortVegPO4 = self%fDissMortVeg * tPMortVeg
!  mortality_flux_of_macrophytes_becoming_dissolved_P_in_sediment
   tPMortVegPO4S = bfRootVeg * tPMortVegPO4
!  mortality_flux_of_macrophytes_becoming_dissolved_P_in_water
   tPMortVegPO4W = tPMortVegPO4 - tPMortVegPO4S
!  total_P_flux_from_macrophytes_module_to_PO4_in_water
   wPBedPO4W = - tPUptVegW + tPExcrVegW + tPMortVegPO4W
!-----------------------------------------------------------------------
!  Update PO4 in sediment
!-----------------------------------------------------------------------
! total_P_flux_from_macrophytes_module_to_pore_water_PO4
   tPBedPO4S = - tPUptVegS + tPExcrVegS + tPMortVegPO4S
!-----------------------------------------------------------------------
!  Update O2 in water
!-----------------------------------------------------------------------
!  O2_production_to_water_due_to_NO3_uptake_by_macrophytes
   tO2UptNO3VegW = O2PerNO3 * molO2molN * bfSubVeg * tNUptNO3VegW
!  correction_of_O2_demand_in_water_at_low_oxygen_conc.
   aCorO2BOD = sO2W / (self%hO2BOD + sO2W)
!  submerged_O2_respiration
   tO2RespVegW = molO2molC * self%cCPerDW * bfSubVeg * tDRespVeg * aCorO2BOD
!  root_O2_respiration
   tO2RespVegS = molO2molC * self%cCPerDW * bfRootVeg * tDRespVeg * afOxySed
!  macrophytes_O2_production
   tO2ProdVeg = molO2molC * self%cCPerDW * tDProdVeg
!  O2_transport_to_roots
   tO2ProdVegS = min (tO2RespVegS,tO2ProdVeg)
!  O2_used_for_macrophytes_production
   tO2ProdVegW = min( tO2ProdVeg - tO2ProdVegS, bfSubVeg * tO2ProdVeg)
!  total_water_O2_flux_in_macrophytes_module
   tO2BedW = tO2ProdVegW - tO2RespVegW + tO2UptNO3VegW
!-----------------------------------------------------------------------
!  Update organic matter in water, DW, N and P
!-----------------------------------------------------------------------
!  mortality_flux_becoming_water_organic matter
   tDMortVegW = self%fTOMWMortVeg * (1.0_rk - bfRootVeg) * tDMortVeg
!  total_DW_flux_from_macrophytes_module_to_water_organic matter
   wDBedPOMW = tDMortVegW
   wDBedPOMW = wDBedTOMW * (1.0_rk - self%fVegDOMW)
   wDBedDOMW = wDBedTOMW * self%fVegDOMW
!  mortality_flux_of_macrophytes_becoming_organic_N
   tNMortVegTOM = tNMortVeg - tNMortVegNH4
!  mortality_flux_of_macrophytes_becoming_organic_N_in_water
   tNMortVegTOMW = self%fTOMWMortVeg * (1.0_rk - bfRootVeg) * tNMortVegTOM
!  total_N_flux_from_macrophytes_module_to_water_organic matter
   wNBedTOMW = tNMortVegTOMW
   wNBedPOMW = wNBedTOMW * (1.0_rk - self%fVegDOMW)
   wNBedDOMW = wNBedTOMW * self%fVegDOMW
!  mortality_flux_of_macrophytes_becoming_organic matter_P
   tPMortVegTOM = tPMortVeg - tPMortVegPO4
!  mortality_flux_of_macrophytes_becoming_organic_P_in_water
   tPMortVegTOMW = self%fTOMWMortVeg * (1.0_rk - bfRootVeg) * tPMortVegTOM
!  total_P_flux_from_macrophytes_module_to_water_organic matter
   wPBedTOMW = tPMortVegTOMW
   wPBedPOMW = wPBedTOMW * (1.0_rk - self%fVegDOMW)
   wPBedDOMW = wPBedTOMW * self%fVegDOMW
!---------------------------------------------------------------------------
!  Update organic matter  in sediment, DW, N and P
!---------------------------------------------------------------------------
!  mortality_flux_becoming_sediment_organic matter
   tDMortVegS = tDMortVeg - tDMortVegW
!  total_DW_flux_from_macrophytes_module_to_sediment_organic matter
   tDBedPOMS = tDMortVegS * (1.0_rk - self%fVegDOMS)
   tDBedDOMS = tDMortVegS * self%fVegDOMS
!  mortality_flux_of_macrophytes_becoming_organic_N_in_sediment
   tNMortVegTOMS = tNMortVegTOM - tNMortVegTOMW
!  total_N_flux_from_macrophytes_module_to_sediment_organic matter
   tNBedPOMS = tNMortVegTOMS* (1.0_rk - self%fVegDOMS)
   tNBedDOMS = tNMortVegTOMS * self%fVegDOMS
!  mortality_flux_of_macrophytes_becoming_organic_P_in_sediment
   tPMortVegTOMS = tPMortVegTOM - tPMortVegTOMW
!  total_P_flux_from_macrophytes_module_to_sediment_organic matter
   tPBedPOMS = tPMortVegTOMS * (1.0_rk - self%fVegDOMS)
   tPBedDOMS = tPMortVegTOMS * self%fVegDOMS
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
   _SET_ODE_BEN_(self%id_sDVeg,tDBedVeg)
   _SET_ODE_BEN_(self%id_sNVeg,tNBedVeg)
   _SET_ODE_BEN_(self%id_sPVeg,tPBedVeg)
!-----------------------------------------------------------------------
!  Output local diagnostic variables
!-----------------------------------------------------------------------
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tO2BedW,tO2BedW)
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
   _SET_BOTTOM_EXCHANGE_(self%id_NH4poolW,wNBedNH4W)
   _SET_BOTTOM_EXCHANGE_(self%id_NO3poolW,wNBedNO3W)
   _SET_BOTTOM_EXCHANGE_(self%id_PO4poolW,wPBedPO4W)
   _SET_BOTTOM_EXCHANGE_(self%id_O2poolW,tO2BedW)
   _SET_BOTTOM_EXCHANGE_(self%id_DPOMpoolW,wDBedPOMW)
   _SET_BOTTOM_EXCHANGE_(self%id_NPOMpoolW,wNBedPOMW)
   _SET_BOTTOM_EXCHANGE_(self%id_PPOMpoolW,wPBedPOMW)
   _SET_BOTTOM_EXCHANGE_(self%id_DDOMpoolW,wDBedDOMW)
   _SET_BOTTOM_EXCHANGE_(self%id_NDOMpoolW,wNBedDOMW)
   _SET_BOTTOM_EXCHANGE_(self%id_PDOMpoolW,wPBedDOMW)
!  in the sediment
   _SET_ODE_BEN_(self%id_NH4poolS,tNBedNH4S)
   _SET_ODE_BEN_(self%id_NO3poolS,tNBedNO3S)
   _SET_ODE_BEN_(self%id_PO4poolS,tPBedPO4S)
   _SET_ODE_BEN_(self%id_DPOMpoolS,tDBedPOMS)
   _SET_ODE_BEN_(self%id_NPOMpoolS,tNBedPOMS)
   _SET_ODE_BEN_(self%id_PPOMpoolS,tPBedPOMS)
   _SET_ODE_BEN_(self%id_DDOMpoolS,tDBedDOMS)
   _SET_ODE_BEN_(self%id_NDOMpoolS,tNBedDOMS)
   _SET_ODE_BEN_(self%id_PDOMpoolS,tPBedDOMS)
!-----------------------------------------------------------------------
!  output diagnostic variables for external links
!-----------------------------------------------------------------------
!  Export diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aDSubVeg,aDSubVeg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aCovVeg,aCovVeg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBedPOMS,tDBedPOMS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_afCovSurfVeg,afCovSurfVeg)
!  for macrophytes light attenuation output
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aDayInitVeg,aDayInitVeg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBedVeg,tDBedVeg)
!  light variables output
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_macroextinction,extc)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aLPAR1Veg,aLPAR1Veg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aLPAR2Veg,aLPAR2Veg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aLLimVeg,aLLimShootVeg)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aNutLimVeg,aNutLimVeg)
!  Output diagnostic variables for modular fluxes
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBedVeg,tDBedVeg*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBedVeg,tNBedVeg*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBedVeg,tPBedVeg*secs_pr_day)
#ifdef _DEVELOPMENT_
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wNBedNH4W,wNBedNH4W/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wNBedNO3W,wNBedNO3W/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wPBedPO4W,wPBedPO4W/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tO2BedW,tO2BedW/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wDBedPOMW,wDBedPOMW/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wNBedPOMW,wNBedPOMW/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wPBedPOMW,wPBedPOMW/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBedNH4S,tNBedNH4S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBedNO3S,tNBedNO3S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBedPO4S,tPBedPO4S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBedPOMSflux,tDBedPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBedPOMS,tNBedPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBedPOMS,tPBedPOMS*secs_pr_day)
#endif
   _FABM_HORIZONTAL_LOOP_END_
!
! Spatial loop end
!
!EOP
!-----------------------------------------------------------------------

    end subroutine do_bottom

!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Get the light extinction coefficient due to macrophytes
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_pclake_macrophytes), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s): Fenjuan Hu
!
! !LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  Enter spatial loop
   _LOOP_BEGIN_


!  Self-shading with explicit contribution from background phytoplankton concentration.
!!!UPDATE, FEH, Dec. 14th, 2016
   _SET_EXTINCTION_(self%cExtSpVeg*self%aDSubVeg)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC
!-----------------------------------------------------------------------
!
!EOP
!-----------------------------------------------------------------------

   end module pclake_macrophytes

!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
