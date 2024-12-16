#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module pclake_fish

! !USES:
   use fabm_types
   use pclake_utility, ONLY:uFunTmBio

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_fish
!  local state variable identifiers
!  id_sDFiJv, zooplanktivorous fish concentration in dry-weight, gDW/m**2
!  id_sPFiJv, zooplanktivorous fish concentration in nitrogen element, gN/m**2
!  id_sNFiJv, zooplanktivorous fish concentration in phosphorus element, gP/m**2
!  id_sDFiAd, benthivoros fish concentration in dry-weight, gDW/m**2
!  id_sPFiAd, benthivoros fish concentration in nitrogen element, gN/m**2
!  id_sNFiAd, benthivoros fish concentration in phosphorus element, gP/m**2
!  id_sDPisc, piscivorous fish concentration in dry-weight, gDW/m**2
   type (type_bottom_state_variable_id)            :: id_sDFiJv,id_sPFiJv,id_sNFiJv
   type (type_bottom_state_variable_id)            :: id_sDFiAd,id_sPFiAd,id_sNFiAd,id_sDPisc
!  Diagnostic variables for local output
!  id_aNPisc, piscivorous fish concentration in nitrogen element, gN/m**2
!  id_aPPisc, piscivorous fish concentration in phosphorus element, gP/m**2
   type (type_horizontal_diagnostic_variable_id)       :: id_aNPisc,id_aPPisc
#ifdef _DEVELOPMENT_
! diagnostic variables for modular fluxes
   type (type_horizontal_diagnostic_variable_id)       :: id_tDFiJv,id_tPFiJv,id_tNFiJv
   type (type_horizontal_diagnostic_variable_id)       :: id_tDFiAd,id_tPFiAd,id_tNFiAd
   type (type_horizontal_diagnostic_variable_id)       :: id_tDPisc,id_tPFishPO4W,id_tNFishNH4W
   type (type_horizontal_diagnostic_variable_id)       :: id_tDFishPOMW,id_tNFishPOMW,id_tPFishPOMW
   type (type_horizontal_diagnostic_variable_id)       :: id_tDFishZoo,id_tNFishZoo,id_tPFishZoo
#endif
!  state dependencies identifiers
   type (type_bottom_state_variable_id)                :: id_DFoodBen,id_NFoodBen,id_PFoodBen
   type (type_state_variable_id)                       :: id_DPOMpoolW,id_PPOMpoolW,id_NPOMpoolW
   type (type_state_variable_id)                       :: id_NH4poolW,id_PO4poolW
   type (type_state_variable_id)                       :: id_DFoodZoo,id_NFoodZoo,id_PFoodZoo
   type (type_state_variable_id)                       :: id_DDOMpoolW,id_PDOMpoolW,id_NDOMpoolW
!  environmental dependencies
   type (type_dependency_id)                :: id_uTm
   type (type_dependency_id)                :: id_dz
   type (type_horizontal_dependency_id)     :: id_sDepthW
   type (type_global_dependency_id)         :: id_Day
!  diagnostic dependencies
   type (type_horizontal_dependency_id)     :: id_aDSubVeg !,id_tDEnvFiAd,id_aDSatFiAd
!  for adult fish assimilation
   type ( type_horizontal_dependency_id)    :: id_aCovVeg
!  Model parameters
!  parameters for fish
   real(rk)      :: kMigrFish,cDFiJvIn,cDFiAdIn
   real(rk)      :: cDPiscIn,kMigrPisc,fDBone
   real(rk)      :: fPBone,cDCarrFish,fDissEgesFish,fDissMortFish
   real(rk)      :: cTmOptFish,cSigTmFish,cDayReprFish,fReprFish
   real(rk)      :: fAgeFish,kDAssFiJv,hDZooFiJv,fDAssFiJv
   real(rk)      :: kDRespFiJv,kMortFiJv
   real(rk)      :: kDRespFiAd,kMortFiAd,cDCarrPiscMax,cDCarrPiscMin
   real(rk)      :: cDCarrPiscBare,cDPhraMinPisc,cCovVegMin
   real(rk)      :: cRelPhraPisc,cRelVegPisc,kDAssPisc,hDVegPisc
   real(rk)      :: hDFishPisc,fDAssPisc,fDissEgesPisc,kDRespPisc
   real(rk)      :: kMortPisc,fDissMortPisc,cTmOptPisc,cSigTmPisc
   real(rk)      :: cPDFishRef,cNDFishRef,cPDPisc,cNDPisc
   real(rk)      :: cDayAgeFish
!  Parameter relating to adult fish assimilation
   real(rk)      :: fDAssFiAd, cRelVegFish, kDAssFiAd, hDBentFiAd
!  minimum state variable values
   real(rk)   :: cDFiJvMin,cDFiAdMin,cDPiscMin
!  dissolved organic fraction from fish
   real(rk)   :: fFisDOMW

   contains

!  Model procedures
   procedure :: initialize
   procedure :: do_bottom

   end type type_pclake_fish

!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
   real(rk)           :: Pi=3.14159265358979_rk
!   Lowest state variable value for fish
   real(rk),parameter :: FishZero=0.00001_rk
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
   class (type_pclake_fish), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!EOP
!-----------------------------------------------------------------------
!BOC
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%kMigrFish,       'kMigrFish',      'd-1',       'fish migration rate',                                                            default=0.001_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDFiJvIn,        'cDFiJvIn',       'gDW m-2',   'external fish density',                                                          default=0.005_rk)
   call self%get_parameter(self%cDFiAdIn,        'cDFiAdIn',       'gDW m-2',   'external fish density',                                                          default=0.005_rk)
   call self%get_parameter(self%cDPiscIn,        'cDPiscIn',       'gDW m-2',   'external Pisc. density',                                                         default=0.001_rk)
   call self%get_parameter(self%kMigrPisc,       'kMigrPisc',      'd-1',       'piscivorous migration rate',                                                     default=0.001_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDBone,          'fDBone',         '[-]',       'fraction of fish C fixed in bones and scales',                                   default=0.35_rk)
   call self%get_parameter(self%fPBone,          'fPBone',         '[-]',       'fraction of fish P fixed in bones and scales',                                   default=0.5_rk)
   call self%get_parameter(self%cDCarrFish,      'cDCarrFish',     'gDW m-2',   'carrying capacity of fish',                                                      default=15.0_rk)
   call self%get_parameter(self%fDissEgesFish,   'fDissEgesFish',  '[-]',       'soluble nutrient fraction of fish egested food',                                 default=0.25_rk)
   call self%get_parameter(self%fDissMortFish,   'fDissMortFish',  '[-]',       'soluble nutrient fraction of dead fish (excl bones and scales)',                 default=0.1_rk)
   call self%get_parameter(self%cTmOptFish,      'cTmOptFish',     'degree C',  'optimal temperature of zoo- and benthivorous fish',                                                   default=25.0_rk)
   call self%get_parameter(self%cSigTmFish,      'cSigTmFish',     'degree C',  'temperature constant of fish (sigma in Gaussian curve)',                         default=10.0_rk)
   call self%get_parameter(self%cDayReprFish,    'cDayReprFish',   '[-]',       'reproduction day of year for fish ',                                             default=120.0_rk)
!  new parameter, fish aging day count: the day young fish become adult fish
   call self%get_parameter(self%cDayAgeFish,     'cDayAgeFish',     '[-]',      'aging day of year for fish ',                                                    default=360.0_rk, maximum=364.0_rk)
   call self%get_parameter(self%fReprFish,       'fReprFish',      '[-]',       'yearly reproduction fraction of benthivorous fish, daily rate',                  default=0.02_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fAgeFish,        'fAgeFish',       '[-]',       'yearly ageing fraction of zooplanktivorous fish, daily rate',                    default=0.5_rk,     scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDAssFiJv,       'kDAssFiJv',      'd-1',       'maximum assimilation rate of zooplanktivorous fish',                             default=0.12_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDZooFiJv,       'hDZooFiJv',      'g m-3',     'half-saturation zooplankton for zooplanktivorous fish predation',                default=1.25_rk)
   call self%get_parameter(self%fDAssFiJv,       'fDAssFiJv',      '[-]',       'C assimilation efficiency of zooplanktivorous fish',                             default=0.4_rk)
   call self%get_parameter(self%kDRespFiJv,      'kDRespFiJv',     'd-1',       'maintenance respiration constant of zooplanktivorous fish',                      default=0.01_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortFiJv,       'kMortFiJv',      'd-1',       'specific mortality of zooplanktivorous fish',                                    default=0.00137_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespFiAd,      'kDRespFiAd',     'd-1',       'maintenance respiration constant of benthivorous fish',                          default=0.004_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortFiAd,       'kMortFiAd',      'd-1',       'specific mortality of benthivorous fish',                                        default=0.00027_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDCarrPiscMax,   'cDCarrPiscMax',  'gDW m-2',   'maximum carrying capacity of piscivorous fish',                                  default=1.2_rk)
   call self%get_parameter(self%cDCarrPiscMin,   'cDCarrPiscMin',  'gDW m-2',   'minimum carrying capacity of piscivorous fish',                                  default=0.1_rk)
   call self%get_parameter(self%cDCarrPiscBare,  'cDCarrPiscBare', 'gDW m-2',   'carrying capacity of piscivorous fish',                                          default=0.1_rk)
   call self%get_parameter(self%cDPhraMinPisc,   'cDPhraMinPisc',  'gDW m-2',   'minimum reed biomass for piscivorous fish',                                      default=50.0_rk)
   call self%get_parameter(self%cCovVegMin,      'cCovVegMin',     '%',         'minimum submerged macrophytes coverage for piscivorous fish',                    default=40.0_rk)
   call self%get_parameter(self%cRelPhraPisc,    'cRelPhraPisc',   'gDW m-2',   'relative piscivorous fish density per reed',                                     default=0.075_rk)
   call self%get_parameter(self%cRelVegPisc,     'cRelVegPisc',    'gDW m-2',   'relative piscivorous fish density per reed if aCovVeg>cCovVegMin',               default=0.03_rk)
   call self%get_parameter(self%kDAssPisc,       'kDAssPisc',      'd-1',       'maximum assimilation rate for piscivorous fish',                                 default=0.025_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDVegPisc,       'hDVegPisc',      'g m-2',     'half-saturation constant for macrophytes on piscivorous fish',                   default=5.0_rk)
   call self%get_parameter(self%hDFishPisc,      'hDFishPisc',     'g m-2',     'half-saturating DFish for piscivorous fish predation',                           default=1.0_rk)
   call self%get_parameter(self%fDAssPisc,       'fDAssPisc',      '[-]',       'C assimilation efficiency of piscivorous fish',                                  default=0.4_rk)
   call self%get_parameter(self%fDissEgesPisc,   'fDissEgesPisc',  '[-]',       'soluble P fraction of fish egested food',                                        default=0.25_rk)
   call self%get_parameter(self%kDRespPisc,      'kDRespPisc',     'd-1',       'respiration constant of piscivorous fish',                                       default=0.005_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortPisc,       'kMortPisc',      'd-1',       'specific mortality of piscivorous fish',                                         default=0.00027_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDissMortPisc,   'fDissMortPisc',  '[-]',       'soluble nutrient fraction of dead piscivorous fish',                             default=0.1_rk)
   call self%get_parameter(self%cTmOptPisc,      'cTmOptPisc',     'degree C',  'optimal temperature for piscivorous fish',                                       default=25.0_rk)
   call self%get_parameter(self%cSigTmPisc,      'cSigTmPisc',     'degree C',  'temperature constant for piscivorous fish (sigma)',                              default=10.0_rk)
   call self%get_parameter(self%cPDFishRef,      'cPDFishRef',     'mgP/mgDW',  'reference P/C ratio of fish',                                                    default=0.022_rk)
   call self%get_parameter(self%cNDFishRef,      'cNDFishRef',     'mgN/mgDW',  'reference N/C ratio of fish',                                                    default=0.1_rk)
   call self%get_parameter(self%cPDPisc,         'cPDPisc',        'mgP/mgDW',  'reference P/C ratio of piscivorous fish',                                        default=0.022_rk)
   call self%get_parameter(self%cNDPisc,         'cNDPisc',        'mgN/mgDW',  'reference N/C ratio of piscivorous fish ',                                       default=0.1_rk)
!  The user defined minimum value for state variables
   call self%get_parameter(self%cDFiJvMin,       'cDFiJvMin',       'gDW m-2',   'minimum zooplanktivorous fish biomass in system',                               default=0.0001_rk)
   call self%get_parameter(self%cDFiAdMin,       'cDFiAdMin',       'gDW m-2',   'minimum benthivorous fish biomass in system',                                   default=0.0001_rk)
   call self%get_parameter(self%cDPiscMin,       'cDPiscMin',       'gDW m-2',   'minimum piscivorous fish biomass in system',                                    default=0.0001_rk)
   call self%get_parameter(self%fFisDOMW,        'fFisDOMW',        '[-]',       'dissolved organic matter fraction from fish',                                          default=0.5_rk)
! Parameters regarding adult fish assimilation
   call self%get_parameter(self%fDAssFiAd,    'fDAssFiAd',    '[-]',      'C assimilation efficiency of adult fish',                                              default=0.4_rk)
   call self%get_parameter(self%cRelVegFish,  'cRelVegFish',  '[-]',      'decrease of fish feeding per macrophytes cover (max. 0.01)',                           default=0.009_rk)
   call self%get_parameter(self%kDAssFiAd,    'kDAssFiAd',    'd-1',      'maximum assimilation rate of adult fish',                                              default=0.06_rk,   scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDBentFiAd,   'hDBentFiAd',   'g m-2',    'half-saturation constant for zoobenthos on adult fish',                                default=2.5_rk)
!  Register local state variable
!  zooplanktivorous fish, transportation is turned off
   call self%register_state_variable(self%id_sDFiJv,'sDFiJv','gDW m-2','zooplanktivorous fish DW',     &
                                    initial_value= 0.5_rk,minimum=self%cDFiJvMin) !,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sDFiJv,'disable_transport',.true.)
   call self%register_state_variable(self%id_sPFiJv,'sPFiJv','gP m-2','zooplanktivorous fish P',     &
                                    initial_value=0.011_rk,minimum=self%cDFiJvMin * self%cPDFishRef) !,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sPFiJv,'disable_transport',.true.)
   call self%register_state_variable(self%id_sNFiJv,'sNFiJv','gN m-2','zooplanktivorous fish N',     &
                                    initial_value=0.05_rk,minimum=self%cDFiJvMin * self%cNDFishRef) !,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_spFiJv,'disable_transport',.true.)
!  benthivoros fish, transportation turned off
   call self%register_state_variable(self%id_sDFiAd,'sDFiAd','gDW m-2','benthivorous fish DW',     &
                                    initial_value=2.0_rk,minimum=self%cDFiAdMin) !,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sDFiAd,'disable_transport',.true.)
   call self%register_state_variable(self%id_sPFiAd,'sPFiAd','gP m-2','benthivorous fish P',     &
                                    initial_value=0.044_rk,minimum=self%cDFiAdMin * self%cPDFishRef) !,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sPFiAd,'disable_transport',.true.)
   call self%register_state_variable(self%id_sNFiAd,'sNFiAd','gN m-2','benthivorous fish N',     &
                                    initial_value=0.2_rk,minimum=self%cDFiAdMin * self%cNDFishRef) !,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sNFiAd,'disable_transport',.true.)
!  piscivorous fish
   call self%register_state_variable(self%id_sDPisc,'sDPisc','gDW m-2','piscivorous fish DW', &
                                    initial_value=0.01_rk,minimum=NearZero) !,no_river_dilution=.TRUE.)
!   call self%set_variable_property(self%id_sDPisc,'disable_transport',.true.)

!  Register diagnostic variables for dependencies in other modules
   call self%register_diagnostic_variable(self%id_aNPisc,    'aNPisc',      'g m-3',     'piscivorous fish nitrogen content',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aPPisc,    'aPPisc',      'g m-3',     'piscivorous fish phosphorus content', output=output_instantaneous)
#ifdef _DEVELOPMENT_
!  register diagnostic variables for modular fluxes
   call self%register_diagnostic_variable(self%id_tDFiJv,     'tDFiJv',     'g m-3 s-1', 'fish_DFiJv_change',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPFiJv,     'tPFiJv',     'g m-3 s-1', 'fish_PFiJv_change',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNFiJv,     'tNFiJv',     'g m-3 s-1', 'fish_NFiJv_change',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDFiAd,     'tDFiAd',     'g m-3 s-1', 'fish_DFiAd_change',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPFiAd,     'tPFiAd',     'g m-3 s-1', 'fish_PFiAd_change',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNFiAd,     'tNFiAd',     'g m-3 s-1', 'fish_NFiAd_change',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDPisc,     'tDPisc',     'g m-3 s-1', 'fish_DPisc_change',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNFishNH4W, 'tNFishNH4W', 'g m-3 s-1', 'fish_NH4W_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPFishPO4W, 'tPFishPO4W', 'g m-3 s-1', 'fish_PO4W_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDFishPOMW, 'tDFishPOMW', 'g m-3 s-1', 'fish_DPOMW_change',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNFishPOMW, 'tNFishPOMW', 'g m-3 s-1', 'fish_NPOMW_change',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPFishPOMW, 'tPFishPOMW', 'g m-3 s-1', 'fish_PPOMW_change',                   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDFishZoo,  'tDFishZoo',  'g m-3 s-1', 'fish_DZoo_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNFishZoo,  'tNFishZoo',  'g m-3 s-1', 'fish_NZoo_change',                    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPFishZoo,  'tPFishZoo',  'g m-3 s-1', 'fish_PZoo_change',                    output=output_instantaneous)
#endif
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNFiJv)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNFiAd)
!   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_aNPisc)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPFiJv)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPFiAd)
!   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_aPPisc)
!  register state variables dependencies
   call self%register_state_dependency(self%id_DPOMpoolW,    'POM_DW_pool_water',     'g m-3', 'POM DW pool in water')
   call self%register_state_dependency(self%id_NPOMpoolW,    'POM_N_pool_water',      'g m-3', 'POM N pool in water')
   call self%register_state_dependency(self%id_PPOMpoolW,    'POM_P_pool_water',      'g m-3', 'POM P pool in water')
   call self%register_state_dependency(self%id_NH4poolW,     'NH4_pool_water',        'g m-3', 'NH4 pool in water')
   call self%register_state_dependency(self%id_PO4poolW,     'PO4_pool_water',        'g m-3', 'PO4 pool in water')
   call self%register_state_dependency(self%id_DFoodZoo,     'zooplankton_D_Food',    'g m-3', 'zooplankton D Food')
   call self%register_state_dependency(self%id_NFoodZoo,     'zooplankton_N_Food',    'g m-3', 'zooplankton N Food')
   call self%register_state_dependency(self%id_PFoodZoo,     'zooplankton_P_Food',    'g m-3', 'zooplankton P Food')
   call self%register_state_dependency(self%id_DDOMpoolW,    'DOM_DW_pool_water',     'g m-3', 'DOM DW in water')
   call self%register_state_dependency(self%id_NDOMpoolW,    'DOM_N_pool_water',      'g m-3', 'DOM N in water')
   call self%register_state_dependency(self%id_PDOMpoolW,    'DOM_P_pool_water',      'g m-3', 'DOM P in water')
   call self%register_state_dependency(self%id_DFoodBen,     'zoobenthos_D_Food',    'g m-3', 'zoobenthos D Food')
   call self%register_state_dependency(self%id_NFoodBen,     'zoobenthos_N_Food',    'g m-3', 'zoobenthos N Food')
   call self%register_state_dependency(self%id_PFoodBen,     'zoobenthos_P_Food',    'g m-3', 'zoobenthos P Food')
   
!  register environmental dependencies
   call self%register_dependency(self%id_uTm,    standard_variables%temperature)
   call self%register_dependency(self%id_Day,    standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)
   call self%register_dependency(self%id_dz,     standard_variables%cell_thickness)
   call self%register_dependency(self%id_aCovVeg,   'macrophytes_coverage',       '[-]',  'macrophytes coverage')
   call self%register_dependency(self%id_aDSubVeg,  'submerged_macrophytes',      'g m-2','submerged macrophytes dry weight')
   
   return

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP

! !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
!  INPUT PARAMETERS:
   class (type_pclake_fish), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!  LOCAL VARIABLES:
!  Carriers for environment dependencies
   real(rk)     :: uTm,Day,sDepthW,dz
!  carriers for local state variables
   real(rk)      :: sDFiJv,sPFiJv,sNFiJv
   real(rk)      :: sDFiAd,sPFiAd,sNFiAd,sDPisc
!  carriers for external link state variables
   real(rk)      :: sDZoo,sNZoo,sPZoo
   real(rk)      :: sDBent,sNBent,sPBent
!  carriers for external link diagnostic variables
   real(rk)      :: aDSubVeg,tDEnvFiAd,aDSatFiAd
!  nutrient ratios variables
   real(rk)      :: rPDFiJv,rNDFiJv,rPDFiAd,rNDFiAd
   real(rk)      :: rPDZoo,rNDZoo,rPDBent,rNDBent
!  status auxiliaries
   real(rk)      :: aDFish,aPFish,aNFish
!  variables for temperature function
   real(rk)      :: uFunTmFish,uFunTmPisc
!  variables for fish (including Jv, juvinile, and Ad, adult)
   real(rk)      :: tDReprFish,tDAgeFish
   real(rk)      :: tPReprFish,tPAgeFish
   real(rk)      :: tNReprFish,tNAgeFish
   real(rk)      :: tDFiJv,tDMigrFiJv
   real(rk)      :: tDAssFiJv,tDRespFiJv,tDMortFiJv,tDConsFiJvPisc
   real(rk)      :: aDSatFiJv,tDEnvFiJv,ukDIncrFiJv,tDConsFiJv
!  variables for zooplanktivorous fish flux_P
   real(rk)      :: wPFiJv
   real(rk)      :: tPFiJv,tPMigrFiJv,tPAssFiJv
   real(rk)      :: tPExcrFiJv,tPMortFiJv,tPConsFiJvPisc
   real(rk)      :: afPAssFiJv,tPConsFiJv,afNAssFiJv,tNConsFiJv
   real(rk)      :: tNFiJv,tNMigrFiJv,tNAssFiJv
   real(rk)      :: tNExcrFiJv,tNMortFiJv,tNConsFiJvPisc
!  variables for benthivorous fish
   real(rk)     :: tDFiAd,tDMigrFiAd,tDRespFiAd,tDMortFiAd
   real(rk)     :: tDConsFiAdPisc
   real(rk)     :: tPFiAd,tPMigrFiAd,tPExcrFiAd,tPMortFiAd
   real(rk)     :: tPConsFiAdPisc
!  assimilation
!  variables for benthivorous fish flux_N
   real(rk)     :: tNFiAd,tNMigrFiAd,tNExcrFiAd,tNMortFiAd
   real(rk)     :: tNConsFiAdPisc
   real(rk)     :: tDConsPisc,tDAssPisc,aDSatPisc,aFunVegPisc
   real(rk)     :: tDEnvPisc,akDIncrPisc,aDCarrPisc
   real(rk)     :: tDMigrPisc,tDRespPisc,tDMortPisc
   real(rk)     :: tDPisc
!  variables for piscivorous fish P process
   real(rk)     :: aPPisc,tPConsPisc,rPDFoodPisc,afPAssPisc,tPAssPisc
   real(rk)     :: tPEgesPisc,tPExcrPisc,tPMortPisc,tPMigrPisc
!  variables for piscivorous fish N process
   real(rk)     :: aNPisc,tNConsPisc,rNDFoodPisc,afNAssPisc,tNAssPisc
   real(rk)     :: tNEgesPisc,tNExcrPisc,tNMortPisc,tNMigrPisc
!  variables for exchange of NH4
   real(rk)     :: tNEgesFiJv,tNMortFishNH4,tNMortFishBot,tNEgesFiJvNH4
   real(rk)     :: tNMortFish,tNEgesPiscNH4,tNMortPiscNH4,tNMortPiscBot
   real(rk)     :: tNFishNH4W
!  variables for exchange of PO4
   real(rk)     :: tPEgesFiJvPO4,tPFishPO4W
   real(rk)     :: tPEgesFiJv,tPMortFish,tPMortFishBot
   real(rk)     :: tPMortFishPO4,tPEgesPiscPO4,tPMortPiscPO4,tPMortPiscBot
!  variables for exchange of organic matter DW
   real(rk)     :: tDEgesFiJv, tDMortFishTOM
   real(rk)     :: tDMortFish,tDMortFishBot,tDEgesPisc,tDMortPiscTOM,tDMortPiscBot
!  variables for exchange of organic matter N
   real(rk)     :: tNEgesFiJvTOM,tNMortFishPOM
   real(rk)     :: tNEgesPiscTOM,tNMortPiscTOM
!  variables for exchange of organic matter P
   real(rk)     :: tPEgesFiJvTOM,tPMortFishTOM
   real(rk)     :: tPEgesPiscTOM,tPMortPiscTOM
!  benthivorous fish assimilation
   real(rk)     :: ukDIncrFiAd
!  Fish organic compoments
   real(rk)     :: tDFishTOMW, tNFishTOMW, tPFishTOMW
   real(rk)     :: tDFishPOMW, tNFishPOMW, tPFishPOMW
   real(rk)     :: tDFishDOMW, tNFishDOMW, tPFishDOMW
!  adult fish assimilation
   real(rk)     :: tDAssFiAd,aFunVegFish
   real(rk)     :: aCovVeg
   real(rk)     :: tDConsFiAd,tNConsFiAd,tPConsFiAd
   real(rk)     :: afNAssFiAd,tNAssFiAd,afPAssFiAd,tPAssFiAd
   real(rk)     :: tDEgesFiAd,tNEgesFiAd,tPEgesFiAd
   real(rk)     :: tNEgesFiAdNH4,tPEgesFiAdPO4,tNEgesFiAdTOM,tPEgesFiAdTOM

#ifdef _DEVELOPMENT_
   integer, save :: n=0
#endif
!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
!   _LOOP_BEGIN_
   _FABM_HORIZONTAL_LOOP_BEGIN_
!  Retrieve current (local) state variable values.
   _GET_HORIZONTAL_(self%id_sDFiJv,sDFiJv)
   _GET_HORIZONTAL_(self%id_sPFiJv,sPFiJv)
   _GET_HORIZONTAL_(self%id_sNFiJv,sNFiJv)
   _GET_HORIZONTAL_(self%id_sDFiAd,sDFiAd)
   _GET_HORIZONTAL_(self%id_sPFiAd,sPFiAd)
   _GET_HORIZONTAL_(self%id_sNFiAd,sNFiAd)
   _GET_HORIZONTAL_(self%id_sDPisc,sDPisc)
!-----------------------------------------------------------------------
!  Retrieve dependencies  value
!-----------------------------------------------------------------------
   _GET_(self%id_DFoodZoo,sDZoo)
   _GET_(self%id_NFoodZoo,sNZoo)
   _GET_(self%id_PFoodZoo,sPZoo)
   _GET_HORIZONTAL_(self%id_DFoodBen,sDBent)
   _GET_HORIZONTAL_(self%id_NFoodBen,sNBent)
   _GET_HORIZONTAL_(self%id_PFoodBen,sPBent)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_GLOBAL_(self%id_Day,Day)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
   _GET_(self%id_dz,dz)
! !retrieve diagnostic dependency
   _GET_HORIZONTAL_(self%id_aDSubVeg,aDSubVeg)
   _GET_HORIZONTAL_(self%id_aCovVeg,aCovVeg)

!-------------------------------------------------------------------------
!  The following section is organized in the order of zooplanktivorous fish, 
!  benthivorous fish, and piscivorous processes. There is a separate section 
!  on zooplankton consumption by zooplanktivorous fish, due to the variables 
!  being dependent on each other. It starts with assimilation of DW to provide 
!  wDAssFiJv for its predation of zooplankton(DW,N,P).The latter process (predation) 
!  provide the variable wDConsFiJv for fish assimilation of N,P
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!  Current local nutrients ratios in zooplankton (check the current state)
!-------------------------------------------------------------------------
!  P/D_ratio_herb.zooplankton
   rPDZoo = sPZoo /(sDZoo+NearZero)
!  N/C_ratio_herb.zooplankton
   rNDZoo = sNZoo/(sDZoo+NearZero)
!-------------------------------------------------------------------------
!  Current local nutrients ratios in zoobenthos(check the current state)
!-------------------------------------------------------------------------
!  P/D_ratio_herb.zooplankton
   rPDBent = sPBent /(sDBent+NearZero)
!  N/C_ratio_herb.zooplankton
   rNDBent = sNBent/(sDBent+NearZero)
!-----------------------------------------------------------------------
!  status auxiliaries---auxiliaries for describing the current status,
!  usually derivatives of state variables
!-----------------------------------------------------------------------
!  P/D_ratio_of_young_fish
   rPDFiJv = sPFiJv /(sDFiJv+NearZero)
!  P/D_ratio_of_adult_fish
   rPDFiAd = sPFiAd /(sDFiAd+NearZero)
!  N/D_ratio_of_young_fish
   rNDFiJv = sNFiJv /(sDFiJv+NearZero)
!  N/D_ratio_of_adult_fish
   rNDFiAd = sNFiAd /(sDFiAd+NearZero)
!  total_fish_biomass
   aDFish = sDFiJv + sDFiAd
!  total_fish_biomass
   aPFish = sPFiJv + sPFiAd
!  total_fish_biomass
   aNFish = sNFiJv + sNFiAd
!-----------------------------------------------------------------------
!  temperature function
!-----------------------------------------------------------------------
!  temp._function_of_fish
   uFunTmFish = uFunTmBio(uTm,self%cSigTmFish,self%cTmOptFish)
!  temp._function_of_Pisc
   uFunTmPisc = uFunTmBio(uTm,self%cSigTmPisc,self%cTmOptPisc)
!-----------------------------------------------------------------------
!  zooplanktivorous fish assimilation_DW
!-----------------------------------------------------------------------
!  food_limitation_function_of_young_fish
   aDSatFiJv = (sDZoo * sDZoo) /(self%hDZooFiJv * self%hDZooFiJv + sDZoo * sDZoo)
!  intrinsic_net_increase_rate_of_fish
   ukDIncrFiJv = (self%kDAssFiJv - self%kDRespFiJv) * uFunTmFish - self%kMortFiJv
!  environmental_correction_of_fish
   tDEnvFiJv = max(0.0_rk,ukDIncrFiJv /(self%cDCarrFish - sDFiAd) * sDFiJv*sDFiJv)
!  assimilation_of_fish
   tDAssFiJv = aDSatFiJv *(self%kDAssFiJv * uFunTmFish * sDFiJv - tDEnvFiJv)
!-----------------------------------------------------------------------
!  zooplankton predated by fish
!-----------------------------------------------------------------------
!  zooplankton_consumption_of_fish
    tDConsFiJv = tDAssFiJv / self%fDAssFiJv
!  (zooplankton)_P_consumption_by_FiJv
   tPConsFiJv = rPDZoo * tDConsFiJv
!  (zooplankton)_N_consumption_by_FiJv
   tNConsFiJv = rNDZoo * tDConsFiJv
!-----------------------------------------------------------------------
!  zooplanktivorous fish assimilation_P
!-----------------------------------------------------------------------
!  P_assim._efficiency_of_FiJv
   afPAssFiJv = min(1.0_rk,self%cPDFishRef / rPDZoo * self%fDAssFiJv)
!  P_assimilation_of_FiJv
   tPAssFiJv = afPAssFiJv * tPConsFiJv
!-----------------------------------------------------------------------
!  zooplanktivorous fish assimilation_N
!-----------------------------------------------------------------------
!  N_assim._efficiency_of_FiJv
   afNAssFiJv = min(1.0_rk,self%cNDFishRef / rNDZoo * self%fDAssFiJv)
!  N_assimilation_of_FiJv
   tNAssFiJv = afNAssFiJv * tNConsFiJv
!-----------------------------------------------------------------------
!  zooplanktivorous fish respiration and excretion
!-----------------------------------------------------------------------
!  respiration_of_fish_DW
   tDRespFiJv = (self%cPDFishRef / rPDFiJv) * self%kDRespFiJv * uFunTmFish * sDFiJv
!  P_excretion_of_FiJv
   tPExcrFiJv = (rPDFiJv / self%cPDFishRef) * self%kDRespFiJv * uFunTmFish * sPFiJv
!  N_excretion_of_FiJv
   tNExcrFiJv = (rNDFiJv / self%cNDFishRef) * self%kDRespFiJv * uFunTmFish * sNFiJv
!-----------------------------------------------------------------------
!  zooplanktivorous fish migration
!-----------------------------------------------------------------------
!  migration_flux of zooplanktivorous fish, DW
   tDMigrFiJv = self%kMigrFish *(self%cDFiJvIn - sDFiJv)
!  net_migration_flux of zooplanktivorous fish,P
   tPMigrFiJv = self%kMigrFish *(self%cPDFishRef * self%cDFiJvIn - sPFiJv)
!  net_migration_flux of zooplanktivorous fish,N
   tNMigrFiJv = self%kMigrFish *(self%cNDFishRef * self%cDFiJvIn - sNFiJv)
!-----------------------------------------------------------------------
!  zooplanktivorous fish mortality
!-----------------------------------------------------------------------
!  fish_mortality_incl._environmental_correction
   tDMortFiJv = self%kMortFiJv * sDFiJv +(1.0_rk - aDSatFiJv) * tDEnvFiJv
!  mortality_of_FiJv_P
   tPMortFiJv = rPDFiJv * tDMortFiJv
!  mortality_of_FiJv_N
   tNMortFiJv = rNDFiJv * tDMortFiJv
!-----------------------------------------------------------------------
!  zooplanktivorous fish egestion
!-----------------------------------------------------------------------
!  egestion_of_fish,zooplanktivorous fish
   tDEgesFiJv = tDConsFiJv - tDAssFiJv
!  egestion_of_FiJv
   tNEgesFiJv = tNConsFiJv - tNAssFiJv
!  egestion_of_FiJv
   tPEgesFiJv = tPConsFiJv - tPAssFiJv
!-----------------------------------------------------------------------
!  benthivorous fish assimilation_DW
!-----------------------------------------------------------------------
!  macrophytes_dependence_of_fish_feeding
   aFunVegFish = max(0.0_rk,1.0_rk - self%cRelVegFish * aCovVeg)
!   for first time step check out
!   aFunVegFish = max(0.0_rk,1.0_rk - self%cRelVegFish * 0.2_rk)
!  food_limitation_function_of_adult_fish
   aDSatFiAd = (aFunVegFish * sDBent) *(aFunVegFish * sDBent) /(self%hDBentFiAd * &
   &self%hDBentFiAd + (aFunVegFish * sDBent) *(aFunVegFish * sDBent))
!  intrinsic_net_increase_rate_of_fish
   ukDIncrFiAd = (self%kDAssFiAd - self%kDRespFiAd) * uFunTmFish - self%kMortFiAd
!  environmental_correction_of_fish,in concentration
   tDEnvFiAd = max(0.0_rk,ukDIncrFiAd /(self%cDCarrFish - sDFiJv) * sDFiAd*sDFiAd)
!  assimilation_of_fish
   tDAssFiAd = aDSatFiAd *(self%kDAssFiAd * uFunTmFish * sDFiAd - tDEnvFiAd)
!-----------------------------------------------------------------------
!  benthivorous fish assimilation_P
!-----------------------------------------------------------------------
!  zoobenthos_consumption_of_fish
   tDConsFiAd = tDAssFiAd / self%fDAssFiAd
!  (zoobenthos)_P_consumption_by_FiAd
   tPConsFiAd = rPDBent * tDConsFiAd
!  P_assim._efficiency_of_FiAd
   afPAssFiAd = min(1.0_rk,self%cPDFishRef / rPDBent * self%fDAssFiAd)
!  P_assimilation_of_FiAd
   tPAssFiAd = afPAssFiAd * tPConsFiAd
!-----------------------------------------------------------------------
!  benthivorous fish assimilation_N
!-----------------------------------------------------------------------
!  Zoobenthos_N_consumption_by_FiAd
   tNConsFiAd = rNDBent * tDConsFiAd
!  N_assimilation_efficiency_of_FiAd
   afNAssFiAd = min(1.0_rk,self%cNDFishRef / rNDBent * self%fDAssFiAd)
!  N_assimilation_of_FiAd
   tNAssFiAd = afNAssFiAd * tNConsFiAd
!-----------------------------------------------------------------------
!  benthivorous fish respiration and excretion
!-----------------------------------------------------------------------
!  respiration_of_fish
   tDRespFiAd = (self%cPDFishRef / rPDFiAd) * self%kDRespFiAd * uFunTmFish * sDFiAd
!  P_excretion_of_FiAd
   tPExcrFiAd = (rPDFiAd / self%cPDFishRef) * self%kDRespFiAd * uFunTmFish * sPFiAd
!  N_excretion_of_FiAd
   tNExcrFiAd = (rNDFiAd / self%cNDFishRef) * self%kDRespFiAd * uFunTmFish * sNFiAd
!-----------------------------------------------------------------------
!  benthivorous fish migration
!-----------------------------------------------------------------------
!  migration_flux of benthivorous fish,DW
   tDMigrFiAd = self%kMigrFish *(self%cDFiAdIn - sDFiAd)
!  net_migration_flux of benthivorous fish, P
   tPMigrFiAd = self%kMigrFish *(self%cPDFishRef * self%cDFiAdIn - sPFiAd)
!  net_migration_flux of benthivorous fish, N
   tNMigrFiAd = self%kMigrFish *(self%cNDFishRef * self%cDFiAdIn - sNFiAd)
!-----------------------------------------------------------------------
!  benthivorous fish mortality
!-----------------------------------------------------------------------
!  fish_mortality_incl._environmental_correction
   tDMortFiAd = self%kMortFiAd * sDFiAd +(1.0_rk - aDSatFiAd) * tDEnvFiAd
!  mortality_of_FiAd
   tPMortFiAd = rPDFiAd * tDMortFiAd
!  mortality_of_FiAd
   tNMortFiAd = rNDFiAd * tDMortFiAd
!-----------------------------------------------------------------------
!  benthivirous fish egestion
!-----------------------------------------------------------------------
!  egestion_of_fish,adult fish
   tDEgesFiAd = tDConsFiAd - tDAssFiAd
!  egestion_of_FiAd
   tPEgesFiAd = tPConsFiAd - tPAssFiAd
!  egestion_of_FiAd
   tNEgesFiAd = tNConsFiAd - tNAssFiAd
!-----------------------------------------------------------------------
!  fish reproduction
!-----------------------------------------------------------------------
!  Reproduction_flux_DW
   if (Day >= self%cDayReprFish .AND. Day <= self%cDayReprFish + 1.0_rk) then
      tDReprFish = self%fReprFish * sDFiAd
   else
      tDReprFish =0.0_rk
   endif
!  Reproduction_flux_P
   tPReprFish = rPDFiAd * tDReprFish
!  Reproduction_flux_N
   tNReprFish = rNDFiAd * tDReprFish
!-----------------------------------------------------------------------
!  fish aging
!-----------------------------------------------------------------------
!  Ageing_DW
   if (Day >= self%cDayAgeFish .AND. Day <= self%cDayAgeFish + 1.0_rk) then
      tDAgeFish = self%fAgeFish * sDFiJv
   else
      tDAgeFish = 0.0_rk
   endif
!  Ageing_P
   tPAgeFish = rPDFiJv * tDAgeFish
!  Ageing_N
   tNAgeFish = rNDFiJv * tDAgeFish
!---------------------------------------------------------------------------
!  Piscivorous fish assimilation 
!---------------------------------------------------------------------------
!  macrophytes_dependence_of_Pisc_growth_rate
   aFunVegPisc = aDSubVeg /(self%hDVegPisc + aDSubVeg + NearZero)
!  food_limitation_function_of_Pisc
   aDSatPisc = aDFish*aDFish /(self%hDFishPisc*self%hDFishPisc + aDFish*aDFish)
!  intrinsic_net_increase_rate_of_Pisc
   akDIncrPisc = (self%kDAssPisc * aFunVegPisc - self%kDRespPisc) * uFunTmPisc - self%kMortPisc
!  Carrying_capacity_of_Pisc_for_lake_without_marsh_zone
   aDCarrPisc = max(self%cDCarrPiscMin,min(self%cDCarrPiscMax,self%cDCarrPiscBare))
!  environmental_correction_of_Pisc
   tDEnvPisc = max(0.0_rk,akDIncrPisc / aDCarrPisc * sDPisc*sDPisc)
!  assimilation_of_Pisc
   tDAssPisc = aDSatPisc *(self%kDAssPisc * aFunVegPisc * uFunTmPisc * sDPisc - tDEnvPisc)
!-----------------------------------------------------------------------
!  Piscivorous fish consumption
!-----------------------------------------------------------------------
!  consumption_of_Pisc
   tDConsPisc = tDAssPisc / self%fDAssPisc
!-----------------------------------------------------------------------
!  zooplanktivorous fish predated by piscivorous fish
!-----------------------------------------------------------------------
!  young_fish_consumption_by_Pisc_DW
   tDConsFiJvPisc = sDFiJv / aDFish * tDConsPisc
!  young_fish_consumption_by_Pisc
   tPConsFiJvPisc = rPDFiJv * tDConsFiJvPisc
!  young_fish_consumption_by_Pisc
   tNConsFiJvPisc = rNDFiJv * tDConsFiJvPisc
!-----------------------------------------------------------------------
!  benthivorous fish predated by piscivorious fish
!-----------------------------------------------------------------------
!  adult_fish_consumption_by_Pisc
   tDConsFiAdPisc = tDConsPisc - tDConsFiJvPisc
!  adult_fish_consumption_by_Pisc
   tPConsFiAdPisc = rPDFiAd * tDConsFiAdPisc
!  adult_fish_consumption_by_Pisc
   tNConsFiAdPisc = rNDFiAd * tDConsFiAdPisc
!-----------------------------------------------------------------------
!  piscivorous fish migration
!-----------------------------------------------------------------------
!  migration_flux
   tDMigrPisc = self%kMigrPisc *(self%cDPiscIn - sDPisc)
!-----------------------------------------------------------------------
!  piscivorous fish respiration
!-----------------------------------------------------------------------
!  respiration_of_Pisc
   tDRespPisc = self%kDRespPisc * uFunTmPisc * sDPisc
!-----------------------------------------------------------------------
!  piscivorous fish mortality
!-----------------------------------------------------------------------
!  mortality_of_Pisc(incl._environmental_correction)
   tDMortPisc = self%kMortPisc * sDPisc +(1.0_rk - aDSatPisc) * tDEnvPisc
!---------------------------------------------------------------------------
!  piscivorous fish P process
!---------------------------------------------------------------------------
!  _Piscivorous_fish
    aPPisc = self%cPDPisc * sDPisc
!  total_P_consumption_by_Pisc
   tPConsPisc = tPConsFiJvPisc + tPConsFiAdPisc
!  average_P/D_ratio_of_Pisc_food
   rPDFoodPisc = tPConsPisc / tDConsPisc
!  P_assim._efficiency_of_Pisc
    afPAssPisc = min(1.0_rk,self%cPDPisc / rPDFoodPisc * self%fDAssPisc)
!  P_assimilation_of_Pisc
   tPAssPisc = afPAssPisc * tPConsPisc
!  respiration_of_Pisc
   tPExcrPisc = self%cPDPisc * tDRespPisc
!  mortality_of_Pisc
   tPMortPisc = self%cPDPisc * tDMortPisc
!  net_migration_flux
   tPMigrPisc = self%kMigrPisc *(self%cPDPisc * self%cDPiscIn - aPPisc)
!-----------------------------------------------------------------------
!  piscivorous fish N process
!-----------------------------------------------------------------------
!  Piscivorous_fish
   aNPisc = self%cNDPisc * sDPisc
!  total_N_consumption_by_Pisc
   tNConsPisc = tNConsFiJvPisc + tNConsFiAdPisc
!  average_N/D_ratio_of_Pisc_food
   rNDFoodPisc = tNConsPisc / tDConsPisc
!  N_assim._efficiency_of_Pisc
    afNAssPisc = min(1.0_rk,self%cNDPisc / rNDFoodPisc * self%fDAssPisc)
!  N_assimilation_of_Pisc
   tNAssPisc = afNAssPisc * tNConsPisc
!  egestion_of_Pisc
   tNEgesPisc = tNConsPisc - tNAssPisc
!  respiration_of_Pisc
   tNExcrPisc = self%cNDPisc * tDRespPisc
!  mortality_of_Pisc
   tNMortPisc = self%cNDPisc * tDMortPisc
!  net_migration_flux
   tNMigrPisc = self%kMigrPisc *(self%cNDPisc * self%cDPiscIn - aNPisc)
!-----------------------------------------------------------------------
!  piscivorous fish egestion
!-----------------------------------------------------------------------
!  egestion_of_Pisc
   tDEgesPisc = tDConsPisc - tDAssPisc
!  egestion_of_Pisc
   tNEgesPisc = tNConsPisc - tNAssPisc
!  egestion_of_Pisc
   tPEgesPisc = tPConsPisc - tPAssPisc
!-----------------------------------------------------------------------
!  total flux of Fish change to state variables
!-----------------------------------------------------------------------
!  total_fish_flux_of_DW_in_Young_fish
   tDFiJv = tDMigrFiJv + tDReprFish + tDAssFiJv - tDRespFiJv - tDMortFiJv - tDConsFiJvPisc - tDAgeFish 
!  total_fish_flux_of_P_in_Young_fish
   tPFiJv = tPMigrFiJv + tPReprFish  + tPAssFiJv - tPExcrFiJv - tPMortFiJv - tPConsFiJvPisc - tPAgeFish 
!  total_fish_flux_of_N_in_Young_fish
   tNFiJv = tNMigrFiJv + tNReprFish + tNAssFiJv - tNExcrFiJv - tNMortFiJv - tNConsFiJvPisc - tNAgeFish 
!  total_fish_flux_of_DW_in_Adult_fish
   tDFiAd = tDMigrFiAd + tDAssFiAd - tDRespFiAd - tDMortFiAd - tDReprFish - tDConsFiAdPisc + tDAgeFish
!  total_fish_flux_of_P_in_Adult_fish
   tPFiAd = tPMigrFiAd + tPAssFiAd - tPExcrFiAd - tPMortFiAd - tPReprFish - tPConsFiAdPisc + tPAgeFish 
!  total_fish_flux_of_N_in_Adult_fish
   tNFiAd = tNMigrFiAd + tNAssFiAd - tNExcrFiAd - tNMortFiAd - tNReprFish - tNConsFiAdPisc + tNAgeFish 
!  total_fish_flux_of_DW_in_predatory_fish
   tDPisc = tDMigrPisc + tDAssPisc - tDRespPisc - tDMortPisc

!=======================================================================
!  fish processes relating to other modules
!=======================================================================
!-----------------------------------------------------------------------
!  Update NH4 in water
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  external state variables change due to adult fish (egestion, from sediment top)
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!  for fish it has t-, units in /m^2
!-----------------------------------------------------------------------
!  NH4_egestion_of_young_fish
   tNEgesFiJvNH4 = self%fDissEgesFish * tNEgesFiJv
!  NH4_egestion_of_adult_fish
   tNEgesFiAdNH4 = self%fDissEgesFish * tNEgesFiAd
!  total fish mortality, N
   tNMortFish = tNMortFiJv + tNMortFiAd
!  part_of_died_fish_N_fixed_in_bones_AND_scales
   tNMortFishBot = self%fDBone * tNMortFish
!  part_of_died_fish_N_becoming_dissolved_N
   tNMortFishNH4 = self%fDissMortFish *(tNMortFish - tNMortFishBot)
!  SRN_egestion_of_Pisc
   tNEgesPiscNH4 = self%fDissEgesPisc * tNEgesPisc
!  part_of_died_Pisc_N_fixed_in_bones_AND_scales
   tNMortPiscBot = self%fDBone * tNMortPisc
!  part_of_died_fish_N_becoming_dissolved_N
   tNMortPiscNH4 = self%fDissMortPisc *(tNMortPisc - tNMortPiscBot)
!  total_fish_flux_of_N_in_ammonium_in_water_in_lake_water
   tNFishNH4W = tNExcrFiJv + tNExcrFiAd + tNEgesFiJvNH4 + tNEgesFiAdNH4 &
   & + tNMortFishNH4 + tNExcrPisc + tNEgesPiscNH4 + tNMortPiscNH4

!-----------------------------------------------------------------------
!  Update PO4 in water
!  for fish it has t-, units in /m^2
!-----------------------------------------------------------------------
!  SRP_egestion_of_adult_fish
   tPEgesFiAdPO4 = self%fDissEgesFish * tPEgesFiAd
!  SRP_egestion_of_young_fish
   tPEgesFiJvPO4 = self%fDissEgesFish * tPEgesFiJv
!  total fish mortality
   tPMortFish = tPMortFiJv + tPMortFiAd
!  part_of_died_fish_P_fixed_in_bones_AND_scales
   tPMortFishBot = self%fPBone * tPMortFish
!  part_of_died_fish_P_becoming_dissolved_P
   tPMortFishPO4 = self%fDissMortFish *(tPMortFish - tPMortFishBot)
!  SRP_egestion_of_Pisc
   tPEgesPiscPO4 = self%fDissEgesPisc * tPEgesPisc
!  part_of_died_Pisc_P_fixed_in_bones_AND_scales
   tPMortPiscBot = self%fPBone * tPMortPisc
!  part_of_died_fish_P_becoming_dissolved_P
   tPMortPiscPO4 = self%fDissMortPisc *(tPMortPisc - tPMortPiscBot)
!  total_fish_flux_of_P_in_SRP_in_water_in_lake_water
   tPFishPO4W= tPExcrFiJv + tPExcrFiAd + tPEgesFiJvPO4 + tPEgesFiAdPO4 + &
   & tPMortFishPO4 + tPExcrPisc + tPEgesPiscPO4 + tPMortPiscPO4

!-----------------------------------------------------------------------
!  Update organic matter DW in water
!  for fish it has t-, units in /m^2
!-----------------------------------------------------------------------
!  bent._fish_mortality
   tDMortFish = tDMortFiJv + tDMortFiAd
!  part_of_died_fish_DW_fixed_in_bones_and_scales
   tDMortFishBot = self%fDBone * tDMortFish
!  part_of_died_fish_DW_becoming_organic matters
   tDMortFishTOM = tDMortFish - tDMortFishBot
!  part_of_died_fish_DW_fixed_in_bones_AND_scales
   tDMortPiscBot = self%fDBone * tDMortPisc
!  part_of_died_Pisc_DW_becoming_organic matters
   tDMortPiscTOM = tDMortPisc - tDMortPiscBot
!  total_fish_flux_of_DW_in_organic matter_in_lake_water
   tDFishTOMW = tDEgesFiJv + tDEgesFiAd + tDMortFishTOM + tDEgesPisc + tDMortPiscTOM
   tDFishPOMW = tDFishTOMW * (1.0_rk - self%fFisDOMW)
   tDFishDOMW = tDFishTOMW * self%fFisDOMW
!-----------------------------------------------------------------------
!  Update organic N in water
!-----------------------------------------------------------------------
!  part_of_dead_Pisc_N_becoming_orgainc_N
   tNMortPiscTOM = tNMortPisc - tNMortPiscBot - tNMortPiscNH4
!  organic_N_egestion_of_Pisc
   tNEgesPiscTOM = tNEgesPisc - tNEgesPiscNH4
!  part_of_dead_fish_NW_becoming_organic matter
   tNMortFishPOM = tNMortFish - tNMortFishBot - tNMortFishNH4
!  organic_N_egestion_of_adult_fish
   tNEgesFiAdTOM = tNEgesFiAd - tNEgesFiAdNH4
!  organic_N_egestion_of_young_fish
   tNEgesFiJvTOM = tNEgesFiJv - tNEgesFiJvNH4
!  total_fish_flux_of_N_in_organics_in_lake_water
   tNFishTOMW = tNEgesFiJvTOM + tNEgesFiAdTOM + tNMortFishPOM + tNEgesPiscTOM + tNMortPiscTOM
   tNFishPOMW = tNFishTOMW * (1.0_rk - self%fFisDOMW)
   tNFishDOMW = tNFishTOMW * self%fFisDOMW
!-----------------------------------------------------------------------
!  Update organic P in water
!-----------------------------------------------------------------------
!  part_of_died_Pisc_P_becoming_organic_P
   tPMortPiscTOM = tPMortPisc - tPMortPiscBot - tPMortPiscPO4
!  organic_P_egestion_of_Pisc
   tPEgesPiscTOM = tPEgesPisc - tPEgesPiscPO4
!  part_of_died_fish_PW_becoming_organic matter
   tPMortFishTOM = tPMortFish - tPMortFishBot - tPMortFishPO4
!  organic_P_egestion_of_young_fish
   tPEgesFiJvTOM = tPEgesFiJv - tPEgesFiJvPO4
!  organic_P_egestion_of_adult_fish
   tPEgesFiAdTOM = tPEgesFiAd - tPEgesFiAdPO4
!  total_fish_flux_of_P_in_organic matter_in_lake_water
   tPFishTOMW = tPEgesFiJvTOM + tPEgesFiAdTOM + tPMortFishTOM + tPEgesPiscTOM + tPMortPiscTOM
   tPFishPOMW = tPFishTOMW * (1.0_rk - self%fFisDOMW)
   tPFishDOMW = tPFishTOMW * self%fFisDOMW
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
   _SET_ODE_BEN_(self%id_sDFiJv,tDFiJv)
   _SET_ODE_BEN_(self%id_sPFiJv,tPFiJv)
   _SET_ODE_BEN_(self%id_sNFiJv,tNFiJv)
   _SET_ODE_BEN_(self%id_sDFiAd,tDFiAd)
   _SET_ODE_BEN_(self%id_sPFiAd,tPFiAd)
   _SET_ODE_BEN_(self%id_sNFiAd,tNFiAd)
   _SET_ODE_BEN_(self%id_sDPisc,tDPisc)   
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
!  update abiotic variables in water
   _SET_BOTTOM_EXCHANGE_(self%id_NH4poolW,  tNFishNH4W)
   _SET_BOTTOM_EXCHANGE_(self%id_PO4poolW,  tPFishPO4W)
   _SET_BOTTOM_EXCHANGE_(self%id_DPOMpoolW, tDFishPOMW)
   _SET_BOTTOM_EXCHANGE_(self%id_NPOMpoolW, tNFishPOMW)
   _SET_BOTTOM_EXCHANGE_(self%id_PPOMpoolW, tPFishPOMW)
   _SET_BOTTOM_EXCHANGE_(self%id_DFoodZoo,  -tDConsFiJv)
   _SET_BOTTOM_EXCHANGE_(self%id_NFoodZoo,  -tNConsFiJv)
   _SET_BOTTOM_EXCHANGE_(self%id_PFoodZoo,  -tPConsFiJv)
   _SET_BOTTOM_EXCHANGE_(self%id_DDOMpoolW, tDFishDOMW)
   _SET_BOTTOM_EXCHANGE_(self%id_NDOMpoolW, tNFishDOMW)
   _SET_BOTTOM_EXCHANGE_(self%id_PDOMpoolW, tPFishDOMW)
   _SET_ODE_BEN_(self%id_DFoodBen,  -tDConsFiAd)
   _SET_ODE_BEN_(self%id_NFoodBen,  -tNConsFiAd)
   _SET_ODE_BEN_(self%id_PFoodBen,  -tPConsFiAd)
!-----------------------------------------------------------------------
!  output diagnostic variables for external links
!-----------------------------------------------------------------------
!!  Export diagnostic variables
!   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aNPisc,aNPisc)
!   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aPPisc,aPPisc)
#ifdef _DEVELOPMENT_
!  output diagnostic variables for modular fluxes
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDFiJv,    tDFiJv*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPFiJv,    tPFiJv*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNFiJv,    tNFiJv*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDPisc,    tDPisc*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDFishZoo, -tDConsFiJv*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNFishZoo, -tNConsFiJv*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPFishZoo, -tPConsFiJv*secs_pr_day)
!  feh: This temperal solution
!  feh: all the variables are connected to tMortFiAd, which has dependency on
!  external variables, so can not be updated at the first time step
   if (n .GE.1000) then
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDFiAd,     tDFiAd*86400.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPFiAd,     tPFiAd*86400.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNFiAd,     tNFiAd*86400.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNFishNH4W, tNFishNH4W*86400.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPFishPO4W, tPFishPO4W*86400.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDFishPOMW, tDFishPOMW*86400.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNFishPOMW, tNFishPOMW*86400.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPFishPOMW, tPFishPOMW*86400.0_rk)
   else
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDFiAd,     0.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPFiAd,     0.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNFiAd,     0.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNFishNH4W, 0.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPFishPO4W, 0.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDFishPOMW, 0.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNFishPOMW, 0.0_rk)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPFishPOMW, 0.0_rk)
   endif
   n=n+1
#endif
! Horizontal loop end
!   _LOOP_END_
   _FABM_HORIZONTAL_LOOP_END_
!
!EOP
!-----------------------------------------------------------------------

   end subroutine do_bottom

!EOC
!-----------------------------------------------------------------------
   end module pclake_fish

!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
