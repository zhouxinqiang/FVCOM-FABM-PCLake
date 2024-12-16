#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module pclake_zooplankton
! !USES:
   use fabm_types
   use pclake_utility, ONLY:uFunTmBio

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_zooplankton
!  local state variable identifiers
!  id_sDZoo,zooplankton concentration in dry-weight, gDW/m**3
!  id_sPZoo,zooplankton concentration in nitrogen element, gN/m**3
!  id_sNZoo,zooplankton concentration in phosphorus element, gP/m**3
   type (type_state_variable_id)            :: id_sDZoo,id_sPZoo,id_sNZoo
#ifdef _DEVELOPMENT_
!  diagnostic variables for dependencies(without output)
!  diagnostic variables for modular fluxes
   type (type_diagnostic_variable_id)       :: id_wDZoo,id_wNZoo,id_wPZoo
   type (type_diagnostic_variable_id)       :: id_wNZooNO3W,id_wPZooPO4W,id_wDZooPOMW
   type (type_diagnostic_variable_id)       :: id_wNZooPOMW,id_wPZooPOMW,id_wSiZooPOMW
   type (type_diagnostic_variable_id)       :: id_wDZooDiatW,id_wNZooDiatW,id_wPZooDiatW
   type (type_diagnostic_variable_id)       :: id_wDZooGrenW,id_wNZooGrenW,id_wPZooGrenW
   type (type_diagnostic_variable_id)       :: id_wDZooBlueW,id_wNZooBlueW,id_wPZooBlueW
#endif
!  state dependencies identifiers
   type (type_state_variable_id)            :: id_DfoodDiat,id_DfoodGren,id_DfoodBlue,id_DPOMpoolW
   type (type_state_variable_id)            :: id_NfoodDiat,id_NfoodGren,id_NfoodBlue,id_NPOMpoolW
   type (type_state_variable_id)            :: id_PfoodDiat,id_PfoodGren,id_PfoodBlue,id_PPOMpoolW,id_SiPOMpoolW
   type (type_state_variable_id)            :: id_DDOMpoolW, id_NDOMpoolW,id_PDOMpoolW,id_SiDOMpoolW
   type (type_state_variable_id)            :: id_NH4poolW,id_NO3poolW,id_PO4poolW
!  environmental dependencies
   type (type_dependency_id)                :: id_uTm ,id_dz
   type (type_horizontal_dependency_id)     :: id_sDepthW
   type (type_global_dependency_id)         :: id_Day
!  Model parameters
!  parameter for temperature
   real(rk)      :: cSigTmZoo,cTmOptZoo
!  parameters for zooplankton
   real(rk)      :: cDCarrZoo,kMortZoo,kDRespZoo,cPrefDiat
   real(rk)      :: cPrefGren,cPrefBlue,cPrefPOM,hFilt
   real(rk)      :: fDAssZoo,cFiltMax
   real(rk)      :: cPDZooRef,cNDZooRef
   real(rk)      :: fDissEgesZoo,fDissMortZoo,cSiDDiat
!  nutrient ratios parameter
   real(rk)   :: cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin
   real(rk)   :: cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax
!  minimum state variable values
   real(rk)   :: cDZooMin
!  dissolved fraction of organic matter
   real(rk)   :: fZooDOMW

   contains

!  Model procedures
   procedure :: initialize
   procedure :: do

   end type type_pclake_zooplankton

!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
   real(rk)           :: Pi=3.14159265358979_rk

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
   class (type_pclake_zooplankton), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cSigTmZoo,     'cSigTmZoo',     'degree C',  'temperature constant for zooplankton (sigma)',                            default=13.0_rk)
   call self%get_parameter(self%cTmOptZoo,     'cTmOptZoo',     'degree C',  'optimum temperature for zooplankton',                                     default=25.0_rk)
   call self%get_parameter(self%cDCarrZoo,     'cDCarrZoo',     'gm-3',      'carrying capacity of zooplankton',                                        default=25.0_rk)
   call self%get_parameter(self%kMortZoo,      'kMortZoo',      'd-1',       'mortality constant for zooplankton',                                      default=0.04_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespZoo,     'kDRespZoo',     'd-1',       'maintenance respiration constant for zooplankton',                        default=0.15_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cPrefDiat,     'cPrefDiat',     '[-]',       'selection factor for diatoms',                                            default=0.75_rk)
   call self%get_parameter(self%cPrefGren,     'cPrefGren',     '[-]',       'selection factor for greens',                                             default=0.75_rk)
   call self%get_parameter(self%cPrefBlue,     'cPrefBlue',     '[-]',       'selection factor for blue-greens',                                        default=0.125_rk)
   call self%get_parameter(self%cPrefPOM,      'cPrefPOM',      '[-]',       'selection factor for particulate organic matter',                         default=0.25_rk)
   call self%get_parameter(self%hFilt,         'hFilt',         'gDW m-3',   'half-saturation constant for food conc. on zooplankton',                  default=1.0_rk)
   call self%get_parameter(self%fDAssZoo,      'fDAssZoo',      '[-]',       'dry weight assimilation efficiency of zooplankton',                       default=0.35_rk)
   call self%get_parameter(self%cFiltMax,      'cFiltMax',      'ltr/mgDW/d','maximum filtering rate',                                                  default=4.5_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cPDZooRef,     'cPDZooRef',     'mgP/mgDW',  'reference P/DW-ratio zooplankton',                                        default=0.01_rk)
   call self%get_parameter(self%cNDZooRef,     'cNDZooRef',     'mgN/mgDW',  'reference N/DW-ratio zooplankton',                                        default=0.07_rk)
   call self%get_parameter(self%fDissEgesZoo,  'fDissEgesZoo',  '[-]',       'soluble nutrient fraction of zooplankton egested food',                   default=0.25_rk)
   call self%get_parameter(self%fDissMortZoo,  'fDissMortZoo',  '[-]',       'soluble nutrient fraction of dead zooplankton',                           default=0.1_rk)
   call self%get_parameter(self%cSiDDiat,      'cSiDDiat',      'mgSi/mgDW', 'Si/DW ratio of diatoms',                                                  default=0.15_rk)
   call self%get_parameter(self%cNDDiatMin,    'cNDDiatMin',    'mgN/mgDW',  'minimum N/DW ratio diatoms',                                              default=0.01_rk)
   call self%get_parameter(self%cPDDiatMin,    'cPDDiatMin',    'mgP/mgDW',  'minimum P/DW ratio diatoms',                                              default=0.0005_rk)
   call self%get_parameter(self%cNDGrenMin,    'cNDGrenMin',    'mgN/mgDW',  'minimum N/DW ratio greens',                                               default=0.02_rk)
   call self%get_parameter(self%cPDGrenMin,    'cPDGrenMin',    'mgP/mgDW',  'minimum P/DW ratio greens',                                               default=0.0015_rk)
   call self%get_parameter(self%cNDBlueMin,    'cNDBlueMin',    'mgN/mgDW',  'minimum N/DW ratio bluegreens',                                           default=0.03_rk)
   call self%get_parameter(self%cPDBlueMin,    'cPDBlueMin',    'mgP/mgDW',  'minimum P/DW ratio bluegreens',                                           default=0.0025_rk)
   call self%get_parameter(self%cNDBlueMax,    'cNDBlueMax',    'mgN/mgDW',  'maximum N/DW ratio bluegreens',                                           default=0.15_rk)
   call self%get_parameter(self%cNDDiatMax,    'cNDDiatMax',    'mgN/mgDW',  'maximum N/DW ratio diatoms',                                              default=0.05_rk)
   call self%get_parameter(self%cNDGrenMax,    'cNDGrenMax',    'mgN/mgDW',  'maximum N/DW ratio greens',                                               default=0.1_rk)
   call self%get_parameter(self%cPDBlueMax,    'cPDBlueMax',    'mgP/mgDW',  'maximum P/DW ratio blue-greens',                                          default=0.025_rk)
   call self%get_parameter(self%cPDDiatMax,    'cPDDiatMax',    'mgP/mgDW',  'maximum P/DW ratio diatoms',                                              default=0.005_rk)
   call self%get_parameter(self%cPDGrenMax,    'cPDGrenMax',    'mgP/mgDW',  'maximum P/DW ratio greens',                                               default=0.015_rk)
!  the user defined minumun value for state variables
   call self%get_parameter(self%cDZooMin,      'cDZooMin',      'gDW/m3',    'minimum zooplankton biomass in system',                                   default=0.00001_rk)
   call self%get_parameter(self%fZooDOMW,   'fZooDOMW',   '[-]',       'dissolved organic fraction from zooplankton',                                   default=0.5_rk)
   
!  Register local state variable
!  zooplankton
   call self%register_state_variable(self%id_sDZoo,'sDZoo','gDW m-3','zooplankton DW',     &
                                    initial_value=0.05_rk,minimum=self%cDZooMin,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPZoo,'sPZoo','gP m-3','zooplankton P',     &
                                    initial_value=0.0005_rk,minimum=self%cDZooMin * self%cPDZooRef,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNZoo,'sNZoo','gN m-3','zooplankton N',     &
                                    initial_value=0.0035_rk,minimum=self%cDZooMin * self%cNDZooRef,no_river_dilution=.TRUE.)
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNZoo)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPZoo)
!  register state variables dependencies
   call self%register_state_dependency(self%id_DfoodDiat,     'diatom_as_food_DW',            'g m-3', 'diatom as food DW')
   call self%register_state_dependency(self%id_DfoodGren,     'green_as_food_DW',             'g m-3', 'green as food DW')
   call self%register_state_dependency(self%id_DfoodBlue,     'blue_as_food_DW',              'g m-3', 'blue as food DW')
   call self%register_state_dependency(self%id_NfoodDiat,     'diatom_as_food_N',             'g m-3', 'diatom as food N')
   call self%register_state_dependency(self%id_NfoodGren,     'green_as_food_N',              'g m-3', 'green as food N')
   call self%register_state_dependency(self%id_NfoodBlue,     'blue_as_food_N',               'g m-3', 'blue as food N')
   call self%register_state_dependency(self%id_PfoodDiat,     'diatom_as_food_P',             'g m-3', 'diatom as food P')
   call self%register_state_dependency(self%id_PfoodGren,     'green_as_food_P',              'g m-3', 'green as food P')
   call self%register_state_dependency(self%id_PfoodBlue,     'blue_as_food_P',               'g m-3', 'blue as food P')
   call self%register_state_dependency(self%id_DPOMpoolW,     'POM_DW_pool_water',            'g m-3', 'POM DW pool in water')
   call self%register_state_dependency(self%id_NPOMpoolW,     'POM_N_pool_water',             'g m-3', 'POM N pool in water')
   call self%register_state_dependency(self%id_PPOMpoolW,     'POM_P_pool_water',             'g m-3', 'POM P pool in water')
   call self%register_state_dependency(self%id_SiPOMpoolW,    'POM_Si_pool_water',            'g m-3', 'POM Si pool water')
   call self%register_state_dependency(self%id_NH4poolW,      'NH4_pool_water',               'g m-3', 'NH4 pool in water')
   call self%register_state_dependency(self%id_NO3poolW,      'NO3_pool_water',               'g m-3', 'NO3 pool in water')
   call self%register_state_dependency(self%id_PO4poolW,      'PO4_pool_water',               'g m-3', 'PO4 pool in water')
   call self%register_state_dependency(self%id_DDOMpoolW,     'DOM_DW_pool_water',            'g m-3', 'DOM DW in water')
   call self%register_state_dependency(self%id_NDOMpoolW,     'DOM_N_pool_water',             'g m-3', 'DOM N in water')
   call self%register_state_dependency(self%id_PDOMpoolW,     'DOM_P_pool_water',             'g m-3', 'DOM P in water')
   call self%register_state_dependency(self%id_SiDOMpoolW,    'DOM_Si_pool_water',            'g m-3', 'DOM Si in water')
#ifdef _DEVELOPMENT_
!  register diagnostic variables for modular fluxes
   call self%register_diagnostic_variable(self%id_wDZoo,      'wDZoo',                   'g m-3 s-1', 'zooplankton_DZoo_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZoo,      'wNZoo',                   'g m-3 s-1', 'zooplankton_NZoo_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZoo,      'wPZoo',                   'g m-3 s-1', 'zooplankton_PZoo_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZooNO3W,  'wNZooNO3W',               'g m-3 s-1', 'zooplankton_NO3W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZooPO4W,  'wPZooPO4W',               'g m-3 s-1', 'zooplankton_PO4W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDZooPOMW,  'wDZooPOMW',               'g m-3 s-1', 'zooplankton_DPOMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZooPOMW,  'wNZooPOMW',               'g m-3 s-1', 'zooplankton_NPOMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZooPOMW,  'wPZooPOMW',               'g m-3 s-1', 'zooplankton_PPOMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wSiZooPOMW, 'wSiZooPOMW',              'g m-3 s-1', 'zooplankton_SiPOMW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDZooDiatW, 'wDZooDiatW',              'g m-3 s-1', 'zooplankton_DDiat_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZooDiatW, 'wNZooDiatW',              'g m-3 s-1', 'zooplankton_NDiat_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZooDiatW, 'wPZooDiatW',              'g m-3 s-1', 'zooplankton_PDiat_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDZooGrenW, 'wDZooGrenW',              'g m-3 s-1', 'zooplankton_DGren_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZooGrenW, 'wNZooGrenW',              'g m-3 s-1', 'zooplankton_NGren_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZooGrenW, 'wPZooGrenW',              'g m-3 s-1', 'zooplankton_PGren_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDZooBlueW, 'wDZooBlueW',              'g m-3 s-1', 'zooplankton_DBlue_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNZooBlueW, 'wNZooBlueW',              'g m-3 s-1', 'zooplankton_NBlue_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPZooBlueW, 'wPZooBlueW',              'g m-3 s-1', 'zooplankton_PBlue_change',  output=output_instantaneous)
#endif
!  register environmental dependencies
   call self%register_dependency(self%id_uTm,    standard_variables%temperature)
   call self%register_dependency(self%id_Day,    standard_variables%number_of_days_since_start_of_the_year)
   call self%register_dependency(self%id_dz,     standard_variables%cell_thickness)
   call self%register_dependency(self%id_sDepthW,standard_variables%bottom_depth)

   return


   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP

! !IROUTINE:
!
! !INTERFACE:
   subroutine do(self,_ARGUMENTS_DO_)
!
!  INPUT PARAMETERS:
   class (type_pclake_zooplankton), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_
!  LOCAL VARIABLES:
!  Carriers for environment dependencies
   real(rk)     :: uTm,Day,dz,sDepthW
!  carriers for local state variables
   real(rk)      :: sDZoo,sNZoo,sPZoo
!  carriers for exteral link state variables
   real(rk)      :: sDDiatW,sDGrenW,sDBlueW,sDPOMW
   real(rk)      :: sNDiatW,sNGrenW,sNBlueW,sNPOMW
   real(rk)      :: sPDiatW,sPGrenW,sPBlueW,sPPOMW
!  nutrient ratios variables
   real(rk)      :: rNDDiatW,rNDGrenW,rNDBlueW,rNDPOMW
   real(rk)      :: rPDDiatW,rPDGrenW,rPDBlueW,rPDPOMW
   real(rk)      :: rPDZoo,rNDZoo
!  variables for temperature function
   real(rk)      :: uFunTmZoo
!  variables of Zooplankton DW fluxes
   real(rk)      :: aDSatZoo,ukDAssTmZoo,wDEnvZoo,ukDIncrZoo
   real(rk)      :: ukDRespTmZoo,oDFoodZoo
   real(rk)      :: oDOMW,oDPhytW
   real(rk)      :: aCorDRespZoo
!  variables for N assimilation
   real(rk)      :: afNAssZoo,wNConsPOMZoo,wDConsPOMZoo,wDConsZoo
   real(rk)      :: wNConsZoo,wNConsPhytZoo,wNConsDiatZoo, wNConsGrenZoo,wNConsBlueZoo
   real(rk)      :: wDConsDiatZoo,wDConsGrenZoo,wDConsBlueZoo
   real(rk)      :: rNDFoodZoo,oNFoodZoo
!  variables for P assimilation
   real(rk)      :: rPDFoodZoo,oPFoodZoo,afPAssZoo,wPConsZoo
   real(rk)      :: wPConsPhytZoo,wPConsPOMZoo,wPConsDiatZoo
   real(rk)      :: wPConsGrenZoo,wPConsBlueZoo
!  variables for Zooplankton dw fluxes
   real(rk)      :: wDZoo,wDAssZoo,wDRespZoo,wDMortZoo
!  variables for Zooplankton N fluxes
   real(rk)      :: wNZoo,wNAssZoo,wNExcrZoo,wNMortZoo
!  variables for Zooplankton P fluxes
   real(rk)      ::wPZoo,wPAssZoo,wPExcrZoo,wPMortZoo
!  variables for exchange of NH4
   real(rk)     :: wNZooNH4W,wNEgesZooNH4,wNEgesZoo,wNMortZooNH4
!  variables for exchange of NH3
   real(rk)     :: wNZooNO3W
!  variables for exchange of PO4
   real(rk)     :: wPZooPO4W,wPEgesZooPO4,wPEgesZoo,wPMortZooPO4
!  variables for exchange of organic DW
   real(rk)     :: wDZooPOMW,wDEgesZoo
!  variables for exchange of organic N
   real(rk)     :: wNZooPOMW,wNEgesZooTOM,wNMortZooTOM
!  variables for exchange of organic P
   real(rk)     :: wPZooPOMW,wPEgesZooTOM,wPMortZooTOM
!  variables for exchange of organic Si
   real(rk)     :: wSiZooPOMW,wSiConsDiatZoo
!  variables for exchange of dissolved organic matter
   real(rk)     :: wDZooTOMW,wNZooTOMW,wPZooTOMW
   real(rk)     :: wSiZooTOMW,wDZooDOMW,wNZooDOMW
   real(rk)     :: wPZooDOMW,wSiZooDOMW
!  variables for exchange of diatoms
   real(rk)     :: wDZooDiatW,wNZooDiatW,wPZooDiatW
!  variables for exchange of green algae
   real(rk)     :: wDZooGrenW,wNZooGrenW,wPZooGrenW
!  variables for exchange of green algae
   real(rk)     :: wDZooBlueW,wNZooBlueW,wPZooBlueW
!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
   _LOOP_BEGIN_
!  Retrieve current (local) state variable values.
   _GET_(self%id_sDZoo,sDZoo)
   _GET_(self%id_sNZoo,sNZoo)
   _GET_(self%id_sPZoo,sPZoo)
!-----------------------------------------------------------------------
!  Retrieve dependencies  
!-----------------------------------------------------------------------
!  Retrieve state dependencies 
   _GET_(self%id_DfoodDiat,sDDiatW)
   _GET_(self%id_DfoodGren,sDGrenW)
   _GET_(self%id_DfoodBlue,sDBlueW)
   _GET_(self%id_DPOMpoolW,sDPOMW)
   _GET_(self%id_NfoodDiat,sNDiatW)
   _GET_(self%id_NfoodGren,sNGrenW)
   _GET_(self%id_NfoodBlue,sNBlueW)
   _GET_(self%id_NPOMpoolW,sNPOMW)
   _GET_(self%id_PfoodDiat,sPDiatW)
   _GET_(self%id_PfoodGren,sPGrenW)
   _GET_(self%id_PfoodBlue,sPBlueW)
   _GET_(self%id_PPOMpoolW,sPPOMW)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_GLOBAL_(self%id_Day,Day)
   _GET_(self%id_dz,dz)
   _GET_HORIZONTAL_(self%id_sDepthW,sDepthW)
!-------------------------------------------------------------------------
!  Current local nutrients ratios in zooplankton (check the current state)
!-------------------------------------------------------------------------
!  P/DW_ratio of food
   rPDDiatW=sPDiatW/(sDDiatw+NearZero)
   rPDGrenW=sPGrenW/(sDGrenW+NearZero)
   rPDBlueW=sPBlueW/(sDBlueW+NearZero)
   rPDPOMW=sPPOMW/(sDPOMW+NearZero)
!  N/DW_ratio of food
   rNDDiatW=sNDiatW/(sDDiatw+NearZero)
   rNDGrenW=sNGrenW/(sDGrenW+NearZero)
   rNDBlueW=sNBlueW/(sDBlueW+NearZero)
   rNDPOMW=sNPOMW/(sDPOMW+NearZero)
   if ( rPDDiatW .GT. self%cPDDiatMax)  then
       rPDDiatW=self%cPDDiatMax
   elseif (rPDDiatW .LT. self%cPDDiatMin)  then
       rPDDiatW = self%cPDDiatMin
   else
       rPDDiatW=rPDDiatW
   endif

   if ( rPDBlueW .GT. self%cPDBlueMax)  then
       rPDBlueW=self%cPDBlueMax
   elseif (rPDBlueW .LT. self%cPDBlueMin)  then
       rPDBlueW = self%cPDBlueMin
   else
       rPDBlueW=rPDBlueW
   endif

   if ( rPDGrenW .GT. self%cPDGrenMax)  then
       rPDGrenW=self%cPDGrenMax
   elseif (rPDGrenW .LT. self%cPDGrenMin)  then
       rPDGrenW = self%cPDGrenMin
   else
       rPDGrenW=rPDGrenW
   endif
!   check for nitrogen nutrient ratios before calculating droop equations
   if ( rNDBlueW .GT. self%cNDBlueMax)  then
       rNDBlueW=self%cNDBlueMax
   elseif (rNDBlueW .LT. self%cNDBlueMin)  then
       rNDBlueW = self%cNDBlueMin
   else
       rNDBlueW =rNDBlueW
   endif

   if ( rNDDiatW .GT. self%cNDDiatMax)  then
       rNDDiatW=self%cNDDiatMax
   elseif (rNDDiatW .LT. self%cNDDiatMin)  then
       rNDDiatW = self%cNDDiatMin
   else
       rNDDiatW=rNDDiatW
   endif


   if ( rNDGrenW .GT. self%cNDGrenMax)  then
       rNDGrenW=self%cNDGrenMax
   elseif (rNDGrenW .LT. self%cNDGrenMin)  then
       rNDGrenW = self%cNDGrenMin
   else
       rNDGrenW =rNDGrenW
   endif
!  P/D_ratio_herb.zooplankton
   rPDZoo = sPZoo /(sDZoo+NearZero)
!  N/C_ratio_herb.zooplankton
   rNDZoo = sNZoo/(sDZoo+NearZero)
!-----------------------------------------------------------------------
!  temperature function
!-----------------------------------------------------------------------
!  temp._function_of_zooplankton
   uFunTmZoo =uFunTmBio(uTm,self%cSigTmZoo,self%cTmOptZoo)
!-----------------------------------------------------------------------
!  zooplankton assimilation DW
!-----------------------------------------------------------------------
!  organic_seston
   oDPhytW= sDDiatW+sDGrenW+sDBlueW
   oDOMW = sDPOMW + oDPhytW
!  food_for_zooplankton
   oDFoodZoo = self%cPrefDiat * sDDiatW + self%cPrefGren * sDGrenW + self%cPrefBlue * sDBlueW + self%cPrefPOM * sDPOMW
!  max._assimilation_rate_of_zooplankton,temp._corrected
   ukDAssTmZoo = self%fDAssZoo * self%cFiltMax * uFunTmZoo * self%hFilt
!  food_saturation_function_of_zooplankton
   aDSatZoo = oDFoodZoo /(self%hFilt + oDOMW)
!  respiration_constant_of_zooplankton
   ukDRespTmZoo = self%kDRespZoo * uFunTmZoo
!  intrinsic_rate_of_increase_of_zooplankton
   ukDIncrZoo = ukDAssTmZoo - ukDRespTmZoo - self%kMortZoo
!  environmental_correction_of_zooplankton
   wDEnvZoo = max(0.0_rk,ukDIncrZoo / self%cDCarrZoo * sDZoo*sDZoo)
!  assimilation_of_zooplankton
   wDAssZoo = aDSatZoo *(ukDAssTmZoo * sDZoo - wDEnvZoo)
!-----------------------------------------------------------------------
!  zooplankton assimilation N
!-----------------------------------------------------------------------
!  Zooplankton_food
   oNFoodZoo = self%cPrefDiat*sNDiatW + self%cPrefGren*sNGrenW + &
               & self%cPrefBlue*sNBlueW + self%cPrefPOM*sNPOMW
!  consumption_of_zooplankton
   wDConsZoo = wDAssZoo / self%fDAssZoo
!  N/C_ratio_of_zooplankton_food
   rNDFoodZoo = oNFoodZoo /(oDFoodZoo+NearZero)
!  DW_diatoms_consumption_by_zooplankton
   wDConsDiatZoo = self%cPrefDiat*sDDiatW / oDFoodZoo * wDConsZoo
!  DW_greens_consumption_by_zooplankton
   wDConsGrenZoo = self%cPrefGren*sDGrenW / oDFoodZoo * wDConsZoo
!  DW_blue-greens_consumption_by_zooplankton
   wDConsBlueZoo = self%cPrefBlue*sDBlueW / oDFoodZoo * wDConsZoo
!  N_diatom_consumption_by_zoopl.
   wNConsDiatZoo = rNDDiatW*wDConsDiatZoo
!  N_green_consumption_by_zoopl.
   wNConsGrenZoo = rNDGrenW*wDConsGrenZoo
!  N_bluegreen_consumption_by_zoopl.
   wNConsBlueZoo = rNDBlueW*wDConsBlueZoo
!  total_N_phytoplankton_consumption_by_zoopl.
   wNConsPhytZoo = wNConsDiatZoo + wNConsGrenZoo + wNConsBlueZoo
!  DW_organic matter_consumption_by_zooplankton
   wDConsPOMZoo = self%cPrefPOM*sDPOMW / oDFoodZoo * wDConsZoo
!  consumption_of_organic_N
   wNConsPOMZoo = rNDPOMW*wDConsPOMZoo
!  total_N_consumption
   wNConsZoo = wNConsPhytZoo + wNConsPOMZoo
!  N_assimilation_efficiency_of_herbivores
   afNAssZoo = min(1.0_rk,self%cNDZooRef / rNDFoodZoo * self%fDAssZoo)
!  assimilation_by_herbivores
   wNAssZoo = afNAssZoo*wNConsZoo
!-----------------------------------------------------------------------
!  zooplankton assimilation P
!-----------------------------------------------------------------------
!  Zooplankton_food
   oPFoodZoo = self%cPrefDiat*sPDiatW + self%cPrefGren*sPGrenW &
               & + self%cPrefBlue*sPBlueW + self%cPrefPOM*sPPOMW
!  P/D_ratio_of_zooplankton_food
   rPDFoodZoo = oPFoodZoo /(oDFoodZoo+NearZero)
!  P_diatom_consumption_by_zoopl.
   wPConsDiatZoo = rPDDiatW * wDConsDiatZoo
!  P_green_consumption_by_zoopl.
   wPConsGrenZoo = rPDGrenW * wDConsGrenZoo
!  P_bluegreen_consumption_by_zoopl.
   wPConsBlueZoo = rPDBlueW * wDConsBlueZoo
!  total_P_phytoplankton_consumption_by_zoopl.
   wPConsPhytZoo = wPConsDiatZoo + wPConsGrenZoo + wPConsBlueZoo
!  consumption_of_organic_P
   wPConsPOMZoo = rPDPOMW * wDConsPOMZoo
!  total_P_consumption
   wPConsZoo = wPConsPhytZoo + wPConsPOMZoo
!  P_assimilation_efficiency_of_herbivores
   afPAssZoo = min(1.0_rk,self%cPDZooRef / rPDFoodZoo * self%fDAssZoo)
!  assimilation_by_herbivores
   wPAssZoo = afPAssZoo * wPConsZoo
!-----------------------------------------------------------------------
!  zooplankton respiration and excretion
!-----------------------------------------------------------------------
!  corr._factor_of_zoopl._respiration_for_P_and_N_content
   aCorDRespZoo = max(self%cPDZooRef / rPDZoo,self%cNDZooRef / rNDZoo)
!  zoopl._respiration_DW
   wDRespZoo = aCorDRespZoo * self%kDRespZoo * uFunTmZoo * sDZoo
!  N_excretion
   wNExcrZoo = rNDZoo / self%cNDZooRef * self%kDRespZoo * uFunTmZoo*sNZoo
!  P_excretion
   wPExcrZoo = rPDZoo / self%cPDZooRef * self%kDRespZoo * uFunTmZoo*sPZoo
!-----------------------------------------------------------------------
!  zooplankton mortality
!-----------------------------------------------------------------------
!  zoopl._mortality,incl._environmental_correction_DW
   wDMortZoo = self%kMortZoo * sDZoo +(1.0_rk - aDSatZoo) * wDEnvZoo
!  zoopl._mortality_N
   wNMortZoo = rNDZoo*wDMortZoo
!  zoopl._mortality_P
   wPMortZoo = rPDZoo * wDMortZoo
!-----------------------------------------------------------------------
!  zooplankton egestion
!-----------------------------------------------------------------------
!  egestion_of_zooplankton
   wDEgesZoo = wDConsZoo - wDAssZoo
!  N_egestion
   wNEgesZoo = wNConsZoo - wNAssZoo
!  P_egestion
    wPEgesZoo = wPConsZoo - wPAssZoo
!-----------------------------------------------------------------------
!  total flux of state variables
!-----------------------------------------------------------------------
!  total_flux_of_DW_in_Herbivorous_zooplankton
   wDZoo = wDAssZoo - wDRespZoo - wDMortZoo
!  total_flux_of_N_in_Herbivorous_zooplankton
   wNZoo = wNAssZoo - wNExcrZoo - wNMortZoo
!  total_flux_of_P_in_Herbivorous_zooplankton
   wPZoo = wPAssZoo - wPExcrZoo - wPMortZoo
!=======================================================================
!  zooplankton processes relating to other modules
!=======================================================================
!-----------------------------------------------------------------------
!  Update NH4 in water
!-----------------------------------------------------------------------
!  soluble_N_egestion
   wNEgesZooNH4 = self%fDissEgesZoo*wNEgesZoo
!  soluble_N_mortality
   wNMortZooNH4 = self%fDissMortZoo*wNMortZoo
!  total_Zoo_flux_of_N_in_ammonium_in_water_in_lake_water
   wNZooNH4W = wNExcrZoo + wNEgesZooNH4 + wNMortZooNH4
!-----------------------------------------------------------------------
!  Update NO3 in water (zooplankton is currently assumed to contribute NH4 only, not NO3)
!-----------------------------------------------------------------------
!  total_Zoo_flux_of_N_in_nitrate_in_water_in_lake_water
   wNZooNO3W = 0.0_rk
!-----------------------------------------------------------------------
!  Update PO4 in water
!-----------------------------------------------------------------------
!  soluble_P_egestion
   wPEgesZooPO4 = self%fDissEgesZoo*wPEgesZoo
!  soluble_P_mortality
   wPMortZooPO4 = self%fDissMortZoo * wPMortZoo
!  total_Zoo_flux_of_P_in_SRP_in_water_in_lake_water
   wPZooPO4W = wPExcrZoo + wPEgesZooPO4 + wPMortZooPO4
!-----------------------------------------------------------------------
!  Update organic DW in water
!-----------------------------------------------------------------------
!  total_Zoo_flux_of_DW_in_organic matter_in_lake_wate
   wDZooTOMW = - wDConsPOMZoo + wDEgesZoo + wDMortZoo
   wDZooPOMW = wDZooTOMW * (1.0_rk - self%fZooDOMW)
   wDZooDOMW = wDZooTOMW * self%fZooDOMW
!-----------------------------------------------------------------------
!  Update organic N in water
!-----------------------------------------------------------------------
!  organic_N_mortality
   wNMortZooTOM = wNMortZoo - wNMortZooNH4
!  organic_N_egestion
   wNEgesZooTOM = wNEgesZoo - wNEgesZooNH4
!  total_Zoo_flux_of_N_in_organic matter_in_lake_water
   wNZooTOMW = - wNConsPOMZoo + wNEgesZooTOM + wNMortZooTOM
   wNZooPOMW = wNZooTOMW * (1.0_rk - self%fZooDOMW)
   wNZooDOMW = wNZooTOMW * self%fZooDOMW
!-----------------------------------------------------------------------
!  Update organic P in water
!-----------------------------------------------------------------------
!  organic_P_mortality
   wPMortZooTOM = wPMortZoo - wPMortZooPO4
!  organic_P_egestion
   wPEgesZooTOM = wPEgesZoo - wPEgesZooPO4
!  total_Zoo_flux_of_P_in_organic matter_in_lake_water
   wPZooTOMW = - wPConsPOMZoo + wPEgesZooTOM + wPMortZooTOM
   wPZooPOMW = wPZooTOMW * (1.0_rk - self%fZooDOMW)
   wPZooDOMW = wPZooTOMW * self%fZooDOMW
!-----------------------------------------------------------------------
!  Update organic Si in water
!-----------------------------------------------------------------------
!  consumption_of_diatoms
   wSiConsDiatZoo = self%cSiDDiat * wDConsDiatZoo
!  total_Zoo_flux_of_silica_in_lake_water_organic matter
   wSiZooPOMW = wSiConsDiatZoo * (1.0_rk - self%fZooDOMW)
   wSiZooDOMW = wSiConsDiatZoo * self%fZooDOMW   
!-----------------------------------------------------------------------
!  Update diatom state variables
!-----------------------------------------------------------------------
!  total_Zoo_flux_of_DW_in_Diatoms_in_lake_water
   wDZooDiatW = - wDConsDiatZoo
!  total_Zoo_flux_of_N_in_Diatoms_in_lake_water
   wNZooDiatW = - wNConsDiatZoo
!  total_Zoo_flux_of_P_in_Diatoms_in_lake_water
   wPZooDiatW = - wPConsDiatZoo
!-----------------------------------------------------------------------
!  Update green algae state variables
!-----------------------------------------------------------------------
!  total_Zoo_flux_of_DW_in_Greens_in_lake_water
   wDZooGrenW = - wDConsGrenZoo
!  total_Zoo_flux_of_N_in_Greens_in_lake_water
   wNZooGrenW = - wNConsGrenZoo
!  total_Zoo_flux_of_P_in_Greens_in_lake_water
   wPZooGrenW = - wPConsGrenZoo
!-----------------------------------------------------------------------
!  Update blue-green algae state variables
!-----------------------------------------------------------------------
!  total_Zoo_flux_of_DW_in_Blue-greens_in_lake_water
   wDZooBlueW = - wDConsBlueZoo
!  total_Zoo_flux_of_N_in_Blue-greens_in_lake_water
   wNZooBlueW = - wNConsBlueZoo
!  total_Zoo_flux_of_P_in_Blue-greens_in_lake_water
   wPZooBlueW = - wPConsBlueZoo
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
   _SET_ODE_(self%id_sDZoo,wDZoo)
   _SET_ODE_(self%id_sNZoo,wNZoo)
   _SET_ODE_(self%id_sPZoo,wPZoo)
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
!  update abiotic variables in water
   _SET_ODE_(self%id_NH4poolW,      wNZooNH4W)
   _SET_ODE_(self%id_NO3poolW,      wNZooNO3W)
   _SET_ODE_(self%id_PO4poolW,      wPZooPO4W)
   _SET_ODE_(self%id_DPOMpoolW,     wDZooPOMW)
   _SET_ODE_(self%id_NPOMpoolW,     wNZooPOMW)
   _SET_ODE_(self%id_PPOMpoolW,     wPZooPOMW)
   _SET_ODE_(self%id_SiPOMpoolW,    wSiZooPOMW)
   _SET_ODE_(self%id_DDOMpoolW,     wDZooDOMW)
   _SET_ODE_(self%id_NDOMpoolW,     wNZooDOMW)
   _SET_ODE_(self%id_PDOMpoolW,     wPZooDOMW)
   _SET_ODE_(self%id_SiDOMpoolW,    wSiZooDOMW)
!  update phytoplankton in water
   _SET_ODE_(self%id_DfoodDiat,     wDZooDiatW)
   _SET_ODE_(self%id_NfoodDiat,     wNZooDiatW)
   _SET_ODE_(self%id_PfoodDiat,     wPZooDiatW)
   _SET_ODE_(self%id_DfoodGren,     wDZooGrenW)
   _SET_ODE_(self%id_NfoodGren,     wNZooGrenW)
   _SET_ODE_(self%id_PfoodGren,     wPZooGrenW)
   _SET_ODE_(self%id_DfoodBlue,     wDZooBlueW)
   _SET_ODE_(self%id_NfoodBlue,     wNZooBlueW)
   _SET_ODE_(self%id_PfoodBlue,     wPZooBlueW)
#ifdef _DEVELOPMENT_
!  output diagnostic variables for modular fluxes
   _SET_DIAGNOSTIC_(self%id_wDZoo,      wDZoo*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wNZoo,      wNZoo*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wPZoo,      wPZoo*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wNZooNO3W,  wNZooNO3W*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wPZooPO4W,  wPZooPO4W*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wDZooPOMW,  wDZooPOMW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wNZooPOMW,  wNZooPOMW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wPZooPOMW,  wPZooPOMW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wSiZooPOMW, wSiZooPOMW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wDZooDiatW, wDZooDiatW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wNZooDiatW, wNZooDiatW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wPZooDiatW, wPZooDiatW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wDZooGrenW, wDZooGrenW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wNZooGrenW, wNZooGrenW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wPZooGrenW, wPZooGrenW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wDZooBlueW, wDZooBlueW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wNZooBlueW, wNZooBlueW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wPZooBlueW, wPZooBlueW*secs_pr_day)
#endif
! Spatial loop end
   _LOOP_END_

!
!EOP
!-----------------------------------------------------------------------

   end subroutine do
!EOC
!-----------------------------------------------------------------------
!
   end module pclake_zooplankton

!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
