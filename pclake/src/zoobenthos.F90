#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module pclake_zoobenthos
! !USES:
   use fabm_types
   use fabm_expressions
   use pclake_utility, ONLY:uFunTmBio
   implicit none
!  default: all is private.
   private
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_zoobenthos
!  local state variable identifiers
!  id_sDBent,zoobenthos concentration in dry-weight, gDW/m**2
!  id_sPBent,zoobenthos concentration in nitrogen element, gN/m**2
!  id_sNBent,zoobenthos concentration in phosphorus element, gP/m**2
   type (type_bottom_state_variable_id) :: id_sDBent,id_sPBent,id_sNBent
!  diagnostic variables for dependencies(without output)
   type (type_horizontal_diagnostic_variable_id)       :: id_tDBenPOMS
#ifdef _DEVELOPMENT_
!  diagnostic variable for modular fluxes
   type (type_horizontal_diagnostic_variable_id)       :: id_tDBenBent,id_tPBenBent,id_tNBenBent
   type (type_horizontal_diagnostic_variable_id)       :: id_tNBenNH4S,id_tNBenNO3S,id_tPBenPO4S
   type (type_horizontal_diagnostic_variable_id)       :: id_tDBenPOMSflux,id_tNBenPOMS,id_tPBenPOMS
   type (type_horizontal_diagnostic_variable_id)       :: id_tSiBenPOMS,id_tDBenDiatS,id_tNBenDiatS
   type (type_horizontal_diagnostic_variable_id)       :: id_tPBenDiatS,id_tDBenGrenS,id_tNBenGrenS
   type (type_horizontal_diagnostic_variable_id)       :: id_tPBenGrenS,id_tDBenBlueS,id_tNBenBlueS
   type (type_horizontal_diagnostic_variable_id)       :: id_tPBenBlueS,id_tDAssFiAd,id_tNAssFiAd
   type (type_horizontal_diagnostic_variable_id)       :: id_tPAssFiAd,id_tNBenNH4W,id_tPBenPO4W
   type (type_horizontal_diagnostic_variable_id)       :: id_tDBenPOMW,id_tNBenPOMW,id_tPBenPOMW
#endif
!  state dependencies identifiers
   type (type_bottom_state_variable_id)            :: id_DfoodDiatS,id_DfoodGrenS,id_DfoodBlueS,id_DPOMpoolS
   type (type_bottom_state_variable_id)            :: id_NfoodDiatS,id_NfoodGrenS,id_NfoodBlueS,id_NPOMpoolS
   type (type_bottom_state_variable_id)            :: id_PfoodDiatS,id_PfoodGrenS,id_PfoodBlueS,id_PPOMpoolS
   type (type_bottom_state_variable_id)            :: id_NH4poolS,id_NO3poolS,id_PO4poolS,id_SiPOMpoolS
   type (type_bottom_state_variable_id)            :: id_DDOMpoolS,id_NDOMpoolS, id_PDOMpoolS,id_SiDOMpoolS
!  environmental dependencies
   type (type_dependency_id)                       :: id_uTm,id_dz
!  Model parameters
   real(rk)           :: cDBentIn,kMigrBent,cDCarrBent,kDAssBent,hDFoodBent,fDAssBent,fDissEgesBent
   real(rk)           :: kDRespBent,kMortBent,fDissMortBent,cTmOptBent,cSigTmBent,cPDBentRef
   real(rk)           :: cNDBentRef,cSiDDiat
!  nutrient ratios parameter
   real(rk)   :: cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin
   real(rk)   :: cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax
!  minimum state variable values
   real(rk)   :: cDBentMin
!  fraction of dissolved organic matter from zoobenthos
   real(rk)   :: fBenDOMS
   contains
!  Module procedures
   procedure :: initialize
   procedure :: do_bottom
   end type type_pclake_zoobenthos
!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero=0.000000000000000000000000000000001_rk
!
!EOP
!-----------------------------------------------------------------------
   contains
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the passive tracer model
!
! !INTERFACE:
   subroutine initialize(self,configunit)
! !INPUT PARAMETERS:
   class (type_pclake_zoobenthos), intent(inout),      target :: self
   integer,                          intent(in)            :: configunit
!EOP
!-----------------------------------------------------------------------
!BOC
!  Store parameter values in derived type
!  NB: all rates must be provided in values per day in .yaml input,
!  and are converted here to values per second.
   call self%get_parameter(self%cDBentIn,     'cDBentIn',     'gDW m-2',  'external zoobenthos density',                                default=0.01_rk)
   call self%get_parameter(self%kMigrBent,    'kMigrBent',    'd-1',      'zoobenthos migration rate',                                  default=0.001_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cDCarrBent,   'cDCarrBent',   'gDW m-2',  'carrying capacity of zoobenthos',                            default=10.0_rk)
   call self%get_parameter(self%kDAssBent,    'kDAssBent',    'd-1',      'maximum assimilation rate',                                  default=0.1_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%hDFoodBent,   'hDFoodBent',   'g m-2',    'half saturation food for zoobenthos',                        default=200.0_rk)
   call self%get_parameter(self%fDAssBent,    'fDAssBent',    '[-]',      'C assimilation efficiency of zoobenthos',                    default=0.3_rk)
   call self%get_parameter(self%fDissEgesBent,'fDissEgesBent','[-]',      'soluble nutrient fraction of by zoobenthos egested food',    default=0.25_rk)
   call self%get_parameter(self%kDRespBent,   'kDRespBent',   'd-1',      'respiration constant of zoobenthos',                         default=0.005_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kMortBent,    'kMortBent',    'd-1',      'mortality constant of zoobenthos',                           default=0.005_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fDissMortBent,'fDissMortBent','[-]',      'soluble P fraction of dead zoobenthos P',                    default=0.1_rk)
   call self%get_parameter(self%cTmOptBent,   'cTmOptBent',   'degreee C','optimal temperature for zoobenthos',                         default=25.0_rk)
   call self%get_parameter(self%cSigTmBent,   'cSigTmBent',   'degreee C','temperature constant for zoobenthos (sigma)',                default=16.0_rk)
   call self%get_parameter(self%cPDBentRef,   'cPDBentRef',   'mgP/mgDW', 'reference P/C ratio of zoobenthos',                          default=0.01_rk)
   call self%get_parameter(self%cNDBentRef,   'cNDBentRef',   'mgN/mgDW', 'reference N/C ratio of zoobenthos',                          default=0.07_rk)
   call self%get_parameter(self%cSiDDiat,     'cSiDDiat',     'mgSi/mgDW','Si/DW ratio of diatoms',                                     default=0.15_rk)
   call self%get_parameter(self%cNDDiatMin,   'cNDDiatMin',   'mgN/mgDW', 'minimum N/DW ratio diatoms',                                default=0.01_rk)
   call self%get_parameter(self%cPDDiatMin,   'cPDDiatMin',   'mgP/mgDW', 'minimum P/DW ratio diatoms',                                default=0.0005_rk)
   call self%get_parameter(self%cNDGrenMin,   'cNDGrenMin',   'mgN/mgDW', 'minimum N/DW ratio greens',                                 default=0.02_rk)
   call self%get_parameter(self%cPDGrenMin,   'cPDGrenMin',   'mgP/mgDW', 'minimum P/DW ratio greens',                                 default=0.0015_rk)
   call self%get_parameter(self%cNDBlueMin,   'cNDBlueMin',   'mgN/mgDW', 'minimum N/DW ratio blue-greens',                             default=0.03_rk)
   call self%get_parameter(self%cPDBlueMin,   'cPDBlueMin',   'mgP/mgDW', 'minimum P/DW ratio blue-greens',                             default=0.0025_rk)
   call self%get_parameter(self%cNDBlueMax,   'cNDBlueMax',   'mgN/mgDW', 'maximum N/DW ratio blue-greens',                                default=0.15_rk)
   call self%get_parameter(self%cNDDiatMax,   'cNDDiatMax',   'mgN/mgDW', 'maximum N/DW ratio diatoms',                                   default=0.005_rk)
   call self%get_parameter(self%cNDGrenMax,   'cNDGrenMax',   'mgN/mgDW', 'maximum N/DW ratio greens',                                    default=0.1_rk)
   call self%get_parameter(self%cPDBlueMax,   'cPDBlueMax',   'mgP/mgDW', 'maximum P/DW ratio blue-greens',                               default=0.025_rk)
   call self%get_parameter(self%cPDDiatMax,   'cPDDiatMax',   'mgP/mgDW', 'maximum P/DW ratio diatoms',                                   default=0.05_rk)
   call self%get_parameter(self%cPDGrenMax,   'cPDGrenMax',   'mgP/mgDW', 'maximum P/DW ratio greens',                                    default=0.015_rk)
!  the user defined minumun value for state variables
   call self%get_parameter(self%cDBentMin,    'cDBentMin',    'gDW/m2',   'minimum zoobenthos concentration in system',                   default=0.00001_rk)
   call self%get_parameter(self%fBenDOMS,  'fBenDOMS',   '[-]',           'dissolved organic fraction from zoobenthos',                   default=0.5_rk)
!  Register local state variable
   call self%register_state_variable(self%id_sDBent,'sDBent','gDW m-2','zoobenthos DW',     &
                                    initial_value=1.0_rk,minimum=self%cDBentMin)
   call self%register_state_variable(self%id_sPBent,'sPBent','gP m-2','zoobenthos P',     &
                                    initial_value=0.1_rk,minimum=self%cDBentMin * self%cPDBentRef)
   call self%register_state_variable(self%id_sNBent,'sNBent','gN m-2','zoobenthos N',     &
                                    initial_value=0.01_rk,minimum=self%cDBentMin * self%cNDBentRef)
!  Register diagnostic variables for dependencies in other modules
   call self%register_diagnostic_variable(self%id_tDBenPOMS,     'tDBenPOMS',    'g m-2 s-1', 'tDBenPOMS',                output=output_none)
#ifdef _DEVELOPMENT_
!  Register diagnostic variables for modular fluxes
   call self%register_diagnostic_variable(self%id_tDBenBent,     'tDBenBent',     'g m-2',    'zoobenthos_DBent_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBenBent,     'tPBenBent',     'g m-2',    'zoobenthos_PBent_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBenBent,     'tNBenBent',     'g m-2',    'zoobenthos_NBent_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBenNH4S,     'tNBenNH4S',     'g m-2',    'zoobenthos_NH4S_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBenNO3S,     'tNBenNO3S',     'g m-2',    'zoobenthos_NO3S_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBenPO4S,     'tPBenPO4S',     'g m-2',    'zoobenthos_PO4S_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDBenPOMSflux, 'tDBenPOMSflux', 'g m-2',    'zoobenthos_DPOMS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBenPOMS,     'tNBenPOMS',     'g m-2',    'zoobenthos_NPOMS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBenPOMS,     'tPBenPOMS',     'g m-2',    'zoobenthos_PPOMS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tSiBenPOMS,    'tSiBenPOMS',    'g m-2',    'zoobenthos_SiPOMS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDBenDiatS,    'tDBenDiatS',    'g m-2',    'zoobenthos_DDiatS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBenDiatS,    'tNBenDiatS',    'g m-2',    'zoobenthos_NDiatS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBenDiatS,    'tPBenDiatS',    'g m-2',    'zoobenthos_PDiatS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDBenGrenS,    'tDBenGrenS',    'g m-2',    'zoobenthos_DGrenS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBenGrenS,    'tNBenGrenS',    'g m-2',    'zoobenthos_NGrenS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBenGrenS,    'tPBenGrenS',    'g m-2',    'zoobenthos_PGrenS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDBenBlueS,    'tDBenBlueS',    'g m-2',    'zoobenthos_DBlueS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBenBlueS,    'tNBenBlueS',    'g m-2',    'zoobenthos_NBlueS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBenBlueS,    'tPBenBlueS',    'g m-2',    'zoobenthos_PBlueS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAssFiAd,     'tDAssFiAd',     'g m-2',    'zoobenthos_DFiAd_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAssFiAd,     'tNAssFiAd',     'g m-2',    'zoobenthos_NFiAd_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAssFiAd,     'tPAssFiAd',     'g m-2',    'zoobenthos_PFiAd_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBenNH4W,     'tNBenNH4W',     'g m-2',    'zoobenthos_NH4W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBenPO4W,     'tPBenPO4W',     'g m-2',    'zoobenthos_PO4W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDBenPOMW,     'tDBenPOMW',     'g m-2',    'zoobenthos_DPOMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNBenPOMW,     'tNBenPOMW',     'g m-2',    'zoobenthos_NPOMw_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPBenPOMW,     'tPBenPOMW',     'g m-2',    'zoobenthos_PPOMw_change',  output=output_instantaneous)
#endif
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_sNBent)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPBent)
!  register state variables dependencies
   call self%register_state_dependency(self%id_DfoodDiatS,   'diatom_as_food_DW',                'g m-2', 'diatom as food DW')
   call self%register_state_dependency(self%id_DfoodGrenS,   'green_as_food_DW',                 'g m-2', 'green as food DW')
   call self%register_state_dependency(self%id_DfoodBlueS,   'blue_as_food_DW',                  'g m-2', 'blue-green as food DW')
   call self%register_state_dependency(self%id_NfoodDiatS,   'diatom_as_food_N',                 'g m-2', 'diatom as food N')
   call self%register_state_dependency(self%id_NfoodGrenS,   'green_as_food_N',                  'g m-2', 'green as food N')
   call self%register_state_dependency(self%id_NfoodBlueS,   'blue_as_food_N',                   'g m-2', 'blue-green as food N')
   call self%register_state_dependency(self%id_PfoodDiatS,   'diatom_as_food_P',                 'g m-2', 'diatom as food P')
   call self%register_state_dependency(self%id_PfoodGrenS,   'green_as_food_P',                  'g m-2', 'green as food P')
   call self%register_state_dependency(self%id_PfoodBlueS,   'blue_as_food_P',                   'g m-2', 'blue-green as food P')
   call self%register_state_dependency(self%id_DPOMpoolS,    'POM_DW_pool_sediment',        'g m-2', 'particulate organic matter DW pool in sediment')
   call self%register_state_dependency(self%id_PPOMpoolS,    'POM_P_pool_sediment',         'g m-2', 'particulate organic matter P pool in sediment')
   call self%register_state_dependency(self%id_NPOMpoolS,    'POM_N_pool_sediment',         'g m-2', 'particulate organic matter N pool in sediment')
   call self%register_state_dependency(self%id_SiPOMpoolS,   'POM_Si_pool_sediment',        'g m-2', 'particulate organic matter Si pool sediment')
   call self%register_state_dependency(self%id_NH4poolS,     'NH4_pool_sediment',                'g m-2', 'NH4 pool in sediment')
   call self%register_state_dependency(self%id_NO3poolS,     'NO3_pool_sediment',                'g m-2', 'NO3 pool in sediment')
   call self%register_state_dependency(self%id_PO4poolS,     'PO4_pool_sediment',                'g m-2', 'PO4 pool in sediment')
   call self%register_state_dependency(self%id_DDOMpoolS,    'DOM_DW_pool_sediment',   'g m-2', 'dissolved organic matter DW in sediment')
   call self%register_state_dependency(self%id_NDOMpoolS,    'DOM_N_pool_sediment',    'g m-2', 'dissolved organic matter N in sediment')
   call self%register_state_dependency(self%id_PDOMpoolS,    'DOM_P_pool_sediment',    'g m-2', 'dissolved organic matter P in sediment')
   call self%register_state_dependency(self%id_SiDOMpoolS,   'DOM_Si_pool_sediment',   'g m-2', 'dissolved organic matter Si in sediment')
!  register environmental dependencies
   call self%register_dependency(self%id_uTm, standard_variables%temperature)
   call self%register_dependency(self%id_dz,  standard_variables%cell_thickness)

   return

   end subroutine initialize
!EOC
!-----------------------------------------------------------------------
!BOP
!
 !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_pclake_zoobenthos), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
! !LOCAL VARIABLES:
!  state variables value carriers
   real(rk)                   :: sDBent,sPBent,sNBent
!  environmental dependencies carriers
   real(rk)                   :: uTm,dz
!  external links variable carriers
   real(rk)         :: sDDiatS,sDGrenS,sDBlueS,sDPOMS
   real(rk)         :: sNDiatS,sNGrenS,sNBlueS,sNPOMS
   real(rk)         :: sPDiatS,sPGrenS,sPBlueS,sPPOMS
!  variables for status auxilaries
   real(rk)          :: aDPhytS,aPPhytS,aNPhytS
   real(rk)          :: rPDBent,rNDBent
   real(rk)          :: rPDBlueS,rPDGrenS,rPDDiatS,rPDPOMS,rPDFoodBent
   real(rk)          :: rNDPOMS,rNDDiatS,rNDGrenS,rNDBlueS,rNDFoodBent
!  variables for temperature functions
   real(rk)          :: uFunTmBent
!  variables for DW fluxes
   real(rk)          :: tDBenBent,tDMigrBent,tDAssBent
   real(rk)          :: tDRespBent,tDMortBent,aDFoodBent
   real(rk)          :: aDSatBent,tDEnvBent,ukDIncrBent,tDConsBent
   real(rk)          :: tDConsBlueBent,tDConsGrenBent,tDConsDiatBent
!  variables for P fluxes
   real(rk)          :: tPBenBent,tPMigrBent,tPAssBent
   real(rk)          :: tPExcrBent,tPMortBent,aPFoodBent,afPAssBent
   real(rk)          :: tPConsBent,tPConsPOMBent,tPConsPhytBent
   real(rk)          :: tPConsDiatBent,tPConsGrenBent,tPConsBlueBent
   real(rk)          :: tDConsPOMBent
!  variables for N fluxes
   real(rk)          :: tNBenBent,tNMigrBent,tNAssBent
   real(rk)          :: tNExcrBent,tNMortBent,aNFoodBent,afNAssBent
   real(rk)          :: tNConsDiatBent,tNConsGrenBent,tNConsBlueBent
   real(rk)          :: tNConsPhytBent,tNConsPOMBent,tNConsBent
!  variables for exchange of NH4S
   real(rk)          :: tNBenNH4S,tNEgesBentNH4,tNEgesBent,tNMortBentNH4
!  variables for exchange of NO3S
   real(rk)          :: tNBenNO3S
!  variables for exchange of PO4S
   real(rk)          :: tPBenPO4S,tPEgesBentPO4,tPEgesBent,tPMortBentPO4
!  variables for exchange for organic matter
   real(rk)          :: tDBenPOMS,tDEgesBent,tNBenPOMS,tNEgesBentTOM
   real(rk)          :: tNMortBentTOM,tPBenPOMS,tPEgesBentTOM,tPMortBentTOM
   real(rk)          :: tSiBenPOMS,tSiConsDiatBent
   real(rk)          :: tDBenTOMS,tNBenTOMS,tPBenTOMS,tSiBenTOMS
   real(rk)          :: tDBenDOMS,tNBenDOMS,tPBenDOMS,tSiBenDOMS
!  variables for exchange for diatom
   real(rk)          :: tDBenDiatS,tNBenDiatS,tPBenDiatS
!  variables for exchange for green algae
   real(rk)          :: tDBenGrenS,tNBenGrenS,tPBenGrenS
!  variables for exchange for green algae
   real(rk)          :: tDBenBlueS,tNBenBlueS,tPBenBlueS

!  Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_
!-----------------------------------------------------------------------
!  Retrieve current (local) state variable values.
!-----------------------------------------------------------------------
   _GET_HORIZONTAL_(self%id_sDBent,sDBent)
   _GET_HORIZONTAL_(self%id_sPBent,sPBent)
   _GET_HORIZONTAL_(self%id_sNBent,sNBent)
!-----------------------------------------------------------------------
!  Retrieve dependencies  value
!-----------------------------------------------------------------------
!  Retrieve state dependencie value
   _GET_HORIZONTAL_(self%id_DPOMpoolS,sDPOMS)
   _GET_HORIZONTAL_(self%id_NPOMpoolS,sNPOMS)
   _GET_HORIZONTAL_(self%id_PPOMpoolS,sPPOMS)
   _GET_HORIZONTAL_(self%id_DfoodDiatS,sDDiatS)
   _GET_HORIZONTAL_(self%id_DfoodGrenS,sDGrenS)
   _GET_HORIZONTAL_(self%id_DfoodBlueS,sDBlueS)
   _GET_HORIZONTAL_(self%id_NfoodDiatS,sNDiatS)
   _GET_HORIZONTAL_(self%id_NfoodGrenS,sNGrenS)
   _GET_HORIZONTAL_(self%id_NfoodBlueS,sNBlueS)
   _GET_HORIZONTAL_(self%id_PfoodDiatS,sPDiatS)
   _GET_HORIZONTAL_(self%id_PfoodGrenS,sPGrenS)
   _GET_HORIZONTAL_(self%id_PfoodBlueS,sPBlueS)

!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_(self%id_dz,dz)
!----------------------------------------------------------------------
!  Current local nutrients ratios in zoobenthos (check the current state)
!----------------------------------------------------------------------
   rPDBent=sPBent/(sDBent+NearZero)
   rNDBent=sNBent/(sDBent+NearZero)
   rPDPOMS=sPPOMS/(sDPOMS+NearZero)
   rNDPOMS=sNPOMS/(sDPOMS+NearZero)
   rPDBlueS=sPBlueS/(sDBlueS+NearZero)
   rPDGrenS=sPGrenS/(sDGrenS+NearZero)
   rPDDiatS=sPDiatS/(sDDiatS+NearZero)
   rNDDiatS=sNBlueS/(sDBlueS+NearZero)
   rNDGrenS=sNGrenS/(sDGrenS+NearZero)
   rNDBlueS=sNDiatS/(sDDiatS+NearZero)
!  check for phosphorus nutrient ratios
   if ( rPDDiatS .GT. self%cPDDiatMax)  then
       rPDDiatS=self%cPDDiatMax
   elseif (rPDDiatS .LT. self%cPDDiatMin)  then
       rPDDiatS = self%cPDDiatMin
   else
       rPDDiatS=rPDDiatS
   endif

   if ( rPDBlueS .GT. self%cPDBlueMax)  then
       rPDBlueS=self%cPDBlueMax
   elseif (rPDBlueS .LT. self%cPDBlueMin)  then
       rPDBlueS = self%cPDBlueMin
   else
       rPDBlueS=rPDBlueS
   endif

   if ( rPDGrenS .GT. self%cPDGrenMax)  then
       rPDGrenS=self%cPDGrenMax
   elseif (rPDGrenS .LT. self%cPDGrenMin)  then
       rPDGrenS = self%cPDGrenMin
   else
       rPDGrenS=rPDGrenS
   endif
!  check for nitrogen nutrient ratios
   if ( rNDBlueS .GT. self%cNDBlueMax)  then
       rNDBlueS=self%cNDBlueMax
   elseif (rNDBlueS .LT. self%cNDBlueMin)  then
       rNDBlueS = self%cNDBlueMin
   else
       rNDBlueS =rNDBlueS
   endif

   if ( rNDDiatS .GT. self%cNDDiatMax)  then
       rNDDiatS=self%cNDDiatMax
   elseif (rNDDiatS .LT. self%cNDDiatMin)  then
       rNDDiatS = self%cNDDiatMin
   else
       rNDDiatS=rNDDiatS
   endif


   if ( rNDGrenS .GT. self%cNDGrenMax)  then
       rNDGrenS=self%cNDGrenMax
   elseif (rNDGrenS .LT. self%cNDGrenMin)  then
       rNDGrenS = self%cNDGrenMin
   else
       rNDGrenS =rNDGrenS
   endif

!  auxilaries for phytoplankton
   aDPhytS=sDDiatS+sDGrenS+sDBlueS
   aPPhytS=sPDiatS+sPGrenS+sPBlueS
   aNPhytS=sNDiatS+sNGrenS+sNBlueS
!-----------------------------------------------------------------------
!  temperature function
!-----------------------------------------------------------------------
!  temp._function_for_zoobenthos
   uFunTmBent = uFunTmBio(uTm,self%cSigTmBent,self%cTmOptBent)
!-----------------------------------------------------------------------
!  zoobenthos migration
!-----------------------------------------------------------------------
!  migration_flux
   tDMigrBent = self%kMigrBent *(self%cDBentIn - sDBent)
!  net_migration_flux
   tPMigrBent = self%kMigrBent *(self%cPDBentRef*self%cDBentIn - sPBent)
!  Net_migration_flux
   tNMigrBent = self%kMigrBent *(self%cNDBentRef*self%cDBentIn - sNBent)
!-----------------------------------------------------------------------
!  zoobenthos assimilation,DW
!-----------------------------------------------------------------------
!  food_for_zoobenthos
   aDFoodBent = sDPOMS + aDPhytS
!  food_limitation_function_of_zoobenthos
   aDSatBent = aDFoodBent /(self%hDFoodBent + aDFoodBent)
!  intrinsic_net_increase_rate_of_zoobenthos
   ukDIncrBent = (self%kDAssBent - self%kDRespBent) * uFunTmBent - self%kMortBent
!  environmental_correction_of_zoobenthos
   tDEnvBent = max(0.0_rk,ukDIncrBent / self%cDCarrBent * sDBent*sDBent)
!  assimilation_of_zoobenthos
   tDAssBent = aDSatBent *(self%kDAssBent * uFunTmBent * sDBent - tDEnvBent)
!  consumption_of_zoobenthos
   tDConsBent = tDAssBent / self%fDAssBent
!  organic_consumption_by_zoobenthos
   tDConsPOMBent = sDPOMS / aDFoodBent * tDConsBent
!  diatoms_consumption_by_zoobenthos
   tDConsDiatBent = sDDiatS / aDFoodBent * tDConsBent
!  greens_consumption_by_zoobenthos
   tDConsGrenBent = sDGrenS / aDFoodBent * tDConsBent
!  blue-greens_consumption_by_zoobenthos
   tDConsBlueBent = sDBlueS / aDFoodBent * tDConsBent
!---------------------------------------------------------------------------
!  zoobenthos assimilation,P
!---------------------------------------------------------------------------
!  food_for_zoobenthos
   aPFoodBent = sPPOMS + aPPhytS
!  average_P/D_ratio_of_zoobenthos_food
   rPDFoodBent = aPFoodBent /(aDFoodBent+NearZero)
!  organic_P_consumption_by_zoobenthos
   tPConsPOMBent = rPDPOMS * tDConsPOMBent
!  diatom_P_consumption_by_zoobenthos
   tPConsDiatBent = rPDDiatS * tDConsDiatBent
!  greens_P_consumption_by_zoobenthos
   tPConsGrenBent = rPDGrenS * tDConsGrenBent
!  blue-greens_P_consumption_by_zoobenthos
   tPConsBlueBent = rPDBlueS * tDConsBlueBent
!  phytoplankton_P_consumption_by_zoobenthos
   tPConsPhytBent = tPConsDiatBent + tPConsGrenBent + tPConsBlueBent
!  total_P_consumption_of_zoobenthos
   tPConsBent = tPConsPOMBent + tPConsPhytBent
!  P_assim._efficiency_of_zoobenthos
   afPAssBent = min(1.0_rk,self%cPDBentRef / rPDFoodBent * self%fDAssBent)
!  P_assimilation_of_zoobenthos
   tPAssBent = afPAssBent * tPConsBent
!---------------------------------------------------------------------------
!  zoobenthos assimilation,N
!---------------------------------------------------------------------------
!  food_for_zoobenthos
   aNFoodBent = sNPOMS + aNPhytS
!  average_N/D_ratio_of_zoobenthos_food
   rNDFoodBent = aNFoodBent /(aDFoodBent+NearZero)
!  organic_N_consumption_by_zoobenthos
   tNConsPOMBent = rNDPOMS * tDConsPOMBent
!  diatom_N_consumption_by_zoobenthos
   tNConsDiatBent = rNDDiatS * tDConsDiatBent
!  greens_N_consumption_by_zoobenthos
   tNConsGrenBent = rNDGrenS * tDConsGrenBent
!  blue-greens_N_consumption_by_zoobenthos
   tNConsBlueBent = rNDBlueS * tDConsBlueBent
!  phytoplankton_N_consumption_by_zoobenthos
   tNConsPhytBent = tNConsDiatBent + tNConsGrenBent + tNConsBlueBent
!  total_N_consumption_of_zoobenthos
   tNConsBent = tNConsPOMBent + tNConsPhytBent
!  N_assim._efficiency_of_zoobenthos
   afNAssBent = min(1.0_rk,self%cNDBentRef / rNDFoodBent * self%fDAssBent)
!  N_assimilation_of_zoobenthos
   tNAssBent = afNAssBent * tNConsBent
!-----------------------------------------------------------------------
!  zoobenthos respiration and excretion, DW,P and N
!-----------------------------------------------------------------------
!  respiration_of_zoobenthos
   tDRespBent = (self%cPDBentRef / rPDBent) * self%kDRespBent * uFunTmBent * sDBent
!  P_excretion_of_zoobenthos
   tPExcrBent = (rPDBent / self%cPDBentRef) * self%kDRespBent * uFunTmBent * sPBent
!  N_excretion_of_zoobenthos
   tNExcrBent = (rNDBent / self%cNDBentRef) * self%kDRespBent * uFunTmBent * sNBent
!-----------------------------------------------------------------------
!  zoobenthos mortality
!-----------------------------------------------------------------------
!  zoobenthos_mortality_incl._environmental_correction
   tDMortBent = self%kMortBent*sDBent +(1.0_rk - aDSatBent) * tDEnvBent
!  mortality_of_zoobenthos
   tPMortBent = rPDBent * tDMortBent
!  mortality_of_zoobenthos
   tNMortBent = rNDBent * tDMortBent
!-----------------------------------------------------------------------
!  zoobenthos egestion
!-----------------------------------------------------------------------
!  egestion_of_zoobenthos_DW
   tDEgesBent = tDConsBent - tDAssBent
!  egestion_of_zoobenthos_N
   tNEgesBent = tNConsBent - tNAssBent

!-----------------------------------------------------------------------
!  total flux of Zoobenthos change to state variables
!-----------------------------------------------------------------------
!  total_flux_of_DW_in_Zoobenthos
   tDBenBent = tDMigrBent + tDAssBent - tDRespBent - tDMortBent !- tDConsFiAd
!  total_flux_of_P_in_Zoobenthos
   tPBenBent = tPMigrBent + tPAssBent - tPExcrBent - tPMortBent !- tPConsFiAd
!  total_flux_of_N_in_Zoobenthos
   tNBenBent = tNMigrBent + tNAssBent - tNExcrBent - tNMortBent ! - tNConsFiAd
!=======================================================================
! zoobenthos part relating to other modules
!=======================================================================
!-----------------------------------------------------------------------
!  Update NH4 in sediment
!-----------------------------------------------------------------------
!  part_of_died_zoobenthos_N_becoming_ammonium-N
   tNMortBentNH4 = self%fDissMortBent*tNMortBent
!  NH4_egestion_of_zoobenthos
   tNEgesBentNH4 = self%fDissEgesBent * tNEgesBent
!  total_flux_of_N_in_Pore_water_ammonium_in_lake_sediment
   tNBenNH4S = tNExcrBent + tNEgesBentNH4 + tNMortBentNH4
!-----------------------------------------------------------------------
!  Update NO3 in sediment
!-----------------------------------------------------------------------
!  total_flux_of_N_in_Pore_water_nitrate_in_lake_sediment
   tNBenNO3S = 0.0_rk
!-----------------------------------------------------------------------
!  Update PO4 in sediment
!-----------------------------------------------------------------------
!  part_of_died_zoobenthos_P_becoming_dissolved_P
   tPMortBentPO4 = self%fDissMortBent * tPMortBent
!  egestion_of_zoobenthos
   tPEgesBent = tPConsBent - tPAssBent
!  SRP_egestion_of_zoobenthos
   tPEgesBentPO4 = self%fDissEgesBent * tPEgesBent
!  total_flux_of_P_in_Pore_water_P_in_lake_sediment
   tPBenPO4S = tPExcrBent + tPEgesBentPO4 + tPMortBentPO4
!-----------------------------------------------------------------------
!  Update organic matter in sediment(DW,N,P)
!-----------------------------------------------------------------------
!  total_flux_of_DW_in_Sediment_organic matter_in_lake
   tDBenTOMS = - tDConsPOMBent + tDEgesBent + tDMortBent
   tDBenPOMS = tDBenTOMS * (1.0_rk -self%fBenDOMS)
   tDBenDOMS = tDBenTOMS * self%fBenDOMS
!  part_of_died_zoobenthos_N_becoming_organic_N
   tNMortBentTOM = (1.0_rk-self%fDissMortBent)*tNMortBent
!  detrital_N_egestion_of_zoobenthos
   tNEgesBentTOM = (1.0_rk - self%fDissEgesBent) * tNEgesBent
!  total_flux_of_N_in_Sediment_N_in_lake_sediment
   tNBenTOMS = - tNConsPOMBent + tNEgesBentTOM + tNMortBentTOM
   tNBenPOMS = tNBenTOMS * (1.0_rk -self%fBenDOMS)
   tNBenDOMS = tNBenTOMS * self%fBenDOMS
!  part_of_died_zoobenthos_P_becoming_organic_P
   tPMortBentTOM = (1.0_rk-self%fDissMortBent)*tPMortBent
!  detrital_P_egestion_of_zoobenthos
   tPEgesBentTOM = (1.0_rk - self%fDissEgesBent) * tPEgesBent
!  total_flux_of_P_in_Sediment_P_in_lake
   tPBenTOMS = - tPConsPOMBent + tPEgesBentTOM + tPMortBentTOM
   tPBenPOMS = tPBenTOMS * (1.0_rk -self%fBenDOMS)
   tPBenDOMS = tPBenTOMS * self%fBenDOMS
!  diatom_consumption_by_zoobenthos
   tSiConsDiatBent = self%cSiDDiat * tDConsDiatBent
!  total_flux_of_silica_in_sediment_organic matter
   tSiBenPOMS = tSiConsDiatBent * (1.0_rk - self%fBenDOMS)
   tSiBenDOMS = tSiConsDiatBent * self%fBenDOMS
!-----------------------------------------------------------------------
!  Update diatom in sediment(DW,N,P)
!-----------------------------------------------------------------------
!  total_flux_of_DW_in_sediment_diatoms_in_lake
   tDBenDiatS = - tDConsDiatBent
!  total_flux_of_N_in_sediment_diatoms_in_lake
   tNBenDiatS = - tNConsDiatBent
!  total_flux_of_P_in_sediment_diatoms_in_lake
   tPBenDiatS = - tPConsDiatBent
!-----------------------------------------------------------------------
!  Update green algae in sediment(DW,N,P)
!-----------------------------------------------------------------------
!  total_flux_of_DW_in_sediment_greens_in_lake
   tDBenGrenS = - tDConsGrenBent
!  total_flux_of_N_in_sediment_greens_in_lake
   tNBenGrenS = - tNConsGrenBent
!  total_flux_of_P_in_sediment_greens_in_lake
   tPBenGrenS = - tPConsGrenBent
!-----------------------------------------------------------------------
!  Update blue algae in sediment(DW,N,P)
!-----------------------------------------------------------------------
!  total_flux_of_DW_in_sediment_blue-greens_in_lake
   tDBenBlueS = - tDConsBlueBent
!  total_flux_of_N_in_sediment_blue-greens_in_lake
   tNBenBlueS = - tNConsBlueBent
!  total_flux_of_P_in_sediment_blue-greens_in_lake
   tPBenBlueS = - tPConsBlueBent
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
!  update state variables
   _SET_ODE_BEN_(self%id_sDBent,tDBenBent)
   _SET_ODE_BEN_(self%id_sPBent,tPBenBent)
   _SET_ODE_BEN_(self%id_sNBent,tNBenBent)
!-----------------------------------------------------------------------
!  Update external links
!-----------------------------------------------------------------------
!  update abiotic variables in sediment
   _SET_ODE_BEN_(self%id_NH4poolS,  tNBenNH4S)
   _SET_ODE_BEN_(self%id_NO3poolS,  tNBenNO3S)
   _SET_ODE_BEN_(self%id_PO4poolS,  tPBenPO4S)
   _SET_ODE_BEN_(self%id_DPOMpoolS, tDBenPOMS)
   _SET_ODE_BEN_(self%id_NPOMpoolS, tNBenPOMS)
   _SET_ODE_BEN_(self%id_PPOMpoolS, tPBenPOMS)
   _SET_ODE_BEN_(self%id_SiPOMpoolS,tSiBenPOMS)
   _SET_ODE_BEN_(self%id_DDOMpoolS, tDBenDOMS)
   _SET_ODE_BEN_(self%id_NDOMpoolS, tNBenDOMS)
   _SET_ODE_BEN_(self%id_PDOMpoolS, tPBenDOMS)
   _SET_ODE_BEN_(self%id_SiDOMpoolS,tSiBenDOMS)
!  update phytoplanktons in sediment
   _SET_ODE_BEN_(self%id_DfoodDiatS, tDBenDiatS)
   _SET_ODE_BEN_(self%id_NfoodDiatS, tNBenDiatS)
   _SET_ODE_BEN_(self%id_PfoodDiatS, tPBenDiatS)
   _SET_ODE_BEN_(self%id_DfoodGrenS, tDBenGrenS)
   _SET_ODE_BEN_(self%id_NfoodGrenS, tNBenGrenS)
   _SET_ODE_BEN_(self%id_PfoodGrenS, tPBenGrenS)
   _SET_ODE_BEN_(self%id_DfoodBlueS, tDBenBlueS)
   _SET_ODE_BEN_(self%id_NfoodBlueS, tNBenBlueS)
   _SET_ODE_BEN_(self%id_PfoodBlueS, tPBenBlueS)

!-----------------------------------------------------------------------
!  output diagnostic variables for external links
!-----------------------------------------------------------------------
!  Export diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBenPOMS,tDBenPOMS)

#ifdef _DEVELOPMENT_
!  update modular fluxes diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBenBent,     tDBenBent*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBenBent,     tPBenBent*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBenBent,     tNBenBent*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBenNH4S,     tNBenNH4S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBenNO3S,     tNBenNO3S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBenPO4S,     tPBenPO4S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBenPOMSflux, tDBenPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBenPOMS,     tNBenPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBenPOMS,     tPBenPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tSiBenPOMS,    tSiBenPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBenDiatS,    tDBenDiatS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBenDiatS,    tNBenDiatS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBenDiatS,    tPBenDiatS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBenGrenS,    tDBenGrenS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBenGrenS,    tNBenGrenS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBenGrenS,    tPBenGrenS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBenBlueS,    tDBenBlueS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBenBlueS,    tNBenBlueS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBenBlueS,    tPBenBlueS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAssFiAd,     tDAssFiAd/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAssFiAd,     tNAssFiAd/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAssFiAd,     tPAssFiAd/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBenNH4W,     tNEgesFiAdNH4/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBenPO4W,     tPEgesFiAdPO4/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDBenPOMW,     tDEgesFiAd/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNBenPOMW,     tNEgesFiAdPOM/dz*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPBenPOMW,     tPEgesFiAdPOM/dz*secs_pr_day)
#endif
!  Spatial loop end
   _FABM_HORIZONTAL_LOOP_END_
   end subroutine do_bottom

!
!EOC
!-----------------------------------------------------------------------
   end module pclake_zoobenthos
!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
