#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module pclake_phytoplankton_sediment
! !USES:
   use fabm_types
   use pclake_utility, ONLY:uFunTmBio

   implicit none

!  default: all is private.
   private

! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_phytoplankton_sediment
!  local state variable identifiers
!  id_sDDiatS,id_sDGrenS,id_sDBlueS: phytoplankton concentration in dry-weight, gDW/m**2
!  id_sPDiatS,id_sPGrenS,id_sPBlueS: phytoplankton concentration in nitrogen element, gN/m**2
!  id_sNDiatS,id_sNGrenS,id_sNBlueS: phytoplankton concentration in phosphorus element, gP/m**2
   type (type_bottom_state_variable_id) :: id_sDDiatS,id_sPDiatS,id_sNDiatS
   type (type_bottom_state_variable_id) :: id_sDGrenS,id_sPGrenS,id_sNGrenS
   type (type_bottom_state_variable_id) :: id_sDBlueS,id_sPBlueS,id_sNBlueS
!  diagnostic variables for local output
!  id_oSiDiatS, diatom concentration in silica element, gSi/m**2
!  id_rPDPhytS, P/D ratio of settled phytoplankton
!  id_rNDPhytS, N/D rartio of settled phytoplankton
   type (type_horizontal_diagnostic_variable_id)       :: id_oSiDiatS
   type (type_horizontal_diagnostic_variable_id)       :: id_rPDPhytS,id_rNDPhytS
!  diagnostic variables used by external dependencies
   type (type_horizontal_diagnostic_variable_id)       :: id_tDPrimPOMS
#ifdef _DEVELOPMENT_
!  diagnostic variables for modular fluxes
   type (type_horizontal_diagnostic_variable_id)       :: id_tDPrimDiatS,id_tNPrimDiatS,id_tPPrimDiatS
   type (type_horizontal_diagnostic_variable_id)       :: id_tDPrimGrenS,id_tNPrimGrenS,id_tPPrimGrenS
   type (type_horizontal_diagnostic_variable_id)       :: id_tDPrimBlueS,id_tNPrimBlueS,id_tPPrimBlueS
   type (type_horizontal_diagnostic_variable_id)       :: id_tPPrimPO4S,id_tNPrimNO3S,id_tNPrimNH4S
   type (type_horizontal_diagnostic_variable_id)       :: id_tNPrimPOMS,id_tPPrimPOMS,id_tSiPrimPOMS,id_tSiExcrDiatS
!  organic fluxes will be named differently due to they are used by
!  external dependencies
   type (type_horizontal_diagnostic_variable_id)       :: id_tDPrimPOMSflux
#endif
!  state dependencies identifiers
   type (type_bottom_state_variable_id) :: id_PO4poolS,id_NO3poolS,id_NH4poolS
   type (type_bottom_state_variable_id) :: id_DPOMpoolS,id_NPOMpoolS,id_PPOMpoolS,id_SiPOMpoolS
   type (type_bottom_state_variable_id) :: id_DDOMpoolS,id_NDOMpoolS,id_PDOMpoolS,id_SiDOMpoolS
   type (type_state_variable_id)        :: id_SiO2poolW
!  environmental dependencies
   type (type_dependency_id)                :: id_uTm,id_dz
!  Model parameters
!  temperature parameters
   real(rk)   :: cSigTmDiat,cTmOptDiat
   real(rk)   :: cSigTmBlue,cTmOptBlue,cSigTmGren,cTmOptGren
!  diatoms related parameters
   real(rk)   :: kMortDiatS,kDRespDiat,kMortGrenS,kDRespGren,kMortBlueS,kDRespBlue
!  nutrient ratios parameter
   real(rk)   :: cNDDiatMin,cPDDiatMin,cNDGrenMin,cPDGrenMin,cNDBlueMin,cPDBlueMin
   real(rk)   :: cNDDiatMax,cPDDiatMax,cNDGrenMax,cPDGrenMax,cNDBlueMax,cPDBlueMax
!  exchange process related parameters
   real(rk)   :: fDissMortPhyt,cSiDDiat
!  minimum state variable values
   real(rk)   :: cDBlueMinS,cDGrenMinS,cDDiatMinS
!  fraction of dissolved organic matter from phytoplankton
   real(rk)   :: fPrimDOMS


   contains

   procedure initialize
   procedure do_bottom

   end type type_pclake_phytoplankton_sediment

!  private data members (API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero=0.000000000000000000000000000000001_rk
!  Lowest phytoplankton value in sediment
   real(rk),parameter :: PhyZeroS=0.0000001_rk
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
! !INPUT PARAMETERS:
   class (type_pclake_phytoplankton_sediment), intent(inout), target :: self
   integer,                     intent(in)            :: configunit

!EOP
!-----------------------------------------------------------------------
!BOC
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
   call self%get_parameter(self%cSigTmDiat,   'cSigTmDiat',   'degree C', 'temperature constant for diatoms (sigma)',                  default=20.0_rk)
   call self%get_parameter(self%cTmOptDiat,   'cTmOptDiat',   'degree C', 'optimal temperature for diatoms',                           default=18.0_rk)
   call self%get_parameter(self%cSigTmBlue,   'cSigTmBlue',   'degree C', 'temperature constant for blue-greens (sigma)',              default=12.0_rk)
   call self%get_parameter(self%cTmOptBlue,   'cTmOptBlue',   'degree C', 'optimal temperature for blue-greens',                       default=25.0_rk)
   call self%get_parameter(self%cSigTmGren,   'cSigTmGren',   'degree C', 'temperature constant for greens (sigma)',                   default=15.0_rk)
   call self%get_parameter(self%cTmOptGren,   'cTmOptGren',   'degree C', 'optimal temperature for greens',                            default=25.0_rk)
   call self%get_parameter(self%kMortDiatS,   'kMortDiatS',   'd-1',      'mortality rate of settled diatoms',                         default=0.05_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespDiat,   'kDRespDiat',   'd-1',      'maintenance respiration constant for diatoms',              default=0.1_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cNDDiatMin,   'cNDDiatMin',   'mgN/mgDW', 'minimum N/DW ratio diatoms',                                default=0.01_rk)
   call self%get_parameter(self%cPDDiatMin,   'cPDDiatMin',   'mgP/mgDW', 'minimum P/DW ratio diatoms',                                default=0.0005_rk)
   call self%get_parameter(self%kMortGrenS,   'kMortGrenS',   'd-1',      'mortality rate of settled greens',                          default=0.05_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespGren,   'kDRespGren',   'd-1',      'maintenance respiration constant for greens',               default=0.075_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cNDGrenMin,   'cNDGrenMin',   'mgN/mgDW', 'minimum N/DW ratio greens',                                 default=0.02_rk)
   call self%get_parameter(self%cPDGrenMin,   'cPDGrenMin',   'mgP/mgDW', 'minimum P/DW ratio greens',                                 default=0.0015_rk)
   call self%get_parameter(self%kMortBlueS,   'kMortBlueS',   'd-1',      'mortality rate of settled blue-greens',                     default=0.2_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kDRespBlue,   'kDRespBlue',   'd-1',      'maintenance respiration constant for blue-greens',          default=0.03_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cNDBlueMin,   'cNDBlueMin',   'mgN/mgDW', 'minimum N/DW ratio blue-greens',                            default=0.03_rk)
   call self%get_parameter(self%cPDBlueMin,   'cPDBlueMin',   'mgP/mgDW', 'minimum P/DW ratio blue-greens',                            default=0.0025_rk)
   call self%get_parameter(self%fDissMortPhyt,'fDissMortPhyt','[-]',      'soluble nutrient fraction of died algae',                   default=0.2_rk)
   call self%get_parameter(self%cSiDDiat,     'cSiDDiat',     'mgSi/mgDW','Si/D ratio of diatoms',                                     default=0.15_rk)
   call self%get_parameter(self%cNDBlueMax,   'cNDBlueMax',   'mgN/mgDW', 'maximum N/DW ratio blue-greens',                            default=0.15_rk)
   call self%get_parameter(self%cNDDiatMax,   'cNDDiatMax',   'mgN/mgDW', 'maximum N/DW ratio diatoms',                                default=0.05_rk)
   call self%get_parameter(self%cNDGrenMax,   'cNDGrenMax',   'mgN/mgDW', 'maximum N/DW ratio greens',                                 default=0.1_rk)
   call self%get_parameter(self%cPDBlueMax,   'cPDBlueMax',   'mgP/mgDW', 'maximum P/DW ratio blue-greens',                            default=0.025_rk)
   call self%get_parameter(self%cPDDiatMax,   'cPDDiatMax',   'mgP/mgDW', 'maximum P/DW ratio diatoms',                                default=0.005_rk)
   call self%get_parameter(self%cPDGrenMax,   'cPDGrenMax',   'mgP/mgDW', 'maximum P/DW ratio greens',                                 default=0.015_rk)
!  the user defined minumun value for state variables
   call self%get_parameter(self%cDBlueMinS,   'cDBlueMinS',   'gDW/m2',   'minimum blue-green algae biomass in system',                default=0.00001_rk)
   call self%get_parameter(self%cDGrenMinS,   'cDGrenMinS',   'gDW/m2',   'minimum green algae biomass in system',                     default=0.00001_rk)
   call self%get_parameter(self%cDDiatMinS,   'cDDiatMinS',   'gDW/m2',   'minimum diatom biomass in system',                          default=0.00001_rk)
   call self%get_parameter(self%fPrimDOMS, 'fPrimDOMS', '[-]',      'fraction of DOM from settled phytoplankton',                      default=0.5_rk)
!  Register local state variable
   call self%register_state_variable(self%id_sDDiatS,'sDDiatS','gDW m-2','diatoms DW',     &
                                    initial_value=0.001_rk,minimum= self%cDDiatMinS)
   call self%register_state_variable(self%id_sPDiatS,'sPDiatS','gP m-2','diatoms P',     &
                                    initial_value=0.00001_rk,minimum= self%cDDiatMinS * self%cPDDiatMin)
   call self%register_state_variable(self%id_sNDiatS,'sNDiatS','gN m-2','diatoms N',     &
                                    initial_value=0.0001_rk,minimum= self%cDDiatMinS * self%cNDDiatMin)
   call self%register_state_variable(self%id_sDGrenS,'sDGrenS','gDW m-2','greens DW',     &
                                    initial_value=0.001_rk,minimum= self%cDGrenMinS)
   call self%register_state_variable(self%id_sPGrenS,'sPGrenS','gP m-2','greens P',     &
                                    initial_value=0.00001_rk,minimum= self%cDGrenMinS * self%cPDGrenMin)
   call self%register_state_variable(self%id_sNGrenS,'sNGrenS','gN m-2','greens N',     &
                                    initial_value=0.0001_rk,minimum=self%cDGrenMinS * self%cNDGrenMin)
   call self%register_state_variable(self%id_sDBlueS,'sDBlueS','gDW m-2','blue-greens DW',     &
                                    initial_value=0.001_rk,minimum=self%cDBlueMinS)
   call self%register_state_variable(self%id_sPBlueS,'sPBlueS','gP m-2','blue-greens P',     &
                                    initial_value=0.00001_rk,minimum=self%cDBlueMinS * self%cPDBlueMin)
   call self%register_state_variable(self%id_sNBlueS,'sNBlueS','gN m-2','blue-greens N',     &
                                    initial_value=0.0001_rk,minimum=self%cDBlueMinS * self%cNDBlueMin)
!  register diagnostic variables
   call self%register_diagnostic_variable(self%id_oSiDiatS,  'oSiDiatS',  'g m-2 s-1','oSiDiatS',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rPDPhytS,  'rPDPhytS',  '[-]',      'rPDPhytS',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDPhytS,  'rNDPhytS',  '[-]',      'rNDPhytS',  output=output_instantaneous)
!  register diagnostic variable for external usage
   call self%register_diagnostic_variable(self%id_tDPrimPOMS,'tDPrimPOMS','gDW m-2 s-1', 'POM from settled phytoplankton',  output=output_none)
#ifdef _DEVELOPMENT_
!  register diagnostic variables for modular fluxes
!  fluxes for local state variables
   call self%register_diagnostic_variable(self%id_tDPrimDiatS,    'tDPrimDiatS',    'g m-2 s-1', 'phytoplankton_sediment_DDiat_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNPrimDiatS,    'tNPrimDiatS',    'g m-2 s-1', 'phytoplankton_sediment_NDiat_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPPrimDiatS,    'tPPrimDiatS',    'g m-2 s-1', 'phytoplankton_sediment_PDiat_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDPrimGrenS,    'tDPrimGrenS',    'g m-2 s-1', 'phytoplankton_sediment_DGren_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNPrimGrenS,    'tNPrimGrenS',    'g m-2 s-1', 'phytoplankton_sediment_NGren_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPPrimGrenS,    'tPPrimGrenS',    'g m-2 s-1', 'phytoplankton_sediment_PGren_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDPrimBlueS,    'tDPrimBlueS',    'g m-2 s-1', 'phytoplankton_sediment_DBlue_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNPrimBlueS,    'tNPrimBlueS',    'g m-2 s-1', 'phytoplankton_sediment_NBlue_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPPrimBlueS,    'tPPrimBlueS',    'g m-2 s-1', 'phytoplankton_sediment_PBlue_change',  output=output_instantaneous)
!  fluxes for other modules, mainly abiotic sediment
   call self%register_diagnostic_variable(self%id_tPPrimPO4S,     'tPPrimPO4S',     'g m-2 s-1', 'phytoplankton_sediment_PO4S_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNPrimNO3S,     'tNPrimNO3S',     'g m-2 s-1', 'phytoplankton_sediment_NO3S_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNPrimNH4S,     'tNPrimNH4S',     'g m-2 s-1', 'phytoplankton_sediment_NH4S_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDPrimPOMSflux, 'tDPrimPOMSflux', 'g m-2 s-1', 'phytoplankton_sediment_DPOMS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNPrimPOMS,     'tNPrimPOMS',     'g m-2 s-1', 'phytoplankton_sediment_NPOMS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPPrimPOMS,     'tPPrimPOMS',     'g m-2 s-1', 'phytoplankton_sediment_PPOMS_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tSiPrimPOMS,    'tSiPrimPOMS',    'g m-2 s-1', 'phytoplankton_sediment_SiPOMS_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tSiExcrDiatS,   'tSiExcrDiatS',   'g m-2 s-1', 'phytoplankton_sediment_SiO2_change',   output=output_instantaneous)
#endif
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNDiatS)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNGrenS)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNBlueS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPDiatS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPGrenS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPBlueS)
!   call self%add_to_aggregate_variable(standard_variables%total_silicate,  self%id_oSiDiatS)
!  register state variables dependencies
   call self%register_state_dependency(self%id_PO4poolS,     'PO4_pool_sediment',                'g m-2', 'PO4 pool in sediment')
   call self%register_state_dependency(self%id_NO3poolS,     'NO3_pool_sediment',                'g m-2', 'NO3 pool in sediment')
   call self%register_state_dependency(self%id_NH4poolS,     'NH4_pool_sediment',                'g m-2', 'NH4 pool in sediment')
   call self%register_state_dependency(self%id_DPOMpoolS,    'POM_DW_pool_sediment',        'g m-2', 'POM DW pool in sediment')
   call self%register_state_dependency(self%id_NPOMpoolS,    'POM_N_pool_sediment',         'g m-2', 'POM N pool in sediment')
   call self%register_state_dependency(self%id_PPOMpoolS,    'POM_P_pool_sediment',         'g m-2', 'POM P pool in sediment')
   call self%register_state_dependency(self%id_SiPOMpoolS,   'POM_Si_pool_sediment',        'g m-2', 'POM Si pool in sediment')
   call self%register_state_dependency(self%id_SiO2poolW,    'SiO2_pool_water',                  'g m-3', 'SiO2 pool in water')
   call self%register_state_dependency(self%id_DDOMpoolS, 'DOM_DW_pool_sediment',   'g m-2', 'DOM DW in sediment')
   call self%register_state_dependency(self%id_NDOMpoolS, 'DOM_N_pool_sediment',    'g m-2', 'DOM N in sediment')
   call self%register_state_dependency(self%id_PDOMpoolS, 'DOM_P_pool_sediment',    'g m-2', 'DOM P in sediment')
   call self%register_state_dependency(self%id_SiDOMpoolS,'DOM_Si_pool_sediment',   'g m-2', 'DOM Si in sediment')
!  register environmental dependencies
   call self%register_dependency(self%id_uTm,standard_variables%temperature)
   call self%register_dependency(self%id_dz,standard_variables%cell_thickness)

   return

   end subroutine initialize
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
! !INPUT PARAMETERS:
   class (type_pclake_phytoplankton_sediment), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
! !LOCAL VARIABLES:
!  state variables value carriers
   real(rk)                   :: sDDiatS,sPDiatS,sNDiatS
   real(rk)                   :: sDGrenS,sPGrenS,sNGrenS
   real(rk)                   :: sDBlueS,sPBlueS,sNBlueS
   real(rk)                   :: oSiDiatS,rPDPhytS,rNDPhytS
!  environmental dependencies carrier
   real(rk)                   :: uTm,dz
!  external links value carrier
   real(rk)                  :: sPO4S,sNO3S,sNH4S,sSiO2W
   real(rk)                  :: sDPOMS,sNPOMS,sPPOMS,sSiPOMS
!  Nutrients rations
   real(rk)   :: rPDDiatS,rNDDiatS
   real(rk)   :: rPDGrenS,rNDGrenS
   real(rk)   :: rPDBlueS,rNDBlueS
!  Temperature variables
   real(rk)   :: uFunTmDiat,uFunTmBlue,uFunTmGren
!  Diatom dry-weight variables
   real(rk)   :: tDPrimDiatS,tDMortDiatS,tDRespDiatS,ukDRespTmDiat
!  Diatom nitrogen variables
   real(rk)   :: tNPrimDiatS,tNMortDiatS,tNExcrDiatS
!  Diatom phosphorus variables
   real(rk)   :: tPPrimDiatS,tPMortDiatS,tPExcrDiatS
!  Green algae dry-weight variables
   real(rk)   :: tDPrimGrenS,tDMortGrenS,tDRespGrenS,ukDRespTmGren
!  Green algae Nitrogen variables
   real(rk)   :: tNPrimGrenS,tNMortGrenS,tNExcrGrenS
!  Green algae phosphorus variables
   real(rk)   :: tPPrimGrenS,tPMortGrenS,tPExcrGrenS
!  Blue-green algae dry-weight variables
   real(rk)   :: tDPrimBlueS,tDMortBlueS,tDRespBlueS,ukDRespTmBlue
!  Blue-green algae nitrogen variables
   real(rk)   :: tNPrimBlueS,tNMortBlueS,tNExcrBlueS
!  Blue-green algae phosphorus variables
   real(rk)   :: tPPrimBlueS,tPMortBlueS,tPExcrBlueS
!  Exchange related variables
   real(rk)   :: tNPrimNH4S,tNExcrPhytS,tNMortPhytNH4S,tNMortPhytS
   real(rk)   :: tNPrimNO3S
   real(rk)   :: tPPrimPO4S,tPExcrPhytS,tPMortPhytPO4S,tPMortPhytS
   real(rk)   :: tDPrimPOMS,tDMortPhytS
   real(rk)   :: tNPrimPOMS,tNMortPhytTOMS
   real(rk)   :: tPPrimPOMS,tPMortPhytTOMS
   real(rk)   :: tSiPrimPOMS,tSiMortDiatS
   real(rk)   :: tSiExcrDiatS,tSiPrimDOMS
   real(rk)   :: tDPrimDOMS,tNPrimDOMS,tPPrimDOMS
!  Enter spatial loops (if any)
   _FABM_HORIZONTAL_LOOP_BEGIN_
!-----------------------------------------------------------------------
!  Retrieve current (local) state variable values.
!-----------------------------------------------------------------------
   _GET_HORIZONTAL_(self%id_sDDiatS,sDDiatS)
   _GET_HORIZONTAL_(self%id_sNDiatS,sNDiatS)
   _GET_HORIZONTAL_(self%id_sPDiatS,sPDiatS)
   _GET_HORIZONTAL_(self%id_sDGrenS,sDGrenS)
   _GET_HORIZONTAL_(self%id_sNGrenS,sNGrenS)
   _GET_HORIZONTAL_(self%id_sPGrenS,sPGrenS)
   _GET_HORIZONTAL_(self%id_sDBlueS,sDBlueS)
   _GET_HORIZONTAL_(self%id_sNBlueS,sNBlueS)
   _GET_HORIZONTAL_(self%id_sPBlueS,sPBlueS)
!-----------------------------------------------------------------------
!  Retrieve dependencies  value
!-----------------------------------------------------------------------
!  Retrieve state dependencie value
   _GET_HORIZONTAL_(self%id_PO4poolS,sPO4S)
   _GET_HORIZONTAL_(self%id_NO3poolS,sNO3S)
   _GET_HORIZONTAL_(self%id_NH4poolS,sNH4S)
   _GET_HORIZONTAL_(self%id_DPOMpoolS,sDPOMS)
   _GET_HORIZONTAL_(self%id_NPOMpoolS,sNPOMS)
   _GET_HORIZONTAL_(self%id_PPOMpoolS,sPPOMS)
   _GET_HORIZONTAL_(self%id_SiPOMpoolS,sSiPOMS)
   _GET_(self%id_SiO2poolW,sSiO2W)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_(self%id_dz,dz)
!-----------------------------------------------------------------------
!  Current local nutrients ratios (check the current state)
!-----------------------------------------------------------------------
!  Local status
!  P/D_ratio_of_Diatom
   rPDDiatS = sPDiatS /(sDDiatS+NearZero)
!  N/D_ratio_of_Diatom
   rNDDiatS = sNDiatS /(sDDiatS+NearZero)
!  P/D_ratio_of_Green_Algae
   rPDGrenS = sPGrenS /(sDGrenS+NearZero)
!  N/D_ratio_of_Green_Algae
   rNDGrenS = sNGrenS /(sDGrenS+NearZero)
!  P/D_ratio_of_Blue_Algae
   rPDBlueS = sPBlueS /(sDBlueS+NearZero)
!  N/D_ratio_of_Blue_Algae
   rNDBlueS = sNBlueS /(sDBlueS+NearZero)
!  Settled diatom s in silica element
   oSiDiatS= self%cSiDDiat*sDDiatS
!  P/D_ratio_of_settled_phytoplankton
   rPDPhytS= (sPDiatS+sPGrenS+sPBlueS)/(sDDiatS+sDGrenS+sDBlueS+NearZero)
   rNDPhytS= (sNDiatS+sNGrenS+sNBlueS)/(sDDiatS+sDGrenS+sDBlueS+NearZero)
!-----------------------------------------------------------------------
!  Temperature functions for pelagic phytoplankton
!-----------------------------------------------------------------------
!  temperature_function_of_Diatom
   uFunTmDiat = uFunTmBio(uTm,self%cSigTmDiat,self%cTmOptDiat)
!  temperature_function_of_Blue_algae
   uFunTmBlue = uFunTmBio(uTm,self%cSigTmBlue,self%cTmOptBlue)
!  temperature_function_of_Green_Algae
   uFunTmGren = uFunTmBio(uTm,self%cSigTmGren,self%cTmOptGren)
!-----------------------------------------------------------------------
!  Dry-weight change of diatoms in sediment
!-----------------------------------------------------------------------
!  temp._corrected_respiration_constant_of_Algae
   ukDRespTmDiat = self%kDRespDiat * uFunTmDiat
!  respiration_of_sediment_Algae
   tDRespDiatS = ukDRespTmDiat * sDDiatS
!  mortality_in_sed.
   tDMortDiatS = self%kMortDiatS * sDDiatS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tDPrimDiatS = - tDMortDiatS - tDRespDiatS
!-----------------------------------------------------------------------
!  Nitrogen change of diatoms in sediment
!-----------------------------------------------------------------------
#ifdef _V509_
   tNExcrDiatS = rNDDiatS /(self%cNDDiatMin + rNDDiatS) * rNDDiatS * tDRespDiatS
   tPExcrDiatS = rPDDiatS /(self%cPDDiatMin + rPDDiatS) * rPDDiatS * tDRespDiatS
   tNExcrGrenS = rNDGrenS /(self%cNDGrenMin + rNDGrenS) * rNDGrenS * tDRespGrenS
   tPExcrGrenS = rPDGrenS /(self%cPDGrenMin + rPDGrenS) * rPDGrenS * tDRespGrenS
   tNExcrBlueS = rNDBlueS /(self%cNDBlueMin + rNDBlueS) * rNDBlueS * tDRespBlueS
   tPExcrBlueS = rPDBlueS /(self%cPDBlueMin + rPDBlueS) * rPDBlueS * tDRespBlueS
#endif
!  N_excretion_of_algae_in_sediment
   tNExcrDiatS = (2.0_rk * rNDDiatS) /(self%cNDDiatMax + rNDDiatS) * rNDDiatS * tDRespDiatS
!  N_mortality_of_algae_in_sediment
   tNMortDiatS = self%kMortDiatS * sNDiatS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tNPrimDiatS = - tNMortDiatS - tNExcrDiatS
!-----------------------------------------------------------------------
!  phosphorus change of diatoms in sediment
!-----------------------------------------------------------------------
!  P_excretion_of_algae_in_sediment 
   tPExcrDiatS = (2.0_rk * rPDDiatS) /(self%cPDDiatMax + rPDDiatS) * rPDDiatS * tDRespDiatS
!  P_mortality_of_algae_in_sediment
   tPMortDiatS = self%kMortDiatS * sPDiatS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tPPrimDiatS = -tPMortDiatS - tPExcrDiatS
!-----------------------------------------------------------------------
!  Dry-weight change of green algae in sediment
!-----------------------------------------------------------------------
!  temp._corrected_respiration_constant_of_Algae
   ukDRespTmGren = self%kDRespGren * uFunTmGren
!  respiration_of_sediment_Algae
   tDRespGrenS = ukDRespTmGren * sDGrenS
!  mortality_in_sed.
   tDMortGrenS = self%kMortGrenS * sDGrenS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tDPrimGrenS = - tDMortGrenS - tDRespGrenS
!-----------------------------------------------------------------------
!  Nitrogen change of green algae in sediment
!-----------------------------------------------------------------------
!  N_excretion_of_algae_in_sediment
   tNExcrGrenS = (2.0_rk * rNDGrenS) /(self%cNDGrenMax + rNDGrenS) * rNDGrenS * tDRespGrenS
!  N_mortality_of_algae_in_sediment
   tNMortGrenS = self%kMortGrenS * sNGrenS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tNPrimGrenS = - tNMortGrenS - tNExcrGrenS
!-----------------------------------------------------------------------
!  Phosphorus change of green algae in sediment
!-----------------------------------------------------------------------
!  P_excretion_of_algae_in_sediment
   tPExcrGrenS = (2.0_rk * rPDGrenS) /(self%cPDGrenMax + rPDGrenS) * rPDGrenS * tDRespGrenS
!  P_mortality_of_algae_in_sediment
   tPMortGrenS = self%kMortGrenS * sPGrenS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tPPrimGrenS = - tPMortGrenS - tPExcrGrenS
!-----------------------------------------------------------------------
!  Dry-weight change of blue-green algae in sediment
!-----------------------------------------------------------------------
!  temp._corrected_respiration_constant_of_Algae
   ukDRespTmBlue = self%kDRespBlue * uFunTmBlue
!  respiration_of_sediment_Algae
   tDRespBlueS = ukDRespTmBlue * sDBlueS
!  mortality_in_sed.
   tDMortBlueS = self%kMortBlueS * sDBlueS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tDPrimBlueS = - tDMortBlueS - tDRespBlueS
!-----------------------------------------------------------------------
!  Nitrogen change of blue-green algae in sediment
!-----------------------------------------------------------------------
!  N_excretion_of_algae_in_sediment
   tNExcrBlueS = (2.0_rk * rNDBlueS) /(self%cNDBlueMax + rNDBlueS) * rNDBlueS * tDRespBlueS
!  N_mortality_of_algae_in_sediment
   tNMortBlueS = self%kMortBlueS * sNBlueS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tNPrimBlueS = - tNMortBlueS - tNExcrBlueS
!-----------------------------------------------------------------------
!  Phosphorus change of blue-green algae in sediment
!-----------------------------------------------------------------------
!  P_excretion_of_algae_in_sediment
   tPExcrBlueS = (rPDBlueS *2.0_rk)/(self%cPDBlueMax + rPDBlueS) * rPDBlueS * tDRespBlueS
!  P_mortality_of_algae_in_sediment
   tPMortBlueS = self%kMortBlueS * sPBlueS
!  total_flux_from_PRIM_module_to_sediment_Algae
   tPPrimBlueS = - tPMortBlueS - tPExcrBlueS
!-----------------------------------------------------------------------
!  exchange to other modules
!-----------------------------------------------------------------------
!  sNH4S exchange
!-----------------------------------------------------------------------
!  total_phytoplankton_mortality
   tNMortPhytS = tNMortDiatS + tNMortGrenS + tNMortBlueS
!  ammonium_flux_from_died_Algae
   tNMortPhytNH4S = self%fDissMortPhyt * tNMortPhytS
!  total_N_excretion_sediment_phytoplankton
   tNExcrPhytS = tNExcrDiatS + tNExcrGrenS + tNExcrBlueS
!  Pore_water_ammonium
   tNPrimNH4S = tNExcrPhytS + tNMortPhytNH4S
!-----------------------------------------------------------------------
!  sNO3S exchange, the phytoplankton in the sediment are considered
!  without uptaking process
!-----------------------------------------------------------------------
!  Pore_water_nitrate
   tNPrimNO3S = 0.0_rk
!-----------------------------------------------------------------------
!  sPO4S exchange
!-----------------------------------------------------------------------
!  total_phytoplankton_mortality
   tPMortPhytS = tPMortDiatS + tPMortGrenS + tPMortBlueS
!  soluble_P_flux_from_died_Algae
   tPMortPhytPO4S = self%fDissMortPhyt * tPMortPhytS
!  total_P_excretion_sediment_phytoplankton
   tPExcrPhytS = tPExcrDiatS + tPExcrGrenS + tPExcrBlueS
!  Pore_water_P
   tPPrimPO4S = tPExcrPhytS + tPMortPhytPO4S
!-----------------------------------------------------------------------
!  sediment organic matter exchange, DW
!-----------------------------------------------------------------------
!  mortality_of_algae_on_bottom
   tDMortPhytS = tDMortDiatS + tDMortGrenS + tDMortBlueS
!  Flux_to_sediment_organic matter
   tDPrimPOMS = tDMortPhytS * (1.0_rk - self%fPrimDOMS)
   tDPrimDOMS = tDMortPhytS * self%fPrimDOMS
!-----------------------------------------------------------------------
!  sediment organic matter exchange, N
!-----------------------------------------------------------------------
!  organic_N_flux_from_died_Algae
   tNMortPhytTOMS = tNMortPhytS - tNMortPhytNH4S
!  Flux_to_sediment_organic matter
   tNPrimPOMS = tNMortPhytTOMS * (1.0_rk - self%fPrimDOMS)
   tNPrimDOMS = tNMortPhytTOMS * self%fPrimDOMS
!-----------------------------------------------------------------------
!  sediment organic matter exchange, P
!-----------------------------------------------------------------------
!  organic_P_flux_from_died_Algae
   tPMortPhytTOMS = tPMortPhytS - tPMortPhytPO4S
!  Flux_to_sediment_organic matter
   tPPrimPOMS = tPMortPhytTOMS * (1.0_rk - self%fPrimDOMS)
   tPPrimDOMS = tPMortPhytTOMS * self%fPrimDOMS
!-----------------------------------------------------------------------
!  sediment organic matter exchange, Si
!-----------------------------------------------------------------------
!  mortality_of_bottom_Algae
   tSiMortDiatS = self%cSiDDiat * tDMortDiatS
!  Sediment_organic matter
   tSiPrimPOMS = tSiMortDiatS* (1.0_rk - self%fPrimDOMS)
   tSiPrimDOMS = tSiMortDiatS * self%fPrimDOMS
!-----------------------------------------------------------------------
!  sSiO2W exchange
!-----------------------------------------------------------------------
!  Si_excretion_of_bottom_Algae
   tSiExcrDiatS = self%cSiDDiat * tDRespDiatS
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
   _SET_ODE_BEN_(self%id_sDDiatS, tDPrimDiatS)
   _SET_ODE_BEN_(self%id_sNDiatS, tNPrimDiatS)
   _SET_ODE_BEN_(self%id_sPDiatS, tPPrimDiatS)
   _SET_ODE_BEN_(self%id_sDGrenS, tDPrimGrenS)
   _SET_ODE_BEN_(self%id_sNGrenS, tNPrimGrenS)
   _SET_ODE_BEN_(self%id_sPGrenS, tPPrimGrenS)
   _SET_ODE_BEN_(self%id_sDBlueS, tDPrimBlueS)
   _SET_ODE_BEN_(self%id_sNBlueS, tNPrimBlueS)
   _SET_ODE_BEN_(self%id_sPBlueS, tPPrimBlueS)
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
   _SET_ODE_BEN_(self%id_PO4poolS,          tPPrimPO4S)
   _SET_ODE_BEN_(self%id_NO3poolS,          tNPrimNO3S)
   _SET_ODE_BEN_(self%id_NH4poolS,          tNPrimNH4S)
   _SET_ODE_BEN_(self%id_DPOMpoolS,         tDPrimPOMS)
   _SET_ODE_BEN_(self%id_NPOMpoolS,         tNPrimPOMS)
   _SET_ODE_BEN_(self%id_PPOMpoolS,         tPPrimPOMS)
   _SET_ODE_BEN_(self%id_SiPOMpoolS,        tSiPrimPOMS)
   _SET_BOTTOM_EXCHANGE_(self%id_SiO2poolW, tSiExcrDiatS)
   _SET_ODE_BEN_(self%id_DDOMpoolS,      tDPrimDOMS)
   _SET_ODE_BEN_(self%id_NDOMpoolS,      tNPrimDOMS)
   _SET_ODE_BEN_(self%id_PDOMpoolS,      tPPrimDOMS)
   _SET_ODE_BEN_(self%id_SiDOMpoolS,     tSiPrimDOMS)
!-----------------------------------------------------------------------
!  Output dependent diagnostic variables for other modules
!-----------------------------------------------------------------------
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_oSiDiatS,   oSiDiatS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rPDPhytS,   rPDPhytS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rNDPhytS,   rNDPhytS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDPrimPOMS, tDPrimPOMS)
#ifdef _DEVELOPMENT_
!  output diagnostic variables for modular fluxes
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDPrimDiatS,    tDPrimDiatS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNPrimDiatS,    tNPrimDiatS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPPrimDiatS,    tPPrimDiatS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDPrimGrenS,    tDPrimGrenS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNPrimGrenS,    tNPrimGrenS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPPrimGrenS,    tPPrimGrenS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDPrimBlueS,    tDPrimBlueS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNPrimBlueS,    tNPrimBlueS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPPrimBlueS,    tPPrimBlueS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPPrimPO4S,     tPPrimPO4S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNPrimNO3S,     tNPrimNO3S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNPrimNH4S,     tNPrimNH4S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDPrimPOMSflux, tDPrimPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNPrimPOMS,     tNPrimPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPPrimPOMS,     tPPrimPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tSiPrimPOMS,    tSiPrimPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tSiExcrDiatS,   tSiExcrDiatS/dz*secs_pr_day)
#endif
   _FABM_HORIZONTAL_LOOP_END_

   end subroutine do_bottom
!-----------------------------------------------------------------------
! Spatial loop end
!-----------------------------------------------------------------------
!
!EOP
!-----------------------------------------------------------------------

   end module pclake_phytoplankton_sediment

!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
