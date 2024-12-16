#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module pclake_abiotic_water
! !USES:
   use fabm_types
   use fabm_standard_variables
   use pclake_utility, ONLY: uFunTmAbio

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_abiotic_water
!  local state variable identifiers
!  sDIMW: inorganic matter concentration, in dry-weight, gDW/m**3
!  sDPOMW, particulate organic matter concentration, in dry-weight, gDW/m**3
!  sNPOMW, particulate organic matter concentration, in nitrogen element, gN/m**3
!  sPPOMW, particulate organic matter concentration, in phosphorus element, gP/m**3
!  sSiPOMW, particulate organic matter concentration, in silicon element, gSi/m**3
!  sDDOMW, dissolved organic matter concentration, in dry-weight, gDW/m**3
!  sNDOMW, dissolved organic matter concentration, in nitrogen element, gN/m**3
!  sPDOMW, dissolved organic matter concentration, in phosphorus element, gP/m**3
!  sSidOMW, dissolved organic matter concentration, in silicon element, gSi/m**3
!  sNH4W,  ammonium concentration, in nitrogen element, gN/m**3
!  sNO3W,  nitrate concentration, in nitrogen element, gN/m**3
!  sPO4W,  phosphate concentration, in phosphorus element, gP/m**3
!  sPAIMW, absorbed phosphorus concentration, in phosphorus element, gP/m**3
!  sSiO2W, silica dioxide concentration, in silica element, gSi/m**3
!  sO2W,   dissolved oxygen, in O2 molecule, gO2/m**3
   type (type_state_variable_id)   :: id_sDIMW,id_sDPOMW,id_sNPOMW,id_sPPOMW,id_sSiPOMW
   type (type_state_variable_id)   :: id_sPO4W,id_sPAIMW,id_sNH4W,id_sNO3W
   type (type_state_variable_id)   :: id_sSiO2W,id_sO2W
!  dissolved organic matter
   type (type_state_variable_id)   :: id_sDDOMW,id_sNDOMW,id_sPDOMW,id_sSiDOMW 
!  diagnostic variables for local output
!  rPDPOMW: P/D ratio of particulate organic matter
!  rNDPOMW: N/D ratio of particulate organic matter
!  rPDDOMW: P/D ratio of dissolved organic matter
!  rNDDOMW: N/D ratio of dissolved organic matter
!  wO2AbioW: abiotic water column oxygen consumption
!  tO2Aer: O2 reaeration rate at the water surface
   type (type_diagnostic_variable_id)  :: id_rPDPOMW,id_rNDPOMW
   type (type_diagnostic_variable_id)  :: id_rPDDOMW,id_rNDDOMW
   type (type_diagnostic_variable_id)           :: id_extIM,id_extPOM
   type (type_horizontal_diagnostic_variable_id):: id_wind,id_tO2Aer
#ifdef _DEVELOPMENT_
!  diagnostic variables for the total fluxes for state variables from abiotic_module
   type (type_diagnostic_variable_id)  :: id_wSiAbioDOMW,id_wDAbioDOMW,id_wPAbioDOMW,id_wNAbioDOMW
   type (type_diagnostic_variable_id)  :: id_wDAbioIMW,id_wPAbioPO4W,id_wPAbioAIMW,id_wNAbioNH4W
   type (type_diagnostic_variable_id)  :: id_wNAbioNO3W,id_wSiAbioSiO2W,id_wO2AbioW
#endif
!  environmental dependencies
   type (type_dependency_id)                :: id_uTm
   type (type_horizontal_dependency_id)     :: id_uVWind
!  Model parameters
!! mineralization process, split into 2 steps
!  from POM to DOM
   real(rk)                   :: cThetaMinPOMW, kDMinPOMW,kNMinPOMW
   real(rk)                   :: kPMinPOMw,kSiMinPOMW
!  from DOM to dissolved nutrients
   real(rk)                   :: cThetaMinDOMW,kDMinDOMW,kNMinDOMW
   real(rk)                   :: kPMinDOMW,kSiMinDOMW
!! nitrification pars
   real(rk)                   :: kNitrW,cThetaNitr,hO2Nitr
   real(rk)                   :: cThetaAer,cCPerDW,hO2BOD,O2PerNH4
   real(rk)                   :: NO3PerC,hNO3Denit,kPSorp,cRelPAdsD
   real(rk)                   :: cRelPAdsFe,fFeDIM,cRelPAdsAl,fAlDIM
   real(rk)                   :: fRedMax,cKPAdsOx
!  sinking parameter
   real(rk)                   :: cVSetIM,cVSetPOM
!  parameter for specific light attenuation coefficient
   real(rk)                   :: cExtSpIM,cExtSpPOM
!  minimum state variable values
   real(rk)   :: cDVegMin, cNVegMin,cPVegMin

   contains
!  Model procedures
   procedure :: initialize
   procedure :: do
   procedure :: get_light_extinction
   procedure :: do_surface

   end type type_pclake_abiotic_water

!  private data members
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
!  ratio of mol.weights, = 32/12 [gO2/gC],
   real(rk),parameter :: molO2molC = 2.6667_rk
!  ratio of mol.weights,32/14 [gO2/gN],
   real(rk),parameter :: molO2molN = 2.2857_rk
!  ratio of mol.weights,14/12 [gN/gC],
   real(rk),parameter :: molNmolC = 1.1667_rk
!
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the inorganic and organic matter model in water
!
! !INTERFACE:

   subroutine initialize(self,configunit)
!

! !INPUT PARAMETERS:
   class (type_pclake_abiotic_water), intent(inout), target  :: self
   integer,                          intent(in)            :: configunit

!EOP
!-----------------------------------------------------------------------------
!BOC

!  Store parameter values in derived type
!  NB: all rates must be provided in values per day in .yaml input,
!  and are converted here to values per second.
   call self%get_parameter(self%cCPerDW,        'cCPerDW',      'gC/gDW', 'C content of organic matter',                                default=0.4_rk)
   call self%get_parameter(self%cExtSpPOM,      'cExtSpPOM',    'm2/gDW', 'specific extinction factor for POM',                         default=0.15_rk)
   call self%get_parameter(self%cExtSpIM,       'cExtSpIM',     'm2/gDW', 'specific extinction inorganic matter',                       default=0.05_rk)
   call self%get_parameter(self%cKPAdsOx,       'cKPAdsOx',     'm3/gP',  'P adsorption affinity at oxidized conditions',               default=0.6_rk)
   call self%get_parameter(self%cRelPAdsAl,     'cRelPAdsAl',   'gP/gAl', 'maximum P adsorption per g Al',                              default=0.134_rk)
   call self%get_parameter(self%cRelPAdsD,      'cRelPAdsD',    'gP/gD',  'maximum P adsorption per g DW',                              default=0.00003_rk)
   call self%get_parameter(self%cRelPAdsFe,     'cRelPAdsFe',   'gP/gFe', 'maximum P adsorption per g Fe',                              default=0.065_rk)
   call self%get_parameter(self%cThetaAer,      'cThetaAer',    '1/e^oC', 'temperature coefficient for reaeration',                     default=1.024_rk)
   call self%get_parameter(self%cThetaNitr,     'cThetaNitr',   '[-]',    'temperature coefficient for nitrification',                  default=1.08_rk)
   call self%get_parameter(self%cVSetPOM,       'cVSetPOM',     'm/day',  'maximum settling rate of POM',                               default=-0.25_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cVSetIM,        'cVSetIM',      'm/day',  'maximum settling rate of iorganic matter',                   default=-1.0_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%fAlDIM,         'fAlDIM',       'gAl/gD', 'Al content of inorganic matter',                             default=0.01_rk)
   call self%get_parameter(self%fFeDIM,         'fFeDIM',       'gFe/gD', 'Fe content of inorganic matter',                             default=0.01_rk)
   call self%get_parameter(self%fRedMax,        'fRedMax',      '[-]',    'maximum reduction factor of P adsorption affinity',          default=0.9_rk)
   call self%get_parameter(self%hNO3Denit,      'hNO3Denit',    'mgN/l',  'quadratic half saturation NO3 conc. for denitrification',    default=2.0_rk)
   call self%get_parameter(self%hO2BOD,         'hO2BOD',       'mgO2/l', 'half saturation oxygen conc. for BOD',                       default=1.0_rk)
   call self%get_parameter(self%hO2Nitr,        'hO2Nitr',      'mgO2/l', 'half saturation oxygen conc. for nitrification',             default=2.0_rk)
   call self%get_parameter(self%kNitrW,         'kNitrW',       'day-1',  'nitrification rate constant in water',                       default=0.1_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kPSorp,         'kPSorp',       'day-1',  'P adsorption rate constant not too high -> model speed day', default=0.05_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%NO3PerC,        'NO3PerC',      'molNO3', 'denitrified NO3 per mol C mineralised',                      default=0.8_rk)
   call self%get_parameter(self%O2PerNH4,       'O2PerNH4',     'molO2',  'used O2 per mol NH4+ nitrified',                             default=2.0_rk)
!  mineralization parameters
!  from POM to DOM
   call self%get_parameter(self%cThetaMinPOMW,  'cThetaMinPOMW','[-]',    'temperature coefficient for mineralization from POM to DOM', default=1.07_rk)
   call self%get_parameter(self%kDMinPOMW,      'kDMinPOMW',    'day-1',  'decomposition constant for POM-DW to DOM-DW',                default=0.01_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kNMinPOMW,      'kNMinPOMW',    'day-1',  'decomposition constant for POM-N to DOM-N',                  default=0.01_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kPMinPOMW,      'kPMinPOMW',    'day-1',  'decomposition constant for POM-P to DOM-P',                  default=0.01_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kSiMinPOMW,     'kSiMinPOMW',   'day-1',  'decomposition constant for POM-Si to DOM-Si',                default=0.01_rk,scale_factor=1.0_rk/secs_pr_day)
!  from DOM to dissolved nutrients
   call self%get_parameter(self%cThetaMinDOMW,  'cThetaMinDOMW',    '[-]', 'temperature coefficient for DOM mineralization',            default=1.07_rk)
   call self%get_parameter(self%kDMinDOMW,      'kDMinDOMW',     'day-1',  'mineralization constant of dissolved organic DW',           default=0.01_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kNMinDOMW,      'kNMinDOMW',     'day-1',  'mineralization constant of dissolved organic N',            default=0.01_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kPMinDOMW,      'kPMinDOMW',     'day-1',  'mineralization constant of dissolved organic P',            default=0.01_rk,scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kSiMinDOMW,     'kSiMinDOMW',    'day-1',  'mineralization constant of dissolved organic Si',           default=0.01_rk,scale_factor=1.0_rk/secs_pr_day)
!  Register local state variable
!  particles, including inorganic matter(sDIM) and organic matter(sDPOMW,sNPOMW,sPPOMW,sSiPOMW) have
!  vertical movement, usually settling(negative values)
   call self%register_state_variable(self%id_sDIMW,'sDIMW','g m-3','inorganic matter',           &
                                    initial_value=5.0_rk,  minimum=_ZERO_, vertical_movement= self%cVSetIM,no_river_dilution=.TRUE.)
!  particulate organic matter
   call self%register_state_variable(self%id_sDPOMW,'sDPOMW','g m-3','particulate organic matter DW',    &
                                    initial_value=1.0_rk,  minimum=_ZERO_, vertical_movement= self%cVSetPOM,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPPOMW,'sPPOMW','g m-3','particulate organic P',     &
                                    initial_value=0.005_rk,minimum=_ZERO_,vertical_movement= self%cVSetPOM,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNPOMW,'sNPOMW','g m-3','particulate organic N',      &
                                    initial_value=0.05_rk ,minimum=_ZERO_,vertical_movement=self%cVSetPOM,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sSiPOMW,'sSiPOMW','g m-3','particulate organic Si',      &
                                    initial_value=0.02_rk, minimum=_ZERO_,vertical_movement= self%cVSetPOM,no_river_dilution=.TRUE.)
!  dissolved organic matter
   call self%register_state_variable(self%id_sDDOMW,'sDDOMW','g m-3','dissolved organic matter DW',     &
                                    initial_value=1.0_rk,minimum=_ZERO_,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNDOMW,'sNDOMW','g m-3','dissolved organic N',     &
                                    initial_value=0.025_rk,minimum=_ZERO_,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPDOMW,'sPDOMW','g m-3','dissolved organic P',     &
                                    initial_value=0.0025_rk,minimum=_ZERO_,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sSiDOMW,'sSiDOMW','g m-3','dissolved organic Si',     &
                                    initial_value=0.0125_rk,minimum=_ZERO_,no_river_dilution=.TRUE.)
!  dissolved nutrients
   call self%register_state_variable(self%id_sPO4W,  'sPO4W',   'g m-3','phosphate',     &
                                    initial_value=0.01_rk, minimum=_ZERO_,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sPAIMW,'sPAIMW','g m-3','IM-adsorbed P',     &
                                    initial_value=0.0_rk,minimum=_ZERO_,vertical_movement= self%cVSetIM,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNH4W,'sNH4W','g m-3','ammonium',     &
                                    initial_value=0.1_rk,minimum=_ZERO_,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sNO3W,'sNO3W','g m-3','nitrate',     &
                                    initial_value=0.1_rk,minimum=_ZERO_,no_river_dilution=.TRUE.)
   call self%register_state_variable(self%id_sSiO2W,'sSiO2W','g m-3','silica dioxide',     &
                                    initial_value=3.0_rk,minimum=_ZERO_,no_river_dilution=.TRUE.)
!  oxygen
   call self%register_state_variable(self%id_sO2W,'sO2W','g m-3','oxygen',     &
                                    initial_value=10.0_rk,minimum=_ZERO_,no_river_dilution=.TRUE.)

!  Register diagnostic variables for dependencies in other modules
   call self%register_diagnostic_variable(self%id_rPDPOMW, 'rPDPOMW',   '[-]',       'particulate organic matter P/D ratio',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDPOMW, 'rNDPOMW',   '[-]',       'particulate organic matter N/D ratio',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rPDDOMW, 'rPDDOMW',   '[-]',       'dissolved organic matter P/D ratio',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDDOMW, 'rNDDOMW',   '[-]',       'dissolved organic matter N/D ratio',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wind,    'wind',      'm/s',       'windspeed',                output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_extIM,   'extIM',     '[-]',       'extIM',                    output=output_instantaneous)!,source=source_get_light_extinction)
   call self%register_diagnostic_variable(self%id_extPOM,  'extPOM',    '[-]',       'extPOM',                   output=output_instantaneous)!,source=source_get_light_extinction)
   call self%register_diagnostic_variable(self%id_tO2Aer,  'tO2Aer',    'g m-2 s-1', 'O2_reareation',            output=output_instantaneous)
#ifdef _DEVELOPMENT_
!  Register diagnostic variables for modular fluxes for state variables
   call self%register_diagnostic_variable(self%id_wDAbioDOMW,  'wDAbioDOMW',   'g m-3 s-1', 'abiotic_water_DDOMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPAbioDOMW,  'wPAbioDOMW',   'g m-3 s-1', 'abiotic_water_PDOMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNAbioDOMW,  'wNAbioDOMW',   'g m-3 s-1', 'abiotic_water_NDOMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wSiAbioDOMW, 'wSiAbioDOMW',  'g m-3 s-1', 'abiotic_water_SiDOMW_change', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wDAbioIMW,   'wDAbioIMW',    'g m-3 s-1', 'abiotic_water_IMW_change',    output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPAbioPO4W,  'wPAbioPO4W',   'g m-3 s-1', 'abiotic_water_PO4W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wPAbioAIMW,  'wPAbioAIMW',   'g m-3 s-1', 'abiotic_water_PAIMW_change',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNAbioNH4W,  'wNAbioNH4W',   'g m-3 s-1', 'abiotic_water_NH4W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wNAbioNO3W,  'wNAbioNO3W',   'g m-3 s-1', 'abiotic_water_NO3W_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wSiAbioSiO2W,'wSiAbioSiO2W', 'g m-3 s-1', 'abiotic_water_SiO2_change',   output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_wO2AbioW,    'wO2AbioW',     'g m-3 s-1', 'abiotic_water_O2_change',     output=output_instantaneous)
#endif
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,              self%id_sNPOMW)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,              self%id_sNH4W)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,              self%id_sNO3W)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,            self%id_sPO4W)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,            self%id_sPAIMW)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,            self%id_sPPOMW)
   call self%add_to_aggregate_variable(standard_variables%total_silicate,              self%id_sSiO2W)
   call self%add_to_aggregate_variable(standard_variables%total_silicate,              self%id_sSiPOMW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totN'),self%id_sNPOMW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totN'),self%id_sNH4W)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totN'),self%id_sNO3W)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totP'),self%id_sPO4W)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totP'),self%id_sPAIMW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totP'),self%id_sPPOMW)
!  add dissolved organic N and P to PCLake totN and totP
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totN'),self%id_sNDOMW)
   call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pclake_totP'),self%id_sPDOMW)

!  register environmental dependencies
   call self%register_dependency(self%id_uTm,   standard_variables%temperature)
   call self%register_dependency(self%id_uVWind,standard_variables%wind_speed)

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
! !INPUT PARAMETERS:
   class (type_pclake_abiotic_water), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_
! !LOCAL VARIABLES:
!  carriers for local state variables
   real(rk)     :: sNH4W,sNO3W,sPO4W,sPAIMW,sSiO2W,sO2W
!  Nutrients ratio value
   real(rk)     :: rPDPOMW,rNDPOMW
   real(rk)     :: rPDDOMW,rNDDOMW
!  carriers for environmental dependencies
   real(rk)     :: uTm
!  local variables for processes
!  Temperature function variables
   real(rk)     :: uFunTmNitr
!  O2 functions variables
   real(rk)     :: aCorO2BOD,aCorO2NitrW,wO2MinDOMW,wO2NitrW
!  adsorption function variables
   real(rk)     :: wPSorpIMW,aPEqIMW,aPIsoAdsW,aPAdsMaxW,aKPAdsW
!  nitrification functions variables
   real(rk)     :: wNNitrW
!  denitrification functions variables
   real(rk)     :: wDDenitW,wNDenitW
!  total fluxes
   real(rk)     :: wNAbioNH4W,wNAbioNO3W,wPAbioPO4W,wPAbioAIMW
   real(rk)     :: wDAbioIMW,wSiAbioSiO2W,wO2AbioW
!-----------------------------------------------------------------------
!  Variable section for organic matter (POM and DOM)
!-----------------------------------------------------------------------
!  state variable carriers
   real(rk)     :: sDPOMW,sNPOMW,sPPOMW,sDIMW,sSiPOMW
   real(rk)     :: sDDOMW, sNDOMW,sPDOMW,sSiDOMW
 !  minerilization step 1, from particulate to dissolved
   real(rk)     :: kNMinPOMW,kPMinPOMW,kSiMinPOMW,uFunTmMinPOMW
   real(rk)     :: wDMinPOMW,wNMinPOMW,wPMinPOMW,wSiMinPOMW
   real(rk)     :: wDAbioPOMW,wNAbioPOMW,wPAbioPOMW,wSiAbioPOMW
!  O2 comsuption during first step
   real(rk)     :: wO2MinPOMW
!  mineralization function variables
   real(rk)     :: wDMinDOMW,wNMinDOMW,wPMinDOMW,wSiMinDOMW,uFunTmMinDOMW  
   real(rk)     :: wDAbioDOMW,wNAbioDOMW,wPAbioDOMW,wSiAbioDOMW
!
!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
   _LOOP_BEGIN_
!  Retrieve current (local) state variable values.
   _GET_(self%id_sDIMW,sDIMW)
   _GET_(self%id_sDPOMW,sDPOMW)
   _GET_(self%id_sPPOMW,sPPOMW)
   _GET_(self%id_sNPOMW,sNPOMW)
   _GET_(self%id_sSiPOMW,sSiPOMW)
   _GET_(self%id_sPO4W,sPO4W)
   _GET_(self%id_sPAIMW,sPAIMW)
   _GET_(self%id_sNH4W,sNH4W)
   _GET_(self%id_sNO3W,sNO3W)
   _GET_(self%id_sO2W,sO2W)
   _GET_(self%id_sSiO2W,sSiO2W)
   _GET_(self%id_sDDOMW,sDDOMW)
   _GET_(self%id_sNDOMW,sNDOMW)
   _GET_(self%id_sPDOMW,sPDOMW)
   _GET_(self%id_sSiDOMW,sSiDOMW)

!  retrieve current environmental dependencies
   _GET_(self%id_uTm,uTm)
!  Nutrients ratio of particulate organic matter
   rPDPOMW=sPPOMW/(sDPOMW+NearZero)
   rNDPOMW= sNPOMW/(sDPOMW+NearZero)
!  Nutrients ratio of dissolved organic matter
   rPDDOMW=sPDOMW/(sDDOMW+NearZero)
   rNDDOMW= sNDOMW/(sDDOMW+NearZero)
!-----------------------------------------------------------------------
!  Temperature functions
!-----------------------------------------------------------------------
!  Temperature_dependence_for_nitrification
   uFunTmNitr = uFunTmAbio(uTm,self%cThetaNitr)
!-----------------------------------------------------------------------
!  mineralization step 1, from POM to DOM
!-----------------------------------------------------------------------
!  temp._function_of_mineralization_in_water for first step min.
   uFunTmMinPOMW = uFunTmAbio(uTm,self%cThetaMinPOMW)
!  decomposition
   wDMinPOMW = self%kDMinPOMW * uFunTmMinPOMW * sDPOMW
!  particulate_OM_P mineralization
   wPMinPOMW = self%kPMinPOMW * uFunTmMinPOMW * sPPOMW
!  particulate_OM_N mineralization
   wNMinPOMW = self%kNMinPOMW * uFunTmMinPOMW * sNPOMW
!  particulate_OM_Si mineralization
   wSiMinPOMW = self%kSiMinPOMW * uFunTmMinPOMW * sSiPOMW
!  correction_of_O2_demand_in_water_at_low_oxygen_conc.
   aCorO2BOD = sO2W / (self%hO2BOD + sO2W)
!  O2 comsuption for step POM to DOM
   wO2MinPOMW = molO2molC * self%cCPerDW * aCorO2BOD * wDMinPOMW
!-----------------------------------------------------------------------
!  mineralization step 2, from DOM to inorganic nutrients
!-----------------------------------------------------------------------
!  temperature_function_of_mineralization_in_water
   uFunTmMinDOMW = uFunTmAbio(uTm,self%cThetaMinDOMW)
!  decomposition
   wDMinDOMW = self%kDMinDOMW * uFunTmMinDOMW * sDDOMW
!  dissolved_OM_P mineralization
   wPMinDOMW = self%kPMinDOMW * uFunTmMinDOMW * sPDOMW
!  dissolved_OM mineralization
   wNMinDOMW = self%kNMinDOMW * uFunTmMinDOMW * sNDOMW
!  dissolved_OM_Si mineralization
   wSiMinDOMW = self%kSiMinDOMW * uFunTmMinDOMW * sSiDOMW
!-----------------------------------------------------------------------
!  Phosphorus adsorption
!-----------------------------------------------------------------------
!  max._P_adsorption_per_g_inorg._matter_in_water
   aPAdsMaxW = self%cRelPAdsD+aCorO2BOD*self%cRelPAdsFe*self%fFeDIM+ &
   & self%cRelPAdsAl*self%fAlDIM
!  P_adsorption_affinity,_corrected_for_redox_conditions
   aKPAdsW = (1.0_rk - self%fRedMax * (1.0_rk-aCorO2BOD)) * self%cKPAdsOx
!  P_adsorption_isotherm_onto_inorg._matter_in_sediment
   aPIsoAdsW = aPAdsMaxW * aKPAdsW * sPO4W / (1.0_rk + aKPAdsW * sPO4W)
!  equilibrium_conc.
   aPEqIMW = aPIsoAdsW * sDIMW
!  adsorption_flux_in_water
   wPSorpIMW = self%kPSorp * (aPEqIMW - sPAIMW)
!-----------------------------------------------------------------------
!  Nitrification functions
!-----------------------------------------------------------------------
   aCorO2NitrW = sO2W*sO2W / (self%hO2Nitr*self%hO2Nitr + sO2W*sO2W)
!  nitrification_flux
   wNNitrW = self%kNitrW * uFunTmNitr * aCorO2NitrW * sNH4W
!-----------------------------------------------------------------------
!  Denitrification functions
!-----------------------------------------------------------------------
!  mineralization_flux_by_denitrification
   wDDenitW = sNO3W*sNO3W/(self%hNO3Denit*self%hNO3Denit+sNO3W*sNO3W)* &
   & (1.0_rk-aCorO2BOD)*(wDMinDOMW+wDMinPOMW)
!  Denitrification_flux
   wNDenitW=self%NO3PerC*molNmolC*self%cCPerDW*wDDenitW
!-----------------------------------------------------------------------
!  O2 dynamics
!-----------------------------------------------------------------------
!  O2_flux_due_to_mineralization_of_organic matter
   wO2MinDOMW = molO2molC * self%cCPerDW * aCorO2BOD * wDMinDOMW
!  O2_flux_due_to_nitrification
   wO2NitrW = self%O2PerNH4 * molO2molN * wNNitrW
!-----------------------------------------------------------------------
!  total abiotic flux for each state variable in water
!  the abiot auxilaries have different meanings of the original one
!-----------------------------------------------------------------------
!  total_abiotic/microbial_DW_inorganic_matter_flux_in_water
   wDAbioIMW=0.0_rk
!  change of particulate organics
!  total_abiotic/microbial_DW particulate OM flux in water
   wDAbioPOMW=-wDMinPOMW
!  total_abiotic/microbial_N particulate OM flux in water
   wNAbioPOMW =-wNMinPOMW
!  total_abiotic/microbial_P particulate OM flux in water
   wPAbioPOMW=-wPMinPOMW
!  total_abiotic/microbial_Si particulate OM flux in water
   wSiAbioPOMW=-wSiMinPOMW
!  change of dissolved detritus
!  total_abiotic/microbial_DW dissoved orignic flux in water
   wDAbioDOMW=wDMinPOMW - wDMinDOMW
!  total_abiotic/microbial_N dissolved organic flux in water
   wNAbioDOMW = wNMinPOMW - wNMinDOMW
!  total_abiotic/microbial_P dissolved organic flux in water
   wPAbioDOMW= wPMinPOMW - wPMinDOMW
!  total_abiotic/microbial_Si dissolved organic flux in water
   wSiAbioDOMW= wSiMinPOMW - wSiMinDOMW
!  total_abiotic/microbial_PO4 flux in water
   wPAbioPO4W=wPMinDOMW-wPSorpIMW
!  total_abiotic/microbial_P absorbed on inorganic matter flux in water
   wPAbioAIMW=wPSorpIMW
!  total_abiotic/microbial_NH4 flux in water
   wNAbioNH4W = wNMinDOMW-wNNitrW
!  total_abiotic/microbial_NO3 flux in water
   wNAbioNO3W = wNNitrW - wNDenitW
!  total_abiotic/microbial_SiO2 flux in water
   wSiAbioSiO2W=wSiMinDOMW
!  total_abiotic/microbial_O2 flux in water
   wO2AbioW = - wO2MinPOMW - wO2MinDOMW - wO2NitrW
!-----------------------------------------------------------------------
!  Update local state variables
!-----------------------------------------------------------------------
!  Set temporal derivatives
   _SET_ODE_(self%id_sDIMW,wDAbioIMW)
   _SET_ODE_(self%id_sPO4W,wPAbioPO4W)
   _SET_ODE_(self%id_sPAIMW,wPAbioAIMW)
   _SET_ODE_(self%id_sNH4W,wNAbioNH4W)
   _SET_ODE_(self%id_sNO3W,wNAbioNO3W)
   _SET_ODE_(self%id_sSiO2W,wSiAbioSiO2W)
   _SET_ODE_(self%id_sO2W,wO2AbioW)
!  update orgainc variables
!  particulate organic matter
   _SET_ODE_(self%id_sDPOMW,wDAbioPOMW)
   _SET_ODE_(self%id_sPPOMW,wPAbioPOMW)
   _SET_ODE_(self%id_sNPOMW,wNAbioPOMW)
   _SET_ODE_(self%id_sSiPOMW,wSiAbioPOMW)
!! Dissolved organic matter
   _SET_ODE_(self%id_sDDOMW,wDAbioDOMW)
   _SET_ODE_(self%id_sPDOMW,wPAbioDOMW)
   _SET_ODE_(self%id_sNDOMW,wNAbioDOMW)
   _SET_ODE_(self%id_sSiDOMW,wSiAbioDOMW)
!-----------------------------------------------------------------------
!  Output local diagnostic variables
!-----------------------------------------------------------------------
   _SET_DIAGNOSTIC_(self%id_rPDPOMW,rPDPOMW)
   _SET_DIAGNOSTIC_(self%id_rNDPOMW,rNDPOMW)
   _SET_DIAGNOSTIC_(self%id_rPDDOMW,rPDDOMW)
   _SET_DIAGNOSTIC_(self%id_rNDDOMW,rNDDOMW)
#ifdef _DEVELOPMENT_
!  output diagnostic variables for modular fluxes
   _SET_DIAGNOSTIC_(self%id_wDAbioDOMW,   wDAbioDOMW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wPAbioDOMW,   wPAbioDOMW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wNAbioDOMW,   wNAbioDOMW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wSiAbioDOMW,  wSiAbioDOMW*secs_pr_day)

   _SET_DIAGNOSTIC_(self%id_wDAbioIMW,    wDAbioIMW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wPAbioPO4W,   wPAbioPO4W*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wPAbioAIMW,   wPAbioAIMW*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wNAbioNH4W,   wNAbioNH4W*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wNAbioNO3W,   wNAbioNO3W*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wSiAbioSiO2W, wSiAbioSiO2W*secs_pr_day)
   _SET_DIAGNOSTIC_(self%id_wO2AbioW,     wO2AbioW*secs_pr_day)
#endif

   _LOOP_END_
!-----------------------------------------------------------------------
!  Spatial loop end
!-----------------------------------------------------------------------
   end subroutine do

!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: 
!
! !INTERFACE:
   subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!
! !INPUT PARAMETERS:
   class (type_pclake_abiotic_water), intent(in) :: self
   _DECLARE_ARGUMENTS_GET_EXTINCTION_
!
! !REVISION HISTORY:
!  Original author(s):Fenjuan Hu
! !LOCAL VARIABLES:
   real(rk) :: sDIMW,sDPOMW
   real(rk) :: extIM,extPOM
!
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Enter spatial loops (if any)
   _LOOP_BEGIN_

   ! Retrieve current (local) state variable values.
   _GET_(self%id_sDIMW,sDIMW)
   _GET_(self%id_sDPOMW,sDPOMW)

   extIM=self%cExtSpIM*sDIMW
   extPOM=self%cExtSpPOM*sDPOMW

   ! Self-shading with explicit contribution from background phytoplankton concentration.
   _SET_EXTINCTION_(extIM+extPOM)

   _SET_DIAGNOSTIC_(self%id_extIM,extIM)
   _SET_DIAGNOSTIC_(self%id_extPOM,extPOM)

   ! Leave spatial loops (if any)
   _LOOP_END_

   end subroutine get_light_extinction
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Surface fluxes for Oxygen reaeration from the air
!
! !INTERFACE:
   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   class (type_pclake_abiotic_water),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_
!  local variables
   real(rk)                   :: uVWind,sO2W,uTm
   real(rk)                   :: tO2Aer,uO2Sat,aFunLemnAer,kAer,uFunTmAer
!EOP
!-----------------------------------------------------------------------
!BOC
   _HORIZONTAL_LOOP_BEGIN_
!  Retrieve environmental dependencies
   _GET_HORIZONTAL_(self%id_uVWind,uVWind)
   _GET_(self%id_uTm,uTm)
!  Retrieve state variable values
   _GET_(self%id_sO2W,sO2W)
!-----------------------------------------------------------------------
!     oxygen_reaeration functions
!-----------------------------------------------------------------------
!  temperature_function_of_reaeration
   uFunTmAer =  uFunTmAbio(uTm,self%cThetaAer)
!  oxygen_saturation_concentration
   uO2Sat = 14.652_rk - 0.41022 * uTm + 0.007991_rk * uTm*uTm - 0.000077774_rk * uTm*uTm*uTm
!  reaeration_coefficient
   kAer=(0.727_rk*((uVWind)**(0.5_rk))-0.371_rk*uVWind+0.0376_rk*uVWind*uVWind)
!  duckweed_function_of_reaeration
   aFunLemnAer = 1.0_rk
!  reaeration_flux_of_O2_into_the_water
   tO2Aer = kAer * uFunTmAer * (uO2Sat - sO2W) * aFunLemnAer
!-----------------------------------------------------------------------
!  update oxygen in reaeration process
!-----------------------------------------------------------------------
!  for PCLake Benchmark
!  convert daily rates to seconds
!  keep this, time unit is in day
   _SET_SURFACE_EXCHANGE_(self%id_sO2W,tO2Aer/secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tO2Aer,tO2Aer)
!   gotm output for 0d input
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_wind,uVWind)
   _HORIZONTAL_LOOP_END_

   end subroutine do_surface
!EOC

!-----------------------------------------------------------------------
   end module pclake_abiotic_water
!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
