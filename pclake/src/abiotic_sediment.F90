#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module pclake_abiotic_sediment

! !USES:
   use fabm_types
   use pclake_utility, ONLY: uFunTmAbio
   implicit none
!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_pclake_abiotic_sediment
!  local state variable identifers
!  sDIMS: inorganic matter concentration, in dry-weight, gDW/m**2
!  sDPOMS, particulate organic matter, in dry-weight, gDW/m**2
!  sNPOMS, particulate organic matter, in nitrogen element, gN/m**2
!  sPPOMS, particulate organic matter, in phosphorus element, gP/m**2
!  sSiPOMS, particulate organic matter, in silicon element, gSi/m**2
!  sDDOMS, dissolved organic matter, in dry-weight, gDW/m**2
!  sNDOMS, dissolved organic matter, in nitrogen element, gN/m**2
!  sPDOMS, dissolved organic matter, in phosphorus element, gP/m**2
!  sSiDOMS, dissolved organic matter, in silicon element, gSi/m**2
!  sNH4S, ammonium concentration, in nitrogen element, gN/m**2
!  sNO3S, nitrate concentration, in nitrogen element, gN/m**2
!  sPO4S, phosphate concentration, in phosphorus element, gP/m**2
!  sPAIMS, absorbed phosphorus concentration, in phosphorus element, gP/m**2
   type (type_bottom_state_variable_id) :: id_sDIMS,id_sDPOMS,id_sNPOMS
   type (type_bottom_state_variable_id) :: id_sPPOMS,id_sSiPOMS
   type (type_bottom_state_variable_id) :: id_sDDOMS,id_sNDOMS
   type (type_bottom_state_variable_id) :: id_sPDOMS,id_sSiDOMS
   type (type_bottom_state_variable_id) :: id_sPO4S,id_sPAIMS,id_sNH4S
   type (type_bottom_state_variable_id) :: id_sNO3S
   type (type_bottom_state_variable_id) :: id_sDHumS,id_sNHumS,id_sPHumS
!  diagnostic variables for local output
!  rPDPOMS: P/D ratio of particulate organic matter
!  rNDPOMS: N/D ratio of particulate organic matter
!  rPDDOMS: P/D ratio of dissolved organic matter
!  rNDDOMS: N/D ratio of dissolved organic matter
   type (type_horizontal_diagnostic_variable_id) :: id_rPDPOMS,id_rNDPOMS
   type (type_horizontal_diagnostic_variable_id) :: id_rPDDOMS,id_rNDDOMS
!  aPEqIMS: equilibrium absorbed phosphorus concentration
   type (type_horizontal_diagnostic_variable_id) :: id_aPEqIMS,id_afOxySed
!  diagnostic variable for external dependencies: 
   type (type_horizontal_diagnostic_variable_id) :: id_tDAbioHumS,id_tDAbioPOMS,id_tDAbioO2S
   type (type_horizontal_diagnostic_variable_id) :: id_tNdifNH4,id_tNdifNO3,id_tPdifPO4
#ifdef _DEVELOPMENT_
!  diagnostic variables for modular fluxes for each module
   type (type_horizontal_diagnostic_variable_id) :: id_tDAbioIMS,id_tPAbioPOMS
   type (type_horizontal_diagnostic_variable_id) :: id_tNAbioPOMS,id_tSiAbioPOMS,id_tNAbioNH4S
   type (type_horizontal_diagnostic_variable_id) :: id_tNAbioNO3S,id_tPAbioPO4S,id_tPAbioAIMS
   type (type_horizontal_diagnostic_variable_id) :: id_tPAbioHumS,id_tNAbioHumS
   type (type_horizontal_diagnostic_variable_id) :: id_tSiAbioSiO2S
!  particulate organic matter and humus fluxes will be named differently as they are used by 
!  external dependencies
   type (type_horizontal_diagnostic_variable_id) :: id_tDAbioHumSflux,id_tDAbioPOMSflux
#endif
!  state dependencies identifiers
!  MinSiO2Sed: mineralization generated SiO2 from sediment
!  O2ConsumpSed: O2 consumption in sediment
!  diff+nut: diffusion fluxes of nutrients between water and sediment
   type (type_state_variable_id) :: id_MinSiO2Sed ,id_O2ConsumpSed
   type (type_state_variable_id) :: id_diffNH4,id_diffNO3,id_diffPO4
!  added dissolved organic matter fluxes dependencies
   type (type_state_variable_id) :: id_diffDDOM,id_diffNDOM
   type (type_state_variable_id) :: id_diffPDOM,id_diffSiDOM
!  environmental dependencies
   type (type_dependency_id)                :: id_uTm,id_dz
   type (type_horizontal_dependency_id)     :: id_depth
!  Model parameters
   real(rk)                   :: cDepthS,cCPerDW,O2PerNH4
!  sediment properties parameters
   real(rk)                   :: bPorS,bPorCorS
!  P-sorption parameters
   real(rk)                   :: kPSorp,cRelPAdsD
   real(rk)                   :: cRelPAdsFe,fFeDIM,cRelPAdsAl,fAlDIM
   real(rk)                   :: fRedMax,cKPAdsOx,kPChemPO4,coPO4Max
!  denitrification parameters
   real(rk)                   :: NO3PerC,hNO3Denit
   real(rk)                   :: fRefrPOMS
   real(rk)                   :: kNitrS,cThetaNitr
!  humus related parameters
   real(rk)                   :: kDMinHum
!  diffusion parameters
   real(rk)                   :: fDepthDifS,cThetaDif,cTurbDifNut
   real(rk)                   :: kNDifNH4,kNDifNO3,kPDifPO4
   real(rk)                   :: kO2Dif,cTurbDifO2
!  mineralization pars from POM to DOM
   real(rk)                   :: cThetaMinPOMS,kDMinPOMS,kNMinPOMS
   real(rk)                   :: kPMinPOMS,kSiMinPOMS
!  mineralization pars from DOM to dissolved nutrients
   real(rk)                   :: cThetaMinDOMS,kDMinDOMS,kNMinDOMS
   real(rk)                   :: kPMinDOMS,kSiMinDOMS
!  added dissolved organic diffusion pars
   real(rk)                   :: kDDifDOM,kNDifDOM,kPDifDOM,kSiDifDOM
!
   contains
!
   procedure initialize
   procedure do_bottom
   end type type_pclake_abiotic_sediment

!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day= 86400.0_rk
   real(rk),parameter :: NearZero=0.000000000000000000000000000000001_rk
!  ratio of mol.weights,=32/12 [gO2/gC],
   real(rk),parameter :: molO2molC=2.6667_rk
!  ratio of mol.weights,32/14 [gO2/gN],
   real(rk),parameter :: molO2molN=2.2857_rk
!  ratio of mol.weights,14/12 [gN/gC],
   real(rk),parameter :: molNmolC=1.1667_rk
!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the pclake_abiotic_sediment module
!
! !INTERFACE:
   subroutine initialize(self,configunit)
! !INPUT PARAMETERS:
   class (type_pclake_abiotic_sediment), intent(inout), target :: self
   integer,                              intent(in)            :: configunit
!  Store parameter values in derived type
!  NB: all rates must be provided in values per day in .yaml input,
!  and are converted here to values per second.
   call self%get_parameter(self%cDepthS,    'cDepthS',     'm',                    'sediment depth',                                           default=0.1_rk)
   call self%get_parameter(self%fRefrPOMS,  'fRefrPOMS',   '[-]',                  'refractory fraction of sediment POM',                      default=0.15_rk)
   call self%get_parameter(self%cCPerDW,    'cCPerDW',     'gC/gDW',               'C content of organic matter',                              default=0.4_rk)
   call self%get_parameter(self%O2PerNH4,   'O2PerNH4',    'mol',                  'O2 used per mol NH4+ nitrified',                           default=2.0_rk)
   call self%get_parameter(self%kNitrS,     'kNitrS',      'd-1',                  'nitrification rate constant',                              default=1.0_rk,    scale_factor =1.0_rk/secs_pr_day)
   call self%get_parameter(self%cThetaNitr, 'cThetaNitr',  '[-]',                  'temperature coefficient for nitrification',                default=1.08_rk)
   call self%get_parameter(self%NO3PerC,    'NO3PerC',     '[-]',                  'NO3 denitrified per mol C mineralised',                    default=0.8_rk)
   call self%get_parameter(self%hNO3Denit,  'hNO3Denit',   'mgN/l',                'quadratic half-sat. NO3 conc. for denitrification',        default=2.0_rk)
   call self%get_parameter(self%kPSorp,     'kPSorp',      'd-1',                  'P sorption rate constant not too high -> model speed',     default=0.05_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cRelPAdsD,  'cRelPAdsD',   'gP/gD',                'max. P adsorption per g DW',                               default=0.00003_rk)
   call self%get_parameter(self%cRelPAdsFe, 'cRelPAdsFe',  'gP/gFe',               'max. P adsorption per g Fe',                               default=0.065_rk)
   call self%get_parameter(self%fFeDIM,     'fFeDIM',      'gFe/gD',               'Fe content of inorg. matter',                              default=0.01_rk)
   call self%get_parameter(self%cRelPAdsAl, 'cRelPAdsAl',  'gP/gAl',               'max. P adsorption per g Al',                               default=0.134_rk)
   call self%get_parameter(self%fAlDIM,     'fAlDIM',      'gAl/gD',               'Al content of inorganic matter',                           default=0.01_rk)
   call self%get_parameter(self%fRedMax,    'fRedMax',     '[-]',                  'max. reduction factor of P adsorption affinity',           default=0.9_rk)
   call self%get_parameter(self%cKPAdsOx,   'cKPAdsOx',    'm3/gP',                'P adsorption affinity at oxidized conditions',             default=0.6_rk)
   call self%get_parameter(self%kPChemPO4,  'kPChemPO4',   'd-1',                  'chem. PO4 loss rate',                                      default=0.03_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%coPO4Max,   'coPO4Max',    'mgP/l',                'max. SRP conc. in pore water',                             default=1.0_rk)
   call self%get_parameter(self%bPorS,      'bPorS',       'm3 water m-3 sediment','sediment porosity',                                        default=0.847947_rk)
   call self%get_parameter(self%cThetaDif,  'cThetaDif',   '[-]',                  'temperature coefficient for diffusion',                    default=1.02_rk)
   call self%get_parameter(self%fDepthDifS, 'fDepthDifS',  '[-]',                  'nutrient diffusion distance as fraction of sediment depth', default=0.5_rk)
   call self%get_parameter(self%kNDifNH4,   'kNDifNH4',    'm2/day',               'molecular NH4 diffusion constant',                         default=0.000112_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cTurbDifNut,'cTurbDifNut', '[-]',                  'bioturbation factor for diffusion',                        default=5.0_rk)
   call self%get_parameter(self%bPorCorS,   'bPorCorS',    'm3 water m-3 sediment','sediment porosity, corrected for tortuosity',              default=0.737275_rk)
   call self%get_parameter(self%kNDifNO3,   'kNDifNO3',    'm2/day',               'molecular NO3 diffusion constant',                         default=0.000086_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kPDifPO4,   'kPDifPO4',    'm2/day',               'molecular PO4 diffusion constant',                         default=0.000072_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kO2Dif,     'kO2Dif',      'm2/day',               'molecular O2 diffusion constant',                          default=0.000026_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%cTurbDifO2, 'cTurbDifO2',  '[-]',                  'bioturbation factor for diffusion',                        default=5.0_rk)
   call self%get_parameter(self%kDMinHum,   'kDMinHum',    'd-1',                  'maximum decomposition constant of humic material (1D-5)',  default=0.00001_rk , scale_factor=1.0_rk/secs_pr_day)
!  organic matter mineralization parameters
!  POM to DOM
   call self%get_parameter(self%cThetaMinPOMS,'cThetaMinPOMS','[-]',           'temperature coeff. for sediment mineralization of POM to DOM',  default=1.07_rk)
   call self%get_parameter(self%kDMinPOMS,  'kDMinPOMS',   'd-1',              'mineralization constant in sediment from POM-DW to DOM-DW',     default=0.002_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kNMinPOMS,  'kNMinPOMS',   'd-1',              'mineralization constant in sediment from POM-N to DOM-N',       default=0.002_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kPMinPOMS,  'kPMinPOMS',   'd-1',              'mineralization constant in sediment from POM-P to DOM-P',       default=0.002_rk,  scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kSiMinPOMS, 'kSiMinPOMS',  'd-1',              'mineralization constant in sediment from POM-Si to DOM-Si',     default=0.002_rk,  scale_factor=1.0_rk/secs_pr_day)
!  added step 2. min. pars
   call self%get_parameter(self%cThetaMinDOMS, 'cThetaMinDOMS',  '[-]',        'exp. temperature constant of sediment mineralization', default=1.07_rk)
   call self%get_parameter(self%kDMinDOMS,  'kDMinDOMS',  'day-1',             'mineralization constant for sediment dissolved organic matter',  default=0.002_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kNMinDOMS,  'kNMinDOMS',  'day-1',             'mineralization constant for sediment dissolved organic N',       default=0.002_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kPMinDOMS,  'kPMinDOMS',  'day-1',             'mineralization constant for sediment dissolvedorganic P',        default=0.002_rk,    scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kSiMinDOMS, 'kSiMinDOMS', 'day-1',             'mineralization constant for sediment dissolved organic Si',      default=0.002_rk,    scale_factor=1.0_rk/secs_pr_day)
!  added diffusion pars
   call self%get_parameter(self%kDDifDOM,    'kDDifDOM',     'm2/day',             'molecular diffusion constant for dissolved organic matter', default=0.000112_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kNDifDOM,    'kNDifDOM',     'm2/day',             'molecular diffusion constant for dissolved organic N',     default=0.000112_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kPDifDOM,    'kPDifDOM',     'm2/day',             'molecular diffusion constant for dissolved organic P',     default=0.000112_rk, scale_factor=1.0_rk/secs_pr_day)
   call self%get_parameter(self%kSiDifDOM,   'kSiDifDOM',    'm2/day',             'molecular diffusion constant for dissolved organic Si',    default=0.000112_rk, scale_factor=1.0_rk/secs_pr_day)
!  Register local state variable
!  Inorganic matter
   call self%register_state_variable(self%id_sDIMS,  'sDIMS',  'g m-2','sediment inorganic matter',  initial_value=39611.3_rk, minimum=_ZERO_)
!  Particulate organic matter
   call self%register_state_variable(self%id_sDPOMS, 'sDPOMS', 'g m-2','sediment particulate organic DW',   initial_value=181.7_rk,   minimum=_ZERO_)
   call self%register_state_variable(self%id_sNPOMS, 'sNPOMS', 'g m-2','sediment particulate organic N',    initial_value=4.54_rk,    minimum=_ZERO_)
   call self%register_state_variable(self%id_sPPOMS, 'sPPOMS', 'g m-2','sediment particulate organic P',    initial_value=0.454_rk,   minimum=_ZERO_)
   call self%register_state_variable(self%id_sSiPOMS,'sSiPOMS','g m-2','sediment particulate organic Si',   initial_value=1.82_rk,    minimum=_ZERO_)
!  Dissolved organic matter
   call self%register_state_variable(self%id_sDDOMS, 'sDDOMS', 'g m-2','sediment dissolved organic matter',   initial_value=0.01_rk,   minimum=_ZERO_)
   call self%register_state_variable(self%id_sNDOMS, 'sNDOMS', 'g m-2','sediment dissolved organic N',  initial_value=0.01_rk,    minimum=_ZERO_)
   call self%register_state_variable(self%id_sPDOMS, 'sPDOMS', 'g m-2','sediment dissolved organic P',  initial_value=0.01_rk,   minimum=_ZERO_)
   call self%register_state_variable(self%id_sSiDOMS,'sSiDOMS','g m-2','sediment dissolved organic Si', initial_value=1.82_rk,    minimum=_ZERO_)
!  Dissolved nutrients
   call self%register_state_variable(self%id_sPO4S,  'sPO4S',   'g m-2','sediment phosphate',    initial_value=0.182_rk,   minimum=_ZERO_)
   call self%register_state_variable(self%id_sPAIMS, 'sPAIMS',  'g m-2','sediment absorbed phosphate',initial_value=17.99_rk,   minimum=_ZERO_)
   call self%register_state_variable(self%id_sNH4S,  'sNH4S',   'g m-2','sediment ammonium',       initial_value=0.02_rk,    minimum=_ZERO_)
   call self%register_state_variable(self%id_sNO3S,  'sNO3S',   'g m-2','sediment nitrate',     initial_value=0.002_rk,   minimum=_ZERO_)
!  Humus
   call self%register_state_variable(self%id_sDHumS,'sDHumS','g m-2','sediment humus DW',        initial_value=3452.34_rk, minimum=_ZERO_)
   call self%register_state_variable(self%id_sNHumS,'sNHumS','g m-2','sediment humus N',         initial_value=172.62_rk,  minimum=_ZERO_)
   call self%register_state_variable(self%id_sPHumS,'sPHumS','g m-2','sediment humus P',         initial_value=17.26_rk,   minimum=_ZERO_)
!  Register contribution of state to global aggregate variables
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNH4S)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNO3S)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNPOMS)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_sNDOMS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPO4S)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPAIMS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPPOMS)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_sPDOMS)
   call self%add_to_aggregate_variable(standard_variables%total_silicate,  self%id_sSiPOMS)
   call self%add_to_aggregate_variable(standard_variables%total_silicate,  self%id_sSiDOMS)
!---------------------------------------------------------------------------------------------------------------
!  register state variables dependencies
!---------------------------------------------------------------------------------------------------------------
!   Register dependencies on external state variables
   call self%register_state_dependency(self%id_O2ConsumpSed, 'oxygen_pool_water',               'g m-3', 'Oxygen pool water')
   call self%register_state_dependency(self%id_MinSiO2Sed,   'SiO2_generated_by_mineralization','g m-3', 'SiO2 generated by mineralization')
   call self%register_state_dependency(self%id_diffNH4,      'NH4_diffusion_flux',              'g m-3', 'NH4 diffusion flux')
   call self%register_state_dependency(self%id_diffNO3,      'NO3_diffusion_flux',              'g m-3', 'NO3 diffusion flux')
   call self%register_state_dependency(self%id_diffPO4,      'PO4_diffusion_flux',              'g m-3', 'PO4 diffusion flux')
   call self%register_state_dependency(self%id_diffDDOM,     'DDOM_diffusion_flux',             'g m-3', 'DDOM diffusion flux')
   call self%register_state_dependency(self%id_diffNDOM,     'NDOM_diffusion_flux',             'g m-3', 'NDOM diffusion flux')
   call self%register_state_dependency(self%id_diffPDOM,     'PDOM_diffusion_flux',             'g m-3', 'PDOM diffusion flux')
   call self%register_state_dependency(self%id_diffSiDOM,    'SiDOM_diffusion_flux',            'g m-3', 'SiDOM diffusion flux')
!  Register diagnostic variables
   call self%register_diagnostic_variable(self%id_afOxySed,       'afOxySed',       '[-]',       'fraction of aerobic sediment',        output=output_instantaneous,domain=domain_bottom)
   call self%register_diagnostic_variable(self%id_aPEqIMS,        'aPEqIMS',        '[-]',       'equilibrium adsorbed PO4',            output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rPDPOMS,        'rPDPOMS',        '[-]',       'particulate organic matter P/D_ratio',output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDPOMS,        'rNDPOMS',        '[-]',       'particulate organic matter N/D_ratio',output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rPDDOMS,        'rPDDOMS',        '[-]',       'dissolved organic matter P/D_ratio',  output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_rNDDOMS,        'rNDDOMS',        '[-]',       'dissolved organic matter N/D_ratio',  output=output_instantaneous)
#ifdef _DEVELOPMENT_
!  Register diagnostic variables for modular fluxes
!  Total fluxes to local state variables
   call self%register_diagnostic_variable(self%id_tDAbioIMS,      'tDAbioIMS',      'g m-2 s-1', 'abiotic_sediment_DIMS_change',        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAbioPOMSflux, 'tDAbioPOMSflux', 'g m-2 s-1', 'abiotic_sediment_DPOMS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAbioPOMS,     'tNAbioPOMS',     'g m-2 s-1', 'abiotic_sediment_NPOMS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAbioPOMS,     'tPAbioPOMS',     'g m-2 s-1', 'abiotic_sediment_PPOMS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tSiAbioPOMS,    'tSiAbioPOMS',    'g m-2 s-1', 'abiotic_sediment_SiPOMS_change',      output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAbioPO4S,     'tPAbioPO4S',     'g m-2 s-1', 'abiotic_sediment_PAIMS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAbioNH4S,     'tNAbioNH4S',     'g m-2 s-1', 'abiotic_sediment_NH4S_change',        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAbioNO3S,     'tNAbioNO3S',     'g m-2 s-1', 'abiotic_sediment_NO3S_change',        output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAbioAIMS,     'tPAbioAIMS',     'g m-2 s-1', 'abiotic_sediment_PAIMS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tDAbioHumSflux, 'tDAbioHumSflux', 'g m-2 s-1', 'abiotic_sediment_DHumS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tNAbioHumS,     'tNAbioHumS',     'g m-2 s-1', 'abiotic_sediment_NHumS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tPAbioHumS,     'tPAbioHumS',     'g m-2 s-1', 'abiotic_sediment_NHumS_change',       output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_tSiAbioSiO2S,   'tSiAbioSiO2S',   'g m-2 s-1', 'abiotic_sediment_SiO2_change',        output=output_instantaneous)
#endif
   call self%register_diagnostic_variable(self%id_tPdifPO4,       'tPdifPO4',       'g m-2 d-1', 'abiotic_sediment_PO4W_diffusion',     output=output_instantaneous,domain=domain_bottom)
   call self%register_diagnostic_variable(self%id_tNdifNH4,       'tNdifNH4',       'g m-2 d-1', 'abiotic_sediment_NH4W_diffusion',     output=output_instantaneous,domain=domain_bottom)
   call self%register_diagnostic_variable(self%id_tNdifNO3,       'tNdifNO3',       'g m-2 d-1', 'abiotic_sediment_NO3W_diffusion',     output=output_instantaneous,domain=domain_bottom)
   call self%register_diagnostic_variable(self%id_tDAbioO2S,      'tDAbioO2S',      'g m-2 d-1', 'abiotic_sediment_O2_consumption',     output=output_instantaneous,domain=domain_bottom)
!  Register diagnostic variable for external dependencies
   call self%register_diagnostic_variable(self%id_tDAbioPOMS,     'tDAbioPOMS',     'g m-2 s-1', 'particulate organic matter change',   output=output_none)
   call self%register_diagnostic_variable(self%id_tDAbioHumS,     'tDAbioHumS',     'g m-2 s-1', 'Humus change',                        output=output_none)
!  Register environmental dependencies
   call self%register_dependency(self%id_uTm,standard_variables%temperature)
   call self%register_dependency(self%id_depth,standard_variables%bottom_depth)
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
! ! INPUT PARAMETERS:
   class (type_pclake_abiotic_sediment), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!
! !LOCAL VARIABLES:
!  carriers for local state variables
   real(rk)                   :: sNH4S,sNO3S,sPO4S,sPAIMS,sDIMS
   real(rk)                   :: sDPOMS,sNPOMS,sPPOMS,sSiPOMS
   real(rk)                   :: sDHumS,sNHumS,sPHumS
!  carriers for dissolved organic matter
   real(rk)                   :: sDDOMS,sNDOMS,sPDOMS,sSiDOMS
   real(rk)                   :: sDDOMW,sNDOMW,sPDOMW,sSiDOMW
!  nutrients ratios
   real(rk)                   :: rPDDOMS,rNDDOMS
   real(rk)                   :: rPDPOMS,rNDPOMS
!  carriers for environmental dependencies
   real(rk)                   :: uTm,depth,dz
!  carriers for diagnostic dependencies
!  in abiotic water column module
   real(rk)                   :: sO2W,sNH4W,sNO3W,sPO4W
!  variables for local processes
   real(rk)                   :: afOxySed, aDepthOxySed
   real(rk)                   :: tSOD,uFunTmNitr,tNNitrS
   real(rk)                   :: oNH4S,oNO3S,oPO4S
   real(rk)                   :: tO2NitrS
   real(rk)                   :: tNAbioNO3S,tNDenitS,tDDenitS
   real(rk)                   :: tNAbioNH4S
   real(rk)                   :: tPAbioPO4S,tPSorpIMS,aPEqIMS
   real(rk)                   :: aPIsoAdsS,aPAdsMaxS,aKPAdsS,tPChemPO4,tPAbioAIMS
   real(rk)                   :: tDAbioIMS
!  variables for diffusion (in the order of appearance)
   real(rk)                   :: aDepthDif,tNDifNH4,tNDifNO3,tPDifPO4
   real(rk)                   :: akO2DifCor,tO2Dif,uFunTmDif
!  Variables for humus
   real(rk)    :: tDMinHumS,tNMinHumS,tPMinHumS
   real(rk)    :: tDAbioHumS,tNAbioHumS,tPAbioHumS
!  variables for mineralization
!  step 1 min. variables (POM mineralization)
   real(rk)                   :: uFunTmMinPOMS
   real(rk)                   :: tDMinPOMS,tNMinPOMS,tPMinPOMS,tSiMinPOMS
!  step 2 min. variables (DOM mineralization)
   real(rk)                   :: uFunTmMinDOMS
   real(rk)                   :: tDMinDOMS,tNMinDOMS, tPMinDOMS,tSiMinDOMS
!  dissolved organic matter diffusion variables
   real(rk)                   :: oDDOMS,oNDOMS,oPDOMS,oSiDOMS
   real(rk)                   :: tDDifDOM,tNDifDOM,tPDifDOM,tSiDifDOM
!  sediment oxygen consumption variables
   real(rk)                   :: tO2MinPOMS, tDMinOxyPOMS
!  sediment oxygen consumption variables for dissolved organic matter
   real(rk)                   :: tDMinOxyDOMS,tO2MinDOMS
!  total sediment oxygen consumpution
   real(rk)                   :: tO2MinTOMS
!  total fluxes for sediment change
   real(rk)                   :: tDAbioPOMS,tNAbioPOMS,tPAbioPOMS,tSiAbioPOMS
!  total fluxes for dissolved organic matter
   real(rk)                   :: tDAbioDOMS,tNAbioDOMS,tPAbioDOMS,tSiAbioDOMS
   
!EOP
!-----------------------------------------------------------------------
!BOC

!  Enter spatial loops
   _FABM_HORIZONTAL_LOOP_BEGIN_
!  Retrieve current (local) state variable values.
   _GET_HORIZONTAL_(self%id_sDIMS,sDIMS)
   _GET_HORIZONTAL_(self%id_sDPOMS,sDPOMS)
   _GET_HORIZONTAL_(self%id_sNPOMS,sNPOMS)
   _GET_HORIZONTAL_(self%id_sPPOMS,sPPOMS)
   _GET_HORIZONTAL_(self%id_sSiPOMS,sSiPOMS)
   _GET_HORIZONTAL_(self%id_sPO4S,sPO4S)
   _GET_HORIZONTAL_(self%id_sPAIMS,sPAIMS)
   _GET_HORIZONTAL_(self%id_sNH4S,sNH4S)
   _GET_HORIZONTAL_(self%id_sNO3S,sNO3S)
!  Retrieve dissolved organic state variables
   _GET_HORIZONTAL_(self%id_sDDOMS,sDDOMS)
   _GET_HORIZONTAL_(self%id_sNDOMS,sNDOMS)
   _GET_HORIZONTAL_(self%id_sPDOMS,sPDOMS)
   _GET_HORIZONTAL_(self%id_sSiDOMS,sSiDOMS)
!  Humus
   _GET_HORIZONTAL_(self%id_sDHumS,sDHumS)
   _GET_HORIZONTAL_(self%id_sNHumS,sNHumS)
   _GET_HORIZONTAL_(self%id_sPHumS,sPHumS)
!  Retrieve dependencis value
!  from abiotic water module
   _GET_(self%id_O2ConsumpSed,sO2W)
   _GET_(self%id_diffNH4,sNH4W)
   _GET_(self%id_diffNO3,sNO3W)
   _GET_(self%id_diffPO4,sPO4W)
!  retrieve dissolved organic matter in the water column
   _GET_(self%id_diffDDOM, sDDOMW)
   _GET_(self%id_diffNDOM, sNDOMW)
   _GET_(self%id_diffPDOM, sPDOMW)
   _GET_(self%id_diffSiDOM,sSiDOMW)
!  retrieve environmental dependencies
   _GET_(self%id_uTm,uTm)
   _GET_HORIZONTAL_(self%id_depth,depth)
   _GET_(self%id_dz,dz)
!-----------------------------------------------------------------------
!  Current local nutrients ratios (check the current state)
!-----------------------------------------------------------------------
   rPDDOMS=sPDOMS/(sDDOMS+NearZero)
   rNDDOMS=sNDOMS/(sDDOMS+NearZero)
   rPDPOMS=sPPOMS/(sDPOMS+NearZero)
   rNDPOMS=sNPOMS/(sDPOMS+NearZero)
!-----------------------------------------------------------------------
!  Temperature functions for sediment abiotic process
!-----------------------------------------------------------------------
!  Temperature_dependence_for_nitrification
   uFunTmNitr=uFunTmAbio(uTm,self%cThetaNitr)
!  temperature_function_of_diffusion
   uFunTmDif= uFunTmAbio(uTm,self%cThetaDif)
!-----------------------------------------------------------------------
!  dissolved nutrients concentration in sediment (converting)
!-----------------------------------------------------------------------
!  conc._dissolved_N-NO3_in_interstitial_water
   oNO3S = max(0.0_rk, sNO3S / self%cDepthS / self%bPorS)
!  conc._dissolved_N-NH4_in_interstitial_water
   oNH4S=max(0.0_rk, sNH4S/self%cDepthS/self%bPorS)
!  conc._dissolved_P_in_interstitial_water
   oPO4S = max(0.0_rk, sPO4S / self%cDepthS / self%bPorS)
!-----------------------------------------------------------------------
!  mineralization from POM to DOM
!-----------------------------------------------------------------------
!  temperature function from POM to DOM
   uFunTmMinPOMS=uFunTmAbio(uTm,self%cThetaMinPOMS)
!  decomposition_of_upper_sediment
   tDMinPOMS=self%kDMinPOMS*uFunTmMinPOMS*sDPOMS
!  mineralization_of_P_in_upper_sediment
   tPMinPOMS=self%kPMinPOMS*uFunTmMinPOMS*sPPOMS
!  mineralization_of_N_in_upper_sediment
   tNMinPOMS=self%kNMinPOMS*uFunTmMinPOMS*sNPOMS
!  mineralization_of_Si_in_upper_sediment
   tSiMinPOMS=self%kSiMinPOMS*uFunTmMinPOMS*sSiPOMS
!-----------------------------------------------------------------------
!  mineralization DOM to dissolved nutrients
!-----------------------------------------------------------------------
!  temp._function_of_mineralization
   uFunTmMinDOMS=uFunTmAbio(uTm,self%cThetaMinDOMS)
!  decomposition_of_upper_sediment
   tDMinDOMS=self%kDMinDOMS*uFunTmMinDOMS*sDDOMS
!  mineralization_of_P_in_upper_sediment
   tPMinDOMS=self%kPMinDOMS*uFunTmMinDOMS*sPDOMS
!  mineralization_of_N_in_upper_sediment
   tNMinDOMS=self%kNMinDOMS*uFunTmMinDOMS*sNDOMS
!  mineralization_of_Si_in_upper_sediment
   tSiMinDOMS=self%kSiMinDOMS*uFunTmMinDOMS*sSiDOMS
!-----------------------------------------------------------------------
!  Diffusion process
!-----------------------------------------------------------------------
!  average_diffusion_distance
   aDepthDif=self%fDepthDifS*self%cDepthS
!  diffusion_flux_of_NH4_from_sediment_to_water
   tNDifNH4=self%kNDifNH4*uFunTmDif*self%cTurbDifNut*self%bPorCorS*(oNH4S-sNH4W)/aDepthDif
!  diffusion_flux_of_NO3_from_sediment_to_water
   tNDifNO3=self%kNDifNO3*uFunTmDif*self%cTurbDifNut*self%bPorCorS*(oNO3S-sNO3W)/aDepthDif
!  diffusion_flux_of_dissolved_P_from_sediment_to_water
   tPDifPO4=self%kPDifPO4*uFunTmDif*self%cTurbDifNut*self%bPorCorS*(oPO4S-sPO4W)/aDepthDif
!  corrected_O2_diffusion_coefficient
   akO2DifCor=self%kO2Dif*uFunTmDif*self%cTurbDifO2*self%bPorCorS
!  O2_diffusion_(water_->_sediment)
   tO2Dif= akO2DifCor*sO2W/aDepthDif
!-----------------------------------------------------------------------
!  diffusion process for dissolved organic matter
!-----------------------------------------------------------------------
!  conc._dissolved_organic matter_in_interstitial_water
   oDDOMS  = sDDOMS  / self%cDepthS / self%bPorS
   oNDOMS  = sNDOMS  / self%cDepthS / self%bPorS
   oPDOMS  = sPDOMS  / self%cDepthS / self%bPorS
   oSiDOMS = sSiDOMS / self%cDepthS / self%bPorS
!  diffusion fluxes of dissolved organic matter
   tDDifDOM =self%kDDifDOM *uFunTmDif*self%cTurbDifNut*self%bPorCorS*(oDDOMS-sDDOMW)/aDepthDif
   tNDifDOM =self%kNDifDOM *uFunTmDif*self%cTurbDifNut*self%bPorCorS*(oNDOMS-sNDOMW)/aDepthDif
   tPDifDOM =self%kPDifDOM *uFunTmDif*self%cTurbDifNut*self%bPorCorS*(oPDOMS-sPDOMW)/aDepthDif
   tSiDifDOM=self%kSiDifDOM*uFunTmDif*self%cTurbDifNut*self%bPorCorS*(oSiDOMS-sSiDOMW)/aDepthDif
!-----------------------------------------------------------------------
!  Oxygen conditions in sediment
!-----------------------------------------------------------------------
!  sediment_oxygen_demand
!   The original equation (before splitting OM in POM and DOM fraction)
!   tSOD=(molO2molC*self%cCPerDW*(1.0_rk-self%fRefrDetS)*tDMinDetS+self%O2PerNH4*molO2molN*self%kNitrS*uFunTmNitr*sNH4S)/self%cDepthS
   tSOD=(molO2molC*self%cCPerDW*((1.0_rk-self%fRefrPOMS)*tDMinPOMS+tDMinDOMS)+self%O2PerNH4*molO2molN*self%kNitrS*uFunTmNitr*sNH4S)/self%cDepthS
!  oxygen_penetration_depth
   aDepthOxySed=(((2.0_rk * sO2W * akO2DifCor / tSOD) )** (0.5_rk))
!  fraction_aerobic_sediment
   afOxySed=min(1.0, max(0.0_rk, aDepthOxySed/self%cDepthS))
!  aerobic_mineralization_of POM
   tDMinOxyPOMS=afOxySed*(1.0_rk-self%fRefrPOMS)*tDMinPOMS
!  aerobic_mineralization of DOM
   tDMinOxyDOMS=afOxySed*tDMinDOMS
!  sediment_oxygen_demand
!  The original equation (before splitting OM in POM and DOM fraction)
!  tO2MinTOMS=molO2molC*self%cCPerDW*tDMinOxyDetS
   tO2MinTOMS=molO2molC*self%cCPerDW*(tDMinOxyPOMS+tDMinOxyDOMS)
!-----------------------------------------------------------------------
!  denitrification flux
!-----------------------------------------------------------------------
!  mineralization_flux_by_denitrification
!  The original equation (before splitting OM in POM and DOM fraction)
!   tDDenitS=oNO3S*oNO3S/(self%hNO3Denit*self%hNO3Denit+oNO3S*oNO3S)*(1.0-afOxySed)*(1.0_rk-self%fRefrDetS)*tDMinDetS
   tDDenitS=oNO3S*oNO3S/(self%hNO3Denit*self%hNO3Denit+oNO3S*oNO3S)*(1.0-afOxySed)*(tDMinDOMS+(1.0_rk-self%fRefrPOMS)*tDMinPOMS)
!  Denitrification_flux
   tNDenitS=self%NO3PerC*molNmolC*self%cCPerDW*tDDenitS
!-----------------------------------------------------------------------
!  nitrification flux
!-----------------------------------------------------------------------
!  nitrification_flux
   tNNitrS=afOxySed*self%kNitrS*uFunTmNitr*sNH4S
!  O2_flux_due_to_nitrification
   tO2NitrS=self%O2PerNH4*molO2molN*tNNitrS
!-----------------------------------------------------------------------
!  absorbed P in sediment, oxygen dependent
!-----------------------------------------------------------------------
!  max._P_adsorption_per_g_inorg._matter_in_sediment
   aPAdsMaxS =self%cRelPAdsD+afOxySed*self%cRelPAdsFe*self%fFeDIM+self%cRelPAdsAl*self%fAlDIM
!  P_adsorption_affinity,_corrected_for_redox_conditions
   aKPAdsS=(1.0_rk-self%fRedMax*(1.0_rk-afOxySed))*self%cKPAdsOx
!  P_adsorption_isotherm_onto_inorg._matter_in_sediment
   aPIsoAdsS=aPAdsMaxS*aKPAdsS*oPO4S/(1.0_rk+aKPAdsS*oPO4S)
!  equilibrium_amount
   aPEqIMS = aPIsoAdsS * sDIMS
!  sorption
   tPSorpIMS=self%kPSorp*(aPEqIMS-sPAIMS)
!  chem._loss_of_dissolved_P_from_pore_water
   tPChemPO4=max( 0.0_rk,self%kPChemPO4*(oPO4S-self%coPO4Max))
!  decomposition_of_upper_sediment_humus
   tDMinHumS = self%kDMinHum * uFunTmMinPOMS * afOxySed * sDHumS
!  mineralization_of_P_in_upper_sediment_humus
   tPMinHumS = self%kDMinHum * uFunTmMinPOMS * afOxySed * sPHumS
!  mineralization_of_N_in_upper_sediment_humus
   tNMinHumS = self%kDMinHum * uFunTmMinPOMS * afOxySed * sNHumS
!-----------------------------------------------------------------------
!  total abiotic flux for each state variable in sediment
!-----------------------------------------------------------------------
!  total_abiotic/microbial_DW_inorganic_matter_flux_in_sediment
   tDAbioIMS=0.0_rk
!  total_abiotic/microbial_DW_POM_flux_in_sediment
   tDAbioPOMS=-tDMinPOMS
!  total_abiotic/microbial_P_POM_flux_in_sediment
   tPAbioPOMS =-tPMinPOMS
!  total_abiotic/microbial_dissolved_P_flux_in_sediment
   tPAbioPO4S= tPMinDOMS + tPMinHumS - tPDifPO4 - tPSorpIMS - tPChemPO4
!  total_abiotic/microbial_P_absorbed_onto_inorganic_matter_flux_in_sediment
   tPAbioAIMS=tPSorpIMS
!  total_abiotic/microbial_N_NH4_flux_in_sediment
   tNAbioNH4S=tNMinDOMS +tNMinHumS -tNDifNH4 -tNNitrS
!  total_abiotic/microbial_N_NO3_flux_in_sediment
   tNAbioNO3S= tNNitrS-tNDenitS-tNDifNO3
!  total_abiotic/microbial_N_POM_flux_in_sediment
   tNAbioPOMS =-tNMinPOMS
!  total_abiotic/microbial_Si_POM_flux_in_sediment
   tSiAbioPOMS =-tSiMinPOMS
!  Humus process
!  total_abiotic/microbial_DW_humus_flux_in_sediment
   tDAbioHumS = self%fRefrPOMS * tDMinPOMS - tDMinHumS
!  total_abiotic/microbial_N_humus_flux_in_sediment
   tNAbioHumS = self%fRefrPOMS * tNMinPOMS - tNMinHumS
!  total_abiotic/microbial_P_humus_flux_in_sediment
   tPAbioHumS = self%fRefrPOMS * tPMinPOMS - tPMinHumS
!  total abiotic fluxes for dissolved organic matter
   tDAbioDOMS  = tDMinPOMS -tDMinDOMS - tDDifDOM
   tNAbioDOMS  = tNMinPOMS -tNMinDOMS - tNDifDOM
   tPAbioDOMS  = tPMinPOMS -tPMinDOMS - tPDifDOM
   tSiAbioDOMS = tSiMinPOMS -tSiMinDOMS - tSiDifDOM
!-----------------------------------------------------------------------
!  update the state variables
!-----------------------------------------------------------------------
   _SET_ODE_BEN_(self%id_sDIMS,tDAbioIMS)
   _SET_ODE_BEN_(self%id_sDPOMS,tDAbioPOMS)
   _SET_ODE_BEN_(self%id_sPPOMS,tPAbioPOMS)
   _SET_ODE_BEN_(self%id_sNPOMS,tNAbioPOMS)
   _SET_ODE_BEN_(self%id_sSiPOMS,tSiAbioPOMS)
   _SET_ODE_BEN_(self%id_sNH4S,tNAbioNH4S)
   _SET_ODE_BEN_(self%id_sNO3S,tNAbioNO3S)
   _SET_ODE_BEN_(self%id_sPO4S,tPAbioPO4S)
   _SET_ODE_BEN_(self%id_sPAIMS,tPAbioAIMS)
   _SET_ODE_BEN_(self%id_sDHumS,tDAbioHumS)
   _SET_ODE_BEN_(self%id_sPHumS,tPAbioHumS)
   _SET_ODE_BEN_(self%id_sNHumS,tNAbioHumS)
!-----------------------------------------------------------------------
!  update the state variables for dissolved organic matter
!-----------------------------------------------------------------------
   _SET_ODE_BEN_(self%id_sDDOMS, tDAbioDOMS)
   _SET_ODE_BEN_(self%id_sPDOMS, tPAbioDOMS)
   _SET_ODE_BEN_(self%id_sNDOMS, tNAbioDOMS)
   _SET_ODE_BEN_(self%id_sSiDOMS,tSiAbioDOMS)
!-----------------------------------------------------------------------
!  Update external state variables
!-----------------------------------------------------------------------
   _SET_BOTTOM_EXCHANGE_(self%id_MinSiO2Sed,tSiMinDOMS)
   _SET_BOTTOM_EXCHANGE_(self%id_diffNH4, tNdifNH4)
   _SET_BOTTOM_EXCHANGE_(self%id_diffNO3,tNdifNO3)
   _SET_BOTTOM_EXCHANGE_(self%id_diffPO4,tPdifPO4)
   _SET_BOTTOM_EXCHANGE_(self%id_O2ConsumpSed,-tO2MinTOMS - tO2NitrS)
!-----------------------------------------------------------------------
!  Update external state variables for dissolved organic matter
!-----------------------------------------------------------------------
   _SET_BOTTOM_EXCHANGE_(self%id_diffDDOM,tDDifDOM)
   _SET_BOTTOM_EXCHANGE_(self%id_diffNDOM,tNDifDOM)
   _SET_BOTTOM_EXCHANGE_(self%id_diffPDOM,tPDifDOM)
   _SET_BOTTOM_EXCHANGE_(self%id_diffSiDOM,tSiDifDOM)
!-----------------------------------------------------------------------
!  Output dependent diagnostic variables for other modules
!-----------------------------------------------------------------------
!  output local diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_afOxySed,afOxySed)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aPEqIMS,aPEqIMS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rPDPOMS,rPDPOMS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rNDPOMS,rNDPOMS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rPDDOMS,rPDDOMS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_rNDDOMS,rNDDOMS)
!  output diagnostic values for external usage
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioPOMS,tDAbioPOMS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioHumS,tDAbioHumS)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPdifPO4,tPdifPO4*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNdifNH4,tNdifNH4*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNdifNO3,tNdifNO3*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioO2S,(-tO2MinTOMS - tO2NitrS)*secs_pr_day)
#ifdef _DEVELOPMENT_
!  output diagnostic variables for modular fluxes
!  total fluxes for local sediment state variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioIMS,tDAbioIMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioPOMSflux,tDAbioPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAbioPOMS,tNAbioPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAbioPOMS,tPAbioPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tSiAbioPOMS,tSiAbioPOMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAbioAIMS,tPAbioAIMS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tDAbioHumSflux,tDAbioHumS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAbioHumS,tNAbioHumS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAbioHumS,tPAbioHumS*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAbioNH4S,tNAbioNH4S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tNAbioNO3S,tNAbioNO3S*secs_pr_day)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPAbioPO4S,tPAbioPO4S*secs_pr_day)
!  total fluxes for abiotic water state variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tSiAbioSiO2S,tSiMinDOMS*secs_pr_day)
#endif


   _FABM_HORIZONTAL_LOOP_END_
!-----------------------------------------------------------------------
! Spatial loop end
!-----------------------------------------------------------------------

   end subroutine do_bottom

!EOC
!-----------------------------------------------------------------------

   end module pclake_abiotic_sediment

!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
