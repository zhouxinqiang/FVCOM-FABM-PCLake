!*******************************************************************************
!*                                                                             *
!* pclake_model_library.F90                                                    *
!*                                                                             *
!*******************************************************************************

module pclake_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use pclake_version
   use pclake_abiotic_water
   use pclake_abiotic_sediment
   use pclake_phytoplankton_water
   use pclake_phytoplankton_sediment
   use pclake_macrophytes
   use pclake_zooplankton
   use pclake_fish
   use pclake_zoobenthos
   use pclake_auxiliary

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: initialize
      procedure :: create
   end type

   type (type_factory),save,target,public :: pclake_model_factory

contains

   subroutine initialize(self)
      class (type_factory), intent(inout) :: self
      call self%register_version('PCLake',git_commit_id//' ('//git_branch_name//' branch)')
   end subroutine initialize

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('abiotic_water');          allocate(type_pclake_abiotic_water::model);
         case ('abiotic_sediment');       allocate(type_pclake_abiotic_sediment::model);
         case ('phytoplankton_water');    allocate(type_pclake_phytoplankton_water::model);
         case ('phytoplankton_sediment'); allocate(type_pclake_phytoplankton_sediment::model);
         case ('macrophytes');            allocate(type_pclake_macrophytes::model);
         case ('zooplankton');            allocate(type_pclake_zooplankton::model);
         case ('fish');                   allocate(type_pclake_fish::model);
         case ('zoobenthos');             allocate(type_pclake_zoobenthos::model);
         case ('auxiliary');              allocate(type_pclake_auxiliary::model);
         case default
            call self%type_base_model_factory%create(name,model)
       end select
   end subroutine create

end module

!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
