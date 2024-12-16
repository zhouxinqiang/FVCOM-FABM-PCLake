! !IROUTINE: Get the vertical movement rates (m/s) for the bio state variables.
! Note that negative values indicate movement towards the bottom, e.g., sinking,
! and positive values indicate movemment towards the surface, e.g., floating.
!
! !INTERFACE:
   subroutine fabm_get_vertical_movement(self _ARGUMENTS_INTERIOR_IN_,velocity)
!
! !INPUT PARAMETERS:
   class (type_model),               intent(inout) :: self
   _DECLARE_ARGUMENTS_INTERIOR_IN_
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk) _DIMENSION_EXT_SLICE_PLUS_1_,intent(out) :: velocity
!
! !LOCAL PARAMETERS:
   type (type_model_list_node), pointer :: node
   integer                              :: i,k
   _DECLARE_INTERIOR_INDICES_
!
!EOP
!-----------------------------------------------------------------------
!BOC
#ifndef NDEBUG
   call check_interior_location(self _ARGUMENTS_INTERIOR_IN_,'fabm_get_vertical_movement')
#  ifdef _FABM_VECTORIZED_DIMENSION_INDEX_
   call check_extents_2d(velocity,loop_stop-loop_start+1,size(self%state_variables),'fabm_get_vertical_movement','velocity','stop-start+1, # interior state variables')
#  else
   call check_extents_1d(velocity,size(self%state_variables),'fabm_get_vertical_movement','velocity','# interior state variables')
#  endif
#endif

   call prefetch_interior(self,self%get_vertical_movement_environment,self%env_int _ARGUMENTS_INTERIOR_IN_)

   ! Now allow models to overwrite with spatially-varying sinking rates - if any.
   node => self%models%first
   do while (associated(node))
      call node%model%get_vertical_movement(self%env_int)
      node => node%next
   end do

   ! Compose total sources-sinks for each state variable, combining model-specific contributions.
   do i=1,size(self%state_variables)
      k = self%state_variables(i)%movement_index
      _UNPACK_TO_PLUS_1_(self%env_int%scratch,k,velocity,i,self%env_int,0.0_rk)
   end do

   call deallocate_prefetch(self,self%get_vertical_movement_environment,self%env_int _ARGUMENTS_INTERIOR_IN_)

   end subroutine fabm_get_vertical_movement
!EOC
