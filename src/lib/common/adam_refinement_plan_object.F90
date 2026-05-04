!< ADAM, refinement plan class definition.
module adam_refinement_plan_object
!< ADAM, refinement plan class definition.
!<
!< A refinement plan is a one-shot transfer object produced by `tree_object%adapt`
!< and consumed by `field_object%adapt` in a single AMR cycle. It carries the
!< field-block lists that the tree computes during refinement/derefinement, keeping
!< them out of the tree's persistent state and decoupling tree_object from
!< field_object.

! third party modules
use :: penf

implicit none
private
public :: refinement_plan_object

type :: refinement_plan_object
   !< Refinement plan class definition.
   integer(I4P)              :: ratio=0_I4P          !< Refinement ratio (2=1D, 4=2D, 8=3D).
   integer(I8P), allocatable :: block_to_refine(:,:) !< Field blocks to be refined [2, n_my_refine].
   integer(I8P), allocatable :: block_refined(:,:)   !< Field refined blocks with Morton code [2, ratio*n_my_refine].
   integer(I8P), allocatable :: block_to_derefine(:) !< Field blocks to be derefined [n_my_derefine].
   integer(I8P), allocatable :: block_derefined(:,:) !< Field derefined blocks with Morton code [2, n_my_derefine/ratio].
   contains
      procedure, pass(self) :: destroy           !< Destroy plan.
      procedure, pass(self) :: has_refinements   !< Return true if there are local blocks to refine.
      procedure, pass(self) :: has_derefinements !< Return true if there are local blocks to derefine.
endtype refinement_plan_object

contains
   pure subroutine destroy(self)
   !< Destroy plan.
   class(refinement_plan_object), intent(inout) :: self !< The refinement plan.

   self%ratio = 0_I4P
   if (allocated(self%block_to_refine  )) deallocate(self%block_to_refine  )
   if (allocated(self%block_refined    )) deallocate(self%block_refined    )
   if (allocated(self%block_to_derefine)) deallocate(self%block_to_derefine)
   if (allocated(self%block_derefined  )) deallocate(self%block_derefined  )
   endsubroutine destroy

   pure function has_refinements(self) result(res)
   !< Return true if there are local blocks to refine.
   class(refinement_plan_object), intent(in) :: self !< The refinement plan.
   logical                                   :: res  !< Result.

   res = allocated(self%block_to_refine) .and. size(self%block_to_refine, dim=2) > 0
   endfunction has_refinements

   pure function has_derefinements(self) result(res)
   !< Return true if there are local blocks to derefine.
   class(refinement_plan_object), intent(in) :: self !< The refinement plan.
   logical                                   :: res  !< Result.

   res = allocated(self%block_to_derefine) .and. size(self%block_to_derefine) > 0
   endfunction has_derefinements
endmodule adam_refinement_plan_object
