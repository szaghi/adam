!< ADAM, WENO class FNL (FNL backend of [[weno_object]]).
module adam_weno_fnl_object
!< ADAM, WENO class FNL (FNL backend of [[weno_object]]).

use adam_mpih_fnl_object
use adam_weno_object
use fundal
use penf

implicit none
save
private
public :: weno_fnl_object

type :: weno_fnl_object
   !< WENO FNL class definition.
   ! ADAM library objects
   type(weno_object), pointer :: weno=>null() !< WENO common handler.
   ! ADAM FNL library objects
   type(mpih_fnl_object) :: mpih !< MPI handler.
   ! device data
   real(R8P),    pointer :: a_gpu(:,:,:)               !< Optimal weights                    [1:2,0:S-1,1:S].
   real(R8P),    pointer :: p_gpu(:,:,:,:)             !< Polinomials coefficients           [1:2,0:S-1,0:S-1,1:S].
   real(R8P),    pointer :: d_gpu(:,:,:,:)             !< Smoothness indicators coefficients [0:S-1,0:S-1,0:S-1,1:S].
   integer(I4P), pointer :: ror_schemes_gpu(:)         !< Scheme (S value) for each ROR step.
   integer(I4P), pointer :: ror_ivar_gpu(:)            !< Index variables to check in ROR.
   integer(I4P), pointer :: ror_stats_gpu(:,:,:,:,:)   !< Scheme (S value) for each ROR step.
   integer(I4P), pointer :: cell_scheme_gpu(:,:,:,:,:) !< Modified order close to solids (GPU variable).
   contains
      ! public methods
      procedure, pass(self) :: initialize !< Initialize class.
endtype weno_fnl_object
contains
   ! public methods
   subroutine initialize(self, weno)
   !< Initialize class.
   class(weno_fnl_object), intent(inout)      :: self !< WENO FNL object.
   type(weno_object),      intent(in), target :: weno !< WENO object.

   call self%mpih%initialize(do_mpi_init=.false.)
   call self%mpih%print_message('weno_fnl_object%initialize start')
   self%weno => weno
   call dev_assign_to_device(dst=self%a_gpu,           src=weno%a          )
   call dev_assign_to_device(dst=self%p_gpu,           src=weno%p          )
   call dev_assign_to_device(dst=self%d_gpu,           src=weno%d          )
   call dev_assign_to_device(dst=self%ror_schemes_gpu, src=weno%ror_schemes)
   call dev_assign_to_device(dst=self%ror_ivar_gpu,    src=weno%ror_ivar   )
   call dev_assign_to_device(dst=self%ror_stats_gpu,   src=weno%ror_stats  )
   call dev_assign_to_device(dst=self%cell_scheme_gpu, src=weno%cell_scheme)
   call self%mpih%print_message('weno_fnl_object%initialize finish')
   endsubroutine initialize
endmodule adam_weno_fnl_object
