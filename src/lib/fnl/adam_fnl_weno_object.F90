!< ADAM, WENO class FNL (FNL backend of [[weno_object]]).
module adam_fnl_weno_object
!< ADAM, WENO class FNL (FNL backend of [[weno_object]]).

! ADAM singleton objects
use :: adam_weno_global,     only : weno
use :: adam_fnl_mpih_global, only : mpih_fnl
! ADAM modules
use :: adam_weno_object
! third party modules
use :: fundal
use :: penf

implicit none
save
private
public :: weno_fnl_object

type :: weno_fnl_object
   !< WENO FNL class definition.
   ! device data
   ! `=> null()` is mandatory: dev_assign_to_device tests `associated(dst)`
   ! before allocating, and associated() on a pointer that was never
   ! nullified is undefined behaviour (nvfortran -Mchkptr traps it as
   ! "Null pointer").
   real(R8P),    pointer :: a_gpu(:,:,:)               => null() !< Optimal weights                    [1:2,0:S-1,1:S].
   real(R8P),    pointer :: p_gpu(:,:,:,:)             => null() !< Polinomials coefficients           [1:2,0:S-1,0:S-1,1:S].
   real(R8P),    pointer :: d_gpu(:,:,:,:)             => null() !< Smoothness indicators coefficients [0:S-1,0:S-1,0:S-1,1:S].
   integer(I4P), pointer :: ror_schemes_gpu(:)         => null() !< Scheme (S value) for each ROR step.
   integer(I4P), pointer :: ror_ivar_gpu(:)            => null() !< Index variables to check in ROR.
   integer(I4P), pointer :: ror_stats_gpu(:,:,:,:,:)   => null() !< Scheme (S value) for each ROR step.
   integer(I4P), pointer :: cell_scheme_gpu(:,:,:,:,:) => null() !< Modified order close to solids (GPU variable).
   contains
      ! public methods
      procedure, pass(self) :: initialize !< Initialize class from weno global singleton.
endtype weno_fnl_object
contains
   ! public methods
   subroutine initialize(self)
   !< Initialize class from program-scope `weno` singleton (adam_weno_global).
   !< Requires `mpih_fnl` (adam_fnl_mpih_global) and `weno` (adam_weno_global) to be
   !< initialized before calling.
   class(weno_fnl_object), intent(inout) :: self !< WENO FNL object.

   call mpih_fnl%print_message('weno_fnl_object%initialize start')
   call dev_assign_to_device(dst=self%a_gpu,           src=weno%a          )
   call dev_assign_to_device(dst=self%p_gpu,           src=weno%p          )
   call dev_assign_to_device(dst=self%d_gpu,           src=weno%d          )
   call dev_assign_to_device(dst=self%ror_schemes_gpu, src=weno%ror_schemes)
   call dev_assign_to_device(dst=self%ror_ivar_gpu,    src=weno%ror_ivar   )
   ! weno%ror_stats is allocated only when enable_ror_stats is set; otherwise
   ! it is unallocated and must not be passed to dev_assign_to_device's
   ! assumed-shape intent(in) src. With ROR stats disabled the GPU pointer
   ! stays null, the correct "nothing to offload" state.
   if (allocated(weno%ror_stats)) call dev_assign_to_device(dst=self%ror_stats_gpu, src=weno%ror_stats)
   call dev_assign_to_device(dst=self%cell_scheme_gpu, src=weno%cell_scheme)
   call mpih_fnl%print_message('weno_fnl_object%initialize finish')
   endsubroutine initialize
endmodule adam_fnl_weno_object
