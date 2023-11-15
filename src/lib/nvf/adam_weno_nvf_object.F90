!< ADAM, WENO class NVF (NVF backend of [[weno_object]]).
module adam_weno_nvf_object
!< ADAM, WENO class NVF (NVF backend of [[weno_object]]).

use adam_memory_nvf_library
use adam_mpih_nvf_object
use adam_weno_object
use finer
use penf

implicit none
save
private
public :: weno_nvf_object

type :: weno_nvf_object
   !< WENO NVF class definition.
   ! ADAM library objects
   type(weno_object), pointer :: weno=>null() !< WENO common handler.
   ! ADAM NVF library objects
   type(mpih_nvf_object) :: mpih !< MPI handler.
   ! device data
   real(R8P),    allocatable, device :: a_gpu(:,:,:)               !< Optimal weights                    [1:2,0:S-1,1:S].
   real(R8P),    allocatable, device :: p_gpu(:,:,:,:)             !< Polinomials coefficients           [1:2,0:S-1,0:S-1,1:S].
   real(R8P),    allocatable, device :: d_gpu(:,:,:,:)             !< Smoothness indicators coefficients [0:S-1,0:S-1,0:S-1,1:S].
   integer(I4P), allocatable, device :: ror_schemes_gpu(:)         !< Scheme (S value) for each ROR step.
   integer(I4P), allocatable, device :: ror_ivar_gpu(:)            !< Index variables to check in ROR.
   integer(I4P), allocatable, device :: ror_stats_gpu(:,:,:,:,:)   !< Scheme (S value) for each ROR step.
   integer(I4P), allocatable, device :: cell_scheme_gpu(:,:,:,:,:) !< Modified order close to solids (GPU variable).
   contains
      ! public methods
      procedure, pass(self) :: initialize !< Initialize class.
endtype weno_nvf_object
contains
   ! public methods
   subroutine initialize(self, weno)
   !< Initialize class.
   class(weno_nvf_object), intent(inout)      :: self !< WENO NVF object.
   type(weno_object),      intent(in), target :: weno !< WENO object.
   character(:), allocatable                  :: msg_ !< Allocating message base.

   call self%mpih%initialize(do_mpi_init=.false.)
   call self%mpih%print_message('weno_nvf_object%initialize start')
   self%weno => weno
   msg_ = self%mpih%myrankstr//'weno_nvf_object%initialize'
   call assign_allocatable_gpu(lhs=self%a_gpu,           rhs=weno%a,           msg=msg_//' a_gpu '          )
   call assign_allocatable_gpu(lhs=self%p_gpu,           rhs=weno%p,           msg=msg_//' p_gpu '          )
   call assign_allocatable_gpu(lhs=self%d_gpu,           rhs=weno%d,           msg=msg_//' d_gpu '          )
   call assign_allocatable_gpu(lhs=self%ror_schemes_gpu, rhs=weno%ror_schemes, msg=msg_//' ror_schemes_gpu ')
   call assign_allocatable_gpu(lhs=self%ror_ivar_gpu,    rhs=weno%ror_ivar,    msg=msg_//' ror_ivar_gpu '   )
   call assign_allocatable_gpu(lhs=self%ror_stats_gpu,   rhs=weno%ror_stats,   msg=msg_//' ror_stats_gpu '  )
   call assign_allocatable_gpu(lhs=self%cell_scheme_gpu, rhs=weno%cell_scheme, msg=msg_//' cell_scheme_gpu ')
   call self%mpih%print_message('weno_nvf_object%initialize finish')
   endsubroutine initialize
endmodule adam_weno_nvf_object
