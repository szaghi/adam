!< ADAM, IB class NVF (NVF backend of [[ib_object]]).
module adam_ib_nvf_object
!< ADAM, IB class NVF (NVF backend of [[ib_object]]).

use adam_memory_nvf_library
use adam_mpih_nvf_object
use adam_ib_object
use penf

implicit none
save
private
public :: ib_nvf_object

type :: ib_nvf_object
   !< IB NVF class definition.
   ! ADAM library objects
   type(ib_object), pointer :: ib=>null() !< IB common handler.
   ! ADAM NVF library objects
   type(mpih_nvf_object) :: mpih !< MPI handler.
   ! device data
   real(R8P), allocatable, device :: q_bcs_vars_gpu(:,:) !< Variables array for immersed boundary on GPU.
   real(R8P), allocatable, device :: phi_gpu(:,:,:,:,:)  !< Distance function on GPU.
   contains
      ! public methods
      procedure, pass(self) :: initialize !< Initialize class.
endtype ib_nvf_object
contains
   ! public methods
   subroutine initialize(self, ib)
   !< Initialize class.
   class(ib_nvf_object), intent(inout)      :: self !< IB NVF object.
   type(ib_object),      intent(in), target :: ib   !< IB object.
   character(:), allocatable                :: msg_ !< Allocating message base.
   character(:), allocatable                :: ms   !< Allocating message.

   call self%mpih%initialize(do_mpi_init=.false.)
   call self%mpih%print_message('ib_nvf_object%initialize start')
   self%ib => ib
   msg_ = self%mpih%myrankstr//'ib_nvf_object%initialize'
   call assign_allocatable_gpu(lhs=self%q_bcs_vars_gpu, rhs=self%ib%q, msg=msg_//' q_bcs_vars_gpu ')
   associate(ngc=>self%ib%grid%ngc, ni=>self%ib%grid%ni, nj=>self%ib%grid%nj, nk=>self%ib%grid%nk, nb=>self%ib%field%nb, &
             solids_number=>self%ib%solids_number)
   if (solids_number>0) then
      ms = msg_//' phi_gpu '
      call alloc_var_gpu(var=self%phi_gpu,ulb=reshape([1,nb,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,solids_number+1],[2,5]),msg=ms)
      self%phi_gpu = -1._R8P
   endif
   endassociate
   endsubroutine initialize
endmodule adam_ib_nvf_object
