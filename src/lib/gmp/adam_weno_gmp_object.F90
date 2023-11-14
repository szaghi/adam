!< ADAM, WENO class GMP (GMP backend of [[weno_object]]).
module adam_weno_gmp_object
!< ADAM, WENO class GMP (GMP backend of [[weno_object]]).

use adam_memory_gmp_library
use adam_mpih_gmp_object
use adam_weno_object
use finer
use penf

implicit none
save
private
public :: weno_gmp_object

type :: weno_gmp_object
   !< WENO GMP class definition.
   ! ADAM library objects
   type(weno_object), pointer :: weno=>null() !< WENO common handler.
   ! ADAM GMP library objects
   type(mpih_gmp_object), pointer :: mpih=>null() !< MPI handler.
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
endtype weno_gmp_object
contains
   ! public methods
   subroutine initialize(self, mpih, weno)
   !< Initialize class.
   class(weno_gmp_object), intent(inout)      :: self !< WENO GMP object.
   type(mpih_gmp_object),  intent(in), target :: mpih !< MPI handler, GMP backend.
   type(weno_object),      intent(in), target :: weno !< WENO object.
   character(:), allocatable                  :: msg_ !< Allocating message base.

   self%mpih => mpih
   call self%mpih%print_message('weno_gmp_object%initialize start')
   self%weno => weno
   msg_ = self%mpih%myrankstr//'weno_gmp_object%initialize'
   call assign_allocatable_gpu(lhs=self%a_gpu,          rhs=weno%a,          omp_dev=self%mpih%mydev,msg=msg_//' a_gpu '          )
   call assign_allocatable_gpu(lhs=self%p_gpu,          rhs=weno%p,          omp_dev=self%mpih%mydev,msg=msg_//' p_gpu '          )
   call assign_allocatable_gpu(lhs=self%d_gpu,          rhs=weno%d,          omp_dev=self%mpih%mydev,msg=msg_//' d_gpu '          )
   call assign_allocatable_gpu(lhs=self%ror_schemes_gpu,rhs=weno%ror_schemes,omp_dev=self%mpih%mydev,msg=msg_//' ror_schemes_gpu ')
   call assign_allocatable_gpu(lhs=self%ror_ivar_gpu,   rhs=weno%ror_ivar,   omp_dev=self%mpih%mydev,msg=msg_//' ror_ivar_gpu '   )
   call assign_allocatable_gpu(lhs=self%ror_stats_gpu,  rhs=weno%ror_stats,  omp_dev=self%mpih%mydev,msg=msg_//' ror_stats_gpu '  )
   call assign_allocatable_gpu(lhs=self%cell_scheme_gpu,rhs=weno%cell_scheme,omp_dev=self%mpih%mydev,msg=msg_//' cell_scheme_gpu ')
   call self%mpih%print_message('weno_gmp_object%initialize finish')
   endsubroutine initialize
endmodule adam_weno_gmp_object
