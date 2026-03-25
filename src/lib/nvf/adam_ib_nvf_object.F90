!< ADAM, IB class NVF (NVF backend of [[ib_object]]).
module adam_ib_nvf_object
!< ADAM, IB class NVF (NVF backend of [[ib_object]]).

use adam_ib_object
use adam_grid_global, only: grid
use adam_ib_nvf_kernels
use adam_field_nvf_object
use adam_memory_nvf_library
use adam_mpih_nvf_object
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
   type(field_nvf_object), pointer :: field_gpu=>null() !< Field NVF handler.
   ! device data
   real(R8P), allocatable, device :: q_bcs_vars_gpu(:,:) !< Variables array for immersed boundary on GPU.
   real(R8P), allocatable, device :: phi_gpu(:,:,:,:,:)  !< Distance function on GPU.
   ! grid/field data replica for easy handling
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P), pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P), pointer :: nv=>null()            !< Number of conservative variables.
   contains
      ! public methods
      procedure, pass(self) :: evolve_eikonal !< Evolve eikonal equation.
      procedure, pass(self) :: initialize     !< Initialize class.
      procedure, pass(self) :: invert_eikonal !< Invert momentum eikonal equation.
endtype ib_nvf_object
contains
   ! public methods
   subroutine evolve_eikonal(self, dq_gpu, q_gpu)
   !< Evolve eikonal equation.
   class(ib_nvf_object), intent(in)            :: self       !< IB.
   real(R8P),            intent(inout), device :: dq_gpu(1:,         &
                                                         1-self%ngc:,&
                                                         1-self%ngc:,&
                                                         1-self%ngc:,&
                                                         1:) !< State variables variations.
   real(R8P),            intent(inout), device :: q_gpu(1:,         &
                                                        1-self%ngc:,&
                                                        1-self%ngc:,&
                                                        1-self%ngc:,&
                                                        1:)  !< Conservative variables.
   integer(I4P)                                :: ib         !< Counter.

   associate(blocks_number=>self%blocks_number, ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv,          &
             dx_gpu=>self%field_gpu%dxyz_gpu(:,1),dy_gpu=>self%field_gpu%dxyz_gpu(:,2),dz_gpu=>self%field_gpu%dxyz_gpu(:,3),&
             solids_number=>self%ib%solids_number, phi_gpu=>self%phi_gpu)
      do ib=1, solids_number
         call compute_eikonal_dq_phi_cuf(ib=ib, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                         dx_gpu=dx_gpu, dy_gpu=dy_gpu, dz_gpu=dz_gpu,                             &
                                         phi_gpu=phi_gpu, dq_gpu=dq_gpu, q_gpu=q_gpu)
         call self%mpih%check_cuda_error(error_code=CUDA_SILENT_ERR, msg='CUDA error in compute eikonal dq phi')
         call evolve_eikonal_q_phi_cuf(ib=ib, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                       phi_gpu=phi_gpu, dq_gpu=dq_gpu, q_gpu=q_gpu)
         call self%mpih%check_cuda_error(error_code=CUDA_SILENT_ERR, msg='CUDA error in evolve eikonal q phi')
      enddo
   endassociate
   endsubroutine evolve_eikonal

   subroutine initialize(self, ib, field_gpu)
   !< Initialize class.
   class(ib_nvf_object),   intent(inout)      :: self      !< IB NVF object.
   type(ib_object),        intent(in), target :: ib        !< IB object.
   type(field_nvf_object), intent(in), target :: field_gpu !< The field.
   character(:), allocatable                  :: msg_      !< Allocating message base.
   character(:), allocatable                  :: ms        !< Allocating message.

   call self%mpih%initialize(do_mpi_init=.false.)
   call self%mpih%print_message('ib_nvf_object%initialize start')
   self%ib => ib
   self%field_gpu => field_gpu
   msg_ = self%mpih%myrankstr//'ib_nvf_object%initialize'
   call assign_allocatable_gpu(lhs=self%q_bcs_vars_gpu, rhs=self%ib%q, msg=msg_//' q_bcs_vars_gpu ')
   associate(ngc=>grid%ngc, ni=>grid%ni, nj=>grid%nj, nk=>grid%nk, nb=>self%ib%field%nb, &
             solids_number=>self%ib%solids_number)
   if (solids_number>0) then
      ms = msg_//' phi_gpu '
      call alloc_var_gpu(var=self%phi_gpu,ulb=reshape([1,nb,1-ngc,ni+ngc,1-ngc,nj+ngc,1-ngc,nk+ngc,1,solids_number+1],[2,5]),msg=ms)
      self%phi_gpu = -1._R8P
   endif
   endassociate
   call associate_adam_data(ib=ib)
   contains
      subroutine associate_adam_data(ib)
      !< Associate objects data to equation for easy handling.
      type(ib_object), intent(in), target :: ib !< IB object.

      self%blocks_number => ib%field%blocks_number
      self%ni            => grid%ni
      self%nj            => grid%nj
      self%nk            => grid%nk
      self%ngc           => grid%ngc
      self%nb            => ib%field%nb
      self%nv            => ib%field%nv
      endsubroutine associate_adam_data
   endsubroutine initialize

   subroutine invert_eikonal(self, q_gpu)
   !< Invert momentum eikonal equation.
   class(ib_nvf_object), intent(in)            :: self       !< IB.
   real(R8P),            intent(inout), device :: q_gpu(1:,         &
                                                        1-self%ngc:,&
                                                        1-self%ngc:,&
                                                        1-self%ngc:,&
                                                        1:) !< Conservative variables.
   integer(I4P)                                 :: ib       !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, ngc=>self%ngc, nv=>self%nv, blocks_number=>self%blocks_number, &
             bcs_type=>self%ib%bc_type, solids_number=>self%ib%solids_number, phi_gpu=>self%phi_gpu)
      do ib=1, solids_number
         call invert_eikonal_q_phi_cuf(BCS_VISCOUS=BCS_VISCOUS, BCS_EULER=BCS_EULER,                            &
                                       ib=ib, ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                       bcs_type=bcs_type(ib), phi_gpu=phi_gpu, q_gpu=q_gpu)
         call self%mpih%check_cuda_error(error_code=CUDA_SILENT_ERR, msg='CUDA error in invert eikonal q phi')
      enddo
   endassociate
   endsubroutine invert_eikonal
endmodule adam_ib_nvf_object
