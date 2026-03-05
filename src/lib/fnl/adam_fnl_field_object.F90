!< ADAM, field class definition, FNL backend.

#include "fundal.H"

module adam_fnl_field_object
!< ADAM, field class definition, FNL backend.

use adam_common_library
use adam_fnl_field_kernels
use adam_fnl_maps_object
use adam_fnl_mpih_object
use fundal
use penf
use mpi

implicit none
save
private
public :: field_fnl_object

type :: field_fnl_object
   !< Field class, FNL backend.
   type(mpih_fnl_object)       :: mpih           !< MPI handler.
   type(maps_fnl_object)       :: maps           !< Maps handler.
   type(field_object), pointer :: field=>null()  !< The field.
   ! device data
   integer(I4P), pointer :: fec_1_6_array_gpu(:)=>null() !< Mapping fec1-26 to fec1-6 for boundaries (GPU).
   ! device data copied from field object
   real(R8P), pointer :: x_cell_gpu(:,:)=>null() !< Cells x coordinates on GPU [nb,1-ngc:ni+ngc].
   real(R8P), pointer :: y_cell_gpu(:,:)=>null() !< Cells y coordinates on GPU [nb,1-ngc:nj+ngc].
   real(R8P), pointer :: z_cell_gpu(:,:)=>null() !< Cells z coordinates on GPU [nb,1-ngc:nk+ngc].
   real(R8P), pointer :: dxyz_gpu(:,:)=>null()   !< Delta cells GPU [nb,3].
   ! grid/field data replica for easy handling
   integer(I4P), pointer :: ngc=>null()           !< Number of ghost cells.
   integer(I4P), pointer :: ni=>null()            !< Number of cells in i direction.
   integer(I4P), pointer :: nj=>null()            !< Number of cells in j direction.
   integer(I4P), pointer :: nk=>null()            !< Number of cells in k direction.
   integer(I4P), pointer :: nb=>null()            !< Total blocks number for MPI.
   integer(I4P), pointer :: blocks_number=>null() !< Actual blocks number.
   integer(I4P), pointer :: nv=>null()            !< Number of variables in q vector.
   contains
      ! public methods
      procedure, pass(self) :: compute_q_gradient     !< Compute maximum gradient module of q element of a block.
      procedure, pass(self) :: copy_cpu_gpu           !< Copy data from (field_object) CPU to (field_fnl_object) GPU.
      ! procedure, pass(self) :: copy_transpose_cpu_gpu !< Transpose data from GPU to CPU.
      ! procedure, pass(self) :: copy_transpose_gpu_cpu !< Transpose data from GPU to CPU.
      procedure, pass(self) :: initialize             !< Initialize field.
      procedure, pass(self) :: update_ghost_local_gpu !< Update ghosts locally.
      procedure, pass(self) :: update_ghost_mpi_gpu   !< Update ghosts MPI.
endtype field_fnl_object

contains
   ! public methods
   subroutine compute_q_gradient(self, b, ivar, q_gpu, gradient)
   !< Compute gradient (module) over q elements.
   class(field_fnl_object), intent(in)  :: self      !< The field.
   integer(I4P),            intent(in)  :: b         !< Block index.
   integer(I4P),            intent(in)  :: ivar      !< Index of q variable.
   real(R8P),               intent(in)  :: q_gpu(1:,                    &
                                                 1-self%field%grid%ngc:,&
                                                 1-self%field%grid%ngc:,&
                                                 1-self%field%grid%ngc:,&
                                                 1:) !< Field component to which apply gradient.
   real(R8P),               intent(out) :: gradient  !< Maximum gradient of q(ivar).

   call compute_q_gradient_dev(b=b, ni=self%field%grid%ni, nj=self%field%grid%nj, nk=self%field%grid%nk, ngc=self%field%grid%ngc, &
                               dx=self%field%dxyz(1,b), dy=self%field%dxyz(2,b), dz=self%field%dxyz(3,b),     &
                               q_gpu=q_gpu, ivar=ivar, gradient=gradient)
   endsubroutine compute_q_gradient

   subroutine copy_cpu_gpu(self, verbose)
   !< Copy data from (field_object) CPU to (field_fnl_object) GPU.
   class(field_fnl_object), intent(inout)        :: self     !< The field.
   logical,                 intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                       :: verbose_ !< Flag to activate verbose mode, local var.
   character(:), allocatable                     :: r        !< My rank stringified.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih%print_message('field_fnl_object%copy_cpu_gpu start')
   r = self%mpih%myrankstr
   call self%maps%copy_cpu_gpu
   call dev_assign_to_device(src=self%field%x_cell, dst=self%x_cell_gpu, ij=[1,2])
   call dev_assign_to_device(src=self%field%y_cell, dst=self%y_cell_gpu, ij=[1,2])
   call dev_assign_to_device(src=self%field%z_cell, dst=self%z_cell_gpu, ij=[1,2])
   call dev_assign_to_device(src=self%field%dxyz,   dst=self%dxyz_gpu,   ij=[1,2])
   if (verbose_) call self%mpih%print_message('field_fnl_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine copy_transpose_cpu_gpu(self, nv, q_cpu, q_t, q_gpu)
   !< Copy transposed data from CPU to GPU.
   !< This routine is called by equation typically passing either q_gpu or q_aux_gpu.
   class(field_fnl_object), intent(inout) :: self      !< The field.
   integer(I4P),            intent(in)    :: nv        !< Number of varibales.
   real(R8P),               intent(in)    :: q_cpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:) !< Conservative variables on CPU.
   real(R8P),               intent(inout) :: q_t(1:,         &
                                                 1-self%ngc:,&
                                                 1-self%ngc:,&
                                                 1-self%ngc:,&
                                                 1:)   !< Transposed conservative variables on CPU.
   real(R8P),               intent(inout) :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:) !< Conservative variables on GPU.
   integer(I4P)                           :: i,j,k,b,v !< Counter.

   associate(blocks_number=>self%blocks_number,ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do v=1, nv
      q_t(b,i,j,k,v) = q_cpu(v,i,j,k,b)
   enddo
   enddo
   enddo
   enddo
   enddo
   call dev_memcpy_to_device(dst=q_gpu, src=q_t)
   endassociate
   endsubroutine copy_transpose_cpu_gpu

   subroutine copy_transpose_gpu_cpu(self, nv, q_gpu, q_t_gpu, q_cpu)
   !< Copy transposed data from GPU to CPU.
   !< This routine is called by equation typically passing either q_gpu or q_aux_gpu.
   class(field_fnl_object), intent(inout) :: self        !< The field.
   integer(I4P),            intent(in)    :: nv          !< Number of varibales.
   real(R8P),               intent(in)    :: q_gpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)   !< Conservative variables on GPU.
   real(R8P),               intent(inout) :: q_t_gpu(1:,         &
                                                     1-self%ngc:,&
                                                     1-self%ngc:,&
                                                     1-self%ngc:,&
                                                     1:) !< Transposed conservative variables on GPU.
   real(R8P),               intent(inout) :: q_cpu(1:,         &
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1-self%ngc:,&
                                                   1:)   !< Conservative variables on CPU.
   integer(I4P)                           :: i,j,k,b,v   !< Counter.

   associate(blocks_number=>self%blocks_number,ni=>self%ni,nj=>self%nj,nk=>self%nk,ngc=>self%ngc)
   !$acc parallel loop independent DEVICEVAR(q_gpu, q_t_gpu)
   !$omp OMPLOOP DEVICEVAR(q_gpu, q_t_gpu)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do v=1, nv
      q_t_gpu(v,i,j,k,b) = q_gpu(b,i,j,k,v)
   enddo
   enddo
   enddo
   enddo
   enddo
   call dev_memcpy_from_device(dst=q_cpu, src=q_t_gpu)
   endassociate
   endsubroutine copy_transpose_gpu_cpu

   subroutine initialize(self, field, nv_aux, q_gpu, verbose)
   !< Initialize field.
   class(field_fnl_object), intent(inout)           :: self             !< The field.
   type(field_object),      intent(in), target      :: field            !< Field variable array.
   integer(I4P),            intent(in),    optional :: nv_aux           !< Number of auxiliary variables.
   real(R8P), pointer,      intent(inout), optional :: q_gpu(:,:,:,:,:) !< Field cell centered variables.
   logical,                 intent(in),    optional :: verbose          !< Flag to activate verbose mode.
   integer(I4P)                                     :: nv_aux_          !< Number of auxiliary variables (local var).
   integer(I4P)                                     :: ierr             !< Error status.

   call self%mpih%initialize
   call self%mpih%print_message('field_fnl_object%initialize start')
   self%field         => field
   self%ngc           => field%ngc
   self%ni            => field%ni
   self%nj            => field%nj
   self%nk            => field%nk
   self%nb            => field%nb
   self%blocks_number => field%blocks_number
   self%nv            => field%nv
   call self%maps%initialize(maps=field%maps)
   associate(nb=>field%nb, ngc=>field%grid%ngc, ni=>field%grid%ni, nj=>field%grid%nj, nk=>field%grid%nk, nv=>field%nv)
      if (present(q_gpu)) &
         call dev_alloc(fptr_dev=q_gpu, ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,nv], lbounds=[1,1-ngc,1-ngc,1-ngc,1], ierr=ierr)
      call dev_assign_to_device(dst=self%fec_1_6_array_gpu, src=FEC_1_6_ARRAY)
      nv_aux_ = self%field%nv ; if (present(nv_aux)) nv_aux_ = max(nv_aux_, nv_aux)
   endassociate
   call self%copy_cpu_gpu
   call self%mpih%print_message('field_fnl_object%initialize finish')
   endsubroutine initialize

   subroutine update_ghost_local_gpu(self, q_gpu)
   !< Update (local) ghost cells.
   class(field_fnl_object), intent(in)    :: self      !< The field.
   real(R8P),               intent(inout) :: q_gpu(1:,                    &
                                                   1-self%field%grid%ngc:,&
                                                   1-self%field%grid%ngc:,&
                                                   1-self%field%grid%ngc:,&
                                                   1:) !< Field component to be updated.
   call update_ghost_local_gpu_dev(l_map_ghost_cell_gpu=self%maps%local_map_ghost_cell_gpu,ngc=self%field%grid%ngc,q_gpu=q_gpu)
   endsubroutine update_ghost_local_gpu

   subroutine update_ghost_mpi_gpu(self, q_gpu, step)
   !< Update ghost cells within other processes.
   class(field_fnl_object), intent(inout)        :: self               !< The field.
   real(R8P),               intent(inout)        :: q_gpu(1:,                    &
                                                          1-self%field%grid%ngc:,&
                                                          1-self%field%grid%ngc:,&
                                                          1-self%field%grid%ngc:,&
                                                          1:)          !< Field component to be updated.
   integer(I4P),            intent(in), optional :: step               !< Step to be perfordmed in asyncronous comp.
   logical                                       :: do_step(3)         !< Steps performed in async comp.
   integer(I4P)                                  :: p                  !< Counter.
   integer(I4P)                                  :: ptr_start, ptr_end !< Counter.
   integer(I4P)                                  :: n_recv, n_send     !< Counter.

   associate(procs_number=>self%mpih%procs_number,                                 &
             error=>self%mpih%error,                                               &
             req_send_recv=>self%mpih%req_send_recv,                               &
             comm_map_send_ptr_ghost=>self%maps%maps%comm_map_send_ptr_ghost,      &
             comm_map_recv_ptr_ghost=>self%maps%maps%comm_map_recv_ptr_ghost,      &
             recv_buffer_ghost_gpu=>self%maps%recv_buffer_ghost_gpu,               &
             send_buffer_ghost_gpu=>self%maps%send_buffer_ghost_gpu,               &
             ngc=>self%field%grid%ngc, q_gpu=>q_gpu)
   do_step = .true.
   if (present(step)) then
      do_step = .false.
      do_step(step) = .true.
   endif

   if (do_step(1)) then
      req_send_recv = MPI_REQUEST_NULL
      call populate_send_buffer_ghost_gpu_dev(ngc=ngc,                                                             &
                                              comm_map_send_ghost_cell_gpu=self%maps%comm_map_send_ghost_cell_gpu, &
                                              send_buffer_ghost_gpu=self%maps%send_buffer_ghost_gpu,               &
                                              q_gpu=q_gpu)
   endif

   if (do_step(2)) then
      ! receive
      do p=0, procs_number - 1_I4P
         ptr_start = comm_map_recv_ptr_ghost(p) + 1
         ptr_end   = comm_map_recv_ptr_ghost(p+1)
         n_recv    = ptr_end - ptr_start + 1
         if (n_recv > 0) then
            call MPI_IRECV(recv_buffer_ghost_gpu(ptr_start), n_recv, MPI_REAL8, p, 100, MPI_COMM_WORLD, &
                           req_send_recv(p), error)
         endif
      enddo
      ! send
      do p=0, procs_number - 1_I4P
         ptr_start = comm_map_send_ptr_ghost(p) + 1
         ptr_end   = comm_map_send_ptr_ghost(p+1)
         n_send    = ptr_end - ptr_start + 1
         if (n_send > 0) then
            call MPI_ISEND(send_buffer_ghost_gpu(ptr_start), n_send, MPI_REAL8, p, 100, MPI_COMM_WORLD, &
                           req_send_recv(p+procs_number), error)
         endif
      enddo
   endif

   if (do_step(3)) then
      call MPI_WAITALL(procs_number * 2, req_send_recv, MPI_STATUSES_IGNORE, error)
      call MPI_Barrier(MPI_COMM_WORLD, error)
      !RIMETTERE SENZA
      call receive_recv_buffer_ghost_gpu_dev(ngc=ngc,                                                             &
                                             comm_map_recv_ghost_cell_gpu=self%maps%comm_map_recv_ghost_cell_gpu, &
                                             recv_buffer_ghost_gpu=self%maps%recv_buffer_ghost_gpu,               &
                                             q_gpu=q_gpu)
   endif
   call MPI_Barrier(MPI_COMM_WORLD, error)
   endassociate
   endsubroutine update_ghost_mpi_gpu
endmodule adam_fnl_field_object
