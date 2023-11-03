!< ADAM, field class definition, NVF backend.
module adam_field_nvf_object
!< ADAM, field class definition, NVF backend.

use adam_common_library
use adam_field_nvf_kernels
use adam_maps_nvf_object
use adam_mpih_nvf_object
use adam_memory_nvf_library
use penf
use mpi
use cudafor

implicit none
save
private
public :: field_nvf_object

type :: field_nvf_object
   !< Field class, NVF backend.
   type(mpih_nvf_object)       :: mpih           !< MPI handler.
   type(maps_nvf_object)       :: maps           !< Maps handler.
   type(field_object), pointer :: field=>null()  !< The field.
   real(R8P), allocatable      :: q_t(:,:,:,:,:) !< Transposed cell centered variables on CPU.
   ! GPU data
   real(R8P),    allocatable, device :: q_t_gpu(:,:,:,:,:)   !< Transposed cell centered variables on GPU.
   integer(I4P), allocatable, device :: fec_1_6_array_gpu(:) !< Mapping fec1-26 to fec1-6 for boundaries (GPU).
   ! GPU data copied from field object
   real(R8P),    allocatable, device :: x_cell_gpu(:,:) !< Cells x coordinates on GPU.
   real(R8P),    allocatable, device :: y_cell_gpu(:,:) !< Cells y coordinates on GPU.
   real(R8P),    allocatable, device :: z_cell_gpu(:,:) !< Cells z coordinates on GPU.
   real(R8P),    allocatable, device :: dxyz_gpu(:,:)   !< Delta cells GPU.
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu           !< Copy data from (field_object) CPU to (field_nvf_object) GPU.
      procedure, pass(self) :: copy_transpose_cpu_gpu !< Transpose data from GPU to CPU.
      procedure, pass(self) :: copy_transpose_gpu_cpu !< Transpose data from GPU to CPU.
      procedure, pass(self) :: initialize             !< Initialize base backend.
      procedure, pass(self) :: update_ghost_local_gpu !< Update ghosts locally.
      procedure, pass(self) :: update_ghost_mpi_gpu   !< Update ghosts MPI.
endtype field_nvf_object

contains
   ! public methods
   subroutine copy_cpu_gpu(self, verbose)
   !< Copy data from (field_object) CPU to (field_nvf_object) GPU.
   class(field_nvf_object), intent(inout)        :: self     !< The base backend.
   logical,                 intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                       :: verbose_ !< Flag to activate verbose mode, local var.
   character(:), allocatable                     :: r        !< My rank stringified.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) call self%mpih%print_message('field_nvf_object%copy_cpu_gpu start')
   r = self%mpih%myrankstr
   call self%maps%copy_cpu_gpu
   call assign_allocatable_gpu(lhs=self%x_cell_gpu, rhs=self%field%x_cell, transposed=.true., &
                               msg=r//'field_nvf_object%copy_cpu_gpu(x_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%y_cell_gpu, rhs=self%field%y_cell, transposed=.true., &
                               msg=r//'field_nvf_object%copy_cpu_gpu(y_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%z_cell_gpu, rhs=self%field%z_cell, transposed=.true., &
                               msg=r//'field_nvf_object%copy_cpu_gpu(z_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=self%dxyz_gpu, rhs=self%field%dxyz, transposed=.true., &
                               msg=r//'field_nvf_object%copy_cpu_gpu(dxyz_gpu) ', verbose=verbose)
   if (verbose_) call self%mpih%print_message('field_nvf_object%copy_cpu_gpu finish')
   endsubroutine copy_cpu_gpu

   subroutine copy_transpose_cpu_gpu(self, nv, q_cpu, q_gpu)
   !< Copy transposed data from CPU to GPU.
   !< This routine is called by equation typically passing either q_gpu or q_aux_gpu.
   class(field_nvf_object), intent(inout)       :: self          !< The equation.
   integer(I4P),            intent(in)          :: nv            !< Number of varibales.
   real(R8P),               intent(in)          :: q_cpu(1:,                    &
                                                         1-self%field%grid%ngc:,&
                                                         1-self%field%grid%ngc:,&
                                                         1-self%field%grid%ngc:,&
                                                         1:)     !< Conservative variables on CPU.
   real(R8P),               intent(out), device :: q_gpu(1:,                    &
                                                        1-self%field%grid%ngc:,&
                                                        1-self%field%grid%ngc:,&
                                                        1-self%field%grid%ngc:,&
                                                        1:)      !< Conservative variables on GPU.
   integer(I4P)                                 :: i, j, k, b, v !< Counter.

   associate(blocks_number=>self%field%blocks_number, &
             ni=>self%field%grid%ni,                  &
             nj=>self%field%grid%nj,                  &
             nk=>self%field%grid%nk,                  &
             ngc=>self%field%grid%ngc,                &
             q_t=>self%q_t)
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
      q_gpu = q_t
   endassociate
   endsubroutine copy_transpose_cpu_gpu

   subroutine copy_transpose_gpu_cpu(self, nv, q_gpu, q_cpu)
   !< Copy transposed data from GPU to CPU.
   !< This routine is called by equation typically passing either q_gpu or q_aux_gpu.
   class(field_nvf_object), intent(inout)      :: self      !< The equation.
   integer(I4P),            intent(in)         :: nv        !< Number of varibales.
   real(R8P),               intent(in), device :: q_gpu(1:,                    &
                                                        1-self%field%grid%ngc:,&
                                                        1-self%field%grid%ngc:,&
                                                        1-self%field%grid%ngc:,&
                                                        1:) !< Conservative variables on GPU.
   real(R8P),               intent(inout)      :: q_cpu(1:,                    &
                                                        1-self%field%grid%ngc:,&
                                                        1-self%field%grid%ngc:,&
                                                        1-self%field%grid%ngc:,&
                                                        1:) !< Conservative variables on CPU.

   associate(blocks_number=>self%field%blocks_number, &
             ni=>self%field%grid%ni,                  &
             nj=>self%field%grid%nj,                  &
             nk=>self%field%grid%nk,                  &
             ngc=>self%field%grid%ngc,                &
             q_t_gpu=>self%q_t_gpu)

      call self%mpih%check_cuda_error(error_code=-15, msg='field_nvf_object%copy_transpose_gpu_cpu start')
      call copy_transpose_gpu_cpu_cuf(ni=ni,nj=nj,nk=nk,ngc=ngc,nv=nv,blocks_number=blocks_number,q_gpu=q_gpu,q_t_gpu=q_t_gpu)
      call self%mpih%check_cuda_error(error_code=-15, msg='field_nvf_object%copy_transpose_gpu_cpu finish')
      ! q_t_gpu has nv_aux variables which can be larger than local nv (i.e., nv or nv_aux)
      ! q_cpu   has local nv variables which is lower than nv_aux
      q_cpu(1:nv,:,:,:,1:blocks_number) = q_t_gpu(1:nv,:,:,:,1:blocks_number)
   endassociate
   endsubroutine copy_transpose_gpu_cpu

   subroutine initialize(self, field, nv_aux, verbose)
   !< Initialize base backend.
   class(field_nvf_object), intent(inout)        :: self    !< The base backend.
   type(field_object),      intent(in), target   :: field   !< Field variable array.
   integer(I4P),            intent(in), optional :: nv_aux  !< Number of auxiliary variables.
   logical,                 intent(in), optional :: verbose !< Flag to activate verbose mode.
   integer(I4P)                                  :: nv_aux_ !< Number of auxiliary variables (local var).

   call self%mpih%initialize
   call self%mpih%print_message('field_nvf_object%initialize start')
   self%field => field
   call self%maps%initialize(maps=field%maps)
   call assign_allocatable_gpu(lhs=self%fec_1_6_array_gpu, &
                               rhsa=    FEC_1_6_ARRAY,     &
                               msg=self%mpih%myrankstr//'field_nvf_object%initialize(fec_1_6_array_gpu) ', verbose=verbose)
   nv_aux_ = self%field%nv ; if (present(nv_aux)) nv_aux_ = max(nv_aux_, nv_aux)
   call alloc_var_cpu(var=self%q_t,                                              &
                      ulb=reshape([1,field%nb,                                   &
                                   1-field%grid%ngc,field%grid%ni+field%grid%ngc,&
                                   1-field%grid%ngc,field%grid%nj+field%grid%ngc,&
                                   1-field%grid%ngc,field%grid%nk+field%grid%ngc,&
                                   1,nv_aux_],[2,5]),                            &
                      msg=self%mpih%myrankstr//'field_nvf_object%initialize alloc_var_cpu(q_t) ', verbose=verbose)
   call alloc_var_gpu(var=self%q_t_gpu,                                          &
                      ulb=reshape([1,nv_aux_,                                    &
                                   1-field%grid%ngc,field%grid%ni+field%grid%ngc,&
                                   1-field%grid%ngc,field%grid%nj+field%grid%ngc,&
                                   1-field%grid%ngc,field%grid%nk+field%grid%ngc,&
                                   1,field%nb],[2,5]),                           &
                      msg=self%mpih%myrankstr//'field_nvf_object%initialize alloc_var_gpu(q_t) ', verbose=verbose)
   call self%copy_cpu_gpu
   call self%mpih%print_message('field_nvf_object%initialize finish')
   endsubroutine initialize

   subroutine update_ghost_local_gpu(self, q_gpu)
   !< Update (local) ghost cells.
   class(field_nvf_object), intent(in)            :: self      !< The base backend.
   real(R8P),               intent(inout), device :: q_gpu(1:,                    &
                                                           1-self%field%grid%ngc:,&
                                                           1-self%field%grid%ngc:,&
                                                           1-self%field%grid%ngc:,&
                                                           1:) !< Field component to be updated.
   call self%mpih%check_cuda_error(error_code=-15, msg='field_nvf_object%update_ghost_local_gpu start')
   call update_ghost_local_gpu_cuf(local_map_ghost_cell_gpu=self%maps%local_map_ghost_cell_gpu,ngc=self%field%grid%ngc,q_gpu=q_gpu)
   call self%mpih%check_cuda_error(error_code=-15, msg='field_nvf_object%update_ghost_local_gpu finish')
   endsubroutine update_ghost_local_gpu

   subroutine update_ghost_mpi_gpu(self, q_gpu, step)
   !< Update ghost cells within other processes.
   class(field_nvf_object), intent(inout)         :: self      !< The base backend.
   real(R8P),               intent(inout), device :: q_gpu(1:,                    &
                                                           1-self%field%grid%ngc:,&
                                                           1-self%field%grid%ngc:,&
                                                           1-self%field%grid%ngc:,&
                                                           1:) !< Field component to be updated.
   integer(I4P),            intent(in), optional  :: step      !< Step to be perfordmed in asyncronous comp.
   logical                                        :: do_step(3)!< Steps performed in async comp.
   integer(I4P)                                   :: p                                      !< Counter.
   integer(I4P)                                   :: ptr_start, ptr_end                     !< Counter.
   integer(I4P)                                   :: n_recv, n_send                         !< Counter.

   call self%mpih%check_cuda_error(error_code=-15, msg='field_nvf_object%update_ghost_mpi_gpu start')
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
      call populate_send_buffer_ghost_gpu_cuf(ngc=ngc,                                                             &
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
      call receive_recv_buffer_ghost_gpu_cuf(ngc=ngc,                                                             &
                                             comm_map_recv_ghost_cell_gpu=self%maps%comm_map_recv_ghost_cell_gpu, &
                                             recv_buffer_ghost_gpu=self%maps%recv_buffer_ghost_gpu,               &
                                             q_gpu=q_gpu)
   endif
   call MPI_Barrier(MPI_COMM_WORLD, error)
   endassociate
   call self%mpih%check_cuda_error(error_code=-15, msg='field_nvf_object%update_ghost_mpi_gpu finish')
   endsubroutine update_ghost_mpi_gpu
endmodule adam_field_nvf_object
