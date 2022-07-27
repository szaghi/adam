!< ADAM, base GPU class definition.
module adam_base_gpu_object
!< ADAM, base GPU class definition: provide methods for GPU backend handling.

use adam_field_object, only : field_object
use adam_mpih_object,  only : mpih_object
use adam_tree_object,  only : tree_object
use adam_memory_gpu_lib
use adam_parameters
use penf
use mpi
use cudafor
use, intrinsic :: iso_fortran_env, only : stderr=>error_unit

implicit none
save
private
public :: base_gpu_object

type :: base_gpu_object
   !< Base GPU class definition.
   !<
   !< Provide methods for GPU backend.
   type(mpih_object)           :: mpih           !< MPI handler.
   type(tree_object),  pointer :: tree=>null()   !< The tree.
   type(field_object), pointer :: field=>null()  !< The field.
   real(R8P), allocatable      :: q_t(:,:,:,:,:) !< Transposed cell centered variables on CPU.
   ! MPI data
   integer(I4P)              :: mydev=0_I4P       !< My GPU rank.
   integer(I4P)              :: local_comm=0_I4P  !< Local communicator.
   integer(I4P), allocatable :: req_send_recv(:)  !< MPI request receive flags.
   ! GPU data
   real(R8P)                         :: memory_avail=0._R8P                  !< Device memory available (Gb).
   real(R8P),    allocatable, device :: q_t_gpu(:,:,:,:,:)                   !< Transposed cell centered variables on GPU.
   integer(I8P), allocatable, device :: local_map_ghost_cell_gpu(:,:)        !< Local map for ghost cells updating, cells order.
   integer(I8P), allocatable, device :: local_map_ghost_fluxes_cell_gpu(:,:) !< Local map for ghost fluxes updating, cells order.
   integer(I8P), allocatable, device :: comm_map_recv_ghost_cell_gpu(:,:)    !< Communication map, `fec` information, cell order.
   integer(I8P), allocatable, device :: comm_map_send_ghost_cell_gpu(:,:)    !< Communication map, `fec` information, cell order.
   real(R8P),    allocatable, device :: send_buffer_ghost_gpu(:)             !< Send buffer of ghost cells.
   real(R8P),    allocatable, device :: recv_buffer_ghost_gpu(:)             !< Receive buffer of ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_crown_gpu(:,:,:)        !< Local map for face BC ghost cells, "crown" order.
   integer(I4P), allocatable, device :: fec_1_6_array_gpu(:)                 !< Mapping fec1-26 to fec1-6 for boundaries (GPU).
   real(R8P),    allocatable, device :: x_cell_gpu(:,:)                      !< Cells x coordinates on GPU.
   real(R8P),    allocatable, device :: y_cell_gpu(:,:)                      !< Cells y coordinates on GPU.
   real(R8P),    allocatable, device :: z_cell_gpu(:,:)                      !< Cells z coordinates on GPU.
   real(R8P),    allocatable, device :: dxyz_gpu(:,:)                        !< Delta cells GPU.
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu           !< Copy data from (field) CPU to (base_gpu) GPU.
      procedure, pass(self) :: copy_transpose_cpu_gpu !< Transpose data from GPU to CPU.
      procedure, pass(self) :: copy_transpose_gpu_cpu !< Transpose data from GPU to CPU.
      procedure, pass(self) :: initialize             !< Initialize base backend.
      procedure, pass(self) :: initialize_gpu         !< Initialize GPU main data.
      procedure, pass(self) :: update_ghost_local_gpu !< Update ghosts locally.
      procedure, pass(self) :: update_ghost_mpi_gpu   !< Update ghosts MPI.
endtype base_gpu_object

contains
   ! public methods
   subroutine copy_cpu_gpu(self, verbose)
   !< Copy data from (field) CPU to (base_gpu) GPU.
   class(base_gpu_object), intent(inout)        :: self     !< The base backend.
   logical,                intent(in), optional :: verbose  !< Flag to activate verbose mode.
   logical                                      :: verbose_ !< Flag to activate verbose mode, local var.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) print '(A)', self%mpih%myrankstr//'base_gpu%copy_cpu_gpu start'
   call assign_allocatable_gpu(lhs=     self%send_buffer_ghost_gpu, &
                               rhs=self%tree%send_buffer_ghost,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(send_buffer_ghost_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=     self%recv_buffer_ghost_gpu, &
                               rhs=self%tree%recv_buffer_ghost,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(recv_buffer_ghost_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%x_cell_gpu, &
                               rhs=self%field%x_cell,     &
                               transposed=.true., msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(x_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%y_cell_gpu, &
                               rhs=self%field%y_cell,     &
                               transposed=.true., msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(y_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%z_cell_gpu, &
                               rhs=self%field%z_cell,     &
                               transposed=.true., msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(z_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%dxyz_gpu, &
                               rhs=self%field%dxyz,     &
                               transposed=.true., msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(dxyz_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=     self%local_map_ghost_cell_gpu, &
                               rhs=self%tree%local_map_ghost_cell,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(local_map_ghost_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=  self%comm_map_send_ghost_cell_gpu, &
                               rhs=self%tree%comm_map_send_ghost_cell, &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(comm_map_send_ghost_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=     self%comm_map_recv_ghost_cell_gpu, &
                               rhs=self%tree%comm_map_recv_ghost_cell,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(comm_map_recv_ghost_cell_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=     self%local_map_bc_crown_gpu, &
                               rhs=self%tree%local_map_bc_crown,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(local_map_bc_crown_gpu) ', verbose=verbose)
   if (verbose_) print '(A)', self%mpih%myrankstr//'base_gpu%copy_cpu_gpu finish'
   endsubroutine copy_cpu_gpu

   subroutine copy_transpose_cpu_gpu(self, nv, q_cpu, q_gpu)
   !< Copy transposed data from CPU to GPU.
   !< This routine is called by equation typically passing either q_gpu or q_aux_gpu.
   class(base_gpu_object), intent(inout)       :: self          !< The equation.
   integer(I4P),           intent(in)          :: nv            !< Number of varibales.
   real(R8P),              intent(in)          :: q_cpu(1:,                    &
                                                        1-self%field%grid%ngc:,&
                                                        1-self%field%grid%ngc:,&
                                                        1-self%field%grid%ngc:,&
                                                        1:)     !< Conservative variables on CPU.
   real(R8P),              intent(out), device :: q_gpu(1:,                    &
                                                       1-self%field%grid%ngc:,&
                                                       1-self%field%grid%ngc:,&
                                                       1-self%field%grid%ngc:,&
                                                       1:)      !< Conservative variables on GPU.
   integer(I4P)                                :: i, j, k, b, v !< Counter.

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
   class(base_gpu_object), intent(inout)      :: self      !< The equation.
   integer(I4P),           intent(in)         :: nv        !< Number of varibales.
   real(R8P),              intent(in), device :: q_gpu(1:,                    &
                                                       1-self%field%grid%ngc:,&
                                                       1-self%field%grid%ngc:,&
                                                       1-self%field%grid%ngc:,&
                                                       1:) !< Conservative variables on GPU.
   real(R8P),              intent(out)        :: q_cpu(1:,                    &
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
      call copy_transpose_gpu_cpu_cuf(ni=ni, nj=nj, nk=nk, ngc=ngc, nv=nv, blocks_number=blocks_number, &
                                      q_gpu=q_gpu, q_t_gpu=q_t_gpu, q_cpu=q_cpu)
   endassociate
   endsubroutine copy_transpose_gpu_cpu

   subroutine initialize(self, tree, field, nv_aux, verbose)
   !< Initialize base backend.
   class(base_gpu_object), intent(inout)        :: self             !< The base backend.
   type(tree_object),      intent(in), target   :: tree             !< The tree.
   type(field_object),     intent(in), target   :: field            !< The field.
   integer(I4P),           intent(in), optional :: nv_aux           !< Number of auxiliary variables.
   logical,                intent(in), optional :: verbose          !< Flag to activate verbose mode.
   integer(I4P)                                 :: nv_aux_          !< Number of auxiliary variables (local var).
   integer(I4P), allocatable                    :: fec_1_6_array(:) !< Mapping fec1-26 to fec1-6 for boundaries.

   print '(A)', self%mpih%myrankstr//'base_gpu%initialize start'
   allocate(fec_1_6_array(26))
   self%tree  => tree
   self%field => field
   allocate(self%req_send_recv(0:self%mpih%procs_number*2-1))
   fec_1_6_array([1,7,9,11,13,19,21,23,25])  = 1
   fec_1_6_array([2,8,10,12,14,20,22,24,26]) = 2
   fec_1_6_array([3,15,17])                  = 3
   fec_1_6_array([4,16,18])                  = 4
   fec_1_6_array([5])                        = 5
   fec_1_6_array([6])                        = 6
   call assign_allocatable_gpu(lhs=self%fec_1_6_array_gpu, &
                               rhs=     fec_1_6_array,     &
                               msg=self%mpih%myrankstr//'base_gpu%initialize(fec_1_6_array_gpu)', verbose=verbose)
   nv_aux_ = self%field%nv
   if (present(nv_aux)) nv_aux_ = max(nv_aux_, nv_aux)
   ! allocate buffers for copy-transposes performed by equation
   allocate(self%q_t(1:field%nb,                                    &
                     1-field%grid%ngc:field%grid%ni+field%grid%ngc, &
                     1-field%grid%ngc:field%grid%nj+field%grid%ngc, &
                     1-field%grid%ngc:field%grid%nk+field%grid%ngc, 1:nv_aux_))
   ! call alloc_var_gpu(var=self%q_t_gpu,                                          &
   !                    ulb=reshape([1,nv_aux_,                                    &
   !                                 1-field%grid%ngc,field%grid%ni+field%grid%ngc,&
   !                                 1-field%grid%ngc,field%grid%nj+field%grid%ngc,&
   !                                 1-field%grid%ngc,field%grid%nk+field%grid%ngc,&
   !                                 1,field%nb],[2,5]),                           &
   !                    msg=self%mpih%myrankstr//'base_gpu%alloc(q_t_gpu) ', verbose=verbose)
   allocate(self%q_t_gpu(1:nv_aux_,                                     &
                         1-field%grid%ngc:field%grid%ni+field%grid%ngc, &
                         1-field%grid%ngc:field%grid%nj+field%grid%ngc, &
                         1-field%grid%ngc:field%grid%nk+field%grid%ngc, 1:field%nb))
   ! copy CPU-to-GPU of base_gpu variables (maps and cells, not q_gpu)
   call self%copy_cpu_gpu
   deallocate(fec_1_6_array)
   print '(A)', self%mpih%myrankstr//'base_gpu%initialize finish'
   endsubroutine initialize

   subroutine initialize_gpu(self)
   !< Initialize GPU main data.
   !< @Note This must be the first routine called.
   class(base_gpu_object), intent(inout) :: self              !< The base backend.
   type(cudadeviceprop)                  :: device_properties !< Device properties.

   call self%mpih%initialize
   print '(A)', self%mpih%myrankstr//'base_gpu%initialize_gpu start'
   call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%local_comm, self%mpih%error)
   call MPI_COMM_RANK(self%local_comm, self%mydev, self%mpih%error)
   self%mpih%error = CudaSetDevice(self%mydev)
   self%mpih%error = cudaGetDeviceProperties(device_properties, self%mydev)
   self%memory_avail = real(device_properties%totalGlobalMem, R8P)/1e9
   call print_device_properties(self, device_properties)
   print '(A)', self%mpih%myrankstr//'base_gpu%initialize_gpu finish'
   endsubroutine initialize_gpu

   subroutine print_device_properties(self, device_properties)
   !< Pretty print device properties.
   class(base_gpu_object), intent(in) :: self               !< The base backend.
   type(cudadeviceprop),   intent(in) :: device_properties  !< Device properties.

   associate(r=>self%mpih%myrankstr)
   print'(A)',r//"total global memory:         "//trim(str(real(device_properties%totalGlobalMem)/1e9           ,.true.))//" Gbytes"
   print'(A)',r//"shared mem per block:        "//trim(str(     device_properties%sharedMemPerBlock             ,.true.))//" bytes"
   print'(A)',r//"regs per block:              "//trim(str(     device_properties%regsPerBlock                  ,.true.))
   print'(A)',r//"warp size:                   "//trim(str(     device_properties%warpSize                      ,.true.))
   print'(A)',r//"max threads per block:       "//trim(str(     device_properties%maxThreadsPerBlock            ,.true.))
   print'(A)',r//"max threads dim:             "//trim(str(     device_properties%maxThreadsDim                 ,.true.))
   print'(A)',r//"clock rate:                  "//trim(str(real(device_properties%clockRate)/1e6                ,.true.))//" GHz"
   print'(A)',r//"total const memory:          "//trim(str(     device_properties%totalConstMem                 ,.true.))//" bytes"
   print'(A)',r//"compute capability revision: "//trim(str(    [device_properties%major,device_properties%minor],.true.))
   print'(A)',r//"multi processor count:       "//trim(str(     device_properties%multiProcessorCount           ,.true.))
   print'(A)',r//"L2 cache size:               "//trim(str(     device_properties%l2CacheSize                   ,.true.))
   print'(A)',r//"max threads per SMP:         "//trim(str(     device_properties%maxThreadsPerMultiProcessor   ,.true.))
   print'(A)',r//"device rank:                 "//trim(str(     self%mydev                                      ,.true.))
   endassociate
   endsubroutine print_device_properties

   subroutine update_ghost_local_gpu(self, q_gpu)
   !< Update (local) ghost cells.
   class(base_gpu_object), intent(in)            :: self      !< The base backend.
   real(R8P),              intent(inout), device :: q_gpu(1:,                    &
                                                          1-self%field%grid%ngc:,&
                                                          1-self%field%grid%ngc:,&
                                                          1-self%field%grid%ngc:,&
                                                          1:) !< Field component to be updated.
   call update_ghost_local_gpu_cuf(local_map_ghost_cell_gpu=self%local_map_ghost_cell_gpu, ngc=self%field%grid%ngc, q_gpu=q_gpu)
   endsubroutine update_ghost_local_gpu

   subroutine update_ghost_mpi_gpu(self, q_gpu, step)
   !< Update ghost cells within other processes.
   class(base_gpu_object), intent(inout)         :: self      !< The base backend.
   real(R8P),              intent(inout), device :: q_gpu(1:,                    &
                                                          1-self%field%grid%ngc:,&
                                                          1-self%field%grid%ngc:,&
                                                          1-self%field%grid%ngc:,&
                                                          1:) !< Field component to be updated.
   integer(I4P),           intent(in), optional  :: step      !< Step to be perfordmed in asyncronous comp.

   call update_ghost_mpi_gpu_cuf(procs_number=self%mpih%procs_number,                            &
                                 req_send_recv=self%field%req_send_recv,                         &
                                 ! comm_map_send_ptr_ghost=self%field%comm_map_send_ptr_ghost,     &
                                 ! comm_map_recv_ptr_ghost=self%field%comm_map_recv_ptr_ghost,     &
                                 comm_map_send_ptr_ghost=self%tree%comm_map_send_ptr_ghost,      &
                                 comm_map_recv_ptr_ghost=self%tree%comm_map_recv_ptr_ghost,      &
                                 comm_map_send_ghost_cell_gpu=self%comm_map_send_ghost_cell_gpu, &
                                 comm_map_recv_ghost_cell_gpu=self%comm_map_recv_ghost_cell_gpu, &
                                 recv_buffer_ghost_gpu=self%recv_buffer_ghost_gpu,               &
                                 send_buffer_ghost_gpu=self%send_buffer_ghost_gpu,               &
                                 ngc=self%field%grid%ngc, q_gpu=q_gpu, step=step)
   endsubroutine update_ghost_mpi_gpu

   ! non TBP CUF methods
   subroutine copy_transpose_gpu_cpu_cuf(ni, nj, nk, ngc, nv, blocks_number, q_gpu, q_t_gpu, q_cpu)
   !< Copy transposed data from GPU to CPU by CUF threads.
   integer(I4P), intent(in)            :: ni            !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj            !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk            !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc           !< Ghost cells number.
   integer(I4P), intent(in)            :: nv            !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number !< Number of blocks.
   real(R8P),    intent(in), device    :: q_gpu(1:,    &
                                                1-ngc:,&
                                                1-ngc:,&
                                                1-ngc:,&
                                                1:)     !< Conservative variables on GPU.
   real(R8P),    intent(inout), device :: q_t_gpu(1:,    &
                                                  1-ngc:,&
                                                  1-ngc:,&
                                                  1-ngc:,&
                                                  1:)   !< Conservative (transposed) variables on GPU.
   real(R8P),    intent(out)           :: q_cpu(1:,    &
                                                1-ngc:,&
                                                1-ngc:,&
                                                1-ngc:,&
                                                1:)     !< Conservative variables on CPU.
   integer(I4P)                        :: i, j, k, b, v !< Counter.
   integer(I4P)                        :: error         !< Error traping flag.
   integer(I4P)                        :: iercuda       !< Error trapping flag for CUDAFortran.

   iercuda = cudaGetLastError()
   if (iercuda /= cudaSuccess) then
      write(stderr, '(A)') 'base_gpu%copy_transpose_gpu_cpu_cuf start, cuda ERROR: '//cudaGetErrorString(iercuda)
      call MPI_Abort(MPI_COMM_WORLD, -15, error)
      stop
   endif

   !$cuf kernel do(4) <<<*,*>>>
   do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
         do i=1-ngc, ni+ngc
            do b=1, blocks_number
               do v=1, nv
                  q_t_gpu(v,i,j,k,b) = q_gpu(b,i,j,k,v)
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   ! q_t_gpu has nv_aux variables which can be larger than local nv (i.e., nv or nv_aux)
   ! q_cpu   has local nv variables which is lower than nv_aux
   q_cpu(1:nv,:,:,:,1:blocks_number) = q_t_gpu(1:nv,:,:,:,1:blocks_number)

   iercuda = cudaGetLastError()
   if (iercuda /= cudaSuccess) then
      write(stderr, '(A)') 'base_gpu%copy_transpose_gpu_cpu_cuf finish, cuda ERROR: '//cudaGetErrorString(iercuda)
      call MPI_Abort(MPI_COMM_WORLD, -15, error)
      stop
   endif
   endsubroutine copy_transpose_gpu_cpu_cuf

   subroutine update_ghost_local_gpu_cuf(ngc, local_map_ghost_cell_gpu, q_gpu)
   !< Update (local) ghost cells.
   integer(I4P), intent(in)                         :: ngc                           !< Ghost cells number.
   integer(I8P), intent(in),    device, allocatable :: local_map_ghost_cell_gpu(:,:) !< Local map of ghost cells.
   real(R8P),    intent(inout), device              :: q_gpu(1:,    &
                                                             1-ngc:,&
                                                             1-ngc:,&
                                                             1-ngc:,1:)              !< Field component to be updated.
   integer(I4P)                                     :: ic, jc, kc, mf, v             !< Counter.
   integer(I4P)                                     :: b_recv                        !< Index of receiving block.
   integer(I4P)                                     :: b_send                        !< Index of sending block.
   integer(I4P)                                     :: i_recv                        !< I recv index.
   integer(I4P)                                     :: j_recv                        !< J recv index.
   integer(I4P)                                     :: k_recv                        !< K recv index.
   integer(I4P)                                     :: i_send                        !< I send index.
   integer(I4P)                                     :: j_send                        !< J send index.
   integer(I4P)                                     :: k_send                        !< K send index.
   integer(I4P)                                     :: one_or_eight                  !< Flag triggering 8 cells mean.
   integer(I4P)                                     :: error                         !< Error traping flag.
   integer(I4P)                                     :: iercuda                       !< Error trapping flag for CUDAFortran.

   iercuda = cudaGetLastError()
   if (iercuda /= cudaSuccess) then
      write(stderr, '(A)') 'base_gpu%update_ghost_local_gpu_cuf start, cuda ERROR: '//cudaGetErrorString(iercuda)
      call MPI_Abort(MPI_COMM_WORLD, -15, error)
      stop
   endif

   if (.not.allocated(local_map_ghost_cell_gpu)) return
   !$cuf kernel do(2) <<<*,*>>>
   do v=1, size(q_gpu, dim=5)
      do mf=1, size(local_map_ghost_cell_gpu, dim=1)
         b_send       = local_map_ghost_cell_gpu(mf,1)
         b_recv       = local_map_ghost_cell_gpu(mf,2)
         i_send       = local_map_ghost_cell_gpu(mf,3)
         j_send       = local_map_ghost_cell_gpu(mf,4)
         k_send       = local_map_ghost_cell_gpu(mf,5)
         i_recv       = local_map_ghost_cell_gpu(mf,6)
         j_recv       = local_map_ghost_cell_gpu(mf,7)
         k_recv       = local_map_ghost_cell_gpu(mf,8)
         one_or_eight = local_map_ghost_cell_gpu(mf,9)
         if (one_or_eight==1) then
            q_gpu(b_recv,i_recv, j_recv, k_recv,v) = q_gpu(b_send, i_send,j_send,k_send,v)
         else
            q_gpu(b_recv,i_recv,j_recv,k_recv,v) = 0._R8P
            do kc=0,1 ; do jc=0,1 ; do ic=0,1
               q_gpu(b_recv,i_recv,j_recv,k_recv,v) = q_gpu(b_recv,i_recv,   j_recv,   k_recv,   v) + &
                                                      q_gpu(b_send,i_send+ic,j_send+jc,k_send+kc,v)
            enddo ; enddo ; enddo
            q_gpu(b_recv,i_recv,j_recv,k_recv,v) = q_gpu(b_recv,i_recv,j_recv,k_recv,v) * 0.125_R8P
         endif
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   iercuda = cudaGetLastError()
   if (iercuda /= cudaSuccess) then
      write(stderr, '(A)') 'base_gpu%update_ghost_local_gpu_cuf finish, cuda ERROR: '//cudaGetErrorString(iercuda)
      call MPI_Abort(MPI_COMM_WORLD, -15, error)
      stop
   endif
   endsubroutine update_ghost_local_gpu_cuf

   subroutine update_ghost_mpi_gpu_cuf(ngc, procs_number, req_send_recv,                           &
                                       comm_map_send_ptr_ghost, comm_map_recv_ptr_ghost,           &
                                       comm_map_recv_ghost_cell_gpu, comm_map_send_ghost_cell_gpu, &
                                       recv_buffer_ghost_gpu, send_buffer_ghost_gpu , q_gpu, step)
   !< Update ghost cells within other processes.
   integer(I4P), intent(in)                         :: ngc                                    !< Ghost cells number.
   integer(I4P),              intent(in)            :: procs_number                           !< Number of MPI processes.
   integer(I4P), allocatable, intent(inout)         :: req_send_recv(:)                       !< MPI request receive flags.
   integer(I4P), allocatable, intent(in)            :: comm_map_send_ptr_ghost(:)             !< Comm map, pntrs list to send.
   integer(I4P), allocatable, intent(in)            :: comm_map_recv_ptr_ghost(:)             !< Comm map, pntrs list to recv.
   integer(I8P), allocatable, intent(in),    device :: comm_map_recv_ghost_cell_gpu(:,:)      !< Comm map, cell information.
   integer(I8P), allocatable, intent(in),    device :: comm_map_send_ghost_cell_gpu(:,:)      !< Comm map, cell information.
   real(R8P),    allocatable, intent(inout), device :: recv_buffer_ghost_gpu(:)               !< Receive buffer of ghost cells.
   real(R8P),    allocatable, intent(inout), device :: send_buffer_ghost_gpu(:)               !< Send buffer of ghost cells.
   real(R8P),                 intent(inout), device :: q_gpu(1:,    &
                                                             1-ngc:,&
                                                             1-ngc:,&
                                                             1-ngc:,&
                                                             1:)                              !< Field component to be updated.
   integer(I4P),              intent(in), optional  :: step                                   !< Step performed in async comp.
   logical                                          :: do_step(3)                             !< Steps performed in async comp.
   integer(I4P)                                     :: i, j, k                                !< Counter.
   integer(I4P)                                     :: ic, jc, kc                             !< Counter.
   integer(I4P)                                     :: fec, mf, rf, sf, n, p, v               !< Counter.
   integer(I4P)                                     :: portion                                !< Portion of fec updated (0=>whole).
   integer(I4P)                                     :: b_send, i_send, j_send, k_send, v_send !< Send indexes.
   integer(I4P)                                     :: b_recv, i_recv, j_recv, k_recv, v_recv !< Receive indexes.
   integer(I4P)                                     :: ptr_start, ptr_end                     !< Counter.
   integer(I4P)                                     :: n_recv, n_send                         !< Counter.
   integer(I4P)                                     :: recv_rank                              !< Rank of receiving block.
   integer(I4P)                                     :: send_ptr, send_ctr                     !< Counter.
   integer(I4P)                                     :: recv_ptr, recv_ctr                     !< Counter.
   integer(I4P)                                     :: c_send, c_recv                         !< Counter.
   integer(I4P)                                     :: one_or_eight                           !< Flag triggering 8 cells mean.
   integer(I4P)                                     :: error                                  !< Error traping flag.
   integer(I4P)                                     :: iercuda                                !< Error trapping flag for CUDAFor.

   iercuda = cudaGetLastError()
   if (iercuda /= cudaSuccess) then
      write(stderr, '(A)') 'base_gpu%update_ghost_mpi_gpu_cuf start, cuda ERROR: '//cudaGetErrorString(iercuda)
      call MPI_Abort(MPI_COMM_WORLD, -15, iercuda)
      stop
   endif

   do_step = .true.
   if (present(step)) then
      do_step = .false.
      do_step(step) = .true.
   endif

   if (do_step(1)) then
      req_send_recv = MPI_REQUEST_NULL

      iercuda = cudaGetLastError()
      if (iercuda /= cudaSuccess) then
         write(stderr, '(A)') 'base_gpu%update_ghost_mpi_gpu_cuf before send buffer fill, cuda ERROR: '//cudaGetErrorString(iercuda)
         call MPI_Abort(MPI_COMM_WORLD, -15, iercuda)
         stop
      endif

      ! populate send buffer
      if(allocated(comm_map_send_ghost_cell_gpu)) then
         !$cuf kernel do(1) <<<*,*>>>
         do sf=1, size(comm_map_send_ghost_cell_gpu, dim=1)
            b_send       = comm_map_send_ghost_cell_gpu(sf,1)
            i_send       = comm_map_send_ghost_cell_gpu(sf,2)
            j_send       = comm_map_send_ghost_cell_gpu(sf,3)
            k_send       = comm_map_send_ghost_cell_gpu(sf,4)
            v_send       = comm_map_send_ghost_cell_gpu(sf,5)
            c_recv       = comm_map_send_ghost_cell_gpu(sf,6)
            one_or_eight = comm_map_send_ghost_cell_gpu(sf,7)
            if (one_or_eight==1) then
               send_buffer_ghost_gpu(c_recv) = q_gpu(b_send,i_send,j_send,k_send,v_send)
            else
               send_buffer_ghost_gpu(c_recv) = 0._R8P
               do kc=0,1 ; do jc=0,1 ; do ic=0,1
                  send_buffer_ghost_gpu(c_recv) = send_buffer_ghost_gpu(c_recv) + &
                                                  q_gpu(b_send,i_send+ic,j_send+jc,k_send+kc,v_send)
               enddo ; enddo ; enddo
               send_buffer_ghost_gpu(c_recv) = send_buffer_ghost_gpu(c_recv) * 0.125_R8P
            endif
         enddo
         !@cuf iercuda=cudaDeviceSynchronize()
      endif

      iercuda = cudaGetLastError()
      if (iercuda /= cudaSuccess) then
         write(stderr, '(A)') 'base_gpu%update_ghost_mpi_gpu_cuf after send buffer fill, cuda ERROR: '//cudaGetErrorString(iercuda)
         call MPI_Abort(MPI_COMM_WORLD, -15, error)
         stop
      endif
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

      if(allocated(comm_map_recv_ghost_cell_gpu)) then
         ! retrive from receive buffer
         !$cuf kernel do(1) <<<*,*>>>
         do rf=1, size(comm_map_recv_ghost_cell_gpu, dim=1)
            c_send = comm_map_recv_ghost_cell_gpu(rf,1)
            b_recv = comm_map_recv_ghost_cell_gpu(rf,2)
            i_recv = comm_map_recv_ghost_cell_gpu(rf,3)
            j_recv = comm_map_recv_ghost_cell_gpu(rf,4)
            k_recv = comm_map_recv_ghost_cell_gpu(rf,5)
            v_recv = comm_map_recv_ghost_cell_gpu(rf,6)
            q_gpu(b_recv,i_recv,j_recv,k_recv,v_recv) = recv_buffer_ghost_gpu(c_send)
         enddo
         !@cuf iercuda=cudaDeviceSynchronize()
      endif
   endif
   call MPI_Barrier(MPI_COMM_WORLD, error)

   iercuda = cudaGetLastError()
   if (iercuda /= cudaSuccess) then
      write(stderr, '(A)') 'base_gpu%update_ghost_mpi_gpu_cuf finish, cuda ERROR: '//cudaGetErrorString(iercuda)
      call MPI_Abort(MPI_COMM_WORLD, -15, error)
      stop
   endif
   endsubroutine update_ghost_mpi_gpu_cuf
endmodule adam_base_gpu_object
