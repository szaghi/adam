!< ADAM, base GPU class definition.
module adam_base_gpu_object
!< ADAM, base GPU class definition: provide update ghosts methods for CPU backend.

use adam_field_object
use adam_parameters
use PENF
use MPI
use CUDAFOR

implicit none
private
public :: base_gpu_object
public :: assign_allocatable_gpu

type :: base_gpu_object
   !< Equation base GPU class definition.
   !<
   !< Provide update ghosts methods for GPU backend.
   type(field_object), pointer :: field=>null() !< The field.
   ! MPI data
   integer(I4P)              :: myrank=0_I4P       !< MPI rank process.
   integer(I4P)              :: procs_number=1_I4P !< Number of MPI processes.
   integer(I4P)              :: error=0_I4P        !< Error traping flag.
   integer(I4P)              :: mydev=0_I4P        !< My GPU rank.
   integer(I4P)              :: local_comm=0_I4P   !< Local communicator.
   integer(I4P), allocatable :: req_send_recv(:)   !< MPI request receive flags.
   real(R8P),    allocatable :: q_t(:,:,:,:,:)     !< Transposed cell centered variables on CPU.
   ! GPU data
   real(R8P),    allocatable, device :: q_t_gpu(:,:,:,:,:)           !< Transposed cell centered variables on GPU.
   integer(I8P), allocatable, device :: local_map_ghost_cell_gpu(:,:)!< Local map for ghost cells updating, cells order.
   integer(I8P), allocatable, device :: local_map_ghost_gpu(:,:)     !< Local map for ghost cells updating, fecs order.
   integer(I8P), allocatable, device :: comm_map_recv_ghost_gpu(:,:) !< Communication map, `fec` information.
   integer(I8P), allocatable, device :: comm_map_send_ghost_gpu(:,:) !< Communication map, `fec` information.
   integer(I8P), allocatable, device :: comm_map_recv_ghost_cell_gpu(:,:) !< Communication map, `fec` information, cell order.
   integer(I8P), allocatable, device :: comm_map_send_ghost_cell_gpu(:,:) !< Communication map, `fec` information, cell order.
   real(R8P),    allocatable, device :: send_buffer_ghost_gpu(:)     !< Send buffer of ghost cells.
   real(R8P),    allocatable, device :: recv_buffer_ghost_gpu(:)     !< Receive buffer of ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_face_gpu(:,:)   !< Local map for face BC ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_edge_gpu(:,:)   !< Local map for edge BC ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_corner_gpu(:,:) !< Local map for corner BC ghost cells.
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu           !< Copy data from CPU to GPU.
      procedure, pass(self) :: copy_transpose_cpu_gpu !< Transpose data from GPU to CPU.
      procedure, pass(self) :: copy_transpose_gpu_cpu !< Transpose data from GPU to CPU.
      procedure, pass(self) :: create_maps_cell       !< Create mapps in cells order form the fecs ordered ones.
      procedure, pass(self) :: destroy                !< Destroy the equation.
      procedure, pass(self) :: initialize             !< Initialize the equation.
      procedure, pass(self) :: update_ghost_local_gpu !< Update ghosts locally.
      procedure, pass(self) :: update_ghost_mpi_gpu   !< Update ghosts MPI.
      ! operators
      generic :: assignment(=) => base_assign_base      !< Overload `=`.
      procedure, pass(lhs), private :: base_assign_base !< Operator `=`.
endtype base_gpu_object

interface assign_allocatable_gpu
   !< Safe assign allocatable arrays (GPU), generic interface.
   module procedure assign_allocatable_I8P_1D_gpu !< Safe assign allocatable arrays, I8P 1D type.
   module procedure assign_allocatable_I8P_2D_gpu !< Safe assign allocatable arrays, I8P 2D type.
   module procedure assign_allocatable_I4P_1D_gpu !< Safe assign allocatable arrays, I4P 1D type.
   module procedure assign_allocatable_I4P_2D_gpu !< Safe assign allocatable arrays, I4P 2D type.
   module procedure assign_allocatable_R8P_1D_gpu !< Safe assign allocatable arrays, R8P 1D type.
   module procedure assign_allocatable_R8P_2D_gpu !< Safe assign allocatable arrays, R8P 2D type.
   module procedure assign_allocatable_R8P_3D_gpu !< Safe assign allocatable arrays, R8P 3D type.
   module procedure assign_allocatable_R8P_4D_gpu !< Safe assign allocatable arrays, R8P 4D type.
   module procedure assign_allocatable_R8P_5D_gpu !< Safe assign allocatable arrays, R8P 5D type.
   module procedure assign_allocatable_R8P_6D_gpu !< Safe assign allocatable arrays, R8P 6D type.
endinterface assign_allocatable_gpu

contains
   ! public methods
   subroutine copy_cpu_gpu(self)
   !< Copy data from CPU to GPU.
   class(base_gpu_object), intent(inout) :: self !< The base backend.

   if (allocated(self%local_map_ghost_gpu    )) deallocate(self%local_map_ghost_gpu    )
   if (allocated(self%comm_map_recv_ghost_gpu)) deallocate(self%comm_map_recv_ghost_gpu)
   if (allocated(self%comm_map_send_ghost_gpu)) deallocate(self%comm_map_send_ghost_gpu)
   if (allocated(self%send_buffer_ghost_gpu  )) deallocate(self%send_buffer_ghost_gpu  )
   if (allocated(self%recv_buffer_ghost_gpu  )) deallocate(self%recv_buffer_ghost_gpu  )
   if (allocated(self%local_map_bc_face_gpu  )) deallocate(self%local_map_bc_face_gpu  )
   if (allocated(self%local_map_bc_corner_gpu)) deallocate(self%local_map_bc_corner_gpu)
   if (allocated(self%local_map_bc_edge_gpu  )) deallocate(self%local_map_bc_edge_gpu  )

   if (allocated(self%field%local_map_ghost    )) then
      self%local_map_ghost_gpu     = self%field%local_map_ghost
   endif
   if (allocated(self%field%comm_map_recv_ghost)) then
      self%comm_map_recv_ghost_gpu = self%field%comm_map_recv_ghost
   endif
   if (allocated(self%field%comm_map_send_ghost)) then
      self%comm_map_send_ghost_gpu = self%field%comm_map_send_ghost
   endif
   if (allocated(self%field%send_buffer_ghost  ).and.size(self%field%send_buffer_ghost)>0) then
      self%send_buffer_ghost_gpu   = self%field%send_buffer_ghost
   endif
   if (allocated(self%field%recv_buffer_ghost  ).and.size(self%field%recv_buffer_ghost)>0) then
      self%recv_buffer_ghost_gpu   = self%field%recv_buffer_ghost
   endif
   if (allocated(self%field%local_map_bc_face)) then
      self%local_map_bc_face_gpu = self%field%local_map_bc_face
   endif
   if (allocated(self%field%local_map_bc_corner)) then
      self%local_map_bc_corner_gpu = self%field%local_map_bc_corner
   endif
   if (allocated(self%field%local_map_bc_edge)) then
      self%local_map_bc_edge_gpu = self%field%local_map_bc_edge
   endif
   call self%create_maps_cell
   endsubroutine copy_cpu_gpu

   subroutine create_maps_cell(self)
   !< Create maps in cells order form the fecs ordered ones.
   class(base_gpu_object), intent(inout) :: self !< The base backend.
   integer(I8P), allocatable             :: local_map_ghost_cell(:,:)     !< Local map ghost cells update, cells order.
   integer(I8P), allocatable             :: comm_map_send_ghost_cell(:,:) !< MPI send map ghost cells update, cells order.
   integer(I8P), allocatable             :: comm_map_recv_ghost_cell(:,:) !< MPI send map ghost cells update, cells order.
   integer(I4P)                          :: c, f, v, n                    !< Counter.
   integer(I4P)                          :: i, j, k                       !< Counter.
   integer(I4P)                          :: iii, jjj, kkk                 !< Counter.
   integer(I4P)                          :: ic, jc, kc                    !< Counter.
   integer(I4P)                          :: fec                           !< Ghost direction, faces/edges/corners.
   integer(I4P)                          :: portion                       !< Portion of fec updated (0=>whole fec).
   integer(I4P)                          :: b_recv                        !< Index of receiving block.
   integer(I4P)                          :: b_send                        !< Index of sending block.
   integer(I4P)                          :: imin                          !< Lower limit of i indexes.
   integer(I4P)                          :: jmin                          !< Lower limit of j indexes.
   integer(I4P)                          :: kmin                          !< Lower limit of j indexes.
   integer(I4P)                          :: imax                          !< Upper limit of i indexes.
   integer(I4P)                          :: jmax                          !< Upper limit of j indexes.
   integer(I4P)                          :: kmax                          !< Upper limit of k indexes.
   integer(I4P)                          :: idelta                        !< Delta offset for ghost-inner cells of i.
   integer(I4P)                          :: jdelta                        !< Delta offset for ghost-inner cells of j.
   integer(I4P)                          :: kdelta                        !< Delta offset for ghost-inner cells of k.
   integer(I4P)                          :: recv_ptr, recv_ctr            !< Counter.
   integer(I4P)                          :: send_ptr, send_ctr            !< Counter.

   if (allocated(self%field%local_map_ghost)) then
      c = 0
      do f=1, size(self%field%local_map_ghost, dim=1)
         portion = self%field%local_map_ghost(f, 4 )
         imin    = self%field%local_map_ghost(f, 5 )
         jmin    = self%field%local_map_ghost(f, 6 )
         kmin    = self%field%local_map_ghost(f, 7 )
         imax    = self%field%local_map_ghost(f, 8 )
         jmax    = self%field%local_map_ghost(f, 9 )
         kmax    = self%field%local_map_ghost(f, 10)
         if (portion>=0) then
            ! receiving from a block with the same refinement or finer than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     c = c + 1
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do kc=0,1 ; do jc=0,1 ; do ic=0,1
                        c = c + 1
                     enddo ; enddo ; enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      allocate(local_map_ghost_cell(1:c,1:9))
      c = 1
      do f=1, size(self%field%local_map_ghost, dim=1)
         b_recv  = self%field%local_map_ghost(f, 1 )
         b_send  = self%field%local_map_ghost(f, 2 )
         portion = self%field%local_map_ghost(f, 4 )
         imin    = self%field%local_map_ghost(f, 5 )
         jmin    = self%field%local_map_ghost(f, 6 )
         kmin    = self%field%local_map_ghost(f, 7 )
         imax    = self%field%local_map_ghost(f, 8 )
         jmax    = self%field%local_map_ghost(f, 9 )
         kmax    = self%field%local_map_ghost(f, 10)
         idelta  = self%field%local_map_ghost(f, 11)
         jdelta  = self%field%local_map_ghost(f, 12)
         kdelta  = self%field%local_map_ghost(f, 13)
         if     (portion==0) then
            ! receiving from a block with the same refinement
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     local_map_ghost_cell(c,1:2) = [b_send, b_recv]
                     local_map_ghost_cell(c,3:5) = [i+idelta, j+jdelta, k+kdelta]
                     local_map_ghost_cell(c,6:8) = [i, j, k]
                     local_map_ghost_cell(c, 9 ) = 1
                     c = c + 1
                  enddo
               enddo
            enddo
         elseif (portion>0) then
            ! receiving from a block finer than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     kkk = 2 * k + kdelta
                     jjj = 2 * j + jdelta
                     iii = 2 * i + idelta
                     local_map_ghost_cell(c,1:2) = [b_send, b_recv]
                     local_map_ghost_cell(c,3:5) = [iii, jjj, kkk]
                     local_map_ghost_cell(c,6:8) = [i, j, k]
                     local_map_ghost_cell(c, 9 ) = 8
                     c = c + 1
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     kkk = 2 * k + kdelta
                     jjj = 2 * j + jdelta
                     iii = 2 * i + idelta
                     do kc=0,1 ; do jc=0,1 ; do ic=0,1
                        local_map_ghost_cell(c,1:2) = [b_send, b_recv]
                        local_map_ghost_cell(c,3:5) = [i, j, k]
                        local_map_ghost_cell(c,6:8) = [iii+ic,  jjj+jc, kkk+kc]
                        local_map_ghost_cell(c, 9 ) = 1
                        c = c + 1
                     enddo ; enddo ; enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      self%local_map_ghost_cell_gpu = local_map_ghost_cell
      deallocate(local_map_ghost_cell)
   endif
   if (allocated(self%field%comm_map_send_ghost)) then
      c = 0
      do f=1, size(self%field%comm_map_send_ghost, dim=1)
         portion = self%field%comm_map_send_ghost(f, 5 )
         imin    = self%field%comm_map_send_ghost(f, 6 )
         jmin    = self%field%comm_map_send_ghost(f, 7 )
         kmin    = self%field%comm_map_send_ghost(f, 8 )
         imax    = self%field%comm_map_send_ghost(f, 9 )
         jmax    = self%field%comm_map_send_ghost(f, 10)
         kmax    = self%field%comm_map_send_ghost(f, 11)
         if (portion>=0_I4P) then
            ! sending to a block at my level or to a block coarser than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     c = c + 1
                  enddo
               enddo
            enddo
         else
            ! sending to a block finer than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do kc=0,1 ; do jc=0,1 ; do ic=0,1
                        c = c + 1
                     enddo ; enddo ; enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      allocate(comm_map_send_ghost_cell(1:c*self%field%nv,1:7))
      c = 1
      do f=1, size(self%field%comm_map_send_ghost, dim=1)
         b_send    = self%field%comm_map_send_ghost(f, 2 )
         portion   = self%field%comm_map_send_ghost(f, 5 )
         imin      = self%field%comm_map_send_ghost(f, 6 )
         jmin      = self%field%comm_map_send_ghost(f, 7 )
         kmin      = self%field%comm_map_send_ghost(f, 8 )
         imax      = self%field%comm_map_send_ghost(f, 9 )
         jmax      = self%field%comm_map_send_ghost(f, 10)
         kmax      = self%field%comm_map_send_ghost(f, 11)
         idelta    = self%field%comm_map_send_ghost(f, 12)
         jdelta    = self%field%comm_map_send_ghost(f, 13)
         kdelta    = self%field%comm_map_send_ghost(f, 14)
         send_ptr  = self%field%comm_map_send_ghost(f, 15)
         if (portion==0_I4P) then
            ! sending to a block at my level
            send_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do v=1,self%field%nv
                        comm_map_send_ghost_cell(c, 1:4) = [b_send,i+idelta,j+jdelta,k+kdelta]
                        comm_map_send_ghost_cell(c,  5 ) = v
                        comm_map_send_ghost_cell(c,  6 ) = send_ptr + send_ctr
                        comm_map_send_ghost_cell(c,  7 ) = 1
                        send_ctr = send_ctr + 1
                        c = c + 1
                     enddo
                  enddo
               enddo
            enddo
         elseif (portion<0_I4P) then ! Beware! This is < 0 because the reference is the receiver
            ! sending to a block finer than me
            send_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do n=1,8
                        do v=1,self%field%nv
                           comm_map_send_ghost_cell(c, 1:4) = [b_send,i,j,k]
                           comm_map_send_ghost_cell(c,  5 ) = v
                           comm_map_send_ghost_cell(c,  6 ) = send_ptr + send_ctr
                           comm_map_send_ghost_cell(c,  7 ) = 1
                           send_ctr = send_ctr + 1
                           c = c + 1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         else
            ! sending to a block coarser than me, loop is over the coarser grid
            send_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     kkk = 2 * k + kdelta
                     jjj = 2 * j + jdelta
                     iii = 2 * i + idelta
                     do v=1,self%field%nv
                        comm_map_send_ghost_cell(c, 1:4) = [b_send,iii,jjj,kkk]
                        comm_map_send_ghost_cell(c,  5 ) = v
                        comm_map_send_ghost_cell(c,  6 ) = send_ptr + send_ctr
                        comm_map_send_ghost_cell(c,  7 ) = 8
                        send_ctr = send_ctr + 1
                        c = c + 1
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      self%comm_map_send_ghost_cell_gpu = comm_map_send_ghost_cell
      deallocate(comm_map_send_ghost_cell)
   endif
   if (allocated(self%field%comm_map_recv_ghost)) then
      c = 0
      do f=1, size(self%field%comm_map_recv_ghost, dim=1)
         portion = self%field%comm_map_recv_ghost(f, 5 )
         imin    = self%field%comm_map_recv_ghost(f, 6 )
         jmin    = self%field%comm_map_recv_ghost(f, 7 )
         kmin    = self%field%comm_map_recv_ghost(f, 8 )
         imax    = self%field%comm_map_recv_ghost(f, 9 )
         jmax    = self%field%comm_map_recv_ghost(f, 10)
         kmax    = self%field%comm_map_recv_ghost(f, 11)
         if (portion>=0_I4P) then
            ! receiving from a block at the same level or finer than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     c = c + 1
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do kc=0,1 ; do jc=0,1 ; do ic=0,1
                        c = c + 1
                     enddo ; enddo ; enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      allocate(comm_map_recv_ghost_cell(1:c*self%field%nv,1:6))
      c = 1
      do f=1, size(self%field%comm_map_recv_ghost, dim=1)
         b_recv   = self%field%comm_map_recv_ghost(f, 1 )
         portion  = self%field%comm_map_recv_ghost(f, 5 )
         imin     = self%field%comm_map_recv_ghost(f, 6 )
         jmin     = self%field%comm_map_recv_ghost(f, 7 )
         kmin     = self%field%comm_map_recv_ghost(f, 8 )
         imax     = self%field%comm_map_recv_ghost(f, 9 )
         jmax     = self%field%comm_map_recv_ghost(f, 10)
         kmax     = self%field%comm_map_recv_ghost(f, 11)
         idelta   = self%field%comm_map_recv_ghost(f, 12)
         jdelta   = self%field%comm_map_recv_ghost(f, 13)
         kdelta   = self%field%comm_map_recv_ghost(f, 14)
         recv_ptr = self%field%comm_map_recv_ghost(f, 15)
         if (portion>=0_I4P) then
            recv_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do v=1, self%field%nv
                        comm_map_recv_ghost_cell(c, 1 ) = recv_ptr + recv_ctr
                        comm_map_recv_ghost_cell(c,2:6) = [b_recv,i,j,k,v]
                        recv_ctr = recv_ctr + 1
                        c = c + 1
                     enddo
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            recv_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     kkk = 2 * k + kdelta
                     jjj = 2 * j + jdelta
                     iii = 2 * i + idelta

                     do kc=0,1 ; do jc=0,1 ; do ic=0,1
                        do v=1, self%field%nv
                           comm_map_recv_ghost_cell(c, 1 ) = recv_ptr + recv_ctr
                           comm_map_recv_ghost_cell(c,2:6) = [b_recv,iii+ic,jjj+jc,kkk+kc,v]
                           recv_ctr = recv_ctr + 1
                           c = c + 1
                        enddo
                     enddo ; enddo ; enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      self%comm_map_recv_ghost_cell_gpu = comm_map_recv_ghost_cell
      deallocate(comm_map_recv_ghost_cell)
   endif
   endsubroutine create_maps_cell

   subroutine destroy(self)
   !< Destroy base backend.
   class(base_gpu_object), intent(inout) :: self  !< The base backend.
   type(base_gpu_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, field, nv_aux)
   !< Initialize base backend.
   class(base_gpu_object), intent(inout)        :: self    !< The base backend.
   integer(I4P),           intent(in), optional :: nv_aux  !< Number of auxiliary variables.
   type(field_object),     intent(in), target   :: field   !< The field.
   integer(I4P)                                 :: nv_aux_ !< Number of auxiliary variables, local variables.

   call self%destroy
   self%field => field
   nv_aux_ = self%field%nv
   if (present(nv_aux)) nv_aux_ = max(nv_aux_, nv_aux)
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   allocate(self%req_send_recv(0:self%procs_number*2-1))
   call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%local_comm, self%error)
   call MPI_COMM_RANK(self%local_comm, self%mydev,self%error)
   self%error = CudaSetDevice(self%mydev)
   allocate(self%q_t(1:field%nb,                                    &
                     1-field%grid%gci:field%grid%ni+field%grid%gci, &
                     1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                     1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nv))
   allocate(self%q_t_gpu(nv_aux_,                                       &
                         1-field%grid%gci:field%grid%ni+field%grid%gci, &
                         1-field%grid%gcj:field%grid%nj+field%grid%gcj, &
                         1-field%grid%gck:field%grid%nk+field%grid%gck, 1:field%nb))
   call self%copy_cpu_gpu
   endsubroutine initialize

   subroutine update_ghost_local_gpu(self, q_gpu)
   !< Update (local) ghost cells.
   class(base_gpu_object), intent(in)            :: self      !< The base backend.
   real(R8P),              intent(inout), device :: q_gpu(1:,                    &
                                                          1-self%field%grid%gci:,&
                                                          1-self%field%grid%gcj:,&
                                                          1-self%field%grid%gck:,&
                                                          1:) !< Field component to be updated.
   call update_ghost_local_gpu_cuf(local_map_ghost_cell_gpu=self%local_map_ghost_cell_gpu, &
                                   gci=self%field%grid%gci, gcj=self%field%grid%gcj, gck=self%field%grid%gck, q_gpu=q_gpu)
   endsubroutine update_ghost_local_gpu

   subroutine update_ghost_mpi_gpu(self, q_gpu, step)
   !< Update ghost cells within other processes.
   class(base_gpu_object), intent(inout)         :: self      !< The base backend.
   real(R8P),              intent(inout), device :: q_gpu(1:,                    &
                                                          1-self%field%grid%gci:,&
                                                          1-self%field%grid%gcj:,&
                                                          1-self%field%grid%gck:,&
                                                          1:) !< Field component to be updated.
   integer(I4P),           intent(in), optional  :: step      !< Step to be perfordmed in asyncronous comp.

   call update_ghost_mpi_gpu_cuf(nv=self%field%nv,                                           &
                                 procs_number=self%field%procs_number,                       &
                                 req_send_recv=self%field%req_send_recv,                     &
                                 comm_map_send_ptr_ghost=self%field%comm_map_send_ptr_ghost, &
                                 comm_map_recv_ptr_ghost=self%field%comm_map_recv_ptr_ghost, &
                                 comm_map_send_ghost_cell_gpu=self%comm_map_send_ghost_cell_gpu, &
                                 comm_map_recv_ghost_cell_gpu=self%comm_map_recv_ghost_cell_gpu, &
                                 ! comm_map_send_ghost_gpu=self%comm_map_send_ghost_gpu,       &
                                 ! comm_map_recv_ghost_gpu=self%comm_map_recv_ghost_gpu,       &
                                 recv_buffer_ghost_gpu=self%recv_buffer_ghost_gpu,           &
                                 send_buffer_ghost_gpu=self%send_buffer_ghost_gpu,           &
                                 gci=self%field%grid%gci, gcj=self%field%grid%gcj, gck=self%field%grid%gck, q_gpu=q_gpu, step=step)
   endsubroutine update_ghost_mpi_gpu

   ! operators
   ! =
   subroutine base_assign_base(lhs, rhs)
   !< Operator `=`.
   class(base_gpu_object), intent(inout) :: lhs !< Left hand side.
   type(base_gpu_object),  intent(in)    :: rhs !< Right hand side.

   lhs%field => rhs%field
   lhs%myrank = rhs%myrank
   lhs%procs_number = rhs%procs_number
   lhs%error = rhs%error
   call assign_allocatable(lhs=lhs%req_send_recv, rhs=rhs%req_send_recv)
   call assign_allocatable(lhs=lhs%q_t,           rhs=rhs%q_t          )
   call assign_allocatable_gpu(lhs=lhs%q_t_gpu,                 rhs=rhs%q_t_gpu                )
   call assign_allocatable_gpu(lhs=lhs%local_map_ghost_gpu    , rhs=rhs%local_map_ghost_gpu    )
   call assign_allocatable_gpu(lhs=lhs%comm_map_recv_ghost_gpu, rhs=rhs%comm_map_recv_ghost_gpu)
   call assign_allocatable_gpu(lhs=lhs%comm_map_send_ghost_gpu, rhs=rhs%comm_map_send_ghost_gpu)
   call assign_allocatable_gpu(lhs=lhs%send_buffer_ghost_gpu  , rhs=rhs%send_buffer_ghost_gpu  )
   call assign_allocatable_gpu(lhs=lhs%recv_buffer_ghost_gpu  , rhs=rhs%recv_buffer_ghost_gpu  )
   endsubroutine base_assign_base

   ! private methods
   subroutine copy_transpose_cpu_gpu(self, q_cpu, q_gpu)
   !< Copy transposed data from CPU to GPU.
   class(base_gpu_object), intent(inout)       :: self          !< The equation.
   real(R8P),              intent(in)          :: q_cpu(1:,                    &
                                                        1-self%field%grid%gci:,&
                                                        1-self%field%grid%gcj:,&
                                                        1-self%field%grid%gck:,&
                                                        1:)     !< Conservative variables on CPU.
   real(R8P),              intent(out), device :: q_gpu(1:,                    &
                                                       1-self%field%grid%gci:,&
                                                       1-self%field%grid%gcj:,&
                                                       1-self%field%grid%gck:,&
                                                       1:)      !< Conservative variables on GPU.
   integer(I4P)                                :: i, j, k, b, v !< Counter.
   associate(blocks_number=>self%field%blocks_number,                                      &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,       &
             gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck, &
             nv=>self%field%nv, q_t=>self%q_t)
      do b=1, blocks_number
         do k=1-gck, nk+gck
            do j=1-gcj, nj+gcj
               do i=1-gci, ni+gci
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
   class(base_gpu_object), intent(inout)      :: self      !< The equation.
   integer(I4P),           intent(in)         :: nv        !< Number of conservative varibales.
   real(R8P),              intent(in), device :: q_gpu(1:,                    &
                                                       1-self%field%grid%gci:,&
                                                       1-self%field%grid%gcj:,&
                                                       1-self%field%grid%gck:,&
                                                       1:) !< Conservative variables on GPU.
   real(R8P),              intent(out)        :: q_cpu(1:,                    &
                                                       1-self%field%grid%gci:,&
                                                       1-self%field%grid%gcj:,&
                                                       1-self%field%grid%gck:,&
                                                       1:) !< Conservative variables on CPU.
   associate(blocks_number=>self%field%blocks_number,                                      &
             ni=>self%field%grid%ni, nj=>self%field%grid%nj, nk=>self%field%grid%nk,       &
             gci=>self%field%grid%gci, gcj=>self%field%grid%gcj, gck=>self%field%grid%gck)
      call copy_transpose_gpu_cpu_cuf(ni=ni, nj=nj, nk=nk, gci=gci, gcj=gcj, gck=gck, nv=nv, &
                                      blocks_number=blocks_number,                           &
                                      q_gpu=q_gpu, q_t_gpu=self%q_t_gpu, q_cpu=q_cpu)
   endassociate
   endsubroutine copy_transpose_gpu_cpu

   ! non TBP CUF methods
   subroutine copy_transpose_gpu_cpu_cuf(ni, nj, nk, gci, gcj, gck, nv, blocks_number, q_gpu, q_t_gpu, q_cpu)
   !< Copy transposed data from GPU to CPU by CUF threads.
   integer(I4P), intent(in)            :: ni            !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj            !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk            !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: gci           !< Ghost grid cells number in I direction.
   integer(I4P), intent(in)            :: gcj           !< Ghost grid cells number in J direction.
   integer(I4P), intent(in)            :: gck           !< Ghost grid cells number in K direction.
   integer(I4P), intent(in)            :: nv            !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number !< Number of blocks.
   real(R8P),    intent(in), device    :: q_gpu(1:,    &
                                                1-gci:,&
                                                1-gcj:,&
                                                1-gck:,&
                                                1:)     !< Conservative variables on GPU.
   real(R8P),    intent(inout), device :: q_t_gpu(1:,    &
                                                  1-gci:,&
                                                  1-gcj:,&
                                                  1-gck:,&
                                                  1:)   !< Conservative (transposed) variables on GPU.
   real(R8P),    intent(out)           :: q_cpu(1:,    &
                                                1-gci:,&
                                                1-gcj:,&
                                                1-gck:,&
                                                1:)     !< Conservative variables on CPU.
   integer(I4P)                        :: i, j, k, b, v !< Counter.
   integer(I4P)                        :: iercuda       !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do k=1-gck, nk+gck
      do j=1-gcj, nj+gcj
         do i=1-gci, ni+gci
            do b=1, blocks_number
               do v=1, nv
                  q_t_gpu(v,i,j,k,b) = q_gpu(b,i,j,k,v)
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   q_cpu = q_t_gpu
   endsubroutine copy_transpose_gpu_cpu_cuf

   subroutine update_ghost_local_gpu_cuf(gci, gcj, gck, local_map_ghost_cell_gpu, q_gpu)
   !< Update (local) ghost cells.
   integer(I4P), intent(in)                         :: gci                           !< Ghost cells number in I direction.
   integer(I4P), intent(in)                         :: gcj                           !< Ghost cells number in J direction.
   integer(I4P), intent(in)                         :: gck                           !< Ghost cells number in K direction.
   integer(I8P), intent(in),    device, allocatable :: local_map_ghost_cell_gpu(:,:) !< Local map of ghost cells.
   real(R8P),    intent(inout), device              :: q_gpu(1:,    &
                                                             1-gci:,&
                                                             1-gcj:,&
                                                             1-gck:,1:)              !< Field component to be updated.
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
   integer(I4P)                                     :: iercuda                       !< Error trapping flag for CUDAFortran.

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
   endsubroutine update_ghost_local_gpu_cuf

   subroutine update_ghost_mpi_gpu_cuf(gci, gcj, gck, nv, procs_number, req_send_recv,   &
                                       comm_map_send_ptr_ghost, comm_map_recv_ptr_ghost, &
                                       ! comm_map_recv_ghost_gpu, comm_map_send_ghost_gpu, &
                                       comm_map_recv_ghost_cell_gpu, comm_map_send_ghost_cell_gpu, &
                                       recv_buffer_ghost_gpu, send_buffer_ghost_gpu , q_gpu, step)
   !< Update ghost cells within other processes.
   integer(I4P), intent(in)                         :: gci                          !< Ghost cells number in I direction.
   integer(I4P), intent(in)                         :: gcj                          !< Ghost cells number in J direction.
   integer(I4P), intent(in)                         :: gck                          !< Ghost cells number in K direction.
   integer(I4P),              intent(in)            :: nv                           !< Number of variables of q.
   integer(I4P),              intent(in)            :: procs_number                 !< Number of MPI processes.
   integer(I4P), allocatable, intent(inout)         :: req_send_recv(:)             !< MPI request receive flags.
   integer(I4P), allocatable, intent(in)            :: comm_map_send_ptr_ghost(:)   !< Communication map, pointers list to send.
   integer(I4P), allocatable, intent(in)            :: comm_map_recv_ptr_ghost(:)   !< Communication map, pointers list to recv.
   ! integer(I8P), allocatable, intent(in),    device :: comm_map_recv_ghost_gpu(:,:) !< Communication map, `fec` information.
   ! integer(I8P), allocatable, intent(in),    device :: comm_map_send_ghost_gpu(:,:) !< Communication map, `fec` information.
   integer(I8P), allocatable, intent(in),    device :: comm_map_recv_ghost_cell_gpu(:,:) !< Communication map, cell information.
   integer(I8P), allocatable, intent(in),    device :: comm_map_send_ghost_cell_gpu(:,:) !< Communication map, cell information.
   real(R8P),    allocatable, intent(inout), device :: recv_buffer_ghost_gpu(:)     !< Receive buffer of ghost cells.
   real(R8P),    allocatable, intent(inout), device :: send_buffer_ghost_gpu(:)     !< Send buffer of ghost cells.
   real(R8P),                 intent(inout), device :: q_gpu(1:,    &
                                                             1-gci:,&
                                                             1-gcj:,&
                                                             1-gck:,&
                                                             1:)                    !< Field component to be updated.
   integer(I4P),              intent(in), optional  :: step                         !< Step to be perfordmed in asyncronous comp.
   logical                                          :: do_step(3)                   !< Steps to be performed in asyncronous comp.
   integer(I4P)                                     :: i, j, k                      !< Counter.
   integer(I4P)                                     :: ic, jc, kc                   !< Counter.
   integer(I4P)                                     :: fec, mf, rf, sf, n, p, v     !< Counter.
   integer(I4P)                                     :: portion                      !< Portion of fec updated (0=>whole fec).
   integer(I4P)                                     :: b_send, i_send, j_send, k_send, v_send !< Send indexes.
   integer(I4P)                                     :: b_recv, i_recv, j_recv, k_recv, v_recv !< Receive indexes.
   integer(I4P)                                     :: ptr_start, ptr_end           !< Counter.
   integer(I4P)                                     :: n_recv, n_send               !< Counter.
   integer(I4P)                                     :: recv_rank                    !< Rank of receiving block.
   integer(I4P)                                     :: error                        !< Error traping flag.
   integer(I4P)                                     :: send_ptr, send_ctr           !< Counter.
   integer(I4P)                                     :: recv_ptr, recv_ctr           !< Counter.
   integer(I4P)                                     :: c_send, c_recv               !< Counter.
   integer(I4P)                                     :: one_or_eight                 !< Flag triggering 8 cells mean.
   integer(I4P)                                     :: iercuda                      !< Error trapping flag for CUDAFortran.

   ! if ((.not.allocated(comm_map_recv_ghost_gpu)).and.&
       ! (.not.allocated(comm_map_send_ghost_gpu))) return
   if ((.not.allocated(comm_map_recv_ghost_cell_gpu)).and.&
       (.not.allocated(comm_map_send_ghost_cell_gpu))) return

   do_step = .true.
   if (present(step)) then
      do_step = .false.
      do_step(step) = .true.
   endif

   if (do_step(1)) then
      req_send_recv = MPI_REQUEST_NULL

      ! populate send buffer
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
      ! do sf=1, size(comm_map_send_ghost_gpu, dim=1)
      !    b_send    = comm_map_send_ghost_gpu(sf, 2 )
      !    send_rank = comm_map_send_ghost_gpu(sf, 3 )
      !    fec       = comm_map_send_ghost_gpu(sf, 4 )
      !    portion   = comm_map_send_ghost_gpu(sf, 5 )
      !    imin      = comm_map_send_ghost_gpu(sf, 6 )
      !    jmin      = comm_map_send_ghost_gpu(sf, 7 )
      !    kmin      = comm_map_send_ghost_gpu(sf, 8 )
      !    imax      = comm_map_send_ghost_gpu(sf, 9 )
      !    jmax      = comm_map_send_ghost_gpu(sf, 10)
      !    kmax      = comm_map_send_ghost_gpu(sf, 11)
      !    idelta    = comm_map_send_ghost_gpu(sf, 12)
      !    jdelta    = comm_map_send_ghost_gpu(sf, 13)
      !    kdelta    = comm_map_send_ghost_gpu(sf, 14)
      !    send_ptr  = comm_map_send_ghost_gpu(sf, 15)
      !    if (portion==0_I4P) then
      !       ! sending to a block at my level
      !       send_ctr = 1
      !       do k=kmin, kmax
      !          do j=jmin, jmax
      !             do i=imin, imax
      !                do v=1,nv
      !                   send_buffer_ghost_gpu(send_ptr + send_ctr) = q_gpu(b_send,i+idelta,j+jdelta,k+kdelta,v)
      !                   send_ctr = send_ctr + 1
      !                enddo
      !             enddo
      !          enddo
      !       enddo
      !    elseif (portion<0_I4P) then ! Beware! This is < 0 because the reference is the receiver
      !       ! sending to a block finer than me
      !       send_ctr = 1
      !       do k=kmin, kmax
      !          do j=jmin, jmax
      !             do i=imin, imax
      !                do n=1,8
      !                   do v=1,nv
      !                      send_buffer_ghost_gpu(send_ptr + send_ctr) = q_gpu(b_send,i,j,k,v)
      !                      send_ctr = send_ctr + 1
      !                   enddo
      !                enddo
      !             enddo
      !          enddo
      !       enddo
      !    else
      !       ! sending to a block coarser than me, loop is over the coarser grid
      !       send_ctr = 1
      !       do k=kmin, kmax
      !          do j=jmin, jmax
      !             do i=imin, imax
      !                kkk = 2 * k + kdelta
      !                jjj = 2 * j + jdelta
      !                iii = 2 * i + idelta
      !                do v=1,nv
      !                   send_buffer_ghost_gpu(send_ptr + send_ctr) =                               &
      !                      (q_gpu(b_send,iii,jjj,  kkk,  v) + q_gpu(b_send,iii+1,jjj,  kkk,  v) +  &
      !                       q_gpu(b_send,iii,jjj+1,kkk,  v) + q_gpu(b_send,iii+1,jjj+1,kkk,  v) +  &
      !                       q_gpu(b_send,iii,jjj,  kkk+1,v) + q_gpu(b_send,iii+1,jjj,  kkk+1,v) +  &
      !                       q_gpu(b_send,iii,jjj+1,kkk+1,v) + q_gpu(b_send,iii+1,jjj+1,kkk+1,v)) / 8._R8P
      !                   send_ctr = send_ctr + 1
      !                enddo
      !             enddo
      !          enddo
      !       enddo
      !    endif
      ! enddo
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
      ! do rf=1, size(comm_map_recv_ghost_gpu, dim=1)
      !    b_recv    = comm_map_recv_ghost_gpu(rf, 1 )
      !    recv_rank = comm_map_recv_ghost_gpu(rf, 3 )
      !    fec       = comm_map_recv_ghost_gpu(rf, 4 )
      !    portion   = comm_map_recv_ghost_gpu(rf, 5 )
      !    imin      = comm_map_recv_ghost_gpu(rf, 6 )
      !    jmin      = comm_map_recv_ghost_gpu(rf, 7 )
      !    kmin      = comm_map_recv_ghost_gpu(rf, 8 )
      !    imax      = comm_map_recv_ghost_gpu(rf, 9 )
      !    jmax      = comm_map_recv_ghost_gpu(rf, 10)
      !    kmax      = comm_map_recv_ghost_gpu(rf, 11)
      !    idelta    = comm_map_recv_ghost_gpu(rf, 12)
      !    jdelta    = comm_map_recv_ghost_gpu(rf, 13)
      !    kdelta    = comm_map_recv_ghost_gpu(rf, 14)
      !    recv_ptr  = comm_map_recv_ghost_gpu(rf, 15)
      !    if (portion==0_I4P) then
      !       recv_ctr = 1
      !       do k=kmin, kmax
      !          do j=jmin, jmax
      !             do i=imin, imax
      !                do v=1, nv
      !                   q_gpu(b_recv,i,j,k,v) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
      !                   recv_ctr = recv_ctr + 1
      !                enddo
      !             enddo
      !          enddo
      !       enddo
      !    elseif (portion>0_I4P) then
      !       ! receiving from a block finer than me
      !       recv_ctr = 1
      !       do k=kmin, kmax
      !          do j=jmin, jmax
      !             do i=imin, imax
      !                do v=1, nv
      !                   q_gpu(b_recv,i,j,k,v) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
      !                   recv_ctr = recv_ctr + 1
      !                enddo
      !             enddo
      !          enddo
      !       enddo
      !    else
      !       ! receiving from a block coarser than me
      !       recv_ctr = 1
      !       do k=kmin, kmax
      !          do j=jmin, jmax
      !             do i=imin, imax
      !                kkk = 2 * k + kdelta
      !                jjj = 2 * j + jdelta
      !                iii = 2 * i + idelta
      !                do v=1, nv
      !                   q_gpu(b_recv,iii,  jjj,  kkk  ,v) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
      !                enddo
      !                do v=1, nv
      !                   q_gpu(b_recv,iii+1,jjj,  kkk  ,v) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
      !                enddo
      !                do v=1, nv
      !                   q_gpu(b_recv,iii,  jjj+1,kkk  ,v) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
      !                enddo
      !                do v=1, nv
      !                   q_gpu(b_recv,iii+1,jjj+1,kkk  ,v) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
      !                enddo
      !                do v=1, nv
      !                   q_gpu(b_recv,iii,  jjj,  kkk+1,v) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
      !                enddo
      !                do v=1, nv
      !                   q_gpu(b_recv,iii+1,jjj,  kkk+1,v) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
      !                enddo
      !                do v=1, nv
      !                   q_gpu(b_recv,iii,  jjj+1,kkk+1,v) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
      !                enddo
      !                do v=1, nv
      !                   q_gpu(b_recv,iii+1,jjj+1,kkk+1,v) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
      !                enddo
      !             enddo
      !          enddo
      !       enddo
      !    endif
      ! enddo
   endif
   endsubroutine update_ghost_mpi_gpu_cuf

   ! assign allocatable interface
   pure subroutine assign_allocatable_I8P_1D_gpu(lhs, rhs)
   !< Safe assign allocatable arrays, I8P 1D type.
   integer(I8P), allocatable, device, intent(inout) :: lhs(:) !< Left hand side.
   integer(I8P), allocatable, device, intent(in)    :: rhs(:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_I8P_1D_gpu

   pure subroutine assign_allocatable_I8P_2D_gpu(lhs, rhs)
   !< Safe assign allocatable arrays, I8P 2D type.
   integer(I8P), allocatable, device, intent(inout) :: lhs(:,:) !< Left hand side.
   integer(I8P), allocatable, device, intent(in)    :: rhs(:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_I8P_2D_gpu

   pure subroutine assign_allocatable_I4P_1D_gpu(lhs, rhs)
   !< Safe assign allocatable arrays, I4P 1D type.
   integer(I4P), allocatable, device, intent(inout) :: lhs(:) !< Left hand side.
   integer(I4P), allocatable, device, intent(in)    :: rhs(:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_I4P_1D_gpu

   pure subroutine assign_allocatable_I4P_2D_gpu(lhs, rhs)
   !< Safe assign allocatable arrays, I4P 2D type.
   integer(I4P), allocatable, device, intent(inout) :: lhs(:,:) !< Left hand side.
   integer(I4P), allocatable, device, intent(in)    :: rhs(:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_I4P_2D_gpu

   pure subroutine assign_allocatable_R8P_1D_gpu(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 1D type.
   real(R8P), allocatable, device, intent(inout) :: lhs(:) !< Left hand side.
   real(R8P), allocatable, device, intent(in)    :: rhs(:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_1D_gpu

   pure subroutine assign_allocatable_R8P_2D_gpu(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 2D type.
   real(R8P), allocatable, device, intent(inout) :: lhs(:,:) !< Left hand side.
   real(R8P), allocatable, device, intent(in)    :: rhs(:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_2D_gpu

   pure subroutine assign_allocatable_R8P_3D_gpu(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 3D type.
   real(R8P), allocatable, device, intent(inout) :: lhs(:,:,:) !< Left hand side.
   real(R8P), allocatable, device, intent(in)    :: rhs(:,:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_3D_gpu

   pure subroutine assign_allocatable_R8P_4D_gpu(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 4D type.
   real(R8P), allocatable, device, intent(inout) :: lhs(:,:,:,:) !< Left hand side.
   real(R8P), allocatable, device, intent(in)    :: rhs(:,:,:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_4D_gpu

   pure subroutine assign_allocatable_R8P_5D_gpu(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 5D type.
   real(R8P), allocatable, device, intent(inout) :: lhs(:,:,:,:,:) !< Left hand side.
   real(R8P), allocatable, device, intent(in)    :: rhs(:,:,:,:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_5D_gpu

   pure subroutine assign_allocatable_R8P_6D_gpu(lhs, rhs)
   !< Safe assign allocatable arrays, R8P 6D type.
   real(R8P), allocatable, device, intent(inout) :: lhs(:,:,:,:,:,:) !< Left hand side.
   real(R8P), allocatable, device, intent(in)    :: rhs(:,:,:,:,:,:) !< Right hand side.

   if (allocated(rhs)) then
      lhs = rhs
   else
      if (allocated(lhs)) deallocate(lhs)
   endif
   endsubroutine assign_allocatable_R8P_6D_gpu
endmodule adam_base_gpu_object
