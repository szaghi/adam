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
   integer(I4P) :: myrank=0_I4P                  !< MPI rank process.
   integer(I4P) :: procs_number=1_I4P            !< Number of MPI processes.
   integer(I4P) :: error=0_I4P                   !< Error traping flag.
   integer(I4P) :: mydev=0_I4P                   !< My GPU rank.
   integer(I4P) :: local_comm=0_I4P              !< Local communicator.
   integer(I4P), allocatable :: req_send_recv(:) !< MPI request receive flags.
   ! GPU data
   integer(I8P), allocatable, device :: local_map_ghost_gpu(:,:)     !< Local map for ghost cells updating.
   integer(I8P), allocatable, device :: comm_map_recv_ghost_gpu(:,:) !< Communication map, `fec` information.
   integer(I8P), allocatable, device :: comm_map_send_ghost_gpu(:,:) !< Communication map, `fec` information.
   real(R8P),    allocatable, device :: send_buffer_ghost_gpu(:)     !< Send buffer of ghost cells.
   real(R8P),    allocatable, device :: recv_buffer_ghost_gpu(:)     !< Receive buffer of ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_face_gpu(:,:)   !< Local map for face BC ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_edge_gpu(:,:)   !< Local map for edge BC ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_corner_gpu(:,:) !< Local map for corner BC ghost cells.
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu           !< Copy data from CPU to GPU.
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
   endsubroutine copy_cpu_gpu

   subroutine destroy(self)
   !< Destroy base backend.
   class(base_gpu_object), intent(inout) :: self  !< The base backend.
   type(base_gpu_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, field)
   !< Initialize base backend.
   class(base_gpu_object), intent(inout)      :: self  !< The base backend.
   type(field_object),     intent(in), target :: field !< The field.

   call self%destroy
   self%field => field
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, self%procs_number, self%error)
   allocate(self%req_send_recv(0:self%procs_number*2-1))
   call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, self%local_comm, self%error)
   call MPI_COMM_RANK(self%local_comm, self%mydev,self%error)
   self%error = CudaSetDevice(self%mydev)
   call self%copy_cpu_gpu
   endsubroutine initialize

   subroutine update_ghost_local_gpu(self, q_gpu)
   !< Update (local) ghost cells.
   class(base_gpu_object), intent(in)            :: self                                !< The base backend.
   real(R8P),              intent(inout), device :: q_gpu(1-self%field%grid%gci:,&
                                                          1-self%field%grid%gcj:,&
                                                          1-self%field%grid%gck:,1:,1:) !< Field component to be updated.
   call update_ghost_local_gpu_cuf(local_map_ghost_gpu=self%local_map_ghost_gpu, &
                                   gci=self%field%grid%gci, gcj=self%field%grid%gcj, gck=self%field%grid%gck, q_gpu=q_gpu)
   endsubroutine update_ghost_local_gpu

   subroutine update_ghost_mpi_gpu(self, q_gpu, step)
   !< Update ghost cells within other processes.
   class(base_gpu_object), intent(inout)         :: self                           !< The field.
   real(R8P),              intent(inout), device :: q_gpu(1-self%field%grid%gci:,&
                                                          1-self%field%grid%gcj:,&
                                                          1-self%field%grid%gck:,&
                                                          1:,1:)                   !< Field component to be updated.
   integer(I4P),           intent(in), optional  :: step                           !< Step to be perfordmed in asyncronous comp.

   call update_ghost_mpi_gpu_cuf(nv=self%field%nv,                                           &
                                 procs_number=self%field%procs_number,                       &
                                 req_send_recv=self%field%req_send_recv,                     &
                                 comm_map_send_ptr_ghost=self%field%comm_map_send_ptr_ghost, &
                                 comm_map_recv_ptr_ghost=self%field%comm_map_recv_ptr_ghost, &
                                 comm_map_send_ghost_gpu=self%comm_map_send_ghost_gpu,       &
                                 comm_map_recv_ghost_gpu=self%comm_map_recv_ghost_gpu,       &
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
   call assign_allocatable_gpu(lhs=lhs%local_map_ghost_gpu    , rhs=rhs%local_map_ghost_gpu    )
   call assign_allocatable_gpu(lhs=lhs%comm_map_recv_ghost_gpu, rhs=rhs%comm_map_recv_ghost_gpu)
   call assign_allocatable_gpu(lhs=lhs%comm_map_send_ghost_gpu, rhs=rhs%comm_map_send_ghost_gpu)
   call assign_allocatable_gpu(lhs=lhs%send_buffer_ghost_gpu  , rhs=rhs%send_buffer_ghost_gpu  )
   call assign_allocatable_gpu(lhs=lhs%recv_buffer_ghost_gpu  , rhs=rhs%recv_buffer_ghost_gpu  )
   endsubroutine base_assign_base

   ! non TBP CUF methods
   subroutine update_ghost_local_gpu_cuf(gci, gcj, gck, local_map_ghost_gpu, q_gpu)
   !< Update (local) ghost cells.
   integer(I4P), intent(in)                         :: gci                      !< Ghost cells number in I direction.
   integer(I4P), intent(in)                         :: gcj                      !< Ghost cells number in J direction.
   integer(I4P), intent(in)                         :: gck                      !< Ghost cells number in K direction.
   integer(I8P), intent(in),    device, allocatable :: local_map_ghost_gpu(:,:) !< Local map of ghost cells.
   real(R8P),    intent(inout), device              :: q_gpu(1-gci:,&
                                                             1-gcj:,&
                                                             1-gck:,1:,1:)      !< Field component to be updated.
   integer(I4P)                                     :: i, j, k, mf              !< Counter.
   integer(I4P)                                     :: iii, jjj, kkk            !< Counter.
   integer(I4P)                                     :: fec                      !< Ghost direction, faces/edges/corners.
   integer(I4P)                                     :: portion                  !< Portion of fec updated (0=>whole fec).
   integer(I4P)                                     :: b_recv                   !< Index of receiving block.
   integer(I4P)                                     :: b_send                   !< Index of sending block.
   integer(I4P)                                     :: imin                     !< Lower limit of i indexes.
   integer(I4P)                                     :: jmin                     !< Lower limit of j indexes.
   integer(I4P)                                     :: kmin                     !< Lower limit of j indexes.
   integer(I4P)                                     :: imax                     !< Upper limit of i indexes.
   integer(I4P)                                     :: jmax                     !< Upper limit of j indexes.
   integer(I4P)                                     :: kmax                     !< Upper limit of k indexes.
   integer(I4P)                                     :: idelta                   !< Delta offset for ghost-inner cells of i.
   integer(I4P)                                     :: jdelta                   !< Delta offset for ghost-inner cells of j.
   integer(I4P)                                     :: kdelta                   !< Delta offset for ghost-inner cells of k.
   integer(I4P)                                     :: iercuda                  !< Error trapping flag for CUDAFortran.

   if (.not.allocated(local_map_ghost_gpu)) return
   !$cuf kernel do(1) <<<*,*>>>
   do mf=1, size(local_map_ghost_gpu, dim=1)
      b_recv  = local_map_ghost_gpu(mf, 1)
      b_send  = local_map_ghost_gpu(mf, 2)
      fec     = local_map_ghost_gpu(mf, 3)
      portion = local_map_ghost_gpu(mf, 4)
      imin    = local_map_ghost_gpu(mf, 5 )
      jmin    = local_map_ghost_gpu(mf, 6 )
      kmin    = local_map_ghost_gpu(mf, 7 )
      imax    = local_map_ghost_gpu(mf, 8 )
      jmax    = local_map_ghost_gpu(mf, 9 )
      kmax    = local_map_ghost_gpu(mf, 10)
      idelta  = local_map_ghost_gpu(mf, 11)
      jdelta  = local_map_ghost_gpu(mf, 12)
      kdelta  = local_map_ghost_gpu(mf, 13)
      if     (portion==0) then
         ! receiving from a block with the same refinement
         do k=kmin, kmax
            do j=jmin, jmax
               do i=imin, imax
                  q_gpu(i,j,k,:,b_recv) = q_gpu(i+idelta,j+jdelta,k+kdelta,:,b_send)
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
                  q_gpu(i,j,k,:,b_recv) = (q_gpu(iii,jjj,  kkk,  :,b_send) + q_gpu(iii+1,jjj,  kkk,  :,b_send) + &
                                           q_gpu(iii,jjj+1,kkk,  :,b_send) + q_gpu(iii+1,jjj+1,kkk,  :,b_send) + &
                                           q_gpu(iii,jjj,  kkk+1,:,b_send) + q_gpu(iii+1,jjj,  kkk+1,:,b_send) + &
                                           q_gpu(iii,jjj+1,kkk+1,:,b_send) + q_gpu(iii+1,jjj+1,kkk+1,:,b_send)) / 8._R8P
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
                  q_gpu(iii,  jjj,  kkk  ,:,b_recv) = q_gpu(i,j,k,:,b_send)
                  q_gpu(iii+1,jjj,  kkk  ,:,b_recv) = q_gpu(i,j,k,:,b_send)
                  q_gpu(iii,  jjj+1,kkk  ,:,b_recv) = q_gpu(i,j,k,:,b_send)
                  q_gpu(iii+1,jjj+1,kkk  ,:,b_recv) = q_gpu(i,j,k,:,b_send)
                  q_gpu(iii,  jjj,  kkk+1,:,b_recv) = q_gpu(i,j,k,:,b_send)
                  q_gpu(iii+1,jjj,  kkk+1,:,b_recv) = q_gpu(i,j,k,:,b_send)
                  q_gpu(iii,  jjj+1,kkk+1,:,b_recv) = q_gpu(i,j,k,:,b_send)
                  q_gpu(iii+1,jjj+1,kkk+1,:,b_recv) = q_gpu(i,j,k,:,b_send)
               enddo
            enddo
         enddo
      endif
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine update_ghost_local_gpu_cuf

   subroutine update_ghost_mpi_gpu_cuf(gci, gcj, gck, nv, procs_number, req_send_recv,   &
                                       comm_map_send_ptr_ghost, comm_map_recv_ptr_ghost, &
                                       comm_map_recv_ghost_gpu, comm_map_send_ghost_gpu, &
                                       recv_buffer_ghost_gpu, send_buffer_ghost_gpu , q_gpu, step)
   !< Update ghost cells within other processes.
   integer(I4P), intent(in)                         :: gci                          !< Ghost cells number in I direction.
   integer(I4P), intent(in)                         :: gcj                          !< Ghost cells number in J direction.
   integer(I4P), intent(in)                         :: gck                          !< Ghost cells number in K direction.
   integer(I4P),              intent(in)            :: nv                           !< Number of variables of q.
   integer(I4P),              intent(in)            :: procs_number                 !< Number of MPI processes.
   integer(I4P),              intent(inout)         :: req_send_recv(:)             !< MPI request receive flags.
   integer(I4P), allocatable, intent(in)            :: comm_map_send_ptr_ghost(:)   !< Communication map, pointers list to send.
   integer(I4P), allocatable, intent(in)            :: comm_map_recv_ptr_ghost(:)   !< Communication map, pointers list to recv.
   integer(I8P), allocatable, intent(in),    device :: comm_map_recv_ghost_gpu(:,:) !< Communication map, `fec` information.
   integer(I8P), allocatable, intent(in),    device :: comm_map_send_ghost_gpu(:,:) !< Communication map, `fec` information.
   real(R8P),                 intent(inout), device :: recv_buffer_ghost_gpu(:)     !< Receive buffer of ghost cells.
   real(R8P),                 intent(inout), device :: send_buffer_ghost_gpu(:)     !< Send buffer of ghost cells.
   real(R8P),                 intent(inout), device :: q_gpu(1-gci:,&
                                                             1-gcj:,&
                                                             1-gck:,&
                                                             1:,1:)                 !< Field component to be updated.
   integer(I4P),              intent(in), optional  :: step                         !< Step to be perfordmed in asyncronous comp.
   logical                                          :: do_step(3)                   !< Steps to be performed in asyncronous comp.
   integer(I4P)                                     :: i, j, k                      !< Counter.
   integer(I4P)                                     :: iii, jjj, kkk                !< Counter.
   integer(I4P)                                     :: fec, mf, rf, sf, n, p, v     !< Counter.
   integer(I4P)                                     :: portion                      !< Portion of fec updated (0=>whole fec).
   integer(I4P)                                     :: b_recv                       !< Index of receiving block.
   integer(I4P)                                     :: b_send                       !< Index of sending block.
   integer(I4P)                                     :: imin                         !< Lower limit of i indexes.
   integer(I4P)                                     :: jmin                         !< Lower limit of j indexes.
   integer(I4P)                                     :: kmin                         !< Lower limit of j indexes.
   integer(I4P)                                     :: imax                         !< Upper limit of i indexes.
   integer(I4P)                                     :: jmax                         !< Upper limit of j indexes.
   integer(I4P)                                     :: kmax                         !< Upper limit of k indexes.
   integer(I4P)                                     :: idelta                       !< Delta offset for ghost-inner cells of i.
   integer(I4P)                                     :: jdelta                       !< Delta offset for ghost-inner cells of j.
   integer(I4P)                                     :: kdelta                       !< Delta offset for ghost-inner cells of k.
   integer(I4P)                                     :: ptr_start, ptr_end           !< Counter.
   integer(I4P)                                     :: n_recv, n_send               !< Counter.
   integer(I4P)                                     :: recv_rank                    !< Rank of receiving block.
   integer(I4P)                                     :: send_rank                    !< Rank of sending block.
   integer(I4P)                                     :: error                        !< Error traping flag.
   integer(I4P)                                     :: send_ptr, send_ctr           !< Counter.
   integer(I4P)                                     :: recv_ptr, recv_ctr           !< Counter.
   integer(I4P)                                     :: iercuda                      !< Error trapping flag for CUDAFortran.

   if ((.not.allocated(comm_map_recv_ghost_gpu)).and.&
       (.not.allocated(comm_map_send_ghost_gpu))) return

   do_step = .true.
   if (present(step)) then
      do_step = .false.
      do_step(step) = .true.
   endif

   if (do_step(1)) then
      req_send_recv = MPI_REQUEST_NULL

      ! populate send buffer
      !$cuf kernel do(1) <<<*,*>>>
      do sf=1, size(comm_map_send_ghost_gpu, dim=1)
         ! b_ghost   =     comm_map_send_ghost(sf, 1) ! block-index
         b_send    = comm_map_send_ghost_gpu(sf, 2) ! neighbor-block-index of block
         send_rank = comm_map_send_ghost_gpu(sf, 3)
         fec       = comm_map_send_ghost_gpu(sf, 4)
         portion   = comm_map_send_ghost_gpu(sf, 5)
         imin      = comm_map_send_ghost_gpu(sf, 6 )
         jmin      = comm_map_send_ghost_gpu(sf, 7 )
         kmin      = comm_map_send_ghost_gpu(sf, 8 )
         imax      = comm_map_send_ghost_gpu(sf, 9 )
         jmax      = comm_map_send_ghost_gpu(sf, 10)
         kmax      = comm_map_send_ghost_gpu(sf, 11)
         idelta    = comm_map_send_ghost_gpu(sf, 12)
         jdelta    = comm_map_send_ghost_gpu(sf, 13)
         kdelta    = comm_map_send_ghost_gpu(sf, 14)
         send_ptr  = comm_map_send_ghost_gpu(sf, 15)
         if (portion==0_I4P) then
            ! sending to a block at my level
            send_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do v=1,nv
                        send_buffer_ghost_gpu(send_ptr + send_ctr) = q_gpu(i+idelta,j+jdelta,k+kdelta,v,b_send)
                        send_ctr = send_ctr + 1
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
                        do v=1,nv
                           send_buffer_ghost_gpu(send_ptr + send_ctr) = q_gpu(i,j,k,v,b_send)
                           send_ctr = send_ctr + 1
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
                     do v=1,nv
                        send_buffer_ghost_gpu(send_ptr + send_ctr) =                               &
                           (q_gpu(iii,jjj,  kkk,  v,b_send) + q_gpu(iii+1,jjj,  kkk,  v,b_send) +  &
                            q_gpu(iii,jjj+1,kkk,  v,b_send) + q_gpu(iii+1,jjj+1,kkk,  v,b_send) +  &
                            q_gpu(iii,jjj,  kkk+1,v,b_send) + q_gpu(iii+1,jjj,  kkk+1,v,b_send) +  &
                            q_gpu(iii,jjj+1,kkk+1,v,b_send) + q_gpu(iii+1,jjj+1,kkk+1,v,b_send)) / 8._R8P
                        send_ctr = send_ctr + 1
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
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
      do rf=1, size(comm_map_recv_ghost_gpu, dim=1)
         b_recv    = comm_map_recv_ghost_gpu(rf, 1 ) ! block-index
         recv_rank = comm_map_recv_ghost_gpu(rf, 3 )
         fec       = comm_map_recv_ghost_gpu(rf, 4 )
         portion   = comm_map_recv_ghost_gpu(rf, 5 )
         imin      = comm_map_recv_ghost_gpu(rf, 6 )
         jmin      = comm_map_recv_ghost_gpu(rf, 7 )
         kmin      = comm_map_recv_ghost_gpu(rf, 8 )
         imax      = comm_map_recv_ghost_gpu(rf, 9 )
         jmax      = comm_map_recv_ghost_gpu(rf, 10)
         kmax      = comm_map_recv_ghost_gpu(rf, 11)
         idelta    = comm_map_recv_ghost_gpu(rf, 12)
         jdelta    = comm_map_recv_ghost_gpu(rf, 13)
         kdelta    = comm_map_recv_ghost_gpu(rf, 14)
         recv_ptr  = comm_map_recv_ghost_gpu(rf, 15)
         if (portion==0_I4P) then
            recv_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do v=1, nv
                        q_gpu(i,j,k,v,b_recv) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                        recv_ctr = recv_ctr + 1
                     enddo
                  enddo
               enddo
            enddo
         elseif (portion>0_I4P) then
            ! receiving from a block finer than me
            recv_ctr = 1
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     do v=1, nv
                        q_gpu(i,j,k,v,b_recv) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr)
                        recv_ctr = recv_ctr + 1
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
                     do v=1, nv
                        q_gpu(iii,  jjj,  kkk  ,v,b_recv) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
                     enddo
                     do v=1, nv
                        q_gpu(iii+1,jjj,  kkk  ,v,b_recv) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
                     enddo
                     do v=1, nv
                        q_gpu(iii,  jjj+1,kkk  ,v,b_recv) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
                     enddo
                     do v=1, nv
                        q_gpu(iii+1,jjj+1,kkk  ,v,b_recv) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
                     enddo
                     do v=1, nv
                        q_gpu(iii,  jjj,  kkk+1,v,b_recv) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
                     enddo
                     do v=1, nv
                        q_gpu(iii+1,jjj,  kkk+1,v,b_recv) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
                     enddo
                     do v=1, nv
                        q_gpu(iii,  jjj+1,kkk+1,v,b_recv) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
                     enddo
                     do v=1, nv
                        q_gpu(iii+1,jjj+1,kkk+1,v,b_recv) = recv_buffer_ghost_gpu(recv_ptr + recv_ctr) ; recv_ctr = recv_ctr + 1
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endsubroutine update_ghost_mpi_gpu_cuf

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
