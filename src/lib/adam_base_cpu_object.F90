!< ADAM, base CPU class definition.
module adam_base_cpu_object
!< ADAM, base CPU class definition: provide update ghosts methods for CPU backend.

use adam_field_object, only : field_object
use PENF
use MPI

implicit none
private
public :: base_cpu_object

type :: base_cpu_object
   !< Equation base CPU class definition.
   !<
   !< Provide update ghosts methods for CPU backend.
   type(field_object), pointer :: field=>null() !< The field.
   ! MPI data, unrelated to field equations
   integer(I4P) :: error=0_I4P  !< Error traping flag.
   integer(I4P) :: myrank=0_I4P !< MPI rank process.
   contains
      ! public methods
      procedure, pass(self) :: destroy            !< Destroy the equation.
      procedure, pass(self) :: initialize         !< Initialize the equation.
      procedure, pass(self) :: update_ghost_local !< Update ghosts locally.
      procedure, pass(self) :: update_ghost_mpi   !< Update ghosts MPI.
      ! operators
      generic :: assignment(=) => eq_assign_eq      !< Overload `=`.
      procedure, pass(lhs), private :: eq_assign_eq !< Operator `=`.
endtype base_cpu_object

contains
   ! public methods
   subroutine destroy(self)
   !< Destroy the equation.
   class(base_cpu_object), intent(inout) :: self  !< The equation.
   type(base_cpu_object)                 :: fresh !< Fresh equation.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, field)
   !< Initialize the equation.
   class(base_cpu_object), intent(inout)      :: self  !< The equation.
   type(field_object),     intent(in), target :: field !< The field.

   call self%destroy
   self%field => field
   call MPI_COMM_RANK(MPI_COMM_WORLD, self%myrank, self%error)
   endsubroutine initialize

   subroutine update_ghost_local(self, q)
   !< Update (local) ghost cells, rank 4.
   class(base_cpu_object), intent(inout) :: self                         !< The equation.
   real(R8P),              intent(inout) :: q(1:,                    &
                                              1-self%field%grid%ngc:,&
                                              1-self%field%grid%ngc:,&
                                              1-self%field%grid%ngc:,1:) !< Field component to be updated.
   integer(I4P)                          :: i, j, k, mf                  !< Counter.
   integer(I4P)                          :: iii, jjj, kkk                !< Counter.
   integer(I4P)                          :: fec                          !< Ghost direction, faces/edges/corners.
   integer(I4P)                          :: portion                      !< Portion of fec updated (0=>whole fec).
   integer(I4P)                          :: b_recv                       !< Index of receiving block.
   integer(I4P)                          :: b_send                       !< Index of sending block.
   integer(I4P)                          :: ijkmin(3)                    !< Lower limit of ijk indexes.
   integer(I4P)                          :: ijkmax(3)                    !< Upper limit of ijk indexes.
   integer(I4P)                          :: ijkdelta(3)                  !< Delta offset for ghost-inner cells mapping.

   if (.not.allocated(self%field%local_map_ghost)) return
   do mf=1, size(self%field%local_map_ghost, dim=1)
      b_recv   = self%field%local_map_ghost(mf, 1)
      b_send   = self%field%local_map_ghost(mf, 2)
      fec      = self%field%local_map_ghost(mf, 3)
      portion  = self%field%local_map_ghost(mf, 4)
      ijkmin   = self%field%local_map_ghost(mf, 5:7)
      ijkmax   = self%field%local_map_ghost(mf, 8:10)
      ijkdelta = self%field%local_map_ghost(mf, 11:13)
      if     (portion==0) then
         ! receiving from a block with the same refinement
         do k=ijkmin(3), ijkmax(3)
            do j=ijkmin(2), ijkmax(2)
               do i=ijkmin(1), ijkmax(1)
                  q(:,i,j,k,b_recv) = q(:,i+ijkdelta(1),j+ijkdelta(2),k+ijkdelta(3),b_send)
               enddo
            enddo
         enddo
      elseif (portion>0) then
         ! receiving from a block finer than me
         do k=ijkmin(3), ijkmax(3)
            do j=ijkmin(2), ijkmax(2)
               do i=ijkmin(1), ijkmax(1)
                  kkk = 2 * k + ijkdelta(3)
                  jjj = 2 * j + ijkdelta(2)
                  iii = 2 * i + ijkdelta(1)
                  q(:,i,j,k,b_recv) = (q(:,iii,jjj,  kkk,  b_send) + q(:,iii+1,jjj,  kkk,  b_send) + &
                                       q(:,iii,jjj+1,kkk,  b_send) + q(:,iii+1,jjj+1,kkk,  b_send) + &
                                       q(:,iii,jjj,  kkk+1,b_send) + q(:,iii+1,jjj,  kkk+1,b_send) + &
                                       q(:,iii,jjj+1,kkk+1,b_send) + q(:,iii+1,jjj+1,kkk+1,b_send)) / 8._R8P
               enddo
            enddo
         enddo
      else
         ! receiving from a block coarser than me
         do k=ijkmin(3), ijkmax(3)
            do j=ijkmin(2), ijkmax(2)
               do i=ijkmin(1), ijkmax(1)
                  kkk = 2 * k + ijkdelta(3)
                  jjj = 2 * j + ijkdelta(2)
                  iii = 2 * i + ijkdelta(1)
                  q(:,iii,  jjj,  kkk  ,b_recv) = q(:,i,j,k,b_send)
                  q(:,iii+1,jjj,  kkk  ,b_recv) = q(:,i,j,k,b_send)
                  q(:,iii,  jjj+1,kkk  ,b_recv) = q(:,i,j,k,b_send)
                  q(:,iii+1,jjj+1,kkk  ,b_recv) = q(:,i,j,k,b_send)
                  q(:,iii,  jjj,  kkk+1,b_recv) = q(:,i,j,k,b_send)
                  q(:,iii+1,jjj,  kkk+1,b_recv) = q(:,i,j,k,b_send)
                  q(:,iii,  jjj+1,kkk+1,b_recv) = q(:,i,j,k,b_send)
                  q(:,iii+1,jjj+1,kkk+1,b_recv) = q(:,i,j,k,b_send)
               enddo
            enddo
         enddo
      endif
   enddo
   endsubroutine update_ghost_local

   subroutine update_ghost_mpi(self, q, step)
   !< Update ghost cells within other processes.
   class(base_cpu_object), intent(inout)        :: self                       !< The field.
   real(R8P),              intent(inout)        :: q(1:,                    &
                                                     1-self%field%grid%ngc:,&
                                                     1-self%field%grid%ngc:,&
                                                     1-self%field%grid%ngc:,&
                                                     1:)                      !< Field component to be updated.
   integer(I4P),           intent(in), optional :: step                       !< Step to be perfordmed in asyncronous comp.
   logical                                      :: do_step(3)                 !< Steps to be performed in asyncronous comp.
   integer(I4P)                                 :: i, j, k                    !< Counter.
   integer(I4P)                                 :: iii, jjj, kkk              !< Counter.
   integer(I4P)                                 :: fec, mf, rf, sf, n, p, v   !< Counter.
   integer(I4P), allocatable                    :: comm_map_send_ctr_ghost(:) !< Comm map, cts to send [procs_number+1].
   integer(I4P), allocatable                    :: comm_map_recv_ctr_ghost(:) !< Comm map, cts to recv [procs_number+1].
   integer(I4P)                                 :: portion                    !< Portion of fec updated (0=>whole fec).
   integer(I4P)                                 :: b_recv                     !< Index of receiving block.
   integer(I4P)                                 :: b_send                     !< Index of sending block.
   integer(I4P)                                 :: delta(3)                   !< Neighbor delta of current fec.
   integer(I4P)                                 :: ijkmin(3)                  !< Lower limit of ijk indexes.
   integer(I4P)                                 :: ijkmax(3)                  !< Upper limit of ijk indexes.
   integer(I4P)                                 :: ijkdelta(3)                !< Delta offset of ghost-inner cells map.
   integer(I4P)                                 :: ptr_start, ptr_end         !< Counter.
   integer(I4P)                                 :: n_recv, n_send             !< Counter.
   integer(I4P)                                 :: recv_rank                  !< Rank of receiving block.
   integer(I4P)                                 :: send_rank                  !< Rank of sending block.

   if ((.not.allocated(self%field%comm_map_recv_ghost)).and.&
       (.not.allocated(self%field%comm_map_send_ghost))) return

   do_step = .true.
   if (present(step)) then
      do_step = .false.
      do_step(step) = .true.
   endif

   if (do_step(1)) then
      self%field%req_send_recv = MPI_REQUEST_NULL
      comm_map_send_ctr_ghost = self%field%comm_map_send_ptr_ghost

      ! populate send buffer
      do sf=1, size(self%field%comm_map_send_ghost, dim=1)
         ! b_ghost   =     comm_map_send_ghost(sf, 1) ! block-index
         b_send    = self%field%comm_map_send_ghost(sf, 2) ! neighbor-block-index of block
         send_rank = self%field%comm_map_send_ghost(sf, 3)
         fec       = self%field%comm_map_send_ghost(sf, 4)
         portion   = self%field%comm_map_send_ghost(sf, 5)
         ijkmin    = self%field%comm_map_send_ghost(sf, 6:8)
         ijkmax    = self%field%comm_map_send_ghost(sf, 9:11)
         ijkdelta  = self%field%comm_map_send_ghost(sf, 12:14)
         if (portion==0_I4P) then
            ! sending to a block at my level
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     do v=1,self%field%nv
                        self%field%send_buffer_ghost(comm_map_send_ctr_ghost(send_rank)+1) = &
                           q(v,i+ijkdelta(1),j+ijkdelta(2),k+ijkdelta(3),b_send)
                        comm_map_send_ctr_ghost(send_rank) = comm_map_send_ctr_ghost(send_rank) + 1
                     enddo
                  enddo
               enddo
            enddo
         elseif (portion<0_I4P) then ! Beware! This is < 0 because the reference is the receiver
            ! sending to a block finer than me
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     do n=1,8
                        do v=1,self%field%nv
                           self%field%send_buffer_ghost(comm_map_send_ctr_ghost(send_rank)+1) = &
                              q(v,i,j,k,b_send)
                           comm_map_send_ctr_ghost(send_rank) = comm_map_send_ctr_ghost(send_rank) + 1
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         else
            ! sending to a block coarser than me, loop is over the coarser grid
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     kkk = 2 * k + ijkdelta(3)
                     jjj = 2 * j + ijkdelta(2)
                     iii = 2 * i + ijkdelta(1)
                     do v=1,self%field%nv
                        self%field%send_buffer_ghost(comm_map_send_ctr_ghost(send_rank)+1) = &
                           (q(v,iii,jjj,  kkk,  b_send) + q(v,iii+1,jjj,  kkk,  b_send) +  &
                            q(v,iii,jjj+1,kkk,  b_send) + q(v,iii+1,jjj+1,kkk,  b_send) +  &
                            q(v,iii,jjj,  kkk+1,b_send) + q(v,iii+1,jjj,  kkk+1,b_send) +  &
                            q(v,iii,jjj+1,kkk+1,b_send) + q(v,iii+1,jjj+1,kkk+1,b_send)) / 8._R8P
                        comm_map_send_ctr_ghost(send_rank) = comm_map_send_ctr_ghost(send_rank) + 1
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
   endif

   if (do_step(2)) then
      ! receive
      do p=0, self%field%procs_number - 1_I4P
         ptr_start = self%field%comm_map_recv_ptr_ghost(p) + 1
         ptr_end   = self%field%comm_map_recv_ptr_ghost(p+1)
         n_recv    = ptr_end - ptr_start + 1
         if (n_recv > 0) then
#ifdef _MPI_
            call MPI_IRECV(self%field%recv_buffer_ghost(ptr_start), n_recv, MPI_REAL8, p, 100, MPI_COMM_WORLD, &
                           self%field%req_send_recv(p), self%field%error)
#endif
         endif
      enddo

      ! send
      do p=0, self%field%procs_number - 1_I4P
         ptr_start = self%field%comm_map_send_ptr_ghost(p) + 1
         ptr_end   = self%field%comm_map_send_ptr_ghost(p+1)
         n_send    = ptr_end - ptr_start + 1
         if (n_send > 0) then
#ifdef _MPI_
            call MPI_ISEND(self%field%send_buffer_ghost(ptr_start), n_send, MPI_REAL8, p, 100, MPI_COMM_WORLD, &
                           self%field%req_send_recv(p+self%field%procs_number), self%field%error)
#endif
         endif
      enddo
   endif

   if (do_step(3)) then
      comm_map_recv_ctr_ghost = self%field%comm_map_recv_ptr_ghost

      call MPI_WAITALL(self%field%procs_number * 2, self%field%req_send_recv, MPI_STATUSES_IGNORE, self%field%error)

      ! retrive from receive buffer
      do rf=1, size(self%field%comm_map_recv_ghost, dim=1)
         b_recv    = self%field%comm_map_recv_ghost(rf, 1) ! block-index
         ! b_recv    =     comm_map_recv_ghost(rf, 2) ! neighbor-block-index of block
         recv_rank = self%field%comm_map_recv_ghost(rf, 3)
         fec       = self%field%comm_map_recv_ghost(rf, 4)
         portion   = self%field%comm_map_recv_ghost(rf, 5)
         ijkmin    = self%field%comm_map_recv_ghost(rf, 6:8)
         ijkmax    = self%field%comm_map_recv_ghost(rf, 9:11)
         ijkdelta  = self%field%comm_map_recv_ghost(rf, 12:14)
         if (portion==0_I4P) then
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     do v=1, self%field%nv
                        q(v,i,j,k,b_recv) = self%field%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                        comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     enddo
                  enddo
               enddo
            enddo
         elseif (portion>0_I4P) then
            ! receiving from a block finer than me
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     do v=1, self%field%nv
                        q(v,i,j,k,b_recv) = self%field%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                        comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     enddo
                  enddo
               enddo
            enddo
         else
            ! receiving from a block coarser than me
            do k=ijkmin(3), ijkmax(3)
               do j=ijkmin(2), ijkmax(2)
                  do i=ijkmin(1), ijkmax(1)
                     kkk = 2 * k + ijkdelta(3)
                     jjj = 2 * j + ijkdelta(2)
                     iii = 2 * i + ijkdelta(1)
                     do v=1, self%field%nv
                        q(v,iii,  jjj,  kkk  ,b_recv) = self%field%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                        comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     enddo
                     do v=1, self%field%nv
                        q(v,iii+1,jjj,  kkk  ,b_recv) = self%field%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                        comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     enddo
                     do v=1, self%field%nv
                        q(v,iii,  jjj+1,kkk  ,b_recv) = self%field%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                        comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     enddo
                     do v=1, self%field%nv
                        q(v,iii+1,jjj+1,kkk  ,b_recv) = self%field%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                        comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     enddo
                     do v=1, self%field%nv
                        q(v,iii,  jjj,  kkk+1,b_recv) = self%field%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                        comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     enddo
                     do v=1, self%field%nv
                        q(v,iii+1,jjj,  kkk+1,b_recv) = self%field%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                        comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     enddo
                     do v=1, self%field%nv
                        q(v,iii,  jjj+1,kkk+1,b_recv) = self%field%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                        comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     enddo
                     do v=1, self%field%nv
                        q(v,iii+1,jjj+1,kkk+1,b_recv) = self%field%recv_buffer_ghost(comm_map_recv_ctr_ghost(recv_rank)+1)
                        comm_map_recv_ctr_ghost(recv_rank) = comm_map_recv_ctr_ghost(recv_rank) + 1
                     enddo
                  enddo
               enddo
            enddo
         endif
      enddo
   endif
   endsubroutine update_ghost_mpi

   ! operators
   ! =
   subroutine eq_assign_eq(lhs, rhs)
   !< Operator `=`.
   class(base_cpu_object), intent(inout) :: lhs !< Left hand side.
   type(base_cpu_object),  intent(in)    :: rhs !< Right hand side.

   lhs%field => rhs%field
   lhs%error  = rhs%error
   lhs%myrank = rhs%myrank
   endsubroutine eq_assign_eq
endmodule adam_base_cpu_object
