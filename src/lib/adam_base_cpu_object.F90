!< ADAM, base CPU class definition.
module adam_base_cpu_object
!< ADAM, base CPU class definition: provide methods for CPU backend handling.

use adam_field_object, only : field_object
use adam_mpih_object,  only : mpih_object
use adam_tree_object,  only : tree_object
use adam_memory_cpu_lib
use penf
use mpi
use, intrinsic :: iso_c_binding, only : C_LONG

implicit none
save
private
public :: base_cpu_object

type :: base_cpu_object
   !< Equation base CPU class definition.
   !<
   !< Provide update ghosts methods for CPU backend.
   type(mpih_object)           :: mpih                !< MPI handler.
   type(tree_object),  pointer :: tree=>null()        !< The tree.
   type(field_object), pointer :: field=>null()       !< The field.
   real(R8P)                   :: memory_avail=0._R8P !< CPU memory available (Gb).
   integer(I4P), allocatable   :: req_send_recv(:)    !< MPI request receive flags.
   contains
      ! public methods
      procedure, pass(self) :: initialize         !< Initialize the equation.
      procedure, pass(self) :: initialize_cpu     !< Initialize CPU main data.
      procedure, pass(self) :: update_ghost_local !< Update ghosts locally.
      procedure, pass(self) :: update_ghost_mpi   !< Update ghosts MPI.
endtype base_cpu_object

contains
   ! public methods
   subroutine initialize(self, tree, field)
   !< Initialize the equation.
   class(base_cpu_object), intent(inout)      :: self  !< The equation.
   type(tree_object),      intent(in), target :: tree  !< The tree.
   type(field_object),     intent(in), target :: field !< The field.

   print '(A)', self%mpih%myrankstr//'base_cpu%initialize start'
   self%tree  => tree
   self%field => field
   allocate(self%req_send_recv(0:self%mpih%procs_number*2-1))
   print '(A)', self%mpih%myrankstr//'base_cpu%initialize finish'
   endsubroutine initialize

   subroutine initialize_cpu(self)
   !< Initialize CPU main data.
   !< @Note This must be the first routine called.
   class(base_cpu_object), intent(inout) :: self !< The base backend.
   integer(C_LONG)                       :: mem_free, mem_total !< CPU memory.

   call self%mpih%initialize
   print '(A)', self%mpih%myrankstr//'base_cpu%initialize_cpu start'
   call cpuMemGetInfo(mem_total, mem_free)
   self%memory_avail = real(mem_total, R8P)/1e9
   print '(A)', self%mpih%myrankstr//'base_cpu%initialize_cpu avail memory: '//trim(str(fm='(F4.1)',n=self%memory_avail))//' [GB]'
   print '(A)', self%mpih%myrankstr//'base_cpu%initialize_cpu finish'
   endsubroutine initialize_cpu

   subroutine update_ghost_local(self, q)
   !< Update (local) ghost cells, rank 4.
   class(base_cpu_object), intent(inout) :: self                         !< The equation.
   real(R8P),              intent(inout) :: q(1:,                    &
                                              1-self%field%grid%ngc:,&
                                              1-self%field%grid%ngc:,&
                                              1-self%field%grid%ngc:,1:) !< Field component to be updated.
   integer(I4P)                          :: mf                           !< Counter.
   integer(I4P)                          :: fec                          !< Ghost direction, faces/edges/corners.
   integer(I4P)                          :: b_recv                       !< Index of receiving block.
   integer(I4P)                          :: b_send                       !< Index of sending block.
   integer(I4P)                          :: ic, jc, kc                   !< Counter.
   integer(I4P)                          :: i_recv                       !< I recv index.
   integer(I4P)                          :: j_recv                       !< J recv index.
   integer(I4P)                          :: k_recv                       !< K recv index.
   integer(I4P)                          :: i_send                       !< I send index.
   integer(I4P)                          :: j_send                       !< J send index.
   integer(I4P)                          :: k_send                       !< K send index.
   integer(I4P)                          :: one_or_eight                 !< Flag triggering 8 cells mean.

   if (.not.allocated(self%tree%local_map_ghost_cell)) return
   do mf=1, size(self%tree%local_map_ghost_cell, dim=1)
      b_send       = self%tree%local_map_ghost_cell(mf,1)
      b_recv       = self%tree%local_map_ghost_cell(mf,2)
      i_send       = self%tree%local_map_ghost_cell(mf,3)
      j_send       = self%tree%local_map_ghost_cell(mf,4)
      k_send       = self%tree%local_map_ghost_cell(mf,5)
      i_recv       = self%tree%local_map_ghost_cell(mf,6)
      j_recv       = self%tree%local_map_ghost_cell(mf,7)
      k_recv       = self%tree%local_map_ghost_cell(mf,8)
      one_or_eight = self%tree%local_map_ghost_cell(mf,9)
      if (one_or_eight==1) then
         q(:,i_recv,j_recv,k_recv,b_recv) = q(:,i_send,j_send,k_send,b_send)
      else
         q(:,i_recv,j_recv,k_recv,b_recv) = 0._R8P
         do kc=0,1 ; do jc=0,1 ; do ic=0,1
            q(:,i_recv,j_recv,k_recv,b_recv) = q(:,i_recv,   j_recv,   k_recv,   b_recv) + &
                                               q(:,i_send+ic,j_send+jc,k_send+kc,b_send)
         enddo ; enddo ; enddo
         q(:,i_recv,j_recv,k_recv,b_recv) = q(:,i_recv,j_recv,k_recv,b_recv) * 0.125_R8P
      endif
   enddo
   endsubroutine update_ghost_local

   subroutine update_ghost_mpi(self, q, step)
   !< Update ghost cells within other processes.
   class(base_cpu_object), intent(inout)        :: self                                   !< The field.
   real(R8P),              intent(inout)        :: q(1:,                    &
                                                     1-self%field%grid%ngc:,&
                                                     1-self%field%grid%ngc:,&
                                                     1-self%field%grid%ngc:,&
                                                     1:)                                  !< Field component to be updated.
   integer(I4P),           intent(in), optional :: step                                   !< Step to be perfordmed in async.
   logical                                      :: do_step(3)                             !< Steps to be performed in async.
   integer(I4P)                                 :: ic, jc, kc                             !< Counter.
   integer(I4P)                                 :: rf, sf, p                              !< Counter.
   integer(I4P)                                 :: b_send, i_send, j_send, k_send, v_send !< Send indexes.
   integer(I4P)                                 :: b_recv, i_recv, j_recv, k_recv, v_recv !< Receive indexes.
   integer(I4P)                                 :: c_send, c_recv                         !< Counter.
   integer(I4P)                                 :: one_or_eight                           !< Flag triggering 8 cells mean.
   integer(I4P)                                 :: ptr_start, ptr_end                     !< Counter.
   integer(I4P)                                 :: n_recv, n_send                         !< Counter.

   associate(comm_map_send_ptr_ghost=>self%tree%comm_map_send_ptr_ghost, &
             comm_map_recv_ptr_ghost=>self%tree%comm_map_recv_ptr_ghost, &
             send_buffer_ghost=>self%tree%send_buffer_ghost,             &
             recv_buffer_ghost=>self%tree%recv_buffer_ghost,             &
             req_send_recv=>self%req_send_recv,                          &
             procs_number=>self%mpih%procs_number,                       &
             error=>self%mpih%error)
   do_step = .true.
   if (present(step)) then
      do_step = .false.
      do_step(step) = .true.
   endif

   if (do_step(1)) then
      req_send_recv = MPI_REQUEST_NULL
      ! populate send buffer
      if (allocated(self%tree%comm_map_send_ghost_cell)) then
         do sf=1, size(self%tree%comm_map_send_ghost_cell, dim=1)
            b_send       = self%tree%comm_map_send_ghost_cell(sf,1)
            i_send       = self%tree%comm_map_send_ghost_cell(sf,2)
            j_send       = self%tree%comm_map_send_ghost_cell(sf,3)
            k_send       = self%tree%comm_map_send_ghost_cell(sf,4)
            v_send       = self%tree%comm_map_send_ghost_cell(sf,5)
            c_recv       = self%tree%comm_map_send_ghost_cell(sf,6)
            one_or_eight = self%tree%comm_map_send_ghost_cell(sf,7)
            if (one_or_eight==1) then
               send_buffer_ghost(c_recv) = q(v_send,i_send,j_send,k_send,b_send)
            else
               send_buffer_ghost(c_recv) = 0._R8P
               do kc=0,1 ; do jc=0,1 ; do ic=0,1
                  send_buffer_ghost(c_recv) = send_buffer_ghost(c_recv) + &
                                              q(v_send,i_send+ic,j_send+jc,k_send+kc,b_send)
               enddo ; enddo ; enddo
               send_buffer_ghost(c_recv) = send_buffer_ghost(c_recv) * 0.125_R8P
            endif
         enddo
      endif
   endif

   if (do_step(2)) then
      ! receive
      do p=0, procs_number - 1_I4P
         ptr_start = comm_map_recv_ptr_ghost(p) + 1
         ptr_end   = comm_map_recv_ptr_ghost(p+1)
         n_recv    = ptr_end - ptr_start + 1
         if (n_recv > 0) then
            call MPI_IRECV(recv_buffer_ghost(ptr_start), n_recv, MPI_REAL8, p, 100, MPI_COMM_WORLD, req_send_recv(p), error)
         endif
      enddo
      ! send
      do p=0, procs_number - 1_I4P
         ptr_start = comm_map_send_ptr_ghost(p) + 1
         ptr_end   = comm_map_send_ptr_ghost(p+1)
         n_send    = ptr_end - ptr_start + 1
         if (n_send > 0) then
            call MPI_ISEND(send_buffer_ghost(ptr_start), n_send, MPI_REAL8, p, 100, MPI_COMM_WORLD, &
                           req_send_recv(p+procs_number), error)
         endif
      enddo
   endif

   if (do_step(3)) then
      call MPI_WAITALL(procs_number * 2, req_send_recv, MPI_STATUSES_IGNORE, error)

      call self%mpih%barrier
      !RIMETTERE SENZA

      if (allocated(self%tree%comm_map_recv_ghost_cell)) then
         ! retrive from receive buffer
         do rf=1, size(self%tree%comm_map_recv_ghost_cell, dim=1)
            c_send = self%tree%comm_map_recv_ghost_cell(rf,1)
            b_recv = self%tree%comm_map_recv_ghost_cell(rf,2)
            i_recv = self%tree%comm_map_recv_ghost_cell(rf,3)
            j_recv = self%tree%comm_map_recv_ghost_cell(rf,4)
            k_recv = self%tree%comm_map_recv_ghost_cell(rf,5)
            v_recv = self%tree%comm_map_recv_ghost_cell(rf,6)
            q(v_recv,i_recv,j_recv,k_recv,b_recv) = recv_buffer_ghost(c_send)
         enddo
      endif
   endif
   call self%mpih%barrier
   endassociate
   endsubroutine update_ghost_mpi
endmodule adam_base_cpu_object
