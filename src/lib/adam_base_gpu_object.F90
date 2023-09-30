!< ADAM, base GPU class definition.
module adam_base_gpu_object
!< ADAM, base GPU class definition: provide methods for GPU backend handling.

use adam_field_object, only : field_object
use adam_mpih_object, only : mpih_object
use adam_memory_gpu_lib
use adam_parameters
use PENF
use MPI
use CUDAFOR
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
   integer(I8P), allocatable, device :: local_map_ghost_gpu(:,:)             !< Local map for ghost cells updating, fecs order.
   integer(I8P), allocatable, device :: comm_map_recv_ghost_gpu(:,:)         !< Communication map, `fec` information.
   integer(I8P), allocatable, device :: comm_map_send_ghost_gpu(:,:)         !< Communication map, `fec` information.
   integer(I8P), allocatable, device :: comm_map_recv_ghost_cell_gpu(:,:)    !< Communication map, `fec` information, cell order.
   integer(I8P), allocatable, device :: comm_map_send_ghost_cell_gpu(:,:)    !< Communication map, `fec` information, cell order.
   real(R8P),    allocatable, device :: send_buffer_ghost_gpu(:)             !< Send buffer of ghost cells.
   real(R8P),    allocatable, device :: recv_buffer_ghost_gpu(:)             !< Receive buffer of ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_face_gpu(:,:)           !< Local map for face BC ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_edge_gpu(:,:)           !< Local map for edge BC ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_corner_gpu(:,:)         !< Local map for corner BC ghost cells.
   integer(I8P), allocatable, device :: local_map_bc_crown_gpu(:,:,:)        !< Local map for face BC ghost cells, "crown" order.
   integer(I4P), allocatable, device :: fec_1_6_array_gpu(:)                 !< Mapping fec1-26 to fec1-6 for boundaries (GPU).
   real(R8P),    allocatable, device :: x_cell_gpu(:,:)                      !< Cells x coordinates on GPU.
   real(R8P),    allocatable, device :: y_cell_gpu(:,:)                      !< Cells y coordinates on GPU.
   real(R8P),    allocatable, device :: z_cell_gpu(:,:)                      !< Cells z coordinates on GPU.
   real(R8P),    allocatable, device :: dxyz_gpu(:,:)                        !< Delta cells GPU.
   contains
      ! public methods
      procedure, pass(self) :: copy_cpu_gpu                  !< Copy data from (field) CPU to (base_gpu) GPU.
      procedure, pass(self) :: copy_transpose_cpu_gpu        !< Transpose data from GPU to CPU.
      procedure, pass(self) :: copy_transpose_gpu_cpu        !< Transpose data from GPU to CPU.
      procedure, pass(self) :: create_maps_cell              !< Create maps in cells order form the fecs ordered ones.
      procedure, pass(self) :: create_maps_fluxes_cell       !< Create fluxes maps in cells order form the fecs ordered ones.
      procedure, pass(self) :: initialize                    !< Initialize base backend.
      procedure, pass(self) :: initialize_gpu                !< Initialize GPU main data.
      procedure, pass(self) :: update_ghost_local_gpu        !< Update ghosts locally.
      procedure, pass(self) :: update_ghost_fluxes_local_gpu !< Update ghosts locally.
      procedure, pass(self) :: update_ghost_mpi_gpu          !< Update ghosts MPI.
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
   call assign_allocatable_gpu(lhs=      self%local_map_ghost_gpu, &
                               rhs=self%field%local_map_ghost,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(local_map_ghost_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%comm_map_recv_ghost_gpu, &
                               rhs=self%field%comm_map_recv_ghost,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(comm_map_recv_ghost_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%comm_map_send_ghost_gpu, &
                               rhs=self%field%comm_map_send_ghost,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(comm_map_send_ghost_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%send_buffer_ghost_gpu, &
                               rhs=self%field%send_buffer_ghost,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(send_buffer_ghost_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%recv_buffer_ghost_gpu, &
                               rhs=self%field%recv_buffer_ghost,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(recv_buffer_ghost_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%local_map_bc_face_gpu, &
                               rhs=self%field%local_map_bc_face,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(local_map_bc_face_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%local_map_bc_corner_gpu, &
                               rhs=self%field%local_map_bc_corner,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(local_map_bc_corner_gpu) ', verbose=verbose)
   call assign_allocatable_gpu(lhs=      self%local_map_bc_edge_gpu, &
                               rhs=self%field%local_map_bc_edge,     &
                               msg=self%mpih%myrankstr//'base_gpu%copy_cpu_gpu(local_map_bc_edge_gpu) ', verbose=verbose)
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
   call self%create_maps_cell(verbose=verbose)
   call self%create_maps_fluxes_cell(verbose=verbose)
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

   subroutine create_maps_cell(self, verbose)
   !< Create maps in cells order form the fecs ordered ones.
   class(base_gpu_object), intent(inout)        :: self                          !< The base backend.
   logical,                intent(in), optional :: verbose                       !< Flag to activate verbose mode.
   logical                                      :: verbose_                      !< Flag to activate verbose mode, local var.
   integer(I8P), allocatable                    :: local_map_ghost_cell(:,:)     !< Local map ghost cells update, cells order.
   integer(I8P), allocatable                    :: comm_map_send_ghost_cell(:,:) !< MPI send map ghost cells update, cells order.
   integer(I8P), allocatable                    :: comm_map_recv_ghost_cell(:,:) !< MPI send map ghost cells update, cells order.
   integer(I8P), allocatable                    :: local_map_bc_crown(:,:,:)     !< Local map for face BC ghost cells, crown order.
   integer(I4P)                                 :: c, f, v, n                    !< Counter.
   integer(I4P)                                 :: i, j, k                       !< Counter.
   integer(I4P)                                 :: iii, jjj, kkk                 !< Counter.
   integer(I4P)                                 :: ic, jc, kc                    !< Counter.
   integer(I4P)                                 :: fec                           !< Ghost direction, faces/edges/corners.
   integer(I4P)                                 :: portion                       !< Portion of fec updated (0=>whole fec).
   integer(I4P)                                 :: b_recv                        !< Index of receiving block.
   integer(I4P)                                 :: b_send                        !< Index of sending block.
   integer(I4P)                                 :: imin                          !< Lower limit of i indexes.
   integer(I4P)                                 :: jmin                          !< Lower limit of j indexes.
   integer(I4P)                                 :: kmin                          !< Lower limit of j indexes.
   integer(I4P)                                 :: imax                          !< Upper limit of i indexes.
   integer(I4P)                                 :: jmax                          !< Upper limit of j indexes.
   integer(I4P)                                 :: kmax                          !< Upper limit of k indexes.
   integer(I4P)                                 :: idelta                        !< Delta offset for ghost-inner cells of i.
   integer(I4P)                                 :: jdelta                        !< Delta offset for ghost-inner cells of j.
   integer(I4P)                                 :: kdelta                        !< Delta offset for ghost-inner cells of k.
   integer(I4P)                                 :: recv_ptr, recv_ctr            !< Counter.
   integer(I4P)                                 :: send_ptr, send_ctr            !< Counter.
   integer(I4P), allocatable                    :: c_crown(:)                    !< Counter.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) print '(A)', self%mpih%myrankstr//'base_gpu%create_maps_cell start'

   if (allocated(self%local_map_ghost_cell_gpu)) deallocate(self%local_map_ghost_cell_gpu)
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
      call assign_allocatable_gpu(lhs=self%local_map_ghost_cell_gpu, &
                                  rhs=     local_map_ghost_cell,     &
                                  msg=self%mpih%myrankstr//'base_gpu%create_maps_cell(local_map_ghost_cell_gpu) ', verbose=verbose)
      deallocate(local_map_ghost_cell)
   endif

   if (allocated(self%comm_map_send_ghost_cell_gpu)) deallocate(self%comm_map_send_ghost_cell_gpu)
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

      ! Nv variables map
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
      call assign_allocatable_gpu(lhs=self%comm_map_send_ghost_cell_gpu,                                               &
                                  rhs=     comm_map_send_ghost_cell,                                                   &
                                  msg=self%mpih%myrankstr//'base_gpu%create_maps_cell(comm_map_send_ghost_cell_gpu) ', &
                                  verbose=verbose)
      deallocate(comm_map_send_ghost_cell)
   endif

   if (allocated(self%comm_map_recv_ghost_cell_gpu)) deallocate(self%comm_map_recv_ghost_cell_gpu)
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

      ! Nv variables map
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
      call assign_allocatable_gpu(lhs=self%comm_map_recv_ghost_cell_gpu,                                               &
                                  rhs=     comm_map_recv_ghost_cell,                                                   &
                                  msg=self%mpih%myrankstr//'base_gpu%create_maps_cell(comm_map_recv_ghost_cell_gpu) ', &
                                  verbose=verbose)
      deallocate(comm_map_recv_ghost_cell)
   endif

   if (allocated(self%local_map_bc_crown_gpu)) deallocate(self%local_map_bc_crown_gpu)
   if (allocated(self%field%local_map_bc_face  ).or.&
       allocated(self%field%local_map_bc_edge  ).or.&
       allocated(self%field%local_map_bc_corner)) then
      c = 0
      if (allocated(self%field%local_map_bc_face  )) c = c + bc_cells_number(self%field%local_map_bc_face  )
      if (allocated(self%field%local_map_bc_edge  )) c = c + bc_cells_number(self%field%local_map_bc_edge  )
      if (allocated(self%field%local_map_bc_corner)) c = c + bc_cells_number(self%field%local_map_bc_corner)
      allocate(local_map_bc_crown(1:c,1:9,1:self%field%grid%ngc))
      allocate(c_crown(1:self%field%grid%ngc))
      local_map_bc_crown = -1
      c_crown = 1
      if (allocated(self%field%local_map_bc_face  )) call populate_local_map_bc_crown(self%field%local_map_bc_face  )
      if (allocated(self%field%local_map_bc_edge  )) call populate_local_map_bc_crown(self%field%local_map_bc_edge  )
      if (allocated(self%field%local_map_bc_corner)) call populate_local_map_bc_crown(self%field%local_map_bc_corner)
      deallocate(c_crown)
      call assign_allocatable_gpu(lhs=self%local_map_bc_crown_gpu, &
                                  rhs=     local_map_bc_crown,     &
                                  msg=self%mpih%myrankstr//'base_gpu%create_maps_cell(local_map_bc_crown_gpu) ', verbose=verbose)
      deallocate(local_map_bc_crown)
   endif

   if (verbose_) print '(A)', self%mpih%myrankstr//'base_gpu%create_maps_cell finish'
   contains
      function bc_cells_number(local_map_bc) result(cells_number)
      !< Return BC cells number.
      integer(I8P), intent(in) :: local_map_bc(:,:) !< Local map for BC ghost cells.
      integer(I4P)             :: f, i, j, k        !< Counter.
      integer(I4P)             :: imin              !< Lower limit of ijk indexes.
      integer(I4P)             :: jmin              !< Lower limit of ijk indexes.
      integer(I4P)             :: kmin              !< Lower limit of ijk indexes.
      integer(I4P)             :: imax              !< Upper limit of ijk indexes.
      integer(I4P)             :: jmax              !< Upper limit of ijk indexes.
      integer(I4P)             :: kmax              !< Upper limit of ijk indexes.
      integer(I4P)             :: cells_number      !< Number of BC cells.

      cells_number = 0
      do f=1, size(local_map_bc, dim=1)
         imin    = local_map_bc(f, 3 )
         jmin    = local_map_bc(f, 4 )
         kmin    = local_map_bc(f, 5 )
         imax    = local_map_bc(f, 6 )
         jmax    = local_map_bc(f, 7 )
         kmax    = local_map_bc(f, 8 )
         do k=kmin, kmax, sign(1, kmax-kmin)
            do j=jmin, jmax, sign(1, jmax-jmin)
               do i=imin, imax, sign(1, imax-imin)
                  cells_number = cells_number + 1
               enddo
            enddo
         enddo
      enddo
      endfunction bc_cells_number

      subroutine populate_local_map_bc_crown(local_map_bc)
      !< Populate map of BC cells in crown order.
      integer(I8P), intent(in) :: local_map_bc(:,:) !< Local map for BC ghost cells.
      integer(I4P)             :: f, i, j, k, b     !< Counter.
      integer(I4P)             :: imin              !< Lower limit of ijk indexes.
      integer(I4P)             :: jmin              !< Lower limit of ijk indexes.
      integer(I4P)             :: kmin              !< Lower limit of ijk indexes.
      integer(I4P)             :: imax              !< Upper limit of ijk indexes.
      integer(I4P)             :: jmax              !< Upper limit of ijk indexes.
      integer(I4P)             :: kmax              !< Upper limit of ijk indexes.
      integer(I4P)             :: idelta            !< IJK delta step for extrapolation.
      integer(I4P)             :: jdelta            !< IJK delta step for extrapolation.
      integer(I4P)             :: kdelta            !< IJK delta step for extrapolation.
      integer(I4P)             :: bc_type           !< Boundary condition type.
      integer(I4P)             :: crown             !< Crown counter.
      integer(I4P)             :: fec               !< Boundary condition fec.

      do f=1, size(local_map_bc, dim=1)
         b       = local_map_bc(f, 1 )
         fec     = local_map_bc(f, 2 )
         imin    = local_map_bc(f, 3 )
         jmin    = local_map_bc(f, 4 )
         kmin    = local_map_bc(f, 5 )
         imax    = local_map_bc(f, 6 )
         jmax    = local_map_bc(f, 7 )
         kmax    = local_map_bc(f, 8 )
         idelta  = local_map_bc(f, 9 )
         jdelta  = local_map_bc(f, 10)
         kdelta  = local_map_bc(f, 11)
         bc_type = local_map_bc(f, 12)
         do k=kmin, kmax, sign(1, kmax-kmin)
            do j=jmin, jmax, sign(1, jmax-jmin)
               do i=imin, imax, sign(1, imax-imin)
                  crown = maxval([abs(i-imin), abs(j-jmin), abs(k-kmin)], mask=[imin/=1, jmin/=1, kmin/=1]) + 1
                  local_map_bc_crown(c_crown(crown), 1, crown) = b
                  local_map_bc_crown(c_crown(crown), 2, crown) = i
                  local_map_bc_crown(c_crown(crown), 3, crown) = j
                  local_map_bc_crown(c_crown(crown), 4, crown) = k
                  local_map_bc_crown(c_crown(crown), 5, crown) = idelta
                  local_map_bc_crown(c_crown(crown), 6, crown) = jdelta
                  local_map_bc_crown(c_crown(crown), 7, crown) = kdelta
                  local_map_bc_crown(c_crown(crown), 8, crown) = bc_type
                  local_map_bc_crown(c_crown(crown), 9, crown) = fec
                  c_crown(crown) = c_crown(crown) + 1
               enddo
            enddo
         enddo
      enddo
      endsubroutine populate_local_map_bc_crown
   endsubroutine create_maps_cell

   subroutine create_maps_fluxes_cell(self, verbose)
   !< Create maps in cells order form the fecs ordered ones.
   class(base_gpu_object), intent(inout)        :: self                             !< The base backend.
   logical,                intent(in), optional :: verbose                          !< Flag to activate verbose mode.
   logical                                      :: verbose_                         !< Flag to activate verbose mode, local var.
   integer(I8P), allocatable                    :: local_map_ghost_fluxes_cell(:,:) !< Local map ghost cells update, cells order.
   integer(I8P), allocatable                    :: comm_map_send_ghost_cell(:,:)    !< MPI send map ghost cells update, cells order.
   integer(I8P), allocatable                    :: comm_map_recv_ghost_cell(:,:)    !< MPI send map ghost cells update, cells order.
   integer(I8P), allocatable                    :: local_map_bc_crown(:,:,:)        !< Local map for face BC ghost cells, crown ord.
   integer(I4P)                                 :: c, f, v, n                       !< Counter.
   integer(I4P)                                 :: i, j, k                          !< Counter.
   integer(I4P)                                 :: iii, jjj, kkk                    !< Counter.
   integer(I4P)                                 :: ic, jc, kc                       !< Counter.
   integer(I4P)                                 :: fec                              !< Ghost direction, faces/edges/corners.
   integer(I4P)                                 :: portion                          !< Portion of fec updated (0=>whole fec).
   integer(I4P)                                 :: b_recv                           !< Index of receiving block.
   integer(I4P)                                 :: b_send                           !< Index of sending block.
   integer(I4P)                                 :: imin                             !< Lower limit of i indexes.
   integer(I4P)                                 :: jmin                             !< Lower limit of j indexes.
   integer(I4P)                                 :: kmin                             !< Lower limit of j indexes.
   integer(I4P)                                 :: imax                             !< Upper limit of i indexes.
   integer(I4P)                                 :: jmax                             !< Upper limit of j indexes.
   integer(I4P)                                 :: kmax                             !< Upper limit of k indexes.
   integer(I4P)                                 :: idelta                           !< Delta offset for ghost-inner cells of i.
   integer(I4P)                                 :: jdelta                           !< Delta offset for ghost-inner cells of j.
   integer(I4P)                                 :: kdelta                           !< Delta offset for ghost-inner cells of k.
   integer(I4P)                                 :: recv_ptr, recv_ctr               !< Counter.
   integer(I4P)                                 :: send_ptr, send_ctr               !< Counter.
   integer(I4P), allocatable                    :: c_crown(:)                       !< Counter.

   verbose_ = .false. ; if (present(verbose)) verbose_ = verbose
   if (verbose_) print '(A)', self%mpih%myrankstr//'base_gpu%create_maps_fluxes_cell start'

   if (allocated(self%local_map_ghost_fluxes_cell_gpu)) deallocate(self%local_map_ghost_fluxes_cell_gpu)
   if (allocated(self%field%local_map_ghost)) then
      c = 0
      do f=1, size(self%field%local_map_ghost, dim=1)
         fec     = self%field%local_map_ghost(f, 3 ) ! recv fec
         portion = self%field%local_map_ghost(f, 4 )
         imin    = self%field%local_map_ghost(f, 5 )
         jmin    = self%field%local_map_ghost(f, 6 )
         kmin    = self%field%local_map_ghost(f, 7 )
         imax    = self%field%local_map_ghost(f, 8 )
         jmax    = self%field%local_map_ghost(f, 9 )
         kmax    = self%field%local_map_ghost(f, 10)
         if (fec <=6 .and. portion>0) then
            ! receiving from a block finer than me
            if(fec == 1) then ; imin = 0                  ; imax = 0                  ; idelta =  self%field%grid%ni   ; endif
            if(fec == 2) then ; imin = self%field%grid%ni ; imax = self%field%grid%ni ; idelta = -2*self%field%grid%ni ; endif
            if(fec == 3) then ; jmin = 0                  ; jmax = 0                  ; jdelta =  self%field%grid%nj   ; endif
            if(fec == 4) then ; jmin = self%field%grid%nj ; jmax = self%field%grid%nj ; jdelta = -2*self%field%grid%nj ; endif
            if(fec == 5) then ; kmin = 0                  ; kmax = 0                  ; kdelta =  self%field%grid%nk   ; endif
            if(fec == 6) then ; kmin = self%field%grid%nk ; kmax = self%field%grid%nk ; kdelta = -2*self%field%grid%nk ; endif
            ! receiving a face-like fec from a block finer than me
            do k=kmin, kmax
               do j=jmin, jmax
                  do i=imin, imax
                     c = c + 1
                  enddo
               enddo
            enddo
         endif
      enddo
      if(allocated(self%local_map_ghost_fluxes_cell_gpu)) deallocate(self%local_map_ghost_fluxes_cell_gpu)
      if(c>0) then
         allocate(local_map_ghost_fluxes_cell(1:c,1:9))
         c = 1
         do f=1, size(self%field%local_map_ghost, dim=1)
            b_recv  = self%field%local_map_ghost(f, 1 )
            b_send  = self%field%local_map_ghost(f, 2 )
            fec     = self%field%local_map_ghost(f, 3 ) ! recv fec
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
            if (fec <=6 .and. portion>0) then
               ! receiving from a block finer than me
               if(fec == 1) then ; imin = 0                  ; imax = 0                  ; idelta =  self%field%grid%ni   ; endif
               if(fec == 2) then ; imin = self%field%grid%ni ; imax = self%field%grid%ni ; idelta = -2*self%field%grid%ni ; endif
               if(fec == 3) then ; jmin = 0                  ; jmax = 0                  ; jdelta =  self%field%grid%nj   ; endif
               if(fec == 4) then ; jmin = self%field%grid%nj ; jmax = self%field%grid%nj ; jdelta = -2*self%field%grid%nj ; endif
               if(fec == 5) then ; kmin = 0                  ; kmax = 0                  ; kdelta =  self%field%grid%nk   ; endif
               if(fec == 6) then ; kmin = self%field%grid%nk ; kmax = self%field%grid%nk ; kdelta = -2*self%field%grid%nk ; endif
               do k=kmin, kmax
                  do j=jmin, jmax
                     do i=imin, imax
                        kkk = 2 * k + kdelta
                        jjj = 2 * j + jdelta
                        iii = 2 * i + idelta
                        local_map_ghost_fluxes_cell(c,1:2) = [b_send, b_recv]
                        local_map_ghost_fluxes_cell(c,3:5) = [iii, jjj, kkk] ! send
                        local_map_ghost_fluxes_cell(c,6:8) = [i, j, k]       ! recv
                        local_map_ghost_fluxes_cell(c, 9 ) = fec             ! recv fec
                        c = c + 1
                     enddo
                  enddo
               enddo
            endif
         enddo
         call assign_allocatable_gpu(lhs=self%local_map_ghost_fluxes_cell_gpu,                                                     &
                                     rhs=     local_map_ghost_fluxes_cell,                                                         &
                                     msg=self%mpih%myrankstr//'base_gpu%create_maps_fluxes_cell(local_map_ghost_fluxes_cell_gpu) ',&
                                     verbose=verbose)
         deallocate(local_map_ghost_fluxes_cell)
      endif
   endif

   if (verbose_) print '(A)', self%mpih%myrankstr//'base_gpu%create_maps_fluxes_cell finish'
   endsubroutine create_maps_fluxes_cell

   subroutine initialize(self, field, nv_aux, verbose)
   !< Initialize base backend.
   class(base_gpu_object), intent(inout)        :: self             !< The base backend.
   type(field_object),     intent(in), target   :: field            !< Field variable array.
   integer(I4P),           intent(in), optional :: nv_aux           !< Number of auxiliary variables.
   logical,                intent(in), optional :: verbose          !< Flag to activate verbose mode.
   integer(I4P)                                 :: nv_aux_          !< Number of auxiliary variables (local var).
   integer(I4P), allocatable                    :: fec_1_6_array(:) !< Mapping fec1-26 to fec1-6 for boundaries.

   print '(A)', self%mpih%myrankstr//'base_gpu%initialize start'
   allocate(fec_1_6_array(26))
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

   subroutine initialize_gpu(self, do_mpi_init)
   !< Initialize GPU main data.
   !< @Note This must be the first routine called before.
   class(base_gpu_object), intent(inout)        :: self              !< The base backend.
   logical,                intent(in), optional :: do_mpi_init       !< Flag to activate MPI init call.
   type(cudadeviceprop)                         :: device_properties !< Device properties.

   call self%mpih%initialize(do_mpi_init=do_mpi_init)
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

   subroutine update_ghost_fluxes_local_gpu(self, flx_gpu, fly_gpu, flz_gpu)
   !< Update (local) ghost cells.
   class(base_gpu_object), intent(in)            :: self      !< The base backend.
   real(R8P),              intent(inout), device :: flx_gpu(1:,                    &
                                                            1-self%field%grid%ngc:,&
                                                            1-self%field%grid%ngc:,&
                                                            1-self%field%grid%ngc:,&
                                                            1:) !< Field component to be updated.
   real(R8P),              intent(inout), device :: fly_gpu(1:,                    &
                                                            1-self%field%grid%ngc:,&
                                                            1-self%field%grid%ngc:,&
                                                            1-self%field%grid%ngc:,&
                                                            1:) !< Field component to be updated.
   real(R8P),              intent(inout), device :: flz_gpu(1:,                    &
                                                            1-self%field%grid%ngc:,&
                                                            1-self%field%grid%ngc:,&
                                                            1-self%field%grid%ngc:,&
                                                            1:) !< Field component to be updated.
   call update_ghost_fluxes_local_gpu_cuf(local_map_ghost_fluxes_cell_gpu=self%local_map_ghost_fluxes_cell_gpu, &
                                          ngc=self%field%grid%ngc, flx_gpu=flx_gpu, fly_gpu=fly_gpu, flz_gpu=flz_gpu)
   endsubroutine update_ghost_fluxes_local_gpu

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
                                 comm_map_send_ptr_ghost=self%field%comm_map_send_ptr_ghost,     &
                                 comm_map_recv_ptr_ghost=self%field%comm_map_recv_ptr_ghost,     &
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

   subroutine update_ghost_fluxes_local_gpu_cuf(ngc, local_map_ghost_fluxes_cell_gpu, flx_gpu, fly_gpu, flz_gpu)
   !< Update (local) ghost fluxes cells.
   integer(I4P), intent(in)                         :: ngc                                  !< Ghost cells number.
   integer(I8P), intent(in),    device, allocatable :: local_map_ghost_fluxes_cell_gpu(:,:) !< Local map of ghost cells.
   real(R8P),    intent(inout), device              :: flx_gpu(1:,    &
                                                               1-ngc:,&
                                                               1-ngc:,&
                                                               1-ngc:,1:)                   !< Field component to be updated.
   real(R8P),    intent(inout), device              :: fly_gpu(1:,    &
                                                               1-ngc:,&
                                                               1-ngc:,&
                                                               1-ngc:,1:)                   !< Field component to be updated.
   real(R8P),    intent(inout), device              :: flz_gpu(1:,    &
                                                               1-ngc:,&
                                                               1-ngc:,&
                                                               1-ngc:,1:)                   !< Field component to be updated.
   integer(I4P)                                     :: ic, jc, kc, mf, v                    !< Counter.
   integer(I4P)                                     :: b_recv                               !< Index of receiving block.
   integer(I4P)                                     :: b_send                               !< Index of sending block.
   integer(I4P)                                     :: i_recv                               !< I recv index.
   integer(I4P)                                     :: j_recv                               !< J recv index.
   integer(I4P)                                     :: k_recv                               !< K recv index.
   integer(I4P)                                     :: i_send                               !< I send index.
   integer(I4P)                                     :: j_send                               !< J send index.
   integer(I4P)                                     :: k_send                               !< K send index.
   integer(I4P)                                     :: one_or_eight                         !< Flag triggering 8 cells mean.
   integer(I4P)                                     :: error                                !< Error traping flag.
   integer(I4P)                                     :: iercuda                              !< Error trapping flag for CUDAFortran.
   integer(I4P)                                     :: fec                                  !< Ghost type.

   iercuda = cudaGetLastError()
   if (iercuda /= cudaSuccess) then
      write(stderr, '(A)') 'base_gpu%update_ghost_fluxes_local_gpu_cuf start, cuda ERROR: '//cudaGetErrorString(iercuda)
      call MPI_Abort(MPI_COMM_WORLD, -15, error)
      stop
   endif

   if (.not.allocated(local_map_ghost_fluxes_cell_gpu)) return
   !$cuf kernel do(2) <<<*,*>>>
   do v=1, size(flx_gpu, dim=5)
      do mf=1, size(local_map_ghost_fluxes_cell_gpu, dim=1)
         b_send       = local_map_ghost_fluxes_cell_gpu(mf,1)
         b_recv       = local_map_ghost_fluxes_cell_gpu(mf,2)
         i_send       = local_map_ghost_fluxes_cell_gpu(mf,3)
         j_send       = local_map_ghost_fluxes_cell_gpu(mf,4)
         k_send       = local_map_ghost_fluxes_cell_gpu(mf,5)
         i_recv       = local_map_ghost_fluxes_cell_gpu(mf,6)
         j_recv       = local_map_ghost_fluxes_cell_gpu(mf,7)
         k_recv       = local_map_ghost_fluxes_cell_gpu(mf,8)
         fec          = local_map_ghost_fluxes_cell_gpu(mf,9)
         if (fec==1 .or. fec==2) then
            flx_gpu(b_recv,i_recv,j_recv,k_recv,v) = 0._R8P
            do kc=0,1 ; do jc=0,1 ; do ic=0,0
               flx_gpu(b_recv,i_recv,j_recv,k_recv,v) = flx_gpu(b_recv,i_recv,   j_recv,   k_recv,   v) + &
                                                        flx_gpu(b_send,i_send+ic,j_send+jc,k_send+kc,v)
            enddo ; enddo ; enddo
            flx_gpu(b_recv,i_recv,j_recv,k_recv,v)    = flx_gpu(b_recv,i_recv,j_recv,k_recv,v) * 0.25_R8P
         endif
         if (fec==3 .or. fec==4) then
            fly_gpu(b_recv,i_recv,j_recv,k_recv,v) = 0._R8P
            do kc=0,1 ; do jc=0,0 ; do ic=0,1
               fly_gpu(b_recv,i_recv,j_recv,k_recv,v) = fly_gpu(b_recv,i_recv,   j_recv,   k_recv,   v) + &
                                                        fly_gpu(b_send,i_send+ic,j_send+jc,k_send+kc,v)
            enddo ; enddo ; enddo
            fly_gpu(b_recv,i_recv,j_recv,k_recv,v)    = fly_gpu(b_recv,i_recv,j_recv,k_recv,v) * 0.25_R8P
         endif
         if (fec==5 .or. fec==6) then
            flz_gpu(b_recv,i_recv,j_recv,k_recv,v) = 0._R8P
            do kc=0,0 ; do jc=0,1 ; do ic=0,1
               flz_gpu(b_recv,i_recv,j_recv,k_recv,v) = flz_gpu(b_recv,i_recv,   j_recv,   k_recv,   v) + &
                                                        flz_gpu(b_send,i_send+ic,j_send+jc,k_send+kc,v)
            enddo ; enddo ; enddo
            flz_gpu(b_recv,i_recv,j_recv,k_recv,v)    = flz_gpu(b_recv,i_recv,j_recv,k_recv,v) * 0.25_R8P
         endif
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   iercuda = cudaGetLastError()
   if (iercuda /= cudaSuccess) then
      write(stderr, '(A)') 'base_gpu%update_ghost_fluxes_local_gpu_cuf finish, cuda ERROR: '//cudaGetErrorString(iercuda)
      call MPI_Abort(MPI_COMM_WORLD, -15, error)
      stop
   endif
   endsubroutine update_ghost_fluxes_local_gpu_cuf

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
