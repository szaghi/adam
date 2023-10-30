!< ADAM, field class NVF kernels (NVF backend of [[field_object]]).
module adam_field_nvf_kernels
!< ADAM, field class NVF kernels (NVF backend of [[field_object]]).

use penf, only : I8P, I4P, R8P
use CUDAFOR

implicit none
private
public :: compute_normL2_residuals_cuf
public :: copy_transpose_gpu_cpu_cuf
public :: populate_send_buffer_ghost_gpu_cuf
public :: receive_recv_buffer_ghost_gpu_cuf
public :: update_ghost_local_gpu_cuf

contains
   subroutine compute_normL2_residuals_cuf(ni, nj, nk, ngc, nv, blocks_number, dq_gpu, norm)
   !< Compute L2 norm of residuals.
   integer(I4P), intent(in)         :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)         :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)         :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)         :: ngc                                !< Ghost grid number.
   integer(I4P), intent(in)         :: nv                                 !< Number of states variables.
   integer(I4P), intent(in)         :: blocks_number                      !< Number of blocks.
   real(R8P),    intent(in), device :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   real(R8P),    intent(inout)      :: norm(1:)                           !< Residuals norm.
   real(R8P)                        :: norm_gpu                           !< Residuals norm, local scalar buffer for reduction.
   integer(I4P)                     :: b, i, j, k, v                      !< Counter.
   integer(I4P)                     :: iercuda                            !< Error trapping flag for CUDAFortran.

   do v=1, nv
      norm_gpu = 0._R8P
      !$cuf kernel do(4) <<<*,*>>> reduce(+:norm_gpu)
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  norm_gpu = norm_gpu + dq_gpu(b, i, j, k, v)**2
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
      norm(v) = norm_gpu
   enddo
   endsubroutine compute_normL2_residuals_cuf

   subroutine copy_transpose_gpu_cpu_cuf(ni, nj, nk, ngc, nv, blocks_number, q_gpu, q_t_gpu)
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
   integer(I4P)                        :: i, j, k, b, v !< Counter.
   integer(I4P)                        :: iercuda       !< Error trapping flag for CUDAFortran.

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
   endsubroutine copy_transpose_gpu_cpu_cuf

   subroutine populate_send_buffer_ghost_gpu_cuf(ngc, comm_map_send_ghost_cell_gpu, send_buffer_ghost_gpu, q_gpu)
   !< Polulate send buffer ghost GPU.
   integer(I4P), intent(in)                         :: ngc                                    !< Ghost cells number.
   integer(I8P), allocatable, intent(in),    device :: comm_map_send_ghost_cell_gpu(:,:)      !< Comm map, cell information.
   real(R8P),    allocatable, intent(inout), device :: send_buffer_ghost_gpu(:)               !< Send buffer of ghost cells.
   real(R8P),                 intent(inout), device :: q_gpu(1:,    &
                                                             1-ngc:,&
                                                             1-ngc:,&
                                                             1-ngc:,&
                                                             1:)                              !< Field component to be updated.
   integer(I4P)                                     :: ic, jc, kc, sf                         !< Counter.
   integer(I4P)                                     :: b_send, i_send, j_send, k_send, v_send !< Send indexes.
   integer(I4P)                                     :: c_recv                                 !< Counter.
   integer(I4P)                                     :: one_or_eight                           !< Flag triggering 8 cells mean.
   integer(I4P)                                     :: iercuda                                !< Error trapping flag for CUDAFor.

   if (allocated(comm_map_send_ghost_cell_gpu)) then
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
   endsubroutine populate_send_buffer_ghost_gpu_cuf

   subroutine receive_recv_buffer_ghost_gpu_cuf(ngc, comm_map_recv_ghost_cell_gpu, recv_buffer_ghost_gpu, q_gpu)
   !< Receive recv buffer ghost GPU.
   integer(I4P), intent(in)                         :: ngc                                    !< Ghost cells number.
   integer(I8P), allocatable, intent(in),    device :: comm_map_recv_ghost_cell_gpu(:,:)      !< Comm map, cell information.
   real(R8P),    allocatable, intent(inout), device :: recv_buffer_ghost_gpu(:)               !< Receive buffer of ghost cells.
   real(R8P),                 intent(inout), device :: q_gpu(1:,    &
                                                             1-ngc:,&
                                                             1-ngc:,&
                                                             1-ngc:,&
                                                             1:)                              !< Field component to be updated.
   integer(I4P)                                     :: rf                                     !< Counter.
   integer(I4P)                                     :: b_recv, i_recv, j_recv, k_recv, v_recv !< Receive indexes.
   integer(I4P)                                     :: c_send                                 !< Counter.
   integer(I4P)                                     :: iercuda                                !< Error trapping flag for CUDAFor.

   if (allocated(comm_map_recv_ghost_cell_gpu)) then
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
   endsubroutine receive_recv_buffer_ghost_gpu_cuf

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
endmodule adam_field_nvf_kernels
