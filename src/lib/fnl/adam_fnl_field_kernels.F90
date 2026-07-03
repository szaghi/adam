!< ADAM, field class FNL kernels (FNL backend of [[field_object]]).

#include "fundal.H"

module adam_fnl_field_kernels
!< ADAM, field class FNL kernels (FNL backend of [[field_object]]).

! ADAM objects
use :: adam_seam_interpolation_library, only : SEAM_FILL_COMPATIBLE, seam_meta_unpack, &
                                               seam_interpolate_compatible, seam_interpolate_tricubic
! third party modules
use :: fundal
use :: penf, only : I8P, I4P, R8P

implicit none
private
public :: compute_q_gradient_dev
public :: compute_normL2_residuals_dev
public :: populate_send_buffer_ghost_gpu_dev
public :: receive_recv_buffer_ghost_gpu_dev
public :: update_ghost_local_gpu_dev

contains
   subroutine compute_q_gradient_dev(b, ni, nj, nk, ngc, dx, dy, dz, q_gpu, ivar, gradient)
   !< Compute gradient of q(ivar).
   integer(I4P), intent(in)  :: b                                 !< Block index.
   integer(I4P), intent(in)  :: ni                                !< Grid cells number in I direction.
   integer(I4P), intent(in)  :: nj                                !< Grid cells number in J direction.
   integer(I4P), intent(in)  :: nk                                !< Grid cells number in K direction.
   integer(I4P), intent(in)  :: ngc                               !< Ghost cells number.
   real(R8P),    intent(in)  :: dx                                !< X space step.
   real(R8P),    intent(in)  :: dy                                !< Y space step.
   real(R8P),    intent(in)  :: dz                                !< Z space step.
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field component to which apply gradient.
   integer(I4P), intent(in)  :: ivar                              !< Ghost cells number.
   real(R8P),    intent(out) :: gradient                          !< Maximum gradient of q.
   real(R8P)                 :: grad                              !< Current gradient of q.
   integer(I4P)              :: i, j, k                           !< Counter.
   real(R8P), parameter      :: tol=1.e-12                        !< Gradient denominator tolerance.

   gradient = 0._R8P
   !$acc parallel loop independent DEVICEVAR(q_gpu) reduction(max:gradient)
   !$omp OMPLOOP DEVICEVAR(q_gpu) reduction(max:gradient)
   do k=1, nk
      do j=1, nj
         do i=1, ni
            grad = sqrt(((q_gpu(b,i+1,j,k,ivar) - q_gpu(b,i-1,j,k,ivar))/(2*dx))**2 + &
                        ((q_gpu(b,i,j+1,k,ivar) - q_gpu(b,i,j-1,k,ivar))/(2*dy))**2 + &
                        ((q_gpu(b,i,j,k+1,ivar) - q_gpu(b,i,j,k-1,ivar))/(2*dz))**2)
            grad = grad/(abs(q_gpu(b,i,j,k,ivar))+tol)
            gradient = max(gradient, grad)
         enddo
      enddo
   enddo
   endsubroutine compute_q_gradient_dev

   subroutine compute_normL2_residuals_dev(ni, nj, nk, ngc, nv, blocks_number, dq_gpu, norm)
   !< Compute L2 norm of residuals.
   integer(I4P), intent(in)    :: ni                                 !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                 !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                 !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                !< Ghost grid number.
   integer(I4P), intent(in)    :: nv                                 !< Number of states variables.
   integer(I4P), intent(in)    :: blocks_number                      !< Number of blocks.
   real(R8P),    intent(in)    :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   real(R8P),    intent(inout) :: norm(1:)                           !< Residuals norm.
   real(R8P)                   :: norm_gpu                           !< Residuals norm, local scalar buffer for reduction.
   integer(I4P)                :: b, i, j, k, v                      !< Counter.

   do v=1, nv
      norm_gpu = 0._R8P
      !$acc parallel loop independent DEVICEVAR(dq_gpu) reduction(+:norm_gpu)
      !$omp OMPLOOP DEVICEVAR(dq_gpu) reduction(+:norm_gpu)
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  norm_gpu = norm_gpu + dq_gpu(b, i, j, k, v)**2
               enddo
            enddo
         enddo
      enddo
      norm(v) = norm_gpu
   enddo
   endsubroutine compute_normL2_residuals_dev

   subroutine populate_send_buffer_ghost_gpu_dev(ngc, seam_ghost_fill, comm_map_send_ghost_cell_gpu, &
                                                 send_buffer_ghost_gpu, q_gpu)
   !< Polulate send buffer ghost GPU.
   integer(I4P),          intent(in)    :: ngc                                    !< Ghost cells number.
   integer(I4P),          intent(in)    :: seam_ghost_fill                        !< Seam ghost-fill regime (flag-4 rows).
   integer(I8P), pointer, intent(in)    :: comm_map_send_ghost_cell_gpu(:,:)      !< Comm map, cell information.
   real(R8P),    pointer, intent(inout) :: send_buffer_ghost_gpu(:)               !< Send buffer of ghost cells.
   real(R8P),             intent(inout) :: q_gpu(1:,    &
                                                 1-ngc:,&
                                                 1-ngc:,&
                                                 1-ngc:,&
                                                 1:)                              !< Field component to be updated.
   integer(I4P)                         :: ic, jc, kc, sf                         !< Counter.
   integer(I4P)                         :: b_send, i_send, j_send, k_send, v_send !< Send indexes.
   integer(I4P)                         :: c_recv                                 !< Counter.
   integer(I4P)                         :: one_or_eight                           !< Flag triggering 8 cells mean.

   if (associated(comm_map_send_ghost_cell_gpu)) then
      !$acc parallel loop independent DEVICEVAR(comm_map_send_ghost_cell_gpu, send_buffer_ghost_gpu, q_gpu)
      !$omp OMPLOOP DEVICEVAR(comm_map_send_ghost_cell_gpu, send_buffer_ghost_gpu, q_gpu)
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
         elseif (one_or_eight==4) then
            ! coarse->fine seam interpolation, pack-side (issue #21 N2, ported #22 F3); metadata in column 8
            send_buffer_ghost_gpu(c_recv) = interp_seam_ghost_gpu_dev(regime=seam_ghost_fill,                          &
                                                                      meta=int(comm_map_send_ghost_cell_gpu(sf,8), I4P), &
                                                                      ngc=ngc, q_gpu=q_gpu, v=v_send, b_send=b_send,   &
                                                                      i_send=i_send, j_send=j_send, k_send=k_send)
         else
            send_buffer_ghost_gpu(c_recv) = 0._R8P
            do kc=0,1 ; do jc=0,1 ; do ic=0,1
               send_buffer_ghost_gpu(c_recv) = send_buffer_ghost_gpu(c_recv) + &
                                               q_gpu(b_send,i_send+ic,j_send+jc,k_send+kc,v_send)
            enddo ; enddo ; enddo
            send_buffer_ghost_gpu(c_recv) = send_buffer_ghost_gpu(c_recv) * 0.125_R8P
         endif
      enddo
   endif
   endsubroutine populate_send_buffer_ghost_gpu_dev

   subroutine receive_recv_buffer_ghost_gpu_dev(ngc, comm_map_recv_ghost_cell_gpu, recv_buffer_ghost_gpu, q_gpu)
   !< Receive recv buffer ghost GPU.
   integer(I4P),          intent(in)    :: ngc                                    !< Ghost cells number.
   integer(I8P), pointer, intent(in)    :: comm_map_recv_ghost_cell_gpu(:,:)      !< Comm map, cell information.
   real(R8P),    pointer, intent(inout) :: recv_buffer_ghost_gpu(:)               !< Receive buffer of ghost cells.
   real(R8P),             intent(inout) :: q_gpu(1:,    &
                                                 1-ngc:,&
                                                 1-ngc:,&
                                                 1-ngc:,&
                                                 1:)                              !< Field component to be updated.
   integer(I4P)                         :: rf                                     !< Counter.
   integer(I4P)                         :: b_recv, i_recv, j_recv, k_recv, v_recv !< Receive indexes.
   integer(I4P)                         :: c_send                                 !< Counter.

   if (associated(comm_map_recv_ghost_cell_gpu)) then
      !$acc parallel loop independent DEVICEVAR(comm_map_recv_ghost_cell_gpu, recv_buffer_ghost_gpu, q_gpu)
      !$omp OMPLOOP DEVICEVAR(comm_map_recv_ghost_cell_gpu, recv_buffer_ghost_gpu, q_gpu)
      do rf=1, size(comm_map_recv_ghost_cell_gpu, dim=1)
         c_send = comm_map_recv_ghost_cell_gpu(rf,1)
         b_recv = comm_map_recv_ghost_cell_gpu(rf,2)
         i_recv = comm_map_recv_ghost_cell_gpu(rf,3)
         j_recv = comm_map_recv_ghost_cell_gpu(rf,4)
         k_recv = comm_map_recv_ghost_cell_gpu(rf,5)
         v_recv = comm_map_recv_ghost_cell_gpu(rf,6)
         q_gpu(b_recv,i_recv,j_recv,k_recv,v_recv) = recv_buffer_ghost_gpu(c_send)
      enddo
   endif
   endsubroutine receive_recv_buffer_ghost_gpu_dev

   subroutine update_ghost_local_gpu_dev(ngc, seam_ghost_fill, l_map_ghost_cell_gpu, q_gpu)
   !< Update (local) ghost cells.
   integer(I4P), intent(in)          :: ngc                       !< Ghost cells number.
   integer(I4P), intent(in)          :: seam_ghost_fill           !< Seam ghost-fill regime (flag-4 rows).
   integer(I8P), intent(in), pointer :: l_map_ghost_cell_gpu(:,:) !< Local map of ghost cells.
   real(R8P),    intent(inout)       :: q_gpu(1:,    &
                                              1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:,1:)          !< Field component to be updated.
   integer(I4P)                      :: ic, jc, kc, mf, v         !< Counter.
   integer(I4P)                      :: b_recv                    !< Index of receiving block.
   integer(I4P)                      :: b_send                    !< Index of sending block.
   integer(I4P)                      :: i_recv                    !< I recv index.
   integer(I4P)                      :: j_recv                    !< J recv index.
   integer(I4P)                      :: k_recv                    !< K recv index.
   integer(I4P)                      :: i_send                    !< I send index.
   integer(I4P)                      :: j_send                    !< J send index.
   integer(I4P)                      :: k_send                    !< K send index.
   integer(I4P)                      :: one_or_eight              !< Flag triggering 8 cells mean.

   if (.not.associated(l_map_ghost_cell_gpu)) return
   !$acc parallel loop independent DEVICEVAR(l_map_ghost_cell_gpu, q_gpu)
   !$omp OMPLOOP DEVICEVAR(l_map_ghost_cell_gpu, q_gpu)
   do v=1, size(q_gpu, dim=5)
      do mf=1, size(l_map_ghost_cell_gpu, dim=1)
         b_send       = l_map_ghost_cell_gpu(mf,1)
         b_recv       = l_map_ghost_cell_gpu(mf,2)
         i_send       = l_map_ghost_cell_gpu(mf,3)
         j_send       = l_map_ghost_cell_gpu(mf,4)
         k_send       = l_map_ghost_cell_gpu(mf,5)
         i_recv       = l_map_ghost_cell_gpu(mf,6)
         j_recv       = l_map_ghost_cell_gpu(mf,7)
         k_recv       = l_map_ghost_cell_gpu(mf,8)
         one_or_eight = l_map_ghost_cell_gpu(mf,9)
         if (one_or_eight==1) then
            q_gpu(b_recv,i_recv, j_recv, k_recv,v) = q_gpu(b_send, i_send,j_send,k_send,v)
         elseif (one_or_eight==4) then
            ! coarse->fine seam interpolation (issue #21 N2, ported #22 F3); metadata in column 10
            q_gpu(b_recv,i_recv,j_recv,k_recv,v) = interp_seam_ghost_gpu_dev(regime=seam_ghost_fill,                &
                                                                             meta=int(l_map_ghost_cell_gpu(mf,10), I4P), &
                                                                             ngc=ngc, q_gpu=q_gpu, v=v,             &
                                                                             b_send=b_send, i_send=i_send,          &
                                                                             j_send=j_send, k_send=k_send)
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
   endsubroutine update_ghost_local_gpu_dev

   ! private procedures
   pure function interp_seam_ghost_gpu_dev(regime, meta, ngc, q_gpu, v, b_send, i_send, j_send, k_send) result(value_)
   !< Device twin of `field_object%interp_seam_ghost` (issue #22 F3): evaluate
   !< the coarse->fine seam ghost interpolant for one (variable, fine ghost)
   !< pair on the FNL layout `q_gpu(b,i,j,k,v)`, calling the single-source
   !< evaluators of `adam_seam_interpolation_library` (device-pinned by
   !< `test_seam_interpolation_fnl`, #22 F2).
   !<
   !< Race discipline (CLAUDE-gpu, #22 F1-bis/F2): the footprint gather
   !< buffers are CONSTANT-BOUND locals of this `acc routine seq` — per-thread
   !< stack storage, the WENO device-kernels precedent — never `private()`
   !< clauses of a contained kernel nor array-section actuals (both are
   !< nvfortran race shapes). The anchor and the shift-clamped positions in
   !< `meta` guarantee the footprint reads REAL donor cells only (G3), so the
   !< gather never sees ghost cells written concurrently by other map rows.
   integer(I4P), intent(in) :: regime           !< Active fill regime (SEAM_FILL_COMPATIBLE | SEAM_FILL_TRICUBIC).
   integer(I4P), intent(in) :: meta             !< Packed interpolation metadata (see seam_meta_pack).
   integer(I4P), intent(in) :: ngc              !< Ghost cells number (bounds of q_gpu).
   real(R8P),    intent(in) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field variables, FNL layout (b,i,j,k,v).
   integer(I4P), intent(in) :: v                !< Variable index.
   integer(I4P), intent(in) :: b_send           !< Donor block index.
   integer(I4P), intent(in) :: i_send           !< Anchor donor cell, i.
   integer(I4P), intent(in) :: j_send           !< Anchor donor cell, j.
   integer(I4P), intent(in) :: k_send           !< Anchor donor cell, k.
   real(R8P)                :: value_           !< Interpolated fine-ghost value.
   real(R8P)                :: fp4(1:4,1:4,1:4) !< Tricubic donor footprint.
   real(R8P)                :: fp3(1:3,1:3,1:3) !< Compatible donor footprint.
   integer(I4P)             :: sub(1:3)         !< Octant sub-position per direction.
   integer(I4P)             :: p4(1:3)          !< Tricubic anchor position per direction.
   integer(I4P)             :: p3(1:3)          !< Compatible anchor position per direction.
   integer(I4P)             :: ii, jj, kk       !< Footprint counters.
   !$acc routine seq

   call seam_meta_unpack(meta=meta, sub=sub, p4=p4, p3=p3)
   if (regime == SEAM_FILL_COMPATIBLE) then
      do kk=1, 3
         do jj=1, 3
            do ii=1, 3
               fp3(ii,jj,kk) = q_gpu(b_send, i_send+ii-p3(1), j_send+jj-p3(2), k_send+kk-p3(3), v)
            enddo
         enddo
      enddo
      value_ = seam_interpolate_compatible(footprint=fp3, sub=sub, anchor_pos=p3)
   else
      do kk=1, 4
         do jj=1, 4
            do ii=1, 4
               fp4(ii,jj,kk) = q_gpu(b_send, i_send+ii-p4(1), j_send+jj-p4(2), k_send+kk-p4(3), v)
            enddo
         enddo
      enddo
      value_ = seam_interpolate_tricubic(footprint=fp4, sub=sub, anchor_pos=p4)
   endif
   endfunction interp_seam_ghost_gpu_dev
endmodule adam_fnl_field_kernels
