!< ADAM, field class FNL kernels (FNL backend of [[field_object]]).

#include "fundal.H"

module adam_fnl_field_kernels
!< ADAM, field class FNL kernels (FNL backend of [[field_object]]).

! ADAM objects
use :: adam_seam_interpolation_library, only : SEAM_FILL_COMPATIBLE, SEAM_FILL_INJECTION, seam_meta_unpack, SEAM_DENSITY_GUARD, &
                                               SEAM_FILL_QUINTIC, SEAM_SUB_COINCIDENT, seam_quintic_weight, &
                                               SEAM_FILL_WENO6, seam_fill_weno6_1d, seam_mp_bound, &
                                               seam_bracket_nodes, seam_bracket_range, &
                                               SEAM_FILL_QUINTIC_MP, seam_fill_quintic_mp_1d, &
                                               seam_interpolate_compatible, seam_interpolate_tricubic, &
                                               SEAM_RESTRICT_MEAN, SEAM_RESTRICT_WENO6, SEAM_RESTRICT_IDENTITY, &
                                               seam_restriction_meta_unpack, seam_restriction_weight, &
                                               seam_restriction_weno6_1d
! third party modules
use :: fundal
use :: penf, only : I8P, I4P, R8P

implicit none
private
public :: compute_q_gradient_dev
public :: compute_normL2_residuals_dev
public :: pack_blocks_transposed_dev
public :: pack_seam_rows_dev
public :: unpack_blocks_transposed_dev
public :: unpack_seam_rows_dev
public :: populate_send_buffer_ghost_gpu_dev
public :: receive_recv_buffer_ghost_gpu_dev
public :: seam_ghosts_to_stages_dev
public :: update_ghost_local_gpu_dev
public :: reduce_ghost_local_gpu_dev
public :: pack_reduce_ghost_mpi_gpu_dev
public :: fold_reduce_ghost_mpi_gpu_dev

contains
   subroutine pack_blocks_transposed_dev(blocks_number, ni, nj, nk, ngc, nv, a_gpu, buf_gpu)
   !< Pack blocks `1:blocks_number` of the device field `a_gpu(b, i, j, k, v)`, ghost cells included, into `buf_gpu` in
   !< the host layout `(v, i, j, k, b)` (issue #75): `buf_gpu` is then one contiguous copy of the host slab
   !< `a(:, :, :, :, 1:blocks_number)`, so a host round trip moves the blocks in use, not the block capacity `nb` (the
   !< leading device index is the block, so the blocks in use are a strided section of `a_gpu`).
   integer(I4P), intent(in)    :: blocks_number                       !< Blocks to pack.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                     !< Block dimensions, ghost cells number.
   integer(I4P), intent(in)    :: nv                                  !< Variables number.
   real(R8P),    intent(in)    :: a_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Field (b, i, j, k, v).
   real(R8P),    intent(inout) :: buf_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Packed blocks (v, i, j, k, b).
   integer(I4P)                :: b, i, j, k, v                       !< Counters.

   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(a_gpu,buf_gpu) &
   !$acc& firstprivate(blocks_number,ni,nj,nk,ngc,nv)
   !$omp OMPLOOP collapse(5) DEVICEPTR(a_gpu,buf_gpu) firstprivate(blocks_number,ni,nj,nk,ngc,nv)
   do b=1, blocks_number
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do v=1, nv
      buf_gpu(v,i,j,k,b) = a_gpu(b,i,j,k,v)
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine pack_blocks_transposed_dev

   subroutine unpack_blocks_transposed_dev(blocks_number, ni, nj, nk, ngc, nv, buf_gpu, a_gpu)
   !< Unpack `buf_gpu`, blocks `1:blocks_number` in the host layout `(v, i, j, k, b)`, into the device field
   !< `a_gpu(b, i, j, k, v)`, ghost cells included: the inverse of `pack_blocks_transposed_dev` (issue #75); the blocks
   !< beyond `blocks_number` are untouched.
   integer(I4P), intent(in)    :: blocks_number                       !< Blocks to unpack.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                     !< Block dimensions, ghost cells number.
   integer(I4P), intent(in)    :: nv                                  !< Variables number.
   real(R8P),    intent(in)    :: buf_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Packed blocks (v, i, j, k, b).
   real(R8P),    intent(inout) :: a_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Field (b, i, j, k, v).
   integer(I4P)                :: b, i, j, k, v                       !< Counters.

   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(buf_gpu,a_gpu) &
   !$acc& firstprivate(blocks_number,ni,nj,nk,ngc,nv)
   !$omp OMPLOOP collapse(5) DEVICEPTR(buf_gpu,a_gpu) firstprivate(blocks_number,ni,nj,nk,ngc,nv)
   do v=1, nv
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, blocks_number
      a_gpu(b,i,j,k,v) = buf_gpu(v,i,j,k,b)
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine unpack_blocks_transposed_dev

   subroutine pack_seam_rows_dev(row_start, row_count, nv, ngc, regime, rows_gpu, q_gpu, buf_gpu, density_var)
   !< Pack the ghost values of the seam send rows `[rank, b, i, j, k, kind, meta]` `row_start..row_start+row_count-1`
   !< from `q_gpu` into `buf_gpu`, `nv` values per row in row order (issues #40, #52), the device twin of
   !< `adam_seam_exchange%pack_seam_rows`: kind 1 copies the cell, kind 4 interpolates the fine ghost around the coarse
   !< anchor cell (`interp_seam_ghost_gpu_dev`, the intra-realm flag-4 fill; injection copies the anchor), kind 8 is the
   !< mean of the 2x2x2 fine cells from the base cell (the intra-realm flag-8 order). The kinds are the values of
   !< `adam_seam_exchange` SEAM_ROW_*.
   integer(I4P), intent(in)    :: row_start                         !< First row.
   integer(I4P), intent(in)    :: row_count                         !< Number of rows.
   integer(I4P), intent(in)    :: nv                                !< Variables number.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)    :: regime                            !< Coarse->fine fill regime (`maps%seam_ghost_fill`).
   integer(I4P), intent(in)    :: rows_gpu(1:,1:)                   !< Rows [rank, b, i, j, k, kind, meta].
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field (b, i, j, k, v).
   real(R8P),    intent(inout) :: buf_gpu(1:)                       !< Packed values.
   integer(I4P), intent(in)    :: density_var                       !< Seam fill density guard (`maps%seam_density_var`).
   integer(I4P)                :: c, v, row                         !< Counters.
   integer(I4P)                :: b, i, j, k                        !< Cell.
   integer(I4P)                :: ic, jc, kc                        !< 2x2x2 counters.
   real(R8P)                   :: total                             !< Restriction sum.

   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(rows_gpu,q_gpu,buf_gpu) &
   !$acc& firstprivate(row_start,row_count,nv,ngc,regime,density_var) private(row,b,i,j,k,total)
   !$omp OMPLOOP collapse(2) DEVICEPTR(rows_gpu,q_gpu,buf_gpu)                                                &
   !$omp& firstprivate(row_start,row_count,nv,ngc,regime,density_var)                                           &
   !$omp& private(row,b,i,j,k,total)
   do c=1, row_count
   do v=1, nv
      row = row_start + c - 1
      b = rows_gpu(row,2) ; i = rows_gpu(row,3) ; j = rows_gpu(row,4) ; k = rows_gpu(row,5)
      if (rows_gpu(row,6) == 4 .and. regime /= SEAM_FILL_INJECTION) then
         buf_gpu(nv*(c-1)+v) = interp_seam_ghost_gpu_dev(regime=regime, meta=rows_gpu(row,7), ngc=ngc, q_gpu=q_gpu, &
                                                         v=v, b_send=b, i_send=i, j_send=j, k_send=k,           &
                                                         density_var=density_var)
      elseif (rows_gpu(row,6) == 8 .and. rows_gpu(row,7) /= 0) then
         ! high-order restriction (#68 C2, C3): the fine owner packed its codes and regime in the meta; 0 is the mean
         buf_gpu(nv*(c-1)+v) = restrict_seam_high_order_gpu_dev(meta=rows_gpu(row,7), ngc=ngc, q_gpu=q_gpu, v=v, &
                                                                b_send=b, i_send=i, j_send=j, k_send=k)
      elseif (rows_gpu(row,6) == 8) then
         total = 0._R8P
         !$acc loop seq
         do kc=0,1 ; do jc=0,1 ; do ic=0,1
            total = total + q_gpu(b,i+ic,j+jc,k+kc,v)
         enddo ; enddo ; enddo
         buf_gpu(nv*(c-1)+v) = total * 0.125_R8P
      else
         buf_gpu(nv*(c-1)+v) = q_gpu(b,i,j,k,v)
      endif
   enddo
   enddo
   endsubroutine pack_seam_rows_dev

   subroutine unpack_seam_rows_dev(row_start, row_count, nv, ngc, rows_gpu, buf_gpu, q_gpu)
   !< Scatter `buf_gpu` (`nv` values per row in row order) into the cells of the cross-rank seam rows
   !< `[rank, b, i, j, k]` `row_start..row_start+row_count-1` of `q_gpu` (issue #40).
   integer(I4P), intent(in)    :: row_start                         !< First row.
   integer(I4P), intent(in)    :: row_count                         !< Number of rows.
   integer(I4P), intent(in)    :: nv                                !< Variables number.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)    :: rows_gpu(1:,1:)                   !< Rows [rank, b, i, j, k].
   real(R8P),    intent(in)    :: buf_gpu(1:)                       !< Packed values.
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field (b, i, j, k, v).
   integer(I4P)                :: c, v, row                         !< Counters.

   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(rows_gpu,q_gpu,buf_gpu) &
   !$acc& firstprivate(row_start,row_count,nv) private(row)
   !$omp OMPLOOP collapse(2) DEVICEPTR(rows_gpu,q_gpu,buf_gpu) firstprivate(row_start,row_count,nv) private(row)
   do c=1, row_count
   do v=1, nv
      row = row_start + c - 1
      q_gpu(rows_gpu(row,2),rows_gpu(row,3),rows_gpu(row,4),rows_gpu(row,5),v) = buf_gpu(nv*(c-1)+v)
   enddo
   enddo
   endsubroutine unpack_seam_rows_dev

   subroutine seam_ghosts_to_stages_dev(row_start, row_count, nv, ngc, ns, cols, rows_gpu, q_gpu, q_rk_gpu)
   !< Copy the seam ghost cells of the rows `row_start..row_start+row_count-1` from `q_gpu` into every stage of
   !< `q_rk_gpu` (issue #113, the device twin of `adam_seam_exchange%seam_ghosts_to_stages`). `cols` are the columns of
   !< the rows holding the ghost cell (b, i, j, k): `[3, 7, 8, 9]` for the local rows, `[2, 3, 4, 5]` for the receive rows.
   integer(I4P), intent(in)    :: row_start                                  !< First row.
   integer(I4P), intent(in)    :: row_count                                  !< Number of rows.
   integer(I4P), intent(in)    :: nv                                         !< Variables number.
   integer(I4P), intent(in)    :: ngc                                        !< Ghost cells number.
   integer(I4P), intent(in)    :: ns                                         !< Stages number.
   integer(I4P), intent(in)    :: cols(4)                                    !< Columns of the ghost cell (b, i, j, k).
   integer(I4P), intent(in)    :: rows_gpu(1:,1:)                            !< Seam rows.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)          !< Committed state (b, i, j, k, v).
   real(R8P),    intent(inout) :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:)    !< Stage buffers (b, i, j, k, v, s).
   integer(I4P)                :: c, v, s, row, cb, ci, cj, ck               !< Counters, row, ghost cell columns.

   cb = cols(1) ; ci = cols(2) ; cj = cols(3) ; ck = cols(4)
   !$acc parallel loop independent gang vector collapse(3) DEVICEVAR(rows_gpu,q_gpu,q_rk_gpu) &
   !$acc& firstprivate(row_start,row_count,nv,ns,cb,ci,cj,ck) private(row)
   !$omp OMPLOOP collapse(3) DEVICEPTR(rows_gpu,q_gpu,q_rk_gpu) firstprivate(row_start,row_count,nv,ns,cb,ci,cj,ck) &
   !$omp& private(row)
   do s=1, ns
   do c=1, row_count
   do v=1, nv
      row = row_start + c - 1
      q_rk_gpu(rows_gpu(row,cb),rows_gpu(row,ci),rows_gpu(row,cj),rows_gpu(row,ck),v,s) = &
         q_gpu(rows_gpu(row,cb),rows_gpu(row,ci),rows_gpu(row,cj),rows_gpu(row,ck),v)
   enddo
   enddo
   enddo
   endsubroutine seam_ghosts_to_stages_dev

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
   !$omp OMPLOOP collapse(3) DEVICEPTR(q_gpu) reduction(max:gradient)
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
      !$acc parallel loop independent collapse(4) DEVICEVAR(dq_gpu) reduction(+:norm_gpu)
      !$omp OMPLOOP collapse(4) DEVICEPTR(dq_gpu) reduction(+:norm_gpu)
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

   subroutine populate_send_buffer_ghost_gpu_dev(ngc, seam_ghost_fill, seam_restriction, ri, rj, rk,              &
                                                 comm_map_send_ghost_cell_gpu, send_buffer_ghost_gpu, q_gpu,     &
                                                 seam_density_var)
   !< Polulate send buffer ghost GPU.
   !<
   !< A fine->coarse row averages the `ri*rj*rk` fine cells of its coarse ghost: 2 per refined axis, 1 on an axis the
   !< tree never refines (quadtree z, issue #46).
   integer(I4P),          intent(in)    :: ngc                                    !< Ghost cells number.
   integer(I4P),          intent(in)    :: seam_ghost_fill                        !< Seam ghost-fill regime (flag-4 rows).
   integer(I4P),          intent(in)    :: seam_restriction                       !< Seam restriction regime (flag-8 rows).
   integer(I4P),          intent(in)    :: ri, rj, rk                             !< Per-axis refinement factor.
   integer(I8P), pointer, intent(in)    :: comm_map_send_ghost_cell_gpu(:,:)      !< Comm map, cell information.
   real(R8P),    pointer, intent(inout) :: send_buffer_ghost_gpu(:)               !< Send buffer of ghost cells.
   integer(I4P),          intent(in)    :: seam_density_var                       !< Seam fill density guard (#101).
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
      !$omp OMPLOOP DEVICEPTR(comm_map_send_ghost_cell_gpu, send_buffer_ghost_gpu, q_gpu)
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
                                                                      i_send=i_send, j_send=j_send, k_send=k_send,     &
                                                                      density_var=seam_density_var)
         elseif (seam_restriction /= SEAM_RESTRICT_MEAN .and. comm_map_send_ghost_cell_gpu(sf,8) /= 0) then
            ! fine->coarse restriction, 6-point Lagrange or WENO, pack-side (#68 C2, C3); meta in column 8
            send_buffer_ghost_gpu(c_recv) = restrict_seam_high_order_gpu_dev(                                         &
                                               meta=int(comm_map_send_ghost_cell_gpu(sf,8), I4P), ngc=ngc, q_gpu=q_gpu, &
                                               v=v_send, b_send=b_send, i_send=i_send, j_send=j_send, k_send=k_send)
         else
            send_buffer_ghost_gpu(c_recv) = 0._R8P
            do kc=0,rk-1 ; do jc=0,rj-1 ; do ic=0,ri-1
               send_buffer_ghost_gpu(c_recv) = send_buffer_ghost_gpu(c_recv) + &
                                               q_gpu(b_send,i_send+ic,j_send+jc,k_send+kc,v_send)
            enddo ; enddo ; enddo
            send_buffer_ghost_gpu(c_recv) = send_buffer_ghost_gpu(c_recv) / real(ri*rj*rk, R8P)
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
      !$omp OMPLOOP DEVICEPTR(comm_map_recv_ghost_cell_gpu, recv_buffer_ghost_gpu, q_gpu)
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

   subroutine update_ghost_local_gpu_dev(ngc, seam_ghost_fill, seam_restriction, ri, rj, rk, l_map_ghost_cell_gpu, q_gpu, &
                                         seam_density_var)
   !< Update (local) ghost cells.
   !<
   !< A fine->coarse row averages the `ri*rj*rk` fine cells of its coarse ghost: 2 per refined axis, 1 on an axis the
   !< tree never refines (quadtree z, issue #46).
   integer(I4P), intent(in)          :: ngc                       !< Ghost cells number.
   integer(I4P), intent(in)          :: seam_ghost_fill           !< Seam ghost-fill regime (flag-4 rows).
   integer(I4P), intent(in)          :: seam_restriction          !< Seam restriction regime (flag-8 rows).
   integer(I4P), intent(in)          :: ri, rj, rk                !< Per-axis refinement factor.
   integer(I8P), intent(in), pointer :: l_map_ghost_cell_gpu(:,:) !< Local map of ghost cells.
   integer(I4P), intent(in)          :: seam_density_var          !< Seam fill density guard (#101).
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
   ! collapse(2) (#68 C3b): without it nvfortran parallelised the variables only (9 threads) and ran the rows in sequence
   ! ("loop carried dependence of q_gpu": the writes go through the map), which made this kernel 98 % of the FNL GPU
   ! time of the PV-5 blast. The rows are independent: one row per receiver cell (dedupe_last_writer, #22 F1), and every
   ! row reads real donor cells only, never a ghost another row writes (the #22 race discipline).
   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(l_map_ghost_cell_gpu, q_gpu)
   !$omp OMPLOOP collapse(2) DEVICEPTR(l_map_ghost_cell_gpu, q_gpu)
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
                                                                             j_send=j_send, k_send=k_send,          &
                                                                             density_var=seam_density_var)
         elseif (seam_restriction /= SEAM_RESTRICT_MEAN .and. l_map_ghost_cell_gpu(mf,10) /= 0) then
            ! fine->coarse restriction, 6-point Lagrange or WENO (#68 C2, C3); meta in column 10
            q_gpu(b_recv,i_recv,j_recv,k_recv,v) = restrict_seam_high_order_gpu_dev(                                 &
                                                      meta=int(l_map_ghost_cell_gpu(mf,10), I4P), ngc=ngc, q_gpu=q_gpu, &
                                                      v=v, b_send=b_send, i_send=i_send, j_send=j_send, k_send=k_send)
         else
            q_gpu(b_recv,i_recv,j_recv,k_recv,v) = 0._R8P
            do kc=0,rk-1 ; do jc=0,rj-1 ; do ic=0,ri-1
               q_gpu(b_recv,i_recv,j_recv,k_recv,v) = q_gpu(b_recv,i_recv,   j_recv,   k_recv,   v) + &
                                                      q_gpu(b_send,i_send+ic,j_send+jc,k_send+kc,v)
            enddo ; enddo ; enddo
            q_gpu(b_recv,i_recv,j_recv,k_recv,v) = q_gpu(b_recv,i_recv,j_recv,k_recv,v) / real(ri*rj*rk, R8P)
         endif
      enddo
   enddo
   endsubroutine update_ghost_local_gpu_dev

   subroutine reduce_ghost_local_gpu_dev(ngc, v_first, v_last, l_map_ghost_cell_gpu, q_gpu, bad_rows)
   !< Fold (local) ghost cells back onto the cells they mirror, the transpose of update_ghost_local_gpu_dev.
   !< Device twin of `field_object%reduce_ghost_local`: same-level rows only (`bad_rows` returns the number of 2:1 rows, in
   !< which case nothing is touched). A donor cell may receive from several ghosts (edges, corners): atomic updates, and the
   !< ghost refill runs in a second kernel, after every fold is complete.
   integer(I4P), intent(in)          :: ngc                       !< Ghost cells number.
   integer(I4P), intent(in)          :: v_first                   !< First variable to reduce.
   integer(I4P), intent(in)          :: v_last                    !< Last variable to reduce.
   integer(I8P), intent(in), pointer :: l_map_ghost_cell_gpu(:,:) !< Local map of ghost cells.
   real(R8P),    intent(inout)       :: q_gpu(1:,    &
                                              1-ngc:,&
                                              1-ngc:,&
                                              1-ngc:,1:)          !< Field component to be reduced.
   integer(I4P), intent(out)         :: bad_rows                  !< Number of 2:1 (non same-level) rows.
   integer(I4P)                      :: mf, v, nrows              !< Counters.
   integer(I4P)                      :: b_recv, i_recv, j_recv, k_recv !< Ghost cell.
   integer(I4P)                      :: b_send, i_send, j_send, k_send !< Donor cell.

   bad_rows = 0_I4P
   if (.not.associated(l_map_ghost_cell_gpu)) return
   nrows = size(l_map_ghost_cell_gpu, dim=1)
   !$acc parallel loop independent gang vector DEVICEVAR(l_map_ghost_cell_gpu) reduction(+:bad_rows)
   !$omp OMPLOOP DEVICEPTR(l_map_ghost_cell_gpu) reduction(+:bad_rows)
   do mf=1, nrows
      if (l_map_ghost_cell_gpu(mf,9) /= 1_I8P) bad_rows = bad_rows + 1_I4P
   enddo
   if (bad_rows > 0_I4P) return
   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(l_map_ghost_cell_gpu, q_gpu) &
   !$acc& private(b_send,i_send,j_send,k_send,b_recv,i_recv,j_recv,k_recv)
   !$omp OMPLOOP collapse(2) DEVICEPTR(l_map_ghost_cell_gpu, q_gpu) &
   !$omp& private(b_send,i_send,j_send,k_send,b_recv,i_recv,j_recv,k_recv)
   do v=v_first, v_last
      do mf=1, nrows
         b_send = l_map_ghost_cell_gpu(mf,1) ; b_recv = l_map_ghost_cell_gpu(mf,2)
         i_send = l_map_ghost_cell_gpu(mf,3) ; j_send = l_map_ghost_cell_gpu(mf,4) ; k_send = l_map_ghost_cell_gpu(mf,5)
         i_recv = l_map_ghost_cell_gpu(mf,6) ; j_recv = l_map_ghost_cell_gpu(mf,7) ; k_recv = l_map_ghost_cell_gpu(mf,8)
         !$acc atomic update
         !$omp atomic update
         q_gpu(b_send,i_send,j_send,k_send,v) = q_gpu(b_send,i_send,j_send,k_send,v) + q_gpu(b_recv,i_recv,j_recv,k_recv,v)
      enddo
   enddo
   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(l_map_ghost_cell_gpu, q_gpu) &
   !$acc& private(b_send,i_send,j_send,k_send,b_recv,i_recv,j_recv,k_recv)
   !$omp OMPLOOP collapse(2) DEVICEPTR(l_map_ghost_cell_gpu, q_gpu) &
   !$omp& private(b_send,i_send,j_send,k_send,b_recv,i_recv,j_recv,k_recv)
   do v=v_first, v_last
      do mf=1, nrows
         b_send = l_map_ghost_cell_gpu(mf,1) ; b_recv = l_map_ghost_cell_gpu(mf,2)
         i_send = l_map_ghost_cell_gpu(mf,3) ; j_send = l_map_ghost_cell_gpu(mf,4) ; k_send = l_map_ghost_cell_gpu(mf,5)
         i_recv = l_map_ghost_cell_gpu(mf,6) ; j_recv = l_map_ghost_cell_gpu(mf,7) ; k_recv = l_map_ghost_cell_gpu(mf,8)
         q_gpu(b_recv,i_recv,j_recv,k_recv,v) = q_gpu(b_send,i_send,j_send,k_send,v)
      enddo
   enddo
   endsubroutine reduce_ghost_local_gpu_dev

   subroutine pack_reduce_ghost_mpi_gpu_dev(ngc, v_first, v_last, comm_map_recv_ghost_cell_gpu, recv_buffer_ghost_gpu, q_gpu)
   !< Pack the received ghost cells of variables v_first:v_last into their receive-buffer slots (transpose of
   !< receive_recv_buffer_ghost_gpu_dev); every other slot is zeroed, so superseded (deduplicated) slots carry nothing back.
   integer(I4P),          intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P),          intent(in)    :: v_first                           !< First variable to reduce.
   integer(I4P),          intent(in)    :: v_last                            !< Last variable to reduce.
   integer(I8P), pointer, intent(in)    :: comm_map_recv_ghost_cell_gpu(:,:) !< Comm map, cell information.
   real(R8P),    pointer, intent(inout) :: recv_buffer_ghost_gpu(:)          !< Receive buffer of ghost cells.
   real(R8P),             intent(in)    :: q_gpu(1:,    &
                                                 1-ngc:,&
                                                 1-ngc:,&
                                                 1-ngc:,&
                                                 1:)                         !< Field component to be reduced.
   integer(I4P)                         :: c, rf, v                          !< Counters.

   if (.not.associated(recv_buffer_ghost_gpu)) return
   !$acc parallel loop independent DEVICEVAR(recv_buffer_ghost_gpu)
   !$omp OMPLOOP DEVICEPTR(recv_buffer_ghost_gpu)
   do c=1, size(recv_buffer_ghost_gpu, dim=1)
      recv_buffer_ghost_gpu(c) = 0._R8P
   enddo
   if (.not.associated(comm_map_recv_ghost_cell_gpu)) return
   !$acc parallel loop independent DEVICEVAR(comm_map_recv_ghost_cell_gpu, recv_buffer_ghost_gpu, q_gpu) private(v)
   !$omp OMPLOOP DEVICEPTR(comm_map_recv_ghost_cell_gpu, recv_buffer_ghost_gpu, q_gpu) private(v)
   do rf=1, size(comm_map_recv_ghost_cell_gpu, dim=1)
      v = int(comm_map_recv_ghost_cell_gpu(rf,6), I4P)
      if (v >= v_first .and. v <= v_last) &
         recv_buffer_ghost_gpu(comm_map_recv_ghost_cell_gpu(rf,1)) = q_gpu(comm_map_recv_ghost_cell_gpu(rf,2), &
                                                                           comm_map_recv_ghost_cell_gpu(rf,3), &
                                                                           comm_map_recv_ghost_cell_gpu(rf,4), &
                                                                           comm_map_recv_ghost_cell_gpu(rf,5), v)
   enddo
   endsubroutine pack_reduce_ghost_mpi_gpu_dev

   subroutine fold_reduce_ghost_mpi_gpu_dev(ngc, v_first, v_last, comm_map_send_ghost_cell_gpu, send_buffer_ghost_gpu, &
                                            q_gpu, bad_rows)
   !< Add the slots sent back by the other ranks to their donor cells (transpose of populate_send_buffer_ghost_gpu_dev),
   !< variables v_first:v_last, same-level rows only (`bad_rows` returns the number of 2:1 rows, in which case nothing is
   !< touched). A donor may own several slots (edges, corners): atomic updates.
   integer(I4P),          intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P),          intent(in)    :: v_first                           !< First variable to reduce.
   integer(I4P),          intent(in)    :: v_last                            !< Last variable to reduce.
   integer(I8P), pointer, intent(in)    :: comm_map_send_ghost_cell_gpu(:,:) !< Comm map, cell information.
   real(R8P),    pointer, intent(in)    :: send_buffer_ghost_gpu(:)          !< Send buffer of ghost cells.
   real(R8P),             intent(inout) :: q_gpu(1:,    &
                                                 1-ngc:,&
                                                 1-ngc:,&
                                                 1-ngc:,&
                                                 1:)                         !< Field component to be reduced.
   integer(I4P),          intent(out)   :: bad_rows                          !< Number of 2:1 (non same-level) rows.
   integer(I4P)                         :: sf, nrows, v                      !< Counters.
   integer(I4P)                         :: b, i, j, k                        !< Donor cell.

   bad_rows = 0_I4P
   if (.not.associated(comm_map_send_ghost_cell_gpu)) return
   nrows = size(comm_map_send_ghost_cell_gpu, dim=1)
   !$acc parallel loop independent gang vector DEVICEVAR(comm_map_send_ghost_cell_gpu) reduction(+:bad_rows)
   !$omp OMPLOOP DEVICEPTR(comm_map_send_ghost_cell_gpu) reduction(+:bad_rows)
   do sf=1, nrows
      if (comm_map_send_ghost_cell_gpu(sf,7) /= 1_I8P) bad_rows = bad_rows + 1_I4P
   enddo
   if (bad_rows > 0_I4P) return
   !$acc parallel loop independent gang vector DEVICEVAR(comm_map_send_ghost_cell_gpu, send_buffer_ghost_gpu, q_gpu) &
   !$acc& private(b,i,j,k,v)
   !$omp OMPLOOP DEVICEPTR(comm_map_send_ghost_cell_gpu, send_buffer_ghost_gpu, q_gpu) private(b,i,j,k,v)
   do sf=1, nrows
      v = int(comm_map_send_ghost_cell_gpu(sf,5), I4P)
      if (v < v_first .or. v > v_last) cycle
      b = int(comm_map_send_ghost_cell_gpu(sf,1), I4P) ; i = int(comm_map_send_ghost_cell_gpu(sf,2), I4P)
      j = int(comm_map_send_ghost_cell_gpu(sf,3), I4P) ; k = int(comm_map_send_ghost_cell_gpu(sf,4), I4P)
      !$acc atomic update
      !$omp atomic update
      q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + send_buffer_ghost_gpu(comm_map_send_ghost_cell_gpu(sf,6))
   enddo
   endsubroutine fold_reduce_ghost_mpi_gpu_dev

   ! private procedures
   pure function interp_seam_ghost_gpu_dev(regime, meta, ngc, q_gpu, v, b_send, i_send, j_send, k_send, density_var) &
      result(value_)
   !< Device twin of `field_object%interp_seam_ghost`: the interpolant of `interp_seam_value_gpu_dev`, with the guard of
   !< issue #101 (the density, followed by the three momenta: the density clamped to the MP bound of the coarse cells'
   !< densities, `seam_mp_bound`, the velocity m/rho to the MP bound of their velocities, the momentum rho times it; a
   !< ghost whose clamped density is still below `SEAM_DENSITY_GUARD` times the anchor's takes in every variable the
   !< anchor cell's value, the injection value). Each thread recomputes the density of its ghost, so the rows of one
   !< ghost decide alike whatever their order; the guard reads only the tricubic donor footprint, real cells. The
   !< footprint buffers are constant-bound locals (see `interp_seam_value_gpu_dev`).
   integer(I4P), intent(in) :: regime           !< Active fill regime (SEAM_FILL_*).
   integer(I4P), intent(in) :: meta             !< Packed interpolation metadata (see seam_meta_pack).
   integer(I4P), intent(in) :: ngc              !< Ghost cells number (bounds of q_gpu).
   real(R8P),    intent(in) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field variables, FNL layout (b,i,j,k,v).
   integer(I4P), intent(in) :: v                !< Variable index.
   integer(I4P), intent(in) :: b_send           !< Donor block index.
   integer(I4P), intent(in) :: i_send           !< Anchor donor cell, i.
   integer(I4P), intent(in) :: j_send           !< Anchor donor cell, j.
   integer(I4P), intent(in) :: k_send           !< Anchor donor cell, k.
   integer(I4P), intent(in) :: density_var      !< Density index for the guard, 0 for none (`maps%seam_density_var`).
   real(R8P)                :: value_           !< Fine-ghost value.
   real(R8P)                :: rho              !< Interpolated (then clamped) density of the ghost.
   real(R8P)                :: fp_rho(1:4,1:4,1:4) !< Tricubic footprint of the density.
   real(R8P)                :: fp_u(1:4,1:4,1:4) !< Tricubic footprint of the velocity component of v.
   real(R8P)                :: lo, hi           !< MP bound.
   real(R8P)                :: u                !< Velocity component of v (clamped).
   integer(I4P)             :: sub(1:3)         !< Octant sub-position per direction.
   integer(I4P)             :: p4(1:3)          !< Tricubic anchor position per direction.
   integer(I4P)             :: p3(1:3)          !< Compatible anchor position per direction (unused).
   integer(I4P)             :: p6(1:3)          !< Quintic anchor position per direction (unused).
   integer(I4P)             :: ii, jj, kk       !< Footprint counters.
   integer(I4P)             :: g_lo(1:3), g_hi(1:3) !< Footprint nodes gathered per direction.
   integer(I4P)             :: b_lo(1:3), b_hi(1:3) !< Bracketing nodes per direction (seam_bracket_nodes).
   integer(I4P)             :: n2(1:3)          !< Bracketing nodes (unused here).
   logical                  :: edge(1:3)        !< Edge flags (unused here).
   !$acc routine seq
   !$omp declare target

   value_ = interp_seam_value_gpu_dev(regime=regime, meta=meta, ngc=ngc, q_gpu=q_gpu, v=v, b_send=b_send, &
                                      i_send=i_send, j_send=j_send, k_send=k_send)
   if (density_var > 0_I4P) then
      call seam_meta_unpack(meta=meta, sub=sub, p4=p4, p3=p3, p6=p6)
      ! the footprint nodes seam_mp_bound reads: all 4 along a refined direction, the anchor along an unrefined one
      g_lo = merge(p4, 1_I4P, sub == SEAM_SUB_COINCIDENT)
      g_hi = merge(p4, 4_I4P, sub == SEAM_SUB_COINCIDENT)
      fp_rho = 0._R8P
      fp_u = 0._R8P
      rho = value_
      if (v /= density_var) rho = interp_seam_value_gpu_dev(regime=regime, meta=meta, ngc=ngc, q_gpu=q_gpu,          &
                                                            v=density_var, b_send=b_send, i_send=i_send, j_send=j_send, &
                                                            k_send=k_send)
      do kk=g_lo(3), g_hi(3)
         do jj=g_lo(2), g_hi(2)
            do ii=g_lo(1), g_hi(1)
               fp_rho(ii,jj,kk) = q_gpu(b_send, i_send+ii-p4(1), j_send+jj-p4(2), k_send+kk-p4(3), density_var)
            enddo
         enddo
      enddo
      call seam_bracket_nodes(sub=sub, anchor_pos=p4, n2=n2, edge=edge, b_lo=b_lo, b_hi=b_hi)
      ! the bracketing cells' range lies inside the MP bound: a value inside it is never clamped (the bound skipped)
      call seam_bracket_range(footprint=fp_rho, anchor_pos=p4, b_lo=b_lo, b_hi=b_hi, lo=lo, hi=hi)
      if (.not.(rho >= lo .and. rho <= hi)) then
         call seam_mp_bound(footprint=fp_rho, sub=sub, anchor_pos=p4, lo=lo, hi=hi)
         ! explicit comparisons, not min/max: a NaN is kept alike by every compiler (gfortran and nvfortran differ)
         if (rho < lo) rho = lo
         if (rho > hi) rho = hi
      endif
      if (.not.(rho > SEAM_DENSITY_GUARD * q_gpu(b_send, i_send, j_send, k_send, density_var))) then
         value_ = q_gpu(b_send, i_send, j_send, k_send, v)
      elseif (v == density_var) then
         value_ = rho
      elseif (v > density_var .and. v <= density_var + 3_I4P) then
         u = value_ / rho
         do kk=b_lo(3), b_hi(3)
            do jj=b_lo(2), b_hi(2)
               do ii=b_lo(1), b_hi(1)
                  fp_u(ii,jj,kk) = q_gpu(b_send, i_send+ii-p4(1), j_send+jj-p4(2), k_send+kk-p4(3), v) / fp_rho(ii,jj,kk)
               enddo
            enddo
         enddo
         call seam_bracket_range(footprint=fp_u, anchor_pos=p4, b_lo=b_lo, b_hi=b_hi, lo=lo, hi=hi)
         if (.not.(u >= lo .and. u <= hi)) then
            do kk=g_lo(3), g_hi(3)
               do jj=g_lo(2), g_hi(2)
                  do ii=g_lo(1), g_hi(1)
                     fp_u(ii,jj,kk) = q_gpu(b_send, i_send+ii-p4(1), j_send+jj-p4(2), k_send+kk-p4(3), v) / &
                                      fp_rho(ii,jj,kk)
                  enddo
               enddo
            enddo
            call seam_mp_bound(footprint=fp_u, sub=sub, anchor_pos=p4, lo=lo, hi=hi)
            if (u < lo) u = lo
            if (u > hi) u = hi
         endif
         value_ = rho * u
      endif
   endif
   endfunction interp_seam_ghost_gpu_dev

   pure function interp_seam_value_gpu_dev(regime, meta, ngc, q_gpu, v, b_send, i_send, j_send, k_send) result(value_)
   !< Device twin of `field_object%interp_seam_value` (issue #22 F3): evaluate
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
   integer(I4P), intent(in) :: regime           !< Active fill regime (SEAM_FILL_*).
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
   integer(I4P)             :: p6(1:3)          !< Quintic anchor position per direction.
   integer(I4P)             :: i0, i1, j0, j1, k0, k1 !< Quintic footprint nodes read per direction.
   real(R8P)                :: sum_x, sum_y     !< Quintic partial sums.
   real(R8P)                :: ux(1:6)          !< WENO: an x line of coarse values.
   real(R8P)                :: uy(1:6)          !< WENO: a y line of x-filled values.
   real(R8P)                :: uz(1:6)          !< WENO: the z line of xy-filled values.
   integer(I4P)             :: ii, jj, kk       !< Footprint counters.
   !$acc routine seq
   !$omp declare target

   call seam_meta_unpack(meta=meta, sub=sub, p4=p4, p3=p3, p6=p6)
   if (regime == SEAM_FILL_WENO6 .or. regime == SEAM_FILL_QUINTIC_MP) then
      ! issue #68 D2: one x line gathered at a time (6 values) and filled line by line by seam_fill_weno6_1d or
      ! seam_fill_quintic_mp_1d, the arithmetic of the host's tensor sweeps; an unrefined direction reads the anchor only
      i0 = 1_I4P ; i1 = 6_I4P ; j0 = 1_I4P ; j1 = 6_I4P ; k0 = 1_I4P ; k1 = 6_I4P
      if (sub(1) == SEAM_SUB_COINCIDENT) then ; i0 = p6(1) ; i1 = p6(1) ; endif
      if (sub(2) == SEAM_SUB_COINCIDENT) then ; j0 = p6(2) ; j1 = p6(2) ; endif
      if (sub(3) == SEAM_SUB_COINCIDENT) then ; k0 = p6(3) ; k1 = p6(3) ; endif
      ux = 0._R8P
      uy = 0._R8P
      uz = 0._R8P
      do kk=k0, k1
         do jj=j0, j1
            do ii=i0, i1
               ux(ii) = q_gpu(b_send, i_send+ii-p6(1), j_send+jj-p6(2), k_send+kk-p6(3), v)
            enddo
            if (regime == SEAM_FILL_WENO6) then
               uy(jj) = seam_fill_weno6_1d(u=ux, sub=sub(1), p=p6(1))
            else
               uy(jj) = seam_fill_quintic_mp_1d(u=ux, sub=sub(1), p=p6(1))
            endif
         enddo
         if (regime == SEAM_FILL_WENO6) then
            uz(kk) = seam_fill_weno6_1d(u=uy, sub=sub(2), p=p6(2))
         else
            uz(kk) = seam_fill_quintic_mp_1d(u=uy, sub=sub(2), p=p6(2))
         endif
      enddo
      if (regime == SEAM_FILL_WENO6) then
         value_ = seam_fill_weno6_1d(u=uz, sub=sub(3), p=p6(3))
      else
         value_ = seam_fill_quintic_mp_1d(u=uz, sub=sub(3), p=p6(3))
      endif
   elseif (regime == SEAM_FILL_QUINTIC) then
      ! issue #68 D1: the tensor quintic summed direction by direction (x innermost) straight from q_gpu, the order of
      ! the host's seam_interpolate_quintic, with no per-thread footprint; an unrefined direction reads the anchor only
      i0 = 1_I4P ; i1 = 6_I4P ; j0 = 1_I4P ; j1 = 6_I4P ; k0 = 1_I4P ; k1 = 6_I4P
      if (sub(1) == SEAM_SUB_COINCIDENT) then ; i0 = p6(1) ; i1 = p6(1) ; endif
      if (sub(2) == SEAM_SUB_COINCIDENT) then ; j0 = p6(2) ; j1 = p6(2) ; endif
      if (sub(3) == SEAM_SUB_COINCIDENT) then ; k0 = p6(3) ; k1 = p6(3) ; endif
      value_ = 0._R8P
      do kk=k0, k1
         sum_y = 0._R8P
         do jj=j0, j1
            sum_x = 0._R8P
            do ii=i0, i1
               sum_x = sum_x + seam_quintic_weight(ii, sub(1), p6(1)) * &
                               q_gpu(b_send, i_send+ii-p6(1), j_send+jj-p6(2), k_send+kk-p6(3), v)
            enddo
            sum_y = sum_y + seam_quintic_weight(jj, sub(2), p6(2)) * sum_x
         enddo
         value_ = value_ + seam_quintic_weight(kk, sub(3), p6(3)) * sum_y
      enddo
   elseif (regime == SEAM_FILL_COMPATIBLE) then
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
   endfunction interp_seam_value_gpu_dev

   pure function restrict_seam_high_order_gpu_dev(meta, ngc, q_gpu, v, b_send, i_send, j_send, k_send) result(value_)
   !< Device twin of `field_object%restrict_seam_high_order` (issue #68 C2, C3): the 6-point Lagrange or the WENO
   !< restriction of one (variable, coarse ghost) pair on the FNL layout `q_gpu(b,i,j,k,v)`, in the same order as the
   !< host (x innermost, then y, then z). The footprint is shifted inside the donor block by the codes of `meta`, so it
   !< reads REAL donor cells only and never a ghost cell another map row writes concurrently (the #22 race discipline).
   !< The Lagrange sums use scalar locals only; the WENO gathers one x line at a time (6 values) and calls
   !< `seam_restriction_weno6_1d` line by line, the arithmetic of the host's `seam_restriction_weno6`.
   integer(I4P), intent(in) :: meta             !< Packed restriction codes and regime (seam_restriction_meta_pack).
   integer(I4P), intent(in) :: ngc              !< Ghost cells number (bounds of q_gpu).
   real(R8P),    intent(in) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field variables, FNL layout (b,i,j,k,v).
   integer(I4P), intent(in) :: v                !< Variable index.
   integer(I4P), intent(in) :: b_send           !< Donor block index.
   integer(I4P), intent(in) :: i_send           !< First fine cell of the restricted pair, i.
   integer(I4P), intent(in) :: j_send           !< First fine cell of the restricted pair, j.
   integer(I4P), intent(in) :: k_send           !< First fine cell of the restricted pair, k.
   real(R8P)                :: value_           !< Restricted coarse-ghost value.
   integer(I4P)             :: code(1:3)        !< Restriction code per direction.
   integer(I4P)             :: regime           !< Restriction regime.
   integer(I4P)             :: ni, nj, nk       !< Footprint nodes per direction (6, or 1 unrefined).
   integer(I4P)             :: si, sj, sk       !< First footprint cell per direction.
   real(R8P)                :: sum_x, sum_y     !< Partial sums.
   real(R8P)                :: ux(1:6)          !< An x line of fine values (WENO).
   real(R8P)                :: uy(1:6)          !< A y line of x-restricted values (WENO).
   real(R8P)                :: uz(1:6)          !< The z line of xy-restricted values (WENO).
   integer(I4P)             :: ii, jj, kk       !< Footprint counters.
   !$acc routine seq
   !$omp declare target

   call seam_restriction_meta_unpack(meta=meta, code=code, regime=regime)
   ni = 6_I4P ; si = i_send - code(1)
   nj = 6_I4P ; sj = j_send - code(2)
   nk = 6_I4P ; sk = k_send - code(3)
   if (code(1) == SEAM_RESTRICT_IDENTITY) then ; ni = 1_I4P ; si = i_send ; endif
   if (code(2) == SEAM_RESTRICT_IDENTITY) then ; nj = 1_I4P ; sj = j_send ; endif
   if (code(3) == SEAM_RESTRICT_IDENTITY) then ; nk = 1_I4P ; sk = k_send ; endif
   if (regime == SEAM_RESTRICT_WENO6) then
      ux = 0._R8P
      uy = 0._R8P
      uz = 0._R8P
      do kk=1, nk
         do jj=1, nj
            do ii=1, ni
               ux(ii) = q_gpu(b_send, si+ii-1, sj+jj-1, sk+kk-1, v)
            enddo
            uy(jj) = seam_restriction_weno6_1d(u=ux, code=code(1))
         enddo
         uz(kk) = seam_restriction_weno6_1d(u=uy, code=code(2))
      enddo
      value_ = seam_restriction_weno6_1d(u=uz, code=code(3))
      return
   endif
   value_ = 0._R8P
   do kk=1, nk
      sum_y = 0._R8P
      do jj=1, nj
         sum_x = 0._R8P
         do ii=1, ni
            sum_x = sum_x + seam_restriction_weight(ii, code(1)) * q_gpu(b_send, si+ii-1, sj+jj-1, sk+kk-1, v)
         enddo
         sum_y = sum_y + seam_restriction_weight(jj, code(2)) * sum_x
      enddo
      value_ = value_ + seam_restriction_weight(kk, code(3)) * sum_y
   enddo
   endfunction restrict_seam_high_order_gpu_dev
endmodule adam_fnl_field_kernels
