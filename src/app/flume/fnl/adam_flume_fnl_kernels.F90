!< ADAM, FLUME FNL device kernels.

#include "fundal.H"

module adam_flume_fnl_kernels
!< ADAM, FLUME FNL device kernels.
!<
!< Module-level kernels (issue #35, D-11/D-12): private work arrays have constant bounds, no array section is ever an
!< actual argument inside a loop body, host scalars are `firstprivate`, device arrays are `DEVICEVAR`/`DEVICEPTR`,
!< and every offload `!$acc` directive is immediately followed by its `!$omp` twin. Device arrays are transposed,
!< `(b, i, j, k, v)`: the block index is stride-1, so it is the innermost loop of every `collapse(4)` nest.

! ADAM classes, libraries, parameters
use :: adam_parameters,           only : FEC_1_6_ARRAY
! FLUME modules
use :: adam_flume_common_library, only : conservative_to_auxiliary, BC_EXTRAPOLATION, BC_INFLOW, BC_WALL_INVISCID, IA_A, &
                                         IA_U, IA_V, IA_W, IQ_RU, NV_AUX, NV_EULER
! third party modules
use :: penf,                      only : I4P, I8P, R8P

implicit none
private
public :: compute_conservation_dev
public :: compute_lambda_max_dev
public :: compute_q_aux_dev
public :: fill_seam_copy_dev
public :: set_boundary_conditions_dev
public :: set_zero_dev

contains
   ! public procedures
   subroutine compute_conservation_dev(ni, nj, nk, ngc, blocks_number, dxyz_gpu, is_null, q_gpu, integrals)
   !< Compute the volume integrals of the conservative variables (interior cells, null directions excluded).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Blocks space steps [nb, 3].
   logical,      intent(in)  :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: integrals(NV_EULER)               !< Volume integrals.
   real(R8P)                 :: wx, wy, wz                        !< Direction weights: 1 active, 0 null.
   real(R8P)                 :: volume                            !< Cell volume.
   real(R8P)                 :: s1, s2, s3, s4, s5                !< Reduction accumulators.
   integer(I4P)              :: b, i, j, k                        !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   s1 = 0._R8P ; s2 = 0._R8P ; s3 = 0._R8P ; s4 = 0._R8P ; s5 = 0._R8P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,wx,wy,wz) private(volume)          &
   !$acc& reduction(+:s1,s2,s3,s4,s5)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,wx,wy,wz) private(volume) &
   !$omp& reduction(+:s1,s2,s3,s4,s5)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      volume = (wx * dxyz_gpu(b,1) + 1._R8P - wx) * (wy * dxyz_gpu(b,2) + 1._R8P - wy) * &
               (wz * dxyz_gpu(b,3) + 1._R8P - wz)
      s1 = s1 + q_gpu(b,i,j,k,1) * volume
      s2 = s2 + q_gpu(b,i,j,k,2) * volume
      s3 = s3 + q_gpu(b,i,j,k,3) * volume
      s4 = s4 + q_gpu(b,i,j,k,4) * volume
      s5 = s5 + q_gpu(b,i,j,k,5) * volume
   enddo
   enddo
   enddo
   enddo
   integrals = [s1, s2, s3, s4, s5]
   endsubroutine compute_conservation_dev

   subroutine compute_lambda_max_dev(ni, nj, nk, ngc, blocks_number, gamma, R, dxyz_gpu, is_null, q_gpu, lambda_max)
   !< Compute `max(sum_d (|u_d| + a) / dx_d)` over the interior cells (null directions excluded).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)  :: gamma                             !< Specific heats ratio.
   real(R8P),    intent(in)  :: R                                 !< Gas constant.
   real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Blocks space steps [nb, 3].
   logical,      intent(in)  :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: lambda_max                        !< Maximum of sum_d (|u_d| + a) / dx_d.
   real(R8P)                 :: wx, wy, wz                        !< Direction weights: 1 active, 0 null.
   real(R8P)                 :: q_(NV_EULER)                      !< Private conservative variables of one cell.
   real(R8P)                 :: qa_(NV_AUX)                       !< Private auxiliary variables of one cell.
   real(R8P)                 :: lambda                            !< Cell value.
   integer(I4P)              :: b, i, j, k, v                     !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   lambda_max = 0._R8P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,gamma,R,wx,wy,wz) private(q_,qa_,lambda) &
   !$acc& reduction(max:lambda_max)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,gamma,R,wx,wy,wz) private(q_,qa_,lambda) &
   !$omp& reduction(max:lambda_max)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do v=1, NV_EULER
         q_(v) = q_gpu(b,i,j,k,v)
      enddo
      call conservative_to_auxiliary(gamma=gamma, R=R, q=q_, qa=qa_)
      lambda = wx * (abs(qa_(IA_U)) + qa_(IA_A)) / dxyz_gpu(b,1) + &
               wy * (abs(qa_(IA_V)) + qa_(IA_A)) / dxyz_gpu(b,2) + &
               wz * (abs(qa_(IA_W)) + qa_(IA_A)) / dxyz_gpu(b,3)
      lambda_max = max(lambda_max, lambda)
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_lambda_max_dev

   subroutine compute_q_aux_dev(ni, nj, nk, ngc, blocks_number, gamma, R, q_gpu, q_aux_gpu)
   !< Compute the auxiliary variables on every cell, ghost cells included.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                                 !< Specific heats ratio.
   real(R8P),    intent(in)    :: R                                     !< Gas constant.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   real(R8P),    intent(inout) :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P)                   :: q_(NV_EULER)                          !< Private conservative variables of one cell.
   real(R8P)                   :: qa_(NV_AUX)                           !< Private auxiliary variables of one cell.
   integer(I4P)                :: b, i, j, k, v                         !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,q_aux_gpu) &
   !$acc& firstprivate(ni,nj,nk,ngc,blocks_number,gamma,R) private(q_,qa_)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu,q_aux_gpu) &
   !$omp& firstprivate(ni,nj,nk,ngc,blocks_number,gamma,R) private(q_,qa_)
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, blocks_number
      !$acc loop seq
      do v=1, NV_EULER
         q_(v) = q_gpu(b,i,j,k,v)
      enddo
      call conservative_to_auxiliary(gamma=gamma, R=R, q=q_, qa=qa_)
      !$acc loop seq
      do v=1, NV_AUX
         q_aux_gpu(b,i,j,k,v) = qa_(v)
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_q_aux_dev

   subroutine fill_seam_copy_dev(row_start, row_count, nv, ngc, rows_gpu, src_gpu, dst_gpu)
   !< Copy a peer realm's interior cells into this realm's seam ghost cells, following the seam map rows.
   integer(I4P), intent(in)    :: row_start, row_count                !< Seam map slice of the peer.
   integer(I4P), intent(in)    :: nv                                  !< Variables number.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: rows_gpu(1:,1:)                     !< Seam map rows.
   real(R8P),    intent(in)    :: src_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Peer active buffer.
   real(R8P),    intent(inout) :: dst_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Own active buffer.
   integer(I4P)                :: c, row, v                           !< Counters.

   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(rows_gpu,src_gpu,dst_gpu) &
   !$acc& firstprivate(row_start,row_count,nv) private(row)
   !$omp OMPLOOP collapse(2) DEVICEPTR(rows_gpu,src_gpu,dst_gpu) &
   !$omp& firstprivate(row_start,row_count,nv) private(row)
   do c=1, row_count
   do v=1, nv
      row = row_start + c - 1
      dst_gpu(rows_gpu(row,3),rows_gpu(row,7),rows_gpu(row,8),rows_gpu(row,9),v) = &
         src_gpu(rows_gpu(row,2),rows_gpu(row,4),rows_gpu(row,5),rows_gpu(row,6),v)
   enddo
   enddo
   endsubroutine fill_seam_copy_dev

   subroutine set_boundary_conditions_dev(ni, nj, nk, ngc, nv, crown, local_map_bc_crown_gpu, q_inflow, q_gpu)
   !< Set boundary conditions on one crown: face ghosts by kind, edge and corner ghosts by extrapolation.
   !<
   !< Crowns are processed in order by the caller, so every source cell is interior or in a lower (already filled)
   !< crown: the rows of one launch are independent.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                                !< Conservative variables number.
   integer(I4P), intent(in)    :: crown                             !< Crown counter.
   integer(I8P), intent(in)    :: local_map_bc_crown_gpu(:,:,:)     !< Boundary crown map (row, field, crown).
   real(R8P),    intent(in)    :: q_inflow(NV_EULER,6)              !< Conservative inflow state of each face.
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   integer(I4P)                :: b, c, i, j, k, v                  !< Counters.
   integer(I4P)                :: idelta, jdelta, kdelta            !< IJK inward step.
   integer(I4P)                :: bc_type                           !< Boundary condition type.
   integer(I4P)                :: fec                               !< Boundary fec (1 to 26).
   integer(I4P)                :: face                              !< Boundary face (1 to 6).
   integer(I4P)                :: iref, jref, kref                  !< Donor indexes.

   !$acc parallel loop independent gang vector DEVICEVAR(local_map_bc_crown_gpu,q_gpu) &
   !$acc& firstprivate(ni,nj,nk,nv,crown,q_inflow)                                    &
   !$acc& private(b,i,j,k,idelta,jdelta,kdelta,bc_type,fec,face,iref,jref,kref)
   !$omp OMPLOOP DEVICEPTR(local_map_bc_crown_gpu,q_gpu) &
   !$omp& firstprivate(ni,nj,nk,nv,crown,q_inflow) &
   !$omp& private(b,i,j,k,idelta,jdelta,kdelta,bc_type,fec,face,iref,jref,kref)
   do c=1, size(local_map_bc_crown_gpu, dim=1)
      b = int(local_map_bc_crown_gpu(c,1,crown), I4P)
      if (b > 0_I4P) then
         i       = int(local_map_bc_crown_gpu(c,2,crown), I4P)
         j       = int(local_map_bc_crown_gpu(c,3,crown), I4P)
         k       = int(local_map_bc_crown_gpu(c,4,crown), I4P)
         idelta  = int(local_map_bc_crown_gpu(c,5,crown), I4P)
         jdelta  = int(local_map_bc_crown_gpu(c,6,crown), I4P)
         kdelta  = int(local_map_bc_crown_gpu(c,7,crown), I4P)
         bc_type = int(local_map_bc_crown_gpu(c,8,crown), I4P)
         fec     = int(local_map_bc_crown_gpu(c,9,crown), I4P)
         iref = i - idelta ; jref = j - jdelta ; kref = k - kdelta
         if (fec <= 6_I4P) then
            face = FEC_1_6_ARRAY(fec)
            if (bc_type == BC_WALL_INVISCID) then
               select case(face)
               case(1)
                  iref = 1_I4P - i
               case(2)
                  iref = 2_I4P * ni + 1_I4P - i
               case(3)
                  jref = 1_I4P - j
               case(4)
                  jref = 2_I4P * nj + 1_I4P - j
               case(5)
                  kref = 1_I4P - k
               case(6)
                  kref = 2_I4P * nk + 1_I4P - k
               endselect
            endif
            if (bc_type == BC_INFLOW) then
               !$acc loop seq
               do v=1, nv
                  q_gpu(b,i,j,k,v) = q_inflow(v,face)
               enddo
            elseif (bc_type == BC_EXTRAPOLATION .or. bc_type == BC_WALL_INVISCID) then
               !$acc loop seq
               do v=1, nv
                  q_gpu(b,i,j,k,v) = q_gpu(b,iref,jref,kref,v)
               enddo
               if (bc_type == BC_WALL_INVISCID) &
                  q_gpu(b,i,j,k,IQ_RU+(face-1)/2) = -q_gpu(b,i,j,k,IQ_RU+(face-1)/2)
            endif
         else
            !$acc loop seq
            do v=1, nv
               q_gpu(b,i,j,k,v) = q_gpu(b,iref,jref,kref,v)
            enddo
         endif
      endif
   enddo
   endsubroutine set_boundary_conditions_dev

   subroutine set_zero_dev(ni, nj, nk, ngc, nv, blocks_number, q_gpu)
   !< Set a device field to zero on every cell, ghost cells included.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                                !< Variables number.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field.
   integer(I4P)                :: b, i, j, k, v                     !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu) firstprivate(ni,nj,nk,ngc,nv,blocks_number)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu) firstprivate(ni,nj,nk,ngc,nv,blocks_number)
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, blocks_number
      !$acc loop seq
      do v=1, nv
         q_gpu(b,i,j,k,v) = 0._R8P
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine set_zero_dev
endmodule adam_flume_fnl_kernels
