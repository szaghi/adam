!< ADAM, FLUME FNL device kernels.

#include "fundal.H"

module adam_flume_fnl_kernels
!< ADAM, FLUME FNL device kernels.
!<
!< Model-independent kernels; the per-model ones (face fluxes, auxiliary variables, signal speed) live in
!< `adam_flume_fnl_<model>_kernels` (issue #41, section 4).
!<
!< Module-level kernels (issue #35, D-11/D-12): private work arrays have constant bounds, no array section is ever an
!< actual argument inside a loop body, host scalars are `firstprivate`, device arrays are `DEVICEVAR`/`DEVICEPTR`,
!< and every offload `!$acc` directive is immediately followed by its `!$omp` twin. Device arrays are transposed,
!< `(b, i, j, k, v)`: the block index is stride-1, so it is the innermost loop of every `collapse(4)` nest.

! ADAM classes, libraries, parameters
use :: adam_parameters,           only : FEC_1_6_ARRAY
! FLUME modules
use :: adam_flume_common_library, only : ib_cut_spacing, seam_skin_cell, BC_EXTRAPOLATION, BC_INFLOW, BC_WALL_INVISCID, &
                                         IQ_RU
! third party modules
use :: penf,                      only : I4P, I8P, R8P

implicit none
private
public :: apply_reflux_face_dev
public :: compute_flux_difference_dev
public :: compute_flux_difference_ib_dev
public :: compute_rk_ssp_residual_dev
public :: fill_seam_copy_dev
public :: pack_seam_skin_dev
public :: set_boundary_conditions_dev

contains
   ! public procedures
   subroutine apply_reflux_face_dev(axis, sgn, b, ni, nj, nk, ngc, nv, nface_cells, scale, delta_gpu, q_gpu)
   !< Add the Berger-Colella correction `scale delta(:, c)` of one register face to the coarse skin cells of block `b`.
   integer(I4P), intent(in)    :: axis, sgn                         !< Face normal axis (1..3) and side (+-1).
   integer(I4P), intent(in)    :: b                                 !< Coarse block.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                                !< Variables number.
   integer(I4P), intent(in)    :: nface_cells                       !< Skin cells number.
   real(R8P),    intent(in)    :: scale                             !< Correction scale, sgn dt / dx_coarse.
   real(R8P),    intent(in)    :: delta_gpu(1:,1:)                  !< Flux mismatch F_coarse - F_fine_sum (nv, cells).
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   integer(I4P)                :: c, i, j, k, v                     !< Counters and skin cell indexes.

   !$acc parallel loop independent gang vector DEVICEVAR(delta_gpu,q_gpu) &
   !$acc& firstprivate(axis,sgn,b,ni,nj,nk,nv,nface_cells,scale) private(i,j,k)
   !$omp OMPLOOP DEVICEPTR(delta_gpu,q_gpu) &
   !$omp& firstprivate(axis,sgn,b,ni,nj,nk,nv,nface_cells,scale) private(i,j,k)
   do c=1, nface_cells
      call seam_skin_cell(axis=axis, sgn=sgn, ni=ni, nj=nj, nk=nk, c=c, i=i, j=j, k=k)
      !$acc loop seq
      do v=1, nv
         q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + scale * delta_gpu(v,c)
      enddo
   enddo
   endsubroutine apply_reflux_face_dev

   subroutine compute_flux_difference_dev(nv, ni, nj, nk, ngc, blocks_number, is_null, dxyz_gpu, flx_f_gpu, fly_f_gpu, &
                                          flz_f_gpu, dq_gpu)
   !< Compute the residuals from the face fluxes, `dq = -sum_d (F_{d,i+1/2} - F_{d,i-1/2}) / dx_d`.
   !<
   !< A null direction weighs zero, and its normal momentum residual is zero (CHASE semantics, issue #35, section 3.4).
   integer(I4P), intent(in)    :: nv                                 !< Conservative variables number.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   logical,      intent(in)    :: is_null(3)                         !< Null directions.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                    !< Blocks space steps [nb, 3].
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)          !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)          !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)          !< Z-face fluxes.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   real(R8P)                   :: wx, wy, wz                         !< Direction weights: 1 active, 0 null.
   logical                     :: nx, ny, nz                         !< Null directions, scalar copies.
   integer(I4P)                :: b, i, j, k, v                      !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1)) ; nx = is_null(1)
   wy = merge(0._R8P, 1._R8P, is_null(2)) ; ny = is_null(2)
   wz = merge(0._R8P, 1._R8P, is_null(3)) ; nz = is_null(3)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,dq_gpu) &
   !$acc& firstprivate(nv,ni,nj,nk,blocks_number,wx,wy,wz,nx,ny,nz)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,dq_gpu) &
   !$omp& firstprivate(nv,ni,nj,nk,blocks_number,wx,wy,wz,nx,ny,nz)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do v=1, nv
         dq_gpu(b,i,j,k,v) = -(wx * (flx_f_gpu(b,i,j,k,v) - flx_f_gpu(b,i-1,j,k,v)) / dxyz_gpu(b,1) + &
                               wy * (fly_f_gpu(b,i,j,k,v) - fly_f_gpu(b,i,j-1,k,v)) / dxyz_gpu(b,2) + &
                               wz * (flz_f_gpu(b,i,j,k,v) - flz_f_gpu(b,i,j,k-1,v)) / dxyz_gpu(b,3))
      enddo
      if (nx) dq_gpu(b,i,j,k,IQ_RU  ) = 0._R8P
      if (ny) dq_gpu(b,i,j,k,IQ_RU+1) = 0._R8P
      if (nz) dq_gpu(b,i,j,k,IQ_RU+2) = 0._R8P
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_flux_difference_dev

   subroutine compute_flux_difference_ib_dev(nv, ni, nj, nk, ngc, blocks_number, is_null, dxyz_gpu, flx_f_gpu, &
                                             fly_f_gpu, flz_f_gpu, phi_gpu, dq_gpu)
   !< Compute the residuals from the face fluxes with immersed solids: the spacing of a fluid cell is cut by the solid
   !< surface (`ib_cut_spacing`, CHASE semantics, issue #35 D-9); device twin of the CPU flux difference with `phi`.
   integer(I4P), intent(in)    :: nv                                  !< Conservative variables number.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                     !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                       !< Actual blocks number.
   logical,      intent(in)    :: is_null(3)                          !< Null directions.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                     !< Blocks space steps [nb, 3].
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)           !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)           !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)           !< Z-face fluxes.
   real(R8P),    intent(in)    :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function [nb, i, j, k, solids+1].
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Residuals.
   real(R8P), parameter        :: IB_EPS=1.e-12_R8P                   !< Guard of the cut spacing (CHASE value).
   real(R8P)                   :: wx, wy, wz                          !< Direction weights: 1 active, 0 null.
   real(R8P)                   :: dx, dy, dz                          !< Cell spacings.
   logical                     :: nx, ny, nz                          !< Null directions, scalar copies.
   integer(I4P)                :: ns                                  !< All-solids summary slot of phi.
   integer(I4P)                :: b, i, j, k, v                       !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1)) ; nx = is_null(1)
   wy = merge(0._R8P, 1._R8P, is_null(2)) ; ny = is_null(2)
   wz = merge(0._R8P, 1._R8P, is_null(3)) ; nz = is_null(3)
   ns = ubound(phi_gpu, dim=5)
   !$acc parallel loop independent gang vector collapse(4)                                    &
   !$acc& DEVICEVAR(dxyz_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,phi_gpu,dq_gpu)                     &
   !$acc& firstprivate(nv,ni,nj,nk,blocks_number,wx,wy,wz,nx,ny,nz,ns) private(dx,dy,dz)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,phi_gpu,dq_gpu) &
   !$omp& firstprivate(nv,ni,nj,nk,blocks_number,wx,wy,wz,nx,ny,nz,ns) private(dx,dy,dz)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      dx = ib_cut_spacing(phi_c=phi_gpu(b,i,j,k,ns), phi_m=phi_gpu(b,i-1,j,k,ns), phi_p=phi_gpu(b,i+1,j,k,ns), &
                          ds=dxyz_gpu(b,1), eps=IB_EPS)
      dy = ib_cut_spacing(phi_c=phi_gpu(b,i,j,k,ns), phi_m=phi_gpu(b,i,j-1,k,ns), phi_p=phi_gpu(b,i,j+1,k,ns), &
                          ds=dxyz_gpu(b,2), eps=IB_EPS)
      dz = ib_cut_spacing(phi_c=phi_gpu(b,i,j,k,ns), phi_m=phi_gpu(b,i,j,k-1,ns), phi_p=phi_gpu(b,i,j,k+1,ns), &
                          ds=dxyz_gpu(b,3), eps=IB_EPS)
      !$acc loop seq
      do v=1, nv
         dq_gpu(b,i,j,k,v) = -(wx * (flx_f_gpu(b,i,j,k,v) - flx_f_gpu(b,i-1,j,k,v)) / dx + &
                               wy * (fly_f_gpu(b,i,j,k,v) - fly_f_gpu(b,i,j-1,k,v)) / dy + &
                               wz * (flz_f_gpu(b,i,j,k,v) - flz_f_gpu(b,i,j,k-1,v)) / dz)
      enddo
      if (nx) dq_gpu(b,i,j,k,IQ_RU  ) = 0._R8P
      if (ny) dq_gpu(b,i,j,k,IQ_RU+1) = 0._R8P
      if (nz) dq_gpu(b,i,j,k,IQ_RU+2) = 0._R8P
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_flux_difference_ib_dev

   subroutine compute_rk_ssp_residual_dev(ni, nj, nk, ngc, nv, blocks_number, nrk, beta_gpu, q_rk_gpu, dq_gpu)
   !< Compute the residual of a strong stability preserving step, `dq = sum_s beta_s dq_s`, from the stored stages.
   !<
   !< Device twin of the residual that the host `rk_object%update_q` returns: the residuals history of the two backends
   !< then reports the same quantity (the effective step residual, not the last stage one).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                         !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                                      !< Variables number.
   integer(I4P), intent(in)    :: blocks_number                           !< Actual blocks number.
   integer(I4P), intent(in)    :: nrk                                     !< Runge-Kutta stages number.
   real(R8P),    intent(in)    :: beta_gpu(1:)                            !< Runge-Kutta beta coefficients.
   real(R8P),    intent(in)    :: q_rk_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:,1:) !< Stored stage residuals.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Step residual.
   real(R8P)                   :: dq                                      !< Private step residual of one cell.
   integer(I4P)                :: b, i, j, k, v, s                        !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(beta_gpu,q_rk_gpu,dq_gpu) &
   !$acc& firstprivate(ni,nj,nk,nv,blocks_number,nrk) private(dq)
   !$omp OMPLOOP collapse(4) DEVICEPTR(beta_gpu,q_rk_gpu,dq_gpu) &
   !$omp& firstprivate(ni,nj,nk,nv,blocks_number,nrk) private(dq)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do v=1, nv
         dq = 0._R8P
         !$acc loop seq
         do s=1, nrk
            dq = dq + beta_gpu(s) * q_rk_gpu(b,i,j,k,v,s)
         enddo
         dq_gpu(b,i,j,k,v) = dq
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_rk_ssp_residual_dev

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

   subroutine pack_seam_skin_dev(fec, b, ni, nj, nk, nv, flx_f_gpu, fly_f_gpu, flz_f_gpu, skin_gpu)
   !< Pack the face fluxes of face `fec` of block `b` into the skin `skin_gpu(v, c)`, `c` running over the two
   !< tangential axes, inner fastest (the flux register order, `accumulate_seam_skin`).
   integer(I4P), intent(in)    :: fec                       !< Face (1..6: -x, +x, -y, +y, -z, +z).
   integer(I4P), intent(in)    :: b                         !< Block.
   integer(I4P), intent(in)    :: ni, nj, nk                !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                        !< Variables number.
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:) !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:) !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:) !< Z-face fluxes.
   real(R8P),    intent(inout) :: skin_gpu(1:,1:)           !< Face skin (nv, inner_n*outer_n).
   integer(I4P)                :: inner, outer, v           !< Counters.
   integer(I4P)                :: inner_n, outer_n          !< Tangential cells numbers.
   integer(I4P)                :: n                         !< Face normal index (0 or n).

   select case(fec)
   case(1_I4P, 2_I4P)
      inner_n = nj ; outer_n = nk ; n = merge(0_I4P, ni, fec == 1_I4P)
   case(3_I4P, 4_I4P)
      inner_n = ni ; outer_n = nk ; n = merge(0_I4P, nj, fec == 3_I4P)
   case default
      inner_n = ni ; outer_n = nj ; n = merge(0_I4P, nk, fec == 5_I4P)
   endselect
   !$acc parallel loop independent gang vector collapse(2) DEVICEVAR(flx_f_gpu,fly_f_gpu,flz_f_gpu,skin_gpu) &
   !$acc& firstprivate(fec,b,nv,n,inner_n,outer_n)
   !$omp OMPLOOP collapse(2) DEVICEPTR(flx_f_gpu,fly_f_gpu,flz_f_gpu,skin_gpu) &
   !$omp& firstprivate(fec,b,nv,n,inner_n,outer_n)
   do outer=1, outer_n
   do inner=1, inner_n
      !$acc loop seq
      do v=1, nv
         select case(fec)
         case(1_I4P, 2_I4P)
            skin_gpu(v,(outer-1)*inner_n+inner) = flx_f_gpu(b,n,inner,outer,v)
         case(3_I4P, 4_I4P)
            skin_gpu(v,(outer-1)*inner_n+inner) = fly_f_gpu(b,inner,n,outer,v)
         case default
            skin_gpu(v,(outer-1)*inner_n+inner) = flz_f_gpu(b,inner,outer,n,v)
         endselect
      enddo
   enddo
   enddo
   endsubroutine pack_seam_skin_dev

   subroutine set_boundary_conditions_dev(ni, nj, nk, ngc, nv, crown, local_map_bc_crown_gpu, q_inflow_gpu, wall_sign_gpu, &
                                          q_gpu)
   !< Set boundary conditions on one crown: face ghosts by kind, edge and corner ghosts by extrapolation.
   !<
   !< Crowns are processed in order by the caller, so every source cell is interior or in a lower (already filled)
   !< crown: the rows of one launch are independent.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                                !< Conservative variables number.
   integer(I4P), intent(in)    :: crown                             !< Crown counter.
   integer(I8P), intent(in)    :: local_map_bc_crown_gpu(:,:,:)     !< Boundary crown map (row, field, crown).
   real(R8P),    intent(in)    :: q_inflow_gpu(1:,1:)               !< Conservative inflow state of each face [nv, 6].
   real(R8P),    intent(in)    :: wall_sign_gpu(1:,1:)              !< Wall mirror sign per variable and direction [nv, 3].
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   integer(I4P)                :: b, c, i, j, k, v                  !< Counters.
   integer(I4P)                :: idelta, jdelta, kdelta            !< IJK inward step.
   integer(I4P)                :: bc_type                           !< Boundary condition type.
   integer(I4P)                :: fec                               !< Boundary fec (1 to 26).
   integer(I4P)                :: face                              !< Boundary face (1 to 6).
   integer(I4P)                :: iref, jref, kref                  !< Donor indexes.

   !$acc parallel loop independent gang vector DEVICEVAR(local_map_bc_crown_gpu,q_inflow_gpu,wall_sign_gpu,q_gpu) &
   !$acc& firstprivate(ni,nj,nk,nv,crown)                                                                   &
   !$acc& private(b,i,j,k,idelta,jdelta,kdelta,bc_type,fec,face,iref,jref,kref)
   !$omp OMPLOOP DEVICEPTR(local_map_bc_crown_gpu,q_inflow_gpu,wall_sign_gpu,q_gpu) &
   !$omp& firstprivate(ni,nj,nk,nv,crown) &
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
                  q_gpu(b,i,j,k,v) = q_inflow_gpu(v,face)
               enddo
            elseif (bc_type == BC_EXTRAPOLATION) then
               !$acc loop seq
               do v=1, nv
                  q_gpu(b,i,j,k,v) = q_gpu(b,iref,jref,kref,v)
               enddo
            elseif (bc_type == BC_WALL_INVISCID) then
               !$acc loop seq
               do v=1, nv
                  q_gpu(b,i,j,k,v) = wall_sign_gpu(v,(face+1)/2) * q_gpu(b,iref,jref,kref,v)
               enddo
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
endmodule adam_flume_fnl_kernels
