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
use :: adam_parameters,           only : FEC_1_6_ARRAY, FEC_TO_DELTA
! FLUME modules
use :: adam_flume_common_library, only : reflux_positivity_theta,                                                 &
                                         ib_cut_spacing, seam_skin_cell, BC_EXTRAPOLATION, BC_INFLOW, BC_WALL_INVISCID, &
                                         BC_WALL_ISOTHERMAL, BC_WALL_NOSLIP, realm_edge_donor, realm_edge_face,           &
                                         wall_noslip_ghost, BC_HYDROSTATIC, hydrostatic_ghost
use :: adam_flume_parameters,     only : IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW
! third party modules
use :: penf,                      only : I4P, I8P, R8P

implicit none
private
public :: apply_reflux_face_dev
public :: compute_flux_difference_dev
public :: add_gravity_sources_dev
public :: add_potential_sources_dev
public :: compute_gravity_correction_faces_dev
public :: subtract_equilibrium_residual_dev
public :: add_gravity_limited_correction_dev
public :: compute_dflux_norms_dev
public :: copy_state_dev
public :: compute_flux_difference_ib_dev
public :: compute_rk_ssp_residual_dev
public :: fill_seam_copy_dev
public :: gather_seam_cells_dev
public :: gather_seam_faces_dev
public :: gather_seam_stencils_dev
public :: pack_density_dev
public :: scatter_seam_cells_dev
public :: scatter_seam_faces_dev
public :: set_boundary_conditions_dev
public :: pack_seam_skins_dev
public :: poison_ghosts_dev
public :: unpack_potential_dev

contains
   ! public procedures
   subroutine apply_reflux_face_dev(axis, sgn, b, ni, nj, nk, ngc, nv, nface_cells, scale, limit, mag_w, psi_w, &
                                    delta_gpu, q_gpu, nclip)
   !< Add the Berger-Colella correction `scale delta(:, c)` of one register face to the coarse skin cells of block `b`.
   !< With `limit` (the positivity limiter, issue #88) each cell receives the largest admissible fraction of its
   !< correction (`reflux_positivity_theta`), exactly the full one whenever admissible; `nclip` counts the others.
   integer(I4P), intent(in)    :: axis, sgn                         !< Face normal axis (1..3) and side (+-1).
   integer(I4P), intent(in)    :: b                                 !< Coarse block.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                                !< Variables number.
   integer(I4P), intent(in)    :: nface_cells                       !< Skin cells number.
   real(R8P),    intent(in)    :: scale                             !< Correction scale, sgn dt / dx_coarse.
   logical,      intent(in)    :: limit                             !< Keep the corrected cells admissible.
   real(R8P),    intent(in)    :: mag_w, psi_w                      !< Magnetic and psi energy weights of the model.
   real(R8P),    intent(in)    :: delta_gpu(1:,1:)                  !< Flux mismatch F_coarse - F_fine_sum (nv, cells).
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   integer(I4P), intent(out)   :: nclip                             !< Cells whose correction was limited.
   integer(I4P)                :: c, i, j, k, v                     !< Counters and skin cell indexes.
   real(R8P)                   :: qc(9), dc(9)                      !< Cell state and full correction (nv <= 9).
   real(R8P)                   :: theta                             !< Admissible fraction of the correction.

   nclip = 0_I4P
   !$acc parallel loop independent gang vector DEVICEVAR(delta_gpu,q_gpu) &
   !$acc& firstprivate(axis,sgn,b,ni,nj,nk,nv,nface_cells,scale,limit,mag_w,psi_w) private(i,j,k,qc,dc,theta) &
   !$acc& reduction(+:nclip)
   !$omp OMPLOOP DEVICEPTR(delta_gpu,q_gpu) &
   !$omp& firstprivate(axis,sgn,b,ni,nj,nk,nv,nface_cells,scale,limit,mag_w,psi_w) private(i,j,k,qc,dc,theta) &
   !$omp& reduction(+:nclip)
   do c=1, nface_cells
      call seam_skin_cell(axis=axis, sgn=sgn, ni=ni, nj=nj, nk=nk, c=c, i=i, j=j, k=k)
      theta = 1._R8P
      if (limit) then
         !$acc loop seq
         do v=1, nv
            qc(v) = q_gpu(b,i,j,k,v)
            dc(v) = scale * delta_gpu(v,c)
         enddo
         theta = reflux_positivity_theta(nv=nv, q=qc(1:nv), d=dc(1:nv), mag_w=mag_w, psi_w=psi_w)
      endif
      if (theta < 1._R8P) then
         nclip = nclip + 1_I4P
         !$acc loop seq
         do v=1, nv
            q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + theta * (scale * delta_gpu(v,c))
         enddo
      else
         !$acc loop seq
         do v=1, nv
            q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + scale * delta_gpu(v,c)
         enddo
      endif
   enddo
   endsubroutine apply_reflux_face_dev

   subroutine copy_state_dev(ni, nj, nk, ngc, nv, blocks_number, src_gpu, dst_gpu)
   !< Copy blocks `1:blocks_number` of a device field, ghost cells included, into another (e.g. an idle RK stage buffer:
   !< the D_flux diagnostic evaluates the residual on a copy of the committed state, issue #68 A2).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                     !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                                  !< Variables number.
   integer(I4P), intent(in)    :: blocks_number                       !< Blocks to copy.
   real(R8P),    intent(in)    :: src_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Source field.
   real(R8P),    intent(inout) :: dst_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Destination field.
   integer(I4P)                :: b, i, j, k, v                       !< Counters.

   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(src_gpu,dst_gpu) &
   !$acc& firstprivate(ni,nj,nk,ngc,nv,blocks_number)
   !$omp OMPLOOP collapse(5) DEVICEPTR(src_gpu,dst_gpu) firstprivate(ni,nj,nk,ngc,nv,blocks_number)
   do v=1, nv
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, blocks_number
      dst_gpu(b,i,j,k,v) = src_gpu(b,i,j,k,v)
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine copy_state_dev

   subroutine compute_dflux_norms_dev(ni, nj, nk, ngc, blocks_number, band, ipsi, ibx, scale, dxyz_gpu, is_null,     &
                                      seam_gpu, flx_gpu, fly_gpu, flz_gpu, corr_gpu, q_gpu, dflux_max, dflux_l1,      &
                                      dflux_seam_max, h_dflux_max, b_max)
   !< Compute the D_flux norms (issue #68 A2), the device twin of the CPU `compute_dflux_norms` loop: D_flux of a cell is
   !< the divergence of the psi face fluxes minus the reflux correction `corr` of the coarse seam skin cells, over the
   !< psi flux per unit B_n (`scale`, c_h^2 or c_h); the norms are its maximum, its integral, its seam-local maximum (as
   !< `compute_divb_norms_dev`), the maximum of `h |D_flux|` (h the smallest active spacing of the block) and `max|B|`.
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                     !< Actual blocks number.
   integer(I4P), intent(in)  :: band                              !< Seam-local band width, in cells.
   integer(I4P), intent(in)  :: ipsi, ibx                         !< Indexes of psi and of B_x (B_y, B_z follow).
   real(R8P),    intent(in)  :: scale                             !< psi flux per unit B_n.
   real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Blocks space steps [nb, 3].
   logical,      intent(in)  :: is_null(3)                        !< Null directions.
   integer(I4P), intent(in)  :: seam_gpu(1:,1:)                   !< Seam faces flags [nb, 6] (0: not a seam face).
   real(R8P),    intent(in)  :: flx_gpu(1:,0:,1:,1:,1:)           !< X-face fluxes [nb, 0:ni, 1:nj, 1:nk, nv].
   real(R8P),    intent(in)  :: fly_gpu(1:,1:,0:,1:,1:)           !< Y-face fluxes [nb, 1:ni, 0:nj, 1:nk, nv].
   real(R8P),    intent(in)  :: flz_gpu(1:,1:,1:,0:,1:)           !< Z-face fluxes [nb, 1:ni, 1:nj, 0:nk, nv].
   real(R8P),    intent(in)  :: corr_gpu(1:,1:,1:,1:)             !< Reflux correction [blocks_number, ni, nj, nk].
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: dflux_max                         !< Maximum of |D_flux|.
   real(R8P),    intent(out) :: dflux_l1                          !< Sum of |D_flux| dV.
   real(R8P),    intent(out) :: dflux_seam_max                    !< Maximum of |D_flux| over the seam-local cells.
   real(R8P),    intent(out) :: h_dflux_max                       !< Maximum of h |D_flux|.
   real(R8P),    intent(out) :: b_max                             !< Maximum of |B|.
   real(R8P)                 :: wx, wy, wz                        !< Direction weights: 1 active, 0 null.
   real(R8P)                 :: div                               !< psi flux divergence of one cell.
   real(R8P)                 :: dflux                             !< |D_flux| of one cell.
   real(R8P)                 :: h                                 !< Smallest active spacing of the block.
   logical                   :: is_seam                           !< Seam-local cell.
   integer(I4P)              :: b, i, j, k                        !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   dflux_max = 0._R8P ; dflux_l1 = 0._R8P ; dflux_seam_max = 0._R8P ; h_dflux_max = 0._R8P ; b_max = 0._R8P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,seam_gpu,flx_gpu,fly_gpu,flz_gpu,corr_gpu,q_gpu)&
   !$acc& firstprivate(ni,nj,nk,blocks_number,band,ipsi,ibx,scale,wx,wy,wz) private(div,dflux,h,is_seam) &
   !$acc& reduction(max:dflux_max,dflux_seam_max,h_dflux_max,b_max) reduction(+:dflux_l1)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,seam_gpu,flx_gpu,fly_gpu,flz_gpu,corr_gpu,q_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,band,ipsi,ibx,scale,wx,wy,wz) private(div,dflux,h,is_seam) &
   !$omp& reduction(max:dflux_max,dflux_seam_max,h_dflux_max,b_max) reduction(+:dflux_l1)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      div = wx * (flx_gpu(b,i,j,k,ipsi) - flx_gpu(b,i-1,j,k,ipsi)) / dxyz_gpu(b,1) + &
            wy * (fly_gpu(b,i,j,k,ipsi) - fly_gpu(b,i,j-1,k,ipsi)) / dxyz_gpu(b,2) + &
            wz * (flz_gpu(b,i,j,k,ipsi) - flz_gpu(b,i,j,k-1,ipsi)) / dxyz_gpu(b,3)
      dflux = abs(div - corr_gpu(b,i,j,k)) / scale
      h = huge(1._R8P)
      if (wx > 0._R8P) h = min(h, dxyz_gpu(b,1))
      if (wy > 0._R8P) h = min(h, dxyz_gpu(b,2))
      if (wz > 0._R8P) h = min(h, dxyz_gpu(b,3))
      dflux_max = max(dflux_max, dflux)
      dflux_l1  = dflux_l1 + dflux * dxyz_gpu(b,1) * dxyz_gpu(b,2) * dxyz_gpu(b,3)
      is_seam = (seam_gpu(b,1) /= 0_I4P .and. i <= band) .or. (seam_gpu(b,2) /= 0_I4P .and. i > ni - band) .or. &
                (seam_gpu(b,3) /= 0_I4P .and. j <= band) .or. (seam_gpu(b,4) /= 0_I4P .and. j > nj - band) .or. &
                (seam_gpu(b,5) /= 0_I4P .and. k <= band) .or. (seam_gpu(b,6) /= 0_I4P .and. k > nk - band)
      if (is_seam) dflux_seam_max = max(dflux_seam_max, dflux)
      h_dflux_max = max(h_dflux_max, h * dflux)
      b_max = max(b_max, sqrt(q_gpu(b,i,j,k,ibx)**2 + q_gpu(b,i,j,k,ibx+1)**2 + q_gpu(b,i,j,k,ibx+2)**2))
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_dflux_norms_dev

   subroutine compute_flux_difference_dev(nv, ni, nj, nk, ngc, blocks_number, is_null, freeze, dxyz_gpu, flx_f_gpu,    &
                                          fly_f_gpu, flz_f_gpu, dq_gpu)
   !< Compute the residuals from the face fluxes, `dq = -sum_d (F_{d,i+1/2} - F_{d,i-1/2}) / dx_d`.
   !<
   !< A null direction weighs zero, and the residual of the variable it freezes, `freeze(d)` (0: none), is zero: the
   !< normal momentum for Euler (CHASE semantics, issue #35, section 3.4), none for MHD (issue #41, M2-P3).
   integer(I4P), intent(in)    :: nv                                 !< Conservative variables number.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   logical,      intent(in)    :: is_null(3)                         !< Null directions.
   integer(I4P), intent(in)    :: freeze(3)                          !< Variable frozen by each null direction (0: none).
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                    !< Blocks space steps [nb, 3].
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)          !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)          !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)          !< Z-face fluxes.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   real(R8P)                   :: wx, wy, wz                         !< Direction weights: 1 active, 0 null.
   integer(I4P)                :: fx, fy, fz                         !< Frozen variables, scalar copies.
   integer(I4P)                :: b, i, j, k, v                      !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1)) ; fx = freeze(1)
   wy = merge(0._R8P, 1._R8P, is_null(2)) ; fy = freeze(2)
   wz = merge(0._R8P, 1._R8P, is_null(3)) ; fz = freeze(3)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,dq_gpu) &
   !$acc& firstprivate(nv,ni,nj,nk,blocks_number,wx,wy,wz,fx,fy,fz)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,dq_gpu) &
   !$omp& firstprivate(nv,ni,nj,nk,blocks_number,wx,wy,wz,fx,fy,fz)
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
      if (fx > 0_I4P) dq_gpu(b,i,j,k,fx) = 0._R8P
      if (fy > 0_I4P) dq_gpu(b,i,j,k,fy) = 0._R8P
      if (fz > 0_I4P) dq_gpu(b,i,j,k,fz) = 0._R8P
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_flux_difference_dev

   subroutine add_gravity_sources_dev(ni, nj, nk, ngc, blocks_number, g, q_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, &
                                      gc_gpu, dq_gpu)
   !< Add the external gravity sources (issue #73, M7-G1; device twin of the CPU `add_gravity_sources`): `rho g` to the
   !< momentum, `1/2 sum_d g_d (F^rho_{d,i+1/2} + F^rho_{d,i-1/2})` plus its fourth-order correction (M7-G2a) to the
   !< energy.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   real(R8P),    intent(in)    :: g(3)                               !< Gravity.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)          !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)          !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)          !< Z-face fluxes.
   real(R8P),    intent(in)    :: gc_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fourth-order correction of the energy source.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Residuals.
   real(R8P)                   :: gx, gy, gz                         !< Gravity, scalar copies.
   integer(I4P)                :: b, i, j, k                         !< Counters.

   gx = g(1) ; gy = g(2) ; gz = g(3)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,gc_gpu,dq_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,gx,gy,gz)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,gc_gpu,dq_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,gx,gy,gz)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      dq_gpu(b,i,j,k,IQ_RU) = dq_gpu(b,i,j,k,IQ_RU) + q_gpu(b,i,j,k,IQ_R) * gx
      dq_gpu(b,i,j,k,IQ_RV) = dq_gpu(b,i,j,k,IQ_RV) + q_gpu(b,i,j,k,IQ_R) * gy
      dq_gpu(b,i,j,k,IQ_RW) = dq_gpu(b,i,j,k,IQ_RW) + q_gpu(b,i,j,k,IQ_R) * gz
      dq_gpu(b,i,j,k,IQ_RE) = dq_gpu(b,i,j,k,IQ_RE) +                                      &
                              0.5_R8P * (gx * (flx_f_gpu(b,i,j,k,IQ_R) + flx_f_gpu(b,i-1,j,k,IQ_R)) + &
                                         gy * (fly_f_gpu(b,i,j,k,IQ_R) + fly_f_gpu(b,i,j-1,k,IQ_R)) + &
                                         gz * (flz_f_gpu(b,i,j,k,IQ_R) + flz_f_gpu(b,i,j,k-1,IQ_R))) + &
                              gc_gpu(b,i,j,k,1)
   enddo
   enddo
   enddo
   enddo
   endsubroutine add_gravity_sources_dev

   subroutine add_potential_sources_dev(ni, nj, nk, ngc, blocks_number, is_null, dxyz_gpu, phi_gpu, q_gpu, flx_f_gpu, &
                                        fly_f_gpu, flz_f_gpu, dq_gpu)
   !< Add the gravity of the cell potential (issue #81, M8-P1; device twin of the CPU `add_potential_sources`): face
   !< gravity `-(Phi_{i+1} - Phi_i) / dx_d`, cell gravity the mean of the two faces; `rho` times the cell gravity to the
   !< momentum, `1/2 sum_d (g_{d,i+1/2} F^rho_{d,i+1/2} + g_{d,i-1/2} F^rho_{d,i-1/2})` to the energy.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                     !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                       !< Actual blocks number.
   logical,      intent(in)    :: is_null(3)                          !< Null directions.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                     !< Blocks space steps [nb, 3].
   real(R8P),    intent(in)    :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Gravitational potential, ghosts filled.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Conservative variables.
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)           !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)           !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)           !< Z-face fluxes.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Residuals.
   logical                     :: ax, ay, az                          !< Active directions, scalar copies.
   real(R8P)                   :: gr, gl                              !< Private gravity of the right and left faces.
   real(R8P)                   :: gmx, gmy, gmz                       !< Private cell gravity.
   real(R8P)                   :: se                                  !< Private energy source.
   integer(I4P)                :: b, i, j, k                          !< Counters.

   ax = .not.is_null(1) ; ay = .not.is_null(2) ; az = .not.is_null(3)
   !$acc parallel loop independent gang vector collapse(4) &
   !$acc& DEVICEVAR(dxyz_gpu,phi_gpu,q_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,dq_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,ax,ay,az) private(gr,gl,gmx,gmy,gmz,se)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,phi_gpu,q_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,dq_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,ax,ay,az) private(gr,gl,gmx,gmy,gmz,se)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      gmx = 0._R8P ; gmy = 0._R8P ; gmz = 0._R8P
      se = 0._R8P
      if (ax) then
         gr = -(phi_gpu(b,i+1,j,k,1) - phi_gpu(b,i,j,k,1)) / dxyz_gpu(b,1)
         gl = -(phi_gpu(b,i,j,k,1) - phi_gpu(b,i-1,j,k,1)) / dxyz_gpu(b,1)
         gmx = 0.5_R8P * (gr + gl)
         se = se + 0.5_R8P * (gr * flx_f_gpu(b,i,j,k,IQ_R) + gl * flx_f_gpu(b,i-1,j,k,IQ_R))
      endif
      if (ay) then
         gr = -(phi_gpu(b,i,j+1,k,1) - phi_gpu(b,i,j,k,1)) / dxyz_gpu(b,2)
         gl = -(phi_gpu(b,i,j,k,1) - phi_gpu(b,i,j-1,k,1)) / dxyz_gpu(b,2)
         gmy = 0.5_R8P * (gr + gl)
         se = se + 0.5_R8P * (gr * fly_f_gpu(b,i,j,k,IQ_R) + gl * fly_f_gpu(b,i,j-1,k,IQ_R))
      endif
      if (az) then
         gr = -(phi_gpu(b,i,j,k+1,1) - phi_gpu(b,i,j,k,1)) / dxyz_gpu(b,3)
         gl = -(phi_gpu(b,i,j,k,1) - phi_gpu(b,i,j,k-1,1)) / dxyz_gpu(b,3)
         gmz = 0.5_R8P * (gr + gl)
         se = se + 0.5_R8P * (gr * flz_f_gpu(b,i,j,k,IQ_R) + gl * flz_f_gpu(b,i,j,k-1,IQ_R))
      endif
      dq_gpu(b,i,j,k,IQ_RU) = dq_gpu(b,i,j,k,IQ_RU) + q_gpu(b,i,j,k,IQ_R) * gmx
      dq_gpu(b,i,j,k,IQ_RV) = dq_gpu(b,i,j,k,IQ_RV) + q_gpu(b,i,j,k,IQ_R) * gmy
      dq_gpu(b,i,j,k,IQ_RW) = dq_gpu(b,i,j,k,IQ_RW) + q_gpu(b,i,j,k,IQ_R) * gmz
      dq_gpu(b,i,j,k,IQ_RE) = dq_gpu(b,i,j,k,IQ_RE) + se
   enddo
   enddo
   enddo
   enddo
   endsubroutine add_potential_sources_dev

   subroutine subtract_equilibrium_residual_dev(nv, ni, nj, nk, ngc, blocks_number, dq_e_gpu, dq_gpu)
   !< Subtract the residual of the hydrostatic equilibrium from the residual of the stage (issue #73, M7-G2b; device
   !< twin of the CPU `subtract_equilibrium_residual`).
   integer(I4P), intent(in)    :: nv                                   !< Variables number.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                      !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                        !< Actual blocks number.
   real(R8P),    intent(in)    :: dq_e_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Equilibrium residual.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Residuals.
   integer(I4P)                :: b, i, j, k, v                        !< Counters.

   !$acc parallel loop independent gang vector collapse(5) DEVICEVAR(dq_e_gpu,dq_gpu) &
   !$acc& firstprivate(nv,ni,nj,nk,blocks_number)
   !$omp OMPLOOP collapse(5) DEVICEPTR(dq_e_gpu,dq_gpu) firstprivate(nv,ni,nj,nk,blocks_number)
   do v=1, nv
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) - dq_e_gpu(b,i,j,k,v)
   enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine subtract_equilibrium_residual_dev

   subroutine compute_gravity_correction_faces_dev(ni, nj, nk, ngc, blocks_number, g, flx_f_gpu, fly_f_gpu, flz_f_gpu, &
                                                   gc_gpu)
   !< Set the fourth-order correction of the gravity energy source from the faces the block holds (issue #73, M7-G2a;
   !< device twin of the CPU `compute_gravity_correction_faces`): `sum_d g_d (-F_{i-3/2} + F_{i-1/2} + F_{i+1/2} -
   !< F_{i+3/2}) / 24`, the faces beyond the block left to `add_gravity_edge_correction_dev`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   real(R8P),    intent(in)    :: g(3)                               !< Gravity.
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)          !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)          !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)          !< Z-face fluxes.
   real(R8P),    intent(inout) :: gc_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Gravity energy correction [..., 1].
   real(R8P)                   :: gx, gy, gz                         !< Gravity, scalar copies.
   real(R8P)                   :: c, t                               !< Correction, one direction's term.
   integer(I4P)                :: b, i, j, k                         !< Counters.

   gx = g(1) ; gy = g(2) ; gz = g(3)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(flx_f_gpu,fly_f_gpu,flz_f_gpu,gc_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,gx,gy,gz) private(c,t)
   !$omp OMPLOOP collapse(4) DEVICEPTR(flx_f_gpu,fly_f_gpu,flz_f_gpu,gc_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,gx,gy,gz) private(c,t)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      c = 0._R8P
      if (gx /= 0._R8P) then
         t = flx_f_gpu(b,i-1,j,k,IQ_R) + flx_f_gpu(b,i,j,k,IQ_R)
         if (i > 1)  t = t - flx_f_gpu(b,i-2,j,k,IQ_R)
         if (i < ni) t = t - flx_f_gpu(b,i+1,j,k,IQ_R)
         c = c + gx * t
      endif
      if (gy /= 0._R8P) then
         t = fly_f_gpu(b,i,j-1,k,IQ_R) + fly_f_gpu(b,i,j,k,IQ_R)
         if (j > 1)  t = t - fly_f_gpu(b,i,j-2,k,IQ_R)
         if (j < nj) t = t - fly_f_gpu(b,i,j+1,k,IQ_R)
         c = c + gy * t
      endif
      if (gz /= 0._R8P) then
         t = flz_f_gpu(b,i,j,k-1,IQ_R) + flz_f_gpu(b,i,j,k,IQ_R)
         if (k > 1)  t = t - flz_f_gpu(b,i,j,k-2,IQ_R)
         if (k < nk) t = t - flz_f_gpu(b,i,j,k+1,IQ_R)
         c = c + gz * t
      endif
      gc_gpu(b,i,j,k,1) = c / 24._R8P
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_gravity_correction_faces_dev

   subroutine add_gravity_limited_correction_dev(ni, nj, nk, ngc, blocks_number, lam_gpu, dq_gpu)
   !< Remove `(1 - Lambda) L` from the gravity energy source (issue #73, M7-G1, option B; device twin of the CPU
   !< `add_gravity_limited_correction`).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                     !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                       !< Actual blocks number.
   real(R8P),    intent(in)    :: lam_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Cell factors (1) and gravity excess (2).
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Residuals.
   integer(I4P)                :: b, i, j, k                          !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(lam_gpu,dq_gpu) firstprivate(ni,nj,nk,blocks_number)
   !$omp OMPLOOP collapse(4) DEVICEPTR(lam_gpu,dq_gpu) firstprivate(ni,nj,nk,blocks_number)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      dq_gpu(b,i,j,k,IQ_RE) = dq_gpu(b,i,j,k,IQ_RE) - (1._R8P - lam_gpu(b,i,j,k,1)) * lam_gpu(b,i,j,k,2)
   enddo
   enddo
   enddo
   enddo
   endsubroutine add_gravity_limited_correction_dev

   subroutine compute_flux_difference_ib_dev(nv, ni, nj, nk, ngc, blocks_number, is_null, freeze, dxyz_gpu, flx_f_gpu, &
                                             fly_f_gpu, flz_f_gpu, phi_gpu, dq_gpu)
   !< Compute the residuals from the face fluxes with immersed solids: the spacing of a fluid cell is cut by the solid
   !< surface (`ib_cut_spacing`, CHASE semantics, issue #35 D-9); device twin of the CPU flux difference with `phi`.
   integer(I4P), intent(in)    :: nv                                  !< Conservative variables number.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                     !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                       !< Actual blocks number.
   logical,      intent(in)    :: is_null(3)                          !< Null directions.
   integer(I4P), intent(in)    :: freeze(3)                           !< Variable frozen by each null direction (0: none).
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                     !< Blocks space steps [nb, 3].
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:)           !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:)           !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:)           !< Z-face fluxes.
   real(R8P),    intent(in)    :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function [nb, i, j, k, solids+1].
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Residuals.
   real(R8P)                   :: wx, wy, wz                          !< Direction weights: 1 active, 0 null.
   real(R8P)                   :: dx, dy, dz                          !< Cell spacings.
   integer(I4P)                :: fx, fy, fz                          !< Frozen variables, scalar copies.
   integer(I4P)                :: ns                                  !< All-solids summary slot of phi.
   integer(I4P)                :: b, i, j, k, v                       !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1)) ; fx = freeze(1)
   wy = merge(0._R8P, 1._R8P, is_null(2)) ; fy = freeze(2)
   wz = merge(0._R8P, 1._R8P, is_null(3)) ; fz = freeze(3)
   ns = ubound(phi_gpu, dim=5)
   !$acc parallel loop independent gang vector collapse(4)                                    &
   !$acc& DEVICEVAR(dxyz_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,phi_gpu,dq_gpu)                     &
   !$acc& firstprivate(nv,ni,nj,nk,blocks_number,wx,wy,wz,fx,fy,fz,ns) private(dx,dy,dz)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,phi_gpu,dq_gpu) &
   !$omp& firstprivate(nv,ni,nj,nk,blocks_number,wx,wy,wz,fx,fy,fz,ns) private(dx,dy,dz)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      dx = ib_cut_spacing(phi_c=phi_gpu(b,i,j,k,ns), phi_m=phi_gpu(b,i-1,j,k,ns), phi_p=phi_gpu(b,i+1,j,k,ns), &
                          ds=dxyz_gpu(b,1))
      dy = ib_cut_spacing(phi_c=phi_gpu(b,i,j,k,ns), phi_m=phi_gpu(b,i,j-1,k,ns), phi_p=phi_gpu(b,i,j+1,k,ns), &
                          ds=dxyz_gpu(b,2))
      dz = ib_cut_spacing(phi_c=phi_gpu(b,i,j,k,ns), phi_m=phi_gpu(b,i,j,k-1,ns), phi_p=phi_gpu(b,i,j,k+1,ns), &
                          ds=dxyz_gpu(b,3))
      !$acc loop seq
      do v=1, nv
         dq_gpu(b,i,j,k,v) = -(wx * (flx_f_gpu(b,i,j,k,v) - flx_f_gpu(b,i-1,j,k,v)) / dx + &
                               wy * (fly_f_gpu(b,i,j,k,v) - fly_f_gpu(b,i,j-1,k,v)) / dy + &
                               wz * (flz_f_gpu(b,i,j,k,v) - flz_f_gpu(b,i,j,k-1,v)) / dz)
      enddo
      if (fx > 0_I4P) dq_gpu(b,i,j,k,fx) = 0._R8P
      if (fy > 0_I4P) dq_gpu(b,i,j,k,fy) = 0._R8P
      if (fz > 0_I4P) dq_gpu(b,i,j,k,fz) = 0._R8P
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_flux_difference_ib_dev

   subroutine poison_ghosts_dev(ni, nj, nk, ngc, nv, blocks_number, q_gpu)
   !< Set every ghost cell to a quiet NaN (verification instrument, `[diagnostics] ghost_poison`, issue #65).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                                !< Variables number.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P)                   :: nan                               !< Quiet NaN.
   integer(I4P)                :: b, i, j, k, v                     !< Counters.

   nan = transfer(int(z'7FF8000000000000', I8P), nan)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu) firstprivate(ni,nj,nk,ngc,nv,nan)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu) firstprivate(ni,nj,nk,ngc,nv,nan)
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, blocks_number
      if (i < 1 .or. i > ni .or. j < 1 .or. j > nj .or. k < 1 .or. k > nk) then
         !$acc loop seq
         do v=1, nv
            q_gpu(b,i,j,k,v) = nan
         enddo
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine poison_ghosts_dev

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

   subroutine pack_seam_skins_dev(n_faces, max_cells, faces_gpu, ni, nj, nk, nv, flx_f_gpu, fly_f_gpu, flz_f_gpu, &
                                  skins_gpu)
   !< Pack the face fluxes of every seam register face `f` of the rank in one launch (issue #90, F-B): face
   !< `faces_gpu(:,f) = [block, fec, first column - 1]` goes to `skins_gpu(:, first column:)`, its cells over the two
   !< tangential axes, inner fastest (x faces: j, k; y faces: i, k; z faces: i, j), the flux register order
   !< (`accumulate_seam_skin`).
   integer(I4P), intent(in)    :: n_faces                   !< Register faces number.
   integer(I4P), intent(in)    :: max_cells                 !< Most skin cells of a face, max(nj nk, ni nk, ni nj).
   integer(I4P), intent(in)    :: faces_gpu(1:,1:)          !< Faces [3, n_faces]: block, fec, first column - 1.
   integer(I4P), intent(in)    :: ni, nj, nk                !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                        !< Variables number.
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:) !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:) !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:) !< Z-face fluxes.
   real(R8P),    intent(inout) :: skins_gpu(1:,1:)          !< Packed skins (nv, columns).
   integer(I4P)                :: f, cell, inner, outer, v  !< Counters.
   integer(I4P)                :: b, fec, col               !< Face block, fec, skin column.
   integer(I4P)                :: inner_n, outer_n          !< Tangential cells numbers.
   integer(I4P)                :: n                         !< Face normal index (0 or n).

   !$acc parallel loop independent gang vector collapse(2) &
   !$acc& DEVICEVAR(faces_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,skins_gpu) &
   !$acc& firstprivate(n_faces,max_cells,ni,nj,nk,nv) private(b,fec,col,inner,outer,inner_n,outer_n,n)
   !$omp OMPLOOP collapse(2) DEVICEPTR(faces_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,skins_gpu) &
   !$omp& firstprivate(n_faces,max_cells,ni,nj,nk,nv) private(b,fec,col,inner,outer,inner_n,outer_n,n)
   do f=1, n_faces
   do cell=1, max_cells
      b = faces_gpu(1,f) ; fec = faces_gpu(2,f)
      select case(fec)
      case(1_I4P, 2_I4P)
         inner_n = nj ; outer_n = nk ; n = merge(0_I4P, ni, fec == 1_I4P)
      case(3_I4P, 4_I4P)
         inner_n = ni ; outer_n = nk ; n = merge(0_I4P, nj, fec == 3_I4P)
      case default
         inner_n = ni ; outer_n = nj ; n = merge(0_I4P, nk, fec == 5_I4P)
      endselect
      if (cell <= inner_n * outer_n) then
         inner = mod(cell - 1_I4P, inner_n) + 1_I4P ; outer = (cell - 1_I4P) / inner_n + 1_I4P
         col = faces_gpu(3,f) + cell
         !$acc loop seq
         do v=1, nv
            select case(fec)
            case(1_I4P, 2_I4P)
               skins_gpu(v,col) = flx_f_gpu(b,n,inner,outer,v)
            case(3_I4P, 4_I4P)
               skins_gpu(v,col) = fly_f_gpu(b,inner,n,outer,v)
            case default
               skins_gpu(v,col) = flz_f_gpu(b,inner,outer,n,v)
            endselect
         enddo
      endif
   enddo
   enddo
   endsubroutine pack_seam_skins_dev

   subroutine pack_density_dev(ni, nj, nk, ngc, blocks_number, q_gpu, buf_gpu)
   !< Pack the interior density of every block into `buf_gpu`, block-major then `k`, `j`, `i` fastest: the order of the
   !< library FFT Poisson solver's `solve_packed` (issue #81, M8-P3), which gathers it on the host.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(inout) :: buf_gpu(1:)                       !< Packed density [blocks_number ni nj nk].
   integer(I4P)                :: b, i, j, k                        !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu,buf_gpu) firstprivate(ni,nj,nk,blocks_number)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu,buf_gpu) firstprivate(ni,nj,nk,blocks_number)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      buf_gpu((((b-1)*nk + (k-1))*nj + (j-1))*ni + i) = q_gpu(b,i,j,k,IQ_R)
   enddo
   enddo
   enddo
   enddo
   endsubroutine pack_density_dev

   subroutine unpack_potential_dev(ni, nj, nk, ngc, blocks_number, buf_gpu, phi_gpu)
   !< Unpack the potential of the active blocks, every ghost layer included, from `buf_gpu` (the host layout
   !< `(i, j, k, b)` flattened, `i` fastest) into `phi_gpu(b, i, j, k, 1)` (issue #81, M8-P3). Only the active blocks
   !< cross the bus: the host and device arrays are sized to the block capacity, not to the blocks in use.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                     !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                       !< Actual blocks number.
   real(R8P),    intent(in)    :: buf_gpu(1:)                         !< Packed potential [(ni+2ngc)(nj+2ngc)(nk+2ngc) nb].
   real(R8P),    intent(inout) :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Potential (b, i, j, k, 1).
   integer(I4P)                :: mi, mj, mk                          !< Cells per block per direction, ghosts included.
   integer(I4P)                :: b, i, j, k                          !< Counters.

   mi = ni + 2 * ngc ; mj = nj + 2 * ngc ; mk = nk + 2 * ngc
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(buf_gpu,phi_gpu) firstprivate(ngc,mi,mj,mk,blocks_number)
   !$omp OMPLOOP collapse(4) DEVICEPTR(buf_gpu,phi_gpu) firstprivate(ngc,mi,mj,mk,blocks_number)
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, blocks_number
      phi_gpu(b,i,j,k,1) = buf_gpu((((b-1)*mk + (k+ngc-1))*mj + (j+ngc-1))*mi + i + ngc)
   enddo
   enddo
   enddo
   enddo
   endsubroutine unpack_potential_dev

   subroutine gather_seam_cells_dev(n, nv, ngc, cells_gpu, a_gpu, out_gpu)
   !< Gather `out_gpu(v, m) = a_gpu(b, i, j, k, v)` at the cells `cells_gpu(:, m) = (i, j, k, b)` of a ghosted field
   !< (issue #50: the seam synchronisation of the positivity limiter works on the host, over the seam cells only).
   integer(I4P), intent(in)    :: n                                 !< Cells number.
   integer(I4P), intent(in)    :: nv                                !< Components gathered (the first `nv`).
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)    :: cells_gpu(1:,1:)                  !< Cells (i, j, k, b) [4, n].
   real(R8P),    intent(in)    :: a_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field.
   real(R8P),    intent(inout) :: out_gpu(1:,1:)                    !< Gathered values [nv, n].
   integer(I4P)                :: m, v                              !< Counters.

   !$acc parallel loop independent gang vector DEVICEVAR(cells_gpu,a_gpu,out_gpu) firstprivate(n,nv)
   !$omp OMPLOOP DEVICEPTR(cells_gpu,a_gpu,out_gpu) firstprivate(n,nv)
   do m=1, n
      !$acc loop seq
      do v=1, nv
         out_gpu(v,m) = a_gpu(cells_gpu(4,m),cells_gpu(1,m),cells_gpu(2,m),cells_gpu(3,m),v)
      enddo
   enddo
   endsubroutine gather_seam_cells_dev

   subroutine gather_seam_stencils_dev(n, nv, ngc, s, cells_gpu, q_gpu, out_gpu)
   !< Gather the three axis stencils of half width `s` around the cells `cells_gpu(:, m) = (i, j, k, b)`:
   !< `out_gpu(v, l + s + 1, d, m)` is the state `l` cells from the cell along `d` (issue #50).
   integer(I4P), intent(in)    :: n                                 !< Cells number.
   integer(I4P), intent(in)    :: nv                                !< Variables number.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)    :: s                                 !< Stencil half width (s <= ngc).
   integer(I4P), intent(in)    :: cells_gpu(1:,1:)                  !< Cells (i, j, k, b) [4, n].
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(inout) :: out_gpu(1:,1:,1:,1:)              !< Stencils [nv, 2 s + 1, 3, n].
   integer(I4P)                :: m, l, v, b, i, j, k               !< Counters.

   !$acc parallel loop independent gang vector DEVICEVAR(cells_gpu,q_gpu,out_gpu) firstprivate(n,nv,s) private(b,i,j,k)
   !$omp OMPLOOP DEVICEPTR(cells_gpu,q_gpu,out_gpu) firstprivate(n,nv,s) private(b,i,j,k)
   do m=1, n
      i = cells_gpu(1,m) ; j = cells_gpu(2,m) ; k = cells_gpu(3,m) ; b = cells_gpu(4,m)
      !$acc loop seq
      do l=-s, s
         !$acc loop seq
         do v=1, nv
            out_gpu(v,l+s+1,1,m) = q_gpu(b,i+l,j,k,v)
            out_gpu(v,l+s+1,2,m) = q_gpu(b,i,j+l,k,v)
            out_gpu(v,l+s+1,3,m) = q_gpu(b,i,j,k+l,v)
         enddo
      enddo
   enddo
   endsubroutine gather_seam_stencils_dev

   subroutine gather_seam_faces_dev(n, nv, faces_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu, out_gpu)
   !< Gather the face fluxes at the faces `faces_gpu(:, m) = (axis, i, j, k, b)`, `(i, j, k)` the face index of the
   !< flux array of `axis` (issue #50).
   integer(I4P), intent(in)    :: n                         !< Faces number.
   integer(I4P), intent(in)    :: nv                        !< Variables number.
   integer(I4P), intent(in)    :: faces_gpu(1:,1:)          !< Faces (axis, i, j, k, b) [5, n].
   real(R8P),    intent(in)    :: flx_f_gpu(1:,0:,1:,1:,1:) !< X-face fluxes.
   real(R8P),    intent(in)    :: fly_f_gpu(1:,1:,0:,1:,1:) !< Y-face fluxes.
   real(R8P),    intent(in)    :: flz_f_gpu(1:,1:,1:,0:,1:) !< Z-face fluxes.
   real(R8P),    intent(inout) :: out_gpu(1:,1:)            !< Gathered fluxes [nv, n].
   integer(I4P)                :: m, v                      !< Counters.

   !$acc parallel loop independent gang vector DEVICEVAR(faces_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,out_gpu) firstprivate(n,nv)
   !$omp OMPLOOP DEVICEPTR(faces_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu,out_gpu) firstprivate(n,nv)
   do m=1, n
      !$acc loop seq
      do v=1, nv
         select case(faces_gpu(1,m))
         case(1)
            out_gpu(v,m) = flx_f_gpu(faces_gpu(5,m),faces_gpu(2,m),faces_gpu(3,m),faces_gpu(4,m),v)
         case(2)
            out_gpu(v,m) = fly_f_gpu(faces_gpu(5,m),faces_gpu(2,m),faces_gpu(3,m),faces_gpu(4,m),v)
         case default
            out_gpu(v,m) = flz_f_gpu(faces_gpu(5,m),faces_gpu(2,m),faces_gpu(3,m),faces_gpu(4,m),v)
         endselect
      enddo
   enddo
   endsubroutine gather_seam_faces_dev

   subroutine scatter_seam_faces_dev(n, nv, faces_gpu, in_gpu, flx_f_gpu, fly_f_gpu, flz_f_gpu)
   !< Set the face fluxes at the faces `faces_gpu(:, m) = (axis, i, j, k, b)` to `in_gpu(:, m)` (issue #50).
   integer(I4P), intent(in)    :: n                         !< Faces number.
   integer(I4P), intent(in)    :: nv                        !< Variables number.
   integer(I4P), intent(in)    :: faces_gpu(1:,1:)          !< Faces (axis, i, j, k, b) [5, n].
   real(R8P),    intent(in)    :: in_gpu(1:,1:)             !< Fluxes [nv, n].
   real(R8P),    intent(inout) :: flx_f_gpu(1:,0:,1:,1:,1:) !< X-face fluxes.
   real(R8P),    intent(inout) :: fly_f_gpu(1:,1:,0:,1:,1:) !< Y-face fluxes.
   real(R8P),    intent(inout) :: flz_f_gpu(1:,1:,1:,0:,1:) !< Z-face fluxes.
   integer(I4P)                :: m, v                      !< Counters.

   !$acc parallel loop independent gang vector DEVICEVAR(faces_gpu,in_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu) firstprivate(n,nv)
   !$omp OMPLOOP DEVICEPTR(faces_gpu,in_gpu,flx_f_gpu,fly_f_gpu,flz_f_gpu) firstprivate(n,nv)
   do m=1, n
      !$acc loop seq
      do v=1, nv
         select case(faces_gpu(1,m))
         case(1)
            flx_f_gpu(faces_gpu(5,m),faces_gpu(2,m),faces_gpu(3,m),faces_gpu(4,m),v) = in_gpu(v,m)
         case(2)
            fly_f_gpu(faces_gpu(5,m),faces_gpu(2,m),faces_gpu(3,m),faces_gpu(4,m),v) = in_gpu(v,m)
         case default
            flz_f_gpu(faces_gpu(5,m),faces_gpu(2,m),faces_gpu(3,m),faces_gpu(4,m),v) = in_gpu(v,m)
         endselect
      enddo
   enddo
   endsubroutine scatter_seam_faces_dev

   subroutine scatter_seam_cells_dev(n, nv, ngc, cells_gpu, in_gpu, a_gpu)
   !< Set `a_gpu(b, i, j, k, v) = in_gpu(v, m)` at the cells `cells_gpu(:, m) = (i, j, k, b)` (issue #50).
   integer(I4P), intent(in)    :: n                                 !< Cells number.
   integer(I4P), intent(in)    :: nv                                !< Components set (the first `nv`).
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)    :: cells_gpu(1:,1:)                  !< Cells (i, j, k, b) [4, n].
   real(R8P),    intent(in)    :: in_gpu(1:,1:)                     !< Values [nv, n].
   real(R8P),    intent(inout) :: a_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field.
   integer(I4P)                :: m, v                              !< Counters.

   !$acc parallel loop independent gang vector DEVICEVAR(cells_gpu,in_gpu,a_gpu) firstprivate(n,nv)
   !$omp OMPLOOP DEVICEPTR(cells_gpu,in_gpu,a_gpu) firstprivate(n,nv)
   do m=1, n
      !$acc loop seq
      do v=1, nv
         a_gpu(cells_gpu(4,m),cells_gpu(1,m),cells_gpu(2,m),cells_gpu(3,m),v) = in_gpu(v,m)
      enddo
   enddo
   endsubroutine scatter_seam_cells_dev

   subroutine set_boundary_conditions_dev(ni, nj, nk, ngc, nv, nv_mhd, nv_entropy, crown, pass, face_kind, gamma, R,    &
                                          psi_energy, &
                                          wall_velocity, wall_temperature, g, has_phi, phi_gpu,             &
                                          local_map_bc_crown_gpu, row_first, row_last, bc_rows_gpu,         &
                                          q_inflow_gpu, wall_sign_gpu, dxyz_gpu, q_gpu)
   !< Set boundary conditions on one crown, for the rows beyond `pass` realm faces (see the CPU
   !< `set_boundary_conditions`): pass 1 the face rows (each face's kind, inward along its normal), pass 2 and 3 the
   !< rows beyond two or three realm faces (the kind of one physical face among them, `realm_edge_face`, on the donor
   !< `realm_edge_donor`, issue #65 P0).
   !<
   !< The caller runs the passes in order, and every donor is interior, filled by the exchange or the seams, or by an
   !< earlier pass: a face row reads the mirrored (walls) or the first interior (extrapolation) cell along its normal,
   !< never another ghost of its pass, so the rows of one launch are independent.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: nv                                !< Conservative variables number.
   integer(I4P), intent(in)    :: nv_mhd                            !< Last slot of the MHD state (5 hydrodynamic, nv MHD).
   integer(I4P), intent(in)    :: nv_entropy                        !< Last entropy slot (6 with dual energy, else nv_mhd).
   integer(I4P), intent(in)    :: crown                             !< Crown counter.
   integer(I4P), intent(in)    :: row_first, row_last               !< The rows of `pass` in `bc_rows_gpu(:,crown)`.
   integer(I4P), intent(in)    :: bc_rows_gpu(:,:)                  !< Crown map rows grouped by pass, per crown.
   integer(I4P), intent(in)    :: pass                              !< Realm faces the rows of this launch lie beyond.
   integer(I4P), intent(in)    :: face_kind(6)                      !< Kind of each realm face, BC_SEAM for a seam.
   real(R8P),    intent(in)    :: gamma, R                          !< Specific heats ratio, gas constant.
   real(R8P),    intent(in)    :: psi_energy                        !< 1 if the total energy holds psi^2 / 2 (EGLM).
   real(R8P),    intent(in)    :: wall_velocity(3,6)                !< No-slip walls: wall velocity of each face.
   real(R8P),    intent(in)    :: wall_temperature(6)               !< Isothermal walls: wall temperature of each face.
   real(R8P),    intent(in)    :: g(3)                              !< Hydrostatic faces: gravity (issue #73).
   logical,      intent(in)    :: has_phi                           !< Hydrostatic faces: a cell potential (issue #81).
   real(R8P),    intent(in)    :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Gravitational potential (read with has_phi).
   integer(I8P), intent(in)    :: local_map_bc_crown_gpu(:,:,:)     !< Boundary crown map (row, field, crown).
   real(R8P),    intent(in)    :: q_inflow_gpu(1:,1:)               !< Conservative inflow state of each face [nv, 6].
   real(R8P),    intent(in)    :: wall_sign_gpu(1:,1:)              !< Wall mirror sign per variable and direction [nv, 3].
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                   !< Cell spacing of each block [nb, 3].
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   integer(I4P)                :: b, c, row, i, j, k, v             !< Counters.
   integer(I4P)                :: idelta, jdelta, kdelta            !< IJK inward step.
   integer(I4P)                :: bc_type                           !< Boundary condition type.
   integer(I4P)                :: fec                               !< Boundary fec (1 to 26).
   integer(I4P)                :: face                              !< Boundary face (1 to 6).
   integer(I4P)                :: iref, jref, kref                  !< Donor indexes.
   real(R8P)                   :: qd(9), qw(9)                      !< No-slip walls: donor and ghost states.
   integer(I4P)                :: i1, j1, k1                        !< Hydrostatic: first interior cell on the normal.
   integer(I4P)                :: dd(3)                             !< Hydrostatic: offset from that cell to the ghost.
   real(R8P)                   :: q1(9), ws(9)                      !< Hydrostatic: its state, the wall signs.
   real(R8P)                   :: gs                                !< Hydrostatic: the potential drop g . s.

   ! the rows of this pass only (issue #90, F-A'): `bc_rows_gpu(row_first:row_last,crown)` lists them; scanning the
   ! whole crown map for them made a launch of a pass with few rows cost as much as one with all of them (the kernel
   ! runs at 255 registers, 2 thread blocks per SM). The pass test below is kept: it is true for every listed row.
   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(local_map_bc_crown_gpu,bc_rows_gpu,q_inflow_gpu,wall_sign_gpu,dxyz_gpu,q_gpu,phi_gpu) &
   !$acc& firstprivate(ni,nj,nk,nv,nv_mhd,nv_entropy,crown,pass,face_kind,gamma,R,psi_energy) &
   !$acc& firstprivate(wall_velocity,wall_temperature,g,has_phi) &
   !$acc& private(b,c,i,j,k,idelta,jdelta,kdelta,bc_type,fec,face,iref,jref,kref,qd,qw,i1,j1,k1,dd,q1,ws,gs)
   !$omp OMPLOOP DEVICEPTR(local_map_bc_crown_gpu,bc_rows_gpu,q_inflow_gpu,wall_sign_gpu,dxyz_gpu,q_gpu,phi_gpu) &
   !$omp& firstprivate(ni,nj,nk,nv,nv_mhd,nv_entropy,crown,pass,face_kind,gamma,R,psi_energy) &
   !$omp& firstprivate(wall_velocity,wall_temperature,g,has_phi) &
   !$omp& private(b,c,i,j,k,idelta,jdelta,kdelta,bc_type,fec,face,iref,jref,kref,qd,qw,i1,j1,k1,dd,q1,ws,gs)
   do row=row_first, row_last
      c = bc_rows_gpu(row,crown)
      b = int(local_map_bc_crown_gpu(c,1,crown), I4P)
      fec = int(local_map_bc_crown_gpu(c,9,crown), I4P)
      if (b > 0_I4P) then
      if (abs(FEC_TO_DELTA(1,fec)) + abs(FEC_TO_DELTA(2,fec)) + abs(FEC_TO_DELTA(3,fec)) == pass) then
         i       = int(local_map_bc_crown_gpu(c,2,crown), I4P)
         j       = int(local_map_bc_crown_gpu(c,3,crown), I4P)
         k       = int(local_map_bc_crown_gpu(c,4,crown), I4P)
         idelta  = int(local_map_bc_crown_gpu(c,5,crown), I4P)
         jdelta  = int(local_map_bc_crown_gpu(c,6,crown), I4P)
         kdelta  = int(local_map_bc_crown_gpu(c,7,crown), I4P)
         bc_type = int(local_map_bc_crown_gpu(c,8,crown), I4P)
         iref = i - idelta ; jref = j - jdelta ; kref = k - kdelta
         if (pass == 1_I4P) then
            face = FEC_1_6_ARRAY(fec)
            ! extrapolation: the first interior cell along the normal, not the previous ghost, which is a row of this
            ! launch (issue #65 P1: the chain raced); the donor then lies beyond block interfaces only
            if (bc_type == BC_EXTRAPOLATION) &
               call realm_edge_donor(face=face, face_kind=BC_EXTRAPOLATION, ni=ni, nj=nj, nk=nk, i=i, j=j, k=k, &
                                     iref=iref, jref=jref, kref=kref)
            if (bc_type == BC_WALL_INVISCID .or. bc_type == BC_WALL_NOSLIP .or. bc_type == BC_WALL_ISOTHERMAL .or. &
                bc_type == BC_HYDROSTATIC) then
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
         else
            ! beyond two or three realm faces: the kind of one physical face on the donor; beyond seams only (outside
            ! the forest) the inward diagonal, as set above
            face = realm_edge_face(fec=fec, face_kind=face_kind)
            bc_type = BC_EXTRAPOLATION
            if (face > 0_I4P) then
               bc_type = face_kind(face)
               call realm_edge_donor(face=face, face_kind=face_kind(face), ni=ni, nj=nj, nk=nk, i=i, j=j, k=k, &
                                     iref=iref, jref=jref, kref=kref)
            endif
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
         elseif (bc_type == BC_WALL_NOSLIP .or. bc_type == BC_WALL_ISOTHERMAL) then
            !$acc loop seq
            do v=1, nv
               qd(v) = q_gpu(b,iref,jref,kref,v)
            enddo
            call wall_noslip_ghost(nv=nv, nv_mhd=nv_mhd, nv_entropy=nv_entropy, d=(face+1)/2, gamma=gamma, R=R,      &
                                   psi_energy=psi_energy,                                                     &
                                   isothermal=bc_type == BC_WALL_ISOTHERMAL, wall_velocity=wall_velocity(:,face), &
                                   wall_temperature=wall_temperature(face), q=qd(1:nv), qg=qw(1:nv))
            !$acc loop seq
            do v=1, nv
               q_gpu(b,i,j,k,v) = qw(v)
            enddo
         elseif (bc_type == BC_HYDROSTATIC) then
            ! the mirror (iref, jref, kref) gives the velocity, the field and psi; the first interior cell along the
            ! normal the density and the temperature (issue #73, M7-G0)
            call realm_edge_donor(face=face, face_kind=BC_EXTRAPOLATION, ni=ni, nj=nj, nk=nk, i=i, j=j, k=k, &
                                  iref=i1, jref=j1, kref=k1)
            !$acc loop seq
            do v=1, nv
               qd(v) = q_gpu(b,iref,jref,kref,v)
               q1(v) = q_gpu(b,i1,j1,k1,v)
               ws(v) = wall_sign_gpu(v,(face+1)/2)
            enddo
            dd = [i - i1, j - j1, k - k1]
            gs = g((face+1)/2) * real(dd((face+1)/2), R8P) * dxyz_gpu(b,(face+1)/2)
            ! with a cell potential (issue #81) its drop, exact whatever the potential: g . s = -(Phi_ghost - Phi_1)
            if (has_phi) gs = gs - (phi_gpu(b,i,j,k,1) - phi_gpu(b,i1,j1,k1,1))
            call hydrostatic_ghost(nv=nv, nv_mhd=nv_mhd, nv_entropy=nv_entropy, gamma=gamma, psi_energy=psi_energy, &
                                   gs=gs,                                                                      &
                                   wall_sign=ws(1:nv), qm=qd(1:nv), q1=q1(1:nv), qg=qw(1:nv))
            !$acc loop seq
            do v=1, nv
               q_gpu(b,i,j,k,v) = qw(v)
            enddo
         endif
      endif
      endif
   enddo
   endsubroutine set_boundary_conditions_dev
endmodule adam_flume_fnl_kernels
