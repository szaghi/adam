   ! MHD FNL kernels shared by the MHD kernel modules (`adam_flume_fnl_mhd_kernels`, `adam_flume_fnl_mhd_glm_kernels`),
   ! which provide `NV_K`, `NV_AUX_K` and the MHD pointwise `conservative_to_auxiliary`, `mhd_fast_speed` by use
   ! association.

   subroutine compute_conservation_dev(ni, nj, nk, ngc, blocks_number, dxyz_gpu, q_gpu, integrals)
   !< Compute the volume integrals of the conservative variables (interior cells), one reduction per variable. The cell
   !< volume includes the null directions (issue #37).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Blocks space steps [nb, 3].
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: integrals(NV_K)                   !< Volume integrals.
   real(R8P)                 :: volume                            !< Cell volume.
   real(R8P)                 :: s                                 !< Reduction accumulator.
   integer(I4P)              :: b, i, j, k, v                     !< Counters.

   do v=1, NV_K
      s = 0._R8P
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu) &
      !$acc& firstprivate(ni,nj,nk,blocks_number,v) private(volume) reduction(+:s)
      !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu) &
      !$omp& firstprivate(ni,nj,nk,blocks_number,v) private(volume) reduction(+:s)
      do k=1, nk
      do j=1, nj
      do i=1, ni
      do b=1, blocks_number
         volume = dxyz_gpu(b,1) * dxyz_gpu(b,2) * dxyz_gpu(b,3)
         s = s + q_gpu(b,i,j,k,v) * volume
      enddo
      enddo
      enddo
      enddo
      integrals(v) = s
   enddo
   endsubroutine compute_conservation_dev

   subroutine compute_lambda_max_dev(ni, nj, nk, ngc, blocks_number, gamma, R, dxyz_gpu, is_null, q_gpu, lambda_max)
   !< Compute `max(sum_d (|u_d| + c_{f,d}) / dx_d)` over the interior cells (null directions excluded), `c_{f,d}` the
   !< fast magnetosonic speed along `d` (issue #41, section 3.5).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)  :: gamma                             !< Specific heats ratio.
   real(R8P),    intent(in)  :: R                                 !< Gas constant.
   real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Blocks space steps [nb, 3].
   logical,      intent(in)  :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: lambda_max                        !< Maximum of sum_d (|u_d| + c_{f,d}) / dx_d.
   real(R8P)                 :: wx, wy, wz                        !< Direction weights: 1 active, 0 null.
   real(R8P)                 :: q_(NV_K)                          !< Private conservative variables of one cell.
   real(R8P)                 :: qa_(NV_AUX_K)                     !< Private auxiliary variables of one cell.
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
      do v=1, NV_K
         q_(v) = q_gpu(b,i,j,k,v)
      enddo
      call conservative_to_auxiliary(gamma=gamma, R=R, q=q_, qa=qa_)
      lambda = wx * (abs(qa_(IA_U)) + mhd_fast_speed(d=1_I4P, qa=qa_)) / dxyz_gpu(b,1) + &
               wy * (abs(qa_(IA_V)) + mhd_fast_speed(d=2_I4P, qa=qa_)) / dxyz_gpu(b,2) + &
               wz * (abs(qa_(IA_W)) + mhd_fast_speed(d=3_I4P, qa=qa_)) / dxyz_gpu(b,3)
      lambda_max = max(lambda_max, lambda)
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_lambda_max_dev

   subroutine apply_floors_dev(ni, nj, nk, ngc, blocks_number, gamma, R, rho_floor, p_floor, q_gpu, floored, nonpositive, &
                               rho_min, p_min)
   !< Apply the MHD positivity floors to the interior cells of a stage state, before its ghost exchange (issue #41,
   !< section 3.8): a cell with `rho < rho_floor` or `p < p_floor` gets `rho = max(rho, rho_floor)` at fixed velocity and
   !< `p = max(p, p_floor)` (internal energy added; the model's `cleaning_energy` adapter adds `psi^2 / 2` under EGLM,
   !< zero otherwise). Interior only: the exchange then carries the floored values to the
   !< ghost cells, and ghost cells that no map fills (never read by the stencils) are not inspected. Returns the floored
   !< cells, the cells with a non-positive density or pressure (counted before flooring) and the minimum density and
   !< pressure. Floors are not conservative: the caller logs the counters.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                                 !< Specific heats ratio.
   real(R8P),    intent(in)    :: R                                     !< Gas constant.
   real(R8P),    intent(in)    :: rho_floor, p_floor                    !< Density and pressure floors.
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   integer(I4P), intent(out)   :: floored                               !< Floored cells.
   integer(I4P), intent(out)   :: nonpositive                           !< Cells with non-positive density or pressure.
   real(R8P),    intent(out)   :: rho_min, p_min                        !< Minimum density and pressure.
   real(R8P)                   :: q_(NV_K)                              !< Private conservative variables of one cell.
   real(R8P)                   :: qa_(NV_AUX_K)                         !< Private auxiliary variables of one cell.
   real(R8P)                   :: rho, p                                !< Density and pressure of one cell.
   integer(I4P)                :: b, i, j, k, v                         !< Counters.

   floored = 0_I4P ; nonpositive = 0_I4P ; rho_min = huge(1._R8P) ; p_min = huge(1._R8P)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu)                       &
   !$acc& firstprivate(ni,nj,nk,ngc,blocks_number,gamma,R,rho_floor,p_floor) private(q_,qa_,rho,p) &
   !$acc& reduction(+:floored,nonpositive) reduction(min:rho_min,p_min)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu) &
   !$omp& firstprivate(ni,nj,nk,ngc,blocks_number,gamma,R,rho_floor,p_floor) private(q_,qa_,rho,p) &
   !$omp& reduction(+:floored,nonpositive) reduction(min:rho_min,p_min)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do v=1, NV_K
         q_(v) = q_gpu(b,i,j,k,v)
      enddo
      call conservative_to_auxiliary(gamma=gamma, R=R, q=q_, qa=qa_)
      rho = qa_(IA_R)
      p   = qa_(IA_P)
      rho_min = min(rho_min, rho)
      p_min   = min(p_min, p)
      if (rho <= 0._R8P .or. p <= 0._R8P) nonpositive = nonpositive + 1_I4P
      if (rho < rho_floor .or. p < p_floor) then
         floored = floored + 1_I4P
         rho = max(rho, rho_floor)
         p   = max(p, p_floor)
         q_(IQ_R)  = rho
         q_(IQ_RU) = rho * qa_(IA_U)
         q_(IQ_RV) = rho * qa_(IA_V)
         q_(IQ_RW) = rho * qa_(IA_W)
         q_(IQ_RE) = p / (gamma - 1._R8P) + 0.5_R8P * mhd_sum3(q_(IQ_RU)**2, q_(IQ_RV)**2, q_(IQ_RW)**2) / rho + &
                     0.5_R8P * mhd_sum3(q_(IQ_BX)**2, q_(IQ_BY)**2, q_(IQ_BZ)**2) + cleaning_energy(q=q_)
         !$acc loop seq
         do v=1, NV_K
            q_gpu(b,i,j,k,v) = q_(v)
         enddo
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine apply_floors_dev

   subroutine compute_speed_max_dev(ni, nj, nk, ngc, blocks_number, gamma, R, is_null, q_gpu, speed_max)
   !< Compute `max(|u_d| + c_{f,d})` over the interior cells and the active directions, the fastest wave speed (the GLM
   !< `c_h` check, issue #41, section 3.5).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)  :: gamma                             !< Specific heats ratio.
   real(R8P),    intent(in)  :: R                                 !< Gas constant.
   logical,      intent(in)  :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: speed_max                         !< Maximum of |u_d| + c_{f,d}.
   real(R8P)                 :: wx, wy, wz                        !< Direction weights: 1 active, 0 null.
   real(R8P)                 :: q_(NV_K)                          !< Private conservative variables of one cell.
   real(R8P)                 :: qa_(NV_AUX_K)                     !< Private auxiliary variables of one cell.
   integer(I4P)              :: b, i, j, k, v                     !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   speed_max = 0._R8P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,gamma,R,wx,wy,wz) private(q_,qa_) reduction(max:speed_max)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,gamma,R,wx,wy,wz) private(q_,qa_) reduction(max:speed_max)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do v=1, NV_K
         q_(v) = q_gpu(b,i,j,k,v)
      enddo
      call conservative_to_auxiliary(gamma=gamma, R=R, q=q_, qa=qa_)
      speed_max = max(speed_max, wx * (abs(qa_(IA_U)) + mhd_fast_speed(d=1_I4P, qa=qa_)), &
                                 wy * (abs(qa_(IA_V)) + mhd_fast_speed(d=2_I4P, qa=qa_)), &
                                 wz * (abs(qa_(IA_W)) + mhd_fast_speed(d=3_I4P, qa=qa_)))
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_speed_max_dev

   subroutine compute_divb_norms_dev(ni, nj, nk, ngc, blocks_number, hs, band, dxyz_gpu, is_null, seam_gpu, q_gpu, &
                                     divb_max, divb_l1, divb_seam_max)
   !< Compute the div(B) norms over the interior cells (issue #41, D-10): `max|div B|`, `sum |div B| dV` and the maximum
   !< over the seam-local cells (within `band` cells of a block face with `seam(b,face) /= 0`, a 2:1 AMR or inter-realm
   !< seam). `div B` is the centred finite difference of the library (half stencil `hs <= HS_MAX`), null directions
   !< weighted zero; the cell volume includes the null directions (issue #37). The caller has refreshed the ghost cells.
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                     !< Actual blocks number.
   integer(I4P), intent(in)  :: hs                                !< Finite difference half stencil, <= min(ngc, HS_MAX).
   integer(I4P), intent(in)  :: band                              !< Seam-local band width, in cells.
   real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Blocks space steps [nb, 3].
   logical,      intent(in)  :: is_null(3)                        !< Null directions.
   integer(I4P), intent(in)  :: seam_gpu(1:,1:)                   !< Seam faces flags [nb, 6] (0: not a seam face).
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: divb_max                          !< Maximum of |div B|.
   real(R8P),    intent(out) :: divb_l1                           !< Sum of |div B| dV.
   real(R8P),    intent(out) :: divb_seam_max                     !< Maximum of |div B| over the seam-local cells.
   integer(I4P), parameter   :: HS_MAX=3_I4P                      !< Maximum half stencil (sixth order).
   real(R8P)                 :: sx(1-HS_MAX:1+HS_MAX)             !< Private Bx stencil along x.
   real(R8P)                 :: sy(1-HS_MAX:1+HS_MAX)             !< Private By stencil along y.
   real(R8P)                 :: sz(1-HS_MAX:1+HS_MAX)             !< Private Bz stencil along z.
   real(R8P)                 :: wx, wy, wz                        !< Direction weights: 1 active, 0 null.
   real(R8P)                 :: dbx, dby, dbz                     !< Derivatives dBx/dx, dBy/dy, dBz/dz.
   real(R8P)                 :: divb                              !< |div B| of one cell.
   logical                   :: is_seam                           !< Seam-local cell.
   integer(I4P)              :: b, i, j, k, m                     !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   divb_max = 0._R8P ; divb_l1 = 0._R8P ; divb_seam_max = 0._R8P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,seam_gpu,q_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,hs,band,wx,wy,wz) private(sx,sy,sz,dbx,dby,dbz,divb,is_seam) &
   !$acc& reduction(max:divb_max,divb_seam_max) reduction(+:divb_l1)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,seam_gpu,q_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,hs,band,wx,wy,wz) private(sx,sy,sz,dbx,dby,dbz,divb,is_seam) &
   !$omp& reduction(max:divb_max,divb_seam_max) reduction(+:divb_l1)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do m=-hs, hs
         sx(1+m) = q_gpu(b,i+m,j,k,IQ_BX)
         sy(1+m) = q_gpu(b,i,j+m,k,IQ_BY)
         sz(1+m) = q_gpu(b,i,j,k+m,IQ_BZ)
      enddo
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,1), q=sx(1-hs:1+hs), dq_ds=dbx)
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,2), q=sy(1-hs:1+hs), dq_ds=dby)
      call compute_derivative1_fd_centered(s=hs, ds=dxyz_gpu(b,3), q=sz(1-hs:1+hs), dq_ds=dbz)
      divb = abs(wx * dbx + wy * dby + wz * dbz)
      divb_max = max(divb_max, divb)
      divb_l1  = divb_l1 + divb * dxyz_gpu(b,1) * dxyz_gpu(b,2) * dxyz_gpu(b,3)
      is_seam = (seam_gpu(b,1) /= 0_I4P .and. i <= band) .or. (seam_gpu(b,2) /= 0_I4P .and. i > ni - band) .or. &
                (seam_gpu(b,3) /= 0_I4P .and. j <= band) .or. (seam_gpu(b,4) /= 0_I4P .and. j > nj - band) .or. &
                (seam_gpu(b,5) /= 0_I4P .and. k <= band) .or. (seam_gpu(b,6) /= 0_I4P .and. k > nk - band)
      if (is_seam) divb_seam_max = max(divb_seam_max, divb)
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_divb_norms_dev
