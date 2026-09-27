   ! MHD CPU loops shared by the MHD kernel modules (`adam_flume_cpu_mhd_kernels`, `adam_flume_cpu_mhd_glm_kernels`),
   ! which provide `NV_AUX_K` and the MHD pointwise `conservative_to_auxiliary`, `mhd_fast_speed` by use association.

   subroutine compute_lambda_max(ni, nj, nk, ngc, blocks_number, gamma, R, dxyz, is_null, q, lambda_max)
   !< Compute `max(sum_d (|u_d| + c_{f,d}) / dx_d)` over the interior cells (null directions excluded), `c_{f,d}` the
   !< fast magnetosonic speed along `d` (issue #41, section 3.5).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc               !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                 !< Actual blocks number.
   real(R8P),    intent(in)  :: gamma                         !< Specific heats ratio.
   real(R8P),    intent(in)  :: R                             !< Gas constant.
   real(R8P),    intent(in)  :: dxyz(1:,1:)                   !< Blocks space steps [3, nb].
   logical,      intent(in)  :: is_null(3)                    !< Null directions.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: lambda_max                    !< Maximum of sum_d (|u_d| + c_{f,d}) / dx_d.
   real(R8P)                 :: qa(NV_AUX_K)                  !< Auxiliary variables of one cell.
   integer(I4P)              :: b, i, j, k                    !< Counters.
   integer(I4P)              :: d                             !< Direction counter.

   lambda_max = 0._R8P
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q) reduction(max:lambda_max)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               call conservative_to_auxiliary(gamma=gamma, R=R, q=q(:,i,j,k,b), qa=qa)
               lambda_max = max(lambda_max, sum([((abs(qa(IA_U+d-1)) + mhd_fast_speed(d=d, qa=qa)) / dxyz(d,b), d=1, 3)], &
                                                mask=.not.is_null))
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_lambda_max

   subroutine apply_floors(ni, nj, nk, ngc, blocks_number, gamma, R, rho_floor, p_floor, q, floored, nonpositive, rho_min, &
                           p_min)
   !< Apply the MHD positivity floors to the interior cells of a stage state, before its ghost exchange (issue #41,
   !< section 3.8): a cell with `rho < rho_floor` or `p < p_floor` gets `rho = max(rho, rho_floor)` at fixed velocity and
   !< `p = max(p, p_floor)` (internal energy added). Interior only: the exchange then carries the floored values to the
   !< ghost cells, and ghost cells that no map fills (never read by the stencils) are not inspected. Returns the floored
   !< cells, the cells with a non-positive density or pressure (counted before flooring) and the minimum density and
   !< pressure. Floors are not conservative: the caller logs the counters.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                             !< Specific heats ratio.
   real(R8P),    intent(in)    :: R                                 !< Gas constant.
   real(R8P),    intent(in)    :: rho_floor, p_floor                !< Density and pressure floors.
   real(R8P),    intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   integer(I4P), intent(out)   :: floored                           !< Floored cells.
   integer(I4P), intent(out)   :: nonpositive                       !< Cells with non-positive density or pressure.
   real(R8P),    intent(out)   :: rho_min, p_min                    !< Minimum density and pressure.
   real(R8P)                   :: qa(NV_AUX_K)                      !< Auxiliary variables of one cell.
   real(R8P)                   :: rho, p                            !< Density and pressure of one cell.
   integer(I4P)                :: b, i, j, k                        !< Counters.

   floored = 0_I4P ; nonpositive = 0_I4P ; rho_min = huge(1._R8P) ; p_min = huge(1._R8P)
   !$omp parallel do collapse(4) default(firstprivate) shared(q) reduction(+:floored,nonpositive) &
   !$omp& reduction(min:rho_min,p_min)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               call conservative_to_auxiliary(gamma=gamma, R=R, q=q(:,i,j,k,b), qa=qa)
               rho = qa(IA_R)
               p   = qa(IA_P)
               rho_min = min(rho_min, rho)
               p_min   = min(p_min, p)
               if (rho <= 0._R8P .or. p <= 0._R8P) nonpositive = nonpositive + 1_I4P
               if (rho < rho_floor .or. p < p_floor) then
                  floored = floored + 1_I4P
                  rho = max(rho, rho_floor)
                  p   = max(p, p_floor)
                  q(IQ_R, i,j,k,b) = rho
                  q(IQ_RU,i,j,k,b) = rho * qa(IA_U)
                  q(IQ_RV,i,j,k,b) = rho * qa(IA_V)
                  q(IQ_RW,i,j,k,b) = rho * qa(IA_W)
                  q(IQ_RE,i,j,k,b) = p / (gamma - 1._R8P)                                                     + &
                                     0.5_R8P * mhd_sum3(q(IQ_RU,i,j,k,b)**2, q(IQ_RV,i,j,k,b)**2, q(IQ_RW,i,j,k,b)**2) / rho + &
                                     0.5_R8P * mhd_sum3(q(IQ_BX,i,j,k,b)**2, q(IQ_BY,i,j,k,b)**2, q(IQ_BZ,i,j,k,b)**2)
               endif
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine apply_floors

   subroutine compute_speed_max(ni, nj, nk, ngc, blocks_number, gamma, R, is_null, q, speed_max)
   !< Compute `max(|u_d| + c_{f,d})` over the interior cells and the active directions, the fastest wave speed (the GLM
   !< `c_h` check, issue #41, section 3.5).
   integer(I4P), intent(in)  :: ni, nj, nk, ngc               !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                 !< Actual blocks number.
   real(R8P),    intent(in)  :: gamma                         !< Specific heats ratio.
   real(R8P),    intent(in)  :: R                             !< Gas constant.
   logical,      intent(in)  :: is_null(3)                    !< Null directions.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: speed_max                     !< Maximum of |u_d| + c_{f,d}.
   real(R8P)                 :: wx, wy, wz                    !< Direction weights: 1 active, 0 null.
   real(R8P)                 :: qa(NV_AUX_K)                  !< Auxiliary variables of one cell.
   integer(I4P)              :: b, i, j, k                    !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   speed_max = 0._R8P
   !$omp parallel do collapse(4) default(firstprivate) shared(q) reduction(max:speed_max)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               call conservative_to_auxiliary(gamma=gamma, R=R, q=q(:,i,j,k,b), qa=qa)
               speed_max = max(speed_max, wx * (abs(qa(IA_U)) + mhd_fast_speed(d=1_I4P, qa=qa)), &
                                          wy * (abs(qa(IA_V)) + mhd_fast_speed(d=2_I4P, qa=qa)), &
                                          wz * (abs(qa(IA_W)) + mhd_fast_speed(d=3_I4P, qa=qa)))
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_speed_max

   subroutine compute_divb_norms(ni, nj, nk, ngc, blocks_number, hs, band, dxyz, is_null, seam, q, divb_max, divb_l1, &
                                 divb_seam_max)
   !< Compute the div(B) norms over the interior cells (issue #41, D-10): `max|div B|`, `sum |div B| dV` and the maximum
   !< over the seam-local cells (within `band` cells of a block face with `seam(b,face) /= 0`, a 2:1 AMR or inter-realm
   !< seam). `div B` is the centred finite difference of the library (half stencil `hs`), null directions weighted zero;
   !< the cell volume includes the null directions (issue #37). The caller has refreshed the ghost cells.
   integer(I4P), intent(in)  :: ni, nj, nk, ngc               !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                 !< Actual blocks number.
   integer(I4P), intent(in)  :: hs                            !< Finite difference half stencil, <= ngc.
   integer(I4P), intent(in)  :: band                          !< Seam-local band width, in cells.
   real(R8P),    intent(in)  :: dxyz(1:,1:)                   !< Blocks space steps [3, nb].
   logical,      intent(in)  :: is_null(3)                    !< Null directions.
   integer(I4P), intent(in)  :: seam(1:,1:)                   !< Seam faces flags [nb, 6] (0: not a seam face).
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: divb_max                      !< Maximum of |div B|.
   real(R8P),    intent(out) :: divb_l1                       !< Sum of |div B| dV.
   real(R8P),    intent(out) :: divb_seam_max                 !< Maximum of |div B| over the seam-local cells.
   real(R8P)                 :: wx, wy, wz                    !< Direction weights: 1 active, 0 null.
   real(R8P)                 :: dbx, dby, dbz                 !< Derivatives dBx/dx, dBy/dy, dBz/dz.
   real(R8P)                 :: divb                          !< |div B| of one cell.
   logical                   :: is_seam                       !< Seam-local cell.
   integer(I4P)              :: b, i, j, k                    !< Counters.

   wx = merge(0._R8P, 1._R8P, is_null(1))
   wy = merge(0._R8P, 1._R8P, is_null(2))
   wz = merge(0._R8P, 1._R8P, is_null(3))
   divb_max = 0._R8P ; divb_l1 = 0._R8P ; divb_seam_max = 0._R8P
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, seam, q) &
   !$omp& reduction(max:divb_max,divb_seam_max) reduction(+:divb_l1)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               call compute_derivative1_fd_centered(s=hs, ds=dxyz(1,b), q=q(IQ_BX,i-hs:i+hs,j,k,b), dq_ds=dbx)
               call compute_derivative1_fd_centered(s=hs, ds=dxyz(2,b), q=q(IQ_BY,i,j-hs:j+hs,k,b), dq_ds=dby)
               call compute_derivative1_fd_centered(s=hs, ds=dxyz(3,b), q=q(IQ_BZ,i,j,k-hs:k+hs,b), dq_ds=dbz)
               divb = abs(wx * dbx + wy * dby + wz * dbz)
               divb_max = max(divb_max, divb)
               divb_l1  = divb_l1 + divb * dxyz(1,b) * dxyz(2,b) * dxyz(3,b)
               is_seam = (seam(b,1) /= 0_I4P .and. i <= band) .or. (seam(b,2) /= 0_I4P .and. i > ni - band) .or. &
                         (seam(b,3) /= 0_I4P .and. j <= band) .or. (seam(b,4) /= 0_I4P .and. j > nj - band) .or. &
                         (seam(b,5) /= 0_I4P .and. k <= band) .or. (seam(b,6) /= 0_I4P .and. k > nk - band)
               if (is_seam) divb_seam_max = max(divb_seam_max, divb)
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_divb_norms
