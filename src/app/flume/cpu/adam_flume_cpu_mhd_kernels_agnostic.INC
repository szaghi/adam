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
