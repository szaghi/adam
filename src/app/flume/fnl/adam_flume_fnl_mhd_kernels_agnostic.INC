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
   !< `p = max(p, p_floor)` (internal energy added). Interior only: the exchange then carries the floored values to the
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
         q_(IQ_RE) = p / (gamma - 1._R8P) + 0.5_R8P * (q_(IQ_RU)**2 + q_(IQ_RV)**2 + q_(IQ_RW)**2) / rho + &
                     0.5_R8P * (q_(IQ_BX)**2 + q_(IQ_BY)**2 + q_(IQ_BZ)**2)
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
