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
