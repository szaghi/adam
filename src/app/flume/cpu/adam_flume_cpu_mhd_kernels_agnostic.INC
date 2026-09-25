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
