   ! Model-agnostic CPU auxiliary-variables and non-finite-count loops: the bodies shared by every per-model kernel
   ! module (issue #41, section 4). Included in the `contains` part of `adam_flume_cpu_<model>_kernels`, which provides, by use
   ! association, `NV_K`, `NV_AUX_K` (compile-time sizes of the local arrays) and the model's pointwise
   ! `conservative_to_auxiliary`.

   subroutine compute_q_aux(ni, nj, nk, ngc, blocks_number, gamma, R, q, q_aux)
   !< Compute the auxiliary variables on every cell, ghost cells included.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                             !< Specific heats ratio.
   real(R8P),    intent(in)    :: R                                 !< Gas constant.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   real(R8P),    intent(inout) :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   integer(I4P)                :: b, i, j, k                        !< Counters.

   !$omp parallel do collapse(4) default(firstprivate) shared(q, q_aux)
   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               call conservative_to_auxiliary(gamma=gamma, R=R, q=q(:,i,j,k,b), qa=q_aux(:,i,j,k,b))
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_q_aux

   subroutine count_nonfinite(ni, nj, nk, ngc, blocks_number, q, n)
   !< Count the non-finite (NaN or infinite) conservative values of the interior cells (issue #45): the exponent field
   !< has all its bits set. Tested on the bit pattern, so that no floating point mode can fold the test away.
   integer(I4P), intent(in)  :: ni, nj, nk, ngc               !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                 !< Actual blocks number.
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   integer(I8P), intent(out) :: n                             !< Non-finite values number.
   integer(I8P), parameter   :: EXPONENT_MASK=2047_I8P        !< The 11 exponent bits of a binary64.
   integer(I4P)              :: b, i, j, k, v                 !< Counters.

   n = 0_I8P
   !$omp parallel do collapse(4) default(firstprivate) shared(q) reduction(+:n)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do v=1, NV_K
                  if (iand(ishft(transfer(q(v,i,j,k,b), 0_I8P), -52), EXPONENT_MASK) == EXPONENT_MASK) n = n + 1_I8P
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine count_nonfinite
