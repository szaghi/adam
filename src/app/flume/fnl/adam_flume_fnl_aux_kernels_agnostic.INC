   ! Model-agnostic FNL auxiliary-variables kernel: the body shared by every per-model kernel module (issue #41,
   ! section 4). Included in the `contains` part of `adam_flume_fnl_<model>_kernels`, which provides, by use
   ! association, `NV_K`, `NV_AUX_K` (compile-time sizes of the private arrays) and the model's pointwise
   ! `conservative_to_auxiliary`.

   subroutine compute_q_aux_dev(ni, nj, nk, ngc, blocks_number, gamma, R, q_gpu, q_aux_gpu)
   !< Compute the auxiliary variables on every cell, ghost cells included.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                                 !< Specific heats ratio.
   real(R8P),    intent(in)    :: R                                     !< Gas constant.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   real(R8P),    intent(inout) :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P)                   :: q_(NV_K)                              !< Private conservative variables of one cell.
   real(R8P)                   :: qa_(NV_AUX_K)                         !< Private auxiliary variables of one cell.
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
      do v=1, NV_K
         q_(v) = q_gpu(b,i,j,k,v)
      enddo
      call conservative_to_auxiliary(gamma=gamma, R=R, q=q_, qa=qa_)
      !$acc loop seq
      do v=1, NV_AUX_K
         q_aux_gpu(b,i,j,k,v) = qa_(v)
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_q_aux_dev
