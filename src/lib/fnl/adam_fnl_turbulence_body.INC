   ! Body of the FNL eddy viscosity kernel of one model: ADAM_SGS_FORMULA names the pointwise routine (see the CPU body).
   integer(I4P), intent(in)    :: ir, iu                                !< Auxiliary indexes of the density and of u.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),    intent(in)    :: c                                     !< Model constant.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                       !< Blocks space steps [nb, 3].
   logical,      intent(in)    :: is_null(3)                            !< Null directions.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: mu_t_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Eddy viscosity, [nb, i, j, k, 1].
   integer(I4P)                :: s1, s2, s3                            !< Steps along x, y, z: 1 active, 0 null.
   real(R8P)                   :: w1, w2, w3                            !< Weights along x, y, z: 1 active, 0 null.
   real(R8P)                   :: rexp                                  !< 1 over the active directions.
   real(R8P)                   :: dl(3)                                 !< Private steps, 0 along null.
   real(R8P)                   :: rdx(3)                                !< Private inverse steps, 0 along null.
   real(R8P)                   :: delta                                 !< Private filter width.
   real(R8P)                   :: g(3,3)                                !< Private velocity gradient.
   real(R8P)                   :: nu                                    !< Private kinematic eddy viscosity.
   integer(I4P)                :: b, i, j, k, a                         !< Counters.

   s1 = merge(0_I4P, 1_I4P, is_null(1)) ; s2 = merge(0_I4P, 1_I4P, is_null(2)) ; s3 = merge(0_I4P, 1_I4P, is_null(3))
   w1 = real(s1, R8P) ; w2 = real(s2, R8P) ; w3 = real(s3, R8P)
   rexp = 1._R8P / real(max(1_I4P, s1 + s2 + s3), R8P)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_aux_gpu,mu_t_gpu)    &
   !$acc& firstprivate(ir,iu,ni,nj,nk,blocks_number,c,s1,s2,s3,w1,w2,w3,rexp) private(dl,rdx,delta,g,nu)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_aux_gpu,mu_t_gpu)                                   &
   !$omp& firstprivate(ir,iu,ni,nj,nk,blocks_number,c,s1,s2,s3,w1,w2,w3,rexp) private(dl,rdx,delta,g,nu)
   do k=1-2*s3, nk+2*s3
   do j=1-2*s2, nj+2*s2
   do i=1-2*s1, ni+2*s1
   do b=1, blocks_number
      dl(1) = w1 * dxyz_gpu(b,1) ; dl(2) = w2 * dxyz_gpu(b,2) ; dl(3) = w3 * dxyz_gpu(b,3)
      rdx(1) = w1 / dxyz_gpu(b,1) ; rdx(2) = w2 / dxyz_gpu(b,2) ; rdx(3) = w3 / dxyz_gpu(b,3)
      delta = (dxyz_gpu(b,1)**s1 * dxyz_gpu(b,2)**s2 * dxyz_gpu(b,3)**s3)**rexp
      !$acc loop seq
      do a=1, 3
         g(a,1) = 0.5_R8P * rdx(1) * (q_aux_gpu(b,i+s1,j,k,iu+a-1) - q_aux_gpu(b,i-s1,j,k,iu+a-1))
         g(a,2) = 0.5_R8P * rdx(2) * (q_aux_gpu(b,i,j+s2,k,iu+a-1) - q_aux_gpu(b,i,j-s2,k,iu+a-1))
         g(a,3) = 0.5_R8P * rdx(3) * (q_aux_gpu(b,i,j,k+s3,iu+a-1) - q_aux_gpu(b,i,j,k-s3,iu+a-1))
      enddo
      call ADAM_SGS_FORMULA(c=c, dl=dl, delta=delta, g=g, nu=nu)
      mu_t_gpu(b,i,j,k,1) = q_aux_gpu(b,i,j,k,ir) * nu
   enddo
   enddo
   enddo
   enddo
