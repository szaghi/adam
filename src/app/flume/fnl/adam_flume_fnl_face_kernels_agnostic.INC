   ! Model-agnostic FNL face-flux kernel: the body shared by every per-model kernel module (issue #41, section 4).
   !
   ! Included in the `contains` part of `adam_flume_fnl_<model>_kernels`, which provides, by use association:
   ! `NV_K`, `NV_AUX_K` (the model's conservative and auxiliary variables numbers, compile-time constants that size
   ! every private array: no model branching and no maximum-model buffers in the kernels), `S_MAX`, the model's
   ! split adapter `face_split_fluxes(gamma, ch, ...)` (one argument list for every model: the GLM cleaning speed `ch` is
   ! ignored by the models without cleaning), `compute_face_flux_back_projection`, and `weno_reconstruct_upwind_dev`.
   ! The host selects the model's module by `select case`, outside the kernels.

   subroutine compute_face_fluxes_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, S, gamma, ch, is_characteristic, &
                                      weno_a_gpu, weno_p_gpu, weno_d_gpu, weno_zeps, q_gpu, q_aux_gpu, fl_gpu)
   !< Compute the WENO face fluxes of direction `d`: face `(i,j,k)` lies between cells `(i,j,k)` and
   !< `(i+di,j+dj,k+dk)`, so the face array starts at index 0 along `d`.
   !<
   !< Per face (device twin of the CPU `compute_face_fluxes`): gather the stencil `m = 1-S ... S` into constant-bound
   !< privates, project and split it, reconstruct every field (the private `v` is the packed stencil of the WENO
   !< primitive), back-project.
   integer(I4P), intent(in)    :: d                                       !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                              !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                         !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                           !< Actual blocks number.
   integer(I4P), intent(in)    :: S                                       !< WENO stencil half-width.
   real(R8P),    intent(in)    :: gamma                                   !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                      !< GLM cleaning speed (models with GLM only).
   logical,      intent(in)    :: is_characteristic                       !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)    :: weno_a_gpu(1:,0:,1:)                    !< WENO optimal weights.
   real(R8P),    intent(in)    :: weno_p_gpu(1:,0:,0:,1:)                 !< WENO polynomials coefficients.
   real(R8P),    intent(in)    :: weno_d_gpu(0:,0:,0:,1:)                 !< WENO smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_zeps                               !< WENO parameter avoiding division by zero.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)         !< Face fluxes of direction `d`.
   real(R8P)                   :: qs(NV_K,1-S_MAX:S_MAX)                  !< Private stencil conservative variables.
   real(R8P)                   :: qas(NV_AUX_K,1-S_MAX:S_MAX)             !< Private stencil auxiliary variables.
   real(R8P)                   :: fsplit(2,1-S_MAX:S_MAX-1,NV_K)          !< Private split fields.
   real(R8P)                   :: er(NV_K,NV_K)                           !< Private right eigenvectors.
   real(R8P)                   :: v(2,2*S_MAX-1)                          !< Private packed stencil of one field.
   real(R8P)                   :: vr_(2)                                  !< Private reconstruction of one field.
   real(R8P)                   :: vr(2,NV_K)                              !< Private reconstructed split fields.
   real(R8P)                   :: flux(NV_K)                              !< Private face flux.
   integer(I4P)                :: b, i, j, k, m, f, w                     !< Counters.

   !$acc parallel loop independent gang vector collapse(4)                                                        &
   !$acc& DEVICEVAR(weno_a_gpu,weno_p_gpu,weno_d_gpu,q_gpu,q_aux_gpu,fl_gpu)                                       &
   !$acc& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,S,gamma,ch,is_characteristic,weno_zeps)                  &
   !$acc& private(qs,qas,fsplit,er,v,vr_,vr,flux)
   !$omp OMPLOOP collapse(4) DEVICEPTR(weno_a_gpu,weno_p_gpu,weno_d_gpu,q_gpu,q_aux_gpu,fl_gpu) &
   !$omp& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,S,gamma,ch,is_characteristic,weno_zeps) &
   !$omp& private(qs,qas,fsplit,er,v,vr_,vr,flux)
   do k=1-dk, nk
   do j=1-dj, nj
   do i=1-di, ni
   do b=1, blocks_number
      !$acc loop seq
      do m=1-S, S
         !$acc loop seq
         do w=1, NV_K
            qs(w,m) = q_gpu(b,i+m*di,j+m*dj,k+m*dk,w)
         enddo
         !$acc loop seq
         do w=1, NV_AUX_K
            qas(w,m) = q_aux_gpu(b,i+m*di,j+m*dj,k+m*dk,w)
         enddo
      enddo
      call face_split_fluxes(gamma=gamma, ch=ch, d=d, S=S, is_characteristic=is_characteristic, qs=qs, qas=qas, &
                             fsplit=fsplit, er=er)
      !$acc loop seq
      do f=1, NV_K
         !$acc loop seq
         do m=1-S, S-1
            v(1,m+S) = fsplit(1,m,f)
            v(2,m+S) = fsplit(2,m,f)
         enddo
         call weno_reconstruct_upwind_dev(S=S, weno_a=weno_a_gpu, weno_p=weno_p_gpu, weno_d=weno_d_gpu, &
                                          weno_zeps=weno_zeps, V=v, VR=vr_)
         vr(1,f) = vr_(1)
         vr(2,f) = vr_(2)
      enddo
      call compute_face_flux_back_projection(is_characteristic=is_characteristic, er=er, vr=vr, flux=flux)
      !$acc loop seq
      do w=1, NV_K
         fl_gpu(b,i,j,k,w) = flux(w)
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_face_fluxes_dev
