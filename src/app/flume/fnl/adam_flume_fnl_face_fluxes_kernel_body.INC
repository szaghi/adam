   ! Body of the FNL face-flux kernel (flux splitting), one source for every stencil width: included by the kernels
   ! `compute_face_fluxes_s<N>_dev` of `adam_flume_fnl_face_kernels_agnostic.INC` (S a named constant, N = 1 ... S_MAX:
   ! the stencil loops of the kernel and of the device routines -gpu=lto inlines into it have constant bounds; measured
   ! faster, see the face kernels file), after their declaration of `S` and their `use` of SETH's WENO-JS kernel of
   ! that width as `weno_width` (issue #139, S3b: SETH's stencil convention, `vp`/`vm` by offset from the face).
   ! `S` is not in the data clauses: a named constant cannot be, and a scalar is firstprivate by default in an
   ! OpenACC parallel and in an OpenMP target region.
   integer(I4P), intent(in)    :: d                                       !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                              !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                         !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                           !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                                   !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                      !< GLM cleaning speed (models with GLM only).
   logical,      intent(in)    :: is_characteristic                       !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)    :: weno_zeps                               !< WENO parameter avoiding division by zero.
   real(R8P),    intent(in)    :: weno_sigma                              !< WENO descaler switch (0 Jiang-Shu, 1 SI).
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)         !< Face fluxes of direction `d`.
   real(R8P)                   :: qs(NV_K,1-S_MAX:S_MAX)                  !< Private stencil conservative variables.
   real(R8P)                   :: qas(NV_AUX_K,1-S_MAX:S_MAX)             !< Private stencil auxiliary variables.
   real(R8P)                   :: fsplit(2,1-S_MAX:S_MAX-1,NV_K)          !< Private split fields.
   real(R8P)                   :: er(NV_K,NV_K)                           !< Private right eigenvectors.
   real(R8P)                   :: mu(NV_K)                                !< Private magnitude of each split field.
   real(R8P)                   :: vp(1-S:S), vm(1-S:S)!< Private stencil of one field, left and right bias.
   real(R8P)                   :: vr(2,NV_K)                              !< Private reconstructed split fields.
   real(R8P)                   :: flux(NV_K)                              !< Private face flux.
   integer(I4P)                :: b, i, j, k, m, f, w                     !< Counters.

   !$acc parallel loop independent gang vector collapse(4)                                                        &
   !$acc& DEVICEVAR(q_gpu,q_aux_gpu,fl_gpu)                                                                        &
   !$acc& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,gamma,ch,is_characteristic,weno_zeps,weno_sigma)         &
   !$acc& private(qs,qas,fsplit,er,mu,vp,vm,vr,flux)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu,q_aux_gpu,fl_gpu)                                  &
   !$omp& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,gamma,ch,is_characteristic,weno_zeps,weno_sigma) &
   !$omp& private(qs,qas,fsplit,er,mu,vp,vm,vr,flux)
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
                             fsplit=fsplit, er=er, mu=mu)
      !$acc loop seq
      do f=1, NV_K
         !$acc loop seq
         do m=1-S, S-1
            vm(m+1) = fsplit(1,m,f)
            vp(m)   = fsplit(2,m,f)
         enddo
         call weno_width(vp=vp, vm=vm, eps=weno_zeps, rmu=weno_sigma / mu(f) + (1._R8P - weno_sigma), &
                         rp=vr(2,f), rm=vr(1,f))
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
