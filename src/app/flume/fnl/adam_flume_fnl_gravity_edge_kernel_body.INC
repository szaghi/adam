   ! Body of the FNL gravity edge-correction kernel (flux splitting), one source for every stencil width: included by
   ! the kernels `add_gravity_edge_correction_s<N>_dev` of `adam_flume_fnl_face_kernels_agnostic.INC` (S a named
   ! constant, N = 1 ... S_MAX), after their declaration of `S` and their `use` of SETH's WENO-JS kernel of that width
   ! as `weno_width` (issue #139, S3b); `S` is not in the data clauses (see adam_flume_fnl_face_fluxes_kernel_body.INC).
   integer(I4P), intent(in)    :: d                                       !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                              !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                         !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                           !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                                   !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                      !< GLM cleaning speed (models with GLM only).
   logical,      intent(in)    :: is_characteristic                       !< Characteristic (or conservative) variables.
   real(R8P),    intent(in)    :: weno_zeps                               !< WENO parameter avoiding division by zero.
   real(R8P),    intent(in)    :: weno_sigma                              !< WENO descaler switch (0 Jiang-Shu, 1 SI).
   integer(I4P), intent(in)    :: ir                                      !< Index of the density (mass) variable.
   real(R8P),    intent(in)    :: gd                                      !< Gravity along `d`.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Auxiliary variables.
   real(R8P),    intent(inout) :: gc_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Gravity energy correction [..., 1].
   integer(I4P)                :: nd, n1, n2                              !< Cells along `d` and across it.
   integer(I4P)                :: fd, cd                                  !< Private face and cell index along `d`.
   integer(I4P)                :: b, i, j, k, m, f, w, t1, t2, side       !< Counters.

   nd = di * ni + dj * nj + dk * nk
   n1 = merge(nj, ni, di == 1_I4P)
   n2 = merge(nj, nk, dk == 1_I4P)
   !$acc parallel loop independent gang vector collapse(3)                                                        &
   !$acc& DEVICEVAR(q_gpu,q_aux_gpu,gc_gpu)                                                                        &
   !$acc& firstprivate(d,di,dj,dk,nd,n1,n2,blocks_number,gamma,ch,is_characteristic,weno_zeps,weno_sigma,ir,gd)   &
   !$acc& private(fd,cd,i,j,k)
   !$omp OMPLOOP collapse(3) DEVICEPTR(q_gpu,q_aux_gpu,gc_gpu)                                  &
   !$omp& firstprivate(d,di,dj,dk,nd,n1,n2,blocks_number,gamma,ch,is_characteristic,weno_zeps,weno_sigma,ir,gd) &
   !$omp& private(fd,cd,i,j,k)
   do t2=1, n2
   do t1=1, n1
   do b=1, blocks_number
      block
      ! locals of the loop body, not private clauses: flang maps a private array of a target region to and from
      ! the host at every launch (a synchronous copy each; #157, the MI300X 2-D step)
      real(R8P) :: qs(NV_K,1-S_MAX:S_MAX)
      real(R8P) :: qas(NV_AUX_K,1-S_MAX:S_MAX)
      real(R8P) :: fsplit(2,1-S_MAX:S_MAX-1,NV_K)
      real(R8P) :: er(NV_K,NV_K)
      real(R8P) :: mu(NV_K)
      real(R8P) :: vp(1-S:S)
      real(R8P) :: vm(1-S:S)
      real(R8P) :: vr(2,NV_K)
      real(R8P) :: flux(NV_K)
      !$acc loop seq
      do side=1, 2
         fd = nd + 1_I4P
         cd = nd
         if (side == 1) then
            fd = -1_I4P
            cd = 1_I4P
         endif
         i = di * fd + (1 - di) * t1
         j = dj * fd + di * t1 + dk * t2
         k = dk * fd + (1 - dk) * t2
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
         i = di * cd + (1 - di) * t1
         j = dj * cd + di * t1 + dk * t2
         k = dk * cd + (1 - dk) * t2
         gc_gpu(b,i,j,k,1) = gc_gpu(b,i,j,k,1) - gd * flux(ir) / 24._R8P
      enddo
      endblock
   enddo
   enddo
   enddo
