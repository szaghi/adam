   ! Body of the FNL face-flux kernel `compute_riemann_face_fluxes_s<N>_dev` / `compute_riemann_face_fluxes_sn_dev` of the scheme
   ! `weno-riemann`, one source for every stencil width: included by the kernels of
   ! `adam_flume_fnl_riemann_face_kernels_agnostic.INC`, after their declaration of `S`.
   ! FNL_WENO_WRATIO_WIDTH (defined by the per-width kernels of `adam_flume_fnl_riemann_face_kernels_agnostic.INC` around
   ! this include): the name of the WENO reconstruction with the weights ratio of their width, S a named constant, whose call
   ! has no S argument. The generic kernel `_sn_dev` does not define it and calls `weno_reconstruct_upwind_wratio_dev` with S.
   ! A macro because one body is shared by all the widths and Fortran cannot pass a procedure name through an include (it
   ! renames no routine, no math, no bound). `S` is not in the data clauses: a named constant cannot be, and a scalar is
   ! firstprivate by default in an OpenACC parallel and in an OpenMP target region.
   integer(I4P), intent(in)    :: d                                       !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                              !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                         !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                           !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                                   !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                      !< GLM cleaning speed (models with GLM only).
   logical,      intent(in)    :: is_characteristic                       !< Characteristic (or primitive) variables.
   real(R8P),    intent(in)    :: weno_a_gpu(1:,0:,1:)                    !< WENO interpolation optimal weights.
   real(R8P),    intent(in)    :: weno_p_gpu(1:,0:,0:,1:)                 !< WENO interpolation polynomials coefficients.
   real(R8P),    intent(in)    :: weno_d_gpu(0:,0:,0:,1:)                 !< WENO smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_zeps                               !< WENO parameter avoiding division by zero.
   real(R8P),    intent(in)    :: weno_sigma                              !< WENO descaler switch (0 Jiang-Shu, 1 SI).
   real(R8P),    intent(in)    :: cc(3)                                   !< Correction coefficients.
   real(R8P),    intent(in)    :: tau                                     !< Sensor threshold.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Conservative variables.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)         !< Face fluxes of direction `d`.
   integer(I4P), intent(out)   :: fallbacks                               !< Faces whose Riemann solver fell back.
   real(R8P)                   :: qs(NV_K,1-S_MAX:S_MAX)                  !< Private stencil conservative variables.
   real(R8P)                   :: qas(NV_AUX_K,1-S_MAX:S_MAX)             !< Private stencil auxiliary variables.
   real(R8P)                   :: fint(2,1-S_MAX:S_MAX-1,NV_K)            !< Private fields to interpolate.
   real(R8P)                   :: er(NV_K,NV_K)                           !< Private right eigenvectors.
   real(R8P)                   :: mu(NV_K)                                !< Private magnitude of each field.
   real(R8P)                   :: v(2,2*S_MAX-1)                          !< Private packed stencil of one field.
   real(R8P)                   :: vr_(2)                                  !< Private interpolation of one field.
   real(R8P)                   :: wr(2)                                   !< Private WENO weights ratio of one field.
   real(R8P)                   :: vr(2,NV_K)                              !< Private interpolated fields.
   real(R8P)                   :: qL(NV_K), qR(NV_K)                      !< Private face states.
   real(R8P)                   :: fr(NV_K)                                !< Private Riemann flux.
   real(R8P)                   :: fc(NV_K,-1:2)                           !< Private physical fluxes of cells -1 ... 2.
   real(R8P)                   :: wmin                                    !< Private minimum WENO weights ratio.
   real(R8P)                   :: sensor                                  !< Private correction switch.
   logical                     :: fallback                                !< Private Riemann solver fallback flag.
   integer(I4P)                :: b, i, j, k, m, f, w                     !< Counters.

   fallbacks = 0_I4P
   !$acc parallel loop independent gang vector collapse(4) reduction(+:fallbacks)                                 &
   !$acc& DEVICEVAR(weno_a_gpu,weno_p_gpu,weno_d_gpu,q_gpu,q_aux_gpu,fl_gpu)                                       &
   !$acc& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,gamma,ch,is_characteristic,weno_zeps,weno_sigma,cc,tau) &
   !$acc& private(qs,qas,fint,er,mu,v,vr_,wr,vr,qL,qR,fr,fc,wmin,sensor,fallback)
   !$omp OMPLOOP collapse(4) reduction(+:fallbacks) &
   !$omp& DEVICEPTR(weno_a_gpu,weno_p_gpu,weno_d_gpu,q_gpu,q_aux_gpu,fl_gpu) &
   !$omp& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,gamma,ch,is_characteristic,weno_zeps,weno_sigma,cc,tau) &
   !$omp& private(qs,qas,fint,er,mu,v,vr_,wr,vr,qL,qR,fr,fc,wmin,sensor,fallback)
   do k=1-dk, nk
   do j=1-dj, nj
   do i=1-di, ni
   do b=1, blocks_number
      !$acc loop seq
      do m=min(1-S,-1), max(S,2)
         !$acc loop seq
         do w=1, NV_K
            qs(w,m) = q_gpu(b,i+m*di,j+m*dj,k+m*dk,w)
         enddo
         !$acc loop seq
         do w=1, NV_AUX_K
            qas(w,m) = q_aux_gpu(b,i+m*di,j+m*dj,k+m*dk,w)
         enddo
      enddo
      call face_interpolation_fields(gamma=gamma, ch=ch, d=d, S=S, is_characteristic=is_characteristic, qs=qs, &
                                     qas=qas, fint=fint, er=er, mu=mu)
      wmin = 1._R8P
      !$acc loop seq
      do f=1, NV_K
         !$acc loop seq
         do m=1-S, S-1
            v(1,m+S) = fint(1,m,f)
            v(2,m+S) = fint(2,m,f)
         enddo
#ifdef FNL_WENO_WRATIO_WIDTH
         call FNL_WENO_WRATIO_WIDTH(weno_a=weno_a_gpu, weno_p=weno_p_gpu, weno_d=weno_d_gpu, &
                                    weno_zeps=weno_zeps, weno_rmu=weno_sigma / mu(f) + (1._R8P - weno_sigma), &
                                    V=v, VR=vr_, WRATIO=wr)
#else
         call weno_reconstruct_upwind_wratio_dev(S=S, weno_a=weno_a_gpu, weno_p=weno_p_gpu, weno_d=weno_d_gpu, &
                                                 weno_zeps=weno_zeps, weno_rmu=weno_sigma / mu(f) + (1._R8P - weno_sigma), &
                                                 V=v, VR=vr_, WRATIO=wr)
#endif
         vr(1,f) = vr_(1)
         vr(2,f) = vr_(2)
         wmin = min(wmin, wr(1), wr(2))
      enddo
      call face_states(gamma=gamma, ch=ch, is_characteristic=is_characteristic, er=er, vr=vr, q0=qs(:,0), q1=qs(:,1), &
                       qL=qL, qR=qR)
      call riemann_flux(gamma=gamma, ch=ch, d=d, qL=qL, qR=qR, f=fr, fallback=fallback)
      if (fallback) fallbacks = fallbacks + 1_I4P
      !$acc loop seq
      do m=-1, 2
         call cell_flux(gamma=gamma, ch=ch, d=d, q=qs(:,m), qa=qas(:,m), f=fc(:,m))
      enddo
      sensor = 0._R8P
      if (wmin >= tau) sensor = 1._R8P
      !$acc loop seq
      do w=1, NV_K
         fl_gpu(b,i,j,k,w) = fr(w) + sensor * (cc(1) * fr(w) + cc(2) * (fc(w,0) + fc(w,1)) + &
                                               cc(3) * (fc(w,-1) + fc(w,2)) - fr(w))
      enddo
   enddo
   enddo
   enddo
   enddo
