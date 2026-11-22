   ! Body of the FLUME CPU routine `compute_riemann_face_fluxes` for ONE WENO stencil width, S a named constant:
   ! included by `compute_riemann_face_fluxes_s1` ... `compute_riemann_face_fluxes_s5` of every per-model kernel
   ! module (through `adam_flume_cpu_riemann_face_kernels_agnostic.INC`), after their declaration of S (issue #139,
   ! S3a: the CPU dispatch structure of the FNL kernels, decision D3).
   integer(I4P),      intent(in)    :: d                                 !< Direction, 1=x, 2=y, 3=z.
   integer(I4P),      intent(in)    :: di, dj, dk                        !< Unit step along `d`.
   integer(I4P),      intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P),      intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),         intent(in)    :: gamma                             !< Specific heats ratio.
   real(R8P),         intent(in)    :: ch                                !< GLM cleaning speed (models with GLM only).
   logical,           intent(in)    :: is_characteristic                 !< Characteristic (or primitive) variables.
   type(weno_object), intent(in)    :: weno                              !< WENO coefficients (with interpolation tables).
   real(R8P),         intent(in)    :: cc(3)                             !< Correction coefficients.
   real(R8P),         intent(in)    :: tau                               !< Sensor threshold.
   real(R8P),         intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   real(R8P),         intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),         intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)       !< Face fluxes of direction `d`.
   integer(I4P),      intent(out)   :: fallbacks                         !< Faces whose Riemann solver fell back.
   real(R8P)                        :: qs(NV_K,1-S_MAX:S_MAX)            !< Stencil conservative variables.
   real(R8P)                        :: qas(NV_AUX_K,1-S_MAX:S_MAX)       !< Stencil auxiliary variables.
   real(R8P)                        :: fint(2,1-S_MAX:S_MAX-1,NV_K)      !< Fields to interpolate.
   real(R8P)                        :: er(NV_K,NV_K)                     !< Right eigenvectors.
   real(R8P)                        :: mu(NV_K)                          !< Magnitude of each field (WENO descaler).
   real(R8P)                        :: vp(1-S:S), vm(1-S:S)              !< Stencil of one field, left and right bias.
   real(R8P)                        :: vr(2,NV_K)                        !< Interpolated fields.
   real(R8P)                        :: wr(2)                             !< WENO weights ratio of one field.
   real(R8P)                        :: wmin                              !< Minimum WENO weights ratio of the face.
   real(R8P)                        :: qL(NV_K), qR(NV_K)                !< Face states.
   real(R8P)                        :: fr(NV_K)                          !< Riemann flux.
   real(R8P)                        :: fc(NV_K,-1:2)                     !< Physical fluxes of the cells -1 ... 2.
   real(R8P)                        :: sensor                            !< Correction switch, 0 or 1.
   logical                          :: fallback                          !< Riemann solver fallback flag of a face.
   integer(I4P)                     :: b, i, j, k, m, f                  !< Counters.

   fallbacks = 0_I4P
   vp(S) = 0._R8P ; vm(1-S) = 0._R8P ! not read by the kernel (SETH's stencil convention)
   !$omp parallel do collapse(4) default(firstprivate) shared(weno, q, q_aux, fl) reduction(+:fallbacks)
   do b=1, blocks_number
      do k=1-dk, nk
         do j=1-dj, nj
            do i=1-di, ni
               do m=min(1-S,-1), max(S,2)
                  qs(:,m)  = q(1:NV_K,i+m*di,j+m*dj,k+m*dk,b)
                  qas(:,m) = q_aux(1:NV_AUX_K,i+m*di,j+m*dj,k+m*dk,b)
               enddo
               call face_interpolation_fields(gamma=gamma, ch=ch, d=d, S=S, is_characteristic=is_characteristic, &
                                              qs=qs, qas=qas, fint=fint, er=er, mu=mu)
               wmin = 1._R8P
               do f=1, NV_K
                  do m=1-S, S-1
                     vp(m)   = fint(2,m,f)
                     vm(m+1) = fint(1,m,f)
                  enddo
                  call weno_width(vp=vp, vm=vm, eps=weno%eps, rmu=weno%sigma / mu(f) + (1._R8P - weno%sigma), &
                                  rp=vr(2,f), rm=vr(1,f), wrp=wr(2), wrm=wr(1))
                  wmin = min(wmin, wr(1), wr(2))
               enddo
               call face_states(gamma=gamma, ch=ch, is_characteristic=is_characteristic, er=er, vr=vr, q0=qs(:,0), &
                                q1=qs(:,1), qL=qL, qR=qR)
               call riemann_flux(gamma=gamma, ch=ch, d=d, qL=qL, qR=qR, f=fr, fallback=fallback)
               if (fallback) fallbacks = fallbacks + 1_I4P
               do m=-1, 2
                  call cell_flux(gamma=gamma, ch=ch, d=d, q=qs(:,m), qa=qas(:,m), f=fc(:,m))
               enddo
               sensor = merge(1._R8P, 0._R8P, wmin >= tau)
               do f=1, NV_K
                  fl(f,i,j,k,b) = fr(f) + sensor * (cc(1) * fr(f) + cc(2) * (fc(f,0) + fc(f,1)) + &
                                                    cc(3) * (fc(f,-1) + fc(f,2)) - fr(f))
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
