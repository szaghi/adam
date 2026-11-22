   ! Body of the FLUME CPU routine `add_gravity_edge_correction` for ONE WENO stencil width, S a named constant:
   ! included by `add_gravity_edge_correction_s1` ... `add_gravity_edge_correction_s5` of every per-model kernel
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
   integer(I4P),      intent(in)    :: ir                                !< Index of the density (mass) variable.
   real(R8P),         intent(in)    :: gd                                !< Gravity along `d`.
   real(R8P),         intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Conservative variables.
   real(R8P),         intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),         intent(inout) :: gc(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Gravity energy correction [1, ...].
   real(R8P)                        :: qs(NV_K,1-S_MAX:S_MAX)            !< Stencil conservative variables.
   real(R8P)                        :: qas(NV_AUX_K,1-S_MAX:S_MAX)       !< Stencil auxiliary variables.
   real(R8P)                        :: fint(2,1-S_MAX:S_MAX-1,NV_K)      !< Fields to interpolate.
   real(R8P)                        :: er(NV_K,NV_K)                     !< Right eigenvectors.
   real(R8P)                        :: mu(NV_K)                          !< Magnitude of each field (WENO descaler).
   real(R8P)                        :: v(2,2*S_MAX-1)                    !< Packed stencil of one field.
   real(R8P)                        :: vr(2,NV_K)                        !< Interpolated fields.
   real(R8P)                        :: wr(2)                             !< WENO weights ratio of one field.
   real(R8P)                        :: wmin                              !< Minimum WENO weights ratio of the face.
   real(R8P)                        :: qL(NV_K), qR(NV_K)                !< Face states.
   real(R8P)                        :: fr(NV_K)                          !< Riemann flux.
   real(R8P)                        :: fc(NV_K,-1:2)                     !< Physical fluxes of the cells -1 ... 2.
   real(R8P)                        :: sensor                            !< Correction switch, 0 or 1.
   real(R8P)                        :: fm                                !< Mass flux of the face.
   logical                          :: fallback                          !< Riemann solver fallback flag of a face.
   integer(I4P)                     :: nd, n1, n2                        !< Cells along `d` and across it.
   integer(I4P)                     :: fd, cd                            !< Face and cell index along `d`.
   integer(I4P)                     :: b, i, j, k, m, f, t1, t2, side    !< Counters.

   nd = di * ni + dj * nj + dk * nk
   n1 = merge(nj, ni, di == 1_I4P)
   n2 = merge(nj, nk, dk == 1_I4P)
   !$omp parallel do collapse(3) default(firstprivate) shared(weno, q, q_aux, gc)
   do b=1, blocks_number
      do t2=1, n2
         do t1=1, n1
            do side=1, 2
               fd = merge(-1_I4P, nd + 1_I4P, side == 1)
               cd = merge(1_I4P, nd, side == 1)
               i = di * fd + (1 - di) * t1
               j = dj * fd + di * t1 + dk * t2
               k = dk * fd + (1 - dk) * t2
               do m=min(1-S,-1), max(S,2)
                  qs(:,m)  = q(1:NV_K,i+m*di,j+m*dj,k+m*dk,b)
                  qas(:,m) = q_aux(1:NV_AUX_K,i+m*di,j+m*dj,k+m*dk,b)
               enddo
               call face_interpolation_fields(gamma=gamma, ch=ch, d=d, S=S, is_characteristic=is_characteristic, &
                                              qs=qs, qas=qas, fint=fint, er=er, mu=mu)
               wmin = 1._R8P
               do f=1, NV_K
                  do m=1-S, S-1
                     v(:,m+S) = fint(:,m,f)
                  enddo
                  call weno_reconstruct_upwind_wratio(S=S, weno_a=weno%a_interp, weno_p=weno%p_interp, weno_d=weno%d, &
                                                      weno_zeps=weno%zeps, weno_rmu=weno%sigma / mu(f) + (1._R8P - weno%sigma), &
                                                      v=v, vr=vr(:,f), wr=wr)
                  wmin = min(wmin, wr(1), wr(2))
               enddo
               call face_states(gamma=gamma, ch=ch, is_characteristic=is_characteristic, er=er, vr=vr, q0=qs(:,0), &
                                q1=qs(:,1), qL=qL, qR=qR)
               call riemann_flux(gamma=gamma, ch=ch, d=d, qL=qL, qR=qR, f=fr, fallback=fallback)
               do m=-1, 2
                  call cell_flux(gamma=gamma, ch=ch, d=d, q=qs(:,m), qa=qas(:,m), f=fc(:,m))
               enddo
               sensor = merge(1._R8P, 0._R8P, wmin >= tau)
               fm = fr(ir) + sensor * (cc(1) * fr(ir) + cc(2) * (fc(ir,0) + fc(ir,1)) + &
                                       cc(3) * (fc(ir,-1) + fc(ir,2)) - fr(ir))
               i = di * cd + (1 - di) * t1
               j = dj * cd + di * t1 + dk * t2
               k = dk * cd + (1 - dk) * t2
               gc(1,i,j,k,b) = gc(1,i,j,k,b) - gd * fm / 24._R8P
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
