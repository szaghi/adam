   ! Body of the FLUME CPU routine `add_gravity_edge_correction` for ONE WENO stencil width, S a named constant:
   ! included by `add_gravity_edge_correction_s1` ... `add_gravity_edge_correction_s5` of every per-model kernel
   ! module (through `adam_flume_cpu_face_kernels_agnostic.INC`), after their declaration of S (issue #139, S3a: the
   ! CPU dispatch structure of the FNL kernels, decision D3).
   integer(I4P),      intent(in)    :: d                                     !< Direction, 1=x, 2=y, 3=z.
   integer(I4P),      intent(in)    :: di, dj, dk                            !< Unit step along `d`.
   integer(I4P),      intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P),      intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),         intent(in)    :: gamma                                 !< Specific heats ratio.
   real(R8P),         intent(in)    :: ch                                    !< GLM cleaning speed (models with GLM only).
   logical,           intent(in)    :: is_characteristic                     !< Characteristic (or conservative) variables.
   type(weno_object), intent(in)    :: weno                                  !< WENO coefficients.
   integer(I4P),      intent(in)    :: ir                                    !< Index of the density (mass) variable.
   real(R8P),         intent(in)    :: gd                                    !< Gravity along `d`.
   real(R8P),         intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)         !< Conservative variables.
   real(R8P),         intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Auxiliary variables.
   real(R8P),         intent(inout) :: gc(1:,1-ngc:,1-ngc:,1-ngc:,1:)        !< Gravity energy correction [1, ...].
   real(R8P)                        :: qs(NV_K,1-S_MAX:S_MAX)                !< Stencil conservative variables.
   real(R8P)                        :: qas(NV_AUX_K,1-S_MAX:S_MAX)           !< Stencil auxiliary variables.
   real(R8P)                        :: fsplit(2,1-S_MAX:S_MAX-1,NV_K)        !< Split fields.
   real(R8P)                        :: er(NV_K,NV_K)                         !< Right eigenvectors.
   real(R8P)                        :: mu(NV_K)                              !< Magnitude of each split field (WENO descaler).
   real(R8P)                        :: v(2,2*S_MAX-1)                        !< Packed stencil of one field.
   real(R8P)                        :: vr(2,NV_K)                            !< Reconstructed split fields.
   real(R8P)                        :: flux(NV_K)                            !< Face flux.
   integer(I4P)                     :: nd, n1, n2                            !< Cells along `d` and across it.
   integer(I4P)                     :: fd, cd                                !< Face and cell index along `d`.
   integer(I4P)                     :: b, i, j, k, m, f, t1, t2, side        !< Counters.

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
               do m=1-S, S
                  qs(:,m)  = q(1:NV_K,i+m*di,j+m*dj,k+m*dk,b)
                  qas(:,m) = q_aux(1:NV_AUX_K,i+m*di,j+m*dj,k+m*dk,b)
               enddo
               call face_split_fluxes(gamma=gamma, ch=ch, d=d, S=S, is_characteristic=is_characteristic, qs=qs, &
                                      qas=qas, fsplit=fsplit, er=er, mu=mu)
               do f=1, NV_K
                  do m=1-S, S-1
                     v(:,m+S) = fsplit(:,m,f)
                  enddo
                  call weno_reconstruct_upwind(S=S, weno_a=weno%a, weno_p=weno%p, weno_d=weno%d, weno_zeps=weno%zeps, &
                                               weno_rmu=weno%sigma / mu(f) + (1._R8P - weno%sigma), v=v, vr=vr(:,f))
               enddo
               call compute_face_flux_back_projection(is_characteristic=is_characteristic, er=er, vr=vr, flux=flux)
               i = di * cd + (1 - di) * t1
               j = dj * cd + di * t1 + dk * t2
               k = dk * cd + (1 - dk) * t2
               gc(1,i,j,k,b) = gc(1,i,j,k,b) - gd * flux(ir) / 24._R8P
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
