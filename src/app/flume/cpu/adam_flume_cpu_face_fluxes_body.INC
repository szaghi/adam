   ! Body of the FLUME CPU routine `compute_face_fluxes` for ONE WENO stencil width, S a named constant: included by
   ! `compute_face_fluxes_s1` ... `compute_face_fluxes_s5` of every per-model kernel module (through
   ! `adam_flume_cpu_face_kernels_agnostic.INC`), after their declaration of S (issue #139, S3a: the CPU dispatch
   ! structure of the FNL kernels, decision D3).
   integer(I4P),      intent(in)    :: d                                     !< Direction, 1=x, 2=y, 3=z.
   integer(I4P),      intent(in)    :: di, dj, dk                            !< Unit step along `d`.
   integer(I4P),      intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P),      intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),         intent(in)    :: gamma                                 !< Specific heats ratio.
   real(R8P),         intent(in)    :: ch                                    !< GLM cleaning speed (models with GLM only).
   logical,           intent(in)    :: is_characteristic                     !< Characteristic (or conservative) variables.
   type(weno_object), intent(in)    :: weno                                  !< WENO coefficients.
   real(R8P),         intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)         !< Conservative variables.
   real(R8P),         intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Auxiliary variables.
   real(R8P),         intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)           !< Face fluxes of direction `d`.
   real(R8P)                        :: qs(NV_K,1-S_MAX:S_MAX)                !< Stencil conservative variables.
   real(R8P)                        :: qas(NV_AUX_K,1-S_MAX:S_MAX)           !< Stencil auxiliary variables.
   real(R8P)                        :: fsplit(2,1-S_MAX:S_MAX-1,NV_K)        !< Split fields.
   real(R8P)                        :: er(NV_K,NV_K)                         !< Right eigenvectors.
   real(R8P)                        :: mu(NV_K)                              !< Magnitude of each split field (WENO descaler).
   real(R8P)                        :: v(2,2*S_MAX-1)                        !< Packed stencil of one field.
   real(R8P)                        :: vr(2,NV_K)                            !< Reconstructed split fields.
   integer(I4P)                     :: b, i, j, k, m, f                      !< Counters.

   !$omp parallel do collapse(4) default(firstprivate) shared(weno, q, q_aux, fl)
   do b=1, blocks_number
      do k=1-dk, nk
         do j=1-dj, nj
            do i=1-di, ni
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
               call compute_face_flux_back_projection(is_characteristic=is_characteristic, er=er, vr=vr, &
                                                      flux=fl(1:NV_K,i,j,k,b))
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
