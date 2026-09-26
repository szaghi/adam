   ! Model-agnostic CPU face-flux loop: the body shared by every per-model kernel module (issue #41, section 4).
   !
   ! Included in the `contains` part of `adam_flume_cpu_<model>_kernels`, which provides, by use association:
   ! `NV_K`, `NV_AUX_K` (the model's conservative and auxiliary variables numbers, compile-time constants that size
   ! every local array: no model branching and no maximum-model buffers in the loops), `S_MAX`, the model's
   ! split adapter `face_split_fluxes(gamma, ch, ...)` (one argument list for every model: the GLM cleaning speed `ch` is
   ! ignored by the models without cleaning), `compute_face_flux_back_projection`, and `weno_object`,
   ! `weno_reconstruct_upwind`. The host selects the model's module by `select case`, outside the loops.

   subroutine compute_face_fluxes(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, weno, &
                                  q, q_aux, fl)
   !< Compute the WENO face fluxes of direction `d`: face `(i,j,k)` lies between cells `(i,j,k)` and
   !< `(i+di,j+dj,k+dk)`, so the face array starts at index 0 along `d`.
   !<
   !< Per face: gather the stencil `m = 1-S ... S` into constant-bound locals, project and split it, reconstruct every
   !< field (the local `v` is the packed stencil of the WENO primitive), back-project.
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
   real(R8P)                        :: v(2,2*S_MAX-1)                        !< Packed stencil of one field.
   real(R8P)                        :: vr(2,NV_K)                            !< Reconstructed split fields.
   integer(I4P)                     :: S                                     !< WENO stencil half-width.
   integer(I4P)                     :: b, i, j, k, m, f                      !< Counters.

   S = weno%S
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
                                      qas=qas, fsplit=fsplit, er=er)
               do f=1, NV_K
                  do m=1-S, S-1
                     v(:,m+S) = fsplit(:,m,f)
                  enddo
                  call weno_reconstruct_upwind(S=S, weno_a=weno%a, weno_p=weno%p, weno_d=weno%d, weno_zeps=weno%zeps, &
                                               v=v, vr=vr(:,f))
               enddo
               call compute_face_flux_back_projection(is_characteristic=is_characteristic, er=er, vr=vr, &
                                                      flux=fl(1:NV_K,i,j,k,b))
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_face_fluxes
