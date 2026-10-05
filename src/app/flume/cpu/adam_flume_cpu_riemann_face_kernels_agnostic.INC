   ! Model- and solver-agnostic CPU face-flux loop of the scheme `weno-riemann` (issue #47): the body shared by every
   ! per-(model, solver) kernel module.
   !
   ! Included in the `contains` part of `adam_flume_cpu_<model>_<solver>_kernels`, which provides, by use association:
   ! `NV_K`, `NV_AUX_K` (compile-time constants sizing every local array), `S_MAX`, `weno_object`,
   ! `weno_reconstruct_upwind_wratio`, and four adapters with one argument list for every model (the GLM cleaning speed
   ! `ch` is ignored by the models without cleaning): `face_interpolation_fields` (the stencil fields to interpolate, in
   ! the WENO upwind layout), `face_states` (the conservative face states, with the first-order fallback),
   ! `riemann_flux` (the solver, with a flag set when it fell back to a simpler solver, e.g. HLLD to HLL) and `cell_flux`
! (the physical flux of a cell). The host selects the module by
   ! `select case` on the model and the solver, outside the loops.

   subroutine compute_riemann_face_fluxes(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, is_characteristic, &
                                          weno, cc, tau, q, q_aux, fl, fallbacks)
   !< Compute the `weno-riemann` face fluxes of direction `d`: face `(i,j,k)` lies between cells `(i,j,k)` and
   !< `(i+di,j+dj,k+dk)`, so the face array starts at index 0 along `d`.
   !<
   !< Per face: gather the stencil `m = min(1-S,-1) ... max(S,2)`, WENO-interpolate its fields to the two face states
   !< (interpolation tables `weno%a_interp`, `weno%p_interp`), take the Riemann flux `F` of the states and add the
   !< correction `F^ = F + s (cc(1) F + cc(2) (f_0 + f_1) + cc(3) (f_-1 + f_2) - F)`, with `f_m` the physical flux of
   !< cell `m` and the sensor `s = 1` where the WENO weights ratio of every field on both sides is `>= tau`, else 0;
   !< `fallbacks` counts the faces whose Riemann solver fell back (HLLD to HLL; always 0 for the other solvers).
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
   real(R8P)                        :: v(2,2*S_MAX-1)                    !< Packed stencil of one field.
   real(R8P)                        :: vr(2,NV_K)                        !< Interpolated fields.
   real(R8P)                        :: wr(2)                             !< WENO weights ratio of one field.
   real(R8P)                        :: wmin                              !< Minimum WENO weights ratio of the face.
   real(R8P)                        :: qL(NV_K), qR(NV_K)                !< Face states.
   real(R8P)                        :: fr(NV_K)                          !< Riemann flux.
   real(R8P)                        :: fc(NV_K,-1:2)                     !< Physical fluxes of the cells -1 ... 2.
   real(R8P)                        :: sensor                            !< Correction switch, 0 or 1.
   logical                          :: fallback                          !< Riemann solver fallback flag of a face.
   integer(I4P)                     :: S                                 !< WENO stencil half-width.
   integer(I4P)                     :: b, i, j, k, m, f                  !< Counters.

   S = weno%S
   fallbacks = 0_I4P
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
   endsubroutine compute_riemann_face_fluxes
