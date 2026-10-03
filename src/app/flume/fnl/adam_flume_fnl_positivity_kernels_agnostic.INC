   ! Model-agnostic FNL device loops of the positivity limiter (issue #47, D-9, M3-P5b): the device twin of
   ! `adam_flume_cpu_positivity_kernels_agnostic.INC` (see it for the method), on the device layout `[nb, i, j, k, nv]`.
   !
   ! Included in the `contains` part of the device kernel modules of the models that accept the limiter (Euler, MHD
   ! without divergence control, MHD with EGLM), which provide by use association `NV_K`, `NV_AUX_K`, `IQ_R`, `S_MAX`,
   ! `POSITIVITY_LIMITER_EPS`, `POSITIVITY_LIMITER_KAPPA` and three device adapters (`!$acc routine seq`):
   ! * `backbone_face_flux(gamma, ch, d, qL, qR, f)`, `internal_energy(q)`: as on the CPU;
   ! * `cell_sources(hs, damping, ds, w, qsx, qsy, qsz, qac, s_lo, s_hi)`: the cell sources from the private stencils of
   !   the three directions (`qs*(:,m)`, `m = -hs ... hs`, gathered by the kernel) and the cell's auxiliary variables.

   subroutine compute_positivity_factors_dev(ni, nj, nk, ngc, blocks_number, gamma, ch, damping, hs, dt, dxyz_gpu, &
                                             is_null, q_gpu, q_aux_gpu, flx_gpu, fly_gpu, flz_gpu, lam_gpu, bad, &
                                             nonfinite)
   !< Compute the cell factor `Lambda` of every interior cell into `lam_gpu(:,:,:,:,1)` and count the cells whose backbone
   !< is not admissible and those with a non-finite high-order face flux (device twin of the CPU
   !< `compute_positivity_factors`).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                        !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                          !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                                  !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                     !< GLM cleaning speed (EGLM only).
   real(R8P),    intent(in)    :: damping                                !< GLM damping rate (EGLM only).
   integer(I4P), intent(in)    :: hs                                     !< Half stencil of the sources, 1 <= hs <= S_MAX.
   real(R8P),    intent(in)    :: dt                                     !< Forward-Euler step of the stage.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                        !< Blocks space steps [nb, 3].
   logical,      intent(in)    :: is_null(3)                             !< Null directions.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Conservative variables (ghosts updated).
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   real(R8P),    intent(in)    :: flx_gpu(1:,0:,1:,1:,1:)                !< X-face high-order fluxes.
   real(R8P),    intent(in)    :: fly_gpu(1:,1:,0:,1:,1:)                !< Y-face high-order fluxes.
   real(R8P),    intent(in)    :: flz_gpu(1:,1:,1:,0:,1:)                !< Z-face high-order fluxes.
   real(R8P),    intent(inout) :: lam_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Cell factors (component 1).
   integer(I4P), intent(out)   :: bad                                    !< Cells with an inadmissible backbone.
   integer(I4P), intent(out)   :: nonfinite                              !< Cells with a non-finite high-order face flux.
   integer(I8P), parameter     :: EXPONENT_MASK=2047_I8P                 !< The 11 exponent bits of a binary64.
   logical                     :: isnf                                   !< Private non-finite face difference flag.
   real(R8P)                   :: qsx(NV_K,-S_MAX:S_MAX)                 !< Private stencil along x.
   real(R8P)                   :: qsy(NV_K,-S_MAX:S_MAX)                 !< Private stencil along y.
   real(R8P)                   :: qsz(NV_K,-S_MAX:S_MAX)                 !< Private stencil along z.
   real(R8P)                   :: qac(NV_AUX_K)                          !< Private cell auxiliary variables.
   real(R8P)                   :: flo(NV_K)                              !< Private backbone face flux.
   real(R8P)                   :: dlt(NV_K,6)                            !< Private antidiffusive differences.
   real(R8P)                   :: qb(NV_K)                               !< Private backbone update.
   real(R8P)                   :: dsrc(NV_K)                             !< Private source difference.
   real(R8P)                   :: s_lo(NV_K), s_hi(NV_K)                 !< Private backbone and high-order sources.
   real(R8P)                   :: dc(NV_K)                               !< Private corner direction.
   real(R8P)                   :: qt(NV_K)                               !< Private corner state.
   real(R8P)                   :: ds(3), w(3)                            !< Private space steps, direction weights.
   real(R8P)                   :: rb, eb, er, ee                         !< Private backbone density, energy, floors.
   real(R8P)                   :: r1, e1, t, lc                          !< Private corner density, energy, factors.
   real(R8P)                   :: dtdx                                   !< Private dt / dx_d.
   integer(I4P)                :: nf                                     !< Private active faces.
   integer(I4P)                :: b, i, j, k, d, f, c, m, v              !< Counters.

   bad = 0_I4P
   nonfinite = 0_I4P
   !$acc parallel loop independent gang vector collapse(4) reduction(+:bad,nonfinite)                                &
   !$acc& DEVICEVAR(dxyz_gpu,q_gpu,q_aux_gpu,flx_gpu,fly_gpu,flz_gpu,lam_gpu)                                         &
   !$acc& firstprivate(ni,nj,nk,blocks_number,gamma,ch,damping,hs,dt,is_null)                                         &
   !$acc& private(qsx,qsy,qsz,qac,flo,dlt,qb,dsrc,s_lo,s_hi,dc,qt,ds,w,rb,eb,er,ee,r1,e1,t,lc,dtdx,nf,isnf)
   !$omp OMPLOOP collapse(4) reduction(+:bad,nonfinite) &
   !$omp& DEVICEPTR(dxyz_gpu,q_gpu,q_aux_gpu,flx_gpu,fly_gpu,flz_gpu,lam_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,gamma,ch,damping,hs,dt,is_null) &
   !$omp& private(qsx,qsy,qsz,qac,flo,dlt,qb,dsrc,s_lo,s_hi,dc,qt,ds,w,rb,eb,er,ee,r1,e1,t,lc,dtdx,nf,isnf)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      !$acc loop seq
      do m=-hs, hs
         !$acc loop seq
         do v=1, NV_K
            qsx(v,m) = q_gpu(b,i+m,j,k,v)
            qsy(v,m) = q_gpu(b,i,j+m,k,v)
            qsz(v,m) = q_gpu(b,i,j,k+m,v)
         enddo
      enddo
      !$acc loop seq
      do v=1, NV_AUX_K
         qac(v) = q_aux_gpu(b,i,j,k,v)
      enddo
      !$acc loop seq
      do d=1, 3
         ds(d) = dxyz_gpu(b,d)
         w(d)  = merge(0._R8P, 1._R8P, is_null(d))
      enddo
      !$acc loop seq
      do v=1, NV_K
         qb(v) = qsx(v,0)
      enddo
      nf = 0
      !$acc loop seq
      do d=1, 3
         if (is_null(d)) cycle
         dtdx = dt / ds(d)
         ! right face i+1/2
         select case(d)
         case(1)
            call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qsx(:,0), qR=qsx(:,1), f=flo)
         case(2)
            call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qsy(:,0), qR=qsy(:,1), f=flo)
         case(3)
            call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qsz(:,0), qR=qsz(:,1), f=flo)
         endselect
         nf = nf + 1
         !$acc loop seq
         do v=1, NV_K
            select case(d)
            case(1)
               dlt(v,nf) = -dtdx * (flx_gpu(b,i,j,k,v) - flo(v))
            case(2)
               dlt(v,nf) = -dtdx * (fly_gpu(b,i,j,k,v) - flo(v))
            case(3)
               dlt(v,nf) = -dtdx * (flz_gpu(b,i,j,k,v) - flo(v))
            endselect
            qb(v) = qb(v) - dtdx * flo(v)
         enddo
         ! left face i-1/2
         select case(d)
         case(1)
            call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qsx(:,-1), qR=qsx(:,0), f=flo)
         case(2)
            call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qsy(:,-1), qR=qsy(:,0), f=flo)
         case(3)
            call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qsz(:,-1), qR=qsz(:,0), f=flo)
         endselect
         nf = nf + 1
         !$acc loop seq
         do v=1, NV_K
            select case(d)
            case(1)
               dlt(v,nf) = dtdx * (flx_gpu(b,i-1,j,k,v) - flo(v))
            case(2)
               dlt(v,nf) = dtdx * (fly_gpu(b,i,j-1,k,v) - flo(v))
            case(3)
               dlt(v,nf) = dtdx * (flz_gpu(b,i,j,k-1,v) - flo(v))
            endselect
            qb(v) = qb(v) + dtdx * flo(v)
         enddo
      enddo
      call cell_sources(hs=hs, damping=damping, ds=ds, w=w, qsx=qsx, qsy=qsy, qsz=qsz, qac=qac, s_lo=s_lo, s_hi=s_hi)
      !$acc loop seq
      do v=1, NV_K
         qb(v)   = qb(v) + dt * s_lo(v)
         dsrc(v) = dt * (s_hi(v) - s_lo(v))
      enddo
      isnf = .false.
      !$acc loop seq
      do f=1, nf
         !$acc loop seq
         do v=1, NV_K
            if (iand(ishft(transfer(dlt(v,f), 0_I8P), -52), EXPONENT_MASK) == EXPONENT_MASK) isnf = .true.
         enddo
      enddo
      rb = qb(IQ_R)
      eb = -1._R8P
      if (rb > 0._R8P) eb = internal_energy(q=qb)
      if (.not.(rb > 0._R8P .and. eb > 0._R8P)) then
         bad = bad + 1_I4P
         lam_gpu(b,i,j,k,1) = 0._R8P
      elseif (isnf) then
         nonfinite = nonfinite + 1_I4P
         lam_gpu(b,i,j,k,1) = 0._R8P
      else
         ! the larger of the absolute and the relative floor: a stage may not take the cell below KAPPA times its backbone
         ! value, else a drained cell (rho ~ EPS, its energy kept) makes the next stage's backbone inadmissible
         er = max(min(POSITIVITY_LIMITER_EPS, 0.5_R8P * rb), POSITIVITY_LIMITER_KAPPA * rb)
         ee = max(min(POSITIVITY_LIMITER_EPS, 0.5_R8P * eb), POSITIVITY_LIMITER_KAPPA * eb)
         lc = 1._R8P
         !$acc loop seq
         do c=0, 2**nf - 1
            !$acc loop seq
            do v=1, NV_K
               dc(v) = dsrc(v)
               !$acc loop seq
               do f=1, nf
                  if (btest(c, f - 1)) dc(v) = dc(v) + dlt(v,f)
               enddo
            enddo
            t = 1._R8P
            r1 = rb + dc(IQ_R)
            if (r1 < er) t = (rb - er) / (rb - r1)
            !$acc loop seq
            do v=1, NV_K
               qt(v) = qb(v) + t * dc(v)
            enddo
            e1 = internal_energy(q=qt)
            if (e1 < ee) t = t * (eb - ee) / (eb - e1)
            lc = min(lc, t)
         enddo
         lam_gpu(b,i,j,k,1) = max(0._R8P, lc)
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_positivity_factors_dev

   subroutine blend_positivity_fluxes_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, q_gpu, lam_gpu, &
                                          fl_gpu, limited)
   !< Blend the face fluxes of direction `d` with the backbone, `theta` the smaller factor of the two cells; faces with
   !< `theta = 1` untouched, faces with `theta = 0` set to the backbone flux (device twin of the CPU
   !< `blend_positivity_fluxes`).
   integer(I4P), intent(in)    :: d                                    !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                           !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                      !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                        !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                                !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                   !< GLM cleaning speed (EGLM only).
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative variables.
   real(R8P),    intent(in)    :: lam_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Cell factors (component 1).
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)      !< Face fluxes of direction `d`.
   integer(I4P), intent(out)   :: limited                              !< Faces with theta < 1.
   real(R8P)                   :: qL(NV_K), qR(NV_K)                   !< Private cell states.
   real(R8P)                   :: flo(NV_K)                            !< Private backbone face flux.
   real(R8P)                   :: theta                                !< Private face factor.
   integer(I4P)                :: b, i, j, k, v                        !< Counters.

   limited = 0_I4P
   !$acc parallel loop independent gang vector collapse(4) reduction(+:limited) DEVICEVAR(q_gpu,lam_gpu,fl_gpu) &
   !$acc& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,gamma,ch) private(qL,qR,flo,theta)
   !$omp OMPLOOP collapse(4) reduction(+:limited) DEVICEPTR(q_gpu,lam_gpu,fl_gpu) &
   !$omp& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,gamma,ch) private(qL,qR,flo,theta)
   do k=1-dk, nk
   do j=1-dj, nj
   do i=1-di, ni
   do b=1, blocks_number
      theta = min(lam_gpu(b,i,j,k,1), lam_gpu(b,i+di,j+dj,k+dk,1), 1._R8P)
      if (theta < 1._R8P) then
         !$acc loop seq
         do v=1, NV_K
            qL(v) = q_gpu(b,i,j,k,v)
            qR(v) = q_gpu(b,i+di,j+dj,k+dk,v)
         enddo
         call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qL, qR=qR, f=flo)
         !$acc loop seq
         do v=1, NV_K
            if (theta > 0._R8P) then
               fl_gpu(b,i,j,k,v) = flo(v) + theta * (fl_gpu(b,i,j,k,v) - flo(v))
            else
               fl_gpu(b,i,j,k,v) = flo(v)
            endif
         enddo
         limited = limited + 1_I4P
      endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine blend_positivity_fluxes_dev
