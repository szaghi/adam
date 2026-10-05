   ! Model-agnostic FNL device loops of the positivity limiter (issue #47, D-9, M3-P5b): the device twin of
   ! `adam_flume_cpu_positivity_kernels_agnostic.INC` (see it for the method), on the device layout `[nb, i, j, k, nv]`.
   !
   ! Included in the `contains` part of the device kernel modules of the models that accept the limiter (Euler, MHD
   ! without divergence control, MHD with EGLM), which provide by use association `NV_K`, `NV_AUX_K`, `IQ_R`, `S_MAX`,
   ! `POSITIVITY_LIMITER_KAPPA` and three device adapters (`!$acc routine seq`):
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
         ! the relative floor: a stage may not take the cell below KAPPA times its backbone value (no absolute floor: issue
         ! #49, the scheme has no unit), else a drained cell (its energy kept) makes the next stage's backbone inadmissible
         er = POSITIVITY_LIMITER_KAPPA * rb
         ee = POSITIVITY_LIMITER_KAPPA * eb
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

   subroutine compute_seam_positivity_factors_host(ncells, hs, gamma, ch, damping, dt, ds, is_null, qs, qa, fh, omask, &
                                                   olo, ohi, lam, dbad, dnonfinite)
   !< Recompute on the host the factor of the cells beside a 2:1 seam with the seam faces' backbone and high-order
   !< fluxes overridden (issue #50, D1; the host twin of the CPU `compute_seam_positivity_factors`, fed by the device
   !< gathers of the seam cells: their axis stencils `qs(:, l + hs + 1, d, n)`, auxiliary variables `qa(:, n)` and face
   !< fluxes `fh(:, fec, n)`, `fec` 1..6 = -x, +x, -y, +y, -z, +z). `dbad`, `dnonfinite`: the change of the counts
   !< against the factors without the overrides, which the device kernel has already counted.
   integer(I4P), intent(in)  :: ncells                     !< Seam cells.
   integer(I4P), intent(in)  :: hs                         !< Stencil half width, 1 <= hs <= S_MAX.
   real(R8P),    intent(in)  :: gamma                      !< Specific heats ratio.
   real(R8P),    intent(in)  :: ch                         !< GLM cleaning speed (EGLM only).
   real(R8P),    intent(in)  :: damping                    !< GLM damping rate (EGLM only).
   real(R8P),    intent(in)  :: dt                         !< Forward-Euler step of the stage.
   real(R8P),    intent(in)  :: ds(1:,1:)                  !< Space steps of the cells' blocks [3, ncells].
   logical,      intent(in)  :: is_null(3)                 !< Null directions.
   real(R8P),    intent(in)  :: qs(1:,1:,1:,1:)            !< Axis stencils [nv, 2 hs + 1, 3, ncells].
   real(R8P),    intent(in)  :: qa(1:,1:)                  !< Auxiliary variables [nv_aux, ncells].
   real(R8P),    intent(in)  :: fh(1:,1:,1:)               !< High-order face fluxes [nv, 6, ncells].
   logical,      intent(in)  :: omask(1:,1:)               !< Overridden faces [6, ncells].
   real(R8P),    intent(in)  :: olo(1:,1:,1:)              !< Override backbone fluxes [nv, 6, ncells].
   real(R8P),    intent(in)  :: ohi(1:,1:,1:)              !< Override high-order fluxes [nv, 6, ncells].
   real(R8P),    intent(out) :: lam(1:)                    !< Cell factors [ncells].
   integer(I4P), intent(out) :: dbad                       !< Change of the inadmissible backbones count.
   integer(I4P), intent(out) :: dnonfinite                 !< Change of the non-finite cells count.
   real(R8P)                 :: qsx(NV_K,-S_MAX:S_MAX)     !< Stencil along x.
   real(R8P)                 :: qsy(NV_K,-S_MAX:S_MAX)     !< Stencil along y.
   real(R8P)                 :: qsz(NV_K,-S_MAX:S_MAX)     !< Stencil along z.
   real(R8P)                 :: qac(NV_AUX_K)              !< Cell auxiliary variables.
   real(R8P)                 :: qc(NV_K), qn(NV_K)         !< Cell and neighbour states.
   real(R8P)                 :: flo(NV_K), fhi(NV_K)       !< Backbone and high-order face fluxes.
   real(R8P)                 :: dlt(NV_K,6), dlt0(NV_K,6)  !< Antidiffusive differences, with/without overrides.
   real(R8P)                 :: qb(NV_K), qb0(NV_K)        !< Backbone updates, with/without overrides.
   real(R8P)                 :: dsrc(NV_K)                 !< Source difference.
   real(R8P)                 :: s_lo(NV_K), s_hi(NV_K)     !< Backbone and high-order sources.
   real(R8P)                 :: w(3)                       !< Direction weights.
   real(R8P)                 :: dtdx, sg, lc, lc0          !< dt / dx_d, face sign; factors with/without overrides.
   integer(I4P)              :: nf, st, st0                !< Active faces; statuses with/without overrides.
   integer(I4P)              :: n, d, s, fec, l            !< Counters.

   dbad = 0_I4P
   dnonfinite = 0_I4P
   do d=1, 3
      w(d) = merge(0._R8P, 1._R8P, is_null(d))
   enddo
   do n=1, ncells
      qsx = 0._R8P ; qsy = 0._R8P ; qsz = 0._R8P
      do l=-hs, hs
         qsx(:,l) = qs(1:NV_K,l+hs+1,1,n)
         qsy(:,l) = qs(1:NV_K,l+hs+1,2,n)
         qsz(:,l) = qs(1:NV_K,l+hs+1,3,n)
      enddo
      qac = qa(1:NV_AUX_K,n)
      qc = qsx(:,0)
      qb = qc
      qb0 = qc
      nf = 0
      do d=1, 3
         if (is_null(d)) cycle
         dtdx = dt / ds(d,n)
         do s=1, 0, -1 ! the right face (s = 1, fec = 2d) first, then the left one, as the device kernel
            fec = 2 * d - 1 + s
            select case(d)
            case(1)
               qn = merge(qsx(:,1), qsx(:,-1), s == 1)
            case(2)
               qn = merge(qsy(:,1), qsy(:,-1), s == 1)
            case default
               qn = merge(qsz(:,1), qsz(:,-1), s == 1)
            endselect
            if (s == 1) then
               call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qc, qR=qn, f=flo)
               sg = -dtdx
            else
               call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qn, qR=qc, f=flo)
               sg = dtdx
            endif
            fhi = fh(1:NV_K,fec,n)
            nf = nf + 1
            dlt0(:,nf) = sg * (fhi - flo)
            qb0 = qb0 + sg * flo
            if (omask(fec,n)) then
               flo = olo(1:NV_K,fec,n)
               fhi = ohi(1:NV_K,fec,n)
            endif
            dlt(:,nf) = sg * (fhi - flo)
            qb = qb + sg * flo
         enddo
      enddo
      call cell_sources(hs=hs, damping=damping, ds=ds(:,n), w=w, qsx=qsx, qsy=qsy, qsz=qsz, qac=qac, s_lo=s_lo, &
                        s_hi=s_hi)
      qb = qb + dt * s_lo
      qb0 = qb0 + dt * s_lo
      dsrc = dt * (s_hi - s_lo)
      call positivity_cell_factor(nf=nf, qb=qb0, dlt=dlt0, dsrc=dsrc, lc=lc0, status=st0)
      call positivity_cell_factor(nf=nf, qb=qb, dlt=dlt, dsrc=dsrc, lc=lc, status=st)
      dbad = dbad + merge(1_I4P, 0_I4P, st == 1_I4P) - merge(1_I4P, 0_I4P, st0 == 1_I4P)
      dnonfinite = dnonfinite + merge(1_I4P, 0_I4P, st == 2_I4P) - merge(1_I4P, 0_I4P, st0 == 2_I4P)
      lam(n) = lc
   enddo
   endsubroutine compute_seam_positivity_factors_host

   subroutine compute_backbone_fluxes_host(n, d, gamma, ch, qL, qR, flo)
   !< Compute on the host the backbone flux of `n` face state pairs of direction `d` (issue #50: the fine faces of a 2:1
   !< seam take the coarse donor cell as their outer state).
   integer(I4P), intent(in)  :: n            !< Faces number.
   integer(I4P), intent(in)  :: d            !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: gamma        !< Specific heats ratio.
   real(R8P),    intent(in)  :: ch           !< GLM cleaning speed (EGLM only).
   real(R8P),    intent(in)  :: qL(1:,1:)    !< Left states [nv, n].
   real(R8P),    intent(in)  :: qR(1:,1:)    !< Right states [nv, n].
   real(R8P),    intent(out) :: flo(1:,1:)   !< Backbone fluxes [nv, n].
   real(R8P)                 :: f(NV_K)      !< Face flux.
   integer(I4P)              :: m            !< Counter.

   do m=1, n
      call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qL(1:NV_K,m), qR=qR(1:NV_K,m), f=f)
      flo(1:NV_K,m) = f
   enddo
   endsubroutine compute_backbone_fluxes_host

   pure subroutine positivity_cell_factor(nf, qb, dlt, dsrc, lc, status)
   !< The factor of one cell from its backbone update, antidiffusive differences and source difference (the decision
   !< of `compute_positivity_factors_dev`, for the host seam pass): `status` 1 inadmissible backbone, 2 non-finite
   !< difference (`lc = 0` both), else 0 with `lc` the largest factor of the corner box.
   integer(I4P), intent(in)  :: nf                      !< Active faces.
   real(R8P),    intent(in)  :: qb(NV_K)                !< Backbone update.
   real(R8P),    intent(in)  :: dlt(NV_K,6)             !< Antidiffusive differences.
   real(R8P),    intent(in)  :: dsrc(NV_K)              !< Source difference.
   real(R8P),    intent(out) :: lc                      !< Cell factor.
   integer(I4P), intent(out) :: status                  !< 0 admissible, 1 inadmissible backbone, 2 non-finite.
   integer(I8P), parameter   :: EXPONENT_MASK=2047_I8P  !< The 11 exponent bits of a binary64.
   real(R8P)                 :: dc(NV_K)                !< Corner direction.
   real(R8P)                 :: rb, eb, er, ee          !< Backbone density, energy, their floors.
   real(R8P)                 :: r1, e1, t               !< Corner density, energy, factor.
   integer(I4P)              :: f, c                    !< Counters.

   lc = 0._R8P
   rb = qb(IQ_R)
   eb = -1._R8P
   if (rb > 0._R8P) eb = internal_energy(q=qb)
   if (.not.(rb > 0._R8P .and. eb > 0._R8P)) then
      status = 1_I4P
      return
   endif
   do f=1, nf
      do c=1, NV_K
         if (iand(ishft(transfer(dlt(c,f), 0_I8P), -52), EXPONENT_MASK) == EXPONENT_MASK) then
            status = 2_I4P
            return
         endif
      enddo
   enddo
   status = 0_I4P
   er = POSITIVITY_LIMITER_KAPPA * rb
   ee = POSITIVITY_LIMITER_KAPPA * eb
   lc = 1._R8P
   do c=0, 2**nf - 1
      dc = dsrc
      do f=1, nf
         if (btest(c, f - 1)) dc = dc + dlt(:,f)
      enddo
      t = 1._R8P
      r1 = rb + dc(IQ_R)
      if (r1 < er) t = (rb - er) / (rb - r1)
      e1 = internal_energy(q=qb + t * dc)
      if (e1 < ee) t = t * (eb - ee) / (eb - e1)
      lc = min(lc, t)
   enddo
   lc = max(0._R8P, lc)
   endsubroutine positivity_cell_factor
