   ! Model-agnostic CPU loops of the positivity limiter (issue #47, D-9, M3-P5a).
   !
   ! Included in the `contains` part of the kernel modules of the models that accept the limiter (Euler, MHD without
   ! divergence control, MHD with EGLM; D-10 refuses it with the mixed GLM), which provide by use association `NV_K`,
   ! `NV_AUX_K`, `IQ_R`, `POSITIVITY_LIMITER_KAPPA`, and three private adapters:
   ! * `backbone_face_flux(gamma, ch, d, qL, qR, f)`: the first-order Lax-Friedrichs flux of two cell states, whose
   !   forward-Euler update is admissible under the CFL condition (Euler: Rusanov; MHD: the Wu speed);
   ! * `internal_energy(q)`: the internal energy per unit volume, `E` minus the kinetic, magnetic and cleaning energies,
   !   a concave function of `q` (positive with the density iff the state is admissible);
   ! * `cell_sources(...)`: the cell sources of the backbone (`s_lo`) and of the high-order update (`s_hi`); zero but for
   !   EGLM (second-order and order-2S centred differences, plus the damping in both).
   !
   ! The limiter (Christlieb et al. 2015, doi:10.1137/140971208, the cell-based form of Xu 2014,
   ! doi:10.1090/S0025-5718-2013-02788-3): per cell, the backbone `q_L = q - dt/dx (F^LLF_{i+1/2} - F^LLF_{i-1/2}) (each
   ! active direction) + dt s_lo` and, per face, the antidiffusive difference `-+dt/dx (F - F^LLF)` of the high-order
   ! flux `F`; the cell factor `Lambda` is the largest value in [0, 1] for which every corner of the box `[0,
   ! Lambda]^faces` (all subsets of the faces at `Lambda`, the source difference `dt (s_hi - s_lo)` always at `Lambda`)
   ! keeps `rho >= eps` and `rho e >= eps`; as both are concave in `q` (the density linear), each corner's factor is a
   ! closed form (Zhang & Shu 2012, doi:10.1016/j.jcp.2011.11.020): first the density along the segment, then the
   ! internal energy through its chord. `eps = POSITIVITY_LIMITER_KAPPA backbone value`: positivity alone keeps the
   ! state admissible but not usable (the planar Sedov blast on HLLC drained its centre cell to `rho = 1e-13` in one stage
   ! under the former absolute floor, the sound speed grew by 1e6 and the next stage's backbone, whose time step is that
   ! of the step's first state, was inadmissible), so a stage may not take a cell below `KAPPA` times its backbone value;
   ! smooth updates differ from the backbone by O(dx), far from it. The floor is relative only (issue #49: FLUME's
   ! equations carry no unit, an absolute floor would make the result depend on the units of the input). A face
   ! takes `theta = min(Lambda)` of its two cells and `F <- F^LLF + theta (F - F^LLF)`, so each forward-Euler step of
   ! size `dt` is admissible; an SSP Runge-Kutta stage is a convex combination of such steps (all of size at most `dt`).

   subroutine compute_positivity_factors(ni, nj, nk, ngc, blocks_number, gamma, ch, damping, hs, dt, dxyz, is_null, q, &
                                         q_aux, flx, fly, flz, lam, bad, nonfinite)
   !< Compute the cell factor `Lambda` of every interior cell into `lam(1,:,:,:,:)` (the other components untouched) and
   !< count the cells whose backbone is not admissible (`Lambda = 0` there: the floors, if enabled, take over) and those
   !< with a non-finite high-order face flux (`Lambda = 0` too: a NaN or infinite flux carries no information, and every
   !< comparison with it is false, so the corner factors would keep it; tested on the bit pattern as `count_nonfinite`).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                              !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                 !< GLM cleaning speed (EGLM only).
   real(R8P),    intent(in)    :: damping                            !< GLM damping rate (EGLM only).
   integer(I4P), intent(in)    :: hs                                 !< Half stencil of the high-order sources.
   real(R8P),    intent(in)    :: dt                                 !< Forward-Euler step of the stage.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                        !< Blocks space steps [3, nb].
   logical,      intent(in)    :: is_null(3)                         !< Null directions.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Conservative variables (ghosts updated).
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   real(R8P),    intent(in)    :: flx(1:,0:,1:,1:,1:)                !< X-face high-order fluxes.
   real(R8P),    intent(in)    :: fly(1:,1:,0:,1:,1:)                !< Y-face high-order fluxes.
   real(R8P),    intent(in)    :: flz(1:,1:,1:,0:,1:)                !< Z-face high-order fluxes.
   real(R8P),    intent(inout) :: lam(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Cell factors (component 1).
   integer(I4P), intent(out)   :: bad                                !< Cells with an inadmissible backbone.
   integer(I4P), intent(out)   :: nonfinite                          !< Cells with a non-finite high-order face flux.
   real(R8P)                   :: qc(NV_K), qn(NV_K)                 !< Cell and neighbour states.
   real(R8P)                   :: flo(NV_K)                          !< Backbone face flux.
   real(R8P)                   :: dlt(NV_K,6)                        !< Antidiffusive difference of each active face.
   real(R8P)                   :: qb(NV_K)                           !< Backbone update.
   real(R8P)                   :: dsrc(NV_K)                         !< Source difference dt (s_hi - s_lo).
   real(R8P)                   :: s_lo(NV_K), s_hi(NV_K)             !< Backbone and high-order sources.
   real(R8P)                   :: lc                                 !< Cell factor.
   real(R8P)                   :: dtdx                               !< dt / dx_d.
   integer(I4P)                :: ids(3,3)                           !< Unit steps of the directions.
   integer(I4P)                :: nf, st                             !< Active faces, factor status.
   integer(I4P)                :: b, i, j, k, d                      !< Counters.

   ids = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   bad = 0_I4P
   nonfinite = 0_I4P
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q, q_aux, flx, fly, flz, lam) &
   !$omp& reduction(+:bad, nonfinite)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               qc = q(1:NV_K,i,j,k,b)
               qb = qc
               nf = 0
               do d=1, 3
                  if (is_null(d)) cycle
                  dtdx = dt / dxyz(d,b)
                  ! right face i+1/2: the backbone flux and the high-order one
                  qn = q(1:NV_K,i+ids(1,d),j+ids(2,d),k+ids(3,d),b)
                  call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qc, qR=qn, f=flo)
                  nf = nf + 1
                  select case(d)
                  case(1)
                     dlt(:,nf) = -dtdx * (flx(1:NV_K,i,j,k,b) - flo)
                  case(2)
                     dlt(:,nf) = -dtdx * (fly(1:NV_K,i,j,k,b) - flo)
                  case(3)
                     dlt(:,nf) = -dtdx * (flz(1:NV_K,i,j,k,b) - flo)
                  endselect
                  qb = qb - dtdx * flo
                  ! left face i-1/2
                  qn = q(1:NV_K,i-ids(1,d),j-ids(2,d),k-ids(3,d),b)
                  call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qn, qR=qc, f=flo)
                  nf = nf + 1
                  select case(d)
                  case(1)
                     dlt(:,nf) = dtdx * (flx(1:NV_K,i-1,j,k,b) - flo)
                  case(2)
                     dlt(:,nf) = dtdx * (fly(1:NV_K,i,j-1,k,b) - flo)
                  case(3)
                     dlt(:,nf) = dtdx * (flz(1:NV_K,i,j,k-1,b) - flo)
                  endselect
                  qb = qb + dtdx * flo
               enddo
               call cell_sources(ngc=ngc, hs=hs, damping=damping, dxyz=dxyz(:,b), is_null=is_null, q=q, q_aux=q_aux, &
                                 i=i, j=j, k=k, b=b, s_lo=s_lo, s_hi=s_hi)
               qb = qb + dt * s_lo
               dsrc = dt * (s_hi - s_lo)
               call positivity_cell_factor(nf=nf, qb=qb, dlt=dlt, dsrc=dsrc, lc=lc, status=st)
               if (st == 1_I4P) bad = bad + 1_I4P
               if (st == 2_I4P) nonfinite = nonfinite + 1_I4P
               lam(1,i,j,k,b) = lc
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_positivity_factors

   subroutine compute_seam_positivity_factors(ngc, gamma, ch, damping, hs, dt, dxyz, is_null, q, q_aux, flx, fly, flz, &
                                              ncells, cell, omask, olo, ohi, lam, dbad, dnonfinite)
   !< Recompute the factor of the cells beside a 2:1 seam with the seam faces' backbone and high-order fluxes overridden
   !< (issue #50, D1: per-stage seam flux synchronisation). Cell `n` is `cell(1:4,n) = (i, j, k, b)`; its face `fec`
   !< (1..6: -x, +x, -y, +y, -z, +z) takes `olo(:,fec,n)` and `ohi(:,fec,n)` where `omask(fec,n)`, the backbone flux of
   !< `q` and the high-order flux of `flx/fly/flz` elsewhere, exactly as `compute_positivity_factors`. `dbad` and
   !< `dnonfinite` return the change of the two counts against the factors without the overrides, which
   !< `compute_positivity_factors` has already counted.
   integer(I4P), intent(in)    :: ngc                                !< Ghost cells number.
   real(R8P),    intent(in)    :: gamma                              !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                                 !< GLM cleaning speed (EGLM only).
   real(R8P),    intent(in)    :: damping                            !< GLM damping rate (EGLM only).
   integer(I4P), intent(in)    :: hs                                 !< Half stencil of the high-order sources.
   real(R8P),    intent(in)    :: dt                                 !< Forward-Euler step of the stage.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                        !< Blocks space steps [3, nb].
   logical,      intent(in)    :: is_null(3)                         !< Null directions.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Conservative variables (ghosts updated).
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   real(R8P),    intent(in)    :: flx(1:,0:,1:,1:,1:)                !< X-face high-order fluxes.
   real(R8P),    intent(in)    :: fly(1:,1:,0:,1:,1:)                !< Y-face high-order fluxes.
   real(R8P),    intent(in)    :: flz(1:,1:,1:,0:,1:)                !< Z-face high-order fluxes.
   integer(I4P), intent(in)    :: ncells                             !< Seam cells.
   integer(I4P), intent(in)    :: cell(1:,1:)                        !< Seam cells (i, j, k, b) [4, ncells].
   logical,      intent(in)    :: omask(1:,1:)                       !< Overridden faces [6, ncells].
   real(R8P),    intent(in)    :: olo(1:,1:,1:)                      !< Override backbone fluxes [nv, 6, ncells].
   real(R8P),    intent(in)    :: ohi(1:,1:,1:)                      !< Override high-order fluxes [nv, 6, ncells].
   real(R8P),    intent(inout) :: lam(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Cell factors (component 1).
   integer(I4P), intent(out)   :: dbad                               !< Change of the inadmissible backbones count.
   integer(I4P), intent(out)   :: dnonfinite                         !< Change of the non-finite cells count.
   real(R8P)                   :: qc(NV_K), qn(NV_K)                 !< Cell and neighbour states.
   real(R8P)                   :: flo(NV_K), fhi(NV_K)               !< Backbone and high-order face fluxes.
   real(R8P)                   :: dlt(NV_K,6), dlt0(NV_K,6)          !< Antidiffusive differences, with/without overrides.
   real(R8P)                   :: qb(NV_K), qb0(NV_K)                !< Backbone updates, with/without overrides.
   real(R8P)                   :: dsrc(NV_K)                         !< Source difference dt (s_hi - s_lo).
   real(R8P)                   :: s_lo(NV_K), s_hi(NV_K)             !< Backbone and high-order sources.
   real(R8P)                   :: dtdx, sg, lc, lc0                  !< dt / dx_d, face sign; factors with/without overrides.
   integer(I4P)                :: ids(3,3)                           !< Unit steps of the directions.
   integer(I4P)                :: nf, st, st0                        !< Active faces; statuses with/without overrides.
   integer(I4P)                :: n, i, j, k, b, d, s, fec           !< Counters.

   ids = reshape([1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   dbad = 0_I4P
   dnonfinite = 0_I4P
   do n=1, ncells
      i = cell(1,n) ; j = cell(2,n) ; k = cell(3,n) ; b = cell(4,n)
      qc = q(1:NV_K,i,j,k,b)
      qb = qc
      qb0 = qc
      nf = 0
      do d=1, 3
         if (is_null(d)) cycle
         dtdx = dt / dxyz(d,b)
         do s=1, 0, -1 ! the right face (s = 1, fec = 2d) first, then the left one, as compute_positivity_factors
            fec = 2 * d - 1 + s
            if (s == 1) then
               qn = q(1:NV_K,i+ids(1,d),j+ids(2,d),k+ids(3,d),b)
               call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qc, qR=qn, f=flo)
               select case(d)
               case(1)
                  fhi = flx(1:NV_K,i,j,k,b)
               case(2)
                  fhi = fly(1:NV_K,i,j,k,b)
               case default
                  fhi = flz(1:NV_K,i,j,k,b)
               endselect
               sg = -dtdx
            else
               qn = q(1:NV_K,i-ids(1,d),j-ids(2,d),k-ids(3,d),b)
               call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=qn, qR=qc, f=flo)
               select case(d)
               case(1)
                  fhi = flx(1:NV_K,i-1,j,k,b)
               case(2)
                  fhi = fly(1:NV_K,i,j-1,k,b)
               case default
                  fhi = flz(1:NV_K,i,j,k-1,b)
               endselect
               sg = dtdx
            endif
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
      call cell_sources(ngc=ngc, hs=hs, damping=damping, dxyz=dxyz(:,b), is_null=is_null, q=q, q_aux=q_aux, &
                        i=i, j=j, k=k, b=b, s_lo=s_lo, s_hi=s_hi)
      qb = qb + dt * s_lo
      qb0 = qb0 + dt * s_lo
      dsrc = dt * (s_hi - s_lo)
      call positivity_cell_factor(nf=nf, qb=qb0, dlt=dlt0, dsrc=dsrc, lc=lc0, status=st0)
      call positivity_cell_factor(nf=nf, qb=qb, dlt=dlt, dsrc=dsrc, lc=lc, status=st)
      dbad = dbad + merge(1_I4P, 0_I4P, st == 1_I4P) - merge(1_I4P, 0_I4P, st0 == 1_I4P)
      dnonfinite = dnonfinite + merge(1_I4P, 0_I4P, st == 2_I4P) - merge(1_I4P, 0_I4P, st0 == 2_I4P)
      lam(1,i,j,k,b) = lc
   enddo
   endsubroutine compute_seam_positivity_factors

   subroutine compute_backbone_fluxes(n, d, gamma, ch, qL, qR, flo)
   !< Compute the backbone flux of `n` face state pairs of direction `d` (issue #50: the fine faces of a 2:1 seam take
   !< the coarse donor cell as their outer state).
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
   endsubroutine compute_backbone_fluxes

   pure subroutine positivity_cell_factor(nf, qb, dlt, dsrc, lc, status)
   !< The factor of one cell from its backbone update `qb`, the antidiffusive differences `dlt(:,1:nf)` of its active
   !< faces and the source difference `dsrc`: `status` 1 (and `lc = 0`) if the backbone is not admissible, 2 (and
   !< `lc = 0`) if a face difference is not finite (a NaN passes every comparison of the corners), else 0 with `lc` the
   !< largest factor of the corner box.
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
   ! the relative floor: a stage may not take the cell below KAPPA times its backbone value (no absolute floor: issue
   ! #49, the scheme has no unit), else a drained cell (its energy kept) makes the next stage's backbone inadmissible
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

   subroutine blend_positivity_fluxes(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, ch, q, lam, fl, limited)
   !< Blend the face fluxes of direction `d` with the backbone: `F <- F^LLF + theta (F - F^LLF)`, `theta` the smaller
   !< factor of the two cells (ghost factors from the exchange; a ghost never written keeps 1, so a physical-boundary
   !< face takes the interior cell's factor). Faces with `theta = 1` are left untouched (bitwise), and counted otherwise;
   !< a face with `theta = 0` takes the backbone flux itself (`0 * NaN` would keep a non-finite high-order flux).
   integer(I4P), intent(in)    :: d                                !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                       !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                  !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                    !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma                            !< Specific heats ratio.
   real(R8P),    intent(in)    :: ch                               !< GLM cleaning speed (EGLM only).
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Conservative variables.
   real(R8P),    intent(in)    :: lam(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Cell factors (component 1).
   real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)      !< Face fluxes of direction `d`.
   integer(I4P), intent(out)   :: limited                          !< Faces with theta < 1.
   real(R8P)                   :: flo(NV_K)                        !< Backbone face flux.
   real(R8P)                   :: theta                            !< Face factor.
   integer(I4P)                :: b, i, j, k                       !< Counters.

   limited = 0_I4P
   !$omp parallel do collapse(4) default(firstprivate) shared(q, lam, fl) reduction(+:limited)
   do b=1, blocks_number
      do k=1-dk, nk
         do j=1-dj, nj
            do i=1-di, ni
               theta = min(lam(1,i,j,k,b), lam(1,i+di,j+dj,k+dk,b), 1._R8P)
               if (theta < 1._R8P) then
                  call backbone_face_flux(gamma=gamma, ch=ch, d=d, qL=q(1:NV_K,i,j,k,b), qR=q(1:NV_K,i+di,j+dj,k+dk,b), &
                                          f=flo)
                  if (theta > 0._R8P) then
                     fl(1:NV_K,i,j,k,b) = flo + theta * (fl(1:NV_K,i,j,k,b) - flo)
                  else
                     fl(1:NV_K,i,j,k,b) = flo
                  endif
                  limited = limited + 1_I4P
               endif
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine blend_positivity_fluxes
