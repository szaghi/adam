   ! Model-agnostic CPU loop of the ghost positivity blend (issue #50, D2).
   !
   ! Included in the `contains` part of every model's kernel module, which provides by use association `NV_K`, `IQ_R`,
   ! `POSITIVITY_LIMITER_KAPPA` and the private adapter `internal_energy(q)` (the internal
   ! energy per unit volume, a concave function of `q`).
   !
   ! Why: a coarse-to-fine seam ghost is a tricubic interpolant of each conservative variable (`[amr] seam_ghost_fill`,
   ! issue #21), and the interpolant of an admissible field need not be admissible. At low beta the internal energy is a
   ! small residual of `E - |B|^2/2`, and on the Balsara-Spicer blast across a 2:1 seam the fine ghosts beside the shock
   ! reach `rho e < 0` (issue #50): the high-order face flux turns non-finite and the limiter's backbone, fed an
   ! inadmissible state, is no longer admissible. The blend `g <- q_in + t (g - q_in)` toward the adjacent interior cell
   ! `q_in` of the same block takes the largest `t` in [0, 1] that keeps the density and the internal energy above the
   ! limiter's floors of `q_in` (the closed form of the limiter's corners: the density along the segment, then the
   ! internal energy through its chord); a non-finite ghost takes `q_in`. Admissible ghosts are not touched (bitwise), nor
   ! is a ghost whose interior cell is itself inadmissible (no admissible anchor). Face ghosts only, in the active
   ! directions: the stencils of the split and Riemann schemes and the EGLM sources read no other ghost.

   subroutine blend_inadmissible_ghosts(ni, nj, nk, ngc, blocks_number, is_null, q, blended)
   !< Blend every inadmissible face ghost of `q` toward the adjacent interior cell and count the ghosts blended.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   logical,      intent(in)    :: is_null(3)                         !< Null directions.
   real(R8P),    intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Conservative variables (ghosts updated).
   integer(I4P), intent(out)   :: blended                            !< Ghosts blended.
   integer(I8P), parameter     :: EXPONENT_MASK=2047_I8P             !< The 11 exponent bits of a binary64.
   logical                     :: isnf                               !< The ghost holds a non-finite value.
   real(R8P)                   :: qg(NV_K), qi(NV_K), dc(NV_K)       !< Ghost, interior anchor, their difference.
   real(R8P)                   :: rb, eb, er, ee, r1, e1, t          !< Anchor density, energy, floors; segment values.
   integer(I4P)                :: lo(3), hi(3), ia(3)                !< Ghost slab bounds, interior index along `d`.
   integer(I4P)                :: b, i, j, k, d, s, g, c             !< Counters.

   blended = 0_I4P
   do d=1, 3
      if (is_null(d)) cycle
      do s=0, 1
         do g=1, ngc
            lo = 1_I4P
            hi = [ni, nj, nk]
            ia = 0_I4P
            ia(d) = merge(1_I4P, hi(d), s == 0)
            lo(d) = merge(1_I4P - g, hi(d) + g, s == 0)
            hi(d) = lo(d)
            !$omp parallel do collapse(4) default(firstprivate) shared(q) reduction(+:blended)
            do b=1, blocks_number
               do k=lo(3), hi(3)
                  do j=lo(2), hi(2)
                     do i=lo(1), hi(1)
                        qg = q(1:NV_K,i,j,k,b)
                        isnf = .false.
                        do c=1, NV_K
                           if (iand(ishft(transfer(qg(c), 0_I8P), -52), EXPONENT_MASK) == EXPONENT_MASK) isnf = .true.
                        enddo
                        if (.not.isnf .and. qg(IQ_R) > 0._R8P) then
                           if (internal_energy(q=qg) > 0._R8P) cycle
                        endif
                        qi = q(1:NV_K,merge(ia(1), i, d == 1),merge(ia(2), j, d == 2),merge(ia(3), k, d == 3),b)
                        rb = qi(IQ_R)
                        if (.not.(rb > 0._R8P)) cycle
                        eb = internal_energy(q=qi)
                        if (.not.(eb > 0._R8P)) cycle
                        blended = blended + 1_I4P
                        if (isnf) then
                           q(1:NV_K,i,j,k,b) = qi
                           cycle
                        endif
                        er = POSITIVITY_LIMITER_KAPPA * rb
                        ee = POSITIVITY_LIMITER_KAPPA * eb
                        dc = qg - qi
                        t = 1._R8P
                        r1 = rb + dc(IQ_R)
                        if (r1 < er) t = (rb - er) / (rb - r1)
                        e1 = internal_energy(q=qi + t * dc)
                        if (e1 < ee) t = t * (eb - ee) / (eb - e1)
                        if (t > 0._R8P) then
                           q(1:NV_K,i,j,k,b) = qi + t * dc
                        else
                           q(1:NV_K,i,j,k,b) = qi
                        endif
                     enddo
                  enddo
               enddo
            enddo
            !$omp end parallel do
         enddo
      enddo
   enddo
   endsubroutine blend_inadmissible_ghosts
