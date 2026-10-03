   ! Model-agnostic FNL device loop of the ghost positivity blend (issue #50, D2): the device twin of
   ! `adam_flume_cpu_ghost_kernels_agnostic.INC` (see it for the method), on the device layout `[nb, i, j, k, nv]`.
   !
   ! Included in the `contains` part of every model's device kernel module, which provides by use association `NV_K`,
   ! `IQ_R`, `POSITIVITY_LIMITER_EPS`, `POSITIVITY_LIMITER_KAPPA` and the device adapter `internal_energy(q)`.

   subroutine blend_inadmissible_ghosts_dev(ni, nj, nk, ngc, blocks_number, is_null, q_gpu, blended)
   !< Blend every inadmissible face ghost of `q_gpu` toward the adjacent interior cell and count the ghosts blended
   !< (device twin of the CPU `blend_inadmissible_ghosts`).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                    !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                      !< Actual blocks number.
   logical,      intent(in)    :: is_null(3)                         !< Null directions.
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables (ghosts updated).
   integer(I4P), intent(out)   :: blended                            !< Ghosts blended.
   integer(I8P), parameter     :: EXPONENT_MASK=2047_I8P             !< The 11 exponent bits of a binary64.
   logical                     :: isnf, ok                           !< Private non-finite ghost, admissible ghost flags.
   real(R8P)                   :: qg(NV_K), qi(NV_K), dc(NV_K)       !< Private ghost, interior anchor, difference.
   real(R8P)                   :: qt(NV_K)                           !< Private segment state.
   real(R8P)                   :: rb, eb, er, ee, r1, e1, t          !< Private anchor values, floors, segment values.
   integer(I4P)                :: lo(3), hi(3), ia(3)                !< Ghost slab bounds, interior index along `d`.
   integer(I4P)                :: il, iu, jl, ju, kl, ku             !< Ghost slab bounds (scalar copies for the launch).
   integer(I4P)                :: ai, aj, ak                         !< Interior index along `d` (0 elsewhere).
   integer(I4P)                :: ii, jj, kk                         !< Private interior anchor indexes.
   integer(I4P)                :: b, i, j, k, d, s, g, v, n          !< Counters.

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
            il = lo(1) ; iu = hi(1) ; jl = lo(2) ; ju = hi(2) ; kl = lo(3) ; ku = hi(3)
            ai = ia(1) ; aj = ia(2) ; ak = ia(3)
            n = 0_I4P
            !$acc parallel loop independent gang vector collapse(4) reduction(+:n) DEVICEVAR(q_gpu) &
            !$acc& firstprivate(d,ai,aj,ak,blocks_number) private(isnf,ok,qg,qi,dc,qt,rb,eb,er,ee,r1,e1,t,ii,jj,kk)
            !$omp OMPLOOP collapse(4) reduction(+:n) DEVICEPTR(q_gpu) &
            !$omp& firstprivate(d,ai,aj,ak,blocks_number) private(isnf,ok,qg,qi,dc,qt,rb,eb,er,ee,r1,e1,t,ii,jj,kk)
            do k=kl, ku
            do j=jl, ju
            do i=il, iu
            do b=1, blocks_number
               ii = merge(ai, i, d == 1)
               jj = merge(aj, j, d == 2)
               kk = merge(ak, k, d == 3)
               isnf = .false.
               !$acc loop seq
               do v=1, NV_K
                  qg(v) = q_gpu(b,i,j,k,v)
                  qi(v) = q_gpu(b,ii,jj,kk,v)
                  if (iand(ishft(transfer(qg(v), 0_I8P), -52), EXPONENT_MASK) == EXPONENT_MASK) isnf = .true.
               enddo
               ok = .false.
               if (.not.isnf .and. qg(IQ_R) > 0._R8P) ok = internal_energy(q=qg) > 0._R8P
               rb = qi(IQ_R)
               eb = -1._R8P
               if (rb > 0._R8P) eb = internal_energy(q=qi)
               if ((.not.ok) .and. rb > 0._R8P .and. eb > 0._R8P) then
                  n = n + 1_I4P
                  t = 0._R8P
                  if (.not.isnf) then
                     er = max(min(POSITIVITY_LIMITER_EPS, 0.5_R8P * rb), POSITIVITY_LIMITER_KAPPA * rb)
                     ee = max(min(POSITIVITY_LIMITER_EPS, 0.5_R8P * eb), POSITIVITY_LIMITER_KAPPA * eb)
                     !$acc loop seq
                     do v=1, NV_K
                        dc(v) = qg(v) - qi(v)
                     enddo
                     t = 1._R8P
                     r1 = rb + dc(IQ_R)
                     if (r1 < er) t = (rb - er) / (rb - r1)
                     !$acc loop seq
                     do v=1, NV_K
                        qt(v) = qi(v) + t * dc(v)
                     enddo
                     e1 = internal_energy(q=qt)
                     if (e1 < ee) t = t * (eb - ee) / (eb - e1)
                  endif
                  !$acc loop seq
                  do v=1, NV_K
                     if (t > 0._R8P) then
                        q_gpu(b,i,j,k,v) = qi(v) + t * dc(v)
                     else
                        q_gpu(b,i,j,k,v) = qi(v)
                     endif
                  enddo
               endif
            enddo
            enddo
            enddo
            enddo
            blended = blended + n
         enddo
      enddo
   enddo
   endsubroutine blend_inadmissible_ghosts_dev
