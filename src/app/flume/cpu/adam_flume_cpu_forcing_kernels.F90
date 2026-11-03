!< ADAM, FLUME CPU kernels of the mass-flow forcing (issue #119, T3.1).

module adam_flume_cpu_forcing_kernels
!< ADAM, FLUME CPU kernels of the mass-flow forcing (issue #119, T3.1).

use :: adam_flume_parameters, only : IA_P, IA_R, IA_U, IA_V, IA_W, IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW
use :: adam_flume_wall_model_library, only : sample_cells, wall_model_solve, wall_model_solve_inverse
use :: penf,                  only : I4P, R8P

implicit none
private
public :: add_uniform_momentum
public :: overwrite_wall_fluxes
public :: close_wall_fluxes
public :: sum_momentum

contains
   subroutine sum_momentum(ni, nj, nk, ngc, blocks_number, axis, dxyz, q, sums)
   !< Sum the momentum along `axis` and the volume over the interior cells; the volume includes the null directions, as the conservation integrals.
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                  !< Actual blocks number.
   integer(I4P), intent(in)  :: axis                           !< Axis of the flow.
   real(R8P),    intent(in)  :: dxyz(1:,1:)                    !< Blocks space steps [3, nb].
   real(R8P),    intent(in)  :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   real(R8P),    intent(out) :: sums(2)                        !< Momentum times volume, volume.
   real(R8P)                 :: s1, s2                         !< Accumulators.
   integer(I4P)              :: b, i, j, k                     !< Counters.

   s1 = 0._R8P ; s2 = 0._R8P
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q) reduction(+:s1,s2)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               s1 = s1 + q(IQ_RU+axis-1,i,j,k,b) * dxyz(1,b) * dxyz(2,b) * dxyz(3,b)
               s2 = s2 + dxyz(1,b) * dxyz(2,b) * dxyz(3,b)
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   sums = [s1, s2]
   endsubroutine sum_momentum

   subroutine add_uniform_momentum(ni, nj, nk, ngc, blocks_number, axis, dm, q)
   !< Add `dm` to the momentum along `axis`, and to the total energy the change of kinetic energy.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                  !< Actual blocks number.
   integer(I4P), intent(in)    :: axis                           !< Axis of the flow.
   real(R8P),    intent(in)    :: dm                             !< Momentum added.
   real(R8P),    intent(inout) :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   integer(I4P)                :: b, i, j, k                     !< Counters.

   !$omp parallel do collapse(4) default(firstprivate) shared(q)
   do b=1, blocks_number
      do k=1, nk
         do j=1, nj
            do i=1, ni
               q(IQ_RE,i,j,k,b) = q(IQ_RE,i,j,k,b) + dm * (q(IQ_RU+axis-1,i,j,k,b) + 0.5_R8P * dm) / q(IQ_R,i,j,k,b)
               q(IQ_RU+axis-1,i,j,k,b) = q(IQ_RU+axis-1,i,j,k,b) + dm
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine add_uniform_momentum
   subroutine overwrite_wall_fluxes(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, R, pr, mu0, tref, omega_mu, faces, tw, hwm, sl, tf, dt, commit, fresh, &
                                    avg, inv, ue, te, dxyz, q_aux, fl)
   !< Overwrite the flux through the modelled wall faces of direction `d` (issue #119, T3.3): no mass, the parallel momentum by the modelled shear, the energy by the heat
   !< flux; the sample is the 3rd cell off the wall. The normal momentum keeps the mirrored inviscid flux, the pressure.
   !< With a filter time `tf` the sample is `(1 - e) a + e s`, `e = min(dt / tf, 1)` (1 when `fresh`), from the averages `a` of the face, which take it when `commit`
   !< (once per step, Kawai & Larsson 2012): `tf` -> 0 is the unfiltered model, bit for bit.
   integer(I4P), intent(in)    :: d                               !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                      !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                 !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                   !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma, R, pr                    !< Specific heats ratio, gas constant, Prandtl number.
   real(R8P),    intent(in)    :: mu0, tref, omega_mu             !< Viscosity law.
   integer(I4P), intent(in)    :: faces(1:,1:)                   !< Modelled block faces [nb, 6].
   real(R8P),    intent(in)    :: tw(6)                           !< Wall temperature per face.
   real(R8P),    intent(in)    :: hwm(6)                          !< Matching height per face, 0 for the 3rd cell.
   logical,      intent(in)    :: sl(6)                           !< Semi-local damping per face.
   real(R8P),    intent(in)    :: tf(6)                           !< Filter time per face, 0 for none.
   real(R8P),    intent(in)    :: dt                              !< Time step.
   logical,      intent(in)    :: commit                          !< The averages take this sample.
   logical,      intent(in)    :: fresh                           !< The averages hold no history.
   real(R8P),    intent(inout) :: avg(1:,1:,1:,1:,1:)             !< Averages of the faces of direction `d` [5, ni or 2, nj or 2, nk or 2, nb].
   logical,      intent(in)    :: inv(6)                          !< The inverse-transformation model per face (T3.5), else the equilibrium one.
   real(R8P),    intent(in)    :: ue(6), te(6)                    !< Inverse model: speed and temperature at the edge of the layer per face.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                    !< Blocks space steps.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)      !< Face fluxes of direction `d`.
   real(R8P)                   :: u(3), us, ts, ps, rs, tau, q, cp, h, s, w, e !< Work.
   integer(I4P)                :: n, side, f, b, i, j, k, l, m, o1, o2, ia, ja, ka !< Counters.

   n = ni * di + nj * dj + nk * dk
   cp = gamma * (R / (gamma - 1._R8P)) ! as the physics object derives it (cv = R / (gamma - 1), cp = gamma cv)
   do side=1, 2
      f = 2_I4P * d - 2_I4P + side
      l = merge(0_I4P, n, side == 1)
      s = merge(1._R8P, -1._R8P, side == 1) ! the flux along +d of a shear pointing to +u at y_min is -tau
      e = 1._R8P
      if (tf(f) > 0._R8P .and. .not.fresh) e = min(dt / tf(f), 1._R8P)
      !$omp parallel do collapse(4) default(firstprivate) shared(faces, dxyz, q_aux, fl, tw, avg)
      do b=1, blocks_number
         do k=merge(l, 1_I4P, dk == 1_I4P), merge(l, nk, dk == 1_I4P)
            do j=merge(l, 1_I4P, dj == 1_I4P), merge(l, nj, dj == 1_I4P)
               do i=merge(l, 1_I4P, di == 1_I4P), merge(l, ni, di == 1_I4P)
                  if (faces(b,f) == 0_I4P) cycle
                  ! the sample: the 3rd cell off the wall, or the matching height between cells m and m + 1 counted from the wall (centre of cell m at
                  ! (m - 1/2) dx), within the block; cell m is i + m di on the low side, i - (m - 1) di on the high side, the face index being l
                  call sample_cells(hwm=hwm(f), dx=dxyz(d,b), n=n, m=m, w=w, h=h)
                  o1 = merge(m, 1_I4P - m, side == 1) ; o2 = merge(m + 1_I4P, -m, side == 1)
                  u = q_aux(IA_U:IA_W,i+di*o1,j+dj*o1,k+dk*o1,b)
                  ps = q_aux(IA_P,i+di*o1,j+dj*o1,k+dk*o1,b)
                  rs = q_aux(IA_R,i+di*o1,j+dj*o1,k+dk*o1,b)
                  if (w > 0._R8P) then
                     u = u + w * (q_aux(IA_U:IA_W,i+di*o2,j+dj*o2,k+dk*o2,b) - u)
                     ps = ps + w * (q_aux(IA_P,i+di*o2,j+dj*o2,k+dk*o2,b) - ps)
                     rs = rs + w * (q_aux(IA_R,i+di*o2,j+dj*o2,k+dk*o2,b) - rs)
                  endif
                  u(d) = 0._R8P
                  if (tf(f) > 0._R8P) then ! the input filter; e = 1 takes the sample exactly
                     ia = merge(side, i, di == 1_I4P) ; ja = merge(side, j, dj == 1_I4P) ; ka = merge(side, k, dk == 1_I4P)
                     u = (1._R8P - e) * avg(1:3,ia,ja,ka,b) + e * u
                     ps = (1._R8P - e) * avg(4,ia,ja,ka,b) + e * ps
                     rs = (1._R8P - e) * avg(5,ia,ja,ka,b) + e * rs
                     if (commit) avg(1:5,ia,ja,ka,b) = [u(1), u(2), u(3), ps, rs]
                  endif
                  us = sqrt(dot_product(u, u))
                  ts = ps / (R * rs)
                  if (inv(f)) then
                     call wall_model_solve_inverse(h=h, us=us, ts=ts, ps=ps, tw=tw(f), ue=ue(f), te=te(f), R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, omega_mu=omega_mu, &
                                                   tau=tau, q=q)
                  else
                     call wall_model_solve(h=h, us=us, ts=ts, ps=ps, tw=tw(f), R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, omega_mu=omega_mu, tau=tau, q=q, semilocal=sl(f))
                  endif
                  if (us > 0._R8P) u = u / us
                  fl(IQ_R,i,j,k,b) = 0._R8P
                  fl(IQ_RU:IQ_RW,i,j,k,b) = merge(fl(IQ_RU:IQ_RW,i,j,k,b), -s * tau * u, [1,2,3] == d)
                  fl(IQ_RE,i,j,k,b) = -s * q
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   enddo
   endsubroutine overwrite_wall_fluxes
   subroutine close_wall_fluxes(d, di, dj, dk, ni, nj, nk, nv, blocks_number, faces, fl)
   !< Zero every flux of direction `d` through the flagged block faces of that direction but the normal momentum (issue #119, the #65 mass leak).
   integer(I4P), intent(in)    :: d                           !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                  !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, nv              !< Grid dimensions, variables.
   integer(I4P), intent(in)    :: blocks_number               !< Actual blocks number.
   integer(I4P), intent(in)    :: faces(1:,1:)                !< Flagged block faces [nb, 6].
   real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:) !< Face fluxes of direction `d`.
   real(R8P)                   :: keep                        !< The normal momentum flux.
   integer(I4P)                :: n, side, f, b, i, j, k, l   !< Counters.

   n = ni * di + nj * dj + nk * dk
   do side=1, 2
      f = 2_I4P * d - 2_I4P + side
      l = merge(0_I4P, n, side == 1)
      !$omp parallel do collapse(4) default(firstprivate) shared(faces, fl)
      do b=1, blocks_number
         do k=merge(l, 1_I4P, dk == 1_I4P), merge(l, nk, dk == 1_I4P)
            do j=merge(l, 1_I4P, dj == 1_I4P), merge(l, nj, dj == 1_I4P)
               do i=merge(l, 1_I4P, di == 1_I4P), merge(l, ni, di == 1_I4P)
                  if (faces(b,f) == 0_I4P) cycle
                  keep = fl(IQ_RU+d-1,i,j,k,b)
                  fl(1:nv,i,j,k,b) = 0._R8P
                  fl(IQ_RU+d-1,i,j,k,b) = keep
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   enddo
   endsubroutine close_wall_fluxes
endmodule adam_flume_cpu_forcing_kernels
