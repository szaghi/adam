!< ADAM, FLUME CPU kernels of the mass-flow forcing (issue #119, T3.1).

module adam_flume_cpu_forcing_kernels
!< ADAM, FLUME CPU kernels of the mass-flow forcing (issue #119, T3.1).

use :: adam_flume_dissipation_library, only : temperature_law
use :: adam_flume_parameters, only : IA_P, IA_R, IA_U, IA_V, IA_W, IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW
use :: adam_flume_wall_model_library, only : sample_cells, wall_eddy_viscosity, wall_model_solve, wall_model_solve_inverse, wall_model_solve_pressure, &
                                         wall_model_solve_convective, WM_CONV_CELLS
use :: penf,                  only : I4P, R8P

implicit none
private
public :: add_uniform_momentum
public :: add_wall_mu_t
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
   subroutine overwrite_wall_fluxes(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, R, pr, mu0, tref, suth, omega_mu, faces, tw, hwm, sl, tf, dt, commit, fresh, &
                                    avg, inv, ue, te, pg, yp, cv, rsf, tp, is_null, dxyz, q_aux, fl)
   !< Overwrite the flux through the modelled wall faces of direction `d` (issue #119, T3.3): no mass, the parallel momentum by the modelled shear, the energy by the heat
   !< flux; the sample is the 3rd cell off the wall. The normal momentum keeps the mirrored inviscid flux, the pressure.
   !< With a filter time `tf` the sample is `(1 - e) a + e s`, `e = min(dt / tf, 1)` (1 when `fresh`), from the averages `a` of the face, which take it when `commit`
   !< (once per step, Kawai & Larsson 2012): `tf` -> 0 is the unfiltered model, bit for bit.
   integer(I4P), intent(in)    :: d                               !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                      !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                 !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                   !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma, R, pr                    !< Specific heats ratio, gas constant, Prandtl number.
   real(R8P),    intent(in)    :: mu0, tref, suth, omega_mu             !< Viscosity law.
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
   logical,      intent(in)    :: pg(6)                           !< The pressure gradient in the equilibrium model per face (issue #121, S4d).
   real(R8P),    intent(in)    :: yp(6)                           !< y_pg of the Bernoulli-consistent pressure gradient per face, 0 for none (#121 P1).
   logical,      intent(in)    :: cv(6)                           !< The model fed with the resolved convection per face (#121 P2).
   logical,      intent(in)    :: rsf(6)                          !< Its resolved shear stress removed per face (#121 P3).
   real(R8P),    intent(in)    :: tp(6)                           !< Filter time of the pressure gradient per face, 0 for none.
   logical,      intent(in)    :: is_null(3)                      !< Null directions.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                    !< Blocks space steps.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)      !< Face fluxes of direction `d`.
   real(R8P)                   :: u(3), us, ts, ps, rs, tau, q, cp, h, s, w, e, px, ep, cc(WM_CONV_CELLS), rc(WM_CONV_CELLS) !< Work.
   real(R8P)                   :: rv(4), ue_, un_, rr, sdr     !< Resolved-stress work: the four averages, u_e, u_eta, rho, sign away from the wall.
   integer(I4P)                :: n, side, f, b, i, j, k, l, m, o1, o2, ia, ja, ka, nc, c, oc, kk, ab, cl(3) !< Counters.

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
                  if (pg(f)) then
                     px = sample_pressure_gradient(i=i, j=j, k=k, b=b, o1=o1, o2=o2, w=w, u=u, us=us)
                     if (tp(f) > 0._R8P) then ! its own exponential filter (the turbulent gradient is as large as tau / h): avg(6)
                        ep = 1._R8P
                        if (.not.fresh) ep = min(dt / tp(f), 1._R8P)
                        ia = merge(side, i, di == 1_I4P) ; ja = merge(side, j, dj == 1_I4P) ; ka = merge(side, k, dk == 1_I4P)
                        px = (1._R8P - ep) * avg(6,ia,ja,ka,b) + ep * px
                        if (commit) avg(6,ia,ja,ka,b) = px
                     endif
                     if (cv(f)) then ! the model fed with the convection of the cells up to the sample (#121 P2), filtered as the gradient: avg(7:)
                        nc = min(m + merge(1_I4P, 0_I4P, w > 0._R8P), WM_CONV_CELLS, n - 1_I4P)
                        cc = 0._R8P
                        do c=1, nc
                           cc(c) = cell_convection(i=i, j=j, k=k, b=b, c=c, u=u, us=us)
                        enddo
                        ia = merge(side, i, di == 1_I4P) ; ja = merge(side, j, dj == 1_I4P) ; ka = merge(side, k, dk == 1_I4P)
                        do c=1, nc
                           cc(c) = (1._R8P - ep) * avg(6+c,ia,ja,ka,b) + ep * cc(c)
                           if (commit) avg(6+c,ia,ja,ka,b) = cc(c)
                        enddo
                        rc = 0._R8P
                        if (rsf(f)) then ! the resolved shear stress of the same cells (#121 P3): rho, rho u_e, rho u_eta, rho u_e u_eta filtered
                           sdr = merge(1._R8P, -1._R8P, side == 1)
                           do c=1, nc
                              oc = merge(c, 1_I4P - c, side == 1) ; cl = [i + di * oc, j + dj * oc, k + dk * oc]
                              rr = q_aux(IA_R,cl(1),cl(2),cl(3),b)
                              ue_ = 0._R8P
                              if (us > 0._R8P) ue_ = dot_product(u, q_aux(IA_U:IA_W,cl(1),cl(2),cl(3),b)) / us
                              un_ = sdr * q_aux(IA_U+d-1,cl(1),cl(2),cl(3),b)
                              rv = [rr, rr * ue_, rr * un_, rr * ue_ * un_]
                              do kk=1, 4
                                 ab = 6_I4P + WM_CONV_CELLS + 4_I4P * (c - 1_I4P) + kk
                                 rv(kk) = (1._R8P - ep) * avg(ab,ia,ja,ka,b) + ep * rv(kk)
                                 if (commit) avg(ab,ia,ja,ka,b) = rv(kk)
                              enddo
                              if (rv(1) > 0._R8P) rc(c) = rv(4) - rv(2) * rv(3) / rv(1)
                           enddo
                        endif
                        call wall_model_solve_convective(h=h, us=us, ts=ts, ps=ps, tw=tw(f), px=px, nc=nc, dyc=dxyz(d,b), cc=cc, rc=rc, ftau=rsf(f), R=R, cp=cp, pr=pr, &
                                                         mu0=mu0, tref=tref, suth=suth, omega_mu=omega_mu, tau=tau, q=q, semilocal=sl(f))
                     else
                        call wall_model_solve_pressure(h=h, us=us, ts=ts, ps=ps, tw=tw(f), px=px, R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, suth=suth, &
                                                       omega_mu=omega_mu, tau=tau, q=q, semilocal=sl(f), ypg=yp(f))
                     endif
                  elseif (inv(f)) then
                     call wall_model_solve_inverse(h=h, us=us, ts=ts, ps=ps, tw=tw(f), ue=ue(f), te=te(f), R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, suth=suth, omega_mu=omega_mu, &
                                                   tau=tau, q=q)
                  else
                     call wall_model_solve(h=h, us=us, ts=ts, ps=ps, tw=tw(f), R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, suth=suth, omega_mu=omega_mu, tau=tau, q=q, semilocal=sl(f))
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
   contains
      function sample_pressure_gradient(i, j, k, b, o1, o2, w, u, us) result(px)
      !< The pressure gradient along the sample velocity `u` (speed `us`, 0 without a speed): per active wall-parallel direction the central difference of
      !< the neighbours of the sample cells o1 and o2, weighted as the sample (issue #121, S4d).
      integer(I4P), intent(in) :: i, j, k, b, o1, o2 !< Face index, block, sample cells.
      real(R8P),    intent(in) :: w, u(3), us        !< Weight of o2, sample velocity and speed.
      real(R8P)                :: px                 !< Gradient along u.
      real(R8P)                :: g, g2              !< Gradient along a direction.
      integer(I4P)             :: t, et(3), c(3)     !< Direction, its unit step, sample cell.

      px = 0._R8P
      if (.not.(us > 0._R8P)) return
      do t=1, 3
         if (t == d .or. is_null(t)) cycle
         et = 0_I4P ; et(t) = 1_I4P
         c = [i + di * o1, j + dj * o1, k + dk * o1]
         g = (q_aux(IA_P,c(1)+et(1),c(2)+et(2),c(3)+et(3),b) - q_aux(IA_P,c(1)-et(1),c(2)-et(2),c(3)-et(3),b)) / (2._R8P * dxyz(t,b))
         if (w > 0._R8P) then
            c = [i + di * o2, j + dj * o2, k + dk * o2]
            g2 = (q_aux(IA_P,c(1)+et(1),c(2)+et(2),c(3)+et(3),b) - q_aux(IA_P,c(1)-et(1),c(2)-et(2),c(3)-et(3),b)) / (2._R8P * dxyz(t,b))
            g = g + w * (g2 - g)
         endif
         px = px + g * u(t)
      enddo
      px = px / us
      endfunction sample_pressure_gradient

      function cell_convection(i, j, k, b, c, u, us) result(cv_)
      !< The convection along the sample velocity `u` (speed `us`) at the cell `c` counted from the wall (issue #121, P2):
      !< `sum_t e_t sum_j d(rho u_t u_j)/dx_j`, `e = u / us`; along each active wall-parallel direction a central difference of the neighbours;
      !< along the normal, in the coordinate away from the wall, `(G_2 + 3 G_1) / (3 dx)` in the first cell (the wall's G = 0 at half a cell)
      !< and a central difference above, `G = rho u_t u_d`.
      integer(I4P), intent(in) :: i, j, k, b, c      !< Face index, block, cell from the wall.
      real(R8P),    intent(in) :: u(3), us           !< Sample velocity and speed.
      real(R8P)                :: cv_                !< Convection.
      integer(I4P)             :: t, jd, oc, cc0(3), et(3) !< Counters, cell.
      real(R8P)                :: sd                 !< Sign of the axis away from the wall.

      cv_ = 0._R8P
      if (.not.(us > 0._R8P)) return
      sd = merge(1._R8P, -1._R8P, side == 1)
      oc = merge(c, 1_I4P - c, side == 1)
      cc0 = [i + di * oc, j + dj * oc, k + dk * oc]
      do t=1, 3
         if (t == d .or. u(t) == 0._R8P) cycle
         do jd=1, 3
            if (jd /= d .and. is_null(jd)) cycle
            et = 0_I4P ; et(jd) = 1_I4P
            if (jd /= d) then
               cv_ = cv_ + u(t) / us * (flux(cc0 + et, b, t, jd) - flux(cc0 - et, b, t, jd)) / (2._R8P * dxyz(jd,b))
            elseif (c == 1_I4P) then
               cv_ = cv_ + u(t) / us * sd * (flux(cc0 + nint(sd) * et, b, t, jd) + 3._R8P * flux(cc0, b, t, jd)) / (3._R8P * dxyz(jd,b))
            else
               cv_ = cv_ + u(t) / us * sd * (flux(cc0 + nint(sd) * et, b, t, jd) - flux(cc0 - nint(sd) * et, b, t, jd)) / (2._R8P * dxyz(jd,b))
            endif
         enddo
      enddo
      endfunction cell_convection

      function flux(cl, b, t, jd) result(g)
      !< `rho u_t u_jd` of the cell `cl` of the block `b`.
      integer(I4P), intent(in) :: cl(3), b, t, jd !< Cell, block, components.
      real(R8P)                :: g            !< Momentum flux.

      g = q_aux(IA_R,cl(1),cl(2),cl(3),b) * q_aux(IA_U+t-1,cl(1),cl(2),cl(3),b) * q_aux(IA_U+jd-1,cl(1),cl(2),cl(3),b)
      endfunction flux
   endsubroutine overwrite_wall_fluxes
   subroutine add_wall_mu_t(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, R, pr, mu0, tref, suth, omega_mu, faces, ev, tw, hwm, sl, inv, ue, te, &
                            is_null, dxyz, q_aux, mu_t)
   !< Raise `mu_t` to the wall model's eddy viscosity in the cells below the matching height of the modelled faces of direction `d` that ask for it
   !< (`wall_model_eddy_viscosity`, issue #119): `wall_eddy_viscosity` with the `tau_w` of the instantaneous sample (never the filtered one, which the ghost
   !< columns do not hold), the damping of the face (semi-local for the inverse model), the cell's density and viscosity. The columns of the 2 ghost layers
   !< along each active tangential direction are computed as well, from their ghost cells: the same inputs as the neighbour's interior column, so the same
   !< value, and the face fluxes of the two blocks agree (the edge columns, which no flux reads, are skipped).
   integer(I4P), intent(in)    :: d                                 !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                        !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma, R, pr                      !< Specific heats ratio, gas constant, Prandtl number.
   real(R8P),    intent(in)    :: mu0, tref, suth, omega_mu               !< Viscosity law.
   integer(I4P), intent(in)    :: faces(1:,1:)                      !< Modelled block faces [nb, 6].
   logical,      intent(in)    :: ev(6)                             !< The eddy viscosity per face.
   real(R8P),    intent(in)    :: tw(6)                             !< Wall temperature per face.
   real(R8P),    intent(in)    :: hwm(6)                            !< Matching height per face, 0 for the 3rd cell.
   logical,      intent(in)    :: sl(6)                             !< Semi-local damping per face.
   logical,      intent(in)    :: inv(6)                            !< The inverse-transformation model per face.
   real(R8P),    intent(in)    :: ue(6), te(6)                      !< Inverse model: edge speed and temperature per face.
   logical,      intent(in)    :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                       !< Blocks space steps.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: mu_t(1-ngc:,1-ngc:,1-ngc:,1:)     !< Eddy viscosity.
   real(R8P)                   :: u(3), us, ts, ps, rs, tau, q, cp, h, w, rhow, muw, rho, mu, mt !< Work.
   integer(I4P)                :: n, side, f, b, i, j, k, l, m, mc, c, o1, o2, oc, gi, gj, gk, ng !< Counters.

   n = ni * di + nj * dj + nk * dk
   cp = gamma * (R / (gamma - 1._R8P))
   gi = merge(0_I4P, 2_I4P, di == 1_I4P .or. is_null(1))
   gj = merge(0_I4P, 2_I4P, dj == 1_I4P .or. is_null(2))
   gk = merge(0_I4P, 2_I4P, dk == 1_I4P .or. is_null(3))
   do side=1, 2
      f = 2_I4P * d - 2_I4P + side
      if (.not.ev(f)) cycle
      l = merge(0_I4P, n, side == 1)
      !$omp parallel do collapse(4) default(firstprivate) shared(faces, dxyz, q_aux, mu_t)
      do b=1, blocks_number
         do k=merge(l, 1_I4P - gk, dk == 1_I4P), merge(l, nk + gk, dk == 1_I4P)
            do j=merge(l, 1_I4P - gj, dj == 1_I4P), merge(l, nj + gj, dj == 1_I4P)
               do i=merge(l, 1_I4P - gi, di == 1_I4P), merge(l, ni + gi, di == 1_I4P)
                  if (faces(b,f) == 0_I4P) cycle
                  ng = 0_I4P
                  if (di == 0_I4P .and. (i < 1_I4P .or. i > ni)) ng = ng + 1_I4P
                  if (dj == 0_I4P .and. (j < 1_I4P .or. j > nj)) ng = ng + 1_I4P
                  if (dk == 0_I4P .and. (k < 1_I4P .or. k > nk)) ng = ng + 1_I4P
                  if (ng > 1_I4P) cycle
                  ! the sample, as overwrite_wall_fluxes takes it (without the filter)
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
                  us = sqrt(dot_product(u, u))
                  if (.not.(us > 0._R8P)) cycle
                  ts = ps / (R * rs)
                  if (inv(f)) then
                     call wall_model_solve_inverse(h=h, us=us, ts=ts, ps=ps, tw=tw(f), ue=ue(f), te=te(f), R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, suth=suth, omega_mu=omega_mu, &
                                                   tau=tau, q=q)
                  else
                     call wall_model_solve(h=h, us=us, ts=ts, ps=ps, tw=tw(f), R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, suth=suth, omega_mu=omega_mu, tau=tau, q=q, semilocal=sl(f))
                  endif
                  if (.not.(tau > 0._R8P)) cycle
                  rhow = ps / (R * tw(f))
                  muw = mu0 * temperature_law(t=tw(f), tref=tref, suth=suth, omega=omega_mu)
                  ! the cells strictly below h: 1 .. m, without m when the sample is its centre (w = 0)
                  mc = m - merge(1_I4P, 0_I4P, w == 0._R8P)
                  do c=1, mc
                     oc = merge(c, 1_I4P - c, side == 1)
                     rho = q_aux(IA_R,i+di*oc,j+dj*oc,k+dk*oc,b)
                     mu = mu0 * temperature_law(t=q_aux(IA_P,i+di*oc,j+dj*oc,k+dk*oc,b) / (R * rho), tref=tref, suth=suth, omega=omega_mu)
                     mt = wall_eddy_viscosity(y=(real(c, R8P) - 0.5_R8P) * dxyz(d,b), rho=rho, mu=mu, tau=tau, rhow=rhow, muw=muw, semilocal=sl(f) .or. inv(f))
                     mu_t(i+di*oc,j+dj*oc,k+dk*oc,b) = max(mu_t(i+di*oc,j+dj*oc,k+dk*oc,b), mt)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end parallel do
   enddo
   endsubroutine add_wall_mu_t

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
