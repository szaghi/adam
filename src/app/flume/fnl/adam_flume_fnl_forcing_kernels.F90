!< ADAM, FLUME FNL kernels of the mass-flow forcing (issue #119, T3.1).

#include "fundal.H"

module adam_flume_fnl_forcing_kernels
!< ADAM, FLUME FNL kernels of the mass-flow forcing (issue #119, T3.1).

use :: adam_flume_dissipation_library, only : temperature_law
use :: adam_flume_parameters, only : IA_P, IA_R, IA_U, IA_V, IA_W, IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW
use :: adam_flume_wall_model_library, only : sample_cells, wall_eddy_viscosity, wall_model_solve, wall_model_solve_inverse, wall_model_solve_pressure
use :: penf,                  only : I4P, R8P

implicit none
private
public :: add_uniform_momentum_dev
public :: add_wall_mu_t_dev
public :: zero_mu_t_dev
public :: overwrite_wall_fluxes_dev
public :: close_wall_fluxes_dev
public :: sum_momentum_dev

contains
   subroutine sum_momentum_dev(ni, nj, nk, ngc, blocks_number, axis, dxyz_gpu, q_gpu, sums)
   !< Sum the momentum along `axis` and the volume over the interior cells.
   integer(I4P), intent(in)  :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)  :: blocks_number                     !< Actual blocks number.
   integer(I4P), intent(in)  :: axis                              !< Axis of the flow.
   real(R8P),    intent(in)  :: dxyz_gpu(1:,1:)                   !< Blocks space steps [nb, 3].
   real(R8P),    intent(in)  :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: sums(2)                           !< Momentum times volume, volume.
   real(R8P)                 :: s1, s2                            !< Reduction accumulators.
   integer(I4P)              :: iv                                !< Momentum slot.
   integer(I4P)              :: b, i, j, k                        !< Counters.

   s1 = 0._R8P ; s2 = 0._R8P ; iv = IQ_RU + axis - 1_I4P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_gpu) &
   !$acc& firstprivate(ni,nj,nk,blocks_number,iv) reduction(+:s1,s2)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_gpu) &
   !$omp& firstprivate(ni,nj,nk,blocks_number,iv) reduction(+:s1,s2)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      s1 = s1 + q_gpu(b,i,j,k,iv) * dxyz_gpu(b,1) * dxyz_gpu(b,2) * dxyz_gpu(b,3)
      s2 = s2 + dxyz_gpu(b,1) * dxyz_gpu(b,2) * dxyz_gpu(b,3)
   enddo
   enddo
   enddo
   enddo
   sums = [s1, s2]
   endsubroutine sum_momentum_dev

   subroutine add_uniform_momentum_dev(ni, nj, nk, ngc, blocks_number, axis, dm, q_gpu)
   !< Add `dm` to the momentum along `axis`, and to the total energy the change of kinetic energy.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   integer(I4P), intent(in)    :: axis                              !< Axis of the flow.
   real(R8P),    intent(in)    :: dm                                !< Momentum added.
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   integer(I4P)                :: iv                                !< Momentum slot.
   integer(I4P)                :: b, i, j, k                        !< Counters.

   iv = IQ_RU + axis - 1_I4P
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(q_gpu) firstprivate(ni,nj,nk,blocks_number,iv,dm)
   !$omp OMPLOOP collapse(4) DEVICEPTR(q_gpu) firstprivate(ni,nj,nk,blocks_number,iv,dm)
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      q_gpu(b,i,j,k,IQ_RE) = q_gpu(b,i,j,k,IQ_RE) + dm * (q_gpu(b,i,j,k,iv) + 0.5_R8P * dm) / q_gpu(b,i,j,k,IQ_R)
      q_gpu(b,i,j,k,iv) = q_gpu(b,i,j,k,iv) + dm
   enddo
   enddo
   enddo
   enddo
   endsubroutine add_uniform_momentum_dev
   subroutine overwrite_wall_fluxes_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, R, pr, mu0, tref, suth, omega_mu, faces_gpu, tw, hwm, sl, tf, dt, commit, &
                                        fresh, avg_gpu, inv, ue, te, pg, tp, is_null, dxyz_gpu, q_aux_gpu, fl_gpu)
   !< Overwrite the flux through the modelled wall faces of direction `d` (issue #119, T3.3): no mass, the parallel momentum by the modelled shear, the energy by the heat
   !< flux; the sample is the 3rd cell off the wall. The normal momentum keeps the mirrored inviscid flux, the pressure.
   !< With a filter time `tf` the sample is `(1 - e) a + e s`, `e = min(dt / tf, 1)` (1 when `fresh`), from the averages `a` of the face, which take it when `commit`.
   integer(I4P), intent(in)    :: d                               !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                      !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                 !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                   !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma, R, pr                    !< Specific heats ratio, gas constant, Prandtl number.
   real(R8P),    intent(in)    :: mu0, tref, suth, omega_mu             !< Viscosity law.
   integer(I4P), intent(in)    :: faces_gpu(1:,1:)                   !< Modelled block faces [nb, 6].
   real(R8P),    intent(in)    :: tw(6)                           !< Wall temperature per face.
   real(R8P),    intent(in)    :: hwm(6)                          !< Matching height per face, 0 for the 3rd cell.
   logical,      intent(in)    :: sl(6)                           !< Semi-local damping per face.
   real(R8P),    intent(in)    :: tf(6)                           !< Filter time per face, 0 for none.
   real(R8P),    intent(in)    :: dt                              !< Time step.
   logical,      intent(in)    :: commit                          !< The averages take this sample.
   logical,      intent(in)    :: fresh                           !< The averages hold no history.
   real(R8P),    intent(inout) :: avg_gpu(1:,1:,1:,1:,1:)         !< Averages of the faces of direction `d` [nb, ni or 2, nj or 2, nk or 2, 5].
   logical,      intent(in)    :: inv(6)                          !< The inverse-transformation model per face (T3.5), else the equilibrium one.
   real(R8P),    intent(in)    :: ue(6), te(6)                    !< Inverse model: speed and temperature at the edge of the layer per face.
   logical,      intent(in)    :: pg(6)                           !< The pressure gradient in the equilibrium model per face (issue #121, S4d).
   real(R8P),    intent(in)    :: tp(6)                           !< Filter time of the pressure gradient per face, 0 for none.
   logical,      intent(in)    :: is_null(3)                      !< Null directions.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                    !< Blocks space steps.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)      !< Face fluxes of direction `d`.
   real(R8P)                   :: u1, u2, u3, us, ts, ps, rs, tau, q, cp, h, s, t1, t2, t3, w, e !< Work.
   integer(I4P)                :: n, side, f, b, i, j, k, l, m, o1, o2, i0, i1, j0, j1, k0, k1, ia, ja, ka !< Counters.
   real(R8P)                   :: twf, hwf, tff, uef, tef            !< Wall temperature, matching height, filter time, edge speed and temperature of the face.
   logical                     :: slf, ivf                           !< Semi-local damping, inverse model of the face.
   logical                     :: pgf, nl1, nl2, nl3                 !< Pressure gradient of the face, null directions.
   real(R8P)                   :: px, g, g2, tpf, ep                 !< Pressure gradient along the sample velocity, along a direction; its filter.
   integer(I4P)                :: t, ei, ej, ek                      !< Wall-parallel direction, its unit step.

   n = ni * di + nj * dj + nk * dk
   cp = gamma * (R / (gamma - 1._R8P)) ! as the physics object derives it (cv = R / (gamma - 1), cp = gamma cv)
   do side=1, 2
      f = 2_I4P * d - 2_I4P + side
      l = merge(0_I4P, n, side == 1)
      s = merge(1._R8P, -1._R8P, side == 1)
      twf = tw(f) ; hwf = hwm(f) ; slf = sl(f) ; tff = tf(f) ; ivf = inv(f) ; uef = ue(f) ; tef = te(f)
      pgf = pg(f) ; nl1 = is_null(1) ; nl2 = is_null(2) ; nl3 = is_null(3) ; tpf = tp(f)
      ep = 1._R8P
      if (tpf > 0._R8P .and. .not.fresh) ep = min(dt / tpf, 1._R8P)
      e = 1._R8P
      if (tff > 0._R8P .and. .not.fresh) e = min(dt / tff, 1._R8P)
      i0 = merge(l, 1_I4P, di == 1_I4P) ; i1 = merge(l, ni, di == 1_I4P)
      j0 = merge(l, 1_I4P, dj == 1_I4P) ; j1 = merge(l, nj, dj == 1_I4P)
      k0 = merge(l, 1_I4P, dk == 1_I4P) ; k1 = merge(l, nk, dk == 1_I4P)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(faces_gpu,dxyz_gpu,q_aux_gpu,fl_gpu,avg_gpu)                                                 &
      !$acc& firstprivate(d,di,dj,dk,f,n,side,s,twf,hwf,slf,tff,e,commit,ivf,uef,tef,pgf,nl1,nl2,nl3,tpf,ep,R,cp,pr,mu0,tref,suth, omega_mu,blocks_number,i0,i1,j0,j1,k0,k1)     &
      !$acc& private(u1,u2,u3,us,ts,ps,rs,tau,q,h,t1,t2,t3,w,m,o1,o2,ia,ja,ka,px,g,g2,t,ei,ej,ek)
      !$omp OMPLOOP collapse(4) DEVICEPTR(faces_gpu,dxyz_gpu,q_aux_gpu,fl_gpu,avg_gpu)                                                                &
      !$omp& firstprivate(d,di,dj,dk,f,n,side,s,twf,hwf,slf,tff,e,commit,ivf,uef,tef,pgf,nl1,nl2,nl3,tpf,ep,R,cp,pr,mu0,tref,suth, omega_mu,blocks_number,i0,i1,j0,j1,k0,k1)     &
      !$omp& private(u1,u2,u3,us,ts,ps,rs,tau,q,h,t1,t2,t3,w,m,o1,o2,ia,ja,ka,px,g,g2,t,ei,ej,ek)
      do k=k0, k1
      do j=j0, j1
      do i=i0, i1
      do b=1, blocks_number
         if (faces_gpu(b,f) /= 0_I4P) then
            call sample_cells(hwm=hwf, dx=dxyz_gpu(b,d), n=n, m=m, w=w, h=h)
            o1 = merge(m, 1_I4P - m, side == 1) ; o2 = merge(m + 1_I4P, -m, side == 1)
            u1 = q_aux_gpu(b,i+di*o1,j+dj*o1,k+dk*o1,IA_U)
            u2 = q_aux_gpu(b,i+di*o1,j+dj*o1,k+dk*o1,IA_V)
            u3 = q_aux_gpu(b,i+di*o1,j+dj*o1,k+dk*o1,IA_W)
            ps = q_aux_gpu(b,i+di*o1,j+dj*o1,k+dk*o1,IA_P)
            rs = q_aux_gpu(b,i+di*o1,j+dj*o1,k+dk*o1,IA_R)
            if (w > 0._R8P) then
               u1 = u1 + w * (q_aux_gpu(b,i+di*o2,j+dj*o2,k+dk*o2,IA_U) - u1)
               u2 = u2 + w * (q_aux_gpu(b,i+di*o2,j+dj*o2,k+dk*o2,IA_V) - u2)
               u3 = u3 + w * (q_aux_gpu(b,i+di*o2,j+dj*o2,k+dk*o2,IA_W) - u3)
               ps = ps + w * (q_aux_gpu(b,i+di*o2,j+dj*o2,k+dk*o2,IA_P) - ps)
               rs = rs + w * (q_aux_gpu(b,i+di*o2,j+dj*o2,k+dk*o2,IA_R) - rs)
            endif
            u1 = u1 * (1._R8P - real(di, R8P)) ; u2 = u2 * (1._R8P - real(dj, R8P)) ; u3 = u3 * (1._R8P - real(dk, R8P))
            if (tff > 0._R8P) then ! the input filter; e = 1 takes the sample exactly
               ia = merge(side, i, di == 1_I4P) ; ja = merge(side, j, dj == 1_I4P) ; ka = merge(side, k, dk == 1_I4P)
               u1 = (1._R8P - e) * avg_gpu(b,ia,ja,ka,1) + e * u1
               u2 = (1._R8P - e) * avg_gpu(b,ia,ja,ka,2) + e * u2
               u3 = (1._R8P - e) * avg_gpu(b,ia,ja,ka,3) + e * u3
               ps = (1._R8P - e) * avg_gpu(b,ia,ja,ka,4) + e * ps
               rs = (1._R8P - e) * avg_gpu(b,ia,ja,ka,5) + e * rs
               if (commit) then
                  avg_gpu(b,ia,ja,ka,1) = u1 ; avg_gpu(b,ia,ja,ka,2) = u2 ; avg_gpu(b,ia,ja,ka,3) = u3
                  avg_gpu(b,ia,ja,ka,4) = ps ; avg_gpu(b,ia,ja,ka,5) = rs
               endif
            endif
            us = sqrt(u1 * u1 + u2 * u2 + u3 * u3)
            ts = ps / (R * rs)
            px = 0._R8P
            if (pgf .and. us > 0._R8P) then ! the pressure gradient along the sample velocity (issue #121, S4d), as the CPU kernel takes it
               do t=1, 3
                  if (t /= d .and. .not.merge(nl1, merge(nl2, nl3, t == 2_I4P), t == 1_I4P)) then
                     ei = merge(1_I4P, 0_I4P, t == 1_I4P) ; ej = merge(1_I4P, 0_I4P, t == 2_I4P) ; ek = merge(1_I4P, 0_I4P, t == 3_I4P)
                     g = (q_aux_gpu(b,i+di*o1+ei,j+dj*o1+ej,k+dk*o1+ek,IA_P) - q_aux_gpu(b,i+di*o1-ei,j+dj*o1-ej,k+dk*o1-ek,IA_P)) / (2._R8P * dxyz_gpu(b,t))
                     if (w > 0._R8P) then
                        g2 = (q_aux_gpu(b,i+di*o2+ei,j+dj*o2+ej,k+dk*o2+ek,IA_P) - q_aux_gpu(b,i+di*o2-ei,j+dj*o2-ej,k+dk*o2-ek,IA_P)) / (2._R8P * dxyz_gpu(b,t))
                        g = g + w * (g2 - g)
                     endif
                     px = px + g * merge(u1, merge(u2, u3, t == 2_I4P), t == 1_I4P)
                  endif
               enddo
               px = px / us
            endif
            if (pgf .and. tpf > 0._R8P) then ! its own exponential filter (the turbulent gradient is as large as tau / h): avg(6)
               ia = merge(side, i, di == 1_I4P) ; ja = merge(side, j, dj == 1_I4P) ; ka = merge(side, k, dk == 1_I4P)
               px = (1._R8P - ep) * avg_gpu(b,ia,ja,ka,6) + ep * px
               if (commit) avg_gpu(b,ia,ja,ka,6) = px
            endif
            if (pgf) then
               call wall_model_solve_pressure(h=h, us=us, ts=ts, ps=ps, tw=twf, px=px, R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, suth=suth, omega_mu=omega_mu, &
                                              tau=tau, q=q, semilocal=slf)
            elseif (ivf) then
               call wall_model_solve_inverse(h=h, us=us, ts=ts, ps=ps, tw=twf, ue=uef, te=tef, R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, suth=suth, omega_mu=omega_mu, tau=tau, q=q)
            else
               call wall_model_solve(h=h, us=us, ts=ts, ps=ps, tw=twf, R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, suth=suth, omega_mu=omega_mu, tau=tau, q=q, semilocal=slf)
            endif
            t1 = 0._R8P ; t2 = 0._R8P ; t3 = 0._R8P
            if (us > 0._R8P) then
               t1 = u1 / us ; t2 = u2 / us ; t3 = u3 / us
            endif
            fl_gpu(b,i,j,k,IQ_R) = 0._R8P
            if (di == 0_I4P) fl_gpu(b,i,j,k,IQ_RU) = -s * tau * t1
            if (dj == 0_I4P) fl_gpu(b,i,j,k,IQ_RV) = -s * tau * t2
            if (dk == 0_I4P) fl_gpu(b,i,j,k,IQ_RW) = -s * tau * t3
            fl_gpu(b,i,j,k,IQ_RE) = -s * q
         endif
      enddo
      enddo
      enddo
      enddo
   enddo
   endsubroutine overwrite_wall_fluxes_dev
   subroutine zero_mu_t_dev(ni, nj, nk, ngc, blocks_number, mu_t_gpu)
   !< Zero the eddy viscosity, ghost cells included: the wall model's eddy viscosity without a subgrid-scale model raises a zero field (issue #119).
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                      !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                        !< Actual blocks number.
   real(R8P),    intent(inout) :: mu_t_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Eddy viscosity, [nb, i, j, k, 1].
   integer(I4P)                :: b, i, j, k                           !< Counters.

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(mu_t_gpu) firstprivate(ni,nj,nk,ngc,blocks_number)
   !$omp OMPLOOP collapse(4) DEVICEPTR(mu_t_gpu) firstprivate(ni,nj,nk,ngc,blocks_number)
   do k=1-ngc, nk+ngc
   do j=1-ngc, nj+ngc
   do i=1-ngc, ni+ngc
   do b=1, blocks_number
      mu_t_gpu(b,i,j,k,1) = 0._R8P
   enddo
   enddo
   enddo
   enddo
   endsubroutine zero_mu_t_dev

   subroutine add_wall_mu_t_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, R, pr, mu0, tref, suth, omega_mu, faces_gpu, ev, tw, hwm, sl, inv, ue, te, &
                                is_null, dxyz_gpu, q_aux_gpu, mu_t_gpu)
   !< Device twin of `add_wall_mu_t` (issue #119, `wall_model_eddy_viscosity`).
   integer(I4P), intent(in)    :: d                                     !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                            !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma, R, pr                          !< Specific heats ratio, gas constant, Prandtl number.
   real(R8P),    intent(in)    :: mu0, tref, suth, omega_mu                   !< Viscosity law.
   integer(I4P), intent(in)    :: faces_gpu(1:,1:)                      !< Modelled block faces [nb, 6].
   logical,      intent(in)    :: ev(6)                                 !< The eddy viscosity per face.
   real(R8P),    intent(in)    :: tw(6)                                 !< Wall temperature per face.
   real(R8P),    intent(in)    :: hwm(6)                                !< Matching height per face, 0 for the 3rd cell.
   logical,      intent(in)    :: sl(6)                                 !< Semi-local damping per face.
   logical,      intent(in)    :: inv(6)                                !< The inverse-transformation model per face.
   real(R8P),    intent(in)    :: ue(6), te(6)                          !< Inverse model: edge speed and temperature per face.
   logical,      intent(in)    :: is_null(3)                            !< Null directions.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                       !< Blocks space steps [nb, 3].
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: mu_t_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Eddy viscosity, [nb, i, j, k, 1].
   real(R8P)                   :: u1, u2, u3, us, ts, ps, rs, tau, q, cp, h, w, rhow, muw, rho, mu, mt, dx !< Work.
   real(R8P)                   :: twf, hwf, uef, tef                    !< Wall temperature, matching height, edge speed and temperature of the face.
   logical                     :: slf, ivf, dlf                         !< Semi-local damping, inverse model, the damping of the eddy viscosity of the face.
   integer(I4P)                :: n, side, f, b, i, j, k, l, m, mc, c, o1, o2, oc, gi, gj, gk, ng, i0, i1, j0, j1, k0, k1 !< Counters.

   n = ni * di + nj * dj + nk * dk
   cp = gamma * (R / (gamma - 1._R8P))
   gi = merge(0_I4P, 2_I4P, di == 1_I4P .or. is_null(1))
   gj = merge(0_I4P, 2_I4P, dj == 1_I4P .or. is_null(2))
   gk = merge(0_I4P, 2_I4P, dk == 1_I4P .or. is_null(3))
   do side=1, 2
      f = 2_I4P * d - 2_I4P + side
      if (.not.ev(f)) cycle
      l = merge(0_I4P, n, side == 1)
      twf = tw(f) ; hwf = hwm(f) ; slf = sl(f) ; ivf = inv(f) ; uef = ue(f) ; tef = te(f) ; dlf = sl(f) .or. inv(f)
      i0 = merge(l, 1_I4P - gi, di == 1_I4P) ; i1 = merge(l, ni + gi, di == 1_I4P)
      j0 = merge(l, 1_I4P - gj, dj == 1_I4P) ; j1 = merge(l, nj + gj, dj == 1_I4P)
      k0 = merge(l, 1_I4P - gk, dk == 1_I4P) ; k1 = merge(l, nk + gk, dk == 1_I4P)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(faces_gpu,dxyz_gpu,q_aux_gpu,mu_t_gpu)                                          &
      !$acc& firstprivate(d,di,dj,dk,ni,nj,nk,f,n,side,twf,hwf,slf,ivf,uef,tef,dlf,R,cp,pr,mu0,tref,suth, omega_mu,blocks_number,i0,i1,j0,j1,k0,k1)       &
      !$acc& private(u1,u2,u3,us,ts,ps,rs,tau,q,h,w,rhow,muw,rho,mu,mt,dx,m,mc,c,o1,o2,oc,ng)
      !$omp OMPLOOP collapse(4) DEVICEPTR(faces_gpu,dxyz_gpu,q_aux_gpu,mu_t_gpu)                                                                    &
      !$omp& firstprivate(d,di,dj,dk,ni,nj,nk,f,n,side,twf,hwf,slf,ivf,uef,tef,dlf,R,cp,pr,mu0,tref,suth, omega_mu,blocks_number,i0,i1,j0,j1,k0,k1)       &
      !$omp& private(u1,u2,u3,us,ts,ps,rs,tau,q,h,w,rhow,muw,rho,mu,mt,dx,m,mc,c,o1,o2,oc,ng)
      do k=k0, k1
      do j=j0, j1
      do i=i0, i1
      do b=1, blocks_number
         ng = 0_I4P
         if (di == 0_I4P .and. (i < 1_I4P .or. i > ni)) ng = ng + 1_I4P
         if (dj == 0_I4P .and. (j < 1_I4P .or. j > nj)) ng = ng + 1_I4P
         if (dk == 0_I4P .and. (k < 1_I4P .or. k > nk)) ng = ng + 1_I4P
         if (faces_gpu(b,f) /= 0_I4P .and. ng <= 1_I4P) then
            dx = dxyz_gpu(b,d)
            call sample_cells(hwm=hwf, dx=dx, n=n, m=m, w=w, h=h)
            o1 = merge(m, 1_I4P - m, side == 1) ; o2 = merge(m + 1_I4P, -m, side == 1)
            u1 = q_aux_gpu(b,i+di*o1,j+dj*o1,k+dk*o1,IA_U)
            u2 = q_aux_gpu(b,i+di*o1,j+dj*o1,k+dk*o1,IA_V)
            u3 = q_aux_gpu(b,i+di*o1,j+dj*o1,k+dk*o1,IA_W)
            ps = q_aux_gpu(b,i+di*o1,j+dj*o1,k+dk*o1,IA_P)
            rs = q_aux_gpu(b,i+di*o1,j+dj*o1,k+dk*o1,IA_R)
            if (w > 0._R8P) then
               u1 = u1 + w * (q_aux_gpu(b,i+di*o2,j+dj*o2,k+dk*o2,IA_U) - u1)
               u2 = u2 + w * (q_aux_gpu(b,i+di*o2,j+dj*o2,k+dk*o2,IA_V) - u2)
               u3 = u3 + w * (q_aux_gpu(b,i+di*o2,j+dj*o2,k+dk*o2,IA_W) - u3)
               ps = ps + w * (q_aux_gpu(b,i+di*o2,j+dj*o2,k+dk*o2,IA_P) - ps)
               rs = rs + w * (q_aux_gpu(b,i+di*o2,j+dj*o2,k+dk*o2,IA_R) - rs)
            endif
            u1 = u1 * (1._R8P - real(di, R8P)) ; u2 = u2 * (1._R8P - real(dj, R8P)) ; u3 = u3 * (1._R8P - real(dk, R8P))
            us = sqrt(u1 * u1 + u2 * u2 + u3 * u3)
            tau = 0._R8P
            if (us > 0._R8P) then
               ts = ps / (R * rs)
               if (ivf) then
                  call wall_model_solve_inverse(h=h, us=us, ts=ts, ps=ps, tw=twf, ue=uef, te=tef, R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, suth=suth, omega_mu=omega_mu, tau=tau, q=q)
               else
                  call wall_model_solve(h=h, us=us, ts=ts, ps=ps, tw=twf, R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, suth=suth, omega_mu=omega_mu, tau=tau, q=q, semilocal=slf)
               endif
            endif
            if (tau > 0._R8P) then
               rhow = ps / (R * twf)
               muw = mu0 * temperature_law(t=twf, tref=tref, suth=suth, omega=omega_mu)
               mc = m - merge(1_I4P, 0_I4P, w == 0._R8P)
               do c=1, mc
                  oc = merge(c, 1_I4P - c, side == 1)
                  rho = q_aux_gpu(b,i+di*oc,j+dj*oc,k+dk*oc,IA_R)
                  mu = mu0 * temperature_law(t=q_aux_gpu(b,i+di*oc,j+dj*oc,k+dk*oc,IA_P) / (R * rho), tref=tref, suth=suth, omega=omega_mu)
                  mt = wall_eddy_viscosity(y=(real(c, R8P) - 0.5_R8P) * dx, rho=rho, mu=mu, tau=tau, rhow=rhow, muw=muw, semilocal=dlf)
                  mu_t_gpu(b,i+di*oc,j+dj*oc,k+dk*oc,1) = max(mu_t_gpu(b,i+di*oc,j+dj*oc,k+dk*oc,1), mt)
               enddo
            endif
         endif
      enddo
      enddo
      enddo
      enddo
   enddo
   endsubroutine add_wall_mu_t_dev

   subroutine close_wall_fluxes_dev(d, di, dj, dk, ni, nj, nk, nv, blocks_number, faces_gpu, fl_gpu)
   !< Zero every flux of direction `d` through the flagged block faces of that direction but the normal momentum (issue #119, the #65 mass leak).
   integer(I4P), intent(in)    :: d                               !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                      !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, nv                  !< Grid dimensions, variables.
   integer(I4P), intent(in)    :: blocks_number                   !< Actual blocks number.
   integer(I4P), intent(in)    :: faces_gpu(1:,1:)                !< Flagged block faces [nb, 6].
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:) !< Face fluxes of direction `d`.
   integer(I4P)                :: n, side, f, b, i, j, k, l, v, vn, i0, i1, j0, j1, k0, k1 !< Counters.

   n = ni * di + nj * dj + nk * dk
   vn = IQ_RU + d - 1_I4P
   do side=1, 2
      f = 2_I4P * d - 2_I4P + side
      l = merge(0_I4P, n, side == 1)
      i0 = merge(l, 1_I4P, di == 1_I4P) ; i1 = merge(l, ni, di == 1_I4P)
      j0 = merge(l, 1_I4P, dj == 1_I4P) ; j1 = merge(l, nj, dj == 1_I4P)
      k0 = merge(l, 1_I4P, dk == 1_I4P) ; k1 = merge(l, nk, dk == 1_I4P)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(faces_gpu,fl_gpu) firstprivate(f,vn,nv,blocks_number,i0,i1,j0,j1,k0,k1)
      !$omp OMPLOOP collapse(4) DEVICEPTR(faces_gpu,fl_gpu) firstprivate(f,vn,nv,blocks_number,i0,i1,j0,j1,k0,k1)
      do k=k0, k1
      do j=j0, j1
      do i=i0, i1
      do b=1, blocks_number
         if (faces_gpu(b,f) /= 0_I4P) then
            do v=1, nv
               if (v /= vn) fl_gpu(b,i,j,k,v) = 0._R8P
            enddo
         endif
      enddo
      enddo
      enddo
      enddo
   enddo
   endsubroutine close_wall_fluxes_dev
endmodule adam_flume_fnl_forcing_kernels
