!< ADAM, FLUME FNL kernels of the mass-flow forcing (issue #119, T3.1).

#include "fundal.H"

module adam_flume_fnl_forcing_kernels
!< ADAM, FLUME FNL kernels of the mass-flow forcing (issue #119, T3.1).

use :: adam_flume_parameters, only : IA_P, IA_R, IA_U, IA_V, IA_W, IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW
use :: adam_flume_wall_model_library, only : sample_cells, wall_model_solve, wall_model_solve_inverse
use :: penf,                  only : I4P, R8P

implicit none
private
public :: add_uniform_momentum_dev
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
   subroutine overwrite_wall_fluxes_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, R, pr, mu0, tref, omega_mu, faces_gpu, tw, hwm, sl, tf, dt, commit, &
                                        fresh, avg_gpu, inv, ue, te, dxyz_gpu, q_aux_gpu, fl_gpu)
   !< Overwrite the flux through the modelled wall faces of direction `d` (issue #119, T3.3): no mass, the parallel momentum by the modelled shear, the energy by the heat
   !< flux; the sample is the 3rd cell off the wall. The normal momentum keeps the mirrored inviscid flux, the pressure.
   !< With a filter time `tf` the sample is `(1 - e) a + e s`, `e = min(dt / tf, 1)` (1 when `fresh`), from the averages `a` of the face, which take it when `commit`.
   integer(I4P), intent(in)    :: d                               !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                      !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                 !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                   !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma, R, pr                    !< Specific heats ratio, gas constant, Prandtl number.
   real(R8P),    intent(in)    :: mu0, tref, omega_mu             !< Viscosity law.
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
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                    !< Blocks space steps.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)      !< Face fluxes of direction `d`.
   real(R8P)                   :: u1, u2, u3, us, ts, ps, rs, tau, q, cp, h, s, t1, t2, t3, w, e !< Work.
   integer(I4P)                :: n, side, f, b, i, j, k, l, m, o1, o2, i0, i1, j0, j1, k0, k1, ia, ja, ka !< Counters.
   real(R8P)                   :: twf, hwf, tff, uef, tef            !< Wall temperature, matching height, filter time, edge speed and temperature of the face.
   logical                     :: slf, ivf                           !< Semi-local damping, inverse model of the face.

   n = ni * di + nj * dj + nk * dk
   cp = gamma * (R / (gamma - 1._R8P)) ! as the physics object derives it (cv = R / (gamma - 1), cp = gamma cv)
   do side=1, 2
      f = 2_I4P * d - 2_I4P + side
      l = merge(0_I4P, n, side == 1)
      s = merge(1._R8P, -1._R8P, side == 1)
      twf = tw(f) ; hwf = hwm(f) ; slf = sl(f) ; tff = tf(f) ; ivf = inv(f) ; uef = ue(f) ; tef = te(f)
      e = 1._R8P
      if (tff > 0._R8P .and. .not.fresh) e = min(dt / tff, 1._R8P)
      i0 = merge(l, 1_I4P, di == 1_I4P) ; i1 = merge(l, ni, di == 1_I4P)
      j0 = merge(l, 1_I4P, dj == 1_I4P) ; j1 = merge(l, nj, dj == 1_I4P)
      k0 = merge(l, 1_I4P, dk == 1_I4P) ; k1 = merge(l, nk, dk == 1_I4P)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(faces_gpu,dxyz_gpu,q_aux_gpu,fl_gpu,avg_gpu)                                                 &
      !$acc& firstprivate(d,di,dj,dk,f,n,side,s,twf,hwf,slf,tff,e,commit,ivf,uef,tef,R,cp,pr,mu0,tref,omega_mu,blocks_number,i0,i1,j0,j1,k0,k1)                              &
      !$acc& private(u1,u2,u3,us,ts,ps,rs,tau,q,h,t1,t2,t3,w,m,o1,o2,ia,ja,ka)
      !$omp OMPLOOP collapse(4) DEVICEPTR(faces_gpu,dxyz_gpu,q_aux_gpu,fl_gpu,avg_gpu)                                                                &
      !$omp& firstprivate(d,di,dj,dk,f,n,side,s,twf,hwf,slf,tff,e,commit,ivf,uef,tef,R,cp,pr,mu0,tref,omega_mu,blocks_number,i0,i1,j0,j1,k0,k1)                              &
      !$omp& private(u1,u2,u3,us,ts,ps,rs,tau,q,h,t1,t2,t3,w,m,o1,o2,ia,ja,ka)
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
            if (ivf) then
               call wall_model_solve_inverse(h=h, us=us, ts=ts, ps=ps, tw=twf, ue=uef, te=tef, R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, omega_mu=omega_mu, tau=tau, q=q)
            else
               call wall_model_solve(h=h, us=us, ts=ts, ps=ps, tw=twf, R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, omega_mu=omega_mu, tau=tau, q=q, semilocal=slf)
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
