!< ADAM, FLUME FNL kernels of the mass-flow forcing (issue #119, T3.1).

#include "fundal.H"

module adam_flume_fnl_forcing_kernels
!< ADAM, FLUME FNL kernels of the mass-flow forcing (issue #119, T3.1).

use :: adam_flume_parameters, only : IA_P, IA_R, IA_U, IA_V, IA_W, IQ_R, IQ_RE, IQ_RU, IQ_RV, IQ_RW
use :: adam_flume_wall_model_library, only : wall_model_solve
use :: penf,                  only : I4P, R8P

implicit none
private
public :: add_uniform_momentum_dev
public :: overwrite_wall_fluxes_dev
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
   subroutine overwrite_wall_fluxes_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, gamma, R, pr, mu0, tref, omega_mu, faces_gpu, tw, dxyz_gpu, q_aux_gpu, fl_gpu)
   !< Overwrite the flux through the modelled wall faces of direction `d` (issue #119, T3.3): no mass, the parallel momentum by the modelled shear, the energy by the heat
   !< flux; the sample is the 3rd cell off the wall. The normal momentum keeps the mirrored inviscid flux, the pressure.
   integer(I4P), intent(in)    :: d                               !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                      !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                 !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                   !< Actual blocks number.
   real(R8P),    intent(in)    :: gamma, R, pr                    !< Specific heats ratio, gas constant, Prandtl number.
   real(R8P),    intent(in)    :: mu0, tref, omega_mu             !< Viscosity law.
   integer(I4P), intent(in)    :: faces_gpu(1:,1:)                   !< Modelled block faces [nb, 6].
   real(R8P),    intent(in)    :: tw(6)                           !< Wall temperature per face.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                    !< Blocks space steps.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)      !< Face fluxes of direction `d`.
   real(R8P)                   :: u1, u2, u3, us, ts, ps, tau, q, cp, h, s, t1, t2, t3 !< Work.
   integer(I4P)                :: n, side, f, b, i, j, k, l, off, i0, i1, j0, j1, k0, k1 !< Counters.
   real(R8P)                   :: twf                                !< Wall temperature of the face.

   n = ni * di + nj * dj + nk * dk
   cp = gamma * R / (gamma - 1._R8P)
   do side=1, 2
      f = 2_I4P * d - 2_I4P + side
      l = merge(0_I4P, n, side == 1)
      s = merge(1._R8P, -1._R8P, side == 1)
      off = merge(3_I4P, -2_I4P, side == 1)
      twf = tw(f)
      i0 = merge(l, 1_I4P, di == 1_I4P) ; i1 = merge(l, ni, di == 1_I4P)
      j0 = merge(l, 1_I4P, dj == 1_I4P) ; j1 = merge(l, nj, dj == 1_I4P)
      k0 = merge(l, 1_I4P, dk == 1_I4P) ; k1 = merge(l, nk, dk == 1_I4P)
      !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(faces_gpu,dxyz_gpu,q_aux_gpu,fl_gpu)      &
      !$acc& firstprivate(d,di,dj,dk,f,off,s,twf,R,cp,pr,mu0,tref,omega_mu,blocks_number,i0,i1,j0,j1,k0,k1) private(u1,u2,u3,us,ts,ps,tau,q,h,t1,t2,t3)
      !$omp OMPLOOP collapse(4) DEVICEPTR(faces_gpu,dxyz_gpu,q_aux_gpu,fl_gpu)                                     &
      !$omp& firstprivate(d,di,dj,dk,f,off,s,twf,R,cp,pr,mu0,tref,omega_mu,blocks_number,i0,i1,j0,j1,k0,k1) private(u1,u2,u3,us,ts,ps,tau,q,h,t1,t2,t3)
      do k=k0, k1
      do j=j0, j1
      do i=i0, i1
      do b=1, blocks_number
         if (faces_gpu(b,f) /= 0_I4P) then
            u1 = q_aux_gpu(b,i+di*off,j+dj*off,k+dk*off,IA_U) * (1._R8P - real(di, R8P))
            u2 = q_aux_gpu(b,i+di*off,j+dj*off,k+dk*off,IA_V) * (1._R8P - real(dj, R8P))
            u3 = q_aux_gpu(b,i+di*off,j+dj*off,k+dk*off,IA_W) * (1._R8P - real(dk, R8P))
            us = sqrt(u1 * u1 + u2 * u2 + u3 * u3)
            ps = q_aux_gpu(b,i+di*off,j+dj*off,k+dk*off,IA_P)
            ts = ps / (R * q_aux_gpu(b,i+di*off,j+dj*off,k+dk*off,IA_R))
            h = 2.5_R8P * dxyz_gpu(b,d)
            call wall_model_solve(h=h, us=us, ts=ts, ps=ps, tw=twf, R=R, cp=cp, pr=pr, mu0=mu0, tref=tref, omega_mu=omega_mu, tau=tau, q=q)
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
endmodule adam_flume_fnl_forcing_kernels
