!< ADAM, FLUME FNL device kernels of the dissipative face fluxes (issue #65, M4).

#include "fundal.H"

module adam_flume_fnl_dissipation_kernels
!< ADAM, FLUME FNL device kernels of the dissipative face fluxes (issue #65, M4).
!<
!< Device twins of `adam_flume_cpu_dissipation_kernels` (same stencils, same summation order, the pointwise flux from the
!< shared `adam_flume_dissipation_library`): see there for the scheme. Device arrays are transposed, `(b, i, j, k, v)`.
!< Separate kernels, never fused with the WENO face kernels (FNL register pressure, issue #47 R-4).

! FLUME modules
use :: adam_flume_dissipation_library, only : compute_dissipative_flux
use :: adam_flume_parameters,          only : IA_T, IA_U, IA_V, IA_W, IQ_RU
! third party modules
use :: penf,                           only : I4P, R8P

implicit none
private
public :: add_dissipative_fluxes_o2_dev
public :: add_dissipative_fluxes_o4_dev

contains
   ! public procedures
   subroutine add_dissipative_fluxes_o2_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, mu0, k0, tref, omega_mu, &
                                            omega_k, dxyz_gpu, is_null, q_aux_gpu, fl_gpu)
   !< Add the 2nd-order dissipative face fluxes of direction `d` to `fl_gpu`.
   integer(I4P), intent(in)    :: d                                     !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                            !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),    intent(in)    :: mu0, k0                               !< Viscosity and conductivity.
   real(R8P),    intent(in)    :: tref, omega_mu, omega_k               !< Temperature laws.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                       !< Blocks space steps [nb, 3].
   logical,      intent(in)    :: is_null(3)                            !< Null directions.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)       !< Face fluxes of direction `d`.
   integer(I4P)                :: iw(4)                                 !< Auxiliary indexes of W = (u, v, w, T).
   integer(I4P)                :: s1, s2, s3                            !< Steps along x, y, z: 1 active, 0 null.
   real(R8P)                   :: w1, w2, w3                            !< Weights along x, y, z: 1 active, 0 null.
   integer(I4P)                :: ot(3)                                 !< Private step along a direction.
   real(R8P)                   :: rdx(3)                                !< Private inverse steps, 0 along null.
   real(R8P)                   :: wl(4), wr(4)                          !< Private left and right cell states.
   real(R8P)                   :: g(4,3)                                !< Private face gradient.
   real(R8P)                   :: wf(4)                                 !< Private face state.
   real(R8P)                   :: h(4)                                  !< Private face flux.
   integer(I4P)                :: b, i, j, k, a, t                      !< Counters.

   iw = [IA_U, IA_V, IA_W, IA_T]
   s1 = merge(0_I4P, 1_I4P, is_null(1)) ; s2 = merge(0_I4P, 1_I4P, is_null(2)) ; s3 = merge(0_I4P, 1_I4P, is_null(3))
   w1 = real(s1, R8P) ; w2 = real(s2, R8P) ; w3 = real(s3, R8P)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_aux_gpu,fl_gpu)                    &
   !$acc& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,mu0,k0,tref,omega_mu,omega_k,iw,s1,s2,s3,w1,w2,w3)       &
   !$acc& private(ot,rdx,wl,wr,g,wf,h)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_aux_gpu,fl_gpu)                                                   &
   !$omp& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,mu0,k0,tref,omega_mu,omega_k,iw,s1,s2,s3,w1,w2,w3)       &
   !$omp& private(ot,rdx,wl,wr,g,wf,h)
   do k=1-dk, nk
   do j=1-dj, nj
   do i=1-di, ni
   do b=1, blocks_number
      rdx(1) = w1 / dxyz_gpu(b,1) ; rdx(2) = w2 / dxyz_gpu(b,2) ; rdx(3) = w3 / dxyz_gpu(b,3)
      !$acc loop seq
      do a=1, 4
         wl(a) = q_aux_gpu(b,i,j,k,iw(a))
         wr(a) = q_aux_gpu(b,i+di,j+dj,k+dk,iw(a))
      enddo
      !$acc loop seq
      do t=1, 3
         ot(1) = 0_I4P ; ot(2) = 0_I4P ; ot(3) = 0_I4P
         if (t == 1) ot(1) = s1
         if (t == 2) ot(2) = s2
         if (t == 3) ot(3) = s3
         !$acc loop seq
         do a=1, 4
            if (t == d) then
               g(a,t) = (wr(a) - wl(a)) * rdx(t)
            else
               g(a,t) = 0.25_R8P * rdx(t) * (q_aux_gpu(b,i+ot(1),j+ot(2),k+ot(3),iw(a))          - &
                                             q_aux_gpu(b,i-ot(1),j-ot(2),k-ot(3),iw(a))          + &
                                             q_aux_gpu(b,i+di+ot(1),j+dj+ot(2),k+dk+ot(3),iw(a)) - &
                                             q_aux_gpu(b,i+di-ot(1),j+dj-ot(2),k+dk-ot(3),iw(a)))
            endif
         enddo
      enddo
      !$acc loop seq
      do a=1, 4
         wf(a) = 0.5_R8P * (wl(a) + wr(a))
      enddo
      call compute_dissipative_flux(d=d, mu0=mu0, k0=k0, tref=tref, omega_mu=omega_mu, omega_k=omega_k, w=wf, g=g, f=h)
      !$acc loop seq
      do a=1, 4
         fl_gpu(b,i,j,k,IQ_RU+a-1) = fl_gpu(b,i,j,k,IQ_RU+a-1) + h(a)
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine add_dissipative_fluxes_o2_dev

   subroutine add_dissipative_fluxes_o4_dev(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, mu0, k0, tref, omega_mu, &
                                            omega_k, dxyz_gpu, is_null, q_aux_gpu, fl_gpu)
   !< Add the 4th-order conservative dissipative face fluxes of direction `d` to `fl_gpu`.
   integer(I4P), intent(in)    :: d                                     !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                            !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                       !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                         !< Actual blocks number.
   real(R8P),    intent(in)    :: mu0, k0                               !< Viscosity and conductivity.
   real(R8P),    intent(in)    :: tref, omega_mu, omega_k               !< Temperature laws.
   real(R8P),    intent(in)    :: dxyz_gpu(1:,1:)                       !< Blocks space steps [nb, 3].
   logical,      intent(in)    :: is_null(3)                            !< Null directions.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl_gpu(1:,1-di:,1-dj:,1-dk:,1:)       !< Face fluxes of direction `d`.
   integer(I4P)                :: iw(4)                                 !< Auxiliary indexes of W = (u, v, w, T).
   integer(I4P)                :: s1, s2, s3                            !< Steps along x, y, z: 1 active, 0 null.
   real(R8P)                   :: w1, w2, w3                            !< Weights along x, y, z: 1 active, 0 null.
   integer(I4P)                :: ot(3)                                 !< Private step along a direction.
   integer(I4P)                :: c(3)                                  !< Private cell index.
   real(R8P)                   :: rdx(3)                                !< Private inverse steps, 0 along null.
   real(R8P)                   :: wm(4,-1:2)                            !< Private states of the cells -1 ... 2.
   real(R8P)                   :: dm(4,-1:2)                            !< Private tangential derivatives there.
   real(R8P)                   :: g(4,3)                                !< Private face (or cell) gradient.
   real(R8P)                   :: gm(4,-1:2)                            !< Private cell fluxes of the correction.
   real(R8P)                   :: wf(4)                                 !< Private face (or cell) state.
   real(R8P)                   :: h(4)                                  !< Private face flux.
   real(R8P)                   :: hc(4)                                 !< Private cell flux.
   integer(I4P)                :: b, i, j, k, a, m, t                   !< Counters.

   iw = [IA_U, IA_V, IA_W, IA_T]
   s1 = merge(0_I4P, 1_I4P, is_null(1)) ; s2 = merge(0_I4P, 1_I4P, is_null(2)) ; s3 = merge(0_I4P, 1_I4P, is_null(3))
   w1 = real(s1, R8P) ; w2 = real(s2, R8P) ; w3 = real(s3, R8P)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(dxyz_gpu,q_aux_gpu,fl_gpu)                    &
   !$acc& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,mu0,k0,tref,omega_mu,omega_k,iw,s1,s2,s3,w1,w2,w3)       &
   !$acc& private(ot,c,rdx,wm,dm,g,gm,wf,h,hc)
   !$omp OMPLOOP collapse(4) DEVICEPTR(dxyz_gpu,q_aux_gpu,fl_gpu)                                                   &
   !$omp& firstprivate(d,di,dj,dk,ni,nj,nk,blocks_number,mu0,k0,tref,omega_mu,omega_k,iw,s1,s2,s3,w1,w2,w3)       &
   !$omp& private(ot,c,rdx,wm,dm,g,gm,wf,h,hc)
   do k=1-dk, nk
   do j=1-dj, nj
   do i=1-di, ni
   do b=1, blocks_number
      rdx(1) = w1 / dxyz_gpu(b,1) ; rdx(2) = w2 / dxyz_gpu(b,2) ; rdx(3) = w3 / dxyz_gpu(b,3)
      !$acc loop seq
      do m=-1, 2
         !$acc loop seq
         do a=1, 4
            wm(a,m) = q_aux_gpu(b,i+m*di,j+m*dj,k+m*dk,iw(a))
         enddo
      enddo
      ! face flux from 4th-order face values and gradients
      !$acc loop seq
      do t=1, 3
         ot(1) = 0_I4P ; ot(2) = 0_I4P ; ot(3) = 0_I4P
         if (t == 1) ot(1) = s1
         if (t == 2) ot(2) = s2
         if (t == 3) ot(3) = s3
         if (t == d) then
            !$acc loop seq
            do a=1, 4
               g(a,t) = (27._R8P * (wm(a,1) - wm(a,0)) - (wm(a,2) - wm(a,-1))) * rdx(t) / 24._R8P
            enddo
         else
            !$acc loop seq
            do m=-1, 2
               c(1) = i + m * di ; c(2) = j + m * dj ; c(3) = k + m * dk
               !$acc loop seq
               do a=1, 4
                  dm(a,m) = (8._R8P * (q_aux_gpu(b,c(1)+ot(1),c(2)+ot(2),c(3)+ot(3),iw(a))        -  &
                                       q_aux_gpu(b,c(1)-ot(1),c(2)-ot(2),c(3)-ot(3),iw(a)))       -  &
                             (q_aux_gpu(b,c(1)+2*ot(1),c(2)+2*ot(2),c(3)+2*ot(3),iw(a))           -  &
                              q_aux_gpu(b,c(1)-2*ot(1),c(2)-2*ot(2),c(3)-2*ot(3),iw(a)))) * rdx(t) / 12._R8P
               enddo
            enddo
            !$acc loop seq
            do a=1, 4
               g(a,t) = (9._R8P * (dm(a,0) + dm(a,1)) - (dm(a,-1) + dm(a,2))) / 16._R8P
            enddo
         endif
      enddo
      !$acc loop seq
      do a=1, 4
         wf(a) = (9._R8P * (wm(a,0) + wm(a,1)) - (wm(a,-1) + wm(a,2))) / 16._R8P
      enddo
      ! as on the CPU: a non-positive face temperature (undershoot across a jump) takes the mean of its two cells
      if (.not.(wf(4) > 0._R8P)) wf(4) = 0.5_R8P * (wm(4,0) + wm(4,1))
      call compute_dissipative_flux(d=d, mu0=mu0, k0=k0, tref=tref, omega_mu=omega_mu, omega_k=omega_k, w=wf, g=g, f=h)
      ! Shu-Osher correction from the cell fluxes, 2nd-order central gradients
      !$acc loop seq
      do m=-1, 2
         c(1) = i + m * di ; c(2) = j + m * dj ; c(3) = k + m * dk
         !$acc loop seq
         do t=1, 3
            ot(1) = 0_I4P ; ot(2) = 0_I4P ; ot(3) = 0_I4P
            if (t == 1) ot(1) = s1
            if (t == 2) ot(2) = s2
            if (t == 3) ot(3) = s3
            !$acc loop seq
            do a=1, 4
               g(a,t) = 0.5_R8P * rdx(t) * (q_aux_gpu(b,c(1)+ot(1),c(2)+ot(2),c(3)+ot(3),iw(a)) - &
                                            q_aux_gpu(b,c(1)-ot(1),c(2)-ot(2),c(3)-ot(3),iw(a)))
            enddo
         enddo
         !$acc loop seq
         do a=1, 4
            wf(a) = wm(a,m)
         enddo
         call compute_dissipative_flux(d=d, mu0=mu0, k0=k0, tref=tref, omega_mu=omega_mu, omega_k=omega_k, w=wf, g=g, &
                                       f=hc)
         !$acc loop seq
         do a=1, 4
            gm(a,m) = hc(a)
         enddo
      enddo
      !$acc loop seq
      do a=1, 4
         h(a) = h(a) - (gm(a,2) - gm(a,1) - gm(a,0) + gm(a,-1)) / 48._R8P
         fl_gpu(b,i,j,k,IQ_RU+a-1) = fl_gpu(b,i,j,k,IQ_RU+a-1) + h(a)
      enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine add_dissipative_fluxes_o4_dev
endmodule adam_flume_fnl_dissipation_kernels
