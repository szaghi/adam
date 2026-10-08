!< ADAM, FLUME CPU kernels of the dissipative face fluxes (issue #65, M4).

module adam_flume_cpu_dissipation_kernels
!< ADAM, FLUME CPU kernels of the dissipative face fluxes (issue #65, M4).
!<
!< The viscous and heat-conduction face fluxes are ADDED to the inviscid face fluxes of direction `d` (`fl`), before the
!< seam accumulation and the flux difference: the AMR reflux and the inter-realm register then see the total flux.
!< They read `W = (u, v, w, T)` from the auxiliary variables, whose indexes are those of every model, so one kernel
!< serves Euler and MHD; the host selects the kernel of the order (2 or 4, `[numerics] dissipative_order`), never a
!< branch inside the loops.
!<
!< Face `(i,j,k)` lies between cells `(i,j,k)` and `(i,j,k) + e_d`; `m` counts the cells along `d` from the left one
!< (`m = 0`). A null direction has a zero step, so its derivatives are differences of a cell with itself, zero exactly,
!< and its ghost cells are never read.
!<
!< Order 2: `h = G(W_f, grad W_f)`, `W_f = (W_0 + W_1) / 2`, normal derivative `(W_1 - W_0) / dx`, tangential ones the
!< mean of the two cells' central differences.
!<
!< Order 4 (Shu & Osher 1989, J. Comput. Phys. 83, 32): the conservative flux of a point-value scheme is
!< `h = G_f - dx^2 / 24 G''_f + O(dx^4)`, with
!<```
!< W_f       = (9 (W_0 + W_1) - (W_-1 + W_2)) / 16
!< dW/dn_f   = (27 (W_1 - W_0) - (W_2 - W_-1)) / (24 dx)
!< dW/dt_f   = (9 (D_0 + D_1) - (D_-1 + D_2)) / 16,  D_m = (8 (W_+1 - W_-1) - (W_+2 - W_-2)) / (12 dt)  at cell m
!< G''_f dx^2 = (G_2 - G_1 - G_0 + G_-1) / 2,        G_m from 2nd-order central gradients at cell m
!<```
!< so `h = G_f - (G_2 - G_1 - G_0 + G_-1) / 48`. Across a strong temperature jump the 4-point face temperature can
!< undershoot below zero (a power law is then NaN): where it is not positive, the face takes the mean of its two cells.
!< Cells `-2 ... 3` along `d` and `-2 ... 2` along the tangents are read: `ngc >= 3`. The correction is built from 2nd-order cell fluxes for every term, linear or not, so for constant
!< `mu` the scheme is a 4th-order `u_xx` on 6 points, not the compact 5-point one; the order is verified by measurement.
!<
!< The Ohmic fluxes (MHD, issue #65 P3) use the same stencils on `W = B` (auxiliary `IA_BX..IA_BZ`), adding to the field
!< and energy fluxes: the kernels `add_resistive_fluxes_o2/o4`, called by the host only when the resistivity is on.

! FLUME modules
use :: adam_flume_dissipation_library, only : compute_dissipative_flux, compute_resistive_flux
use :: adam_flume_parameters,          only : IA_BX, IA_BY, IA_BZ, IA_T, IA_U, IA_V, IA_W, IQ_BX, IQ_RE, IQ_RU
! third party modules
use :: penf,                           only : I4P, R8P

implicit none
private
public :: add_dissipative_fluxes_o2
public :: add_dissipative_fluxes_o4
public :: add_resistive_fluxes_o2
public :: add_resistive_fluxes_o4

integer(I4P), parameter :: IW(4)=[IA_U, IA_V, IA_W, IA_T]  !< Auxiliary indexes of W = (u, v, w, T).
integer(I4P), parameter :: IB(3)=[IA_BX, IA_BY, IA_BZ]     !< Auxiliary indexes of the magnetic field.

contains
   ! public procedures
   subroutine add_dissipative_fluxes_o2(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, mu0, k0, tref, omega_mu, omega_k, &
                                        dxyz, is_null, q_aux, fl)
   !< Add the 2nd-order dissipative face fluxes of direction `d` to `fl`.
   integer(I4P), intent(in)    :: d                                 !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                        !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: mu0, k0                           !< Viscosity and conductivity.
   real(R8P),    intent(in)    :: tref, omega_mu, omega_k           !< Temperature laws.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                       !< Blocks space steps [3, nb].
   logical,      intent(in)    :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)       !< Face fluxes of direction `d`.
   integer(I4P)                :: o(3,3)                            !< Steps o(:,t) along each direction.
   real(R8P)                   :: w0(4), w1(4)                      !< Left and right cell states.
   real(R8P)                   :: g(4,3)                            !< Face gradient.
   real(R8P)                   :: h(4)                              !< Face flux.
   real(R8P)                   :: rdx(3)                            !< Inverse steps, 0 along null directions.
   integer(I4P)                :: c(3)                              !< Right cell index.
   integer(I4P)                :: b, i, j, k, t                     !< Counters.

   o = steps(is_null)
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q_aux, fl)
   do b=1, blocks_number
      do k=1-dk, nk
         do j=1-dj, nj
            do i=1-di, ni
               rdx = inverse_steps(dxyz(:,b), is_null)
               c = [i+di, j+dj, k+dk]
               w0 = q_aux(IW,i,j,k,b)
               w1 = q_aux(IW,c(1),c(2),c(3),b)
               do t=1, 3
                  if (t == d) then
                     g(:,t) = (w1 - w0) * rdx(t)
                  else
                     g(:,t) = 0.25_R8P * rdx(t) * (q_aux(IW,i+o(1,t),j+o(2,t),k+o(3,t),b)          - &
                                                   q_aux(IW,i-o(1,t),j-o(2,t),k-o(3,t),b)          + &
                                                   q_aux(IW,c(1)+o(1,t),c(2)+o(2,t),c(3)+o(3,t),b) - &
                                                   q_aux(IW,c(1)-o(1,t),c(2)-o(2,t),c(3)-o(3,t),b))
                  endif
               enddo
               call compute_dissipative_flux(d=d, mu0=mu0, k0=k0, tref=tref, omega_mu=omega_mu, omega_k=omega_k, &
                                             w=0.5_R8P * (w0 + w1), g=g, f=h)
               fl(IQ_RU:IQ_RU+3,i,j,k,b) = fl(IQ_RU:IQ_RU+3,i,j,k,b) + h
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine add_dissipative_fluxes_o2

   subroutine add_dissipative_fluxes_o4(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, mu0, k0, tref, omega_mu, omega_k, &
                                        dxyz, is_null, q_aux, fl)
   !< Add the 4th-order conservative dissipative face fluxes of direction `d` to `fl`.
   integer(I4P), intent(in)    :: d                                 !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                        !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: mu0, k0                           !< Viscosity and conductivity.
   real(R8P),    intent(in)    :: tref, omega_mu, omega_k           !< Temperature laws.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                       !< Blocks space steps [3, nb].
   logical,      intent(in)    :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)       !< Face fluxes of direction `d`.
   integer(I4P)                :: o(3,3)                            !< Steps o(:,t) along each direction.
   real(R8P)                   :: wm(4,-1:2)                        !< States of the cells -1 ... 2 along `d`.
   real(R8P)                   :: dm(4,-1:2)                        !< 4th-order tangential derivative at those cells.
   real(R8P)                   :: g(4,3)                            !< Face (or cell) gradient.
   real(R8P)                   :: gm(4,-1:2)                        !< Cell fluxes of the correction.
   real(R8P)                   :: wf(4)                             !< Face state.
   real(R8P)                   :: h(4)                              !< Face flux.
   real(R8P)                   :: rdx(3)                            !< Inverse steps, 0 along null directions.
   integer(I4P)                :: c(3)                              !< Cell index.
   integer(I4P)                :: b, i, j, k, m, t                  !< Counters.

   o = steps(is_null)
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q_aux, fl)
   do b=1, blocks_number
      do k=1-dk, nk
         do j=1-dj, nj
            do i=1-di, ni
               rdx = inverse_steps(dxyz(:,b), is_null)
               do m=-1, 2
                  wm(:,m) = q_aux(IW,i+m*di,j+m*dj,k+m*dk,b)
               enddo
               ! face flux from 4th-order face values and gradients
               do t=1, 3
                  if (t == d) then
                     g(:,t) = (27._R8P * (wm(:,1) - wm(:,0)) - (wm(:,2) - wm(:,-1))) * rdx(t) / 24._R8P
                  else
                     do m=-1, 2
                        c = [i+m*di, j+m*dj, k+m*dk]
                        dm(:,m) = (8._R8P * (q_aux(IW,c(1)+o(1,t),c(2)+o(2,t),c(3)+o(3,t),b)        -  &
                                             q_aux(IW,c(1)-o(1,t),c(2)-o(2,t),c(3)-o(3,t),b))       -  &
                                   (q_aux(IW,c(1)+2*o(1,t),c(2)+2*o(2,t),c(3)+2*o(3,t),b)           -  &
                                    q_aux(IW,c(1)-2*o(1,t),c(2)-2*o(2,t),c(3)-2*o(3,t),b))) * rdx(t) / 12._R8P
                     enddo
                     g(:,t) = (9._R8P * (dm(:,0) + dm(:,1)) - (dm(:,-1) + dm(:,2))) / 16._R8P
                  endif
               enddo
               wf = (9._R8P * (wm(:,0) + wm(:,1)) - (wm(:,-1) + wm(:,2))) / 16._R8P
               ! across a strong temperature jump the 4-point face value undershoots, possibly below zero, where a
               ! power law (T / tref)^omega is NaN: there the face temperature is the mean of its two cells
               if (.not.(wf(4) > 0._R8P)) wf(4) = 0.5_R8P * (wm(4,0) + wm(4,1))
               call compute_dissipative_flux(d=d, mu0=mu0, k0=k0, tref=tref, omega_mu=omega_mu, omega_k=omega_k, &
                                             w=wf, g=g, f=h)
               ! Shu-Osher correction from the cell fluxes, 2nd-order central gradients
               do m=-1, 2
                  c = [i+m*di, j+m*dj, k+m*dk]
                  do t=1, 3
                     g(:,t) = 0.5_R8P * rdx(t) * (q_aux(IW,c(1)+o(1,t),c(2)+o(2,t),c(3)+o(3,t),b) - &
                                                  q_aux(IW,c(1)-o(1,t),c(2)-o(2,t),c(3)-o(3,t),b))
                  enddo
                  call compute_dissipative_flux(d=d, mu0=mu0, k0=k0, tref=tref, omega_mu=omega_mu, omega_k=omega_k, &
                                                w=wm(:,m), g=g, f=gm(:,m))
               enddo
               h = h - (gm(:,2) - gm(:,1) - gm(:,0) + gm(:,-1)) / 48._R8P
               fl(IQ_RU:IQ_RU+3,i,j,k,b) = fl(IQ_RU:IQ_RU+3,i,j,k,b) + h
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine add_dissipative_fluxes_o4

   subroutine add_resistive_fluxes_o2(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, eta, dxyz, is_null, q_aux, fl)
   !< Add the 2nd-order Ohmic face fluxes of direction `d` to `fl` (MHD).
   integer(I4P), intent(in)    :: d                                 !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                        !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: eta                               !< Magnetic diffusivity.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                       !< Blocks space steps [3, nb].
   logical,      intent(in)    :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)       !< Face fluxes of direction `d`.
   integer(I4P)                :: o(3,3)                            !< Steps o(:,t) along each direction.
   real(R8P)                   :: b0(3), b1(3)                      !< Left and right cell fields.
   real(R8P)                   :: g(3,3)                            !< Face gradient.
   real(R8P)                   :: h(4)                              !< Face flux.
   real(R8P)                   :: rdx(3)                            !< Inverse steps, 0 along null directions.
   integer(I4P)                :: c(3)                              !< Right cell index.
   integer(I4P)                :: b, i, j, k, t                     !< Counters.

   o = steps(is_null)
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q_aux, fl)
   do b=1, blocks_number
      do k=1-dk, nk
         do j=1-dj, nj
            do i=1-di, ni
               rdx = inverse_steps(dxyz(:,b), is_null)
               c = [i+di, j+dj, k+dk]
               b0 = q_aux(IB,i,j,k,b)
               b1 = q_aux(IB,c(1),c(2),c(3),b)
               do t=1, 3
                  if (t == d) then
                     g(:,t) = (b1 - b0) * rdx(t)
                  else
                     g(:,t) = 0.25_R8P * rdx(t) * (q_aux(IB,i+o(1,t),j+o(2,t),k+o(3,t),b)          - &
                                                   q_aux(IB,i-o(1,t),j-o(2,t),k-o(3,t),b)          + &
                                                   q_aux(IB,c(1)+o(1,t),c(2)+o(2,t),c(3)+o(3,t),b) - &
                                                   q_aux(IB,c(1)-o(1,t),c(2)-o(2,t),c(3)-o(3,t),b))
                  endif
               enddo
               call compute_resistive_flux(d=d, eta=eta, b=0.5_R8P * (b0 + b1), g=g, f=h)
               fl(IQ_BX:IQ_BX+2,i,j,k,b) = fl(IQ_BX:IQ_BX+2,i,j,k,b) + h(1:3)
               fl(IQ_RE,i,j,k,b) = fl(IQ_RE,i,j,k,b) + h(4)
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine add_resistive_fluxes_o2

   subroutine add_resistive_fluxes_o4(d, di, dj, dk, ni, nj, nk, ngc, blocks_number, eta, dxyz, is_null, q_aux, fl)
   !< Add the 4th-order conservative Ohmic face fluxes of direction `d` to `fl` (MHD): the stencils of
   !< `add_dissipative_fluxes_o4` on `W = B`, with the Shu-Osher correction.
   integer(I4P), intent(in)    :: d                                 !< Direction, 1=x, 2=y, 3=z.
   integer(I4P), intent(in)    :: di, dj, dk                        !< Unit step along `d`.
   integer(I4P), intent(in)    :: ni, nj, nk, ngc                   !< Grid dimensions.
   integer(I4P), intent(in)    :: blocks_number                     !< Actual blocks number.
   real(R8P),    intent(in)    :: eta                               !< Magnetic diffusivity.
   real(R8P),    intent(in)    :: dxyz(1:,1:)                       !< Blocks space steps [3, nb].
   logical,      intent(in)    :: is_null(3)                        !< Null directions.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)       !< Face fluxes of direction `d`.
   integer(I4P)                :: o(3,3)                            !< Steps o(:,t) along each direction.
   real(R8P)                   :: bm(3,-1:2)                        !< Fields of the cells -1 ... 2 along `d`.
   real(R8P)                   :: dm(3,-1:2)                        !< 4th-order tangential derivative at those cells.
   real(R8P)                   :: g(3,3)                            !< Face (or cell) gradient.
   real(R8P)                   :: gm(4,-1:2)                        !< Cell fluxes of the correction.
   real(R8P)                   :: h(4)                              !< Face flux.
   real(R8P)                   :: rdx(3)                            !< Inverse steps, 0 along null directions.
   integer(I4P)                :: c(3)                              !< Cell index.
   integer(I4P)                :: b, i, j, k, m, t                  !< Counters.

   o = steps(is_null)
   !$omp parallel do collapse(4) default(firstprivate) shared(dxyz, q_aux, fl)
   do b=1, blocks_number
      do k=1-dk, nk
         do j=1-dj, nj
            do i=1-di, ni
               rdx = inverse_steps(dxyz(:,b), is_null)
               do m=-1, 2
                  bm(:,m) = q_aux(IB,i+m*di,j+m*dj,k+m*dk,b)
               enddo
               ! face flux from 4th-order face values and gradients
               do t=1, 3
                  if (t == d) then
                     g(:,t) = (27._R8P * (bm(:,1) - bm(:,0)) - (bm(:,2) - bm(:,-1))) * rdx(t) / 24._R8P
                  else
                     do m=-1, 2
                        c = [i+m*di, j+m*dj, k+m*dk]
                        dm(:,m) = (8._R8P * (q_aux(IB,c(1)+o(1,t),c(2)+o(2,t),c(3)+o(3,t),b)        -  &
                                             q_aux(IB,c(1)-o(1,t),c(2)-o(2,t),c(3)-o(3,t),b))       -  &
                                   (q_aux(IB,c(1)+2*o(1,t),c(2)+2*o(2,t),c(3)+2*o(3,t),b)           -  &
                                    q_aux(IB,c(1)-2*o(1,t),c(2)-2*o(2,t),c(3)-2*o(3,t),b))) * rdx(t) / 12._R8P
                     enddo
                     g(:,t) = (9._R8P * (dm(:,0) + dm(:,1)) - (dm(:,-1) + dm(:,2))) / 16._R8P
                  endif
               enddo
               call compute_resistive_flux(d=d, eta=eta,                                                       &
                                           b=(9._R8P * (bm(:,0) + bm(:,1)) - (bm(:,-1) + bm(:,2))) / 16._R8P, &
                                           g=g, f=h)
               ! Shu-Osher correction from the cell fluxes, 2nd-order central gradients
               do m=-1, 2
                  c = [i+m*di, j+m*dj, k+m*dk]
                  do t=1, 3
                     g(:,t) = 0.5_R8P * rdx(t) * (q_aux(IB,c(1)+o(1,t),c(2)+o(2,t),c(3)+o(3,t),b) - &
                                                  q_aux(IB,c(1)-o(1,t),c(2)-o(2,t),c(3)-o(3,t),b))
                  enddo
                  call compute_resistive_flux(d=d, eta=eta, b=bm(:,m), g=g, f=gm(:,m))
               enddo
               h = h - (gm(:,2) - gm(:,1) - gm(:,0) + gm(:,-1)) / 48._R8P
               fl(IQ_BX:IQ_BX+2,i,j,k,b) = fl(IQ_BX:IQ_BX+2,i,j,k,b) + h(1:3)
               fl(IQ_RE,i,j,k,b) = fl(IQ_RE,i,j,k,b) + h(4)
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine add_resistive_fluxes_o4

   ! private procedures
   pure function steps(is_null) result(o)
   !< Return the unit steps along each direction, `o(:,t) = e_t`, zero along the null directions.
   logical, intent(in) :: is_null(3) !< Null directions.
   integer(I4P)        :: o(3,3)     !< Steps.
   integer(I4P)        :: t          !< Counter.

   o = 0_I4P
   do t=1, 3
      if (.not.is_null(t)) o(t,t) = 1_I4P
   enddo
   endfunction steps

   pure function inverse_steps(ds, is_null) result(rdx)
   !< Return the inverse space steps, zero along the null directions.
   real(R8P), intent(in) :: ds(3)      !< Space steps.
   logical,   intent(in) :: is_null(3) !< Null directions.
   real(R8P)             :: rdx(3)     !< Inverse steps.

   rdx = merge(0._R8P, 1._R8P / ds, is_null)
   endfunction inverse_steps
endmodule adam_flume_cpu_dissipation_kernels
