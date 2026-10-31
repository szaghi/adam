!< Measure the central flux of the hybrid scheme (issue #119, T2): order, kinetic energy identity, pressure equilibrium.
program test_kep
!< Measure the central flux of the hybrid scheme (issue #119, T2): order, kinetic energy identity, pressure equilibrium.
!<
!< On periodic lines along x, y and z, with `compute_kep_flux`: the flux difference of a smooth field against its exact derivative, at n = 32, 64,
!< 128, order >= 5.5 on every component; `sum(u . d(rho u) - |u|^2 / 2 d rho)` on a random field of uniform pressure, the discrete kinetic energy
!< identity, within 1e-12 of its scale; dp/dt of a uniform velocity and pressure across a 100:1 step, within 1e-12. Stops with an error otherwise.

use :: adam_flume_hybrid_library, only : compute_kep_flux
use :: penf,                      only : I4P, R8P

implicit none
real(R8P), parameter :: GAMMA = 1.4_R8P
real(R8P), parameter :: TWO_PI = 8._R8P * atan(1._R8P)
real(R8P), parameter :: A(3) = [0.3_R8P, 0.2_R8P, -0.1_R8P] !< Transverse and normal velocity amplitudes.
integer(I4P), parameter :: M = 64_I4P !< Cells of the random line.
real(R8P)            :: err(5,3,3), order
integer(I4P)         :: d, l, n, v
logical              :: ok

ok = .true.
do d=1, 3
   do l=1, 3
      call smooth_error(d, 32_I4P * 2_I4P**(l-1_I4P), err(:,d,l))
   enddo
   do v=1, 5
      order = log(err(v,d,2) / err(v,d,3)) / log(2._R8P)
      print '(A,I1,A,I1,A,3ES11.3,A,F6.3)', 'd ', d, ' v ', v, ' errors', err(v,d,:), ' order', order
      ok = ok .and. order >= 5.5_R8P
   enddo
enddo
call identities(ok)
if (.not.ok) error stop 'test_kep: FAIL'
print '(A)', 'test_kep: PASS'

contains
   subroutine smooth_error(d, n, e)
   !< Max error of the flux difference against the derivative of the exact flux, the latter by a 4th-order difference with a step 100 times smaller.
   integer(I4P), intent(in)  :: d, n
   real(R8P),    intent(out) :: e(5)
   real(R8P), allocatable    :: st(:,:), fl(:,:), ex(:,:)
   real(R8P)                 :: x, dx, w(5,-2:3), f(5), h, up(5)
   integer(I4P)              :: i, m

   dx = TWO_PI / n
   allocate(st(5,-2:n+3), fl(5,0:n), ex(5,n))
   do i=-2, n+3
      st(:,i) = prim(real(i, R8P) * dx)
   enddo
   do i=0, n
      do m=-2, 3
         w(:,m) = st(:,i+m)
      enddo
      call compute_kep_flux(d, GAMMA, w, f)
      fl(:,i) = f
   enddo
   h = 1.e-2_R8P * dx
   do i=1, n
      x = real(i, R8P) * dx
      up = 8._R8P * (physical(prim(x + h), d) - physical(prim(x - h), d)) - physical(prim(x + 2._R8P * h), d) + physical(prim(x - 2._R8P * h), d)
      ex(:,i) = up / (12._R8P * h)
   enddo
   e = 0._R8P
   do i=1, n
      e = max(e, abs((fl(:,i) - fl(:,i-1)) / dx - ex(:,i)))
   enddo
   endsubroutine smooth_error

   pure function prim(x) result(w)
   !< A smooth state on the line.
   real(R8P), intent(in) :: x
   real(R8P)             :: w(5)

   w(1) = 1._R8P + 0.2_R8P * sin(x)
   w(2) = A(1) * cos(x) + 0.1_R8P ; w(3) = A(2) * sin(2._R8P * x) ; w(4) = A(3) * cos(x)
   w(5) = 1._R8P + 0.1_R8P * cos(x)
   endfunction prim

   pure function physical(w, d) result(f)
   !< The Euler flux of direction d.
   real(R8P),    intent(in) :: w(5)
   integer(I4P), intent(in) :: d
   real(R8P)                :: f(5), e

   e = w(5) / (GAMMA - 1._R8P) + 0.5_R8P * w(1) * sum(w(2:4)**2)
   f(1) = w(1) * w(1+d)
   f(2:4) = w(1) * w(1+d) * w(2:4)
   f(1+d) = f(1+d) + w(5)
   f(5) = w(1+d) * (e + w(5))
   endfunction physical

   subroutine identities(ok)
   !< The kinetic energy identity and the pressure equilibrium, along x.
   logical, intent(inout)   :: ok
   real(R8P)                :: st(5,1-3:M+3), dq(5,M), r(5,1-3:M+3), ke, scale, dp
   integer(I4P)             :: i

   call random_seed() ; call random_number(r)
   st = r ; st(1,:) = 0.5_R8P + st(1,:) ; st(2:4,:) = st(2:4,:) - 0.5_R8P ; st(5,:) = 1._R8P
   call residual(st, dq)
   ke = 0._R8P ; scale = 0._R8P
   do i=1, M
      ke = ke + sum(st(2:4,i) * dq(2:4,i)) - 0.5_R8P * sum(st(2:4,i)**2) * dq(1,i)
      scale = scale + abs(st(2,i) * dq(2,i))
   enddo
   print '(A,ES10.2,A)', 'kinetic energy identity, relative ', abs(ke) / scale, ' (bound 1e-12)'
   ok = ok .and. abs(ke) / scale < 1.e-12_R8P
   st = 0._R8P ; st(2,:) = 1._R8P ; st(5,:) = 1._R8P
   do i=1-3, M+3
      st(1,i) = merge(100._R8P, 1._R8P, i > M/4 .and. i <= 3*M/4)
   enddo
   call residual(st, dq)
   dp = 0._R8P
   do i=1, M
      dp = max(dp, abs((GAMMA - 1._R8P) * (dq(5,i) - st(2,i) * dq(2,i) + 0.5_R8P * st(2,i)**2 * dq(1,i))))
   enddo
   print '(A,ES10.2,A)', 'pressure equilibrium, largest dp/dt across 100:1 ', dp, ' (bound 1e-12)'
   ok = ok .and. dp < 1.e-12_R8P

   endsubroutine identities

   subroutine residual(st, dq)
   !< Flux difference on the periodic line along x, unit step.
   real(R8P), intent(inout) :: st(5,1-3:M+3)
   real(R8P), intent(out)   :: dq(5,M)
   real(R8P)                :: fl(5,0:M), w(5,-2:3), f(5)
   integer(I4P)             :: i, c

   st(:,1-3:0) = st(:,M-2:M) ; st(:,M+1:M+3) = st(:,1:3)
   do i=0, M
      do c=-2, 3
         w(:,c) = st(:,i+c)
      enddo
      call compute_kep_flux(1_I4P, GAMMA, w, f)
      fl(:,i) = f
   enddo
   do i=1, M
      dq(:,i) = -(fl(:,i) - fl(:,i-1))
   enddo
   endsubroutine residual
endprogram test_kep
