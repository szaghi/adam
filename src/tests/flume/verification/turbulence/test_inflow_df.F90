!< The digital-filter inflow generator of FLUME (issue #119, T3.7), on its plane alone: leg df1 of check.sh.

program test_inflow_df
!< On a synthetic profile (U = eta^(1/7), stresses decaying to 0 at the layer edge), the generator must give, over a long
!< sample:
!< 1. the prescribed Reynolds stresses: the sums over the layer of <u'u'>, <v'v'>, <w'w'>, <u'v'> within 3 % of the target
!<    sums, every node within 10 %, the mean fluctuations below 0.03 of their rms;
!< 2. the temporal correlation exp(-pi t / (2 T)) of the Xie-Castro recursion at the lags T / 2 and T, within 0.03;
!< 3. the spanwise integral length of the filtered signal (one-sided integral of its correlation) within 10 % of the one of
!<    exp(-pi r / (2 L)), 2 L / pi (the shape is (1 + pi r / n) exp(-pi r / n), the integral the same);
!< 4. a periodic span: the correlation across the seam (nodes nz and 1) equal to the interior one at lag 1 within 0.02;
!< 5. a bitwise restart: a copy reloaded from the saved state continues identically.

use :: adam_flume_inflow_df_object, only : flume_inflow_df_object, profile_mean
use :: adam_mpih_global,            only : mpih
use :: penf,                        only : I4P, I8P, R8P

implicit none
type(flume_inflow_df_object) :: df, dg                              !< Generators.
real(R8P), parameter         :: PI=acos(-1._R8P)
real(R8P), allocatable       :: s1(:,:), s2(:,:), s3(:,:), s4(:,:)  !< Sums of u'u', v'v', w'w', u'v' per node [ny, 1].
real(R8P), allocatable       :: m1(:), m2(:), m3(:)                 !< Sums of u', v', w' per node.
real(R8P), allocatable       :: hist(:,:)                           !< psi1 history at a few nodes [lags, nodes].
real(R8P)                    :: target(4), got(4), err, rho_t(2), lag_t(2), cz(0:40), cs, ci, lz, dt, f(4)
integer(I4P), parameter      :: NSTEP=40000_I4P, NWARM=200_I4P
integer(I4P)                 :: it, j, k, n, m, mlag(2)
logical                      :: ok

call mpih%initialize(do_mpi_init=.true., verbose=.false.)
ok = .true.
call make(df)
call df%setup(emin=[0._R8P, 0._R8P, 0._R8P], emax=[1._R8P, 2._R8P, 2._R8P], periodic=[.false., .false., .true.])
dt = 0.05_R8P * df%t_scale
mlag = [10_I4P, 20_I4P] ! T / 2 and T
allocate(s1(df%ny,1), s2(df%ny,1), s3(df%ny,1), s4(df%ny,1), m1(df%ny), m2(df%ny), m3(df%ny), source=0._R8P)
allocate(hist(0:NSTEP, 4), source=0._R8P)
cz = 0._R8P ; cs = 0._R8P ; ci = 0._R8P
n = 0_I4P
do it=1, NWARM
   call df%advance(dt=dt)
enddo
do it=1, NSTEP
   call df%advance(dt=dt)
   do k=1, df%nz
      do j=1, df%ny
         f = df%fluctuation(j, k)
         s1(j,1) = s1(j,1) + f(1)**2 ; s2(j,1) = s2(j,1) + f(2)**2 ; s3(j,1) = s3(j,1) + f(3)**2 ; s4(j,1) = s4(j,1) + f(1) * f(2)
         m1(j) = m1(j) + f(1) ; m2(j) = m2(j) + f(2) ; m3(j) = m3(j) + f(3)
      enddo
      ! spanwise correlation of psi3 on the node row j = 5, periodic lags
      do m=0, 40
         cz(m) = cz(m) + df%psi(3,5,k) * df%psi(3,5,modulo(k+m-1, df%nz)+1)
      enddo
   enddo
   cs = cs + sum(df%psi(3,:,df%nz) * df%psi(3,:,1))                ! across the seam, every node row
   ci = ci + sum(df%psi(3,:,1:df%nz-1) * df%psi(3,:,2:df%nz)) / (df%nz - 1) ! interior, lag 1, every node row
   hist(it,:) = [df%psi(1,3,3), df%psi(1,10,17), df%psi(1,20,31), df%psi(1,30,9)]
   n = n + df%nz
enddo
! 1. stresses
do j=1, 4
   target(j) = sum(stress(j))
enddo
got = [sum(s1) / n, sum(s2) / n, sum(s3) / n, sum(s4) / n]
do j=1, 4
   err = abs(got(j) - target(j)) / abs(target(j))
   print '(A,I1,A,ES12.4,A,ES12.4,A,F7.4,A)', 'stress ', j, ': sum over the layer ', got(j), ' target ', target(j), &
                                           ' relative ', err, ' (bound 0.03)'
   ok = ok .and. err < 0.03_R8P
enddo
err = 0._R8P
do j=1, df%ny
   if (df%node_mean(3,j) > 0.2_R8P * maxval(df%node_mean(3,:))) then
      err = max(err, abs(s1(j,1) / n - df%node_mean(3,j)**2) / df%node_mean(3,j)**2)
      err = max(err, abs(s3(j,1) / n - df%node_mean(6,j)**2) / df%node_mean(6,j)**2)
   endif
enddo
print '(A,F7.4,A)', 'stresses u''u'', w''w'' node by node (where u_rms > 0.2 max): worst relative ', err, ' (bound 0.10)'
ok = ok .and. err < 0.10_R8P
err = 0._R8P
do j=1, df%ny
   if (df%node_mean(3,j) > 0._R8P) err = max(err, abs(m1(j) / n) / df%node_mean(3,j), abs(m3(j) / n) / df%node_mean(6,j))
enddo
print '(A,F7.4,A)', 'mean fluctuations over their rms: worst ', err, ' (bound 0.03)'
ok = ok .and. err < 0.03_R8P
! 2. temporal correlation
do m=1, 2
   rho_t(m) = 0._R8P
   do j=1, 4
      rho_t(m) = rho_t(m) + sum(hist(1:NSTEP-mlag(m),j) * hist(1+mlag(m):NSTEP,j)) / sum(hist(1:NSTEP,j)**2)
   enddo
   rho_t(m) = rho_t(m) / 4._R8P
   lag_t(m) = real(mlag(m), R8P) * dt
   err = abs(rho_t(m) - exp(-PI * lag_t(m) / (2._R8P * df%t_scale)))
   print '(A,F5.2,A,F7.4,A,F7.4,A,F7.4,A)', 'temporal correlation at t/T = ', lag_t(m) / df%t_scale, ': ', rho_t(m), &
                                            ' expected ', exp(-PI * lag_t(m) / (2._R8P * df%t_scale)), ' diff ', err, ' (bound 0.03)'
   ok = ok .and. err < 0.03_R8P
enddo
! 3. spanwise integral length
cz = cz / cz(0)
lz = (0.5_R8P * cz(0) + sum(cz(1:df%nz/2))) * df%dz ! half the periodic span: the lags beyond it repeat the ones below
err = abs(lz - 2._R8P * df%length(3) * df%delta / PI) / (2._R8P * df%length(3) * df%delta / PI)
print '(A,F7.4,A,F7.4,A,F7.4,A)', 'spanwise integral length ', lz, ' expected 2 L / pi = ', 2._R8P * df%length(3) * df%delta / PI, &
                                  ' relative ', err, ' (bound 0.10)'
ok = ok .and. err < 0.10_R8P
! 4. periodic seam
cs = cs / (NSTEP * df%ny) ; ci = ci / (NSTEP * df%ny)
print '(A,F7.4,A,F7.4,A,F7.4,A)', 'lag-1 correlation across the seam ', cs, ' interior ', ci, ' diff ', abs(cs - ci), ' (bound 0.02)'
ok = ok .and. abs(cs - ci) < 0.02_R8P
! 5. restart
call df%save_state(basename='test_inflow_df')
call make(dg)
call dg%setup(emin=[0._R8P, 0._R8P, 0._R8P], emax=[1._R8P, 2._R8P, 2._R8P], periodic=[.false., .false., .true.])
call dg%load_state(basename='test_inflow_df')
do it=1, 50
   call df%advance(dt=dt)
   call dg%advance(dt=dt)
enddo
print '(A,L1,A)', 'restart: identical after 50 steps ', all(df%psi == dg%psi) .and. all(df%q_plane == dg%q_plane), ' (bitwise)'
ok = ok .and. all(df%psi == dg%psi) .and. all(df%q_plane == dg%q_plane)
open(newunit=it, file='test_inflow_df.dff') ; close(it, status='delete')
if (ok) then
   print '(A)', 'test_inflow_df: PASSED'
else
   print '(A)', 'test_inflow_df: FAILED'
endif
call mpih%finalize
if (.not.ok) error stop 1
contains
   subroutine make(g)
   !< A generator on the synthetic profile: Mach 2.5 free stream, delta 1, integral lengths 1 (time), 0.3, 0.25.
   type(flume_inflow_df_object), intent(inout) :: g !< Generator.
   integer(I4P)                                 :: r !< Row.
   real(R8P)                                    :: eta !< y / delta.

   g%is_active = .true.
   g%face = 1_I4P ; g%axis_n = 1_I4P ; g%axis_w = 2_I4P ; g%axis_s = 3_I4P ; g%sgn = 1._R8P
   g%gamma = 1.4_R8P ; g%r_inf = 1._R8P ; g%u_inf = 1._R8P ; g%p_inf = 1._R8P / (1.4_R8P * 2.5_R8P**2)
   g%delta = 1._R8P ; g%spacing = 0.05_R8P ; g%length = [1._R8P, 0.3_R8P, 0.25_R8P] ; g%seed = 7_I8P
   if (allocated(g%prof)) deallocate(g%prof)
   allocate(g%prof(7, 21))
   do r=1, 21
      eta = real(r - 1, R8P) / 20._R8P
      g%prof(:,r) = [eta, max(eta, 1.e-3_R8P)**(1._R8P / 7._R8P), 1._R8P + 0.6_R8P * (1._R8P - eta**2), &
                     0.010_R8P * (1._R8P - eta), 0.004_R8P * (1._R8P - eta), 0.006_R8P * (1._R8P - eta), -0.003_R8P * (1._R8P - eta)]
   enddo
   endsubroutine make

   function stress(c) result(r)
   !< The target stress `c` (1 uu, 2 vv, 3 ww, 4 uv) of every node, from the table itself (not the generator's factors).
   integer(I4P), intent(in) :: c        !< Component.
   real(R8P)                :: r(df%ny) !< Stresses.
   real(R8P)                :: u, t, rr(4) !< Table values.
   integer(I4P)             :: jj       !< Node.

   do jj=1, df%ny
      call profile_mean(prof=df%prof, eta=(real(jj, R8P) - 0.5_R8P) * df%dy / df%delta, u=u, t=t, r=rr)
      r(jj) = rr(c) * df%u_inf**2
   enddo
   endfunction stress
endprogram test_inflow_df
