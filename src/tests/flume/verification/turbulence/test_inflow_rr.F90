!< The recycling-rescaling inflow of FLUME (issue #121, S5), on a synthetic block: leg tb3 (a) of check.sh.

program test_inflow_rr
!< One block of 24 x 8 x 8 cells (spacing 0.1, three ghost layers), the face x_min, the station 1.2 (12 cells), the span shift half
!< the span; a synthetic state whose density and velocity vary along y and z (and with the layer). Checks:
!< 1. the blend and the damping: W(0) = 0, W(0.2) = 1/2, W(1) = 1, the damping 1/2 at eta = 2;
!< 2. the plane: the cells listed are the 3 layers ending at the station, ghost layer g taking cell 13 - g, each (layer, j, k) once;
!<    the fluctuations have zero spanwise means (density, and the Favre velocity: sum rho u'' = 0) per layer and row;
!< 3. end to end with no rescaling (beta = delta ratio = 1): the ghost (1 - g, j, k) of the rows built from a crown map is the table's
!<    mean plus the damped fluctuations of cell (13 - g, j, k + 4) computed here from the synthetic state, to round-off;
!< 4. the rescaling: with a plane whose fluctuations are linear in y, beta = 1.5 and delta ratio = 2, the velocity fluctuation is
!<    beta ((1 - W) 1.5 + 2 W) times the plane's at y, damped; a source above the plane gives the table's mean exactly.

use :: adam_flume_inflow_rr_object, only : flume_inflow_rr_object, rr_blend, rr_damping
use :: adam_flume_inflow_df_object, only : profile_mean
use :: adam_mpih_global,            only : mpih
use :: penf,                        only : I4P, I8P, R8P

implicit none
integer(I4P), parameter      :: NI=24_I4P, NJ=8_I4P, NK=8_I4P, NG=3_I4P
real(R8P),    parameter      :: DX=0.1_R8P, PI=acos(-1._R8P)
type(flume_inflow_rr_object) :: rr                                   !< Inflow.
real(R8P)                    :: x_cell(1-NG:NI+NG,1), y_cell(1-NG:NJ+NG,1), z_cell(1-NG:NK+NG,1) !< Cell centres.
real(R8P)                    :: dxyz(3,1)                            !< Spacing.
real(R8P)                    :: q(5,1-NG:NI+NG,1-NG:NJ+NG,1-NG:NK+NG) !< Synthetic state.
real(R8P), allocatable       :: q_rec(:,:)                           !< Plane cells.
integer(I8P), allocatable    :: crown(:,:,:)                         !< Crown map.
real(R8P)                    :: err, e1, e2, rhot, srho, sru, um, tm, r4(4), rm, f(4), qe(5), vel(3), w, d, y
integer(I4P)                 :: i, j, k, g, m, c, ks, cnt(NG,NJ,NK)
logical                      :: ok, ok1

call mpih%initialize(do_mpi_init=.true., verbose=.false.)
ok = .true.
! 1. blend and damping
ok1 = abs(rr_blend(0._R8P)) < 1.e-15_R8P .and. abs(rr_blend(0.2_R8P) - 0.5_R8P) < 1.e-15_R8P .and. &
      abs(rr_blend(1._R8P) - 1._R8P) < 1.e-15_R8P .and. abs(rr_damping(2._R8P) - 0.5_R8P) < 1.e-15_R8P
print '(A,L1)', 'blend W(0) 0, W(0.2) 1/2, W(1) 1, damping(2) 1/2: ', ok1
ok = ok .and. ok1
! the block and the synthetic state
do i=1-NG, NI+NG
   x_cell(i,1) = (real(i, R8P) - 0.5_R8P) * DX
enddo
do j=1-NG, NJ+NG
   y_cell(j,1) = (real(j, R8P) - 0.5_R8P) * DX
enddo
do k=1-NG, NK+NG
   z_cell(k,1) = (real(k, R8P) - 0.5_R8P) * DX
enddo
dxyz(:,1) = DX
do k=1-NG, NK+NG
   do j=1-NG, NJ+NG
      do i=1-NG, NI+NG
         rhot = 1._R8P + 0.2_R8P * y_cell(j,1) + 0.05_R8P * sin(2._R8P * PI * real(k, R8P) / NK) + 0.01_R8P * i
         vel = [0.8_R8P + 0.3_R8P * y_cell(j,1) + 0.07_R8P * cos(2._R8P * PI * real(k, R8P) / NK + 0.3_R8P * j) + 0.002_R8P * i, &
                0.03_R8P * sin(2._R8P * PI * real(2 * k, R8P) / NK) * y_cell(j,1), 0.02_R8P * cos(2._R8P * PI * real(k, R8P) / NK + j)]
         q(:,i,j,k) = [rhot, rhot * vel, 2.5_R8P + 0.5_R8P * rhot * sum(vel * vel)]
      enddo
   enddo
enddo
call make(rr, beta=1._R8P, ratio=1._R8P)
call rr%setup(emin=[0._R8P, 0._R8P, 0._R8P], emax=[NI * DX, NJ * DX, NK * DX], periodic=[.false., .false., .true.], ngc=NG)
! 2. the plane
call rr%list_recycle_cells(nb=1_I4P, nc=[NI, NJ, NK], ngc=NG, x_cell=x_cell, y_cell=y_cell, z_cell=z_cell, dxyz=dxyz)
cnt = 0_I4P
ok1 = rr%nrec == NG * NJ * NK
do m=1, rr%nrec
   associate(cc=>rr%rec_cells(:,m), ix=>rr%rec_idx(:,m))
   ok1 = ok1 .and. cc(1) == 13_I4P - ix(1) .and. cc(2) == ix(2) .and. cc(3) == ix(3) .and. cc(4) == 1_I4P
   cnt(ix(1),ix(2),ix(3)) = cnt(ix(1),ix(2),ix(3)) + 1_I4P
   endassociate
enddo
ok1 = ok1 .and. all(cnt == 1_I4P)
allocate(q_rec(5, rr%nrec))
do m=1, rr%nrec
   q_rec(:,m) = q(:,rr%rec_cells(1,m),rr%rec_cells(2,m),rr%rec_cells(3,m))
enddo
call rr%build_plane(q_rec=q_rec)
e1 = 0._R8P ; e2 = 0._R8P
do g=1, NG
   do j=1, NJ
      e1 = max(e1, abs(sum(rr%fluct(1,g,j,:))))
      srho = 0._R8P ; sru = 0._R8P
      do k=1, NK
         srho = srho + q(1,13-g,j,k)
         sru = sru + q(1,13-g,j,k) * rr%fluct(2,g,j,k)
      enddo
      e2 = max(e2, abs(sru) / srho)
   enddo
enddo
ok1 = ok1 .and. e1 < 1.e-14_R8P .and. e2 < 1.e-14_R8P
print '(A,I0,A,ES9.2,A,ES9.2,A,L1)', 'plane: ', rr%nrec, ' cells, layer g <- cell 13 - g, each once; sum rho'' ', e1, &
      ', sum rho u'''' / sum rho ', e2, ': ', ok1
ok = ok .and. ok1
! 3. end to end, no rescaling: rows from a crown map of the face's three ghost layers
allocate(crown(NJ*NK, 9, NG), source=0_I8P)
do g=1, NG
   c = 0_I4P
   do k=1, NK
      do j=1, NJ
         c = c + 1_I4P
         crown(c,:,g) = [1_I8P, int(1 - g, I8P), int(j, I8P), int(k, I8P), 1_I8P, 0_I8P, 0_I8P, 10_I8P, 1_I8P]
      enddo
   enddo
enddo
call rr%build_rows(crown_map=crown, x_cell=x_cell, y_cell=y_cell, z_cell=z_cell)
err = 0._R8P
ok1 = rr%nrows == NG * NJ * NK
do m=1, rr%nrows
   associate(rw=>rr%rows(:,m))
   g = 1_I4P - rw(2) ; j = rw(3) ; k = rw(4)
   ks = modulo(k - 1_I4P + NK / 2_I4P, NK) + 1_I4P
   y = y_cell(j,1)
   call profile_mean(prof=rr%prof, eta=y, u=um, t=tm, r=r4)
   rm = rr%r_inf / tm
   d = rr_damping(y)
   ! the fluctuations of cell (13 - g, j, ks), from the synthetic state
   srho = 0._R8P ; vel = 0._R8P
   do i=1, NK
      srho = srho + q(1,13-g,j,i)
      vel = vel + q(2:4,13-g,j,i)
   enddo
   f(1) = q(1,13-g,j,ks) - srho / NK
   f(2:4) = q(2:4,13-g,j,ks) / q(1,13-g,j,ks) - vel / srho
   f = f * d
   vel = f(2:4) ; vel(1) = vel(1) + um * rr%u_inf
   qe = [rm + f(1), (rm + f(1)) * vel, rr%p_inf / (rr%gamma - 1._R8P) + 0.5_R8P * (rm + f(1)) * sum(vel * vel)]
   err = max(err, maxval(abs(rr%q_rows(:,m) - qe)))
   endassociate
enddo
ok1 = ok1 .and. err < 1.e-13_R8P
print '(A,I0,A,ES9.2,A,L1)', 'ghost rows ', rr%nrows, ', ghost (1-g, j, k) against the mean plus the damped fluctuations of cell (13-g, j, k+4): ', &
      err, ': ', ok1
ok = ok .and. ok1
! 4. the rescaling: fluctuations linear in y (rho' 0, u'' = y, v'' = 2 y, w'' = -y)
call make(rr, beta=1.5_R8P, ratio=2._R8P)
call rr%setup(emin=[0._R8P, 0._R8P, 0._R8P], emax=[NI * DX, NJ * DX, NK * DX], periodic=[.false., .false., .true.], ngc=NG)
do j=1, NJ
   rr%fluct(1,:,j,:) = 0._R8P
   rr%fluct(2,:,j,:) = y_cell(j,1)
   rr%fluct(3,:,j,:) = 2._R8P * y_cell(j,1)
   rr%fluct(4,:,j,:) = -y_cell(j,1)
enddo
err = 0._R8P
do j=1, 3 ! heights 0.05, 0.15, 0.25: sources 2 y up to 0.5, within the plane (centres 0.05 to 0.75)
   y = y_cell(j,1)
   call profile_mean(prof=rr%prof, eta=y, u=um, t=tm, r=r4)
   rm = rr%r_inf / tm
   w = rr_blend(y) ; d = rr_damping(y)
   vel = 1.5_R8P * ((1._R8P - w) * 1.5_R8P + w * 2._R8P) * y * d * [1._R8P, 2._R8P, -1._R8P]
   vel(1) = vel(1) + um * rr%u_inf
   qe = [rm, rm * vel, rr%p_inf / (rr%gamma - 1._R8P) + 0.5_R8P * rm * sum(vel * vel)]
   err = max(err, maxval(abs(rr%ghost_state(g=1_I4P, y=y, kz=1_I4P) - qe)))
enddo
y = 0.45_R8P ! the outer source 0.9 lies above the plane's last centre, 0.75: the table's mean exactly
call profile_mean(prof=rr%prof, eta=y, u=um, t=tm, r=r4)
rm = rr%r_inf / tm
qe = [rm, rm * um * rr%u_inf, 0._R8P, 0._R8P, rr%p_inf / (rr%gamma - 1._R8P) + 0.5_R8P * rm * (um * rr%u_inf)**2]
e1 = maxval(abs(rr%ghost_state(g=1_I4P, y=y, kz=1_I4P) - qe))
ok1 = err < 1.e-13_R8P .and. e1 == 0._R8P
print '(A,ES9.2,A,ES9.2,A,L1)', 'rescaling beta 1.5, delta ratio 2 on linear fluctuations: ', err, '; a source above the plane: ', e1, ': ', ok1
ok = ok .and. ok1
if (ok) then
   print '(A)', 'test_inflow_rr: PASS'
else
   print '(A)', 'test_inflow_rr: FAIL'
endif
call mpih%finalize
if (.not.ok) error stop 1
contains
   subroutine make(g, beta, ratio)
   !< An inflow on a synthetic profile: Mach 2.5 free stream, delta 1, the station 1.2, shift half the span.
   type(flume_inflow_rr_object), intent(inout) :: g     !< Inflow.
   real(R8P),                    intent(in)    :: beta  !< Viscous length ratio.
   real(R8P),                    intent(in)    :: ratio !< Thickness ratio.
   integer(I4P)                                :: r     !< Row.
   real(R8P)                                   :: eta   !< y / delta.

   g%is_active = .true.
   g%face = 1_I4P ; g%axis_n = 1_I4P ; g%axis_w = 2_I4P ; g%axis_s = 3_I4P ; g%sgn = 1._R8P
   g%gamma = 1.4_R8P ; g%r_inf = 1._R8P ; g%u_inf = 1._R8P ; g%p_inf = 1._R8P
   g%delta = 1._R8P ; g%spacing = DX ; g%x_r = 1.2_R8P ; g%beta = beta ; g%delta_ratio = ratio ; g%shift = 0.5_R8P
   if (allocated(g%prof)) deallocate(g%prof)
   allocate(g%prof(7, 21))
   do r=1, 21
      eta = real(r - 1, R8P) / 20._R8P
      g%prof(:,r) = [eta, max(eta, 1.e-3_R8P)**(1._R8P / 7._R8P), 1._R8P + 0.6_R8P * (1._R8P - eta**2), 0._R8P, 0._R8P, 0._R8P, 0._R8P]
   enddo
   endsubroutine make
endprogram test_inflow_rr
