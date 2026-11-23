!< Unit test KM-1 of the local Kolmogorov resolution marker `r = h / eta` (issue #133, phase 2).
program test_flume_kolmogorov_marker
!< Unit test KM-1 of the local Kolmogorov resolution marker `r = h / eta` (issue #133, phase 2).
!<
!< **Why this test exists.** The marker (`[amr_marker_N] mode = 7`) refines a block where its largest
!< `r = h (2 S_d:S_d)^(1/4) / nu^(1/2)` exceeds `refine_tol`. A wrong factor in the strain rate (the full gradient instead
!< of its symmetric part, the trace left in, a missing 2), a wrong spacing (the null axis in the mean, the arithmetic
!< instead of the geometric mean) or a viscosity law applied the wrong way still refines somewhere, and only shifts the
!< grid; here each one fails a closed form.
!<
!< **What it pins**, through `kolmogorov_block_max` (the estimator the FLUME marker calls on every block, with the
!< gradient marker's stencil) and the library's parser:
!< 1. the exact Couette flow `u = (S y, 0, 0)`: `r^2 = h^2 |S| / nu`, `h = (dx dy dz)^(1/3)`, on cubic and anisotropic
!<    cells, with every face closed (one-sided differences at the block ends) and every face open (ghosts filled with
!<    the exact field): a linear field is differentiated exactly, relative error below 1e-12;
!< 2. the same shear along another pair of axes, `u = (0, 0, S x)`, and on a 2-D block (null z): `h = (dx dy)^(1/2)`,
!<    the null extent left out;
!< 3. the viscosity law: `nu = mu0 (T / T_ref)^omega / rho` (a power law, T and rho away from 1);
!< 4. the deviatoric part: a solid-body rotation and a pure dilatation give `r = 0` (vorticity and compression are not
!<    dissipation; `S_d:S_d` is then a difference of round-offs whose fourth root `r` takes, so `r` is compared with the
!<    shear scale `h (S / nu)^(1/2)` within 1e-3), and a parabolic profile (open faces) takes its maximum at the cells next to the block ends;
!< 5. the hysteresis rule `0 < derefine_tol < refine_tol / 2` (`kolmogorov_tolerances_admissible`): refused at
!<    `refine_tol / 2` and above, at 0 and below; accepted above 1 (no Loehner cap), and parsed by
!<    `amr_object%load_from_file` from an INI source with `mode = 7`.

use :: adam_amr_object,               only : amr_object, AMR_KOLMOGOROV, kolmogorov_tolerances_admissible
use :: adam_flume_kolmogorov_library, only : kolmogorov_block_max, kolmogorov_spacing
use :: adam_flume_parameters,         only : IA_P, IA_R, IA_T, IA_U, IA_V, IA_W, NV_AUX
use :: adam_globals,                  only : mpih
use :: finer,                         only : file_ini
use :: penf,                          only : I4P, R8P, str

implicit none

integer(I4P), parameter :: NI=6_I4P, NJ=8_I4P, NK=4_I4P !< Block cells.
integer(I4P), parameter :: NGC=3_I4P                     !< Ghost cells.
real(R8P),    parameter :: TOL=1.e-12_R8P                !< Relative tolerance of the closed forms.
real(R8P),    parameter :: TOL_ZERO=1.e-3_R8P            !< Tolerance of a zero r, relative to the shear scale.
real(R8P),    parameter :: SHEAR=2.5_R8P                 !< Shear rate.
real(R8P),    parameter :: MU0=1.e-3_R8P                 !< Viscosity.
character(1), parameter :: NL=new_line('a')            !< New line.
logical,      parameter :: ALL_ACTIVE(3)=.true.          !< 3-D block.
logical,      parameter :: CLOSED(2,3)=.false.           !< Every face on the realm boundary.
logical,      parameter :: OPEN_(2,3)=.true.             !< Every face internal.
real(R8P)               :: qa(NV_AUX,1-NGC:NI+NGC,1-NGC:NJ+NGC,1-NGC:NK+NGC) !< Auxiliary variables of the block.
real(R8P)               :: dxyz(3)                       !< Cell spacings.
real(R8P)               :: h                             !< Spacing of the estimator.
real(R8P)               :: r, r_exact                    !< Estimator and its closed form.
real(R8P)               :: nu                            !< Kinematic viscosity.
type(amr_object)        :: amr                           !< Parser under test.
type(file_ini)          :: ini                           !< INI source.
integer(I4P)            :: t                             !< Counter.
logical                 :: test_passed                   !< Aggregate pass flag.

test_passed = .true.
call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)

! 1. Couette, cubic and anisotropic cells, closed and open faces
do t=1, 2
   if (t == 1) dxyz = [0.1_R8P, 0.1_R8P, 0.1_R8P]
   if (t == 2) dxyz = [0.3_R8P, 0.02_R8P, 0.07_R8P]
   call set_state(dxyz=dxyz, kind='couette-y', rho=1._R8P, temp=1._R8P)
   h = (dxyz(1) * dxyz(2) * dxyz(3))**(1._R8P / 3._R8P)
   r_exact = h * sqrt(SHEAR / MU0)
   call check(label='couette u(y), cells '//cells(dxyz)//', faces closed', &
              r=block_max(active=ALL_ACTIVE, is_open=CLOSED, omega=0._R8P), r_exact=r_exact)
   call check(label='couette u(y), cells '//cells(dxyz)//', faces open',   &
              r=block_max(active=ALL_ACTIVE, is_open=OPEN_, omega=0._R8P), r_exact=r_exact)
enddo

! 2. shear along another pair of axes; a 2-D block (null z): the null extent is not in h
dxyz = [0.3_R8P, 0.02_R8P, 0.07_R8P]
call set_state(dxyz=dxyz, kind='couette-x', rho=1._R8P, temp=1._R8P)
call check(label='shear w(x), anisotropic cells', r=block_max(active=ALL_ACTIVE, is_open=CLOSED, omega=0._R8P), &
           r_exact=(dxyz(1) * dxyz(2) * dxyz(3))**(1._R8P / 3._R8P) * sqrt(SHEAR / MU0))
dxyz = [0.3_R8P, 0.02_R8P, 100._R8P]
call set_state(dxyz=dxyz, kind='couette-y', rho=1._R8P, temp=1._R8P)
h = kolmogorov_spacing(dxyz=dxyz, active=[.true., .true., .false.])
if (abs(h - sqrt(dxyz(1) * dxyz(2))) > TOL * h) then
   print '(A)', 'FAILED: 2-D spacing '//trim(str(h))//', expected sqrt(dx dy) = '//trim(str(sqrt(dxyz(1) * dxyz(2))))
   test_passed = .false.
endif
call check(label='couette u(y), 2-D block (null z, dz = 100)', &
           r=block_max(active=[.true., .true., .false.], is_open=CLOSED, omega=0._R8P), r_exact=h * sqrt(SHEAR / MU0))

! 3. the viscosity law: power law at T = 3, T_ref = 1 (block_max), rho = 0.4
dxyz = [0.3_R8P, 0.02_R8P, 0.07_R8P]
call set_state(dxyz=dxyz, kind='couette-y', rho=0.4_R8P, temp=3._R8P)
nu = MU0 * 3._R8P**0.7_R8P / 0.4_R8P
call check(label='couette u(y), power law 0.7, T = 3, rho = 0.4', &
           r=block_max(active=ALL_ACTIVE, is_open=OPEN_, omega=0.7_R8P), &
           r_exact=(dxyz(1) * dxyz(2) * dxyz(3))**(1._R8P / 3._R8P) * sqrt(SHEAR / nu))

! 4. the deviatoric strain rate: rotation and dilatation dissipate nothing; a parabola peaks next to the block ends
! (S_d:S_d is then a difference of round-offs, and r takes its fourth root: compared with the shear scale at TOL_ZERO)
h = (dxyz(1) * dxyz(2) * dxyz(3))**(1._R8P / 3._R8P)
call set_state(dxyz=dxyz, kind='rotation', rho=1._R8P, temp=1._R8P)
call check(label='solid-body rotation', r=block_max(active=ALL_ACTIVE, is_open=CLOSED, omega=0._R8P), r_exact=0._R8P, &
           scale=h * sqrt(SHEAR / MU0))
call set_state(dxyz=dxyz, kind='dilatation', rho=1._R8P, temp=1._R8P)
call check(label='pure dilatation', r=block_max(active=ALL_ACTIVE, is_open=CLOSED, omega=0._R8P), r_exact=0._R8P, &
           scale=h * sqrt(SHEAR / MU0))
call set_state(dxyz=dxyz, kind='parabola', rho=1._R8P, temp=1._R8P)
! u = SHEAR y^2 / 2 with y the cell centre from the block centre: |du/dy| = SHEAR |y|, largest at the end cells
call check(label='parabola u(y), faces open', r=block_max(active=ALL_ACTIVE, is_open=OPEN_, omega=0._R8P), &
           r_exact=(dxyz(1) * dxyz(2) * dxyz(3))**(1._R8P / 3._R8P) * &
                   sqrt(SHEAR * (0.5_R8P * real(NJ - 1, R8P) * dxyz(2)) / MU0))

! 5. the hysteresis rule, and the parser
call admissible(refine_tol=4._R8P, derefine_tol=1.9_R8P, expected=.true.)
call admissible(refine_tol=4._R8P, derefine_tol=2._R8P, expected=.false.)
call admissible(refine_tol=4._R8P, derefine_tol=3._R8P, expected=.false.)
call admissible(refine_tol=4._R8P, derefine_tol=0._R8P, expected=.false.)
call admissible(refine_tol=4._R8P, derefine_tol=-1._R8P, expected=.false.)
call admissible(refine_tol=0.5_R8P, derefine_tol=0.2_R8P, expected=.true.)
call admissible(refine_tol=40._R8P, derefine_tol=15._R8P, expected=.true.)
call ini%load(source='[amr]'//NL//'frequency = 10'//NL//'[amr_marker_1]'//NL//'mode = 7'//NL//'delta_type = max'//NL// &
                     'delta_fine = 0.0'//NL//'delta_coarse = 0.0'//NL//'refine_tol = 4.0'//NL//'derefine_tol = 1.5'//NL)
call amr%load_from_file(file_parameters=ini)
if (amr%markers_number /= 1_I4P) then
   print '(A)', 'FAILED: parser, markers number '//trim(str(amr%markers_number))
   test_passed = .false.
elseif (amr%markers(1)%mode /= AMR_KOLMOGOROV .or. amr%markers(1)%refine_tol /= 4._R8P .or. &
        amr%markers(1)%derefine_tol /= 1.5_R8P) then
   print '(A)', 'FAILED: parser, mode '//trim(str(amr%markers(1)%mode))//', refine_tol '// &
                trim(str(amr%markers(1)%refine_tol))//', derefine_tol '//trim(str(amr%markers(1)%derefine_tol))
   test_passed = .false.
else
   print '(A)', 'ok: parser, mode 7 with refine_tol 4 and derefine_tol 1.5'
endif

if (test_passed) then
   print '(A)', 'TEST PASSED: Kolmogorov resolution marker (KM-1)'
else
   print '(A)', 'TEST FAILED: Kolmogorov resolution marker (KM-1)'
endif
call mpih%finalize
if (.not.test_passed) error stop 1
contains
   subroutine set_state(dxyz, kind, rho, temp)
   !< Fill every cell of the block, ghosts included, with an analytic velocity field, uniform density and temperature.
   real(R8P),    intent(in) :: dxyz(3) !< Cell spacings.
   character(*), intent(in) :: kind    !< Velocity field.
   real(R8P),    intent(in) :: rho     !< Density.
   real(R8P),    intent(in) :: temp    !< Temperature.
   real(R8P)                :: x(3)    !< Cell centre, from the block centre.
   integer(I4P)             :: i, j, k !< Counters.

   qa = 0._R8P
   do k=1-NGC, NK+NGC
      do j=1-NGC, NJ+NGC
         do i=1-NGC, NI+NGC
            x = ([real(i, R8P), real(j, R8P), real(k, R8P)] - 0.5_R8P * ([NI, NJ, NK] + 1)) * dxyz
            qa(IA_R,i,j,k) = rho
            qa(IA_T,i,j,k) = temp
            qa(IA_P,i,j,k) = rho * temp
            select case(kind)
            case('couette-y')  ; qa(IA_U,i,j,k) = SHEAR * x(2) + 0.3_R8P
            case('couette-x')  ; qa(IA_W,i,j,k) = SHEAR * x(1) - 1.1_R8P
            case('rotation')   ; qa(IA_U,i,j,k) = -SHEAR * x(2) ; qa(IA_V,i,j,k) = SHEAR * x(1)
            case('dilatation') ; qa(IA_U:IA_W,i,j,k) = SHEAR * x
            case('parabola')   ; qa(IA_U,i,j,k) = 0.5_R8P * SHEAR * x(2)**2
            endselect
         enddo
      enddo
   enddo
   endsubroutine set_state

   function block_max(active, is_open, omega) result(r)
   !< The estimator of the block under the current state, T_ref = 1, no Sutherland term.
   logical,   intent(in) :: active(3)    !< Active axes.
   logical,   intent(in) :: is_open(2,3) !< Internal faces.
   real(R8P), intent(in) :: omega        !< Exponent of the viscosity law.
   real(R8P)             :: r            !< Largest r of the block.

   r = kolmogorov_block_max(nv=NV_AUX, ni=NI, nj=NJ, nk=NK, ngc=NGC, qa=qa, dxyz=dxyz, active=active, is_open=is_open, &
                            mu0=MU0, tref=1._R8P, suth=0._R8P, omega_mu=omega)
   endfunction block_max

   subroutine check(label, r, r_exact, scale)
   !< Compare the estimator with its closed form: relative error below TOL, or with `scale` (a zero closed form) the
   !< error over the scale below TOL_ZERO.
   character(*), intent(in)           :: label   !< Case.
   real(R8P),    intent(in)           :: r       !< Estimator.
   real(R8P),    intent(in)           :: r_exact !< Closed form.
   real(R8P),    intent(in), optional :: scale   !< Scale of a zero closed form.
   real(R8P)                          :: err     !< Error.
   real(R8P)                          :: tol_    !< Tolerance.

   if (present(scale)) then
      err = abs(r - r_exact) / scale ; tol_ = TOL_ZERO
   else
      err = abs(r - r_exact) / abs(r_exact) ; tol_ = TOL
   endif
   if (err > tol_) then
      print '(A)', 'FAILED: '//label//': r = '//trim(str(r))//', exact '//trim(str(r_exact))//', error '//trim(str(err))
      test_passed = .false.
   else
      print '(A)', 'ok: '//label//': r = '//trim(str(r))//' (exact '//trim(str(r_exact))//')'
   endif
   endsubroutine check

   subroutine admissible(refine_tol, derefine_tol, expected)
   !< Check the hysteresis rule on one pair.
   real(R8P), intent(in) :: refine_tol, derefine_tol !< Thresholds.
   logical,   intent(in) :: expected                 !< Admissible.

   if (kolmogorov_tolerances_admissible(refine_tol=refine_tol, derefine_tol=derefine_tol) .neqv. expected) then
      print '(A)', 'FAILED: thresholds refine '//trim(str(refine_tol))//', derefine '//trim(str(derefine_tol))// &
                   merge(' refused, expected accepted', ' accepted, expected refused', expected)
      test_passed = .false.
   else
      print '(A)', 'ok: thresholds refine '//trim(str(refine_tol))//', derefine '//trim(str(derefine_tol))// &
                   merge(' accepted', ' refused ', expected)
   endif
   endsubroutine admissible

   function cells(dxyz) result(s)
   !< Return the spacings as text.
   real(R8P), intent(in)     :: dxyz(3) !< Spacings.
   character(:), allocatable :: s       !< Text.

   s = trim(str(dxyz(1), compact=.true.))//' x '//trim(str(dxyz(2), compact=.true.))//' x '// &
       trim(str(dxyz(3), compact=.true.))
   endfunction cells
endprogram test_flume_kolmogorov_marker
