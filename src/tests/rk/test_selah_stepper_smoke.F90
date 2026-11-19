!< Smoke test of the SELAH dependency: a SELAH stepper takes one step on an ADAM-shaped field (SELAH M7, ADAM #132 0.3).
program test_selah_stepper_smoke
!< Smoke test of the SELAH dependency: a SELAH stepper takes one step on an ADAM-shaped field (SELAH M7, ADAM #132 0.3).
!<
!< **Why this test exists.** ADAM compiles SELAH from source with its own flags (the host build in the CPU modes, the
!< OpenACC build in the FNL modes, `-DDEV_OAC -DCOMPILER_NVF`). This program is the first ADAM target that reaches
!< SELAH: it proves that the dependency is fetched, compiled and linked on each backend, before any ADAM object uses it.
!<
!< **What it does.** A field shaped as ADAM's `q(nv, 1-ngc:ni+ngc, 1-ngc:nj+ngc, 1-ngc:nk+ngc, nb)`, a pool of 3 blocks
!< of which 2 are active, allocated through FUNDAL (device memory in the FNL build). One SSP-33 step of y' = -y, in the
!< Butcher form and in the two-register form, the residual written on the device into SELAH's residual register.
!< Checks: every interior cell of the active blocks equals y0 (1 - dt + dt^2/2 - dt^3/6) to 4 ulp; every other cell (the
!< ghosts, the inactive block) is bitwise untouched.
#include "fundal.H"
use fundal, only : dev_alloc, dev_free, dev_memcpy_to_device, dev_memcpy_from_device
#if defined DEV_OAC
use fundal, only : dev_init
#endif
use selah,  only : stepper, box, box_range, FORM_BUTCHER, FORM_SHU_OSHER_2R, SCHEME_SSP_33, R_P
use penf,   only : I4P, R8P

implicit none

integer(I4P), parameter :: NV=5, NI=4, NJ=4, NK=4, NGC=2, NB=3, BN=2 !< Field shape, active blocks.
real(R8P),    parameter :: DT=0.1_R8P                                  !< Time step.
logical                 :: ok                                          !< Overall result.

if (R_P /= R8P) error stop 'test_selah_stepper_smoke: SELAH must be built with R8P in ADAM'
#if defined DEV_OAC
call dev_init(require_device=.true.)
print '(A)', 'FNL build: the field and the registers are device memory'
#else
print '(A)', 'CPU build: FUNDAL host memory'
#endif
ok = .true.
call one_step(form=FORM_BUTCHER,      label='ssp-33 butcher', ok=ok)
call one_step(form=FORM_SHU_OSHER_2R, label='ssp-33 two-register', ok=ok)
if (.not.ok) error stop 'test_selah_stepper_smoke: failed'
print '(A)', 'test_selah_stepper_smoke: ok'

contains
   subroutine one_step(form, label, ok)
   !< One SSP-33 step of y' = -y in the given storage form, checked cell by cell.
   character(*), intent(in)    :: form  !< Storage form.
   character(*), intent(in)    :: label !< Label of the report.
   logical,      intent(inout) :: ok    !< Overall result.
   type(stepper)               :: st    !< Stepper.
   type(box)                   :: bx    !< Interior of the active blocks.
   real(R8P), pointer          :: q_d(:,:,:,:,:) => null() !< Field, device memory in the FNL build.
   real(R8P), pointer          :: u(:,:,:,:,:)             !< Stage state.
   real(R8P), pointer          :: r(:,:,:,:,:)             !< Residual register.
   real(R8P), allocatable      :: q0(:,:,:,:,:), q1(:,:,:,:,:) !< Host copies before and after.
   real(R8P)                   :: g, err                   !< Amplification factor, largest relative error.
   integer(I4P)                :: lb(5), ub(5), lo(5), hi(5), k, ierr, v, i, j, kk, b, nbad !< Bounds, counters.

   lb = [1, 1-NGC, 1-NGC, 1-NGC, 1] ; ub = [NV, NI+NGC, NJ+NGC, NK+NGC, NB]
   lo = [1, 1, 1, 1, 1]             ; hi = [NV, NI, NJ, NK, BN]
   bx = box_range(lb=lb, lo=lo, hi=hi, components=.true.)
   allocate(q0(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3),lb(4):ub(4),lb(5):ub(5)))
   do b=lb(5), ub(5) ; do kk=lb(4), ub(4) ; do j=lb(3), ub(3) ; do i=lb(2), ub(2) ; do v=lb(1), ub(1)
      q0(v,i,j,kk,b) = 1._R8P + 0.01_R8P * v + 0.1_R8P * i - 0.05_R8P * j + 0.02_R8P * kk + 0.3_R8P * b
   enddo ; enddo ; enddo ; enddo ; enddo
   call dev_alloc(fptr_dev=q_d, ubounds=ub, lbounds=lb, ierr=ierr)
   if (ierr /= 0) error stop 'test_selah_stepper_smoke: dev_alloc failed'
   call dev_memcpy_to_device(dst=q_d, src=q0)
   call st%initialize(scheme=SCHEME_SSP_33, form=form, lb=lb, ub=ub)
   call st%begin_step(bx, q_d)
   do k=1, st%stages()
      call st%begin_stage(k, bx, DT, q_d)
      if (st%stage_state_is_q(k)) then
         u => q_d
      else
         u => st%stage_state()
      endif
      r => st%rhs_register()
      call residual(u=u, r=r)
      call st%end_stage(k, bx, DT, q_d)
   enddo
   call st%finish_step(bx, DT, q_d)
   allocate(q1, mold=q0)
   call dev_memcpy_from_device(dst=q1, src=q_d)
   g = 1._R8P - DT + DT**2 / 2._R8P - DT**3 / 6._R8P
   err = 0._R8P ; nbad = 0
   do b=lb(5), ub(5) ; do kk=lb(4), ub(4) ; do j=lb(3), ub(3) ; do i=lb(2), ub(2) ; do v=lb(1), ub(1)
      if (all([v,i,j,kk,b] >= lo) .and. all([v,i,j,kk,b] <= hi)) then
         err = max(err, abs(q1(v,i,j,kk,b) - g * q0(v,i,j,kk,b)) / abs(g * q0(v,i,j,kk,b)))
      elseif (q1(v,i,j,kk,b) /= q0(v,i,j,kk,b)) then
         nbad = nbad + 1
      endif
   enddo ; enddo ; enddo ; enddo ; enddo
   print '(A,ES10.3,A,I0)', '  '//label//': largest relative error ', err, ', cells outside the box changed: ', nbad
   if (err > 4._R8P * epsilon(1._R8P) .or. nbad > 0) then
      print '(A)', '  '//label//': FAILED'
      ok = .false.
   endif
   call st%destroy
   call dev_free(q_d)
   endsubroutine one_step

   subroutine residual(u, r)
   !< The residual of y' = -y over the interior of the active blocks, on the device in the FNL build.
   real(R8P), pointer, intent(in) :: u(:,:,:,:,:) !< Stage state.
   real(R8P), pointer, intent(in) :: r(:,:,:,:,:) !< Residual register.
   integer(I4P)                   :: v, i, j, k, b !< Counters.

   !$acc parallel loop collapse(5) deviceptr(u, r)
   do b=1, BN
      do k=1, NK
         do j=1, NJ
            do i=1, NI
               do v=1, NV
                  r(v,i,j,k,b) = -u(v,i,j,k,b)
               enddo
            enddo
         enddo
      enddo
   enddo
   endsubroutine residual
endprogram test_selah_stepper_smoke
