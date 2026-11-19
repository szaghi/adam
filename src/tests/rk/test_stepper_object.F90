!< Unit test of the SELAH adapter `stepper_object` (SELAH M7, ADAM #132 box 1.2).
program test_stepper_object
!< Unit test of the SELAH adapter `stepper_object` (SELAH M7, ADAM #132 box 1.2).
!<
!< **Why this test exists.** `stepper_object` is ADAM's only door to SELAH: the scheme names of `[runge_kutta]` become a
!< SELAH scheme and storage form, ADAM's field and immersed-boundary mask become SELAH's box and mask, and the registers
!< are sized and owned here. A wrong mapping or a wrong box is a silent change of the time integrator, so each is pinned
!< here on its own, with no grid, tree or MPI run.
!<
!< **What it pins.**
!<   1. The mapping of every ADAM scheme name (`stepper_map`), and an unknown name maps to nothing.
!<   2. The stored fields (`stepper_stored_fields`) equal SELAH's `registers_needed`: 2 for the two-register forms, s + 1
!<      for the Butcher forms, 1 for the Yoshida splitting; every Runge-Kutta scheme has at least 2, so the D_flux
!<      diagnostic, which borrows registers 1 and 2 between steps, is accepted for every scheme (#128 D1: ADAM's
!<      low-storage schemes stored 1 stage while `[diagnostics] divb_flux` was accepted from 2 stages up).
!<   3. One step of y' = -y on an ADAM-shaped field (nv, ghosts, 2 of 3 pool blocks active) for every Runge-Kutta
!<      scheme, unmasked and with half the cells solid: the active cells equal the scheme's amplification factor (its
!<      Butcher tableau, evaluated in scalar arithmetic) to 4 ulp; solid cells, ghosts and the inactive block bitwise
!<      untouched; in the Butcher forms the step residual is -y times the factor's derivative part where active and
!<      exactly 0 in the solid cells (#128 D4), and the masked cells of every stage state are q at the stage start.
!<   4. Writing garbage into the scratch registers between two steps leaves the second step bitwise unchanged
!<      (#132 rule 3).
!<   5. The stage weights and stage times are the tableau's b and c.
!< The refusal of an unknown scheme (`mpih%error_stop`) is checked by test_stepper_object_xfail.
use adam_field_object,   only : field_object
use adam_grid_object,    only : grid_object
use adam_mpih_global,    only : mpih
use adam_stepper_object, only : stepper_object, stepper_map, stepper_stored_fields, &
                                RK_1, RK_2, RK_3, RK_SSP_11, RK_SSP_22, RK_SSP_33, RK_SSP_54, RK_YOSHIDA
use selah,               only : butcher, butcher_tableau, registers_needed, FORM_BUTCHER, FORM_SHU_OSHER_2R,        &
                                FORM_SPLITTING, SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, SCHEME_SSP_54,          &
                                SCHEME_YOSHIDA_4
use penf,                only : I4P, R8P

implicit none

character(19), parameter :: NAMES(8)  = [character(19) :: RK_1, RK_2, RK_3, RK_SSP_11, RK_SSP_22, RK_SSP_33, &
                                                          RK_SSP_54, RK_YOSHIDA] !< ADAM schemes.
character(9),  parameter :: SELAH(8)  = [character(9) :: SCHEME_EULER, SCHEME_SSP_22, SCHEME_SSP_33, SCHEME_EULER, &
                                                         SCHEME_SSP_22, SCHEME_SSP_33, SCHEME_SSP_54, SCHEME_YOSHIDA_4]
                                                                                 !< Expected SELAH schemes.
character(12), parameter :: FORMS(8)  = [character(12) :: FORM_SHU_OSHER_2R, FORM_SHU_OSHER_2R, FORM_SHU_OSHER_2R, &
                                                          FORM_BUTCHER, FORM_BUTCHER, FORM_BUTCHER, FORM_BUTCHER,    &
                                                          FORM_SPLITTING]       !< Expected forms.
integer(I4P),  parameter :: FIELDS(8) = [2, 2, 2, 2, 3, 4, 6, 1]                 !< Expected stored fields.
integer(I4P),  parameter :: NV=3, NI=4, NJ=4, NK=4, NGC=2, NB=3, BN=2           !< Field shape, active blocks.
real(R8P),     parameter :: DT=0.1_R8P                                           !< Time step.
character(:), allocatable :: sch, frm !< Mapped scheme and form.
integer(I4P)              :: m        !< Counter.
integer(I4P)              :: failures !< Failed checks.

mpih%myrankstr = ''
failures = 0
do m=1, size(NAMES)
   call stepper_map(scheme=NAMES(m), selah_scheme=sch, form=frm)
   call check(sch == trim(SELAH(m)) .and. frm == trim(FORMS(m)), 'map '//trim(NAMES(m))//' -> '//sch//', '//frm)
   call check(stepper_stored_fields(NAMES(m)) == FIELDS(m) .and. &
              stepper_stored_fields(NAMES(m)) == registers_needed(trim(SELAH(m)), trim(FORMS(m))), &
              'stored fields of '//trim(NAMES(m)))
   if (m <= 7) call check(stepper_stored_fields(NAMES(m)) >= 2, 'D_flux scratch (2 registers) for '//trim(NAMES(m)))
enddo
call stepper_map(scheme='runge-kutta-4', selah_scheme=sch, form=frm)
call check(len(sch) == 0 .and. len(frm) == 0, 'unknown scheme maps to nothing')
call check(stepper_stored_fields('runge-kutta-4') == -1, 'unknown scheme: -1 stored fields')
do m=1, 7
   call one_step(scheme=NAMES(m), tableau=butcher(trim(SELAH(m))), masked=.false.)
   call one_step(scheme=NAMES(m), tableau=butcher(trim(SELAH(m))), masked=.true.)
enddo
if (failures > 0) then
   print '(A,I0,A)', 'test_stepper_object: ', failures, ' check(s) FAILED'
   error stop 1
endif
print '(A)', 'test_stepper_object: ok'

contains
   subroutine check(ok, label)
   !< Report one check.
   logical,      intent(in) :: ok    !< Result.
   character(*), intent(in) :: label !< Label.

   if (ok) then
      print '(A)', '  ok   '//label
   else
      print '(A)', '  FAIL '//label
      failures = failures + 1
   endif
   endsubroutine check

   subroutine one_step(scheme, tableau, masked)
   !< One step of y' = -y on an ADAM-shaped field, checked cell by cell; then the scratch check.
   character(*),          intent(in) :: scheme  !< ADAM scheme.
   type(butcher_tableau), intent(in) :: tableau !< Its tableau.
   logical,               intent(in) :: masked  !< Half the cells solid.
   type(field_object)                :: field   !< Field.
   type(grid_object)                 :: grid    !< Grid.
   type(stepper_object)              :: ts      !< Stepper.
   real(R8P), allocatable            :: q0(:,:,:,:,:), q(:,:,:,:,:), dq(:,:,:,:,:), q2(:,:,:,:,:) !< Fields.
   real(R8P), allocatable            :: phi(:,:,:,:,:) !< IB distance (2 components, the last one the summary).
   real(R8P), pointer                :: r(:,:,:,:,:)   !< A register.
   real(R8P)                         :: g, d           !< Amplification factor, residual factor.
   real(R8P)                         :: y(tableau%stages) !< Scalar stage values.
   real(R8P)                         :: err, derr      !< Largest relative errors.
   integer(I4P)                      :: nbad, nmask    !< Cells changed outside the active set; masked-state misses.
   integer(I4P)                      :: s, j, v, i, jj, kk, b !< Counters.
   character(:), allocatable         :: label          !< Label.
   logical                           :: active         !< Cell in the box and fluid.

   label = trim(scheme)//merge(' masked  ', ' unmasked', masked)
   field%nv = NV ; field%nb = NB ; field%blocks_number = BN
   grid%ni = NI ; grid%nj = NJ ; grid%nk = NK ; grid%ngc = NGC
   allocate(q0(1:NV,1-NGC:NI+NGC,1-NGC:NJ+NGC,1-NGC:NK+NGC,1:NB))
   allocate(phi(1:2,1-NGC:NI+NGC,1-NGC:NJ+NGC,1-NGC:NK+NGC,1:NB))
   do b=1, NB ; do kk=1-NGC, NK+NGC ; do jj=1-NGC, NJ+NGC ; do i=1-NGC, NI+NGC
      do v=1, NV
         q0(v,i,jj,kk,b) = 1._R8P + 0.01_R8P * v + 0.1_R8P * i - 0.05_R8P * jj + 0.02_R8P * kk + 0.3_R8P * b
      enddo
      phi(:,i,jj,kk,b) = merge(1._R8P, -1._R8P, masked .and. mod(i + jj + kk, 2) == 0) ! solid where >= 0
   enddo ; enddo ; enddo ; enddo
   q = q0
   dq = q0 ! a value the Butcher step residual must overwrite in the box
   call ts%initialize(field=field, grid=grid, scheme=scheme)
   call take_step(ts, field, q, dq, phi, masked, nmask)
   ! the scalar reference: the tableau on y' = -y
   do s=1, tableau%stages
      y(s) = 1._R8P
      do j=1, s - 1
         y(s) = y(s) + DT * tableau%a(s,j) * (-y(j))
      enddo
   enddo
   g = 1._R8P ; d = 0._R8P
   do s=1, tableau%stages
      d = d + tableau%b(s) * (-y(s))
   enddo
   g = 1._R8P + DT * d
   err = 0._R8P ; derr = 0._R8P ; nbad = 0
   do b=1, NB ; do kk=1-NGC, NK+NGC ; do jj=1-NGC, NJ+NGC ; do i=1-NGC, NI+NGC ; do v=1, NV
      active = b <= BN .and. all([i, jj, kk] >= 1) .and. i <= NI .and. jj <= NJ .and. kk <= NK
      if (active .and. masked) active = phi(2,i,jj,kk,b) < 0._R8P
      if (active) then
         err = max(err, abs(q(v,i,jj,kk,b) - g * q0(v,i,jj,kk,b)) / abs(g * q0(v,i,jj,kk,b)))
         if (ts%is_butcher()) derr = max(derr, abs(dq(v,i,jj,kk,b) - d * q0(v,i,jj,kk,b)) / abs(d * q0(v,i,jj,kk,b)))
      else
         if (q(v,i,jj,kk,b) /= q0(v,i,jj,kk,b)) nbad = nbad + 1
         if (ts%is_butcher() .and. b <= BN .and. all([i, jj, kk] >= 1) .and. i <= NI .and. jj <= NJ .and. kk <= NK) then
            if (dq(v,i,jj,kk,b) /= 0._R8P) nbad = nbad + 1 ! masked: the step residual is 0 (#128 D4)
         endif
      endif
   enddo ; enddo ; enddo ; enddo ; enddo
   call check(err <= 4._R8P * epsilon(1._R8P) .and. derr <= 8._R8P * epsilon(1._R8P) .and. nbad == 0 .and. nmask == 0, &
              label//': step (q rel '//trim(rs(err))//', dq rel '//trim(rs(derr))//', cells outside changed '// &
              trim(is(nbad))//', masked stage-state misses '//trim(is(nmask))//')')
   ! stage weights and times
   do s=1, ts%stages()
      if (ts%stage_weight(s) /= tableau%b(s) .or. ts%stage_time(s) /= tableau%c(s)) nbad = nbad + 1
   enddo
   call check(nbad == 0 .and. ts%stages() == tableau%stages, label//': stage weights and times are b and c')
   ! rule 3: garbage in the scratch registers between two steps changes nothing
   q2 = q
   call take_step(ts, field, q2, dq, phi, masked, nmask)
   do j=1, 2
      r => ts%scratch(j)
      r = huge(1._R8P)
   enddo
   call take_step(ts, field, q, dq, phi, masked, nmask)
   call check(all(q == q2), label//': scratch registers written between steps, next step bitwise unchanged')
   call ts%destroy
   endsubroutine one_step

   subroutine take_step(ts, field, q, dq, phi, masked, nmask)
   !< One step of y' = -y through the adapter, the residual written into SELAH's register.
   type(stepper_object), intent(inout) :: ts                                      !< Stepper.
   type(field_object),   intent(in)    :: field                                   !< Field.
   real(R8P), target,    intent(inout) :: q(1:,1-NGC:,1-NGC:,1-NGC:,1:)           !< State.
   real(R8P),            intent(inout) :: dq(1:,1-NGC:,1-NGC:,1-NGC:,1:)          !< Step residual.
   real(R8P),            intent(in)    :: phi(1:,1-NGC:,1-NGC:,1-NGC:,1:)         !< IB distance.
   logical,              intent(in)    :: masked                                  !< Pass the mask.
   integer(I4P),         intent(out)   :: nmask                                   !< Masked stage-state misses.
   real(R8P), pointer                  :: u(:,:,:,:,:), r(:,:,:,:,:)              !< Stage state, residual.
   integer(I4P)                        :: s, i, jj, kk, b                         !< Counters.

   nmask = 0
   if (masked) then
      call ts%begin_step(field=field, q=q, phi=phi)
   else
      call ts%begin_step(field=field, q=q)
   endif
   do s=1, ts%stages()
      if (masked) then
         call ts%begin_stage(k=s, field=field, dt=DT, q=q, phi=phi)
      else
         call ts%begin_stage(k=s, field=field, dt=DT, q=q)
      endif
      if (ts%state_is_q(s)) then
         u(1:,1-NGC:,1-NGC:,1-NGC:,1:) => q
      else
         u => ts%stage_state()
         if (masked) then
            do b=1, BN ; do kk=1, NK ; do jj=1, NJ ; do i=1, NI
               if (phi(2,i,jj,kk,b) >= 0._R8P .and. any(u(:,i,jj,kk,b) /= q(:,i,jj,kk,b))) nmask = nmask + 1
            enddo ; enddo ; enddo ; enddo
         endif
      endif
      r => ts%rhs_register()
      r(:,1:NI,1:NJ,1:NK,1:BN) = -u(:,1:NI,1:NJ,1:NK,1:BN)
      if (masked) then
         call ts%end_stage(k=s, field=field, dt=DT, q=q, phi=phi)
      else
         call ts%end_stage(k=s, field=field, dt=DT, q=q)
      endif
   enddo
   if (masked) then
      call ts%finish_step(field=field, dt=DT, q=q, dq=dq, phi=phi)
   else
      call ts%finish_step(field=field, dt=DT, q=q, dq=dq)
   endif
   endsubroutine take_step

   function rs(x) result(s)
   !< A real as a short string.
   real(R8P), intent(in) :: x !< Value.
   character(12)         :: s !< String.

   write(s, '(ES10.3)') x
   endfunction rs

   function is(n) result(s)
   !< An integer as a short string.
   integer(I4P), intent(in) :: n !< Value.
   character(12)            :: s !< String.

   write(s, '(I0)') n
   endfunction is
endprogram test_stepper_object
