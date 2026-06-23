!< Unit test of the PRISM-CPU Berger-Colella reflux correction TBP.
program test_prism_reflux_apply
!< Unit test of `prism_cpu_object%apply_reflux_to_stage_forest` — the realm-side
!< body that holds the actual Berger-Colella reflux arithmetic for issue #13
!< Phase A.
!<
!< **Why this test exists.** The flux register's corrective path is, as of
!< develop, a no-op in every *reachable* PRISM configuration: inter-realm seams
!< are same-resolution mirrors (`F_coarse - F_fine_sum` is round-off zero) and
!< there is no intra-realm AMR registration pass yet. The register accumulators
!< therefore never carry an asymmetric flux in any regression case, so the
!< correction formula is exercised by *no* end-to-end test. This unit test
!< pins the formula directly by hand-feeding an asymmetric
!< `(F_coarse - F_fine_sum)` and asserting the exact per-cell q update.
!<
!< **What it pins** (verbatim from `apply_reflux_to_stage_forest`):
!<```
!< q_rk(:, i_c, j_c, k_c, coarse_block, stage) +=
!<     sgn * ark(stage) * (dt / dxyz(axis, coarse_block)) * (F_coarse(:,c,1) - F_fine_sum(:,c,1))
!<```
!< with `(axis, sgn)` from `face_axis_sign(coarse_face)`, the end-of-step gate
!< `stage == rk%nrk`, and the seam slice at the normal-axis extreme cell
!< (`i = ni` for X_MAX). It also asserts the interior is untouched and that an
!< earlier (non-final) stage is a no-op.
!<
!< This is a pure-arithmetic check: no INI, no solve, no IO — the realm's
!< `rk`/`grid`/`field` members are populated by hand to the minimum the TBP
!< reads.

use adam_prism_cpu_object, only : prism_cpu_object
use adam_flux_register_object, only : flux_register_object, SEAM_KIND_INTER_REALM
use adam_maps_object, only : FACE_X_MAX, face_axis_sign
use adam_globals, only : mpih
use penf, only : I4P, R8P

implicit none

type(prism_cpu_object)      :: realm        !< Realm under test (coarse side of the seam).
type(flux_register_object)  :: reg          !< Forest flux register holding one seam face.
integer(I4P), parameter     :: nv = 5_I4P   !< State-vector width.
integer(I4P), parameter     :: ni = 4_I4P   !< Cells, i.
integer(I4P), parameter     :: nj = 3_I4P   !< Cells, j.
integer(I4P), parameter     :: nk = 2_I4P   !< Cells, k.
integer(I4P), parameter     :: nb = 1_I4P   !< Blocks.
integer(I4P), parameter     :: nrk = 3_I4P  !< RK stages.
integer(I4P), parameter     :: cblk = 1_I4P !< Coarse block index.
real(R8P),    parameter     :: dt = 0.25_R8P                  !< Time step.
real(R8P),    parameter     :: dx = 0.05_R8P                  !< Coarse-side dx on the seam-normal axis.
real(R8P),    parameter     :: ark_final = 0.6_R8P            !< ark(nrk): the end-of-step accumulation weight.
real(R8P),    parameter     :: tol = 1.0e-14_R8P              !< Round-off comparison tolerance.
integer(I4P)                :: nface_cells                   !< Coarse cells on the seam face (= nj*nk for x-normal).
integer(I4P)                :: axis, sgn                      !< Decoded seam-face axis/sign.
integer(I4P)                :: i, j, k, v, c, c0              !< Counters.
real(R8P), allocatable      :: q_before(:,:,:,:,:,:)          !< Snapshot of q_rk(final stage) before the correction.
real(R8P), allocatable      :: expected(:)                    !< Expected per-(v,c) delta.
real(R8P)                   :: delta, got, want               !< Scratch for the assertions.
logical                     :: test_passed                    !< Aggregate pass flag.

test_passed = .true.
call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)

! --- Build the realm's minimal state by hand (only what the TBP reads). ---
realm%realm_index = 1_I4P
realm%rk%nrk = nrk
allocate(realm%rk%ark(1:nrk))
realm%rk%ark = [0.3_R8P, 0.45_R8P, ark_final]            ! arbitrary; only ark(nrk) is used by the end-of-step gate
realm%adam%grid%ni = ni
realm%adam%grid%nj = nj
realm%adam%grid%nk = nk
allocate(realm%adam%field%dxyz(1:3, 1:nb))
realm%adam%field%dxyz = -1.0_R8P                          ! poison: any axis but the seam-normal must be irrelevant
realm%adam%field%dxyz(1, cblk) = dx                       ! x-normal seam → axis 1 is the one read

! q_rk shape (nv, i, j, k, b, stage); the TBP writes interior cells (no ghosts needed for the seam slice).
allocate(realm%rk%q_rk(1:nv, 1:ni, 1:nj, 1:nk, 1:nb, 1:nrk))
call seed_q_rk(realm%rk%q_rk)
q_before = realm%rk%q_rk

! --- Register one X_MAX seam face on the coarse realm, with asymmetric fluxes. ---
call face_axis_sign(FACE_X_MAX, axis, sgn)               ! → axis=1, sgn=+1
nface_cells = nj * nk                                    ! x-normal face: tangential cells are (j,k)
call reg%initialize(nfaces=1_I4P)
call reg%register_face(face_index=1_I4P, seam_kind=SEAM_KIND_INTER_REALM, &
                       coarse_realm=1_I4P, coarse_block=cblk, coarse_face=FACE_X_MAX, &
                       fine_realm=2_I4P, nface_cells=nface_cells, nv=nv, n_stages=1_I4P)
! Asymmetric coarse vs fine flux → non-zero (F_coarse - F_fine_sum). Distinct per (v,c)
! so a transposed/mis-walked index would be caught.
allocate(expected(1:nv*nface_cells))
call feed_asymmetric_fluxes(reg, expected)

! --- Gate check: a non-final stage must be a no-op. ---
call realm%apply_reflux_to_stage_forest(stage=1_I4P, dt=dt, flux_register=reg)
if (maxval(abs(realm%rk%q_rk - q_before)) > tol) then
   write(*,'(A)') 'FAIL: non-final stage (1) modified q_rk — end-of-step gate broken'
   test_passed = .false.
else
   write(*,'(A)') 'PASS: non-final stage is a no-op (end-of-step gate holds)'
endif

! --- Apply the real correction at the final stage. ---
call realm%apply_reflux_to_stage_forest(stage=nrk, dt=dt, flux_register=reg)

! --- Assertion 1: the seam slice (i = ni for X_MAX, sgn>0) moved by exactly the formula. ---
block
   real(R8P) :: scale_
   logical   :: seam_ok
   seam_ok = .true.
   scale_ = real(sgn, R8P) * ark_final * dt / dx
   do c = 1_I4P, nface_cells
      c0 = c - 1_I4P
      j = 1_I4P + mod(c0, nj)
      k = 1_I4P + c0 / nj
      do v = 1_I4P, nv
         want = q_before(v, ni, j, k, cblk, nrk) + scale_ * expected((c-1)*nv + v)
         got  = realm%rk%q_rk(v, ni, j, k, cblk, nrk)
         if (abs(got - want) > tol * max(1.0_R8P, abs(want))) seam_ok = .false.
      enddo
   enddo
   if (seam_ok) then
      write(*,'(A)') 'PASS: seam-cell correction matches sgn*ark*dt/dx*(F_coarse - F_fine_sum)'
   else
      write(*,'(A)') 'FAIL: seam-cell correction does not match the Berger-Colella formula'
      test_passed = .false.
   endif
endblock

! --- Assertion 2: every interior cell (i /= ni) is untouched at the final stage. ---
block
   logical :: interior_ok
   interior_ok = .true.
   do k = 1_I4P, nk ; do j = 1_I4P, nj ; do i = 1_I4P, ni-1_I4P ; do v = 1_I4P, nv
      delta = abs(realm%rk%q_rk(v, i, j, k, cblk, nrk) - q_before(v, i, j, k, cblk, nrk))
      if (delta > tol) interior_ok = .false.
   enddo ; enddo ; enddo ; enddo
   if (interior_ok) then
      write(*,'(A)') 'PASS: interior cells (i /= ni) untouched by the seam correction'
   else
      write(*,'(A)') 'FAIL: interior cells were modified — seam slice addressing is wrong'
      test_passed = .false.
   endif
endblock

! --- Assertion 3: other RK stages (1, 2) untouched (write target was stage=nrk only). ---
block
   logical :: stages_ok
   stages_ok = .true.
   if (maxval(abs(realm%rk%q_rk(:,:,:,:,:,1) - q_before(:,:,:,:,:,1))) > tol) stages_ok = .false.
   if (maxval(abs(realm%rk%q_rk(:,:,:,:,:,2) - q_before(:,:,:,:,:,2))) > tol) stages_ok = .false.
   if (stages_ok) then
      write(*,'(A)') 'PASS: non-target RK stages (1,2) untouched'
   else
      write(*,'(A)') 'FAIL: a non-target RK stage was modified'
      test_passed = .false.
   endif
endblock

if (test_passed) then
   write(*,'(A)') 'TEST PASSED: prism reflux apply'
else
   write(*,'(A)') 'TEST FAILED: prism reflux apply'
endif

call reg%destroy
call mpih%finalize
if (.not. test_passed) error stop 1

contains
   subroutine seed_q_rk(q)
   !< Seed q_rk with a deterministic, cell-distinct pattern so any spurious
   !< write (wrong stage, wrong cell, wrong component) is detectable.
   real(R8P), intent(out) :: q(:,:,:,:,:,:) !< RK stage buffer (nv,i,j,k,b,stage).
   integer(I4P)           :: vv, ii, jj, kk, bb, ss !< Counters.

   do ss = 1, size(q,6) ; do bb = 1, size(q,5) ; do kk = 1, size(q,4)
      do jj = 1, size(q,3) ; do ii = 1, size(q,2) ; do vv = 1, size(q,1)
         q(vv,ii,jj,kk,bb,ss) = real(vv,R8P) + 0.1_R8P*ii + 0.01_R8P*jj + 0.001_R8P*kk + 10.0_R8P*ss
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   endsubroutine seed_q_rk

   subroutine feed_asymmetric_fluxes(register, exp_delta)
   !< Push distinct coarse and fine fluxes into the one registered face so that
   !< `(F_coarse - F_fine_sum)((c-1)*nv+v)` is a known, per-(v,c)-distinct value.
   type(flux_register_object), intent(inout) :: register      !< Register (face 1 already allocated).
   real(R8P),                  intent(out)    :: exp_delta(:)  !< Expected (F_coarse - F_fine_sum), flat (c-major, v-minor).
   real(R8P), allocatable :: fc(:,:), ff(:,:)                  !< Coarse / fine flux skins (nv, nface_cells).
   integer(I4P)           :: cc, vv                            !< Counters.

   allocate(fc(1:nv, 1:nface_cells), ff(1:nv, 1:nface_cells))
   do cc = 1, nface_cells ; do vv = 1, nv
      fc(vv,cc) = 2.0_R8P*vv + 0.5_R8P*cc        ! coarse-side flux
      ff(vv,cc) =        vv - 0.25_R8P*cc        ! fine-side summed flux (different → asymmetric)
      exp_delta((cc-1)*nv + vv) = fc(vv,cc) - ff(vv,cc)
   enddo ; enddo
   call register%accumulate_coarse_flux(face_index=1_I4P, stage=1_I4P, flux_face=fc)
   call register%accumulate_fine_flux  (face_index=1_I4P, stage=1_I4P, flux_face=ff)
   endsubroutine feed_asymmetric_fluxes
endprogram test_prism_reflux_apply
