!< Unit test of the PRISM-CPU shared-induction-flux correction (issue #13 Phase B).
program test_prism_induction_flux_share
!< Unit test of `prism_cpu_object%apply_induction_flux_sharing_forest` — the
!< realm-side body that makes the magnetic-field (induction) face flux *canonical*
!< across an inter-realm seam, the cell-centered analogue of constrained-transport
!< edge-E sharing (issue #13 Phase B).
!<
!< **Why this exists, and why it is shaped like the Phase A reflux test.**
!< PRISM has no edge-centered electric field and no UCT upwinding coefficients;
!< its induction update is a finite-volume flux divergence of a face-centered
!< Maxwell flux, `dq(B) = -(flx_f(i) - flx_f(i-1))/dx - ...`. So the route to
!< algebraic ∇·B = 0 at a seam is NOT a separate edge register (issue §4.3) but
!< a new *consumption mode* of the existing `flux_register_object`: the B-rows
!< (VAR_BX/BY/BZ) of the per-face flux are already accumulated there by Phase A's
!< `accumulate_seam_fluxes_fv`. Phase B reads them back as a single canonical
!< value used by BOTH sides of the seam, so the discrete ∇·B telescopes across
!< the seam exactly as it does across an interior face.
!<
!< Like reflux, the correction is round-off zero on a same-resolution mirror seam
!< (both sides reconstruct the identical face flux); it does real work only at a
!< non-conforming seam or intra-realm AMR jump, neither of which exists in the
!< codebase yet. This test therefore pins the *mechanism* on a hand-fed asymmetric
!< register, with no INI / solve / IO.
!<
!< **What it pins:**
!<   - canonical value = ½(F_coarse_B + F_fine_sum_B) per seam cell, per B component;
!<   - after the call, the coarse-side seam `flx_f(VAR_B*, ...)` equals that canonical
!<     value (the write-back happened on the B rows only);
!<   - non-B rows of `flx_f` at the seam are untouched (D, cleaning scalars, currents);
!<   - interior face fluxes (i /= seam) are untouched;
!<   - a SYMMETRIC register (F_coarse == F_fine_sum) is a no-op (mirror-seam property).

use adam_prism_cpu_object, only : prism_cpu_object
use adam_flux_register_object, only : flux_register_object, SEAM_KIND_INTER_REALM
use adam_prism_physics_object, only : VAR_BX, VAR_BY, VAR_BZ
use adam_maps_object, only : FACE_X_MAX, face_axis_sign
use adam_globals, only : mpih
use penf, only : I4P, R8P

implicit none

type(prism_cpu_object)      :: realm        !< Realm under test (coarse side of the seam).
type(flux_register_object)  :: reg          !< Flux register holding one seam face with B-row data.
integer(I4P), parameter     :: nv  = 9_I4P  !< State-vector width (DX..DZ,BX..BZ,Jx..Jz; no cleaning).
integer(I4P), parameter     :: nvc = 6_I4P  !< Conservative width (D and B).
integer(I4P), parameter     :: ni = 4_I4P   !< Cells, i.
integer(I4P), parameter     :: nj = 3_I4P   !< Cells, j.
integer(I4P), parameter     :: nk = 2_I4P   !< Cells, k.
integer(I4P), parameter     :: nb = 1_I4P   !< Blocks.
integer(I4P), parameter     :: ngc = 3_I4P  !< Ghost cells (flx_f is allocated 0:ni etc., ngc only affects q).
integer(I4P), parameter     :: cblk = 1_I4P !< Coarse block index.
integer(I4P), parameter     :: fec_xmax = 2_I4P !< +x face code in the maps fec 1..6 numbering.
real(R8P),    parameter     :: tol = 1.0e-14_R8P              !< Round-off comparison tolerance.
integer(I4P)                :: nface_cells                   !< Coarse cells on the seam face (= nj*nk for x-normal).
integer(I4P)                :: axis, sgn                      !< Decoded seam-face axis/sign.
integer(I4P)                :: i, j, k, v, c                  !< Counters.
integer(I4P), parameter     :: bvars(3) = [VAR_BX, VAR_BY, VAR_BZ] !< The induction (B) rows.
real(R8P), allocatable      :: flx_before(:,:,:,:,:)          !< Snapshot of flx_f before the call.
real(R8P), allocatable      :: canon(:,:)                     !< Expected canonical B-flux (3, nface_cells).
logical                     :: test_passed                    !< Aggregate pass flag.

test_passed = .true.
call mpih%initialize(do_mpi_init=.true., do_device_init=.false.)

! --- Minimal realm state the TBP reads. ---
realm%realm_index = 1_I4P
realm%rk%nrk = 3_I4P
realm%adam%grid%ni = ni
realm%adam%grid%nj = nj
realm%adam%grid%nk = nk
! NB: nv_c on prism_common_object is a POINTER (associated to physics%nv_c during
! init); the Phase-B TBP does not read it (it loops the fixed B rows), so we leave
! it null here rather than dereference an unassociated pointer. The TBP also walks
! flux_register%face(:) directly, gated on coarse_realm == realm_index, so it never
! consults inter_realm_face_register_index — no map needs to be built for this test.

! flx_f shape mirrors the residual routine: (nv, 0:ni, 1:nj, 1:nk, nb).
allocate(realm%flx_f(1:nv, 0:ni, 1:nj, 1:nk, 1:nb))
call seed_flx(realm%flx_f)
flx_before = realm%flx_f

! --- One +x seam face, with asymmetric B-row fluxes coarse vs fine. ---
call face_axis_sign(FACE_X_MAX, axis, sgn)               ! axis=1, sgn=+1
nface_cells = nj * nk
call reg%initialize(nfaces=1_I4P)
call reg%register_face(face_index=1_I4P, seam_kind=SEAM_KIND_INTER_REALM, &
                       coarse_realm=1_I4P, coarse_block=cblk, coarse_face=FACE_X_MAX, &
                       fine_realm=2_I4P, nface_cells=nface_cells, nv=nv, n_stages=1_I4P)
allocate(canon(1:3, 1:nface_cells))
call feed_asymmetric_B_rows(reg, canon)

! --- Apply the shared-induction-flux write-back at the final stage. ---
call realm%apply_induction_flux_sharing_forest(stage=realm%rk%nrk, flux_register=reg)

! --- Assertion 1: coarse-side seam flx_f(VAR_B*, ni, j, k) == canonical value. ---
block
   logical :: bok
   integer(I4P) :: vb
   bok = .true.
   do c = 1_I4P, nface_cells
      j = 1_I4P + mod(c-1_I4P, nj)
      k = 1_I4P + (c-1_I4P) / nj
      do vb = 1_I4P, 3_I4P
         if (abs(realm%flx_f(bvars(vb), ni, j, k, cblk) - canon(vb, c)) > &
             tol * max(1.0_R8P, abs(canon(vb, c)))) bok = .false.
      enddo
   enddo
   if (bok) then
      write(*,'(A)') 'PASS: seam B-flux replaced by canonical 1/2(F_coarse + F_fine_sum)'
   else
      write(*,'(A)') 'FAIL: seam B-flux is not the canonical shared value'
      test_passed = .false.
   endif
endblock

! --- Assertion 2: non-B rows at the seam face are untouched. ---
block
   logical :: nb_ok
   nb_ok = .true.
   do k = 1_I4P, nk ; do j = 1_I4P, nj ; do v = 1_I4P, nv
      if (v == VAR_BX .or. v == VAR_BY .or. v == VAR_BZ) cycle
      if (abs(realm%flx_f(v, ni, j, k, cblk) - flx_before(v, ni, j, k, cblk)) > tol) nb_ok = .false.
   enddo ; enddo ; enddo
   if (nb_ok) then
      write(*,'(A)') 'PASS: non-induction rows (D, currents) untouched at the seam'
   else
      write(*,'(A)') 'FAIL: a non-B row was modified at the seam face'
      test_passed = .false.
   endif
endblock

! --- Assertion 3: interior face fluxes (i /= ni) untouched. ---
block
   logical :: int_ok
   int_ok = .true.
   do k = 1_I4P, nk ; do j = 1_I4P, nj ; do i = 0_I4P, ni-1_I4P ; do v = 1_I4P, nv
      if (abs(realm%flx_f(v, i, j, k, cblk) - flx_before(v, i, j, k, cblk)) > tol) int_ok = .false.
   enddo ; enddo ; enddo ; enddo
   if (int_ok) then
      write(*,'(A)') 'PASS: interior face fluxes (i /= ni) untouched'
   else
      write(*,'(A)') 'FAIL: an interior face flux was modified'
      test_passed = .false.
   endif
endblock

! --- Assertion 4: the mirror-seam invariant — write-back is IDEMPOTENT. ---
! At a same-resolution mirror seam both sides reconstruct the identical face
! flux, so the register's canonical B-flux equals what is already in flx_f at
! the seam. Feeding the register exactly the seeded seam values must therefore
! leave flx_f unchanged. (This is the real no-op property; a register carrying
! a *different* symmetric value would, correctly, overwrite the seam — that is
! not the mirror-seam case.)
block
   logical :: noop_ok
   realm%flx_f = flx_before
   call feed_mirror_B_rows(reg, flx_before)
   call realm%apply_induction_flux_sharing_forest(stage=realm%rk%nrk, flux_register=reg)
   noop_ok = maxval(abs(realm%flx_f - flx_before)) <= tol
   if (noop_ok) then
      write(*,'(A)') 'PASS: mirror-seam register (canonical == existing seam flux) is idempotent'
   else
      write(*,'(A)') 'FAIL: idempotent write-back changed the flux (mirror-seam invariant broken)'
      test_passed = .false.
   endif
endblock

if (test_passed) then
   write(*,'(A)') 'TEST PASSED: prism induction flux share'
else
   write(*,'(A)') 'TEST FAILED: prism induction flux share'
endif

call reg%destroy
call mpih%finalize
if (.not. test_passed) error stop 1

contains
   subroutine seed_flx(f)
   !< Seed flx_f with a deterministic, cell-and-component-distinct pattern so any
   !< spurious write (wrong row, wrong cell) is detectable.
   real(R8P), intent(out) :: f(:,0:,:,:,:) !< Face flux (nv, 0:ni, nj, nk, nb).
   integer(I4P)           :: vv, ii, jj, kk, bb !< Counters.

   do bb = 1, size(f,5) ; do kk = 1, size(f,4) ; do jj = 1, size(f,3)
      do ii = 0, ni ; do vv = 1, size(f,1)
         f(vv,ii,jj,kk,bb) = real(vv,R8P) + 0.1_R8P*ii + 0.01_R8P*jj + 0.001_R8P*kk
      enddo ; enddo
   enddo ; enddo ; enddo
   endsubroutine seed_flx

   subroutine feed_asymmetric_B_rows(register, expected_canon)
   !< Push distinct coarse and fine B-row fluxes into the one registered face;
   !< return the expected canonical value = mean of the two.
   type(flux_register_object), intent(inout) :: register            !< Register (face 1 allocated).
   real(R8P),                  intent(out)    :: expected_canon(:,:) !< Expected canonical B flux (3, nface_cells).
   real(R8P), allocatable :: fc(:,:), ff(:,:)                        !< Coarse / fine flux skins (nv, nface_cells).
   integer(I4P)           :: cc, vb                                  !< Counters.

   allocate(fc(1:nv, 1:nface_cells), ff(1:nv, 1:nface_cells))
   fc = 0._R8P ; ff = 0._R8P
   do cc = 1, nface_cells
      do vb = 1, 3
         fc(bvars(vb), cc) = 2.0_R8P*vb + 0.5_R8P*cc        ! coarse-side B flux
         ff(bvars(vb), cc) =        vb - 0.25_R8P*cc        ! fine-side summed B flux (different)
         expected_canon(vb, cc) = 0.5_R8P*(fc(bvars(vb),cc) + ff(bvars(vb),cc))
      enddo
   enddo
   call register%accumulate_coarse_flux(face_index=1_I4P, stage=1_I4P, flux_face=fc)
   call register%accumulate_fine_flux  (face_index=1_I4P, stage=1_I4P, flux_face=ff)
   endsubroutine feed_asymmetric_B_rows

   subroutine feed_mirror_B_rows(register, flx_snap)
   !< Reset accumulators and push, on BOTH sides, exactly the B-row flux already
   !< present at the seam in `flx_snap` (the +x face, i = ni). This is the true
   !< same-resolution mirror seam: both sides reconstruct the identical face flux,
   !< so the canonical mean equals the existing seam value and the write-back is
   !< idempotent.
   type(flux_register_object), intent(inout) :: register      !< Register (face 1 allocated).
   real(R8P),                  intent(in)     :: flx_snap(:,0:,:,:,:) !< Snapshot of flx_f (nv, 0:ni, nj, nk, nb).
   real(R8P), allocatable :: fs(:,:)                            !< Mirror flux skin (nv, nface_cells).
   integer(I4P)           :: cc, jj, kk, vb                     !< Counters.

   call register%reset
   allocate(fs(1:nv, 1:nface_cells))
   fs = 0._R8P
   cc = 0_I4P
   do kk = 1, nk
      do jj = 1, nj
         cc = cc + 1_I4P
         do vb = 1, 3
            fs(bvars(vb), cc) = flx_snap(bvars(vb), ni, jj, kk, cblk)  ! the +x seam value already in flx_f
         enddo
      enddo
   enddo
   call register%accumulate_coarse_flux(face_index=1_I4P, stage=1_I4P, flux_face=fs)
   call register%accumulate_fine_flux  (face_index=1_I4P, stage=1_I4P, flux_face=fs)
   endsubroutine feed_mirror_B_rows
endprogram test_prism_induction_flux_share
