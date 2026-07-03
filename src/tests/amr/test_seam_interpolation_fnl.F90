!< Device unit test of the seam coarse->fine ghost-fill weight tables (issue #22 F2).

#include "fundal.H"

module test_seam_interpolation_fnl_kernels
!< Device kernels for the seam-interpolation device test: MODULE-LEVEL
!< subroutines (the `adam_fnl_field_kernels` production shape — NOT contained
!< kernels, whose private arrays are untrustworthy on nvfortran, see
!< CLAUDE-gpu / #22 F1-bis). Each kernel gathers its case data SCALAR-WISE
!< from shared device arrays into CONSTANT-BOUND private arrays, then calls
!< the single-source `!$acc routine seq` evaluators with whole-array actuals:
!< no array sections cross a device call boundary (nvfortran materializes
!< such sections into copy-in temporaries — the F1-bis race shape — and a
!< host-side `contiguous` repack of a device pointer segfaults outright).
!< This is EXACTLY the pattern the F3 fill kernels must use.
use adam_seam_interpolation_library, only : SEAM_W_TRICUBIC, SEAM_W_COMPATIBLE, SEAM_W_QUADRATIC, &
                                            seam_interpolate_tricubic, seam_interpolate_compatible, &
                                            seam_interpolate_quadratic, &
                                            seam_tricubic_centered_pos, seam_compatible_centered_pos, &
                                            seam_shift_anchor_pos, seam_meta_pack, seam_meta_unpack
use penf, only : I4P, R8P

implicit none
private
public :: dump_tables_dev, eval_tricubic_dev, eval_compatible_dev, eval_quadratic_dev, &
          meta_roundtrip_dev, shift_clamp_dev

contains
   subroutine dump_tables_dev(w4_dev, w3_dev, wq_dev)
   !< Copy the module parameter weight tables into device arrays, element-wise
   !< from device code: pins the GA3 parameter-residency risk in isolation.
   real(R8P), intent(inout) :: w4_dev(1:,1:,1:) !< Device copy of SEAM_W_TRICUBIC.
   real(R8P), intent(inout) :: w3_dev(1:,1:,1:) !< Device copy of SEAM_W_COMPATIBLE.
   real(R8P), intent(inout) :: wq_dev(1:,1:)    !< Device copy of SEAM_W_QUADRATIC.
   integer(I4P) :: i, s, p

   !$acc parallel loop independent collapse(3) gang vector &
   !$acc& DEVICEVAR(w4_dev)
   do p=1, 4
      do s=1, 2
         do i=1, 4
            w4_dev(i,s,p) = SEAM_W_TRICUBIC(i,s,p)
         enddo
      enddo
   enddo
   !$acc parallel loop independent collapse(3) gang vector &
   !$acc& DEVICEVAR(w3_dev)
   do p=1, 3
      do s=1, 2
         do i=1, 3
            w3_dev(i,s,p) = SEAM_W_COMPATIBLE(i,s,p)
         enddo
      enddo
   enddo
   !$acc parallel loop independent collapse(2) gang vector &
   !$acc& DEVICEVAR(wq_dev)
   do s=1, 2
      do i=1, 3
         wq_dev(i,s) = SEAM_W_QUADRATIC(i,s)
      enddo
   enddo
   endsubroutine dump_tables_dev

   subroutine eval_tricubic_dev(ncase, fp_dev, sub_dev, pos_dev, got_dev)
   !< Evaluate the tricubic interpolant on device, one thread per case:
   !< scalar gather into constant-bound private arrays, whole-array call.
   integer(I4P), intent(in)    :: ncase               !< Number of cases.
   real(R8P),    intent(in)    :: fp_dev(1:,1:,1:,1:) !< Footprints (4,4,4,ncase).
   integer(I4P), intent(in)    :: sub_dev(1:,1:)      !< Octant sub-positions (3,ncase).
   integer(I4P), intent(in)    :: pos_dev(1:,1:)      !< Anchor positions (3,ncase).
   real(R8P),    intent(inout) :: got_dev(1:)         !< Interpolated values (ncase).
   real(R8P)    :: fp(1:4,1:4,1:4) !< Private footprint gather buffer.
   integer(I4P) :: sub(1:3)        !< Private octant sub-position.
   integer(I4P) :: pos(1:3)        !< Private anchor position.
   integer(I4P) :: ic, i, j, k, d

   !$acc parallel loop independent gang vector              &
   !$acc& DEVICEVAR(fp_dev, sub_dev, pos_dev, got_dev)      &
   !$acc& firstprivate(ncase) private(fp, sub, pos)
   do ic=1, ncase
      !$acc loop seq
      do d=1, 3
         sub(d) = sub_dev(d,ic)
         pos(d) = pos_dev(d,ic)
      enddo
      !$acc loop seq
      do k=1, 4
         do j=1, 4
            do i=1, 4
               fp(i,j,k) = fp_dev(i,j,k,ic)
            enddo
         enddo
      enddo
      got_dev(ic) = seam_interpolate_tricubic(footprint=fp, sub=sub, anchor_pos=pos)
   enddo
   endsubroutine eval_tricubic_dev

   subroutine eval_compatible_dev(ncase, fp_dev, sub_dev, pos_dev, got_dev)
   !< Evaluate the restriction-compatible interpolant on device, one thread per case.
   integer(I4P), intent(in)    :: ncase               !< Number of cases.
   real(R8P),    intent(in)    :: fp_dev(1:,1:,1:,1:) !< Footprints (3,3,3,ncase).
   integer(I4P), intent(in)    :: sub_dev(1:,1:)      !< Octant sub-positions (3,ncase).
   integer(I4P), intent(in)    :: pos_dev(1:,1:)      !< Anchor positions (3,ncase).
   real(R8P),    intent(inout) :: got_dev(1:)         !< Interpolated values (ncase).
   real(R8P)    :: fp(1:3,1:3,1:3) !< Private footprint gather buffer.
   integer(I4P) :: sub(1:3)        !< Private octant sub-position.
   integer(I4P) :: pos(1:3)        !< Private anchor position.
   integer(I4P) :: ic, i, j, k, d

   !$acc parallel loop independent gang vector              &
   !$acc& DEVICEVAR(fp_dev, sub_dev, pos_dev, got_dev)      &
   !$acc& firstprivate(ncase) private(fp, sub, pos)
   do ic=1, ncase
      !$acc loop seq
      do d=1, 3
         sub(d) = sub_dev(d,ic)
         pos(d) = pos_dev(d,ic)
      enddo
      !$acc loop seq
      do k=1, 3
         do j=1, 3
            do i=1, 3
               fp(i,j,k) = fp_dev(i,j,k,ic)
            enddo
         enddo
      enddo
      got_dev(ic) = seam_interpolate_compatible(footprint=fp, sub=sub, anchor_pos=pos)
   enddo
   endsubroutine eval_compatible_dev

   subroutine eval_quadratic_dev(ncase, fp_dev, sub_dev, got_dev)
   !< Evaluate the centered quadratic interpolant on device, one thread per case.
   integer(I4P), intent(in)    :: ncase               !< Number of cases.
   real(R8P),    intent(in)    :: fp_dev(1:,1:,1:,1:) !< Footprints (3,3,3,ncase).
   integer(I4P), intent(in)    :: sub_dev(1:,1:)      !< Octant sub-positions (3,ncase).
   real(R8P),    intent(inout) :: got_dev(1:)         !< Interpolated values (ncase).
   real(R8P)    :: fp(1:3,1:3,1:3) !< Private footprint gather buffer.
   integer(I4P) :: sub(1:3)        !< Private octant sub-position.
   integer(I4P) :: ic, i, j, k, d

   !$acc parallel loop independent gang vector &
   !$acc& DEVICEVAR(fp_dev, sub_dev, got_dev)  &
   !$acc& firstprivate(ncase) private(fp, sub)
   do ic=1, ncase
      !$acc loop seq
      do d=1, 3
         sub(d) = sub_dev(d,ic)
      enddo
      !$acc loop seq
      do k=1, 3
         do j=1, 3
            do i=1, 3
               fp(i,j,k) = fp_dev(i,j,k,ic)
            enddo
         enddo
      enddo
      got_dev(ic) = seam_interpolate_quadratic(footprint=fp, sub=sub)
   enddo
   endsubroutine eval_quadratic_dev

   subroutine meta_roundtrip_dev(ncase, sub_dev, p4_dev, p3_dev, meta_dev, osub_dev, op4_dev, op3_dev)
   !< Pack then unpack the per-ghost metadata on device (ishft/ibits in device
   !< code): each thread round-trips its own case through the packed integer.
   integer(I4P), intent(in)    :: ncase           !< Number of cases.
   integer(I4P), intent(in)    :: sub_dev(1:,1:)  !< Octant sub-positions in (3,ncase).
   integer(I4P), intent(in)    :: p4_dev(1:,1:)   !< Tricubic anchor positions in (3,ncase).
   integer(I4P), intent(in)    :: p3_dev(1:,1:)   !< Compatible anchor positions in (3,ncase).
   integer(I4P), intent(inout) :: meta_dev(1:)    !< Packed metadata out (ncase).
   integer(I4P), intent(inout) :: osub_dev(1:,1:) !< Unpacked sub-positions out (3,ncase).
   integer(I4P), intent(inout) :: op4_dev(1:,1:)  !< Unpacked tricubic positions out (3,ncase).
   integer(I4P), intent(inout) :: op3_dev(1:,1:)  !< Unpacked compatible positions out (3,ncase).
   integer(I4P) :: sub(1:3)  !< Private sub-position (in).
   integer(I4P) :: p4(1:3)   !< Private tricubic position (in).
   integer(I4P) :: p3(1:3)   !< Private compatible position (in).
   integer(I4P) :: usub(1:3) !< Private sub-position (unpacked).
   integer(I4P) :: up4(1:3)  !< Private tricubic position (unpacked).
   integer(I4P) :: up3(1:3)  !< Private compatible position (unpacked).
   integer(I4P) :: ic, d

   !$acc parallel loop independent gang vector                                     &
   !$acc& DEVICEVAR(sub_dev, p4_dev, p3_dev, meta_dev, osub_dev, op4_dev, op3_dev) &
   !$acc& firstprivate(ncase) private(sub, p4, p3, usub, up4, up3)
   do ic=1, ncase
      !$acc loop seq
      do d=1, 3
         sub(d) = sub_dev(d,ic)
         p4(d)  = p4_dev(d,ic)
         p3(d)  = p3_dev(d,ic)
      enddo
      meta_dev(ic) = seam_meta_pack(sub=sub, p4=p4, p3=p3)
      call seam_meta_unpack(meta=meta_dev(ic), sub=usub, p4=up4, p3=up3)
      !$acc loop seq
      do d=1, 3
         osub_dev(d,ic) = usub(d)
         op4_dev(d,ic)  = up4(d)
         op3_dev(d,ic)  = up3(d)
      enddo
   enddo
   endsubroutine meta_roundtrip_dev

   subroutine shift_clamp_dev(ncase, anchor_dev, ncell_dev, subv_dev, fpn_dev, pc_dev, p_dev)
   !< Compute the centered anchor position and the shift-inward clamp on
   !< device, one thread per case (scalar arguments only: no gather needed).
   integer(I4P), intent(in)    :: ncase          !< Number of cases.
   integer(I4P), intent(in)    :: anchor_dev(1:) !< Donor anchor cell indexes (ncase).
   integer(I4P), intent(in)    :: ncell_dev(1:)  !< Donor block interior cell counts (ncase).
   integer(I4P), intent(in)    :: subv_dev(1:)   !< Octant sub-positions (ncase).
   integer(I4P), intent(in)    :: fpn_dev(1:)    !< Footprint widths, 4 tricubic | 3 compatible (ncase).
   integer(I4P), intent(inout) :: pc_dev(1:)     !< Centered anchor positions out (ncase).
   integer(I4P), intent(inout) :: p_dev(1:)      !< Clamped anchor positions out (ncase).
   integer(I4P) :: ic

   !$acc parallel loop independent gang vector                                &
   !$acc& DEVICEVAR(anchor_dev, ncell_dev, subv_dev, fpn_dev, pc_dev, p_dev)  &
   !$acc& firstprivate(ncase)
   do ic=1, ncase
      if (fpn_dev(ic) == 4_I4P) then
         pc_dev(ic) = seam_tricubic_centered_pos(sub=subv_dev(ic))
      else
         pc_dev(ic) = seam_compatible_centered_pos(sub=subv_dev(ic))
      endif
      p_dev(ic) = seam_shift_anchor_pos(anchor=anchor_dev(ic), n_cells=ncell_dev(ic), &
                                        p_centered=pc_dev(ic), footprint_n=fpn_dev(ic))
   enddo
   endsubroutine shift_clamp_dev
endmodule test_seam_interpolation_fnl_kernels

program test_seam_interpolation_fnl
!< Device twin of `test_seam_interpolation` (issue #22 F2, GA3): runs the SAME
!< assertion matrix as the CPU test, but every `got` value is computed on the
!< GPU by the single-source `!$acc routine seq` procedures of
!< `adam_seam_interpolation_library`, called from OpenACC parallel regions.
!<
!< **What this front-loads (the GA3 risks, before the F3 fill kernels):**
!<   0. Parameter-table device residency: `SEAM_W_*` are module PARAMETERS
!<      read inside device code — block 0 copies them back element-wise and
!<      demands BIT-EXACT equality with the host tables (they are binary
!<      rationals: any deviation is a residency/materialization bug, not
!<      round-off).
!<   1.-8. The full CPU matrix on device: partition of unity, tricubic
!<      monomial exactness (q=4), compatible linears (q=2), exact R o P = I
!<      on random data, the incompatibility theorem, the quadratic safety set,
!<      the metadata pack/unpack round trip (`ishft`/`ibits` IN DEVICE CODE),
!<      and the shift-inward clamp.
!<
!< **Race discipline (CLAUDE-gpu, #22 F1-bis):** the kernels live in a MODULE
!< (see `test_seam_interpolation_fnl_kernels` above) and gather into
!< constant-bound private arrays — the F3 production pattern. The case data
!< are VARIED (monomials, random fields) and every case has an exact analytic
!< answer, so a privatization failure cannot hide behind value coincidence.
!<
!< Memory management follows the `test_fdv_fnl_*` precedent: FUNDAL
!< dev_alloc/dev_memcpy/dev_free with raw device pointers seen via DEVICEVAR.

use adam_seam_interpolation_library, only : seam_tricubic_centered_pos, seam_compatible_centered_pos, &
                                            seam_meta_pack, &
                                            SEAM_W_TRICUBIC, SEAM_W_COMPATIBLE, SEAM_W_QUADRATIC
use fundal
use penf, only : I4P, R8P
use test_seam_interpolation_fnl_kernels

implicit none

real(R8P), parameter :: tol_exact = 1.e-14_R8P !< Tolerance for identities exact up to accumulation round-off.
real(R8P), parameter :: eta(1:2) = [-0.25_R8P, +0.25_R8P] !< Octant sub-position coordinates.
logical              :: test_passed             !< Aggregate pass flag.
integer(I4P)         :: ierr                    !< Device allocation error flag.

test_passed = .true.

! initialize fundal device environment
myhos   = dev_get_host_num()
devtype = dev_get_device_type()
call dev_set_device_num(0)
mydev = dev_get_device_num()

table_residency: block
   real(R8P), pointer :: w4_dev(:,:,:), w3_dev(:,:,:), wq_dev(:,:)
   real(R8P)          :: w4(1:4,1:2,1:4), w3(1:3,1:2,1:3), wq(1:3,1:2)
   logical            :: ok

   call dev_alloc(fptr_dev=w4_dev, lbounds=[1,1,1], ubounds=[4,2,4], ierr=ierr)
   call dev_alloc(fptr_dev=w3_dev, lbounds=[1,1,1], ubounds=[3,2,3], ierr=ierr)
   call dev_alloc(fptr_dev=wq_dev, lbounds=[1,1],   ubounds=[3,2],   ierr=ierr)
   call dump_tables_dev(w4_dev=w4_dev, w3_dev=w3_dev, wq_dev=wq_dev)
   call dev_memcpy_from_device(dst=w4, src=w4_dev)
   call dev_memcpy_from_device(dst=w3, src=w3_dev)
   call dev_memcpy_from_device(dst=wq, src=wq_dev)
   call dev_free(w4_dev, mydev)
   call dev_free(w3_dev, mydev)
   call dev_free(wq_dev, mydev)
   ! binary rationals: device copies must be BIT-EXACT, no tolerance
   ok = all(w4 == SEAM_W_TRICUBIC) .and. all(w3 == SEAM_W_COMPATIBLE) .and. all(wq == SEAM_W_QUADRATIC)
   if (ok) then
      write(*,'(A)') 'PASS: weight tables bit-exact on device (parameter residency, GA3 core risk)'
   else
      write(*,'(A)') 'FAIL: a device-read weight table differs from the host parameter'
      test_passed = .false.
   endif
endblock table_residency

partition_of_unity: block
   real(R8P),    allocatable :: fp4(:,:,:,:), fp3(:,:,:,:), got(:)
   integer(I4P), allocatable :: sub(:,:), pos(:,:)
   real(R8P),    pointer     :: fp_dev(:,:,:,:), got_dev(:)
   integer(I4P), pointer     :: sub_dev(:,:), pos_dev(:,:)
   integer(I4P)              :: sx, sy, sz, px, py, pz, ic, nc
   logical                   :: ok

   ok = .true.

   ! tricubic: all octants x all anchor states, constant footprint
   nc = 2*2*2 * 4*4*4
   allocate(fp4(1:4,1:4,1:4,1:nc), sub(1:3,1:nc), pos(1:3,1:nc), got(1:nc))
   ic = 0
   do sz=1,2 ; do sy=1,2 ; do sx=1,2
      do pz=1,4 ; do py=1,4 ; do px=1,4
         ic = ic + 1
         fp4(:,:,:,ic) = 1._R8P
         sub(:,ic) = [sx,sy,sz]
         pos(:,ic) = [px,py,pz]
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   call dev_alloc(fptr_dev=fp_dev,  lbounds=[1,1,1,1], ubounds=[4,4,4,nc], ierr=ierr)
   call dev_alloc(fptr_dev=sub_dev, lbounds=[1,1],     ubounds=[3,nc],     ierr=ierr)
   call dev_alloc(fptr_dev=pos_dev, lbounds=[1,1],     ubounds=[3,nc],     ierr=ierr)
   call dev_alloc(fptr_dev=got_dev, lbounds=[1],       ubounds=[nc],       ierr=ierr)
   call dev_memcpy_to_device(dst=fp_dev,  src=fp4)
   call dev_memcpy_to_device(dst=sub_dev, src=sub)
   call dev_memcpy_to_device(dst=pos_dev, src=pos)
   call eval_tricubic_dev(ncase=nc, fp_dev=fp_dev, sub_dev=sub_dev, pos_dev=pos_dev, got_dev=got_dev)
   call dev_memcpy_from_device(dst=got, src=got_dev)
   if (any(abs(got - 1._R8P) > tol_exact)) ok = .false.
   call dev_free(fp_dev, mydev) ; call dev_free(got_dev, mydev)
   deallocate(fp4, got)

   ! compatible: all octants x all anchor states (reuse sub/pos device arrays, first 8*27 slots)
   nc = 2*2*2 * 3*3*3
   allocate(fp3(1:3,1:3,1:3,1:nc), got(1:nc))
   ic = 0
   do sz=1,2 ; do sy=1,2 ; do sx=1,2
      do pz=1,3 ; do py=1,3 ; do px=1,3
         ic = ic + 1
         fp3(:,:,:,ic) = 1._R8P
         sub(:,ic) = [sx,sy,sz]
         pos(:,ic) = [px,py,pz]
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   call dev_alloc(fptr_dev=fp_dev,  lbounds=[1,1,1,1], ubounds=[3,3,3,nc], ierr=ierr)
   call dev_alloc(fptr_dev=got_dev, lbounds=[1],       ubounds=[nc],       ierr=ierr)
   call dev_memcpy_to_device(dst=fp_dev,  src=fp3)
   call dev_memcpy_to_device(dst=sub_dev, src=sub)
   call dev_memcpy_to_device(dst=pos_dev, src=pos)
   call eval_compatible_dev(ncase=nc, fp_dev=fp_dev, sub_dev=sub_dev, pos_dev=pos_dev, got_dev=got_dev)
   call dev_memcpy_from_device(dst=got, src=got_dev)
   if (any(abs(got(1:nc) - 1._R8P) > tol_exact)) ok = .false.

   ! quadratic: all octants, centered (reuse fp_dev first 8 slots; the compatible
   ! enumeration has pos innermost, so sub(:,1:8) must be refilled with the octants)
   ic = 0
   do sz=1,2 ; do sy=1,2 ; do sx=1,2
      ic = ic + 1
      sub(:,ic) = [sx,sy,sz]
   enddo ; enddo ; enddo
   call dev_memcpy_to_device(dst=sub_dev, src=sub)
   call eval_quadratic_dev(ncase=8_I4P, fp_dev=fp_dev, sub_dev=sub_dev, got_dev=got_dev)
   call dev_memcpy_from_device(dst=got, src=got_dev)
   if (any(abs(got(1:8) - 1._R8P) > tol_exact)) ok = .false.
   call dev_free(fp_dev, mydev) ; call dev_free(got_dev, mydev)
   call dev_free(sub_dev, mydev) ; call dev_free(pos_dev, mydev)
   deallocate(fp3, got, sub, pos)

   if (ok) then
      write(*,'(A)') 'PASS: partition of unity (constants reproduced by every regime/state, on device)'
   else
      write(*,'(A)') 'FAIL: a weight set does not sum to 1 on device'
      test_passed = .false.
   endif
endblock partition_of_unity

tricubic_monomials: block
   real(R8P),    allocatable :: fp4(:,:,:,:), got(:), want(:)
   integer(I4P), allocatable :: sub(:,:), pos(:,:)
   real(R8P),    pointer     :: fp_dev(:,:,:,:), got_dev(:)
   integer(I4P), pointer     :: sub_dev(:,:), pos_dev(:,:)
   real(R8P)                 :: x, y, z
   integer(I4P)              :: a, b, c, sx, sy, sz, px, py, pz, i, j, k, ic, nc
   logical                   :: ok

   nc = 4*4*4 * 2*2*2 * 4*4*4
   allocate(fp4(1:4,1:4,1:4,1:nc), sub(1:3,1:nc), pos(1:3,1:nc), got(1:nc), want(1:nc))
   ic = 0
   do a=0,3 ; do b=0,3 ; do c=0,3
      do sz=1,2 ; do sy=1,2 ; do sx=1,2
         do pz=1,4 ; do py=1,4 ; do px=1,4
            ic = ic + 1
            do k=1,4 ; do j=1,4 ; do i=1,4
               x = real(i - px, R8P) ; y = real(j - py, R8P) ; z = real(k - pz, R8P)
               fp4(i,j,k,ic) = x**a * y**b * z**c
            enddo ; enddo ; enddo
            sub(:,ic)  = [sx,sy,sz]
            pos(:,ic)  = [px,py,pz]
            want(ic)   = eta(sx)**a * eta(sy)**b * eta(sz)**c
         enddo ; enddo ; enddo
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   call dev_alloc(fptr_dev=fp_dev,  lbounds=[1,1,1,1], ubounds=[4,4,4,nc], ierr=ierr)
   call dev_alloc(fptr_dev=sub_dev, lbounds=[1,1],     ubounds=[3,nc],     ierr=ierr)
   call dev_alloc(fptr_dev=pos_dev, lbounds=[1,1],     ubounds=[3,nc],     ierr=ierr)
   call dev_alloc(fptr_dev=got_dev, lbounds=[1],       ubounds=[nc],       ierr=ierr)
   call dev_memcpy_to_device(dst=fp_dev,  src=fp4)
   call dev_memcpy_to_device(dst=sub_dev, src=sub)
   call dev_memcpy_to_device(dst=pos_dev, src=pos)
   call eval_tricubic_dev(ncase=nc, fp_dev=fp_dev, sub_dev=sub_dev, pos_dev=pos_dev, got_dev=got_dev)
   call dev_memcpy_from_device(dst=got, src=got_dev)
   call dev_free(fp_dev, mydev)  ; call dev_free(got_dev, mydev)
   call dev_free(sub_dev, mydev) ; call dev_free(pos_dev, mydev)
   ok = .true.
   do ic=1, nc
      if (abs(got(ic) - want(ic)) > tol_exact*max(1._R8P, abs(want(ic)))) ok = .false.
   enddo
   deallocate(fp4, sub, pos, got, want)
   if (ok) then
      write(*,'(A)') 'PASS: tricubic exact on all per-direction monomials deg <= 3, all octants, all shift states (q=4, on device)'
   else
      write(*,'(A)') 'FAIL: tricubic missed a degree <= 3 monomial on device'
      test_passed = .false.
   endif
endblock tricubic_monomials

compatible_linears: block
   real(R8P),    allocatable :: fp3(:,:,:,:), got(:), want(:)
   integer(I4P), allocatable :: sub(:,:), pos(:,:)
   real(R8P),    pointer     :: fp_dev(:,:,:,:), got_dev(:)
   integer(I4P), pointer     :: sub_dev(:,:), pos_dev(:,:)
   real(R8P)                 :: x, y, z
   integer(I4P)              :: a, b, c, sx, sy, sz, px, py, pz, i, j, k, ic, nc
   logical                   :: ok

   nc = 2*2*2 * 2*2*2 * 3*3*3
   allocate(fp3(1:3,1:3,1:3,1:nc), sub(1:3,1:nc), pos(1:3,1:nc), got(1:nc), want(1:nc))
   ic = 0
   do a=0,1 ; do b=0,1 ; do c=0,1
      do sz=1,2 ; do sy=1,2 ; do sx=1,2
         do pz=1,3 ; do py=1,3 ; do px=1,3
            ic = ic + 1
            do k=1,3 ; do j=1,3 ; do i=1,3
               x = real(i - px, R8P) ; y = real(j - py, R8P) ; z = real(k - pz, R8P)
               fp3(i,j,k,ic) = x**a * y**b * z**c
            enddo ; enddo ; enddo
            sub(:,ic)  = [sx,sy,sz]
            pos(:,ic)  = [px,py,pz]
            want(ic)   = eta(sx)**a * eta(sy)**b * eta(sz)**c
         enddo ; enddo ; enddo
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   call dev_alloc(fptr_dev=fp_dev,  lbounds=[1,1,1,1], ubounds=[3,3,3,nc], ierr=ierr)
   call dev_alloc(fptr_dev=sub_dev, lbounds=[1,1],     ubounds=[3,nc],     ierr=ierr)
   call dev_alloc(fptr_dev=pos_dev, lbounds=[1,1],     ubounds=[3,nc],     ierr=ierr)
   call dev_alloc(fptr_dev=got_dev, lbounds=[1],       ubounds=[nc],       ierr=ierr)
   call dev_memcpy_to_device(dst=fp_dev,  src=fp3)
   call dev_memcpy_to_device(dst=sub_dev, src=sub)
   call dev_memcpy_to_device(dst=pos_dev, src=pos)
   call eval_compatible_dev(ncase=nc, fp_dev=fp_dev, sub_dev=sub_dev, pos_dev=pos_dev, got_dev=got_dev)
   call dev_memcpy_from_device(dst=got, src=got_dev)
   call dev_free(fp_dev, mydev)  ; call dev_free(got_dev, mydev)
   call dev_free(sub_dev, mydev) ; call dev_free(pos_dev, mydev)
   ok = .true.
   do ic=1, nc
      if (abs(got(ic) - want(ic)) > tol_exact*max(1._R8P, abs(want(ic)))) ok = .false.
   enddo
   deallocate(fp3, sub, pos, got, want)
   if (ok) then
      write(*,'(A)') 'PASS: compatible arm exact on linears, all octants, all shift states (q=2, on device)'
   else
      write(*,'(A)') 'FAIL: compatible arm missed a linear monomial on device'
      test_passed = .false.
   endif
endblock compatible_linears

compatible_round_trip: block
   real(R8P),    allocatable :: fp3(:,:,:,:), got(:), anchor(:)
   integer(I4P), allocatable :: sub(:,:), pos(:,:)
   real(R8P),    pointer     :: fp_dev(:,:,:,:), got_dev(:)
   integer(I4P), pointer     :: sub_dev(:,:), pos_dev(:,:)
   real(R8P)                 :: fp_state(1:3,1:3,1:3), mean8
   integer(I4P)              :: sx, sy, sz, px, py, pz, is, ic, nc
   logical                   :: ok

   nc = 3*3*3 * 2*2*2 ! 27 anchor states x 8 octants, one random field per state
   allocate(fp3(1:3,1:3,1:3,1:nc), sub(1:3,1:nc), pos(1:3,1:nc), got(1:nc), anchor(1:nc/8))
   call random_seed()
   ic = 0 ; is = 0
   do pz=1,3 ; do py=1,3 ; do px=1,3
      is = is + 1
      call random_number(fp_state)
      anchor(is) = fp_state(px,py,pz)
      do sz=1,2 ; do sy=1,2 ; do sx=1,2
         ic = ic + 1
         fp3(:,:,:,ic) = fp_state
         sub(:,ic) = [sx,sy,sz]
         pos(:,ic) = [px,py,pz]
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   call dev_alloc(fptr_dev=fp_dev,  lbounds=[1,1,1,1], ubounds=[3,3,3,nc], ierr=ierr)
   call dev_alloc(fptr_dev=sub_dev, lbounds=[1,1],     ubounds=[3,nc],     ierr=ierr)
   call dev_alloc(fptr_dev=pos_dev, lbounds=[1,1],     ubounds=[3,nc],     ierr=ierr)
   call dev_alloc(fptr_dev=got_dev, lbounds=[1],       ubounds=[nc],       ierr=ierr)
   call dev_memcpy_to_device(dst=fp_dev,  src=fp3)
   call dev_memcpy_to_device(dst=sub_dev, src=sub)
   call dev_memcpy_to_device(dst=pos_dev, src=pos)
   call eval_compatible_dev(ncase=nc, fp_dev=fp_dev, sub_dev=sub_dev, pos_dev=pos_dev, got_dev=got_dev)
   call dev_memcpy_from_device(dst=got, src=got_dev)
   call dev_free(fp_dev, mydev)  ; call dev_free(got_dev, mydev)
   call dev_free(sub_dev, mydev) ; call dev_free(pos_dev, mydev)
   ok = .true.
   do is=1, nc/8
      mean8 = sum(got(8*(is-1)+1:8*is)) * 0.125_R8P
      if (abs(mean8 - anchor(is)) > tol_exact) ok = .false.
   enddo
   deallocate(fp3, sub, pos, got, anchor)
   if (ok) then
      write(*,'(A)') 'PASS: compatible arm satisfies exact R o P = I on random data (8-octant mean = anchor, on device)'
   else
      write(*,'(A)') 'FAIL: compatible arm broke the restriction round trip on device'
      test_passed = .false.
   endif
endblock compatible_round_trip

incompatibility_theorem: block
   real(R8P),    allocatable :: fp4(:,:,:,:), got(:)
   integer(I4P), allocatable :: sub(:,:), pos(:,:)
   real(R8P),    pointer     :: fp_dev(:,:,:,:), got_dev(:)
   integer(I4P), pointer     :: sub_dev(:,:), pos_dev(:,:)
   real(R8P)                 :: mean8, x
   integer(I4P)              :: sx, sy, sz, px, py, pz, i, j, k, ic
   logical                   :: ok

   ! f = x^2, centered per-sub anchor states (the production kernel convention)
   allocate(fp4(1:4,1:4,1:4,1:8), sub(1:3,1:8), pos(1:3,1:8), got(1:8))
   ic = 0
   do sz=1,2 ; do sy=1,2 ; do sx=1,2
      ic = ic + 1
      px = seam_tricubic_centered_pos(sx)
      py = seam_tricubic_centered_pos(sy)
      pz = seam_tricubic_centered_pos(sz)
      do k=1,4 ; do j=1,4 ; do i=1,4
         x = real(i - px, R8P)
         fp4(i,j,k,ic) = x*x
      enddo ; enddo ; enddo
      sub(:,ic) = [sx,sy,sz]
      pos(:,ic) = [px,py,pz]
   enddo ; enddo ; enddo
   call dev_alloc(fptr_dev=fp_dev,  lbounds=[1,1,1,1], ubounds=[4,4,4,8], ierr=ierr)
   call dev_alloc(fptr_dev=sub_dev, lbounds=[1,1],     ubounds=[3,8],     ierr=ierr)
   call dev_alloc(fptr_dev=pos_dev, lbounds=[1,1],     ubounds=[3,8],     ierr=ierr)
   call dev_alloc(fptr_dev=got_dev, lbounds=[1],       ubounds=[8],       ierr=ierr)
   call dev_memcpy_to_device(dst=fp_dev,  src=fp4)
   call dev_memcpy_to_device(dst=sub_dev, src=sub)
   call dev_memcpy_to_device(dst=pos_dev, src=pos)
   call eval_tricubic_dev(ncase=8_I4P, fp_dev=fp_dev, sub_dev=sub_dev, pos_dev=pos_dev, got_dev=got_dev)
   call dev_memcpy_from_device(dst=got, src=got_dev)
   call dev_free(fp_dev, mydev)  ; call dev_free(got_dev, mydev)
   call dev_free(sub_dev, mydev) ; call dev_free(pos_dev, mydev)
   ! correct interpolation: each octant exactly eta^2 = 1/16; mean defect = H^2/16
   ok = all(abs(got - 0.0625_R8P) <= tol_exact)
   mean8 = sum(got) * 0.125_R8P
   if (abs(mean8 - 0.0625_R8P) > tol_exact) ok = .false.
   deallocate(fp4, sub, pos, got)
   if (ok) then
      write(*,'(A)') 'PASS: incompatibility theorem pinned on device (tricubic octants of x^2 exact = 1/16; mean defect = H^2/16)'
   else
      write(*,'(A)') 'FAIL: tricubic x^2 octant values or mean defect differ from the theorem on device'
      test_passed = .false.
   endif
endblock incompatibility_theorem

quadratic_safety: block
   real(R8P),    allocatable :: fp3(:,:,:,:), got(:), want(:)
   integer(I4P), allocatable :: sub(:,:)
   real(R8P),    pointer     :: fp_dev(:,:,:,:), got_dev(:)
   integer(I4P), pointer     :: sub_dev(:,:)
   real(R8P)                 :: x, y, z
   integer(I4P)              :: a, b, c, sx, sy, sz, i, j, k, ic, nc
   logical                   :: ok

   nc = 3*3*3 * 2*2*2
   allocate(fp3(1:3,1:3,1:3,1:nc), sub(1:3,1:nc), got(1:nc), want(1:nc))
   ic = 0
   do a=0,2 ; do b=0,2 ; do c=0,2
      do sz=1,2 ; do sy=1,2 ; do sx=1,2
         ic = ic + 1
         do k=1,3 ; do j=1,3 ; do i=1,3
            x = real(i - 2, R8P) ; y = real(j - 2, R8P) ; z = real(k - 2, R8P)
            fp3(i,j,k,ic) = x**a * y**b * z**c
         enddo ; enddo ; enddo
         sub(:,ic) = [sx,sy,sz]
         want(ic)  = eta(sx)**a * eta(sy)**b * eta(sz)**c
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   call dev_alloc(fptr_dev=fp_dev,  lbounds=[1,1,1,1], ubounds=[3,3,3,nc], ierr=ierr)
   call dev_alloc(fptr_dev=sub_dev, lbounds=[1,1],     ubounds=[3,nc],     ierr=ierr)
   call dev_alloc(fptr_dev=got_dev, lbounds=[1],       ubounds=[nc],       ierr=ierr)
   call dev_memcpy_to_device(dst=fp_dev,  src=fp3)
   call dev_memcpy_to_device(dst=sub_dev, src=sub)
   call eval_quadratic_dev(ncase=nc, fp_dev=fp_dev, sub_dev=sub_dev, got_dev=got_dev)
   call dev_memcpy_from_device(dst=got, src=got_dev)
   call dev_free(fp_dev, mydev) ; call dev_free(got_dev, mydev) ; call dev_free(sub_dev, mydev)
   ok = .true.
   do ic=1, nc
      if (abs(got(ic) - want(ic)) > tol_exact*max(1._R8P, abs(want(ic)))) ok = .false.
   enddo
   deallocate(fp3, sub, got, want)
   if (ok) then
      write(*,'(A)') 'PASS: quadratic safety set exact on per-direction monomials deg <= 2, centered (q=3, on device)'
   else
      write(*,'(A)') 'FAIL: quadratic safety set missed a degree <= 2 monomial on device'
      test_passed = .false.
   endif
endblock quadratic_safety

metadata_round_trip: block
   integer(I4P), allocatable :: sub(:,:), p4(:,:), p3(:,:), meta(:), osub(:,:), op4(:,:), op3(:,:)
   integer(I4P), pointer     :: sub_dev(:,:), p4_dev(:,:), p3_dev(:,:), meta_dev(:)
   integer(I4P), pointer     :: osub_dev(:,:), op4_dev(:,:), op3_dev(:,:)
   integer(I4P)              :: s1, s2, s3, a1, a2, a3, b1, b2, b3, ic, nc
   logical                   :: ok

   nc = 2*2*2 * 4*4*4 * 3*3*3
   allocate(sub(1:3,1:nc), p4(1:3,1:nc), p3(1:3,1:nc), meta(1:nc))
   allocate(osub(1:3,1:nc), op4(1:3,1:nc), op3(1:3,1:nc))
   ic = 0
   do s3=1,2 ; do s2=1,2 ; do s1=1,2
      do a3=1,4 ; do a2=1,4 ; do a1=1,4
         do b3=1,3 ; do b2=1,3 ; do b1=1,3
            ic = ic + 1
            sub(:,ic) = [s1,s2,s3] ; p4(:,ic) = [a1,a2,a3] ; p3(:,ic) = [b1,b2,b3]
         enddo ; enddo ; enddo
      enddo ; enddo ; enddo
   enddo ; enddo ; enddo
   call dev_alloc(fptr_dev=sub_dev,  lbounds=[1,1], ubounds=[3,nc], ierr=ierr)
   call dev_alloc(fptr_dev=p4_dev,   lbounds=[1,1], ubounds=[3,nc], ierr=ierr)
   call dev_alloc(fptr_dev=p3_dev,   lbounds=[1,1], ubounds=[3,nc], ierr=ierr)
   call dev_alloc(fptr_dev=meta_dev, lbounds=[1],   ubounds=[nc],   ierr=ierr)
   call dev_alloc(fptr_dev=osub_dev, lbounds=[1,1], ubounds=[3,nc], ierr=ierr)
   call dev_alloc(fptr_dev=op4_dev,  lbounds=[1,1], ubounds=[3,nc], ierr=ierr)
   call dev_alloc(fptr_dev=op3_dev,  lbounds=[1,1], ubounds=[3,nc], ierr=ierr)
   call dev_memcpy_to_device(dst=sub_dev, src=sub)
   call dev_memcpy_to_device(dst=p4_dev,  src=p4)
   call dev_memcpy_to_device(dst=p3_dev,  src=p3)
   call meta_roundtrip_dev(ncase=nc, sub_dev=sub_dev, p4_dev=p4_dev, p3_dev=p3_dev, meta_dev=meta_dev, &
                           osub_dev=osub_dev, op4_dev=op4_dev, op3_dev=op3_dev)
   call dev_memcpy_from_device(dst=meta, src=meta_dev)
   call dev_memcpy_from_device(dst=osub, src=osub_dev)
   call dev_memcpy_from_device(dst=op4,  src=op4_dev)
   call dev_memcpy_from_device(dst=op3,  src=op3_dev)
   call dev_free(sub_dev, mydev)  ; call dev_free(p4_dev, mydev)  ; call dev_free(p3_dev, mydev)
   call dev_free(meta_dev, mydev) ; call dev_free(osub_dev, mydev)
   call dev_free(op4_dev, mydev)  ; call dev_free(op3_dev, mydev)
   ok = .true.
   do ic=1, nc
      ! device pack must agree with host pack (same bit layout across compilation targets)
      if (meta(ic) /= seam_meta_pack(sub=sub(:,ic), p4=p4(:,ic), p3=p3(:,ic))) ok = .false.
      ! device unpack must invert device pack
      if (any(osub(:,ic) /= sub(:,ic)) .or. any(op4(:,ic) /= p4(:,ic)) .or. any(op3(:,ic) /= p3(:,ic))) ok = .false.
   enddo
   deallocate(sub, p4, p3, meta, osub, op4, op3)
   if (ok) then
      write(*,'(A)') 'PASS: metadata pack/unpack round trip on device, device pack = host pack (ishft/ibits in device code)'
   else
      write(*,'(A)') 'FAIL: device metadata pack/unpack corrupted a field or diverged from host packing'
      test_passed = .false.
   endif
endblock metadata_round_trip

shift_clamp: block
   integer(I4P), allocatable :: anchor(:), ncell(:), subv(:), fpn(:), pc(:), p(:)
   integer(I4P), pointer     :: anchor_dev(:), ncell_dev(:), subv_dev(:), fpn_dev(:), pc_dev(:), p_dev(:)
   integer(I4P)              :: a, n, sub_, f, ic, nc
   logical                   :: ok

   nc = (4 + 16) * 2 * 2 ! n in {4,16} x anchors x subs x {tricubic, compatible}
   allocate(anchor(1:nc), ncell(1:nc), subv(1:nc), fpn(1:nc), pc(1:nc), p(1:nc))
   ic = 0
   do n=4,16,12 ! smallest legal block and the E0 block size
      do a=1,n
         do sub_=1,2
            do f=3,4
               ic = ic + 1
               anchor(ic) = a ; ncell(ic) = n ; subv(ic) = sub_ ; fpn(ic) = f
            enddo
         enddo
      enddo
   enddo
   call dev_alloc(fptr_dev=anchor_dev, lbounds=[1], ubounds=[nc], ierr=ierr)
   call dev_alloc(fptr_dev=ncell_dev,  lbounds=[1], ubounds=[nc], ierr=ierr)
   call dev_alloc(fptr_dev=subv_dev,   lbounds=[1], ubounds=[nc], ierr=ierr)
   call dev_alloc(fptr_dev=fpn_dev,    lbounds=[1], ubounds=[nc], ierr=ierr)
   call dev_alloc(fptr_dev=pc_dev,     lbounds=[1], ubounds=[nc], ierr=ierr)
   call dev_alloc(fptr_dev=p_dev,      lbounds=[1], ubounds=[nc], ierr=ierr)
   call dev_memcpy_to_device(dst=anchor_dev, src=anchor)
   call dev_memcpy_to_device(dst=ncell_dev,  src=ncell)
   call dev_memcpy_to_device(dst=subv_dev,   src=subv)
   call dev_memcpy_to_device(dst=fpn_dev,    src=fpn)
   call shift_clamp_dev(ncase=nc, anchor_dev=anchor_dev, ncell_dev=ncell_dev, subv_dev=subv_dev, &
                        fpn_dev=fpn_dev, pc_dev=pc_dev, p_dev=p_dev)
   call dev_memcpy_from_device(dst=pc, src=pc_dev)
   call dev_memcpy_from_device(dst=p,  src=p_dev)
   call dev_free(anchor_dev, mydev) ; call dev_free(ncell_dev, mydev) ; call dev_free(subv_dev, mydev)
   call dev_free(fpn_dev, mydev)    ; call dev_free(pc_dev, mydev)    ; call dev_free(p_dev, mydev)
   ok = .true.
   do ic=1, nc
      ! device centered positions must match the host functions
      if (fpn(ic) == 4) then
         if (pc(ic) /= seam_tricubic_centered_pos(subv(ic))) ok = .false.
      else
         if (pc(ic) /= seam_compatible_centered_pos(subv(ic))) ok = .false.
      endif
      ! clamp invariants (same assertions as the CPU test)
      if (p(ic) < 1 .or. p(ic) > fpn(ic)) ok = .false.
      if (anchor(ic) + 1 - p(ic) < 1 .or. anchor(ic) + fpn(ic) - p(ic) > ncell(ic)) ok = .false.
      if (anchor(ic) + 1 - pc(ic) >= 1 .and. anchor(ic) + fpn(ic) - pc(ic) <= ncell(ic) .and. p(ic) /= pc(ic)) ok = .false.
   enddo
   deallocate(anchor, ncell, subv, fpn, pc, p)
   if (ok) then
      write(*,'(A)') 'PASS: shift-inward clamp on device keeps every footprint inside real cells and stays centered when possible'
   else
      write(*,'(A)') 'FAIL: device shift-inward clamp produced an out-of-block footprint or diverged from host'
      test_passed = .false.
   endif
endblock shift_clamp

if (test_passed) then
   write(*,'(A)') 'TEST PASSED: seam interpolation weight tables (FNL device backend)'
else
   write(*,'(A)') 'TEST FAILED: seam interpolation weight tables (FNL device backend)'
endif
if (.not. test_passed) error stop 1
endprogram test_seam_interpolation_fnl
