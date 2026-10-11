!< ADAM seam-interpolation library: coarse->fine ghost-fill weight tables (issue #21 N0).
module adam_seam_interpolation_library
!< Weight tables and pure evaluators for the 2:1 coarse->fine seam ghost fill
!< (issue #21, plan amendment v2, milestone N0).
!<
!< **Geometry.** Uniform Cartesian 2:1 refinement, cell-centered point values.
!< The 8 fine ghost cells overlying a coarse donor (anchor) cell sit at octant
!< sub-positions eta = -1/4 (sub=1) and eta = +1/4 (sub=2) of the coarse
!< spacing, per direction. Weights are 1D vectors composed in tensor product
!< (w = wx*wy*wz) by the evaluators; the fill kernel (N2) reads the same tables.
!<
!< **Anchor position ("shift state", G3/G4).** The anchor cell occupies
!< position `p` of the 1D footprint; node offsets relative to the anchor are
!< `([1..n] - p)`. Centered states: tricubic p=2 for sub=2, p=3 for sub=1;
!< compatible p=2 for both. Off-centered states implement the shift-inward
!< rule at tangential block boundaries (donors are real cells, always).
!<
!< **The two regimes (G1 as amended, #21 issuecomment-4869444520).**
!< Exact restriction-compatibility R o P = I (R = the existing 8-cell-mean
!< restriction, `one_or_eight=8`) is PROVABLY incompatible with reproduction
!< of quadratics for point values: any P exact on quadratics fills the octants
!< of f = x^2 with the exact values (all H^2/16), whose 8-mean is H^2/16 /= 0
!< = f(anchor). Hence:
!< - `SEAM_FILL_TRICUBIC` (default candidate): plain tensor tricubic, exact for
!<   per-direction degree <= 3, ghost accuracy q = 4, mean defect O(H^2).
!< - `SEAM_FILL_COMPATIBLE` (diagnostic arm): the unique symmetric 3-node pair
!<   with exact R o P = I for arbitrary data, exact for linears only (q = 2).
!< The quadratic set is the G3 safety valve (centered only; must never be
!< reached in production maps -- reaching it is a maps-build bug).
!<
!< All weights are integer/2^k rationals: EXACT in binary floating point.
!< Partition of unity therefore holds to 0 ulp in every table.
!<
!< **Unrefined directions (issue #46).** A quadtree refines x and y only (a binary
!< tree x only): along an unrefined direction the fine ghost and the coarse donor
!< share the cell, the sub-position is "coincident" (sub=3, eta=0) and the 1D
!< weight is the identity at the anchor node, in every table. The evaluators skip
!< the zero-weight footprint nodes, so the unread donors (beyond the cells of a
!< block one cell thick, nk = 1) never enter the sum; on refined directions no
!< tricubic weight is zero and the octree results are unchanged bit for bit.
!<
!< **Device availability (issue #22 F2).** Every public procedure carries
!< `!$acc routine seq`: single-source host/device, callable from OpenACC
!< parallel regions on the FNL backend (weight tables are parameters,
!< materialized by the compiler on device). Inert comments on non-OpenACC
!< builds. Pinned on device by `test_seam_interpolation_fnl`.
use penf, only : I4P, R8P

implicit none
private
public :: SEAM_FILL_INJECTION, SEAM_FILL_COMPATIBLE, SEAM_FILL_TRICUBIC
public :: SEAM_W_TRICUBIC, SEAM_W_COMPATIBLE, SEAM_W_QUADRATIC
public :: seam_tricubic_centered_pos, seam_compatible_centered_pos
public :: seam_interpolate_tricubic, seam_interpolate_compatible, seam_interpolate_quadratic
public :: seam_shift_anchor_pos, seam_meta_pack, seam_meta_unpack
public :: SEAM_SUB_COINCIDENT
public :: SEAM_RESTRICT_MEAN, SEAM_RESTRICT_LAGRANGE6, SEAM_W_LAGRANGE6, SEAM_RESTRICT_IDENTITY
public :: seam_restriction_code, seam_restriction_meta_pack, seam_restriction_meta_unpack

integer(I4P), parameter :: SEAM_FILL_INJECTION  = 0_I4P !< Regime: 0th-order injection (legacy behavior).
integer(I4P), parameter :: SEAM_FILL_COMPATIBLE = 1_I4P !< Regime: restriction-compatible fill, q=2 (diagnostic arm).
integer(I4P), parameter :: SEAM_FILL_TRICUBIC   = 2_I4P !< Regime: plain tensor tricubic, q=4 (default candidate).
integer(I4P), parameter :: SEAM_SUB_COINCIDENT  = 3_I4P !< Sub-position of an unrefined direction (eta = 0, issue #46).

real(R8P), parameter :: SEAM_W_TRICUBIC(1:4,1:3,1:4) = reshape([ &
   ! p=1 (anchor at node 1; offsets 0..3)
   195._R8P, -117._R8P,   65._R8P,  -15._R8P,  & ! sub=1, eta=-1/4
    77._R8P,   77._R8P,  -33._R8P,    7._R8P,  & ! sub=2, eta=+1/4
   128._R8P,    0._R8P,    0._R8P,    0._R8P,  & ! sub=3, eta=0 (unrefined direction)
   ! p=2 (offsets -1..2; centered for sub=2)
    15._R8P,  135._R8P,  -27._R8P,    5._R8P,  & ! sub=1
    -7._R8P,  105._R8P,   35._R8P,   -5._R8P,  & ! sub=2
     0._R8P,  128._R8P,    0._R8P,    0._R8P,  & ! sub=3
   ! p=3 (offsets -2..1; centered for sub=1)
    -5._R8P,   35._R8P,  105._R8P,   -7._R8P,  & ! sub=1
     5._R8P,  -27._R8P,  135._R8P,   15._R8P,  & ! sub=2
     0._R8P,    0._R8P,  128._R8P,    0._R8P,  & ! sub=3
   ! p=4 (anchor at node 4; offsets -3..0)
     7._R8P,  -33._R8P,   77._R8P,   77._R8P,  & ! sub=1
   -15._R8P,   65._R8P, -117._R8P,  195._R8P,  & ! sub=2
     0._R8P,    0._R8P,    0._R8P,  128._R8P], & ! sub=3
   [4,3,4]) / 128._R8P
   !< 1D cubic Lagrange weights at eta = -+1/4 and the identity at eta = 0: (node 1:4, sub 1:3, anchor position p 1:4).
   !< Exact for per-direction polynomials of degree <= 3 in every state.

real(R8P), parameter :: SEAM_W_COMPATIBLE(1:3,1:3,1:3) = reshape([ &
   ! p=1 (offsets 0..2): linear interp/extrapolation pair, mean = e_1
   10._R8P,  -2._R8P,   0._R8P,  & ! sub=1
    6._R8P,   2._R8P,   0._R8P,  & ! sub=2
    8._R8P,   0._R8P,   0._R8P,  & ! sub=3 (unrefined direction)
   ! p=2 (offsets -1..1): the unique symmetric compatible pair, mean = e_2
    1._R8P,   8._R8P,  -1._R8P,  & ! sub=1
   -1._R8P,   8._R8P,   1._R8P,  & ! sub=2
    0._R8P,   8._R8P,   0._R8P,  & ! sub=3
   ! p=3 (offsets -2..0): mirror of p=1, mean = e_3
    0._R8P,   2._R8P,   6._R8P,  & ! sub=1
    0._R8P,  -2._R8P,  10._R8P,  & ! sub=2
    0._R8P,   0._R8P,   8._R8P], & ! sub=3
   [3,3,3]) / 8._R8P
   !< 1D restriction-compatible weights (node 1:3, sub 1:3, anchor position p 1:3; sub=3 the identity):
   !< exact for linears; per direction (w(sub=1)+w(sub=2))/2 = e_p exactly, so the
   !< 8-octant mean of the tensor interpolant equals the anchor value for ARBITRARY
   !< data (exact R o P = I). Maximal-order compatible fill (q = 2, see theorem above).

real(R8P), parameter :: SEAM_W_QUADRATIC(1:3,1:3) = reshape([ &
    5._R8P,  30._R8P,  -3._R8P,  & ! sub=1, eta=-1/4
   -3._R8P,  30._R8P,   5._R8P,  & ! sub=2, eta=+1/4
    0._R8P,  32._R8P,   0._R8P], & ! sub=3, eta=0 (unrefined direction)
   [3,3]) / 32._R8P
   !< 1D centered quadratic weights (node 1:3 at offsets -1..1, sub 1:2): the G3
   !< safety valve, exact for per-direction degree <= 2 (q = 3). Centered only.

! Fine->coarse restriction of the coarse seam ghosts (issue #68 P1/C1): `[amr] seam_restriction`.
integer(I4P), parameter :: SEAM_RESTRICT_MEAN      = 0_I4P !< Regime: the mean of the r(1) r(2) r(3) fine cells, q = 2.
integer(I4P), parameter :: SEAM_RESTRICT_LAGRANGE6 = 1_I4P !< Regime: tensor 6-point Lagrange at the coarse centre, q = 6.
integer(I4P), parameter :: SEAM_RESTRICT_IDENTITY  = 5_I4P !< Per-direction code of an unrefined direction (issue #46).

real(R8P), parameter :: SEAM_W_LAGRANGE6(1:6,0:5) = reshape([ &
    63._R8P,  315._R8P, -210._R8P,  126._R8P,  -45._R8P,    7._R8P,  & ! code 0: the coarse centre between nodes 1, 2
    -7._R8P,  105._R8P,  210._R8P,  -70._R8P,   21._R8P,   -3._R8P,  & ! code 1: between nodes 2, 3
     3._R8P,  -25._R8P,  150._R8P,  150._R8P,  -25._R8P,    3._R8P,  & ! code 2: between nodes 3, 4 (centred)
    -3._R8P,   21._R8P,  -70._R8P,  210._R8P,  105._R8P,   -7._R8P,  & ! code 3: between nodes 4, 5
     7._R8P,  -45._R8P,  126._R8P, -210._R8P,  315._R8P,   63._R8P,  & ! code 4: between nodes 5, 6
   256._R8P,    0._R8P,    0._R8P,    0._R8P,    0._R8P,    0._R8P], & ! code 5: identity (unrefined direction)
   [6,6]) / 256._R8P
   !< 1D weights of the 6-point Lagrange restriction to the coarse centre, which sits midway between the two fine cells
   !< it covers (`first`, `first + 1`): the 6 nodes are fine cells `first - code .. first - code + 5`, shifted inside
   !< the donor block at its edge (code /= 2, the shift-inward footprint of #68 B3), exact for per-direction degree
   !< <= 5. Integer/2^8 rationals (computed with exact fractions): exact in binary floating point.

contains
   pure function seam_tricubic_centered_pos(sub) result(p)
   !< Centered anchor position for the tricubic footprint: nodes -1..2 for
   !< eta=+1/4 (p=2), nodes -2..1 for eta=-1/4 (p=3).
   integer(I4P), intent(in) :: sub !< Octant sub-position (1 => eta=-1/4, 2 => eta=+1/4, 3 => eta=0, unrefined).
   integer(I4P)             :: p   !< Anchor position in the 4-node footprint.
   !$acc routine seq
   !$omp declare target

   p = 4_I4P - min(sub, 2_I4P)
   endfunction seam_tricubic_centered_pos

   pure function seam_compatible_centered_pos(sub) result(p)
   !< Centered anchor position for the compatible footprint (p=2 for both subs).
   integer(I4P), intent(in) :: sub !< Octant sub-position (unused: symmetric footprint).
   integer(I4P)             :: p   !< Anchor position in the 3-node footprint.
   !$acc routine seq
   !$omp declare target

   p = 2_I4P + 0_I4P*sub
   endfunction seam_compatible_centered_pos

   pure function seam_interpolate_tricubic(footprint, sub, anchor_pos) result(value_)
   !< Tensor tricubic interpolation on a 4x4x4 coarse footprint. `footprint(i,j,k)`
   !< holds the donor values at node offsets `(i - anchor_pos(1), j - anchor_pos(2),
   !< k - anchor_pos(3))` in coarse-cell units; the target is the octant center
   !< `(eta(sub(1)), eta(sub(2)), eta(sub(3)))`, eta = (-1)^sub / 4.
   real(R8P),    intent(in) :: footprint(1:4,1:4,1:4) !< Coarse donor values.
   integer(I4P), intent(in) :: sub(1:3)               !< Octant sub-position per direction (1|2|3).
   integer(I4P), intent(in) :: anchor_pos(1:3)        !< Anchor position per direction (1..4).
   real(R8P)                :: value_                 !< Interpolated fine-ghost value.
   real(R8P)                :: wyz                    !< Tangential weight product.
   integer(I4P)             :: i, j, k                !< Footprint counters.
   !$acc routine seq
   !$omp declare target

   value_ = 0._R8P
   do k = 1, 4
      do j = 1, 4
         wyz = SEAM_W_TRICUBIC(j,sub(2),anchor_pos(2)) * SEAM_W_TRICUBIC(k,sub(3),anchor_pos(3))
         if (wyz == 0._R8P) cycle ! an unrefined direction: the donor is not read
         do i = 1, 4
            if (SEAM_W_TRICUBIC(i,sub(1),anchor_pos(1)) == 0._R8P) cycle
            value_ = value_ + SEAM_W_TRICUBIC(i,sub(1),anchor_pos(1)) * wyz * footprint(i,j,k)
         enddo
      enddo
   enddo
   endfunction seam_interpolate_tricubic

   pure function seam_interpolate_compatible(footprint, sub, anchor_pos) result(value_)
   !< Tensor restriction-compatible interpolation on a 3x3x3 coarse footprint;
   !< same index conventions as the tricubic evaluator.
   real(R8P),    intent(in) :: footprint(1:3,1:3,1:3) !< Coarse donor values.
   integer(I4P), intent(in) :: sub(1:3)               !< Octant sub-position per direction (1|2).
   integer(I4P), intent(in) :: anchor_pos(1:3)        !< Anchor position per direction (1..3).
   real(R8P)                :: value_                 !< Interpolated fine-ghost value.
   real(R8P)                :: wyz                    !< Tangential weight product.
   integer(I4P)             :: i, j, k                !< Footprint counters.
   !$acc routine seq
   !$omp declare target

   value_ = 0._R8P
   do k = 1, 3
      do j = 1, 3
         wyz = SEAM_W_COMPATIBLE(j,sub(2),anchor_pos(2)) * SEAM_W_COMPATIBLE(k,sub(3),anchor_pos(3))
         if (wyz == 0._R8P) cycle ! a zero weight or an unrefined direction: the donor is not read
         do i = 1, 3
            if (SEAM_W_COMPATIBLE(i,sub(1),anchor_pos(1)) == 0._R8P) cycle
            value_ = value_ + SEAM_W_COMPATIBLE(i,sub(1),anchor_pos(1)) * wyz * footprint(i,j,k)
         enddo
      enddo
   enddo
   endfunction seam_interpolate_compatible

   pure function seam_shift_anchor_pos(anchor, n_cells, p_centered, footprint_n) result(p)
   !< Anchor position after the G3 shift-inward rule: clamp the centered
   !< position so the footprint node offsets `([1..footprint_n] - p)` around
   !< `anchor` stay inside the donor block's REAL cells `[1..n_cells]`.
   !< Valid for `n_cells >= footprint_n` and `anchor` in `[1..n_cells]`
   !< (both guaranteed at maps-build; the maps pass asserts `n_cells >= 4`).
   integer(I4P), intent(in) :: anchor      !< Donor anchor cell index (interior, 1-based).
   integer(I4P), intent(in) :: n_cells     !< Donor block interior cell count in this direction.
   integer(I4P), intent(in) :: p_centered  !< Centered anchor position for the active sub-position.
   integer(I4P), intent(in) :: footprint_n !< Footprint width (4 tricubic, 3 compatible).
   integer(I4P)             :: p           !< Shifted anchor position (1..footprint_n).
   !$acc routine seq
   !$omp declare target

   p = min(max(p_centered, anchor + footprint_n - n_cells), anchor)
   endfunction seam_shift_anchor_pos

   pure function seam_meta_pack(sub, p4, p3) result(meta)
   !< Pack the per-fine-ghost interpolation metadata into one integer
   !< (map-row payload, G4): bits 0-2 = octant (sub-1 per direction),
   !< bits 3-8 = tricubic anchor positions - 1 (2 bits each),
   !< bits 9-14 = compatible anchor positions - 1 (2 bits each),
   !< bits 15-17 = unrefined (coincident, sub=3) flag per direction, issue #46
   !< (0 on every octree row: the octree payload is unchanged).
   !< Regime-independent: the kernel unpacks the field matching the active
   !< fill regime. NOTE: 0 is a VALID packed value — consumers must gate on
   !< the map's flag column, never on `meta /= 0`.
   integer(I4P), intent(in) :: sub(1:3) !< Octant sub-position per direction (1|2|3).
   integer(I4P), intent(in) :: p4(1:3)  !< Tricubic anchor position per direction (1..4).
   integer(I4P), intent(in) :: p3(1:3)  !< Compatible anchor position per direction (1..3).
   integer(I4P)             :: meta     !< Packed metadata.
   integer(I4P)             :: d        !< Direction counter.
   !$acc routine seq
   !$omp declare target

   meta = ishft(p4(1) - 1,  3) + ishft(p4(2) - 1,  5) + ishft(p4(3) - 1,  7) + &
          ishft(p3(1) - 1,  9) + ishft(p3(2) - 1, 11) + ishft(p3(3) - 1, 13)
   do d=1, 3
      if (sub(d) == SEAM_SUB_COINCIDENT) then
         meta = meta + ishft(1_I4P, 14 + d)
      else
         meta = meta + ishft(sub(d) - 1, d - 1)
      endif
   enddo
   endfunction seam_meta_pack

   pure subroutine seam_meta_unpack(meta, sub, p4, p3)
   !< Unpack the per-fine-ghost interpolation metadata (inverse of
   !< `seam_meta_pack`).
   integer(I4P), intent(in)  :: meta     !< Packed metadata.
   integer(I4P), intent(out) :: sub(1:3) !< Octant sub-position per direction (1|2|3).
   integer(I4P), intent(out) :: p4(1:3)  !< Tricubic anchor position per direction (1..4).
   integer(I4P), intent(out) :: p3(1:3)  !< Compatible anchor position per direction (1..3).
   integer(I4P)              :: d        !< Direction counter.
   !$acc routine seq
   !$omp declare target

   do d=1, 3
      sub(d) = 1_I4P + ibits(meta, d - 1,      1)
      if (ibits(meta, 14 + d, 1) == 1_I4P) sub(d) = SEAM_SUB_COINCIDENT
      p4(d)  = 1_I4P + ibits(meta, 3 + 2*(d-1), 2)
      p3(d)  = 1_I4P + ibits(meta, 9 + 2*(d-1), 2)
   enddo
   endsubroutine seam_meta_unpack

   pure function seam_interpolate_quadratic(footprint, sub) result(value_)
   !< Tensor centered quadratic interpolation on a 3x3x3 coarse footprint
   !< (safety valve: centered anchor only, node offsets -1..1 per direction).
   real(R8P),    intent(in) :: footprint(1:3,1:3,1:3) !< Coarse donor values.
   integer(I4P), intent(in) :: sub(1:3)               !< Octant sub-position per direction (1|2).
   real(R8P)                :: value_                 !< Interpolated fine-ghost value.
   real(R8P)                :: wyz                    !< Tangential weight product.
   integer(I4P)             :: i, j, k                !< Footprint counters.
   !$acc routine seq
   !$omp declare target

   value_ = 0._R8P
   do k = 1, 3
      do j = 1, 3
         wyz = SEAM_W_QUADRATIC(j,sub(2)) * SEAM_W_QUADRATIC(k,sub(3))
         do i = 1, 3
            value_ = value_ + SEAM_W_QUADRATIC(i,sub(1)) * wyz * footprint(i,j,k)
         enddo
      enddo
   enddo
   endfunction seam_interpolate_quadratic

   pure function seam_restriction_code(first, n_cells, ratio) result(code)
   !< Return the per-direction code of the 6-point restriction (issue #68 C1): the position of the fine pair
   !< (`first`, `first + 1`) in the 6-node footprint, the footprint kept inside the donor block's cells 1..n_cells
   !< (shift-inward); `SEAM_RESTRICT_IDENTITY` on an unrefined direction (ratio 1); -1 when the fine pair is not inside
   !< the donor block (the caller keeps the mean there).
   integer(I4P), intent(in) :: first   !< First of the two fine cells covered by the coarse cell (donor index).
   integer(I4P), intent(in) :: n_cells !< Donor block cells in this direction (at least 6 on a refined one).
   integer(I4P), intent(in) :: ratio   !< Refinement ratio in this direction (1 or 2).
   integer(I4P)             :: code    !< Restriction code (0..4, SEAM_RESTRICT_IDENTITY, or -1).
   !$acc routine seq
   !$omp declare target

   if (ratio == 1_I4P) then
      code = SEAM_RESTRICT_IDENTITY
   elseif (first < 1_I4P .or. first + 1_I4P > n_cells) then
      code = -1_I4P
   else
      code = first - min(max(first - 2_I4P, 1_I4P), n_cells - 5_I4P)
   endif
   endfunction seam_restriction_code

   pure function seam_restriction_meta_pack(code) result(meta)
   !< Pack the per-direction restriction codes (0..5) into one integer, 1 + code(1) + 6 code(2) + 36 code(3) (1..216):
   !< 0 is left free for the rows that keep the mean (the map's meta default), so (0, 0, 0) is not mistaken for it.
   integer(I4P), intent(in) :: code(1:3) !< Restriction code per direction.
   integer(I4P)             :: meta      !< Packed metadata.
   !$acc routine seq
   !$omp declare target

   meta = 1_I4P + code(1) + 6_I4P * code(2) + 36_I4P * code(3)
   endfunction seam_restriction_meta_pack

   pure subroutine seam_restriction_meta_unpack(meta, code)
   !< Unpack the per-direction restriction codes (inverse of `seam_restriction_meta_pack`).
   integer(I4P), intent(in)  :: meta      !< Packed metadata.
   integer(I4P), intent(out) :: code(1:3) !< Restriction code per direction.
   !$acc routine seq
   !$omp declare target

   code(1) = mod(meta - 1_I4P, 6_I4P)
   code(2) = mod((meta - 1_I4P) / 6_I4P, 6_I4P)
   code(3) = (meta - 1_I4P) / 36_I4P
   endsubroutine seam_restriction_meta_unpack
endmodule adam_seam_interpolation_library
