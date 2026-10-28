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
public :: SEAM_FILL_INJECTION, SEAM_FILL_COMPATIBLE, SEAM_FILL_TRICUBIC, SEAM_FILL_QUINTIC
public :: SEAM_DENSITY_GUARD, seam_mp_bound, seam_bracket_nodes, seam_bracket_range
public :: SEAM_W_QUINTIC, seam_quintic_centered_pos, seam_interpolate_quintic, seam_quintic_weight
public :: SEAM_FILL_WENO6, seam_fill_weno6_1d, seam_interpolate_weno6
public :: SEAM_FILL_QUINTIC_MP, seam_fill_quintic_mp_1d, seam_interpolate_quintic_mp
public :: SEAM_W_TRICUBIC, SEAM_W_COMPATIBLE, SEAM_W_QUADRATIC
public :: seam_tricubic_centered_pos, seam_compatible_centered_pos
public :: seam_interpolate_tricubic, seam_interpolate_compatible, seam_interpolate_quadratic
public :: seam_shift_anchor_pos, seam_meta_pack, seam_meta_unpack
public :: SEAM_SUB_COINCIDENT
public :: SEAM_RESTRICT_MEAN, SEAM_RESTRICT_LAGRANGE6, SEAM_RESTRICT_WENO6, SEAM_W_LAGRANGE6, SEAM_RESTRICT_IDENTITY
public :: seam_restriction_code, seam_restriction_meta_pack, seam_restriction_meta_unpack, seam_restriction_weight
public :: seam_restriction_weno6_1d, seam_restriction_weno6

integer(I4P), parameter :: SEAM_FILL_INJECTION  = 0_I4P !< Regime: 0th-order injection (legacy behavior).
integer(I4P), parameter :: SEAM_FILL_COMPATIBLE = 1_I4P !< Regime: restriction-compatible fill, q=2 (diagnostic arm).
integer(I4P), parameter :: SEAM_FILL_TRICUBIC   = 2_I4P !< Regime: plain tensor tricubic, q=4 (default candidate).
integer(I4P), parameter :: SEAM_FILL_QUINTIC    = 3_I4P !< Regime: tensor 6-point quintic, q=6 (issue #68 D1).
integer(I4P), parameter :: SEAM_FILL_WENO6      = 4_I4P !< Regime: WENO on the quintic footprint (issue #68 D2).
integer(I4P), parameter :: SEAM_FILL_QUINTIC_MP = 5_I4P !< Regime: quintic clamped to an MP bound (issue #68 D2).
real(R8P),    parameter :: SEAM_DENSITY_GUARD = 0.1_R8P !< The density guard of the coarse->fine fill (issue #101): a fine
                                                        !< ghost whose interpolated density is below this fraction of its
                                                        !< anchor's (the coarse cell holding it) takes the anchor's state.
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

real(R8P), parameter :: SEAM_W_QUINTIC(1:6,1:3,1:6) = reshape([ &
   ! p=1 (anchor at node 1; offsets 0..5)
   13923._R8P, -13923._R8P,  15470._R8P, -10710._R8P,   4095._R8P,   -663._R8P,  & ! sub=1, eta=-1/4
    4389._R8P,   7315._R8P,  -6270._R8P,   3990._R8P,  -1463._R8P,    231._R8P,  & ! sub=2, eta=+1/4
    8192._R8P,      0._R8P,      0._R8P,      0._R8P,      0._R8P,      0._R8P,  & ! sub=3, eta=0 (unrefined)
   ! p=2 (offsets -1..4)
     663._R8P,   9945._R8P,  -3978._R8P,   2210._R8P,   -765._R8P,    117._R8P,  & ! sub=1
    -231._R8P,   5775._R8P,   3850._R8P,  -1650._R8P,    525._R8P,    -77._R8P,  & ! sub=2
       0._R8P,   8192._R8P,      0._R8P,      0._R8P,      0._R8P,      0._R8P,  & ! sub=3
   ! p=3 (offsets -2..3; centred for sub=2)
    -117._R8P,   1365._R8P,   8190._R8P,  -1638._R8P,    455._R8P,    -63._R8P,  & ! sub=1
      77._R8P,   -693._R8P,   6930._R8P,   2310._R8P,   -495._R8P,     63._R8P,  & ! sub=2
       0._R8P,      0._R8P,   8192._R8P,      0._R8P,      0._R8P,      0._R8P,  & ! sub=3
   ! p=4 (offsets -3..2; centred for sub=1)
      63._R8P,   -495._R8P,   2310._R8P,   6930._R8P,   -693._R8P,     77._R8P,  & ! sub=1
     -63._R8P,    455._R8P,  -1638._R8P,   8190._R8P,   1365._R8P,   -117._R8P,  & ! sub=2
       0._R8P,      0._R8P,      0._R8P,   8192._R8P,      0._R8P,      0._R8P,  & ! sub=3
   ! p=5 (offsets -4..1)
     -77._R8P,    525._R8P,  -1650._R8P,   3850._R8P,   5775._R8P,   -231._R8P,  & ! sub=1
     117._R8P,   -765._R8P,   2210._R8P,  -3978._R8P,   9945._R8P,    663._R8P,  & ! sub=2
       0._R8P,      0._R8P,      0._R8P,      0._R8P,   8192._R8P,      0._R8P,  & ! sub=3
   ! p=6 (anchor at node 6; offsets -5..0)
     231._R8P,  -1463._R8P,   3990._R8P,  -6270._R8P,   7315._R8P,   4389._R8P,  & ! sub=1
    -663._R8P,   4095._R8P, -10710._R8P,  15470._R8P, -13923._R8P,  13923._R8P,  & ! sub=2
       0._R8P,      0._R8P,      0._R8P,      0._R8P,      0._R8P,   8192._R8P], & ! sub=3
   [6,3,6]) / 8192._R8P
   !< 1D quintic Lagrange weights at eta = -+1/4 and the identity at eta = 0 (node 1:6, sub 1:3, anchor position p 1:6),
   !< issue #68 D1. Exact for per-direction polynomials of degree <= 5 in every state; integer/2^13 rationals (computed
   !< with exact fractions), exact in binary floating point; no weight of a refined direction is zero.

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
integer(I4P), parameter :: SEAM_RESTRICT_WENO6     = 2_I4P !< Regime: WENO over the 4-point sub-stencils, mean if troubled.
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

! The WENO restriction (issue #68 C3): per direction a nonlinear convex combination of the cubic interpolants of the
! 4-point sub-stencils of the 6 nodes (sub-stencil s = nodes s..s+3, s = 1..3) that contain the coarse centre, i.e. the
! interval between nodes code+1 and code+2. Sub-stencils that would extrapolate are never used: across a jump they
! import the far side's state. The linear weights reproduce the Lagrange interpolant of the candidates' union, so the
! design order is 6 for code 2 (three candidates, the 6-point Lagrange of `SEAM_W_LAGRANGE6`), 5 for codes 1 and 3
! (two candidates, the 5-point Lagrange) and 4 for codes 0 and 4 (one candidate, at a block edge).
real(R8P), parameter :: SEAM_W_CUBIC_MID(1:4,1:3) = reshape([ &
    5._R8P,  15._R8P,  -5._R8P,   1._R8P,  & ! the midpoint of the sub-stencil's interval 1 (nodes 1, 2)
   -1._R8P,   9._R8P,   9._R8P,  -1._R8P,  & ! interval 2 (nodes 2, 3), centred
    1._R8P,  -5._R8P,  15._R8P,   5._R8P], & ! interval 3 (nodes 3, 4)
   [4,3]) / 16._R8P
   !< Cubic Lagrange weights at the midpoint of each interval of a 4-point sub-stencil (node 1:4, interval 1:3).

real(R8P), parameter :: SEAM_D_WENO6(1:3,0:4) = reshape([ &
   16._R8P,   0._R8P,   0._R8P,  & ! code 0: sub-stencil 1 only
   10._R8P,   6._R8P,   0._R8P,  & ! code 1: sub-stencils 1, 2 -> 5-point Lagrange on nodes 1..5
    3._R8P,  10._R8P,   3._R8P,  & ! code 2: sub-stencils 1, 2, 3 -> 6-point Lagrange (centred)
    0._R8P,   6._R8P,  10._R8P,  & ! code 3: sub-stencils 2, 3 -> 5-point Lagrange on nodes 2..6
    0._R8P,   0._R8P,  16._R8P], & ! code 4: sub-stencil 3 only
   [3,5]) / 16._R8P
   !< Linear weights of the sub-stencils (sub-stencil 1:3, code 0:4); zero marks a sub-stencil that is not a candidate.

real(R8P), parameter :: SEAM_WENO6_JUMP_RATIO = 2._R8P
   !< Jump detector of the WENO restriction: an interval whose undivided difference exceeds this ratio times its
   !< neighbours' is a jump. Smooth data never trips it once resolved (the neighbours average to the interval's
   !< difference up to O(h^3)); a candidate containing a jump is dropped, and with no candidate left the restriction
   !< falls back to the mean of the fine pair (the troubled-cell fallback).
real(R8P), parameter :: SEAM_WENO6_EPS = 1.e-40_R8P !< Guard of the WENO-Z weights against 0/0 on constant data.

! The WENO coarse->fine fill (issue #68 D2): the target, anchor node p + eta (eta = -+1/4), sits at t4/4 in node units,
! t4 = 4 p -+ 1 (3..25, target index j = (t4 - 1)/2 = 1..12); its candidates are the sub-stencils containing the
! interval m = clamp(floor(t4/4), 1, 5) (at the block edge, j = 1 or 12, the one candidate extrapolates a quarter cell,
! the tricubic's own one-sided cubic).
real(R8P), parameter :: SEAM_W_CUBIC_Q(1:4,1:8) = reshape([ &
   195._R8P, -117._R8P,   65._R8P,  -15._R8P,  & ! r = 1: a quarter node spacing before node 1 (block edge)
    77._R8P,   77._R8P,  -33._R8P,    7._R8P,  & ! r = 2: node 1 + 1/4
    15._R8P,  135._R8P,  -27._R8P,    5._R8P,  & ! r = 3: node 1 + 3/4
    -7._R8P,  105._R8P,   35._R8P,   -5._R8P,  & ! r = 4: node 2 + 1/4
    -5._R8P,   35._R8P,  105._R8P,   -7._R8P,  & ! r = 5: node 2 + 3/4
     5._R8P,  -27._R8P,  135._R8P,   15._R8P,  & ! r = 6: node 3 + 1/4
     7._R8P,  -33._R8P,   77._R8P,   77._R8P,  & ! r = 7: node 3 + 3/4
   -15._R8P,   65._R8P, -117._R8P,  195._R8P], & ! r = 8: node 4 + 1/4 (block edge)
   [4,8]) / 128._R8P
   !< Cubic Lagrange weights of a 4-point sub-stencil at its quarter positions (node 1:4, r = (t4 - 4 s + 3)/2 = 1:8).

real(R8P), parameter :: SEAM_D_FILL_WENO6(1:3,1:12) = reshape([ &
   320._R8P,   0._R8P,   0._R8P,  & ! j = 1, t = 3/4 (block edge): sub-stencil 1 only
   320._R8P,   0._R8P,   0._R8P,  & ! j = 2, t = 5/4
   320._R8P,   0._R8P,   0._R8P,  & ! j = 3, t = 7/4
   220._R8P, 100._R8P,   0._R8P,  & ! j = 4, t = 9/4: sub-stencils 1, 2 -> 5-point Lagrange
   180._R8P, 140._R8P,   0._R8P,  & ! j = 5, t = 11/4
    77._R8P, 198._R8P,  45._R8P,  & ! j = 6, t = 13/4: sub-stencils 1, 2, 3 -> 6-point Lagrange (centred, eta = +1/4)
    45._R8P, 198._R8P,  77._R8P,  & ! j = 7, t = 15/4: (centred, eta = -1/4)
     0._R8P, 140._R8P, 180._R8P,  & ! j = 8, t = 17/4
     0._R8P, 100._R8P, 220._R8P,  & ! j = 9, t = 19/4
     0._R8P,   0._R8P, 320._R8P,  & ! j = 10, t = 21/4
     0._R8P,   0._R8P, 320._R8P,  & ! j = 11, t = 23/4
     0._R8P,   0._R8P, 320._R8P], & ! j = 12, t = 25/4 (block edge)
   [3,12]) / 320._R8P
   !< Linear weights of the fill's sub-stencils (sub-stencil 1:3, target 1:12), all positive (computed with exact
   !< fractions); zero marks a sub-stencil that does not contain the target's interval.

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

   pure function seam_quintic_centered_pos(sub) result(p)
   !< Centered anchor position for the quintic footprint (issue #68 D1): nodes -2..3 for eta=+1/4 (p=3), nodes -3..2
   !< for eta=-1/4 (p=4); 3 on an unrefined direction (sub=3, only the anchor is read).
   integer(I4P), intent(in) :: sub !< Octant sub-position (1 => eta=-1/4, 2 => eta=+1/4, 3 => eta=0, unrefined).
   integer(I4P)             :: p   !< Anchor position in the 6-node footprint.
   !$acc routine seq
   !$omp declare target

   p = 5_I4P - min(sub, 2_I4P)
   endfunction seam_quintic_centered_pos

   pure function seam_quintic_weight(node, sub, p) result(w)
   !< Return the 1D weight `SEAM_W_QUINTIC(node, sub, p)` (issue #68 D1): the device kernels read the table through this
   !< routine, inside the module that owns it.
   integer(I4P), intent(in) :: node !< Footprint node (1..6).
   integer(I4P), intent(in) :: sub  !< Octant sub-position (1|2|3).
   integer(I4P), intent(in) :: p    !< Anchor position (1..6).
   real(R8P)                :: w    !< Weight.
   !$acc routine seq
   !$omp declare target

   w = SEAM_W_QUINTIC(node, sub, p)
   endfunction seam_quintic_weight

   pure function seam_interpolate_quintic(footprint, sub, anchor_pos) result(value_)
   !< Tensor quintic interpolation on a 6x6x6 coarse footprint (issue #68 D1), same index conventions as the tricubic
   !< evaluator. Summed direction by direction (x innermost, then y, then z), the order of the FNL twin; an unrefined
   !< direction (sub=3) reads the anchor node only, so the other footprint entries along it are not referenced.
   real(R8P),    intent(in) :: footprint(1:6,1:6,1:6) !< Coarse donor values.
   integer(I4P), intent(in) :: sub(1:3)               !< Octant sub-position per direction (1|2|3).
   integer(I4P), intent(in) :: anchor_pos(1:3)        !< Anchor position per direction (1..6).
   real(R8P)                :: value_                 !< Interpolated fine-ghost value.
   integer(I4P)             :: lo(1:3), hi(1:3)       !< Footprint nodes read per direction.
   real(R8P)                :: sum_x, sum_y           !< Partial sums.
   integer(I4P)             :: i, j, k                !< Footprint counters.
   !$acc routine seq
   !$omp declare target

   do i = 1, 3
      lo(i) = 1_I4P
      hi(i) = 6_I4P
      if (sub(i) == SEAM_SUB_COINCIDENT) then
         lo(i) = anchor_pos(i)
         hi(i) = anchor_pos(i)
      endif
   enddo
   value_ = 0._R8P
   do k = lo(3), hi(3)
      sum_y = 0._R8P
      do j = lo(2), hi(2)
         sum_x = 0._R8P
         do i = lo(1), hi(1)
            sum_x = sum_x + SEAM_W_QUINTIC(i,sub(1),anchor_pos(1)) * footprint(i,j,k)
         enddo
         sum_y = sum_y + SEAM_W_QUINTIC(j,sub(2),anchor_pos(2)) * sum_x
      enddo
      value_ = value_ + SEAM_W_QUINTIC(k,sub(3),anchor_pos(3)) * sum_y
   enddo
   endfunction seam_interpolate_quintic

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

   pure function seam_meta_pack(sub, p4, p3, p6) result(meta)
   !< Pack the per-fine-ghost interpolation metadata into one integer
   !< (map-row payload, G4): bits 0-2 = octant (sub-1 per direction),
   !< bits 3-8 = tricubic anchor positions - 1 (2 bits each),
   !< bits 9-14 = compatible anchor positions - 1 (2 bits each),
   !< bits 15-17 = unrefined (coincident, sub=3) flag per direction, issue #46
   !< (0 on every octree row: the octree payload is unchanged),
   !< bits 18-26 = quintic anchor positions - 1 (3 bits each), issue #68 D1.
   !< Regime-independent: the kernel unpacks the field matching the active
   !< fill regime. NOTE: 0 is a VALID packed value — consumers must gate on
   !< the map's flag column, never on `meta /= 0`.
   integer(I4P), intent(in) :: sub(1:3) !< Octant sub-position per direction (1|2|3).
   integer(I4P), intent(in) :: p4(1:3)  !< Tricubic anchor position per direction (1..4).
   integer(I4P), intent(in) :: p3(1:3)  !< Compatible anchor position per direction (1..3).
   integer(I4P), intent(in) :: p6(1:3)  !< Quintic anchor position per direction (1..6).
   integer(I4P)             :: meta     !< Packed metadata.
   integer(I4P)             :: d        !< Direction counter.
   !$acc routine seq
   !$omp declare target

   meta = ishft(p4(1) - 1,  3) + ishft(p4(2) - 1,  5) + ishft(p4(3) - 1,  7) + &
          ishft(p3(1) - 1,  9) + ishft(p3(2) - 1, 11) + ishft(p3(3) - 1, 13) + &
          ishft(p6(1) - 1, 18) + ishft(p6(2) - 1, 21) + ishft(p6(3) - 1, 24)
   do d=1, 3
      if (sub(d) == SEAM_SUB_COINCIDENT) then
         meta = meta + ishft(1_I4P, 14 + d)
      else
         meta = meta + ishft(sub(d) - 1, d - 1)
      endif
   enddo
   endfunction seam_meta_pack

   pure subroutine seam_meta_unpack(meta, sub, p4, p3, p6)
   !< Unpack the per-fine-ghost interpolation metadata (inverse of
   !< `seam_meta_pack`).
   integer(I4P), intent(in)  :: meta     !< Packed metadata.
   integer(I4P), intent(out) :: sub(1:3) !< Octant sub-position per direction (1|2|3).
   integer(I4P), intent(out) :: p4(1:3)  !< Tricubic anchor position per direction (1..4).
   integer(I4P), intent(out) :: p3(1:3)  !< Compatible anchor position per direction (1..3).
   integer(I4P), intent(out) :: p6(1:3)  !< Quintic anchor position per direction (1..6).
   integer(I4P)              :: d        !< Direction counter.
   !$acc routine seq
   !$omp declare target

   do d=1, 3
      sub(d) = 1_I4P + ibits(meta, d - 1,      1)
      if (ibits(meta, 14 + d, 1) == 1_I4P) sub(d) = SEAM_SUB_COINCIDENT
      p4(d)  = 1_I4P + ibits(meta, 3 + 2*(d-1), 2)
      p3(d)  = 1_I4P + ibits(meta, 9 + 2*(d-1), 2)
      p6(d)  = 1_I4P + ibits(meta, 18 + 3*(d-1), 3)
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

   pure function seam_restriction_meta_pack(code, regime) result(meta)
   !< Pack the per-direction restriction codes (0..5) and the regime into one integer,
   !< 1 + code(1) + 6 code(2) + 36 code(3) + 216 (regime - SEAM_RESTRICT_LAGRANGE6): 1..216 under the 6-point Lagrange,
   !< 217..432 under the WENO restriction (#68 C3), so the rows carry their regime to every evaluator (inter-realm
   !< packers included); 0 is left free for the rows that keep the mean (the map's meta default).
   integer(I4P), intent(in) :: code(1:3) !< Restriction code per direction.
   integer(I4P), intent(in) :: regime    !< Restriction regime (SEAM_RESTRICT_LAGRANGE6 or SEAM_RESTRICT_WENO6).
   integer(I4P)             :: meta      !< Packed metadata.
   !$acc routine seq
   !$omp declare target

   meta = 1_I4P + code(1) + 6_I4P * code(2) + 36_I4P * code(3) + 216_I4P * (regime - SEAM_RESTRICT_LAGRANGE6)
   endfunction seam_restriction_meta_pack

   pure subroutine seam_restriction_meta_unpack(meta, code, regime)
   !< Unpack the per-direction restriction codes and the regime (inverse of `seam_restriction_meta_pack`).
   integer(I4P), intent(in)  :: meta      !< Packed metadata (/= 0).
   integer(I4P), intent(out) :: code(1:3) !< Restriction code per direction.
   integer(I4P), intent(out) :: regime    !< Restriction regime.
   integer(I4P)              :: m         !< Codes part of the meta.
   !$acc routine seq
   !$omp declare target

   regime  = SEAM_RESTRICT_LAGRANGE6 + (meta - 1_I4P) / 216_I4P
   m       = mod(meta - 1_I4P, 216_I4P)
   code(1) = mod(m, 6_I4P)
   code(2) = mod(m / 6_I4P, 6_I4P)
   code(3) = m / 36_I4P
   endsubroutine seam_restriction_meta_unpack

   pure subroutine seam_weno6_jumps(u, jump)
   !< Flag the jumps of 6 values (issue #68 C3): interval i (nodes i, i+1) is a jump when its undivided difference
   !< exceeds `SEAM_WENO6_JUMP_RATIO` times the larger of its neighbours' (at the footprint ends, the inner neighbour
   !< and its linear extrapolation). Smooth data never trips it once resolved.
   real(R8P), intent(in)  :: u(1:6)    !< Values at the footprint nodes.
   logical,   intent(out) :: jump(1:5) !< Jump flag per interval.
   real(R8P)              :: d(1:5)    !< Undivided first differences.
   real(R8P)              :: ref       !< Neighbours' difference.
   integer(I4P)           :: i         !< Counter.
   !$acc routine seq
   !$omp declare target

   do i=1, 5
      d(i) = u(i+1) - u(i)
   enddo
   do i=1, 5
      if (i == 1) then
         ref = max(abs(d(2)), abs(2._R8P * d(2) - d(3)))
      elseif (i == 5) then
         ref = max(abs(d(4)), abs(2._R8P * d(4) - d(3)))
      else
         ref = max(abs(d(i-1)), abs(d(i+1)))
      endif
      jump(i) = abs(d(i)) > SEAM_WENO6_JUMP_RATIO * ref
   enddo
   endsubroutine seam_weno6_jumps

   pure function seam_weno6_core(u, dlin, wsub, fallback) result(value_)
   !< The WENO kernel shared by the restriction and the fill (issue #68 C3, D2): a convex combination of the cubic
   !< interpolants of the candidate 4-point sub-stencils (sub-stencil s = nodes s..s+3 of the 6) at one target.
   !< 1. Jumps: interval i (nodes i, i+1) is a jump when its undivided difference exceeds `SEAM_WENO6_JUMP_RATIO` times
   !<    the larger of its neighbours' (at the footprint ends, the inner neighbour and its linear extrapolation).
   !< 2. Candidates: the sub-stencils with a positive linear weight `dlin(s)` and free of jumps; none left: `fallback`
   !<    (the troubled-cell value, chosen by the caller).
   !< 3. WENO-Z weights (Borges et al. 2008, p = 2) on the Jiang-Shu indicators of each cubic over its central
   !<    interval, beta = c^2 + 13/12 a^2 + 781/720 t^2 (c, a, t its central first, mean second and third differences),
   !<    with tau the highest undivided difference of the candidates' union squared, O(h^8..h^12) on smooth data, so the
   !<    weights reach the linear ones fast enough to keep the design order at critical points.
   real(R8P), intent(in) :: u(1:6)        !< Values at the footprint nodes.
   real(R8P), intent(in) :: dlin(1:3)     !< Linear weight per sub-stencil (0: not a candidate).
   real(R8P), intent(in) :: wsub(1:4,1:3) !< Cubic weights of each candidate sub-stencil at the target.
   real(R8P), intent(in) :: fallback      !< Value returned when no candidate is left.
   real(R8P)             :: value_        !< WENO value.
   real(R8P)             :: w(1:6)        !< Work array of the highest difference.
   real(R8P)             :: tau           !< WENO-Z global indicator.
   real(R8P)             :: c, a, t       !< Central first, mean second and third differences of a sub-stencil.
   real(R8P)             :: alpha         !< Unnormalised weight.
   real(R8P)             :: asum          !< Sum of the unnormalised weights.
   logical               :: jump(1:5)     !< Jump flag per interval.
   integer(I4P)          :: i, m, s       !< Counters.
   integer(I4P)          :: lo, hi        !< First and last candidate.
   !$acc routine seq
   !$omp declare target

   call seam_weno6_jumps(u=u, jump=jump)
   lo = 0_I4P
   hi = 0_I4P
   do s=1, 3 ! the candidates are contiguous: a sub-stencil between two jump-free ones is jump-free
      if (dlin(s) > 0._R8P .and. .not.(jump(s) .or. jump(s+1) .or. jump(s+2))) then
         if (lo == 0_I4P) lo = s
         hi = s
      endif
   enddo
   if (lo == 0_I4P) then
      value_ = fallback
      return
   endif
   if (lo == hi) then
      value_ = wsub(1,lo) * u(lo) + wsub(2,lo) * u(lo+1) + wsub(3,lo) * u(lo+2) + wsub(4,lo) * u(lo+3)
      return
   endif
   do i=1, hi - lo + 4 ! nodes lo..hi+3
      w(i) = u(lo+i-1)
   enddo
   do m=1, hi - lo + 3
      do i=1, hi - lo + 4 - m
         w(i) = w(i+1) - w(i)
      enddo
   enddo
   tau = w(1) * w(1)
   value_ = 0._R8P
   asum = 0._R8P
   do s=lo, hi
      c = u(s+2) - u(s+1)
      a = 0.5_R8P * (u(s+3) - u(s+2) - u(s+1) + u(s))
      t = u(s+3) - 3._R8P * u(s+2) + 3._R8P * u(s+1) - u(s)
      alpha = dlin(s) * (1._R8P + (tau / (c * c + 13._R8P / 12._R8P * a * a + &
                                          781._R8P / 720._R8P * t * t + SEAM_WENO6_EPS))**2)
      value_ = value_ + alpha * (wsub(1,s) * u(s) + wsub(2,s) * u(s+1) + wsub(3,s) * u(s+2) + wsub(4,s) * u(s+3))
      asum = asum + alpha
   enddo
   value_ = value_ / asum
   endfunction seam_weno6_core

   pure function seam_restriction_weno6_1d(u, code) result(value_)
   !< Return the WENO restriction of 6 fine values `u` to the coarse centre between nodes code+1 and code+2 (#68 C3):
   !< `seam_weno6_core` on the sub-stencils of `SEAM_D_WENO6(:,code)` with their midpoint cubics, the troubled-cell
   !< fallback the mean of the fine pair. `SEAM_RESTRICT_IDENTITY` (an unrefined direction) returns u(1).
   real(R8P),    intent(in) :: u(1:6)        !< Fine values (footprint nodes).
   integer(I4P), intent(in) :: code          !< Restriction code (0..4, or SEAM_RESTRICT_IDENTITY).
   real(R8P)                :: value_        !< Restricted value.
   real(R8P)                :: wsub(1:4,1:3) !< Cubic weights of the candidates at the coarse centre.
   integer(I4P)             :: s             !< Sub-stencil counter.
   !$acc routine seq
   !$omp declare target

   if (code == SEAM_RESTRICT_IDENTITY) then
      value_ = u(1)
      return
   endif
   wsub = 0._R8P
   do s=1, 3
      if (SEAM_D_WENO6(s,code) > 0._R8P) wsub(:,s) = SEAM_W_CUBIC_MID(:,code-s+2)
   enddo
   value_ = seam_weno6_core(u=u, dlin=SEAM_D_WENO6(:,code), wsub=wsub, fallback=0.5_R8P * (u(code+1) + u(code+2)))
   endfunction seam_restriction_weno6_1d

   pure function seam_fill_weno6_1d(u, sub, p) result(value_)
   !< Return the WENO coarse->fine fill of 6 coarse values `u` at the fine ghost's centre, anchor node `p` + eta(sub)
   !< (issue #68 D2): `seam_weno6_core` on the sub-stencils that contain the target's interval (none extrapolates but
   !< at the block edge, where the one candidate is the tricubic's own one-sided cubic), linear weights
   !< `SEAM_D_FILL_WENO6` reproducing the Lagrange interpolant of their union (order 6 at the centred anchor
   !< positions, 5 one cell from the block edge, 4 at it); the troubled-cell fallback is the linear interpolation on
   !< the target's interval clamped to its two values. An unrefined direction (sub = 3) returns u(p).
   real(R8P),    intent(in) :: u(1:6)        !< Coarse values (footprint nodes).
   integer(I4P), intent(in) :: sub           !< Octant sub-position (1 => eta=-1/4, 2 => eta=+1/4, 3 => unrefined).
   integer(I4P), intent(in) :: p             !< Anchor position (1..6).
   real(R8P)                :: value_        !< Fine-ghost value.
   real(R8P)                :: wsub(1:4,1:3) !< Cubic weights of the candidates at the target.
   real(R8P)                :: fallback      !< Clamped linear interpolation on the target's interval.
   integer(I4P)             :: t4            !< Target position in quarter units of the node spacing (3..25).
   integer(I4P)             :: j             !< Target index (1..12).
   integer(I4P)             :: m             !< Interval holding the target (1..5, the block edge clamped).
   integer(I4P)             :: s             !< Sub-stencil counter.
   !$acc routine seq
   !$omp declare target

   if (sub == SEAM_SUB_COINCIDENT) then
      value_ = u(p)
      return
   endif
   t4 = 4_I4P * p + 2_I4P * sub - 3_I4P
   j = (t4 - 1_I4P) / 2_I4P
   m = min(max(t4 / 4_I4P, 1_I4P), 5_I4P)
   wsub = 0._R8P
   do s=1, 3
      if (SEAM_D_FILL_WENO6(s,j) > 0._R8P) wsub(:,s) = SEAM_W_CUBIC_Q(:,(t4 - 4_I4P * s + 3_I4P) / 2_I4P)
   enddo
   fallback = u(m) + 0.25_R8P * real(t4 - 4_I4P * m, R8P) * (u(m+1) - u(m))
   fallback = min(max(fallback, min(u(m), u(m+1))), max(u(m), u(m+1)))
   value_ = seam_weno6_core(u=u, dlin=SEAM_D_FILL_WENO6(:,j), wsub=wsub, fallback=fallback)
   endfunction seam_fill_weno6_1d

   pure function seam_fill_quintic_mp_1d(u, sub, p) result(value_)
   !< Return the monotonicity-preserving quintic coarse->fine fill of 6 coarse values `u` at the fine ghost's centre,
   !< anchor node `p` + eta(sub) (issue #68 D2): the quintic of `SEAM_W_QUINTIC` clamped to a Suresh-Huynh (1997) MP
   !< bound. The target lies at fraction th of the interval m (nodes m, m+1), linear value L = u(m) + th (u(m+1) - u(m));
   !< dM is the M4 minmod of the two second differences around the interval (shifted inward at the block edge), zero
   !< across a jump; the bound is [min(u(m), u(m+1), Q), max(...)] with Q = L + th (th - 1)/2 dM, widened by
   !< |th (th - 1)/2 dM|: on smooth data it contains the quintic (which differs from Q at O(h^3)), across a jump it
   !< collapses to the interval's values. At the block edge (the target a quarter cell outside the footprint, th < 0 or
   !< > 1) L is clamped to the interval and the bound widened by the extrapolation's linear and curvature terms only when
   !< the three edge-side intervals hold no jump (`seam_weno6_jumps`). An unrefined direction (sub = 3) returns u(p).
   real(R8P),    intent(in) :: u(1:6)     !< Coarse values (footprint nodes).
   integer(I4P), intent(in) :: sub        !< Octant sub-position (1 => eta=-1/4, 2 => eta=+1/4, 3 => unrefined).
   integer(I4P), intent(in) :: p          !< Anchor position (1..6).
   real(R8P)                :: value_     !< Fine-ghost value.
   logical                  :: jump(1:5)  !< Jump flag per interval.
   real(R8P)                :: th         !< Target fraction of its interval.
   real(R8P)                :: c          !< Curvature coefficient th (th - 1)/2.
   real(R8P)                :: dk, dk1    !< Second differences around the interval.
   real(R8P)                :: dm         !< M4 minmod curvature.
   real(R8P)                :: lin        !< Linear value on the interval.
   real(R8P)                :: q          !< Curvature-corrected linear value.
   real(R8P)                :: widen      !< Bound widening.
   real(R8P)                :: d2e, d2i   !< Edge: the second differences at the edge-side and inner nodes.
   real(R8P)                :: out        !< Edge: the extrapolation distance (in interval units).
   real(R8P)                :: lo, hi     !< MP bound.
   integer(I4P)             :: t4, m, k, n !< Target (quarter units), interval, curvature nodes, counter.
   !$acc routine seq
   !$omp declare target

   if (sub == SEAM_SUB_COINCIDENT) then
      value_ = u(p)
      return
   endif
   value_ = 0._R8P
   do n=1, 6
      value_ = value_ + SEAM_W_QUINTIC(n,sub,p) * u(n)
   enddo
   t4 = 4_I4P * p + 2_I4P * sub - 3_I4P
   m = min(max(t4 / 4_I4P, 1_I4P), 5_I4P)
   th = 0.25_R8P * real(t4 - 4_I4P * m, R8P)
   c = 0.5_R8P * th * (th - 1._R8P)
   k = min(max(m, 2_I4P), 4_I4P)
   dk  = u(k-1) - 2._R8P * u(k)   + u(k+1)
   dk1 = u(k)   - 2._R8P * u(k+1) + u(k+2)
   dm = seam_minmod4(4._R8P * dk - dk1, 4._R8P * dk1 - dk, dk, dk1)
   lin = u(m) + th * (u(m+1) - u(m))
   q = lin + c * dm
   widen = abs(c * dm)
   if (th < 0._R8P .or. th > 1._R8P) then ! the block edge
      out = max(-th, th - 1._R8P)
      lin = min(max(lin, min(u(m), u(m+1))), max(u(m), u(m+1)))
      q = lin
      call seam_weno6_jumps(u=u, jump=jump)
      widen = 0._R8P
      if (th < 0._R8P) then
         if (.not.(jump(1) .or. jump(2) .or. jump(3))) then
            d2e = dk ; d2i = dk1
            widen = out * abs(u(m+1) - u(m)) + &
                    (abs(c) + out) * max(abs(d2e), abs(seam_minmod4(d2e, d2i, d2e, d2i)))
         endif
      else
         if (.not.(jump(3) .or. jump(4) .or. jump(5))) then
            d2e = dk1 ; d2i = dk
            widen = out * abs(u(m+1) - u(m)) + &
                    (abs(c) + out) * max(abs(d2e), abs(seam_minmod4(d2e, d2i, d2e, d2i)))
         endif
      endif
   endif
   lo = min(u(m), u(m+1), q) - widen
   hi = max(u(m), u(m+1), q) + widen
   value_ = min(max(value_, lo), hi)
   endfunction seam_fill_quintic_mp_1d

   pure function seam_minmod4(a, b, c, d) result(mm)
   !< Return the minmod of four values: the one of least magnitude if all share a sign, zero otherwise.
   real(R8P), intent(in) :: a, b, c, d !< Values.
   real(R8P)             :: mm         !< Minmod.
   !$acc routine seq
   !$omp declare target

   if (a > 0._R8P .and. b > 0._R8P .and. c > 0._R8P .and. d > 0._R8P) then
      mm = min(a, b, c, d)
   elseif (a < 0._R8P .and. b < 0._R8P .and. c < 0._R8P .and. d < 0._R8P) then
      mm = max(a, b, c, d)
   else
      mm = 0._R8P
   endif
   endfunction seam_minmod4

   pure function seam_restriction_weno6(code, u) result(value_)
   !< Return the WENO restriction of the fine footprint `u` (#68 C3), direction by direction: x on every (j, k) line,
   !< then y, then z (the same order on every backend). An unrefined direction (code `SEAM_RESTRICT_IDENTITY`) reads
   !< index 1 only, the other entries of `u` along it are not referenced. A troubled line falls back to the mean of its
   !< fine pair, so a footprint troubled on every line returns the mean of the fine cells.
   integer(I4P), intent(in) :: code(1:3)     !< Restriction code per direction.
   real(R8P),    intent(in) :: u(1:6,1:6,1:6) !< Fine footprint (i, j, k).
   real(R8P)                :: value_        !< Restricted value.
   real(R8P)                :: uy(1:6)       !< A y line of x-restricted values.
   real(R8P)                :: uz(1:6)       !< The z line of xy-restricted values.
   integer(I4P)             :: n(1:3)        !< Footprint nodes per direction.
   integer(I4P)             :: jj, kk        !< Counters.
   !$acc routine seq
   !$omp declare target

   n = 6_I4P
   where (code == SEAM_RESTRICT_IDENTITY) n = 1_I4P
   uy = 0._R8P
   uz = 0._R8P
   do kk=1, n(3)
      do jj=1, n(2)
         uy(jj) = seam_restriction_weno6_1d(u=u(:,jj,kk), code=code(1))
      enddo
      uz(kk) = seam_restriction_weno6_1d(u=uy, code=code(2))
   enddo
   value_ = seam_restriction_weno6_1d(u=uz, code=code(3))
   endfunction seam_restriction_weno6

   pure function seam_interpolate_weno6(footprint, sub, anchor_pos) result(value_)
   !< WENO coarse->fine fill on a 6x6x6 coarse footprint (issue #68 D2), same index conventions as the quintic
   !< evaluator, direction by direction: x on every (j, k) line, then y, then z (the order of the FNL twin). An
   !< unrefined direction reads the anchor node only, so the other footprint entries along it are not referenced.
   real(R8P),    intent(in) :: footprint(1:6,1:6,1:6) !< Coarse donor values.
   integer(I4P), intent(in) :: sub(1:3)               !< Octant sub-position per direction (1|2|3).
   integer(I4P), intent(in) :: anchor_pos(1:3)        !< Anchor position per direction (1..6).
   real(R8P)                :: value_                 !< Interpolated fine-ghost value.
   real(R8P)                :: uy(1:6)                !< A y line of x-filled values.
   real(R8P)                :: uz(1:6)                !< The z line of xy-filled values.
   integer(I4P)             :: lo(1:3), hi(1:3)       !< Footprint nodes read per direction.
   integer(I4P)             :: i, j, k                !< Counters.
   !$acc routine seq
   !$omp declare target

   do i = 1, 3
      lo(i) = 1_I4P
      hi(i) = 6_I4P
      if (sub(i) == SEAM_SUB_COINCIDENT) then
         lo(i) = anchor_pos(i)
         hi(i) = anchor_pos(i)
      endif
   enddo
   uy = 0._R8P
   uz = 0._R8P
   do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         uy(j) = seam_fill_weno6_1d(u=footprint(:,j,k), sub=sub(1), p=anchor_pos(1))
      enddo
      uz(k) = seam_fill_weno6_1d(u=uy, sub=sub(2), p=anchor_pos(2))
   enddo
   value_ = seam_fill_weno6_1d(u=uz, sub=sub(3), p=anchor_pos(3))
   endfunction seam_interpolate_weno6

   pure function seam_interpolate_quintic_mp(footprint, sub, anchor_pos) result(value_)
   !< Monotonicity-preserving quintic coarse->fine fill on a 6x6x6 coarse footprint (issue #68 D2), direction by
   !< direction (x on every (j, k) line, then y, then z: the order of the FNL twin); an unrefined direction reads the
   !< anchor node only.
   real(R8P),    intent(in) :: footprint(1:6,1:6,1:6) !< Coarse donor values.
   integer(I4P), intent(in) :: sub(1:3)               !< Octant sub-position per direction (1|2|3).
   integer(I4P), intent(in) :: anchor_pos(1:3)        !< Anchor position per direction (1..6).
   real(R8P)                :: value_                 !< Interpolated fine-ghost value.
   real(R8P)                :: uy(1:6)                !< A y line of x-filled values.
   real(R8P)                :: uz(1:6)                !< The z line of xy-filled values.
   integer(I4P)             :: lo(1:3), hi(1:3)       !< Footprint nodes read per direction.
   integer(I4P)             :: i, j, k                !< Counters.
   !$acc routine seq
   !$omp declare target

   do i = 1, 3
      lo(i) = 1_I4P
      hi(i) = 6_I4P
      if (sub(i) == SEAM_SUB_COINCIDENT) then
         lo(i) = anchor_pos(i)
         hi(i) = anchor_pos(i)
      endif
   enddo
   uy = 0._R8P
   uz = 0._R8P
   do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         uy(j) = seam_fill_quintic_mp_1d(u=footprint(:,j,k), sub=sub(1), p=anchor_pos(1))
      enddo
      uz(k) = seam_fill_quintic_mp_1d(u=uy, sub=sub(2), p=anchor_pos(2))
   enddo
   value_ = seam_fill_quintic_mp_1d(u=uz, sub=sub(3), p=anchor_pos(3))
   endfunction seam_interpolate_quintic_mp

   pure function seam_restriction_weight(node, code) result(w)
   !< Return the 1D weight `SEAM_W_LAGRANGE6(node, code)` (issue #68 C2): the device kernels read the table through this
   !< routine, inside the module that owns it, as they read the coarse->fine tables through the evaluators.
   integer(I4P), intent(in) :: node !< Footprint node (1..6).
   integer(I4P), intent(in) :: code !< Restriction code (0..5).
   real(R8P)                :: w    !< Weight.
   !$acc routine seq
   !$omp declare target

   w = SEAM_W_LAGRANGE6(node, code)
   endfunction seam_restriction_weight

   pure subroutine seam_mp_bound(footprint, sub, anchor_pos, lo, hi)
   !< The monotonicity-preserving bound of a coarse->fine fill (issue #101), from the 4x4x4 tricubic footprint of the
   !< ghost (`footprint(i,j,k)` at node offsets `(i,j,k) - anchor_pos`, the conventions of `seam_interpolate_tricubic`),
   !< built direction by direction as `seam_fill_quintic_mp_1d` (Suresh-Huynh 1997). Along each refined direction the
   !< ghost centre lies at eta = (-1)^sub/4 from the anchor:
   !< - inside the footprint, it is a quarter of the way from the anchor to its neighbour on the ghost's side: those two
   !<   nodes bracket it, and the widening is |th (th - 1)/2 dM| (th = 1/4, 3/32 |dM|), dM the M4 minmod of the second
   !<   differences at the two nodes (the one available when the other is the footprint's edge node);
   !< - at the footprint's edge (the anchor is the donor block's last real cell on that side), it is a quarter of a cell
   !<   beyond the anchor: the anchor and its inner neighbour bracket the interval it extrapolates from, and the widening
   !<   is the extrapolation's, (1/4) |slope| + (|c| + 1/4) |d2|, c = th (th - 1)/2 at th = 5/4, d2 the second difference
   !<   at the inner neighbour.
   !< The bound is the range of the bracketing cells (2 per refined direction, the tensor product) widened by the sum over
   !< the directions of the largest widening among the direction's bracketing lines (the curvature of a smooth field
   !< changes from line to line: the anchor line alone misses saddles), plus 64 ulps of the range's magnitude (round-off
   !< of a flat field). The caller clamps the interpolant to it: smooth data, which may graze it near an inflection,
   !< moves by its O(h^3) excess at most; an interpolant overshooting a jump is brought back to the coarse cells' range.
   real(R8P),    intent(in)  :: footprint(1:4,1:4,1:4) !< Coarse donor values.
   integer(I4P), intent(in)  :: sub(1:3)               !< Octant sub-position per direction (1|2|3).
   integer(I4P), intent(in)  :: anchor_pos(1:3)        !< Anchor position per direction (1..4).
   real(R8P),    intent(out) :: lo                     !< Lower bound.
   real(R8P),    intent(out) :: hi                     !< Upper bound.
   integer(I4P)              :: n2(1:3)                !< The other bracketing node per direction (the anchor's if none).
   integer(I4P)              :: b_lo(1:3), b_hi(1:3)   !< Bracketing node range per direction.
   logical                   :: edge(1:3)              !< The ghost extrapolates beyond the anchor along the direction.
   integer(I4P)              :: base(1:3)              !< A bracketing line (its nodes along the direction vary).
   integer(I4P)              :: lm(1:3), lc(1:3), lp(1:3) !< Nodes of a second difference.
   real(R8P)                 :: d2(1:2)                !< Second differences at the anchor's node (1) and at n2 (2).
   logical                   :: has(1:2)               !< Second difference available.
   real(R8P)                 :: dm                     !< M4 minmod curvature.
   real(R8P)                 :: w_line                 !< Widening of one bracketing line.
   real(R8P)                 :: w_dir                  !< Widening of one direction (largest over its lines).
   real(R8P)                 :: widen                  !< Bound widening.
   integer(I4P)              :: d, i, j, k, m          !< Counters.
   !$acc routine seq
   !$omp declare target

   call seam_bracket_nodes(sub=sub, anchor_pos=anchor_pos, n2=n2, edge=edge, b_lo=b_lo, b_hi=b_hi)
   call seam_bracket_range(footprint=footprint, anchor_pos=anchor_pos, b_lo=b_lo, b_hi=b_hi, lo=lo, hi=hi)
   widen = 0._R8P
   do d=1, 3
      if (sub(d) == SEAM_SUB_COINCIDENT) cycle
      has(1) = anchor_pos(d) >= 2_I4P .and. anchor_pos(d) <= 3_I4P
      has(2) = n2(d) /= anchor_pos(d) .and. n2(d) >= 2_I4P .and. n2(d) <= 3_I4P
      w_dir = 0._R8P
      ! the bracketing lines along d: every combination of the other directions' bracketing nodes
      do k=merge(anchor_pos(3), b_lo(3), d == 3), merge(anchor_pos(3), b_hi(3), d == 3)
         do j=merge(anchor_pos(2), b_lo(2), d == 2), merge(anchor_pos(2), b_hi(2), d == 2)
            do i=merge(anchor_pos(1), b_lo(1), d == 1), merge(anchor_pos(1), b_hi(1), d == 1)
               base = [i, j, k]
               d2 = 0._R8P
               do m=1, 2
                  if (.not.has(m)) cycle
                  lc = base
                  lc(d) = merge(anchor_pos(d), n2(d), m == 1)
                  lm = lc
                  lm(d) = lc(d) - 1_I4P
                  lp = lc
                  lp(d) = lc(d) + 1_I4P
                  d2(m) = footprint(lm(1),lm(2),lm(3)) - 2._R8P * footprint(lc(1),lc(2),lc(3)) + &
                          footprint(lp(1),lp(2),lp(3))
               enddo
               if (edge(d)) then
                  lc = base
                  lc(d) = anchor_pos(d)
                  lm = base
                  lm(d) = n2(d)
                  w_line = 0.25_R8P * abs(footprint(lc(1),lc(2),lc(3)) - footprint(lm(1),lm(2),lm(3))) + &
                           0.40625_R8P * abs(d2(2)) ! |c| + 1/4 = 5/32 + 1/4 at th = 5/4
               else
                  if (has(1) .and. has(2)) then
                     dm = seam_minmod4(4._R8P * d2(1) - d2(2), 4._R8P * d2(2) - d2(1), d2(1), d2(2))
                  elseif (has(1)) then
                     dm = d2(1)
                  else
                     dm = d2(2)
                  endif
                  w_line = 0.09375_R8P * abs(dm) ! |th (th - 1)/2| = 3/32 at th = 1/4
               endif
               if (.not.(w_line <= w_dir)) w_dir = w_line ! a NaN widening propagates, alike everywhere
            enddo
         enddo
      enddo
      widen = widen + w_dir
   enddo
   widen = widen + 64._R8P * epsilon(1._R8P) * max(abs(lo), abs(hi))
   lo = lo - widen
   hi = hi + widen
   endsubroutine seam_mp_bound

   pure subroutine seam_bracket_nodes(sub, anchor_pos, n2, edge, b_lo, b_hi)
   !< The coarse nodes that bracket a coarse->fine ghost in its 4x4x4 tricubic footprint (issue #101, `seam_mp_bound`):
   !< per refined direction, the anchor and its neighbour on the ghost's side (eta = (-1)^sub/4), or, when that
   !< neighbour is outside the footprint (the donor block's edge: the ghost extrapolates), the anchor's inner neighbour;
   !< the anchor alone along an unrefined direction.
   integer(I4P), intent(in)  :: sub(1:3)        !< Octant sub-position per direction (1|2|3).
   integer(I4P), intent(in)  :: anchor_pos(1:3) !< Anchor position per direction (1..4).
   integer(I4P), intent(out) :: n2(1:3)         !< The other bracketing node per direction (the anchor's if none).
   logical,      intent(out) :: edge(1:3)       !< The ghost extrapolates beyond the anchor along the direction.
   integer(I4P), intent(out) :: b_lo(1:3)       !< First bracketing node per direction.
   integer(I4P), intent(out) :: b_hi(1:3)       !< Last bracketing node per direction.
   integer(I4P)              :: d, m, s         !< Counter, node, side.
   !$acc routine seq
   !$omp declare target

   do d=1, 3
      n2(d) = anchor_pos(d)
      edge(d) = .false.
      if (sub(d) /= SEAM_SUB_COINCIDENT) then
         s = merge(1_I4P, -1_I4P, sub(d) == 2_I4P)
         m = anchor_pos(d) + s
         if (m >= 1_I4P .and. m <= 4_I4P) then
            n2(d) = m
         else
            n2(d) = anchor_pos(d) - s
            edge(d) = .true.
         endif
      endif
      b_lo(d) = min(anchor_pos(d), n2(d))
      b_hi(d) = max(anchor_pos(d), n2(d))
   enddo
   endsubroutine seam_bracket_nodes

   pure subroutine seam_bracket_range(footprint, anchor_pos, b_lo, b_hi, lo, hi)
   !< The range of the bracketing coarse cells of a ghost (`seam_bracket_nodes`): the unwidened part of `seam_mp_bound`,
   !< contained in it. A value inside it is never clamped, so a caller may skip the full bound (issue #101). Explicit
   !< comparisons: a NaN is handled alike by every compiler (min/max differ).
   real(R8P),    intent(in)  :: footprint(1:4,1:4,1:4) !< Coarse donor values.
   integer(I4P), intent(in)  :: anchor_pos(1:3)        !< Anchor position per direction (1..4).
   integer(I4P), intent(in)  :: b_lo(1:3)              !< First bracketing node per direction.
   integer(I4P), intent(in)  :: b_hi(1:3)              !< Last bracketing node per direction.
   real(R8P),    intent(out) :: lo                     !< Lower end.
   real(R8P),    intent(out) :: hi                     !< Upper end.
   integer(I4P)              :: i, j, k                !< Counters.
   !$acc routine seq
   !$omp declare target

   lo = footprint(anchor_pos(1), anchor_pos(2), anchor_pos(3))
   hi = lo
   do k=b_lo(3), b_hi(3)
      do j=b_lo(2), b_hi(2)
         do i=b_lo(1), b_hi(1)
            if (footprint(i,j,k) < lo) lo = footprint(i,j,k)
            if (footprint(i,j,k) > hi) hi = footprint(i,j,k)
         enddo
      enddo
   enddo
   endsubroutine seam_bracket_range
endmodule adam_seam_interpolation_library
