# FDV coefficient generator — `fdv_coefficients.py`

Offline, exact-rational source of truth for the finite-difference (FD) and
finite-volume (FV) stencil coefficients hard-coded in
[`src/lib/common/adam_fdv_operators_library.F90`](../src/lib/common/adam_fdv_operators_library.F90).

The Fortran module ships the coefficients as `parameter` tables (compile-time
constants, `!$acc declare copyin` to the device). This script is **not** part of
the build: it derives those same coefficients from first principles in exact
rational arithmetic (`sympy.Rational`) and checks them bit-for-bit, so the tables
can be audited, regenerated, and extended to higher order without hand-computing
rationals.

> **Design decision.** Coefficients are generated *offline* and committed as
> Fortran `parameter` source — they are **not** computed at program init. An
> earlier alternative (host-side runtime generator + `!$acc declare create` +
> `update device`) was prototyped and validated on nvfortran (Phase 0) but
> deliberately not adopted: committed `parameter` tables keep the runtime cost at
> zero, the device-copy path unchanged, and every coefficient diff-auditable.

---

## Usage

```bash
# from the repository root; requires sympy (1.12+)
python3 scripts/fdv_coefficients.py
```

The script prints one line per `(family, s)` pair and a final verdict. Exit code
`0` = all generated tables match the committed Fortran bit-for-bit; `1` = at
least one mismatch (the offending rows are printed with `generated` vs
`committed`).

```
[OK ] FV1_CC s=3: generated=[37/60, -2/15, 1/60] committed=[37/60, -2/15, 1/60]
...
VALIDATION PASSED: all generated tables match committed Fortran bit-for-bit
```

`s` is the **half stencil length**; the classical order of accuracy is `2*s` for
the centered families. `S_MAX = 5` mirrors the current Fortran limit (10th order).

### Roles

- **Validation gate** (current use): regression-check the committed tables before
  any rebuild work. Run it after touching the Fortran coefficient blocks.
- **Generation** (future use): the per-family functions return exact `Rational`
  weights for *any* `s`, so they are the source for raising `S_MAX` and for the
  new SOTA pieces (higher-derivative FV tables, deconvolution) that have no
  committed counterpart yet.

---

## What it generates, and the theory behind each family

All families share one idea: find weights `w_j` on a set of integer grid nodes
(in units of the spacing `h`, set to 1) that match a target linear functional
(a derivative or a face value) exactly for all monomials up to the stencil's
degree. This is **moment matching** — solving a small linear system in exact
rationals. At these sizes the system is trivially solvable; the classical
conditioning worry that motivates Fornberg's recurrence does not arise because we
work in exact arithmetic, not floating point.

### 1. Finite-difference centered derivatives — `FD0..FD6_CC`

Pointwise nodal data `q_i = q(x_i)`. For the `n`-th derivative at `x_i`, weights
`c_m` on nodes `m = -s..s` satisfy

```
sum_{m=-s}^{s} c_m * m^k = n! * delta_{k,n}   for k = 0..2s
```

so that `d^n q/ds^n |_i ≈ (1/h^n) * sum_m c_m q_{i+m}` is accurate to `O(h^{2s})`
(odd `n`) or `O(h^{2s})` (even `n`, with the symmetric stencil). Storage exploits
parity:

- **odd `n` (FD1, FD3, FD5)** — antisymmetric (`c_{-m} = -c_m`, `c_0 = 0`); store
  `[c_1..c_s]`, apply to pairs `q(1+m) - q(1-m)`.
- **even `n` (FD2, FD4, FD6)** — symmetric (`c_{-m} = c_m`); store `[c_0, c_1..c_s]`,
  apply `c_0 q(1)` plus pairs `c_m (q(1+m) + q(1-m))`.
- **FD0** — the special case `n = 0` evaluated at the **face** `i+1/2`: a Lagrange
  *reconstruction* of the face value from point samples, symmetric about the face.

When `n > 2s` the stencil is too small to represent the derivative; the table
stores an all-zero row (e.g. FD3 at `s=1`, FD5 at `s≤2`). The generator reproduces
these zero rows.

### 2. Finite-volume centered reconstruction — `FV1_CC`

Here the data are **cell averages** `⟨q⟩_i = (1/h) ∫_{cell i} q dx`, not point
values — the distinction is the heart of high-order FV. Weights `a_j` on cell
offsets `j = 1-s..s` reconstruct the face value:

```
sum_j a_j * avg_k(j) = (1/2)^k   for k = 0..2s-1,
   where  avg_k(j) = (1/(k+1)) * [ (j+1/2)^{k+1} - (j-1/2)^{k+1} ]
```

`avg_k(j)` is the cell-average of the monomial `x^k` over the unit cell centered at
`j`; the right-hand side is the point value of `x^k` at the face `x = +1/2`. This
is exactly the **primitive-function reconstruction**: interpolating the indefinite
integral `F = ∫q` (whose face point-values are exact running sums of the cell
averages) and differentiating once. Solving the system *is* differentiating that
interpolant. The face derivative gives the conservative FV first derivative
`dq/ds|_i ≈ (q_{i+1/2} − q_{i-1/2})/h`, whose telescoping guarantees discrete
conservation.

### 3. Finite-volume upwind reconstruction — `FV1_UR`, `FV1_UL`

Same moment match as `FV1_CC` but on a **one-sided** cell set. The right-upwind
reconstruction of the face `i+1/2` uses the `s` cells *at and behind* the face,
offsets `0, -1, ..., -(s-1)` (the upwind side for a right-moving wave). `FV1_UL` is
the mirror image (`FV1_UL = reverse(FV1_UR)`).

### 4. (Future) cell-average ↔ point-value deconvolution

Not yet in the Fortran module; the gap that caps high-order FV of nonlinear/product
terms at 2nd order. To 4th order,

```
⟨q⟩_i = q_i + (h^2/24) ∇^2 q_i + O(h^4),     q_i = ⟨q⟩_i − (h^2/24) ∇^2 ⟨q⟩_i + O(h^4)
```

with `∇^2` an undivided second-difference Laplacian. The same moment-matching
machinery generates the arbitrary-order `2p` conversion as a linear combination of
even undivided differences.

---

## Storage convention (matching the Fortran)

- Tables are stored on the **half stencil** (length `s`, or `s+1` for the even-FD
  `[c_0, c_1..s]` layout), exploiting (anti)symmetry.
- Some families (notably `FD0_CC`) store **un-reduced** numerators over a common
  per-row divisor (e.g. `150/256` rather than the reduced `75/128`) for
  readability. The validator compares normalized `Rational`s, so reduced and
  un-reduced forms compare equal — but when *emitting* new Fortran source, preserve
  the un-reduced common-divisor house style rather than auto-reducing.

---

## Theoretical references

- **B. Fornberg**, *Calculation of weights in finite difference formulas*, SIAM
  Review 40(3), 685–691 (1998); and *Generation of finite difference formulas on
  arbitrarily spaced grids*, Math. Comp. 51, 699–706 (1988). The stable recurrence
  for FD weights and the moment-matching formulation used here.
- **P. Colella, M. R. Dorr, J. A. F. Hittinger, D. F. Martin**, *High-order,
  finite-volume methods in mapped coordinates*, J. Comput. Phys. 230(8), 2952–2976
  (2011). Canonical high-order FV framework: face-averaged flux reconstruction,
  the cell-average↔point-value relation `⟨q⟩ = q + (h²/24)∇²q`, and the
  face-gradient (diffusion) stencils.
- **A. Mignone**, *High-order conservative reconstruction schemes for finite
  volume methods in cylindrical and spherical coordinates*, J. Comput. Phys. 270,
  784–814 (2014). The primitive-function reconstruction stated cleanly (interpolate
  the indefinite integral, differentiate); the arbitrary-order avg↔point
  conversion (Eq. 52).
- **C.-W. Shu**, *High order weighted essentially non-oscillatory schemes for
  convection dominated problems*, SIAM Review 51(1), 82–126 (2009). Reference for
  the reconstruction-as-primitive viewpoint and upwind-biased FV stencils.

---

## See also

- [`src/lib/common/adam_fdv_operators_library.F90`](../src/lib/common/adam_fdv_operators_library.F90)
  — the consuming Fortran module (CPU/common backend).
- [`src/lib/fnl/adam_fnl_fdv_operators_library.F90`](../src/lib/fnl/adam_fnl_fdv_operators_library.F90)
  — FNL (OpenACC) device backend.
- `src/tests/fdv_operators/` — the trigonometric order-of-accuracy tests that gate
  the *new* (uncommitted) higher-derivative and deconvolution tables, which have no
  bit-for-bit reference.
