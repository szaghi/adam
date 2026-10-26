# 1. The Sod shock tube

The first FLUME tutorial reads the Sod input key by key, then changes the reconstruction order and measures what it
buys against the exact solution. It continues from the [first run](/start/first-run): the same build, the same input.

## The input, section by section

`src/tests/flume/regression/sod-x/input.ini`, the sections that define the problem (the full reference is
[FLUME input](/applications/flume/input) and [library input](/reference/library-input)):

```ini
[grid]
ni     = 50          ; cells per block along x
nj     = 4
nk     = 4
ngc    = 3           ; ghost cells: at least the WENO stencil half-width
emin_x = 0.0
emax_x = 1.0
null_x = .false.
null_y = .true.      ; y and z are null axes: a 1-D problem
null_z = .true.

[amr]
max_level      = 2
ratio          = 8   ; octree
iu_ref_levels  = 2   ; two uniform refinements of the root block: 64 blocks, 4 along x, 200 cells
markers_number = 0   ; no adaptive refinement

[runge_kutta]
scheme = runge-kutta-ssp-33

[weno]
scheme = weno-u-5    ; upwind WENO, 5th order

[physics]
physical_model = euler
cp             = 1040.004
cv             = 742.86          ; gamma = cp/cv = 1.4

[numerics]
scheme_space             = weno            ; WENO flux splitting
reconstruction_variables = characteristic

[initial_conditions]
type           = riemann-problem
regions_number = 2               ; [initial_conditions_region_1] (x < 0.5) and _2 (x > 0.5)

[time]
it_max   = -1                    ; run to time_max
time_max = 0.2
CFL      = 0.5
```

The two regions hold the left state $(\rho, u, p) = (1, 0, 1)$ and the right state $(0.125, 0, 0.1)$; all six faces
are `extrapolation`.

## Change the reconstruction order

Copy the input into three directories and change `[weno] scheme` (and `ngc` where the stencil needs it):

| `[weno] scheme` | `ngc` | $L_1(\rho)$ against the exact solution, 200 cells, $t = 0.2$ |
|---|---|---|
| `weno-u-3` | 3 | $4.585912\cdot10^{-3}$ |
| `weno-u-5` (the input) | 3 | $3.244140\cdot10^{-3}$ |
| `weno-u-7` | 4 | $2.823079\cdot10^{-3}$ |

Each row is a run on 2 ranks checked with

```bash
../../exe/.regression-venv/bin/python ../../src/tests/flume/verification/sod/sod_oracle.py input.ini .
```

The error falls with the order but not as fast as the formal order suggests: on a problem made of a shock, a contact
and a rarefaction the $L_1$ error is dominated by the discontinuities, where every scheme drops to first order.

**`ngc` is checked.** `weno-u-7` reads 4 cells on each side of a face. With `ngc = 3` the run stops at initialization:

```text
[mpi-00000]error stop : [grid].(ngc)=+3 is smaller than the WENO stencil half-width +4
```

## The conservation history

`sod-x-conservation_history.dat` integrates the conservative variables over the domain at every step. Between step 0
and step 174 the mass stays at $0.5625$ (to the last digit), and the momentum grows from 0 to
$0.18 = (p_L - p_R)\,t = (1 - 0.1)\times0.2$: the extrapolated ends carry the pressure difference, the only flux the
waves have not yet reached.

::: tip What you learned
- An input is a set of sections: the grid and its blocks, the AMR tree, the schemes, the physics, the initial
  condition, the time.
- Null axes make a 1-D problem out of the 3-D machinery.
- The WENO order is one key; the ghost width must follow it, and FLUME refuses a run where it does not.
- The exact solution and the conservation history are the two checks you can always run.
:::

Next: [2. The Lax problem](./lax).
