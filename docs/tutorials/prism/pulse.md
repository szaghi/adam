# 1. An electromagnetic pulse in a periodic box

The first PRISM tutorial runs a Gaussian pulse once around a periodic box: the build, the input, the energy history,
and how far the pulse comes back to where it started.

## Build

```bash
fobis build --mode prism-cpu-gnu         # exe/adam_prism_cpu; links FFTW (the t0 elliptic solves)
```

## The input

`src/tests/prism/regression/periodic-pulse-blocks/input.ini`, the keys that define the problem (the full reference is
[PRISM input](/applications/prism/input)):

```ini
[physics]
physical_model = EM              ; Maxwell's equations in D and B, SI units

[numerics]
scheme_space          = fd_centered
scheme_time           = runge_kutta
constrained_transport = no
divergence_correction = no

[fdv]
fdv_scheme = fd
fdv_order  = 6                   ; 6th-order centred derivatives

[runge_kutta]
scheme = runge-kutta-ssp-54

[grid]
ni     = 12                      ; 12^3 cells per block
emin_x = -0.1                    ; the box [-0.1, 0.1]^3 (metres)
emax_x = 0.1

[amr]
iu_ref_levels = 2                ; 64 blocks: 48^3 cells, dx = 0.2/48

[initial_conditions]
type            = gaussian_pulse
pulse_direction = +x
r0              = 0.07           ; centre of the packet along x
sigma           = 0.015          ; width
k0              = 0.0            ; no carrier: a pure Gaussian
B0              = 10.0

[bc_x_min]
type = periodic                  ; all six faces periodic

[time]
it_max   = 192
CFL      = 0.5
```

The pulse is $E_y = c_0B_0\,e^{-(x - r_0)^2/2\sigma^2}$, $B_z = E_y/c_0$, travelling along $+x$
([initial conditions](/applications/prism/initial-conditions)). The time step is $\Delta t = \Delta x/(4c_0)$, so 192
steps are exactly one crossing of the box: the pulse leaves through $+x$, comes
back through $-x$ and returns to its start.

## Run it

```bash
mkdir -p run/pulse && cp src/tests/prism/regression/periodic-pulse-blocks/input.ini run/pulse/
cd run/pulse && mpirun -np 2 ../../exe/adam_prism_cpu input.ini > log.txt
```

It writes the fields every 48 steps (`periodic_pulse_blocks-000000000.xdmf` … `-000000192.xdmf`), the restart files,
`-residuals.dat`, `-divergence_history.dat` and `-energy_history.dat`.

## What to look at

**Energy** (`periodic_pulse_blocks-energy_history.dat`, every 10 steps; electric and magnetic energy in joules):

<<< @/snippets/tutorials/prism-pulse-energy.txt

The two energies stay equal (a travelling wave) and their sum drifts by $1.2\cdot10^{-6}$ relative over the period:
the centred scheme is non-dissipative, and the drift is the Runge–Kutta's.

**Divergence**: the last row of `-divergence_history.dat` is zero for $\nabla\cdot\mathbf D$ and $\nabla\cdot\mathbf B$:
on a uniform grid the centred divergence of the centred curl vanishes identically
([numerics](/applications/prism/numerics)).

**The return.** `src/tests/tutorials/compare_saves.py` compares two saves of a run, block by block (it needs h5py: the
venv of the [first run](/start/first-run#check-it)):

```bash
../../exe/.regression-venv/bin/python ../../src/tests/tutorials/compare_saves.py . periodic_pulse_blocks 0 192
```

<<< @/snippets/tutorials/prism-pulse-return.txt

A 4% difference after one crossing, and most of it is not the scheme: the initial pulse is set as a Gaussian about
`r0 = 0.07` that is **not wrapped across the periodic faces**. The last cells before the face $x = 0.1$, $1.86\sigma$
from the centre, hold 18% of the peak, and the cells beyond it, the periodic images of $x = -0.1$, almost nothing: the
field jumps there, the jump travels with the pulse and dominates the return error, and refining the grid does not
remove it. The same run with the pulse centred in the box (`r0 = 0.0`, as the [next page](./faces) does), where
nothing is cut, shows the scheme's own error:

<<< @/snippets/tutorials/prism-pulse-centred.txt

::: tip What you learned
- PRISM's input names the model, the space and time schemes, the grid, the initial field and the faces.
- The energy and divergence histories are the first checks: energy flat, divergence zero on a uniform grid.
- The 4% return error here is the unwrapped initial cut, not the scheme; centre the pulse to see the scheme's own error.
:::

Next: [2. Reflecting and absorbing faces](./faces).
