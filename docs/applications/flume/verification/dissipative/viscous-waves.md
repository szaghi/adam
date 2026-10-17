# Viscous waves

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| VV-1, VV-2 | Dissipative terms (M4) | M4 ([#65](https://github.com/szaghi/adam/issues/65)) | CPU, FNL | yes | `shear-viscous` |

## Case and scope

A small sine wave on a gas at rest is one Fourier mode of the linearised compressible Navier–Stokes equations. The
oracle (`viscous/waves.py`) integrates that mode exactly, $\exp(Mt)$ of the $3\times3$ system for $(\rho', u_n, T')$
plus the decoupled transverse velocity, from the Fourier coefficients of the initial state, so it holds every mode the
start excites (the acoustic pair and the entropy mode), not only a decay rate. The error is the RMS over the cells of
the point values (FLUME is a finite-difference scheme).

- **VV-1, the shear wave:** the transverse velocity $u_t = A\sin(k\,\mathbf n\cdot\mathbf x)\,e^{-\nu k^2 t}$; the
  error is measured on $u_t$. It pins the order of the viscous flux.
- **VV-2, the right-running acoustic wave** with a conductivity as well, so $\mu$ and $k$ act together; the error is
  measured on $\rho'$.

What the case catches: a dissipative flux that is consistent and still wrong in its **order** (a missing correction
term), its **coefficient** (a missing 4/3, a conductivity without $R$) or its **direction** (a cross derivative read
with the wrong stencil, exercised by the diagonal leg).

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Navier–Stokes, ideal gas, $\gamma = 1.4$, $R = 1$, constant $\mu$ and $k$ | [models](../../models#compressible-euler) |
| Dissipative fluxes | `dissipative_order = 4` (conservative fourth-order flux with the Shu–Osher correction) or `2` (compact central) | `add_dissipative_fluxes`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `add_dissipative_fluxes_o4`, `add_dissipative_fluxes_o2`, [`adam_flume_cpu_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_dissipation_kernels.F90); `compute_dissipative_flux`, [`adam_flume_dissipation_library.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_dissipation_library.F90); [numerics](../../numerics#dissipative-fluxes-navier–stokes) |
| Inviscid fluxes | WENO5 (`weno-u-5`) flux splitting on characteristic variables (`scheme_space = weno`) | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Time | SSP Runge–Kutta, 5 stages, 4th order (`runge-kutta-ssp-54`), CFL 0.2, the diffusive limit included | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90); `compute_lambda_max_dissipative`, [`adam_flume_cpu_euler_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_euler_kernels.F90) |
| Boundaries | periodic in $x$ and $y$, `extrapolation` on the null $z$ | `set_boundary_conditions`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| FNL backend | the same fluxes as device kernels | `add_dissipative_fluxes_dev`, [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90); `add_dissipative_fluxes_o4_dev`, [`adam_flume_fnl_dissipation_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_dissipation_kernels.F90) |

## Initial and boundary conditions

| Item | Value |
|---|---|
| Domain | $[0, 1]^2$, $z$ null |
| Grid | quadtree, $4 \times 4$ blocks, $N/4$ cells per block along the active axes, `ngc = 3`. A wave along an axis keeps the other in-plane axis: active and thin for the shear wave (4 blocks of 6 cells across it, because FLUME freezes the momentum normal to a null direction, which would freeze the shear), null (1 cell) for the acoustic wave |
| Ranks | 2 |
| Initial condition | [`sine-wave`](../../initial-conditions#sine-wave-a-small-sine-wave-on-a-uniform-state-verification), $\rho_0 = p_0 = 1$, $u = 0$; wavelength 1 along an axis, $1/\sqrt2$ along the diagonal |
| Boundaries | `periodic` on the $x$ and $y$ faces, `extrapolation` on the $z$ faces |

| Leg | Mode | Direction | `dissipative_order` | $A$ | $\mu$, $k$ | $N$ | $t$ |
|---|---|---|---|---|---|---|---|
| VV-1 shear x | shear | $x$ | 4 | $10^{-4}$ | 0.01, 0 | 32 / 64 / 128 | 0.5 |
| VV-1 shear y | shear | $y$ | 4 | $10^{-4}$ | 0.01, 0 | 32 / 64 / 128 | 0.5 |
| VV-1 shear x, order 2 | shear | $x$ | 2 | $10^{-4}$ | 0.01, 0 | 32 / 64 / 128 | 0.5 |
| VV-1 shear diagonal | shear | 45° (both directions active: the cross terms $\partial_i u_d$ and the tangential derivatives) | 4 | $10^{-4}$ | 0.01, 0 | 24 / 48 / 96 | 0.1 |
| VV-2 acoustic | acoustic, and its $-A$ twin | $x$ | 4 | $10^{-5}$ | 0.01, 0.02 | 32 / 64 / 128 | 0.5 |

A sound wave steepens at $O(A^2)$: at $A = 10^{-4}$ that residue ($10^{-8}$) flattened the error at every resolution.
VV-2 therefore runs the twin with amplitude $-A$ and judges the odd part $(f(+A) - f(-A))/2$, whose residue is
$O(A^3)$, at $A = 10^{-5}$.

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/viscous/waves.py` | `make` writes each input; `oracle` integrates the linear mode exactly and judges the ladder |
| `src/tests/flume/verification/viscous/check.sh --legs vv1,vv2` | runs the ladders (`--np N`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/regression/shear-viscous/input.ini` | the frozen diagonal shear wave, $48^2$, $t = 0.05$ |
| `make_doc_figures.py` (`fig_viscous`, `fig_viscous_fields`, `fig_viscous_ladders`) | the figures, from the CPU (and, for the ladders, FNL) work directories of `viscous/check.sh` |

```bash
cd src/tests/flume/verification/viscous && ./check.sh --legs vv1,vv2             # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --legs vv1,vv2       # FNL
```

## Results

![VV-1 to VV-5: convergence, and Becker's shock against the exact profile](/flume/viscous.png)

*Left: the Navier–Stokes ladders (VV-1 to VV-4), each normalised by its coarsest error, against the order-2 and
order-4 slopes; middle: the Ohmic ladders (VV-5); right: Becker's shock (VV-4).*

![VV-1, VV-2 and VV-5: absolute errors on CPU and FNL, and the Ohmic decay](/flume/viscous-ladders.png)

*Absolute RMS errors of the VV-1 and VV-2 legs (left) and of the VV-5 legs (centre), CPU (filled, solid) and FNL
(open, dotted), which coincide; right, the Ohmic decay: the measured amplitude of the mode over the exact one.*

![VV-1 and VV-8 in 2-D: the shear wave and its error on a uniform grid; the shear and Ohmic errors with the centre refined 2:1](/flume/viscous-fields.png)

*Top row: the diagonal shear wave on the uniform $96^2$ grid at $t = 0.1$ and its pointwise error, smooth and spread
over the domain. The bottom row belongs to [VV-8](./seam-accuracy).*

Oracle: the observed order of the finest pair $\ge 3.8$ at order 4, $\ge 1.9$ at order 2.

| Leg | N | RMS error (CPU) | Order (CPU) | FNL |
|---|---|---|---|---|
| shear x, order 4 | 32 / 64 / 128 | $3.64\cdot10^{-10}$ / $2.29\cdot10^{-11}$ / $1.43\cdot10^{-12}$ | +3.99, +4.00 | the same to the printed digits |
| shear y, order 4 | 32 / 64 / 128 | identical to x, bit for bit | +3.99, +4.00 | the same |
| shear x, order 2 | 32 / 64 / 128 | $3.68\cdot10^{-8}$ / $9.20\cdot10^{-9}$ / $2.30\cdot10^{-9}$ | +2.00, +2.00 | the same |
| shear diagonal, order 4 | 24 / 48 / 96 | $6.72\cdot10^{-10}$ / $1.16\cdot10^{-11}$ / $2.43\cdot10^{-13}$ | +5.86, +5.57 | the same |
| acoustic, order 4 | 32 / 64 / 128 | $8.20\cdot10^{-11}$ / $1.69\cdot10^{-12}$ / $2.18\cdot10^{-14}$ | +5.60, +6.28 | $8.20\cdot10^{-11}$ / $1.69\cdot10^{-12}$ / $2.21\cdot10^{-14}$, +5.60, +6.26 (at 32 cells $8.196853\cdot10^{-11}$ against the CPU's $8.196891\cdot10^{-11}$) |

Every leg passes on both backends (2 ranks). On FNL the numbers are the same to the printed digits except where the
round-off of the run itself shows: the acoustic odd part, whose $\pm A$ difference amplifies the round-off of a
$10^{-5}$ wave.

Related cases: the same oracle drives [Ohmic decay and the Alfvén wave](./ohmic-alfven); the diagonal waves cross 2:1
seams in [VV-7](./seam-conservation) and [VV-8](./seam-accuracy); the run-by-run backend comparison is in
[CPU and FNL agree](./cpu-fnl-agreement).

## Findings

- The measured decay rate of the shear wave is $\nu k^2$ to $2\cdot10^{-8}$ relative at 128 cells.
- The diagonal wave converges faster than its design order: on the diagonal the leading truncation errors of the $x$
  and $y$ fluxes appear to cancel.
- On the acoustic wave the error is dominated by the inviscid WENO part, not by the viscous one. VV-2 therefore checks
  that $\mu$ and $k$ enter with the right coefficients (the measured amplitude equals the exact one to
  $6\cdot10^{-10}$ relative), and VV-1 pins the order of the dissipative flux.
- The measured acoustic rate, $0.29576$, differs by 4% from the first-order Kirchhoff–Stokes rate
  $k^2/(2\rho_0)\,(\tfrac43\mu + (\gamma-1)k/c_p) = 0.30831$: here $\nu k/a = 0.053$, and the exact linear oracle
  carries that correction.
- #85 lists absolute shear and acoustic ladders (RMS error against $N$ with reference slopes) as a figure still to
  extend `fig_viscous` with; the published figure shows the ladders normalised by their coarsest error.
