# Stratified atmospheres and the `hydrostatic` boundary

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| GV (G0) | External gravity | M7 ([#73](https://github.com/szaghi/adam/issues/73)) | CPU, FNL | yes | — |

## Case and scope

An atmosphere at rest under a constant gravity $\mathbf g$, in hydrostatic equilibrium, closed along $\mathbf g$ by
`hydrostatic` faces. Every later gravity test starts from such an atmosphere: G0 checks that the initial state and its
ghosts are what the formulas say, before any force acts on them. Each case runs one negligible step with the ghosts
poisoned, so the checkpoint holds the initial atmosphere and its ghosts filled in the order of a Runge–Kutta stage.

What the case catches:
- **a wrong initial atmosphere:** every interior cell against the analytic isothermal or polytropic profile, to
  $10^{-13}$;
- **a wrong `hydrostatic` boundary rule:** every ghost beyond a `hydrostatic` face against the rule (the isothermal
  continuation of the boundary cell), with zero velocity and the field unchanged; on oblique gravity the realm edges
  and corners too;
- **a ghost read before it is filled:** the ghosts are poisoned with NaN, and a NaN in any checked cell fails;
- **an oracle that cannot tell the boundary rule from a copy:** a negative control with `extrapolation` faces must be
  rejected;
- **a wrong conversion of the gravity input** by the `[reference]` layer (NV-5): the components of $\mathbf g$ convert
  as an acceleration, `froude` does not.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Equations | compressible Euler, ideal gas, $\gamma = 1.4$; `iso2d-mhd`: ideal MHD with EGLM cleaning (`glm_ch = 3`, `glm_alpha = 0.18`, `glm_damping_length = 1`) | [models](../../models#external-gravity), [EGLM](../../models#divergence-control-eglm) |
| Gravity input | `[gravity]` components `g_x`, `g_y`, `g_z`, or the number form `froude` + `axis` ($1/\mathrm{Fr}^2$ along the signed axis) | `load_from_file`, [`adam_flume_gravity_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_gravity_object.F90) |
| Atmosphere | isothermal $\rho = \rho_1 e^{\mathbf g\cdot\mathbf x/T_1}$, $p = \rho T_1$, $T_1 = p_1/\rho_1$; polytropic $\rho = \rho_1 s^{1/(G-1)}$, $p = p_1 s^{G/(G-1)}$, $s = 1 + \frac{G-1}{G}\frac{\rho_1}{p_1}\mathbf g\cdot\mathbf x$; region 1 at the origin | `atmosphere_state`, [`adam_flume_ic_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_ic_object.F90) |
| `hydrostatic` boundary | ghost density and pressure continue the isothermal equilibrium of the first interior cell along the normal, $\rho_g = \rho_1 e^{g_n s/T_1}$, $p_g = \rho_g T_1$ ($s$ the signed distance from that cell to the ghost); velocity mirrored as at `wall-inviscid` (normal component odd), field and $\psi$ even (Zingale et al. 2002) | `hydrostatic_ghost`, [`adam_flume_bc_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/common/adam_flume_bc_object.F90); `set_boundary_conditions` (`hydrostatic_ghost_cell`), [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); `set_boundary_conditions_dev`, [`adam_flume_fnl_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_kernels.F90) |
| Ghost poisoning | `[diagnostics] ghost_poison = .true.`: NaN in every ghost before the fill | `poison_ghosts`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90) |
| Space | WENO5 (`weno-u-5`) flux splitting on characteristic variables (`weno`) | `compute_residuals_weno`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); [numerics](../../numerics#scheme-space-weno-flux-splitting) |
| Time | SSP Runge–Kutta, 3 stages (`runge-kutta-ssp-33`), one step at CFL $10^{-30}$ | [`adam_rk_object.F90`](https://github.com/szaghi/adam/blob/develop/src/lib/common/adam_rk_object.F90) |
| AMR (`tilted2d`) | static refinement box, 2:1 seams (`reflux = .true.`) | [numerics](../../numerics#adaptive-mesh-refinement) |
| FNL backend | the same kernels on the device | [`adam_flume_fnl_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_object.F90) |

The G0 cases do not exercise the gravity sources (they came with G1): one step at CFL $10^{-30}$ leaves the state
unchanged to round-off.

## Initial and boundary conditions

| Case | Atmosphere | Gravity | Faces | Grid |
|---|---|---|---|---|
| `iso2d` | isothermal | $(0, -1)$ | x `periodic`, y `hydrostatic` | $[0, 1]^2$ quadtree, 4×4 blocks of $8^2$ cells ($32^2$) |
| `poly2d` | polytropic, $G = 1.2$ | $(0, -1)$ | as `iso2d` | as `iso2d` |
| `tilted2d` | isothermal | $(-0.6, -0.8)$ | every face `hydrostatic` (the realm corners too) | as `iso2d`, the quadrant $[0, 0.5]^2$ refined one level (2:1 seams; 28 blocks) |
| `tilted3d` | isothermal | $(-0.48, -0.64, -0.6)$ | every face `hydrostatic` | $[0, 1]^3$ octree, 4×4×4 blocks of $8^3$ cells ($32^3$) |
| `iso2d-mhd` | isothermal, EGLM, $\mathbf B = (0.5, 0.7, 0.2)$ (oblique to $\mathbf g$) | $(0, -1)$ | as `iso2d` | as `iso2d` |
| `froude2d` | isothermal | `froude = 1`, `axis = -y` | as `iso2d` | as `iso2d` |

| Item | Value |
|---|---|
| Region 1 (at the origin) | $(\rho_1, \mathbf u, p_1) = (1, \mathbf 0, 1)$ |
| Initial condition | [`isothermal-atmosphere`, `polytropic-atmosphere`](../../initial-conditions#isothermal-atmosphere-polytropic-atmosphere-hydrostatic-equilibria-under-gravity) |
| Grid | two uniform levels (`iu_ref_levels = 2`; with one level every block stayed on rank 0, so an np 2 run decomposed nothing), `ngc = 4`; 2-D cases `extrapolation` on the null z faces |
| Ranks | 1 and 2 |
| Steps | one, at CFL $10^{-30}$ |
| Negative control | `iso2d` with `extrapolation` on the y faces, np 1, judged as if they were `hydrostatic` |
| NV-5 | `scaling/reference.sh --inputs` on `iso2d`, `iso2d-mhd`, `froude2d` |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/gravity/make_atmosphere.py` | writes the input of each case (`<out.ini> <case>`) |
| `src/tests/flume/verification/gravity/check.sh` | runs the cases (`--np`, `--cases`, `--no-force`, `--reference` for NV-5, `--keep`; `FLUME_EXE` selects the backend) |
| `src/tests/flume/verification/gravity/atmosphere_oracle.py` | `profile`: the interior and the ghosts against the closed form, tolerance $10^{-13}$ |
| `src/tests/flume/verification/scaling/reference.sh` | NV-5, the `[reference]` layer on the gravity inputs |
| `src/tests/flume/verification/make_doc_runs.sh` (group `gravity`) | re-runs `iso2d`, `poly2d` with their checkpoints |
| `src/tests/flume/verification/make_doc_figures.py` (`fig_gravity_atmospheres`) | draws the figure |

```bash
cd src/tests/flume/verification/gravity
./check.sh --no-force                                          # G0 only, CPU, np 1 and 2
./check.sh --no-force --reference                              # with NV-5
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh --no-force   # FNL
```

## Results

![GV (G0): the iso2d and poly2d atmospheres with their hydrostatic ghosts, and their relative error](/flume/gravity-atmospheres.png)

*Left, the two atmospheres along y: interior cells and the ghosts beyond the y faces against the analytic profiles;
right, the relative error of every one of those cells.*

Left, every interior cell (dots) and every ghost beyond the y faces (crosses, shaded). Beyond a polytropic face the
ghosts leave the polytropic curve on purpose: the `hydrostatic` rule extends the boundary cell isothermally, and they
are checked against that extension. Right, the relative error of every one of those cells: zero or one unit in the
last place, ghosts included.

Every case PASS, CPU and FNL, np 1 and 2:

| Check | Cells (per run) | Bound | CPU | FNL |
|---|---|---|---|---|
| interior $\rho$, $p$ | 1 024 (2-D), 32 768 (3-D) | $10^{-13}$ | $2.2\cdot10^{-16}$ | $6.7\cdot10^{-16}$ |
| hydrostatic ghosts $\rho$, $p$ | 192 (`iso2d`) to 22 104 (`tilted3d`, edges and corners) | $10^{-13}$ | $\le 8.9\cdot10^{-16}$ | $\le 8.9\cdot10^{-16}$ |
| hydrostatic ghosts $\mathbf u$, $\mathbf B$ | as above | $10^{-13}$ | velocity 0 to $10^{-32}$, field uniform | the same |
| negative control | 192 | must exceed $10^{-13}$ | $9.8\cdot10^{-2}$ (rejected, as it must be) | the same |
| NV-5 | `iso2d`, `iso2d-mhd`, `froude2d` | bitwise | bitwise; conversions exact and complete | bitwise |

The ghost counts date from G0, when the inputs had `ngc = 3`: 192 = 2 faces × 32 columns × 3 rows. G2a raised the
inputs to `ngc = 4` ([energy order](./energy-order)), so the current `iso2d` checks 256 hydrostatic ghosts.

Related cases: [ghost probe](../infrastructure/ghost-probe) checks the exchange of the interface and seam ghosts that
this case counts but does not check; [scaling covariance](../infrastructure/scaling-covariance) holds the
`[reference]` layer of NV-5; [the gravity force](./force) and [the well-balanced mode](./well-balanced) start from
these atmospheres.

## Findings

- Asserted: the interior against the analytic atmosphere; every ghost beyond a `hydrostatic` face against the boundary
  rule. For the isothermal atmosphere the rule is the analytic profile, so the realm edges and corners are checked
  too; for the polytropic one, the ghosts beyond one face only (the rule is per face). The velocity in those ghosts
  must be zero and the field uniform. The periodic ghosts are checked against the analytic profile at the wrapped
  position.
- Boundary ghosts beside their own block (their donor along the normal is an interface or seam ghost, not exact on an
  exponential) are counted, not checked; the exchange is the job of the [ghost probe](../infrastructure/ghost-probe)
  (GP).
- The negative control proves the oracle discriminates: `extrapolation` ghosts copied from the boundary cell miss the
  hydrostatic rule by $9.8\cdot10^{-2}$.
- NV-5: the gravity components convert as an acceleration ($u_0^2/L_0$), `froude` does not, and the dimensional run
  equals the base one bit for bit.
- Two uniform levels are needed to spread the blocks over the ranks: with one level every block stayed on rank 0, so
  an np 2 run decomposed nothing.
- **Found on the way** (pre-existing, not gravity): on FNL an MHD run died when a rank owned no blocks
  (`dev_assign_to_device: device memory not allocated`), with or without gravity and with the executable before G0;
  the Euler model and the CPU backend ran. The G0 cases give every rank blocks. Fixed and asserted since, see
  [ranks without blocks](../infrastructure/ranks-without-blocks) ([#82](https://github.com/szaghi/adam/issues/82)).
