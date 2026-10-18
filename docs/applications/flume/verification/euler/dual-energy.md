# Dual energy: the advected entropy of euler-de

| ID | Group | Milestone | Backends | Figure | Regression |
|---|---|---|---|---|---|
| DE (de1 to de4) | Euler | M9 ([#91](https://github.com/szaghi/adam/issues/91)), P2, P3 | CPU, FNL | yes | — |

## Case and scope

In a cold supersonic flow the internal energy $e = E - |\rho\mathbf u|^2/(2\rho)$ is a small difference of large
numbers, and the truncation error of the total-energy scheme can drain it to zero and below: the runtime-AMR collapse of
[#88](https://github.com/szaghi/adam/issues/88) died that way. `physical_model = euler-de` carries the advected entropy
$\rho s$, $s = p/\rho^\gamma$, and uses it as a floor of the internal energy in the cold cells. These legs check that it
changes nothing where it must not, and that it keeps a shock heated; the collapse it rescues is
[SV-17](../self-gravity/polytrope).

What the case catches:
- **a model that perturbs a warm flow:** where no cell is cold the floor never acts, and the Euler fields must be those
  of `euler` bitwise;
- **a floor that removes shock heat:** the entropy does not see a shock, so a scheme that resets the energy of a cold
  shocked cell from it loses the shock (the first form of D-M9-1 did, on this very case);
- **a silent input:** a switch outside $(0, 1)$;
- **a backend that disagrees:** the FNL kernels (`adam_flume_fnl_euler_de_kernels`, `sync_dual_energy_dev`) against
  the CPU ones on the same run.

## Models and numerics

| Item | Choice | Code |
|---|---|---|
| Model | `euler-de`: the Euler equations plus $\partial_t(\rho s) + \nabla\cdot(\rho s\,\mathbf u) = 0$, `nv = 6` | [models](../../numerics#dual-energy); `adam_flume_euler_de_library`, `adam_flume_euler_library_body.INC` |
| Space, time | as [Sod](./sod) (V1): WENO5 characteristic flux splitting, SSP-33, CFL 0.5; the entropy is the sixth characteristic field (eigenvalue $u_n$) | `adam_flume_cpu_euler_de_kernels.F90` |
| Synchronisation | once per step: in a cold cell ($e < \eta E$, $\eta$ = `[physics] dual_energy_eta`, default 0.05) $e = \max(e, e_s)$, $e_s = \rho s\,\rho^{\gamma-1}/(\gamma-1)$; then everywhere $\rho s = (\gamma-1)\,e\,\rho^{1-\gamma}$ | `sync_dual_energy`, [`adam_flume_cpu_object.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/cpu/adam_flume_cpu_object.F90); FNL: `sync_dual_energy_dev`, [`adam_flume_fnl_euler_de_kernels.F90`](https://github.com/szaghi/adam/blob/develop/src/app/flume/fnl/adam_flume_fnl_euler_de_kernels.F90) |

## Initial and boundary conditions

| Leg | Setup |
|---|---|
| de1 | the `sod-x` regression input ([Sod](./sod)), 50 cells on $[0, 1]$, $t = 0.2$, with `euler` and with `euler-de` |
| de2 | the planar Noh problem in a cold gas: $\rho = 1$, $u = \pm 1$, $p = 10^{-6}$ (Mach 845, $\gamma = 1.4$) colliding at $x = 0.5$, the `sod-x` grid and numerics, `extrapolation` boundaries, $t = 0.6$. Exact: shocks moving out at $(\gamma-1)/2 = 0.2$, behind them $\rho = (\gamma+1)/(\gamma-1) = 6$ and $p = \rho_0 u_0 (u_0 + 0.2) = 1.2$ |
| de3 | the de1 input with `dual_energy_eta` = 0 and 1 |
| de4 | the de2 `euler-de` runs of CPU and FNL |

## Input files

| File | Role |
|---|---|
| `src/tests/flume/verification/dualenergy/make_de.py` | `sod` and `noh` write the inputs (`--model`, `--eta`) |
| `src/tests/flume/verification/dualenergy/check.sh` | runs the legs (`--legs de1,de2,de3,de4`, default all but de4; `FLUME_EXE` selects the backend; de4 compares the work directories of both) |
| `src/tests/flume/verification/dualenergy/de_oracle.py` | the bitwise comparison (de1) and the Noh bounds (de2) |
| `make_doc_runs.sh <runs> dualenergy`, `make_doc_figures.py --only dual-energy` | the figure (2 runs, a minute) |

```bash
cd src/tests/flume/verification/dualenergy && ./check.sh                              # CPU
FLUME_EXE=$PWD/../../../../../exe/adam_flume_fnl ./check.sh                         # FNL
./check.sh --legs de4                                                                 # CPU against FNL
```

## Results

![The cold Noh shock with and without dual energy](/flume/dual-energy-noh.png)

*The cold planar Noh shock at $t = 0.6$, 50 cells. Both models put the shocks at $0.5 \pm 0.12$, the density at 6 and the
pressure at 1.2 behind them (the dip at the centre is the wall heating of the Noh problem, the same in both). In the cold
streams the pressure stays at $10^{-6}$; `euler-de` lifts the two cells at each shock foot to $1.4\cdot10^{-5}$, the
advected entropy diffused across the foot.*

| Leg | Asserted | Bound | CPU | FNL |
|---|---|---|---|---|
| de1 | `r`, `ru`, `rv`, `rw`, `rE` of `euler-de` against `euler` | bitwise | bitwise equal | bitwise equal |
| de2 | post-shock pressure (median within 0.08 of the centre), `euler` / `euler-de` | 2 % of 1.2 | $1.0\cdot10^{-4}$ / $1.2\cdot10^{-4}$ | $1.0\cdot10^{-4}$ / $1.2\cdot10^{-4}$ |
| de2 | cold-stream pressure (beyond 0.2 of the centre), largest relative error, `euler` / `euler-de` | $10^{-3}$ | $1.3\cdot10^{-8}$ / $2.6\cdot10^{-5}$ | $1.2\cdot10^{-8}$ / $2.6\cdot10^{-5}$ |
| de2 | largest density difference of the two models | $2\cdot10^{-2}$ | $6.9\cdot10^{-3}$ | $6.9\cdot10^{-3}$ |
| de3 | `dual_energy_eta` = 0, 1 | refused | refused | refused |
| de4 | per field, FNL against CPU relative to its largest magnitude, `euler-de` against the `euler` pair (`rs` against their largest) | $\le 2\times$ the `euler` pair's, or $10^{-10}$ | `ru` $1.5\cdot10^{-9}$ (`euler` $1.0\cdot10^{-9}$), `rs` $5.2\cdot10^{-10}$, the rest $\le 6.1\cdot10^{-10}$ | (the same comparison) |

Related cases: [Sod](./sod) is the warm baseline of de1; the collapse dual energy rescues is
[SV-17](../self-gravity/polytrope).

## Findings

- **Where no cell is cold, `euler-de` is `euler`.** In Sod every cell has $e \ge 0.05 E$, the floor never acts, and the
  entropy only follows the energy: the five Euler fields are bitwise those of `euler`.
- **A floor, not a reset.** The first form of the synchronisation (D-M9-1 as approved: the energy of a cold cell reset
  from its entropy) took the shock's heat away in the first steps of this case, when every cell is cold, and the run
  went non-finite at step 3. Numerical dissipation only raises the entropy, so the drained cells are the ones below
  their adiabat: lifting only those keeps every shock heated (the post-shock pressure of the two models agrees to
  $3\cdot10^{-5}$), and the energy the floor adds is the only change to the total energy.
- **The price is a little preheating at a cold shock's foot.** The advected entropy diffuses forward across the foot,
  so the floor lifts the one or two cold cells there, to $10^{-5}$ in a $10^{-6}$ stream ahead of a 1.2 shock. A switch
  at Enzo's $10^{-3}$ would lift the same cells (their $e/E$ is below it as well): it is the floor, not the switch.
- **The FNL port runs the same rule.** On a Mach 845 shock the two backends already differ by $10^{-9}$ with `euler`
  (the device arithmetic); `euler-de` differs by at most 1.5 times that, and its entropy by half of it (de4, the
  relative-to-baseline rule of [AV](../runtime-amr/accuracy) av8: an absolute $10^{-10}$ bound fails `euler` itself).
- **The switch is 0.05, not Enzo's $10^{-3}$** (D-M9-3, revised): in the #88 collapse the coldest cells sit at
  $e/E \approx 10^{-2}$ two steps before they go negative, so $10^{-3}$ never acts; 0.05 (Cholla uses 0.035) rescues the
  collapse ([SV-17](../self-gravity/polytrope)).
