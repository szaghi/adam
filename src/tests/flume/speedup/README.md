# FLUME speedup matrix

**Status: being rebuilt** ([issue #133](https://github.com/szaghi/adam/issues/133), the plan's only source of truth:
every step is ticked there with its commit, every phase ends with a comment carrying its results). The previous
matrix (16 shock-over-a-body cases, `matrix.conf`, `run.sh`, `analyze.py`) was designed before the
`profile-couette3d` campaign, never ran, and was removed in phase 0. This page is a skeleton until phase 3 (the
harness) and phase 8 (the results).

## What it will measure

Three **legs** (flows), each on a uniform grid and on AMR grids refined to the same finest spacing, on every
**backend** (`cpu`, `fnl-nvf`, `fnl-amd`):

1. **2-D**: the Woodward-Colella Mach 3 forward-facing step at WENO-9, a 3-realm forest: (a) uniform fine grid,
   (b) the runtime grid frozen, (c) runtime AMR with a density-gradient marker;
2. **3-D**: shock diffraction over a sphere, immersed boundary, octree: (a) uniform, (b) runtime AMR with a
   density-gradient marker, (c) the runtime grid frozen;
3. **3-D viscous, turbulent**: the periodic Coleman-Kim-Moser channel: (a) uniform, (b) runtime AMR with a
   turbulence-resolution marker `r = h/η`, (c) the runtime grid frozen.

The questions: how much faster a GPU node is than a CPU node on each flow; how much runtime AMR buys at equal finest
resolution, and where the AMR overhead goes (the AMR grid itself against the regrids). The case table, the
measurement protocol (warm-up, settle, interleaved 10- and 110-step windows), the ratios and the hypotheses H1-H4 are
in the issue.

## Decisions

Owner, 2026-11-20:

- **D1, equal finest h.** Every AMR case refines to the uniform case's level with the same block size: same solution
  where the AMR grid is finest, same global dt. The AMR speedup is measured, not estimated, and the per-step ratio is
  the time-to-solution ratio.
- **D2, a fixed grid is the runtime grid frozen.** The frozen case restarts from the dump that starts the runtime
  case's window with `[amr] frequency = 0`: same grid, same state.
- **D3, one warm-up per leg.** Restart dumps are per rank (`<basename>-proc<rank>.fbd`); a repartition tool (phase 1)
  makes a dump written on N ranks readable on M. Each leg warms up once, on the fastest backend, and every backend and
  rank count starts its windows from the same bits.
- **D4, the channel marker is the local Kolmogorov resolution** `r = h/η` (phase 2).
- **D5, words.** *Leg* = a flow of the matrix (as `--legs` in the verification scripts); *case* = one grid of a leg
  (1a, 1b, ...); *backend* = one build made by `build.sh`. `profile-couette3d` keeps calling builds "legs".

Owner, 2026-11-22 (the re-analysis against `0bcdf870`, comment on #133):

1. **Leg 1 stays 2-D** (D1 holds on every backend); a **3-D demonstration step at h = 1/320 on 4 A30** is added for
   the figure of #142.
2. The paper commit is tagged after the step and #142 Part A.
3. The immersed-boundary rectangle distance bug goes to its own issue.

Revised order: the step forest redesign first (aligned blocks, `auto` seams, WENO-9, the feasibility gate widened),
then phase 1 (repartition, with the `.tnd` rewrite), phase 2 (a Jeans-style marker, AMR mode 7), the harness and the
legs.

## What is here now

| Path | What |
|---|---|
| `build.sh` | builds one backend from a clean tree into `exe/<backend>/`, records the toolchain in `builds/<backend>/` |
| `ini_set.py` | rewrites `section.key=value` options of an input in place (stdlib) |
| `templates/sphere-legacy.ini` | leg 2's starting point, the only 3-D sphere input FLUME had; **stale, does not run**, regenerated in phase 3 |
| `jobs/medooza-{cpu,fnl}.slurm` | placeholders until phase 3 (they keep the shape of a medooza job and stop) |
| `profile-couette3d/` | the profiling case (outside the matrix; shares the executables of `build.sh`) |

Medooza: batch only, at most one CPU and one FNL job at a time, runs in `/ifs/hpc/home/szaghi/adam-tests`, never in
the repository. WSL: dry runs and correctness only, never timings.
