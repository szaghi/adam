# rmf-2realm-stagesync regression case

Two-realm x-split rmf with **β (stage-coincident)** inter-realm seam
coupling (PRD #18). The geometry is byte-identical to `../rmf-2realm/`;
the only difference is the manifest's `coupling_cadence` flag on `face_1`:

```ini
[forest.topology.face_1]
...
coupling_cadence = stage_coincident    ; β: per-substage seam fill (issue #18)
```

Under β the forest fills inter-realm seam ghost cells **at every RK
substage** inside the K loop (Phase 2 of `forest_object%evolve_one_step`)
rather than once per global timestep. This gives the seam the same
temporal order as the per-realm interior, restoring bit-equivalence to a
monolithic single-realm run on the union grid — the strong oracle this
case enforces.

## Files

```
rmf-2realm-stagesync/
   input.ini       # forest manifest, with coupling_cadence = stage_coincident on face_1
   realm_1.ini     # complete PRISM INI for realm 1 (x ∈ [-0.06, 0],   ni=8, nj=16, nk=16)
   realm_2.ini     # complete PRISM INI for realm 2 (x ∈ [ 0,  0.06], ni=8, nj=16, nk=16)
   golden/
      cpu/         # CPU-backend β golden (digest.txt + per-realm residuals.dat)
      fnl/         # FNL-backend β golden (digest.txt + per-realm residuals.dat)
```

The realm INIs are byte-identical to `../rmf-2realm/realm_{1,2}.ini`
except for the header comments — same `runge-kutta-ssp-54` (K=5), same
physics, same spatial operator. Only the manifest distinguishes this
case from `../rmf-2realm/`.

## β admissibility contract

β is admissible on a seam iff both endpoint realms agree on:

1. ODE solver scheme: `numerics%scheme_time` AND `rk%scheme`
2. Stage count: `rk%nrk`
3. Physics layout: `physics%nv` (and the variable-index assignments)

Spatial operator (`numerics%scheme_space`) and grid resolution MAY
differ — that is precisely the use case β exists to serve (accuracy-
driven spatial decomposition: e.g., WENO near shocks, FD-centered in
smooth regions, with same time-integrator throughout).

Enforced at forest init by `forest_object%check_beta_admissibility`
(`src/lib/common/adam_forest_object.F90`). Any disagreement triggers
`mpih%error_stop` with a precise diagnostic naming the offending
face_pair and the specific descriptor field that mismatched. There is no
silent downgrade to α: if the user asked for β and the manifest is
inadmissible, the user wants to know loudly.

## Cross-config equivalence oracle (the load-bearing invariant)

Under same-K + same-physics + same-ODE-solver (the β admissibility
conditions), the multi-realm digest MUST match the monolithic
single-realm `rmf` digest within `rtol=1e-06, atol=1e-3` for every
FIELD variable.

The 9 per-block-metadata mismatches (`dxdydz`, `origin`,
`time_iteration`) are downgraded to `SKIP_METADATA` by `digest.py
compare` automatically — these scale with block count and the 2-realm
config has 2× the blocks of the 1-realm config.

Enforced continuously by `run.sh` and `run-fnl-local.sh`: after the
normal own-golden comparison, the harness runs a SECOND `digest.py
compare` against `rmf/golden/<backend>/digest.txt`. Mismatch fails the
case with a β-specific diagnostic; success logs an explicit oracle
confirmation:

```
>> [rmf-2realm-stagesync/cpu] β cross-config oracle: matches single-realm rmf golden
PASS [rmf-2realm-stagesync/cpu]
```

This is the canonical correctness check for β. If a future refactor to
the orchestrator, the seam-fill TBP, or the per-realm integrator path
breaks bit-equivalence with single-realm rmf, this oracle catches it
on every CI run.

## Why α gives this up

PRD #16 deliberately changed the seam coupling cadence to end-of-step
to admit asymmetric per-realm K — a flexibility that requires giving up
exactly this bit-equivalence (peer ghosts mid-step become first-order
stale, Berger-Oliger 1984; AMReX `FillCoarsePatch` convention). β
recovers the equivalence by paying the K-1 extra inter-realm exchanges
per step under the homogeneous-realm regime, and is opt-in per seam in
the manifest. α remains the default-reliable path for heterogeneous-K
forests.

## Capturing the goldens

Both backends are captured on host `adam` (gfortran 15.2.0 + OpenMPI
5.0.10 for CPU; nvfortran 26.1 + nvhpc OpenMPI for FNL, RTX 4070 ×2).
To recapture:

```bash
# CPU
REGRESSION_RUN_GOLDENLESS=1 ./run.sh cpu
cp rmf-2realm-stagesync/work-cpu/digest.txt                   rmf-2realm-stagesync/golden/cpu/
cp rmf-2realm-stagesync/work-cpu/rmf_2realm_r*-residuals.dat  rmf-2realm-stagesync/golden/cpu/

# FNL
REGRESSION_RUN_GOLDENLESS=1 ./run-fnl-local.sh
cp rmf-2realm-stagesync/work-fnl/digest.txt                   rmf-2realm-stagesync/golden/fnl/
cp rmf-2realm-stagesync/work-fnl/rmf_2realm_r*-residuals.dat  rmf-2realm-stagesync/golden/fnl/
```

The cross-config β oracle MUST pass at capture time (`>> β cross-config
oracle: matches single-realm rmf golden`); if it fails, β has a bug —
debug rather than committing the captured artefacts.

## Sibling cases

See `../rmf-2realm/README.md` for the full table of sibling configurations:

- `../rmf-2realm/` — α, symmetric K (own golden only)
- `../rmf-2realm-asymK/` — α, asymmetric K (own golden only)
- `./` — β, symmetric K (own golden + cross-config oracle vs single-realm)
