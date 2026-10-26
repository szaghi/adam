# 4. The positivity limiter

The blast of the previous page without the refined patch: run it without the limiter and see it fail, then with it
and see it reach the final time.

## The case

`src/tests/flume/regression/blast-limiter/input.ini` is the same Balsara–Spicer blast on the uniform $64\times64$ grid
(a quadtree: `ratio = 4`, `nk = 1`). Two changes to the copied input run it to the end instead of the regression's 100
steps, and select the limiter:

```ini
[numerics]
positivity_limiter = none     ; then: cell

[time]
it_max   = -1                 ; was 100: run to time_max
time_max = 0.01
```

The `[mhd]` floors are off (`rho_floor = 0.0`, `p_floor = 0.0`), so nothing hides a negative pressure.

## Without the limiter

```bash
mpirun -np 2 ../../exe/adam_flume_cpu input.ini > log.txt
```

```text
[mpi-00000]error stop : +4 cells with a non-positive density or pressure at step +3 (min rho +0.72590331193952229E+000, min p -0.15378822578668403E-002); the [mhd] floors rho_floor, p_floor are disabled
```

At $\beta = 2.5\cdot10^{-4}$ the pressure is a tiny difference of large energies, and the high-order update makes it
negative in the third step. FLUME stops rather than continue with a non-physical state.

## With the limiter

`positivity_limiter = cell` blends every face flux of a cell with the first-order Lax–Friedrichs flux, just enough to
keep the cell's density and pressure above a fraction (`POSITIVITY_LIMITER_KAPPA` = 0.1) of their first-order update.
The same input with `cell`:

```text
+193 +0.99908507877250665E-002 +8 ...       (the last row of blast-limiter-residuals.dat)
```

It reaches $t = 0.01$ in 193 steps (the row shows the time at the start of the last step). The limiter changes only
the cells where the high-order update would leave the admissible set; elsewhere the scheme is unchanged. The
[positivity limiter page](/applications/flume/verification/mhd/positivity-limiter) runs this blast at $128^2$ with
figures, its 3-D and AMR versions, and the near-vacuum tests.

::: tip What you learned
- FLUME stops on a non-positive density or pressure unless floors are set; it never continues silently.
- `positivity_limiter = cell` keeps the update admissible by blending with the first-order flux, cell by cell.
- The regression cases are small, complete inputs: copy one, change two keys, and you have a new experiment.
:::

Previous: [3. A blast on a static AMR grid](./blast-amr). This is the end of the FLUME track; the
[build-your-own-solver](/tutorials/build-your-own/core) track shows how such an application is built.
