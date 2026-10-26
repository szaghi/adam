# 2. The Lax problem

The Sod input becomes the Lax problem by changing two states and the final time; then the face flux changes from WENO
flux splitting to a Riemann solver. Both runs are checked against the exact solution.

## From Sod to Lax

Copy `src/tests/flume/regression/sod-x/input.ini` and change these keys:

```ini
[IO]
output_basename = lax

[initial_conditions_region_1]   ; x < 0.5
r = 0.445
u = 0.698
p = 3.528

[initial_conditions_region_2]   ; x > 0.5
r = 0.5
u = 0.0
p = 0.571

[time]
time_max = 0.13
```

Unlike Sod, the left state moves ($u = 0.698$), and its contact is stronger, so the dissipation of the face flux shows
in the error. Run and check it:

```bash
mpirun -np 2 ../../exe/adam_flume_cpu input.ini > log.txt
../../exe/.regression-venv/bin/python ../../src/tests/flume/verification/sod/sod_oracle.py input.ini .
```

```text
/lax: axis x, it 245, cells 200, L1(rho) = 1.063839e-02
```

## A Riemann solver at the faces

`scheme_space = weno` splits the flux and reconstructs the two halves. `weno-riemann` instead interpolates the face
states with WENO, solves a Riemann problem at every face, and adds a 6th-order flux correction gated by the WENO
weights. Change `scheme_space` and add the three keys of the Riemann flux:

```ini
[numerics]
scheme_space             = weno-riemann
reconstruction_variables = characteristic
riemann_solver           = hllc      ; llf, hll or hllc
flux_correction          = 6th
flux_correction_sensor   = weno
reflux                   = .true.
```

```text
/lax: axis x, it 245, cells 200, L1(rho) = 1.080670e-02
```

Both numbers are the ones of the [Lax verification page](/applications/flume/verification/euler/lax) ($1.064\cdot10^{-2}$
splitting, $1.081\cdot10^{-2}$ HLLC). Lax is the one problem where splitting beats every Riemann solver, by 1.6%; the
page shows why, and gives LLF and HLL too.

::: tip What you learned
- A Riemann problem is two region sections and a final time.
- The face flux is a `[numerics]` choice: `weno` (flux splitting) or `weno-riemann` with a solver, a correction and its
  sensor.
- The verification pages record the numbers your runs must reproduce.
:::

Previous: [1. Sod](./sod). Next: [3. A blast on a static AMR grid](./blast-amr).
