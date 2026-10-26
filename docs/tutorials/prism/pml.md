# 3. A perfectly matched layer

The Silver–Müller faces of the previous page let the pulse out to first order. A perfectly matched layer (PML) absorbs
it inside a layer next to the face, with a damping that grows towards the boundary. This page adds one on the $x$
faces and runs it to the same physical time.

## The changes

Start from the Silver–Müller input of the previous page (pulse at $r_0 = 0$, `Silver_Muller` on the $x$ faces) and
add a `[PML]` section:

```ini
[PML]
PML_type       = classic
width          = 0.04          ; metres: about 10 cells of 0.0042 at each x face
x_minus_layer  = .true.
x_plus_layer   = .true.
y_minus_layer  = .false.
y_plus_layer   = .false.
z_minus_layer  = .false.
z_plus_layer   = .false.
gamma_max      = 7.8E10        ; peak damping rate (1/s)
gamma_exponent = 2.0           ; quadratic profile

[time]
it_max   = -1                  ; run to time_max, not a step count
time_max = 6.6712814328639916E-10
```

Each `PML_type` has its own parameters (`classic` takes `gamma_max` and `gamma_exponent`; the list and the profiles are
on the [PML page](/applications/prism/pml)). The PML needs centred FD and an SSP Runge–Kutta scheme, which this input
has.

**Why `time_max`.** The PML caps the time step at $C_\gamma/\gamma_{\max}$ (`C_gamma` = 0.1 by default). The log shows

```text
[mpi-00000]      C_gamma: +0.10000000000000001E+000
[mpi-00000]time step:    0.12820512820512821E-011
```

$1.28\cdot10^{-12}$ s, against $3.47\cdot10^{-12}$ s without the layer: 192 steps would reach only $t = 2.5\cdot10^{-10}$
s. The `time_max` above is the time the 192 steps of the previous pages reached, and the run takes 521 steps.

## The result

Total energy $W_D + W_B$ from `-energy_history.dat`:

| step | 10 | 110 | 210 | 310 | 410 | 521 |
|---|---|---|---|---|---|---|
| $t$ (s) | $1.28\cdot10^{-11}$ | $1.41\cdot10^{-10}$ | $2.69\cdot10^{-10}$ | $3.97\cdot10^{-10}$ | $5.26\cdot10^{-10}$ | $6.67\cdot10^{-10}$ |
| $W_D + W_B$ (J) | 84628.44 | 84491.46 | 35867.34 | 46.86 | 0.256 | $9.9\cdot10^{-3}$ |

At the same final time the PML leaves $1.2\cdot10^{-7}$ of the energy, against $2.3\cdot10^{-3}$ with the Silver–Müller
faces alone: four orders of magnitude less.

**One caveat.** The PML profile is computed per rank, from the blocks the rank holds: a layer deeper than one block,
whose blocks sit on different ranks, gives results that depend on the decomposition
([#108](https://github.com/szaghi/adam/issues/108)). Here the layer (0.04) lies inside the outermost layer of blocks
(0.05 wide), so every rank that holds part of it sees its whole width. Keep the layer within one block until #108 is
fixed.

::: tip What you learned
- A PML is a `[PML]` section: a type, a width, the faces, the type's parameters.
- The PML caps the time step; compare runs at equal physical time, not equal step count.
- The energy history measures the absorption: here $10^{-7}$ against $10^{-3}$ for the first-order condition.
:::

Previous: [2. Reflecting and absorbing faces](./faces). Next: [4. The same on the GPU](./gpu).
