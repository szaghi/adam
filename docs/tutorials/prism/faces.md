# 2. Reflecting and absorbing faces

The pulse of the first tutorial, centred in the box, against two kinds of $x$ faces: a perfect conductor that
reflects it and a Silver–Müller face that lets it out. The energy history tells them apart.

## The changes

Start from `src/tests/prism/regression/periodic-pulse-blocks/input.ini` and change:

```ini
[initial_conditions]
r0 = 0.0             ; the packet starts at the centre, away from the x faces

[bc_x_min]
type = PEC           ; then: Silver_Muller
[bc_x_max]
type = PEC           ; then: Silver_Muller
```

The $y$ and $z$ faces stay periodic. The two faces of an axis must be both periodic or neither; the accepted kinds and
their ghost relations are on the [boundary conditions](/applications/prism/boundary-conditions) page (`Silver_Muller`
is spelled with an underscore).

## PEC: the pulse reflects

A perfect electric conductor sets $\mathbf E_t = 0$ and $B_n = 0$ on the face, by odd reflection into the ghosts. The
packet's centre reaches $+x$ at step 96; the pulse reflects and travels back. Total energy, $W_D + W_B$ from `-energy_history.dat`:

| step | 10 | 50 | 100 | 150 | 192 |
|---|---|---|---|---|---|
| $W_D + W_B$ (J) | 84628.44 | 84628.44 | 84628.43 | 84628.43 | 84628.42 |

The box is closed: the energy is conserved to $2\cdot10^{-7}$ over the run.

## Silver–Müller: the pulse leaves

The first-order Silver–Müller condition imposes the impedance of an outgoing plane wave, $|\mathbf E_t| = c_0|\mathbf B_t|$,
on the face ghosts: a wave arriving along the normal passes through.

| step | 10 | 50 | 100 | 150 | 192 |
|---|---|---|---|---|---|
| $W_D + W_B$ (J) | 84628.44 | 84628.16 | 29704.97 | 375.70 | 191.57 |

The pulse crosses the face between steps 50 and 150; what remains at step 192, $2.3\cdot10^{-3}$ of the initial energy,
is what the first-order condition does not let out. The next page replaces it with a perfectly matched layer.

::: tip What you learned
- A face kind is one key per face; `periodic` pairs the two faces of an axis.
- PEC reflects and conserves the energy; Silver–Müller transmits a normally incident wave to first order.
- The energy history measures what a boundary lets out.
:::

Previous: [1. A pulse in a periodic box](./pulse). Next: [3. A perfectly matched layer](./pml).
