# 3. A blast on a static AMR grid

A 2-D blast on a grid with a refined patch: the AMR keys, what the refined grid looks like in the output, and the check
that matters across 2:1 faces, conservation.

## The case

`src/tests/flume/regression/blast-amr-limiter/input.ini` is the Balsara–Spicer magnetised blast: $\rho = 1$, $p = 0.1$
outside, $p = 1000$ inside $r < 0.1$ around $(0.5, 0.5)$, $|\mathbf B| = 100/\sqrt{4\pi}$ at 45° (ambient
$\beta = 2.5\cdot10^{-4}$), ideal MHD with EGLM cleaning, periodic in $x$ and $y$. It is a hard problem on purpose:
the low $\beta$ drives the pressure negative unless the positivity limiter is on, which is the subject of
[the next page](./limiter). Here it is the grid that matters.

## The AMR keys

(Excerpt: the committed input also sets the marker's `delta_type`, `delta_fine`, `delta_coarse`, `box_zmin` and
`box_zmax`, which the parser requires.)

```ini
[grid]
ni = 16
nj = 16
nk = 4             ; z is a null axis (null_z = .true.)

[amr]
max_level      = 3
ratio          = 8  ; octree
iu_ref_levels  = 2  ; two uniform levels: 64 blocks, 64 x 64 cells
markers_number = 1  ; one marker refines further
frequency      = 0  ; no runtime regridding

[amr_marker_1]
mode         = 1               ; geometric marker
geo_type     = primitive-box
box_xmin     = 0.55
box_ymin     = 0.55
box_xmax     = 0.7
box_ymax     = 0.7
target_level = 3               ; refine the blocks in the box up to level 3

[initial_conditions]
amr_iterations = 1             ; set the IC, refine, set the IC again on the new grid
```

A block is refined when its centroid lies in the box; the box selects the base block $[0.5, 0.75]^2$, so two 2:1
faces cross the blast centre. Every refined axis needs an even block cell count of at least `2 ngc`
([AMR cookbook](/cookbook/amr)).

## Run it

```bash
mpirun -np 2 ../../exe/adam_flume_cpu input.ini > log.txt
```

The case stops after `it_max = 100` steps, at $t = 2.604\cdot10^{-3}$. Counting the blocks of the field files, step 0
and step 100 both hold 92 blocks with two spacings, $\Delta x = 1/64$ and $1/128$: the 4 blocks of the box column
(the null $z$ axis keeps 4 copies) became 32.

In ParaView, colour by `r` and switch the representation to *Surface With Edges*: the refined patch shows as blocks of
half the cell size.

## The check: conservation across the 2:1 faces

`blast-amr-limiter-conservation_history.dat`, first and last rows, the mass column:

```text
VARIABLES="it" "time" "int_r" ...
0   0.0000000000000000E+000 1.0000000000000000E+000
100 2.6041666666666683E-003 9.9999999999999967E-001
```

On a periodic box the mass is conserved exactly by the scheme, and it is here to $3\cdot10^{-16}$ although the blast
crosses the 2:1 faces: at the end of each step the reflux replaces the coarse face flux with the sum of the fine ones,
so both sides exchange the same flux. The [AMR conservation page](/applications/flume/verification/euler/amr-conservation)
checks the same property on a dedicated case, and the [positivity limiter page](/applications/flume/verification/mhd/positivity-limiter)
(PV-5) runs the same blast to $t = 0.01$ on a different refined region, $[0.25, 0.75]^2$.

::: tip What you learned
- Static refinement is `max_level`, a marker section and `amr_iterations`; `frequency = 0` keeps the grid fixed.
- The output holds every block at its own spacing; counting blocks and spacings is the quickest check of the grid.
- Across 2:1 faces the reflux keeps the scheme conservative; the conservation history shows it.
:::

Previous: [2. Lax](./lax). Next: [4. The positivity limiter](./limiter).
