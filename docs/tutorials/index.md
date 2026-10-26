# Tutorials

Three linear tracks. Each page builds on the previous one, every step is a run you can repeat, and every number on
the pages comes from such a run. Start with the [first run](/start/first-run) if you have not built ADAM yet.

| Track | Pages | What you do |
|---|---|---|
| **FLUME** | [1. Sod](./flume/sod) · [2. Lax](./flume/lax) · [3. A blast on a static AMR grid](./flume/blast-amr) · [4. The positivity limiter](./flume/limiter) | run compressible Euler and MHD cases from the regression inputs, change one key at a time, check against exact solutions and conservation |
| **PRISM** | [1. A pulse in a periodic box](./prism/pulse) · [2. Reflecting and absorbing faces](./prism/faces) · [3. A perfectly matched layer](./prism/pml) · [4. The same on the GPU](./prism/gpu) | run an electromagnetic pulse, compare boundary treatments through the energy history, run it on the GPU backend |
| **Build your own solver** | [1. A realm and the core layer](./build-your-own/core) · [2. The library's numerics](./build-your-own/numerics) · [3. AMR and the staged layer](./build-your-own/amr) · [4. The GPU twin](./build-your-own/gpu) | read RIPPLE, the scalar advection app written for this track, and the realm contract it implements, layer by layer |

The runs on these pages were made on a WSL2 workstation (gfortran 16.2 and Open MPI 5.0.11 for the CPU, nvfortran 26.9
for the GPU), which is for correctness only: no page quotes a timing.

The [cookbook](/cookbook/) answers single tasks; the application pages ([FLUME](/applications/flume/),
[PRISM](/applications/prism/)) and the [references](/reference/library-input) give the full picture.
