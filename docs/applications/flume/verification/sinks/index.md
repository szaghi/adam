# Sink particles

Sink particles in FLUME: point masses that act through self-gravity only ([models](../../models#sink-particles),
[numerics](../../numerics#sink-particles), input [`[sinks]`](../../input#sinks-optional-sink-particles)). The case
verifies them at rest and on a Kepler binary, their creation, accretion at Shu's rate and mergers, the input contract,
restart, CPU against FNL and runtime AMR. They need `[gravity] self_gravity = .true.`, so every case runs on the
self-gravity executables (`exe/adam_flume_cpu_fft`, `exe/adam_flume_fnl_fft`) in the box of
[the polytrope](../self-gravity/polytrope).

| Case | Type | What it verifies |
|---|---|---|
| [Sink particles](./sink-particles) | P+C | point masses through the Poisson source: rest, Kepler binary, the input contract, restart, creation, accretion at Shu's rate, mergers, CPU against FNL, runtime AMR |

The figure is regenerated with `make_doc_runs.sh <dir> sinks` (needs `exe/adam_flume_cpu_fft`) and
`make_doc_figures.py --only sink-particles` (`src/tests/flume/verification/`).

## How it was built

Sink particles are milestone M12 ([#94](https://github.com/szaghi/adam/issues/94)). The case is the leg group SI
(si1 to si11). The [verification catalog](../) places the group in implementation order.
