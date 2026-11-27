# #161: the FLUME CPU channel slowdown, located and fixed

Issue [#161](https://github.com/szaghi/adam/issues/161), steps 1.1-1.3 and G1. medooza prod-gn CPU nodes (2 x 32 cores,
gfortran 16.2.0 + Open MPI 5.0.11, `iac_gnu`), 64 ranks, the 3a windows of #160 (the uniform channel, 6.29 M cells,
512 blocks of 32 x 24 x 16, restarted from the dumps of `i160/work`, repartitioned for 64 ranks). Every comparison is
an A/B on one node: the executables alternated run by run, the order rotated each repetition, 3 pairs of (10, 110)-step
runs; the step cost is (loop(110) - loop(10)) / 100 (`logs/steps.py`), the median over the pairs. The "wait" column is
the same difference of the step-end `MPI_ALLREDUCE` wait (the forest's "regrid" timer since b8641a59), per step.

## The probes

Only six commits between the endpoints (b516b303, the #133 night's last, and 5b200ffa) touch the CPU path's sources
(`src` outside `src/tests`, `src/app/prism` and `src/lib/fnl`): 0f751be1, b8641a59, 4a029436, 44cfdfa8, 82c70653; the
later ones (7bb5b134, f5826118, 144add34) touch only the FNL library and app. Each was built from a `git archive` of its
commit (`fobos`, `fobos.d`, `src/app`, `src/lib`, `build.sh`), with the dependencies of `i160/deps.tgz` (no dependency
commit since 2026-11-22), by `build.sh cpu --varset iac_gnu` (`logs/i161-build.slurm`, job 1370554). All of them were
timed instead of bisected: the same answer in one round.

| Label | Commit | Link flags | sha256 (16) |
|---|---|---|---|
| A | the #133 matrix's CPU executable | `-O2` | `5b9bd6873c1468a6` |
| b516 | b516b303 | `-O2` | `5b9bd6873c1468a6` |
| b516ftz | b516b303, `template-gnu` lflags patched as 2a67b16c | `-O2 -mdaz-ftz` | `44dd6ea3b5ab7745` |
| 0f75ftz | 0f751be1, idem | `-O2 -mdaz-ftz` | `f1fca2ad9041452c` |
| b864 | b8641a59 | `-O2 -mdaz-ftz` | `51c6c205131ce344` |
| 4a02 | 4a029436 | `-O2 -mdaz-ftz` | `2fdb1d1a6d64677b` |
| 44cf | 44cfdfa8 | `-O2 -mdaz-ftz` | `470163253e324fa7` |
| 82c7 | 82c70653 (= 5b200ffa on the CPU path) | `-O2 -mdaz-ftz` | `b408a22080269248` |

**Provenance (1.1).** b516b303 rebuilds to the #133 matrix's CPU executable byte for byte (`5b9bd687`): the "src/
modified" tree of its `toolchain.md` held b516b303's sources. 82c7 rebuilds to the executable #160 measured as B
(`b408a220`): the builds are reproducible.

## Step costs (s per step, 3a, 64 ranks)

| Job (node) | Executable | Pairs | Median | vs reference | Wait (s/step) |
|---|---|---|---|---|---|
| t11, 1370555 (gnode20) | A | 3.501 / 3.107 / 3.087 | 3.107 | reference | 0.053 / 0.003 / 0.002 |
| | b516 | 3.106 / 3.202 / 3.100 | 3.106 | -0.0 % | 0.004 / 0.004 / 0.004 |
| | 82c7 | 3.758 / 3.729 / 3.694 | 3.729 | **+20.0 %** | 0.131 / 0.143 / 0.155 |
| t12a, 1370556 (gnode03) | b516ftz | 3.094 / 3.219 / 3.123 | 3.123 | reference | 0.003 / 0.002 / 0.003 |
| | 0f75ftz | 3.128 / 3.253 / 3.334 | 3.253 | +4.2 % | 0.006 / 0.003 / 0.023 |
| | b864 | 3.660 / 3.412 / 3.724 | 3.660 | **+17.2 %** (+12.5 % over 0f75ftz) | 0.146 / 0.079 / 0.171 |
| t12b, 1370557 (gnode04) | b516ftz | 3.099 / 3.246 / 3.208 | 3.208 | reference | 0.001 / 0.006 / 0.025 |
| | 4a02 | 3.424 / 3.644 / 3.940 | 3.644 | +13.6 % | 0.075 / 0.130 / 0.187 |
| | 44cf | 3.664 / 3.442 / 3.601 | 3.601 | +12.2 % | 0.122 / 0.068 / 0.115 |
| | 82c7 | 3.539 / 3.483 / 3.303 | 3.483 | +8.6 % | 0.075 / 0.085 / 0.039 |

- **1.1 (endpoints):** b516b303 equals A (the same binary, medians 0.0 % apart); 5b200ffa's CPU code is 20 % above.
  The #160 job (1370519) measured 8.3 % on the same node: its A pairs (3.185-3.263) ran slower than today's.
- **1.2 (the commit): b8641a59** (`perf(field)`: the ghost exchange drops its two `MPI_Barrier`, #156 L1). Every build
  without the barriers is 9-20 % above the reference, every build with them within noise of it. The criterion "the
  parent within 2 %" is missed on the median (0f75ftz +4.2 %) through one disturbed pair (3.334, with the step-end
  wait of a disturbed run); the pair ratios are 1.011 / 1.011 / 1.068, and 0f751be1 changes only `refine_uniform`,
  which runs at initialization before the restart replaces the layout: the time loop is the same code in both.
- **The signature:** without the barriers the wait at the step-end allreduce is 0.04-0.17 s per step, with them
  0.001-0.006. Out of noise in every pair.

## Where the time goes (1.3)

Two variants of 82c7 with per-rank timers in `update_ghost_mpi` (`logs/instrument.py`: pack, post, `MPI_WAITALL`,
barriers, unpack; min / mean / max over the ranks at the end of the run): i0 as is, i1 with b8641a59's two barriers
restored. Job 1370596 (`logs/i161-var.slurm`), the same A/B against b516ftz.

| 110-step run, mean over the ranks (s) | i0 (no barrier) | i1 (two barriers) |
|---|---|---|
| step cost, pairs (median) | 3.502 / 3.490 / 3.173 (3.490) | 3.177 / 3.432 / 3.256 (3.256) |
| b516ftz in the same job | 3.298 / 3.148 / 3.388 (3.298) | |
| `MPI_WAITALL` (+ barriers), per run | 43.9 / 52.3 / 15.4 | 8.5 / 31.3 / 15.0 |
| the ranks' own work (steps - waitall - barriers) | 334.8 / 323.2 / 337.2 | 340.0 / 340.6 / 340.1 |
| pack + unpack | 5.9 / 6.2 / 6.6 | 8.3 / 8.1 / 8.3 |

The ranks do not compute more without the barriers: the time is spent waiting in `MPI_WAITALL`, 44-52 s per run on the
typical runs against 8-15 s with the barriers (one i0 run was fast, one i1 run disturbed). It is not load imbalance,
which the barriers would pay at each exchange. The exchange packs, posts, waits and unpacks with no overlap; its
messages (a block face, hundreds of kB per neighbour) are far above the shared-memory transport's eager limit (Open MPI
5.0.11 on these nodes: `btl_sm_eager_limit` 4096, `max_send_size` 32768, single copy by CMA or XPMEM), so each goes by
rendezvous, which advances only while both ranks are inside MPI. With the barrier every rank enters the transfer
together; without it a rank waits for peers still computing, and the waits chain.

## The step does not show it (G1)

The 1a windows (the 3-realm step, 4.13 M cells, 64 ranks; the dumps of `i160/work`), 2 pairs alternated on one node
(job 1370620, `logs/i161-ab3.slurm`, main input `step.ini`):

| Executable | Pairs | Median | vs 0f75ftz | Wait (s/step) |
|---|---|---|---|---|
| 0f75ftz | 1.245 / 1.268 | 1.256 | reference | 0.004 / 0.004 |
| b864 | 1.269 / 1.272 | 1.271 | +1.1 % | 0.004 / 0.004 |
| 82c7 | 1.248 / 1.307 | 1.278 | +1.7 % | 0.004 / 0.004 |

b8641a59 leaves the step within 2 % of its parent, and the step-end wait does not move (0.004 s per step in every
build): the barriers do not matter there. A likely reason (an estimate, not measured): the step's blocks are 2-D
(768-12288 cells, faces of tens of cells), so its ghost messages are one to two orders of magnitude smaller than the
channel's (faces of 32 x 24 cells, 3-D), close to or below the eager limit, and the exchange does not depend on both
ranks being in MPI at once. (Job 1370601 ran the same 1a A/B with the wrong main input, `input.ini`: every run stopped
at the input check, rc 1, an error of the job script; its 2a part is valid.)

The #133 README held the same signal and read it as spread: "Skipping the barriers on the CPU brings nothing (one node:
1150 against 974 ms without the body, 1367 against 1266 with it; single runs, inside the CPU node's spread or worse)".

## The fix (D1, 2.1)

D1 = A, on the CPU path only: `update_ghost_mpi` (`src/lib/common/adam_field_object.F90`, the CPU exchange) calls one
`MPI_Barrier` before its posts, so every rank enters the transfer together; the FNL exchange `update_ghost_mpi_gpu`
keeps none (its gain, 0.4-1.3 %, was measured on GPUs). Two candidates, built from HEAD e8ac01f5 (= 5b200ffa on the CPU
path; job 1370599, `logs/i161-fix.slurm`): fixA (one barrier, before the posts) and fixB (b8641a59's CPU hunk
reverted: two barriers, after the wait and after the unpack).

| Case, job | Executable | Pairs | Median | vs reference | Wait (s/step) |
|---|---|---|---|---|---|
| 3a, 1370599 | b516ftz | 3.132 / 3.154 / 3.226 | 3.154 | reference | 0.004 / 0.004 / 0.023 |
| | **fixA** | 3.211 / 3.155 / 3.142 | **3.155** | **+0.0 %** | 0.024 / 0.012 / 0.006 |
| | fixB | 3.248 / 3.251 / 3.143 | 3.248 | +3.0 % | 0.004 / 0.026 / 0.003 |
| | 82c7 (HEAD) | 3.518 / 3.604 / 3.899 | 3.604 | +14.3 % | 0.101 / 0.100 / 0.189 |
| 1a, 1370620 | 0f75ftz | 1.245 / 1.268 | 1.256 | reference | 0.004 |
| | fixA | 1.257 / 1.284 | 1.270 | +1.1 % | 0.004 |
| | 82c7 | 1.248 / 1.307 | 1.278 | +1.7 % | 0.004 |
| 2a, 1370601 (the #133 matrix's 2a dumps) | 82c7 | 8.068 / 7.792 | 7.930 | reference | 0.443 / 0.211 |
| | fixA | 7.701 / 7.669 | 7.685 | -3.1 % | 0.042 / 0.036 |

fixA restores the #133 cost of the channel (3.155 against 3.154) and is no slower than fixB: it is the fix. The sphere
was slowed too (the step-end wait 0.2-0.4 s per step at HEAD) and gains 3.1 %; the step does not move.

**Gate (job 1370600, `logs/i161-gate-cpu.slurm`), base = HEAD, fix = HEAD + fixA, clean builds of flume-cpu-gnu(-fft),
prism-cpu-gnu, ripple-cpu-gnu and the RIPPLE criterion test:** FLUME CPU suite FAIL 0 SKIP 0; PRISM CPU suite fails only
rmf-hyperbolic (the known pre-existing mismatch), SKIP 10; RIPPLE oracles PASSED; DV-1 (empty ranks: the exchange is
collective again) PASSED at np 2 and 4 (10 legs each); the A/B of every FLUME and PRISM case base against fix: 91
bitwise, pic-periodic and pic-slab-pulse not completing on the CPU with either executable (as in the #157 night's A/B).
The comment in the fix was reworded after these jobs; the code they built is the committed one.

## The channel's ratios after the fix (2.2)

| Machine, job | Runs | 3a CPU node (s per step) | Spread | GPU node (s) | GPU node / CPU node |
|---|---|---|---|---|---|
| medooza, 1370614 (gnode16), `medooza-2.2/` | `run.sh windows --cases 3a`, 64 ranks, the fix (sha256 `ece2fda4`) | **3.201** | 2.3 % | 0.2219 (C4) | **14.4** (13.9 with 0.230, before C4) |
| SpaceHPC, 185266 (gpu-015), `spacehpc/` | 4 windows alternated, 96 ranks: fix1, head1, fix2, head2 | **1.478 / 1.466** (fix), 1.532 / 1.466 (HEAD) | 4.5 / 1.2 %, 7.1 / 6.0 % | 0.0609 (C4) | **24.2** (23.1 with 0.0638, before C4) |

- medooza: 1.0 % above the #133 build's 3.17 (gate: within 3 %); the step-end wait 0.0031 s per step (0.1288 in #160's
  HEAD run, 3.59 s). The channel's factor 16.2 of #160 becomes 14.4.
- SpaceHPC (the EPYC node): HEAD runs as fast as the fix in the same job (step-end wait 0.043-0.051 against
  0.030-0.034): the barriers' removal did not cost the EPYC node what it cost medooza's Xeon nodes, and the fix is
  neutral there. The 3a cost of #160 (1.564 on gpu-005) against 1.47 here is the spread between nodes, not the fix;
  the ratio 24.2 (25.7 in #160) is this job's.
- The tree of the SpaceHPC build is `adam-i160` (07a414ba, HEAD on the CPU path) with the fixed field object, hence
  "src/ modified" in `spacehpc/toolchain.md` (executable sha256 `25c8c2aa`).
