# #161: the FLUME CPU channel slowdown, located (phase 1)

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
