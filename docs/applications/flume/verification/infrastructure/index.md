# Infrastructure

The legs of this group verify the machinery under every FLUME case rather than a flow: what the ghost cells hold on
every fill path, whether the scheme carries absolute constants that break the scaling of the equations, and whether a
rank that owns no block changes the result. Each oracle is an exact identity (a linear field reproduced, a rescaled run
equal bit for bit, np 2 equal to np 1 bit for bit), so the bounds sit at round-off or at zero.

| Case | ID | What it verifies |
|---|---|---|
| [Ghost probe](./ghost-probe) | GP | every face and edge ghost, on every fill path (copy, MPI, 2:1, inter-realm seams, physical boundaries), holds the value a stencil must read; a poison exposes reads before the fill |
| [Scaling covariance and the reference layer](./scaling-covariance) | NV-0, NV-4, NV-5, NV-6, NV-8, NV-9 | rescaled inputs reproduce the base run bit for bit; dimensional input, restart and output through `[reference]` |
| [Ranks without blocks](./ranks-without-blocks) | DV-1 | a rank that owns no block runs every model and option, and np 2 equals np 1 bit for bit |
| [Non-finite state guard](./nonfinite) | — (`nonfinite`) | a NaN planted in the Sod state stops the run at step 1 on the guard's message, on every rank, with the first non-finite value located |

The two figures of the group are regenerated from kept runs:

```bash
cd src/tests/flume/verification
(cd ghosts && ./check.sh --cases step --np 2 --keep)      # the "after" run; the "before" run is ghosts/work-before-p0-np2-step
./make_doc_runs.sh <dir> ranks
python3 make_doc_figures.py --out ../../../../docs/public/flume --runs <dir> --only ghosts ranks
```
