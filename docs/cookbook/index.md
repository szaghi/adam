# Cookbook

Recipes: a problem, the input or code that solves it (taken from a committed input or source), what the run
produces, and where the reference is.

| Cookbook | Recipes |
|----------|---------|
| [Input and output](/cookbook/io) | the files a run writes, ParaView (and the WSL file-locking fix), reading checkpoints in Python, restart, slices |
| [Adaptive mesh refinement](/cookbook/amr) | which marker in which application, a static refined region, FLUME runtime regridding, Jeans refinement, pitfalls |
| [GPU (FNL)](/cookbook/gpu) | building FNL modes, checked device allocation and the block budget, the `buf=` rule, OpenACC/OpenMP kernel twins, WSL2 settings |
| [MPI and clusters](/cookbook/clusters) | running on several ranks, ranks without blocks, a Slurm job, GPU binding, workstation against cluster settings |
| [Regression testing](/cookbook/regression) | running the suites, what a golden is, adding a case, the lint, the refusals, verification legs |
| [Forest](/cookbook/forest) | the driver pattern, multi-realm forests (manifests, seams, cadence, 2:1 coupling), restart, forests on FNL |
