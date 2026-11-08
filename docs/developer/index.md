# Developer guide

This section is for people who extend ADAM: write a new application on the library, add a kernel or a boundary
condition to FLUME, or change the library itself. Running the solvers needs none of it; start from
[Concepts](/concepts/architecture) and the application pages instead.

| Page | Content |
|------|---------|
| [Library objects](./library) | The common-library objects, the realm/forest contract an application implements, the module summary |
| [FNL library (GPU)](./fnl) | The device layer: memory layout, FUNDAL, wrapper objects and kernel modules |
| [FLAIL elliptic solver](./flail) | The FFT and multigrid Poisson solvers, their boundary kinds, keys and measured solve times |
| [FLUME source layout](./flume-source) | Every source file of FLUME's common directory and of its CPU and FNL backends, and their build modes |
| [FLUME ghost fill and AMR gravity](./flume-ghost-fill) | How FLUME fills its ghost cells and keeps self-gravity's energy exchange exact on AMR grids |

The step-by-step path to a new application is the [build-your-own tutorial](/tutorials/build-your-own/core). The
generated API reference of every module is under [API](/api/).
