# ADAM FNL Library — `src/lib/fnl`

The FNL backend is ADAM's **GPU acceleration layer**, built on top of the [FUNDAL](https://github.com/szaghi/FUNDAL) GPU memory management library. The same sources compile as OpenACC with nvfortran (`*-fnl-nvf` modes) or as OpenMP target offloading with amdflang (`*-fnl-omp-amd` modes, `DEV_OMP`): every kernel carries an `!$acc` directive and its `!$omp` twin. It follows a consistent two-tier pattern for every subsystem: a **wrapper object** that holds the device-resident arrays, paired with a **kernel module** that contains the device subroutines implementing the actual computation. The wrappers are separate types, not extensions of the common classes (the Poisson solvers excepted): each is a component of the application's FNL type, initialised from the realm's common objects passed as arguments.

**Documentation:** <https://szaghi.github.io/adam/concepts/library-fnl>

## Module summary

| Module | Role | Extends |
|--------|------|---------|
| `adam_fnl_library` | Aggregate entry point | — |
| `adam_fnl_field_object` | GPU field wrapper: device geometry, maps, ghost exchange | — |
| `adam_fnl_field_kernels` | Gradient, L2 norm, ghost-cell pack/unpack | — |
| `adam_fnl_weno_object` | GPU WENO coefficient mirror + ROR tables | — |
| `adam_fnl_weno_kernels` | Upwind WENO reconstruction (`!$acc routine seq`) | — |
| `adam_fnl_rk_object` | GPU RK stage storage and update dispatch | — |
| `adam_fnl_rk_kernels` | Stage assign, accumulate, low-storage, final update | — |
| `adam_fnl_ib_object` | GPU distance field + eikonal BC wrapper | — |
| `adam_fnl_ib_kernels` | Eikonal evolution, sphere distance, momentum inversion | — |
| `adam_fnl_maps_object` | GPU communication index tables + MPI buffer staging | — |
| `adam_fnl_mpih_object` | FUNDAL MPI handler alias | — |
| `adam_fnl_poisson_fft_object` | FFT Poisson solver, cuFFT on each rank's GPU (the slab transform; not re-exported by `adam_fnl_library`) | `poisson_fft_object` |
| `adam_fnl_poisson_mg_object` | Multigrid Poisson solver on AMR, the FCG iteration and the V-cycle on each rank's GPU (not re-exported by `adam_fnl_library`) | `poisson_mg_object` |
| `adam_fnl_fdv_operators_library` | Device-callable FD/FV spatial operators | — |
| `adam_fnl_mpih_global` | Singleton — `mpih_fnl` | — |

---

ADAM is dual-licensed under MIT or LGPLv3 (`LICENSE.mit.md`, `LICENSE.lgpl3.md`).
