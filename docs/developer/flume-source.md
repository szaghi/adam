# FLUME source layout

FLUME has one backend-independent directory and two backends, under `src/app/flume/`. Every type of `common/` is
aggregated into `flume_common_object`, the base type both backends extend. The user documentation (models, numerics,
input, verification) starts at [FLUME](/applications/flume/).

| Directory | Content | Build modes |
|-----------|---------|-------------|
| [`common/`](#common) | Input sections, pointwise physics libraries, `flume_common_object` | built by every mode |
| [`cpu/`](#cpu-backend) | MPI (+ OpenMP) backend, `flume_cpu_object` | `flume-cpu-gnu*` |
| [`fnl/`](#fnl-backend) | GPU backend on the [FNL library](./fnl), `flume_fnl_object` | `flume-fnl-nvf*`, `flume-fnl-omp-amd*` |

## Common

The pointwise libraries (`*_library.F90`) are `pure`, take explicit-size arguments and are tagged `!$acc routine seq`
+ `!$omp declare target`, so the CPU loops and the GPU kernels call the same source.

| File | Content |
|------|---------|
| `adam_flume_parameters.F90` | Variable indexes (conservative `IQ_*`, auxiliary `IA_*`), model ids, accepted option values, `strip_control` |
| `adam_flume_physics_object.F90` | `[physics]` and `[eos]`: physical model (`euler`, `euler-de`, `mhd-ideal`), passive scalars, `gamma` or `cp`, `cv`, the equation of state (ideal, isothermal, barotropic), variable counts; primitive-to-conservative dispatch |
| `adam_flume_reference_object.F90` | `[reference]`: dimensional input converted to code units on the host (issue #49) |
| `adam_flume_gravity_object.F90` | `[gravity]`: external gravity (issue #73) and the self-gravity keys, the Poisson solver and its boundaries (issue #81) |
| `adam_flume_eos_tables_library.F90` | Barotropic equation-of-state tables, continuous piecewise polytropes (issue #100) |
| `adam_flume_dissipation_object.F90` | `[physics]` dissipative coefficients (issue #65): viscosity, conductivity, resistivity, each as a coefficient or its number, and the power-law viscosity |
| `adam_flume_cooling_object.F90` | `[cooling]` radiative cooling (issue #93): the node table, the particle mass, the time-step limiter, the precomputed exponents and evolution function, the radiated energy |
| `adam_flume_cooling_library.F90` | Townsend's exact integration shared by both backends: the cooling function, the evolution function and its inverse, the cooled internal energy of a cell |
| `adam_flume_sinks_object.F90` | `[sinks]` sink particles (issue #94): the replicated list (masses, positions, velocities, accelerations), the restart file and the sinks history; the deposit and the gather live in the common object |
| `adam_flume_dissipation_library.F90` | Pointwise dissipative physics shared by both backends: the viscous and heat-conduction flux of a direction, the Ohmic flux (MHD), the diffusivity of the time-step limit |
| `adam_flume_mhd_object.F90` | `[mhd]`: divergence control (`glm`, `eglm`, `none`), GLM speed, damping and check, div(B) monitor, positivity floors |
| `adam_flume_numerics_object.F90` | `[numerics]`: `scheme_space` (`weno`, `weno-riemann`, `hybrid-kep-weno` and its `hybrid_*` keys), `reconstruction_variables`, Riemann solver, flux correction and sensor, `reflux` and its policy and cadence, `positivity_limiter`, `dissipative_order` |
| `adam_flume_euler_library.F90` | Pointwise Euler physics: conversions, fluxes, Roe average and eigenvectors, flux splitting, face states, LLF/HLL/HLLC Riemann solvers; its body `adam_flume_euler_library_body.INC` is shared with the dual-energy library |
| `adam_flume_euler_de_library.F90`, `adam_flume_euler_ps_library.F90`, `adam_flume_euler_de_ps_library.F90` | The Euler physics with dual energy (issue #91), passive scalars (issue #92), or both |
| `adam_flume_euler_iso_library.F90`, `adam_flume_euler_baro_library.F90` | Isothermal and barotropic Euler physics (issue #100), one body `adam_flume_euler_iso_library_body.INC` behind the EOS hooks |
| `adam_flume_mhd_iso_library.F90`, `adam_flume_mhd_baro_library.F90` | Isothermal and barotropic ideal MHD physics and Riemann solvers (issue #100), one body `adam_flume_mhd_iso_library_body.INC` |
| `adam_flume_mhd_library.F90` | Pointwise ideal MHD physics: conversions, fluxes, fast speed, Roe–Balsara eigensystem (with and without GLM), flux splitting |
| `adam_flume_mhd_riemann_library.F90` | MHD face states and Riemann solvers (LLF, HLL, HLLD; exact GLM subsystem) of `weno-riemann` |
| `adam_flume_bc_object.F90` | `[bc_*]`: boundary condition types (extrapolation, inflow, inviscid, modelled, no-slip and isothermal walls, periodic, hydrostatic, digital-filter inflow, non-reflecting, shock generator), inflow states, wall velocities and temperatures; `wall_noslip_ghost` (the no-slip ghost state); `realm_edge_face`, `realm_edge_donor` (the face condition and the donor of a ghost at a realm edge or corner, shared by the CPU and FNL backends) |
| `adam_flume_ic_object.F90` | `[initial_conditions]`: 29 initial conditions (uniform, channel, boundary layer, Taylor–Green vortex, isotropic turbulence, isentropic vortex, Riemann regions, Shu–Osher, rotated Riemann, GLM pulse, div(B) peak, MHD linear wave, CPAW, magnetised vortex, Orszag–Tang, rotor, field loop, linear, polynomial, sine wave, Couette, Becker shock, isothermal and polytropic atmospheres, Rayleigh–Taylor, polytrope, singular isothermal sphere, Bonnor–Ebert, Boss–Bodenheimer); init-time AMR passes |
| `adam_flume_hybrid_library.F90` | Pointwise central flux of `hybrid-kep-weno` (issue #119): the 6th-order kinetic-energy- and pressure-equilibrium-preserving face flux (Pirozzoli 2011, Shima et al. 2021) |
| `adam_flume_wall_model_library.F90` | Equilibrium wall model of `wall-modelled` faces (issue #119): wall shear stress and heat flux from a sample in the outer layer (Kawai & Larsson 2012) |
| `adam_flume_inflow_df_object.F90` | `inflow-digital-filter` faces (issue #119): mean profile table plus fluctuations with prescribed Reynolds stresses and length scales (Klein, Sadiki & Janicka 2003; Touber & Sandham 2009; Xie & Castro 2008) |
| `adam_flume_forcing_object.F90` | `[forcing]` (issue #119): `mass-flow` forcing of a periodic channel, holding the bulk momentum along an axis; the forcing history |
| `adam_flume_seam_sync_object.F90` | Per-stage seam flux synchronisation data of the positivity limiter at 2:1 seams (issue #50) |
| `adam_flume_time_object.F90` | `[time]`: CFL, iteration and time limits |
| `adam_flume_diagnostics_object.F90` | `[diagnostics]`: conservation (and div(B)) history cadence |
| `adam_flume_common_object.F90` | `flume_common_object`: initialization and cross-section checks, AMR markers (box, gradient, Löhner, solid, Jeans) and runtime regrids, gravity and self-gravity sources (the solve by the library's FLAIL), sink particles' deposit and gather, seam flux accumulation for reflux, immersed-boundary spacing, restart, slices, auxiliary and MHD derived fields |
| `adam_flume_common_library.F90` | Barrel re-export of all common modules |

## CPU backend

The CPU backend (`flume-cpu`) is the MPI (+ optional OpenMP) implementation, extending `flume_common_object`.

| File | Description |
|------|-------------|
| `adam_flume_cpu.F90` | Entry point (program) |
| `adam_flume_cpu_object.F90` | `flume_cpu_object`: ghost update and boundary conditions, WENO residual (with the immersed-boundary variant), Runge-Kutta integration, reflux accumulation and application, output |
| `adam_flume_cpu_<model>[_<solver>]_kernels.F90` | One module per physical model (`euler`, `euler_de`, `euler_ps`, `euler_de_ps`, `euler_iso`, `euler_baro`, `mhd`, `mhd_glm`, `mhd_eglm`, `mhd_iso`, `mhd_glm_iso`, `mhd_baro`, `mhd_glm_baro`) and per Riemann solver of `weno-riemann` (`llf`, `hll`, `hllc` for Euler, `hlld` for MHD): the face fluxes of that model, so no kernel branches on the model |
| `adam_flume_cpu_*_kernels_agnostic.INC` | The model-agnostic kernel bodies (face, Riemann face, ghost, auxiliary, positivity, MHD) the model modules include |
| `adam_flume_cpu_dissipation_kernels.F90` | Viscous, heat-conduction and Ohmic (MHD) face fluxes of order 2 and 4 (issue #65), model-agnostic, added to the inviscid face fluxes before the seam accumulation |
| `adam_flume_cpu_hybrid_kernels.F90` | Shock sensor and central-flux overwrite of `hybrid-kep-weno` (issue #119) |
| `adam_flume_cpu_forcing_kernels.F90` | Mass-flow forcing of `[forcing]` (issue #119) |

### CPU build modes

```bash
fobis build --mode flume-cpu-gnu          # MPI
fobis build --mode flume-cpu-gnu-omp      # MPI + OpenMP
fobis build --mode flume-cpu-gnu-debug    # MPI, runtime checks
fobis build --mode flume-cpu-gnu-fft      # MPI + FFTW (FFT Poisson solver of self-gravity), exe adam_flume_cpu_fft
fobis build --mode flume-cpu-gnu-fft-debug  # the same, runtime checks
```

## FNL backend

The FNL backend (`flume-fnl`) is the GPU implementation on the FNL library: OpenACC directives, each with its OpenMP
offload twin. The fields live on the device. The host copies are refreshed for output, restart and the coarse-fine
seam fluxes of the reflux register; the seam face fluxes of a rank go to the host register in one copy per stage.

| File | Description |
|------|-------------|
| `adam_flume_fnl.F90` | Entry point (program) |
| `adam_flume_fnl_object.F90` | `flume_fnl_object`: device data, ghost update, residual, Runge-Kutta integration, reflux, output |
| `adam_flume_fnl_dissipation_kernels.F90` | Device twins of the CPU viscous, heat and Ohmic face-flux kernels (issue #65) |
| `adam_flume_fnl_kernels.F90` | Device kernels: face fluxes and flux difference (with the immersed-boundary variant), boundary conditions, time-step and conservation reductions, auxiliary fields, seam skin packing and reflux application, radiative cooling and the cooling time (issue #93) |
| `adam_flume_fnl_<model>[_<solver>]_kernels.F90` | Device twins of the CPU per-model and per-Riemann-solver kernel modules (same model and solver names) |
| `adam_flume_fnl_*_kernels_agnostic.INC` | The model-agnostic device kernel bodies the model modules include |
| `adam_flume_fnl_*_kernel_body.INC` | Per-stencil-width bodies of the face-flux and gravity edge-correction kernels, flux splitting and `weno-riemann` (`face_fluxes`, `gravity_edge`, `riemann_face_fluxes`, `riemann_gravity_edge`), included by the agnostic kernel files |
| `adam_flume_fnl_hybrid_kernels.F90` | Device twins of the CPU `hybrid-kep-weno` kernels (issue #119) |
| `adam_flume_fnl_forcing_kernels.F90` | Device twins of the CPU mass-flow forcing kernels (issue #119) |
| `adam_flume_fnl_library.F90` | Barrel re-export of all FNL modules |

### FNL build modes

```bash
fobis build --mode flume-fnl-nvf --varset local_nvf          # NVIDIA, OpenACC (verified)
fobis build --mode flume-fnl-nvf-debug --varset local_nvf    # NVIDIA, OpenACC, runtime checks
fobis build --mode flume-fnl-nvf-fft --varset local_nvf      # NVIDIA, OpenACC + cuFFT (FFT Poisson solver of self-gravity), exe adam_flume_fnl_fft
fobis build --mode flume-fnl-nvf-fft-debug --varset local_nvf  # the same, runtime checks
fobis build --mode flume-fnl-omp-amd                          # AMD, OpenMP offload
fobis build --mode flume-fnl-omp-amd-debug                    # AMD, OpenMP offload, runtime checks
```

The AMD build runs on an MI300X (amdflang). On the `profile-couette3d` case (size `large`, WENO-5, 20 steps, one GPU)
its loop takes 50.7 s, down from 632.2 s before the WENO port of commit `e4caa8ea` (12.5x;
[profile-couette3d README](https://github.com/szaghi/adam/tree/develop/src/tests/flume/speedup/profile-couette3d)).
Its results have not been compared with the FNL regression goldens.
