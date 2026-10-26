# 1. A realm and the core layer

This track builds a solver on the ADAM library by reading one: **RIPPLE** (`src/app/ripple`), the scalar linear
advection $u_t + \mathbf a\cdot\nabla u = 0$ on a periodic box, written for this tutorial
([#59](https://github.com/szaghi/adam/issues/59)). It is small on purpose: every chapter quotes whole procedures from
the committed sources, and CI runs its checks on every push.

## The contract in three layers

An application is a **realm**: a type that extends the library's `realm_object` and implements the type-bound
procedures the **forest** (the time loop, `forest_object`) calls. Which of them a configuration reaches depends on the
grid and the number of realms:

| Layer | Needed when | Procedures |
|---|---|---|
| **Core** | always (one realm, no 2:1 face: the fast path) | `initialize_forest`, `compute_local_dt_forest`, `advance_one_step_forest`, `post_step_forest`, `is_done_forest`, `finalize_forest` |
| **Staged** | any 2:1 AMR face, or more than one realm | `stages_per_step_forest`, `open_step_forest`, `begin_stage_forest` (sets `stage_active = k`), `end_stage_forest` (fills the flux register), `close_step_forest` (resets `stage_active = 0`), `apply_reflux_to_stage_forest` |
| **Seams** | forests of several realms | `fill_seam_from_peer_forest`, `pack_seam_cells_forest`, `unpack_seam_cells_forest`; `coupling_descriptor_forest` for the cadence |

`realm_object` is concrete: each of these procedures has a default that stops with a named `error_stop` when it is
called (a few hooks, such as `finalize_mpi_forest`, have working defaults). A procedure a configuration needs and the
application did not override shows up at its first call, at startup or in the first step. This chapter is the core
layer; [chapter 3](./amr) adds the staged one; RIPPLE has no seams layer (see [chapter 4](./gpu#the-seams-layer)).

## The files

| File | Lines | Content |
|---|---|---|
| `common/adam_ripple_common_object.F90` | 353 | the realm: state, input, initial condition, AMR, outputs, the backend-independent procedures |
| `cpu/adam_ripple_cpu_object.F90` | 263 | the CPU realm: residuals, Runge–Kutta, flux register, reflux |
| `cpu/adam_ripple_cpu.F90` | 25 | the driver |
| `fnl/adam_ripple_fnl_object.F90`, `fnl/adam_ripple_fnl_kernels.F90`, `fnl/adam_ripple_fnl.F90` | | the GPU twin ([chapter 4](./gpu)) |

## The driver

The whole program, `src/app/ripple/cpu/adam_ripple_cpu.F90`:

```fortran
!< RIPPLE, scalar linear advection on the ADAM library, CPU backend.
program adam_ripple_cpu
!< RIPPLE, scalar linear advection on the ADAM library, CPU backend.
!<
!< Usage: `adam_ripple_cpu [input.ini]`. One realm, driven by the library's forest: `forest%simulate` initializes the
!< realm, loops over the steps (dt, advance, post-step, termination) and finalizes it.

! ADAM classes, libraries, parameters
use :: adam_forest_object,     only : forest_object
use :: adam_realm_object,      only : realm_pointer
! RIPPLE modules
use :: adam_ripple_cpu_object, only : ripple_cpu_object

implicit none

type(ripple_cpu_object), target :: realm(1)          !< The realm.
type(realm_pointer)             :: realm_pointers(1) !< A pointer to it: what the forest takes (issue #96).
type(forest_object)             :: forest            !< Orchestrator.
character(999)                  :: input_file_name   !< Input file name.

input_file_name = 'input.ini'
if (command_argument_count() > 0) call get_command_argument(1, input_file_name)
realm_pointers(1)%ptr => realm(1)
call forest%simulate(realm=realm_pointers, filename=trim(adjustl(input_file_name)))
endprogram adam_ripple_cpu
```

The forest takes an array of `realm_pointer`, not an array of realms: a realm is polymorphic, and an array of
`class(realm_object)` cannot be indexed safely. `forest%simulate` initializes the realm, loops (dt, advance, post-step,
termination) and finalizes. FLUME's driver adds one branch, `simulate_from_manifest` for a forest manifest
([forest cookbook](/cookbook/forest)).

## The realm

`src/app/ripple/common/adam_ripple_common_object.F90`, the type:

```fortran
type, extends(realm_object) :: ripple_common_object
   !< RIPPLE common realm.
   real(R8P), allocatable :: q(:,:,:,:,:)            !< State, point values u [1, 1-ngc:ni+ngc, ..., nb].
   real(R8P), allocatable :: dq(:,:,:,:,:)           !< Residuals, same shape.
   real(R8P)              :: a(3)=0._R8P             !< Advection velocity.
   real(R8P)              :: cfl=0.5_R8P             !< CFL number.
   real(R8P)              :: time=0._R8P             !< Current time.
   real(R8P)              :: time_max=0._R8P         !< Final time.
   real(R8P)              :: dt=0._R8P               !< Current time step.
   integer(I4P)           :: it=0_I4P                !< Current iteration.
   integer(I4P)           :: it_max=0_I4P            !< Iterations number; <= 0: run to time_max.
   integer(I4P)           :: history_unit=0_I4P      !< Unit of the history file (rank 0).
   contains
      ! ... the type-bound procedures (lines 50-63 of the source)
endtype ripple_common_object
```

`realm_object` already owns the grid, the tree, the field (block geometry), the maps, the WENO and Runge–Kutta objects
(`self%adam%grid`, `self%adam%field`, `self%weno`, `self%rk`, …). The application owns its state: `q` has the library's
5-D layout `(nv, i, j, k, b)` with the ghost cells, here `nv = 1`.

The CPU realm extends it again, `ripple_cpu_object` in `cpu/adam_ripple_cpu_object.F90`, adding the face fluxes and
the procedures that touch `q` on the host.

## Initialization

The first lines of `initialize_ripple` (common) build the library objects; the rest reads the `[ripple]` and `[time]`
keys, calls `check_input`, allocates `q` and `dq` and refines uniformly:

```fortran
   call self%io%initialize(filename=trim(filename), verbose=.true.)
   ! the blocks budget counts the block-sized fields per block: q, dq, the three face fluxes, the stored stages
   call self%io%file_parameters%get(section_name='runge_kutta', option_name='scheme', val=rk_scheme, error=error)
   stages = rk_stored_stages_number(scheme=rk_scheme)
   if (error > 0 .or. stages < 0_I4P) call mpih%error_stop(msg=': failed to load [runge_kutta].(scheme)')
   call self%realm_object%initialize(filename=filename, memory_avail=memory_avail, nv=1_I4P, &
                                     fields_number=5_I4P + stages, verbose=.true.)
   ! every face periodic; the tree and the maps were built before the faces were known: rebuild them (issue #103)
   call self%adam%grid%set_bc_type(bc_type=[BC_PERIODIC, BC_PERIODIC, BC_PERIODIC, BC_PERIODIC, BC_PERIODIC, BC_PERIODIC])
   call self%adam%tree%rebuild_neighborhood(grid=self%adam%grid)
   call self%adam%make_comm_local_maps_ghost_bc
```

Three rules of the library show here:
- **`fields_number`** is the number of block-sized scalar arrays the application allocates per block. The library
  divides the memory by it to decide how many blocks a rank can hold; under-count it and the run over-commits memory
  ([GPU cookbook](/cookbook/gpu)).
- **The face kinds come after `realm_object%initialize`**, which builds the tree and the maps before they are known;
  `set_bc_type` must be followed by `rebuild_neighborhood` and `make_comm_local_maps_ghost_bc`, or a periodic face of
  an unrefined grid stays a boundary row ([#103](https://github.com/szaghi/adam/issues/103)).
- **Refuse what you cannot honour.** `check_input` stops on every library key RIPPLE would otherwise ignore: a non-SSP
  scheme, a centred WENO, `weights = si`, `ngc` below the WENO stencil, runtime AMR, solids.

The CPU realm's `initialize_forest` then starts MPI, allocates the fluxes, refines, sets the initial condition and
writes the first output:

```fortran
   subroutine initialize_forest(self, filename, realms_number, memory_avail, nv, verbose)
   !< Initialize the realm: common data, face fluxes, initial AMR, initial condition, first output.
   class(ripple_cpu_object), intent(inout)        :: self          !< The realm.
   character(*),             intent(in)           :: filename      !< Input file name.
   integer(I4P),             intent(in), optional :: realms_number !< Unused: RIPPLE is single-realm.
   real(R8P),                intent(in), optional :: memory_avail  !< Unused: the budget comes from the MPI handler.
   integer(I4P),             intent(in), optional :: nv            !< Unused: nv = 1.
   logical,                  intent(in), optional :: verbose       !< Unused.

   call mpih%initialize(do_mpi_init=.true., verbose=.true.)
   call self%initialize_ripple(filename=filename, memory_avail=mpih%memory_avail)
   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, nb=>self%nb)
      allocate(self%flx(1,0:ni,1:nj,1:nk,nb), self%fly(1,1:ni,0:nj,1:nk,nb), self%flz(1,1:ni,1:nj,0:nk,nb), source=0._R8P)
   endassociate
   if (self%amr%markers_number > 0) call self%amr_update
   call self%set_initial_conditions
   call self%save_output(fields=.true.)
   endsubroutine initialize_forest
```

## The time step, the step, the end

The forest asks every realm for its local time step and takes the minimum over the ranks:

```fortran
   subroutine compute_local_dt_forest(self, dt_local)
   !< Return the local CFL time step, `CFL / max_b sum_d |a_d| / dx_d` (the forest takes the minimum over the ranks).
   class(ripple_common_object), intent(in)  :: self     !< The realm.
   real(R8P),                   intent(out) :: dt_local !< Local time step.
   integer(I4P)                             :: b        !< Counter.

   dt_local = huge(1._R8P)
   do b=1, self%blocks_number
      dt_local = min(dt_local, self%cfl / max(sum(abs(self%a) / self%adam%field%dxyz(:,b)), tiny(1._R8P)))
   enddo
   endsubroutine compute_local_dt_forest
```

On the fast path (one realm, no 2:1 face) the forest hands the whole step to the realm:

```fortran
   subroutine advance_one_step_forest(self, dt)
   !< Advance one whole SSP Runge-Kutta step (fast path: no 2:1 face).
   class(ripple_cpu_object), intent(inout) :: self !< The realm.
   real(R8P),                intent(in)    :: dt   !< Time step from the forest.
   integer(I4P)                            :: s    !< Stage counter.

   call self%advance_time(dt=dt)
   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   do s=1, self%rk%nrk
      call self%rk%compute_stage(field=self%adam%field, s=s, dt=self%dt)
      call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,s))
      call self%rk%assign_stage(field=self%adam%field, s=s, q=self%dq)
   enddo
   call self%rk%update_q(field=self%adam%field, dt=self%dt, q=self%q, dq=self%dq)
   self%time = self%time + self%dt
   endsubroutine advance_one_step_forest
```

The Runge–Kutta object holds the stages; the application supplies only the residuals of a state
([chapter 2](./numerics)). After the step the forest calls `post_step_forest` (RIPPLE writes a history row, and the
fields every `[IO] it_save` steps and at the end), then `is_done_forest`:

```fortran
   subroutine is_done_forest(self, done)
   !< Return true when `it_max` steps are done or, with `it_max <= 0`, when `time_max` is reached.
   class(ripple_common_object), intent(in)  :: self !< The realm.
   logical,                     intent(out) :: done !< True if the realm is done.

   if (self%it_max > 0_I4P) then
      done = self%it >= self%it_max
   else
      done = self%time >= self%time_max * (1._R8P - 10._R8P * epsilon(1._R8P))
   endif
   endsubroutine is_done_forest
```

and, once every realm is done, `finalize_forest` (RIPPLE closes its history file; the forest finalizes MPI once).

## Build and run

```bash
fobis build --mode ripple-cpu-gnu        # exe/adam_ripple_cpu
mkdir -p run/ripple && cp src/tests/ripple/periodic.ini run/ripple/input.ini
cd run/ripple && mpirun -np 2 ../../exe/adam_ripple_cpu input.ini > log.txt
```

`periodic.ini` advects $u_0 = \sin 2\pi x\,\sin 2\pi y$ at $\mathbf a = (1, 0.5)$ across the periodic unit square to
$t = 1$, on 4×4 blocks of 8×8 cells. The run writes `periodic-history.dat` (step, time, total of $u\,dV$) and the field
files of step 0 and of the last step (`periodic-000000000.xdmf`, `periodic-000000120.xdmf` and the per-rank `.h5`).
`periodic-history.dat`, first and last rows:

```text
# it time total(u dV)
        0  0.000000000000000E+000  3.122502256758253E-017
      120  9.999999999999989E-001 -1.040834085586084E-017
```

120 steps to $t = 1$; the total of $u\,dV$ is zero for this profile, and stays zero to round-off: the scheme is
conservative.

::: tip What you learned
- An application is a realm type extending `realm_object`, plus a driver that hands it to `forest%simulate`.
- The core layer is six procedures: initialize, local dt, one step, post-step, done, finalize.
- The library owns the grid and the numerics objects; the application owns its state and declares how many
  block-sized fields it allocates.
:::

Next: [2. The library's numerics](./numerics).
