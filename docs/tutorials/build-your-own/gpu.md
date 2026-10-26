# 4. The GPU twin

RIPPLE's FNL backend runs the same realm on an NVIDIA GPU with OpenACC. It is a **sibling** of the CPU realm, not an
extension of it: `ripple_fnl_object` extends the same common realm, keeps the state on the device, and replaces the
procedures that touch `q`. The forest does not know which backend it drives.

## The driver

`src/app/ripple/fnl/adam_ripple_fnl.F90` is the CPU driver of [chapter 1](./core) with one type changed:

```fortran
use :: adam_ripple_fnl_object, only : ripple_fnl_object
type(ripple_fnl_object), target :: realm(1)          !< The realm.
```

(the lines that differ, besides the header comment and the `program`/`endprogram` names; the rest of the file is the
CPU driver's)

## The device state

`src/app/ripple/fnl/adam_ripple_fnl_object.F90`, the type:

```fortran
type, extends(ripple_common_object) :: ripple_fnl_object
   !< RIPPLE FNL realm.
   type(field_fnl_object) :: field_fnl                !< Field helper: coordinates, maps, ghost exchange.
   type(rk_fnl_object)    :: rk_fnl                   !< Runge-Kutta helper: the stages on the device.
   type(weno_fnl_object)  :: weno_fnl                 !< WENO helper: the coefficients on the device.
   real(R8P), pointer     :: q_gpu(:,:,:,:,:)=>null() !< State [nb, 1-ngc:ni+ngc, ..., 1].
   real(R8P), pointer     :: dq_gpu(:,:,:,:,:)=>null()  !< Residuals, same shape.
   real(R8P), pointer     :: flx_gpu(:,:,:,:,:)=>null() !< X-face fluxes [nb, 0:ni, 1:nj, 1:nk, 1].
   real(R8P), pointer     :: fly_gpu(:,:,:,:,:)=>null() !< Y-face fluxes [nb, 1:ni, 0:nj, 1:nk, 1].
   real(R8P), pointer     :: flz_gpu(:,:,:,:,:)=>null() !< Z-face fluxes [nb, 1:ni, 1:nj, 0:nk, 1].
   real(R8P), pointer     :: face_gpu(:,:)=>null()      !< One face's skin or reflux mismatch [1, largest face].
   real(R8P), allocatable :: face(:,:)                  !< Its host copy, the same extent.
   contains
      ! ... the type-bound procedures (lines 50-68 of the source)
endtype ripple_fnl_object
```

Three things change from the CPU realm:
- **The library's FNL helpers.** `field_fnl`, `rk_fnl` and `weno_fnl` are the device side of the library objects the
  common realm already owns: the ghost exchange (GPU-direct MPI), the Runge–Kutta stages, the WENO coefficients. They
  are components of the realm, initialized from the common objects by argument.
- **The layout is transposed**, `(b, i, j, k, v)` against the host `(v, i, j, k, b)`: the block index is stride-1, so
  it is the innermost loop of every kernel and adjacent threads touch adjacent memory.
- **Device arrays are FUNDAL pointers**, allocated by `dev_alloc`, each call checked:

```fortran
   call dev_alloc(fptr_dev=self%q_gpu, lbounds=[1,1-ngc,1-ngc,1-ngc,1], ubounds=[nb,ni+ngc,nj+ngc,nk+ngc,1], &
                  init_value=0._R8P, ierr=ierr) ; call check('q_gpu')
```

`check` stops the run with the array's name when `ierr /= 0` (code 101 is a device out of memory). The host copy of a
device buffer has the device buffer's extent: a copy between arrays of different extents goes through a host temporary
and fails on the device ([GPU cookbook](/cookbook/gpu)).

## Initialization

`initialize_forest`, whole:

```fortran
   subroutine initialize_forest(self, filename, realms_number, memory_avail, nv, verbose)
   !< Initialize the realm: MPI and the device (once per process), the common data within the device budget, the FNL
   !< helpers, the device data; initial AMR and initial condition on the host, then the state to the device.
   class(ripple_fnl_object), intent(inout)        :: self          !< The realm.
   character(*),             intent(in)           :: filename      !< Input file name.
   integer(I4P),             intent(in), optional :: realms_number !< Unused: RIPPLE is single-realm.
   real(R8P),                intent(in), optional :: memory_avail  !< Unused: the budget comes from the device.
   integer(I4P),             intent(in), optional :: nv            !< Unused: nv = 1.
   logical,                  intent(in), optional :: verbose       !< Unused.
   logical                                        :: is_mpi        !< MPI already initialized.

   if (.not.mpih_fnl_is_initialized) then
      call MPI_INITIALIZED(is_mpi, mpih_fnl%error)
      call mpih_fnl%initialize(do_mpi_init=.not.is_mpi, do_device_init=.true., verbose=.true.)
      mpih_fnl_is_initialized = .true.
      mpih_fnl_device_sharers = count_device_sharers()
   endif
   ! the blocks budget is the device memory, shared by the ranks bound to the same device
   call self%initialize_ripple(filename=filename, memory_avail=real(mpih_fnl%dev_memory_total, R8P) / 1e9_R8P / &
                                                               real(mpih_fnl_device_sharers, R8P))
   call self%field_fnl%initialize(grid=self%adam%grid, field=self%adam%field, maps=self%adam%maps, verbose=.true.)
   call self%rk_fnl%initialize(grid=self%adam%grid, field=self%adam%field, rk=self%rk)
   call self%weno_fnl%initialize(weno=self%weno)
   call self%allocate_gpu
   if (self%amr%markers_number > 0) call self%amr_update
   call self%set_initial_conditions
   call self%field_fnl%copy_cpu_gpu(field=self%adam%field, maps=self%adam%maps, verbose=.true.)
   call self%copy_state(to_device=.true.)
   call self%save_output(fields=.true.)
   endsubroutine initialize_forest
```

The block budget is the device memory divided by the ranks that share the device; the initial AMR and the initial
condition run on the host, then the maps and the state go to the device once. The backend also overrides
`finalize_mpi_forest`, which the forest calls once after every realm is finalized, to close `mpih_fnl`.

## A kernel

The face fluxes of chapter 2 on the device, `compute_face_fluxes_dev` in `fnl/adam_ripple_fnl_kernels.F90`, the loop
nest:

```fortran
   ap = max(a, 0._R8P)
   am = min(a, 0._R8P)
   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(weno_a_gpu,weno_p_gpu,weno_d_gpu,q_gpu,fl_gpu) &
   !$acc& firstprivate(di,dj,dk,ni,nj,nk,blocks_number,S,ap,am,weno_zeps) private(v,vr)
   do k=1-dk, nk
   do j=1-dj, nj
   do i=1-di, ni
   do b=1, blocks_number
      !$acc loop seq
      do m=1-S, S-1
         v(2,m+S) = ap * q_gpu(b,i+m*di,j+m*dj,k+m*dk,1)
         v(1,m+S) = am * q_gpu(b,i+(m+1)*di,j+(m+1)*dj,k+(m+1)*dk,1)
      enddo
      call weno_reconstruct_upwind_dev(S=S, weno_a=weno_a_gpu, weno_p=weno_p_gpu, weno_d=weno_d_gpu, &
                                       weno_zeps=weno_zeps, weno_rmu=1._R8P, V=v, VR=vr)
      fl_gpu(b,i,j,k,1) = vr(1) + vr(2)
   enddo
   enddo
   enddo
   enddo
```

It is the CPU loop with the block index innermost: host scalars are `firstprivate`, device arrays `DEVICEVAR` (FUNDAL's
macro for `deviceptr`), and the WENO primitive is the library's device twin, `weno_reconstruct_upwind_dev`. RIPPLE's
kernels carry OpenACC directives only (NVIDIA); FLUME's kernels also carry an `!$omp` twin of each directive for the
AMD OpenMP-offload build ([GPU cookbook](/cookbook/gpu)).

## Build and check: oracle r4

```bash
fobis build --mode ripple-fnl-nvf --varset local_nvf    # exe/adam_ripple_fnl
./src/tests/ripple/check.sh --legs r4                   # after r1 and r2: it compares against their CPU runs
```

```text
16 blocks, largest |u_1 - u_2| = 1.665e-15
PASS [r4 periodic]
28 blocks, largest |u_1 - u_2| = 4.241e-14
PASS [r4 amr]
total of u dV: 482 rows, largest drift 1.388e-17 of max|u0| times the volume
PASS [r4 conservation]
RIPPLE oracles (r4) PASSED
```

The FNL runs of `periodic.ini` and `amr.ini` agree with the CPU runs to round-off, and the device reflux conserves as
the host one does. These runs were made on a WSL2 workstation, which checks correctness only: never time a GPU run
there.

## The seams layer

RIPPLE has one realm, so it stops at the staged layer. A forest of several realms glued by seams adds the third layer:
`fill_seam_from_peer_forest` (or, across ranks, `pack_seam_cells_forest` and `unpack_seam_cells_forest`) fills a
realm's seam ghosts from its peer, `coupling_descriptor_forest` reports what the forest checks before it lets two
realms exchange stage data, and on FNL `after_topology_build_forest` moves the seam maps to the device once the forest
has built them. FLUME implements all of them (`src/app/flume/cpu/adam_flume_cpu_object.F90`,
`src/app/flume/fnl/adam_flume_fnl_object.F90`); the [forest cookbook](/cookbook/forest) shows the inputs, and
[R11](/cookbook/forest#r11-forest-on-fnl) the FNL side.

::: tip What you learned
- A GPU backend is a sibling realm: same common realm, device state, the library's FNL helpers as components.
- Device arrays are transposed with the block index stride-1; every `dev_alloc` is checked.
- A kernel is the CPU loop with OpenACC directives and the library's device primitives; the oracle compares the two
  backends to round-off.
:::

Previous: [3. AMR and the staged layer](./amr). This is the end of the build-your-own track.
