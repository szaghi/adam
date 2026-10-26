# 3. AMR and the staged layer

A refined patch puts 2:1 faces in the grid. The fast path of [chapter 1](./core) no longer applies: the forest drives
each Runge–Kutta stage itself, collects the face fluxes on both sides of every 2:1 face in a **flux register**, and
after the step corrects the coarse cells so that both sides exchange the same flux (the **reflux**). This chapter adds
the procedures of that **staged layer**.

## A refined patch

`src/tests/ripple/amr.ini` is `periodic.ini` with one geometric marker (excerpt: the committed input also sets the
marker's `delta_type`, `delta_fine`, `delta_coarse`, `box_zmin` and `box_zmax`, which the parser requires):

```ini
[amr]
max_level      = 3
markers_number = 1

[amr_marker_1]
mode         = 1             ; AMR_GEO (geometric marker)
geo_type     = primitive-box
box_xmin     = 0.25
box_ymin     = 0.25
box_xmax     = 0.75
box_ymax     = 0.75
target_level = 3
```

The four central blocks of the 4×4 base grid are refined once: the profile crosses 2:1 faces, edges and corners for
the whole run. RIPPLE's `amr_update` (common object) marks the blocks whose centroid lies in the box and calls the
library's `adam%amr_update`, which refines them and moves `q` to the new blocks; refinement happens at initialization
only (`check_input` refuses `[amr] frequency > 0`).

## Which path the forest takes

`forest_object%evolve_one_step` (`src/lib/common/adam_forest_object.F90`) decides per step:

```fortran
   if (int(size(realm), I4P) == 1_I4P .and. self%flux_register%nfaces == 0_I4P) then
      call realm(1)%ptr%advance_one_step_forest(dt=dt)
   else
      do is = 1_I4P, int(size(realm), I4P)
         call realm(is)%ptr%open_step_forest(dt=dt)
      enddo
```

One realm and no 2:1 face: the fast path. Anything else, a single 2:1 face included, takes the staged loop:
`open_step_forest` once, then for every stage `begin_stage_forest` and `end_stage_forest`, then `close_step_forest`,
then the reflux through `apply_reflux_to_stage_forest`. The number of stages comes from `stages_per_step_forest`, which
RIPPLE returns from its SSP scheme.

## The staged procedures

`src/app/ripple/cpu/adam_ripple_cpu_object.F90`, the whole staged layer:

```fortran
   subroutine open_step_forest(self, dt)
   !< Open a step: count it, cap dt, copy q into every stage.
   class(ripple_cpu_object), intent(inout) :: self !< The realm.
   real(R8P),                intent(in)    :: dt   !< Time step from the forest.

   call self%advance_time(dt=dt)
   call self%rk%initialize_stages(field=self%adam%field, q=self%q)
   endsubroutine open_step_forest

   subroutine begin_stage_forest(self, k, K_total, dt, realm)
   !< Begin stage `k`: publish it (`stage_active`) and build its state.
   class(ripple_cpu_object), intent(inout)                   :: self     !< The realm.
   integer(I4P),             intent(in)                      :: k        !< Stage index.
   integer(I4P),             intent(in)                      :: K_total  !< Unused: the forest's stage count.
   real(R8P),                intent(in)                      :: dt       !< Unused: the capped dt is self%dt.
   type(realm_pointer),      intent(inout), optional, target :: realm(:) !< Unused.

   self%stage_active = k
   call self%rk%compute_stage(field=self%adam%field, s=k, dt=self%dt)
   endsubroutine begin_stage_forest

   subroutine end_stage_forest(self, k, K_total, dt, realm, flux_register)
   !< End stage `k`: its residuals, the seam fluxes to the register, the stage assignment.
   class(ripple_cpu_object),    intent(inout)                   :: self          !< The realm.
   integer(I4P),                intent(in)                      :: k             !< Stage index.
   integer(I4P),                intent(in)                      :: K_total       !< Unused.
   real(R8P),                   intent(in)                      :: dt            !< Unused.
   type(realm_pointer),         intent(inout), optional, target :: realm(:)      !< Unused.
   class(flux_register_object), intent(inout), optional         :: flux_register !< Forest's flux register.

   call self%compute_residuals(q=self%rk%q_rk(:,:,:,:,:,k), s=k, flux_register=flux_register)
   call self%rk%assign_stage(field=self%adam%field, s=k, q=self%dq)
   endsubroutine end_stage_forest

   subroutine close_step_forest(self, dt)
   !< Close a step: commit q, advance the time, clear the active stage.
   class(ripple_cpu_object), intent(inout) :: self !< The realm.
   real(R8P),                intent(in)    :: dt   !< Unused: the capped dt is self%dt.

   call self%rk%update_q(field=self%adam%field, dt=self%dt, q=self%q, dq=self%dq)
   self%time = self%time + self%dt
   self%stage_active = 0_I4P
   endsubroutine close_step_forest

   subroutine apply_reflux_to_stage_forest(self, stage, dt, flux_register)
   !< Reflux the committed step (the forest calls it after `close_step_forest`, with the last stage).
   class(ripple_cpu_object),    intent(inout) :: self          !< The realm.
   integer(I4P),                intent(in)    :: stage         !< Stage: the correction acts at the last one.
   real(R8P),                   intent(in)    :: dt            !< Unused: the capped dt is self%dt.
   class(flux_register_object), intent(in)    :: flux_register !< Forest's flux register.

   if (stage == self%rk%nrk) call self%apply_reflux(flux_register=flux_register)
   endsubroutine apply_reflux_to_stage_forest
```

The step is the fast path's loop, split at the stage boundaries. Two rules come with it:
- **`stage_active`** is set to `k` when a stage begins and back to 0 when the step closes. In a forest the peers of a
  realm read it to choose between its stage buffer and its committed state; leaving it stale corrupts their reads.
- **The flux register** reaches `compute_residuals` through `end_stage_forest`. There, `accumulate_seam_fluxes` packs
  the face fluxes of every block face the maps mark as a 2:1 face and hands them to `accumulate_seam_skin` (common
  object), weighted by the stage's SSP coefficient `rk%beta(s)`: the register then holds
  $\sum_s \beta_s F_s$, the flux of the whole step. A coarse face is stored as it is; a fine face is restricted (2×2
  average, 2×1 on a quadtree) into its quadrant of the coarse face.

## The reflux

`apply_reflux` (common object), whole:

```fortran
   subroutine apply_reflux(self, flux_register)
   !< Berger-Colella reflux on the host `q`: every coarse seam cell this rank owns receives
   !< `sgn dt / dx (F_coarse - F_fine)`, the coarse flux the step used replaced by the restricted fine one, so that the
   !< total of u dV is conserved across the 2:1 faces to round-off.
   class(ripple_common_object), intent(inout) :: self          !< The realm.
   class(flux_register_object), intent(in)    :: flux_register !< Forest's flux register.
   integer(I4P)                               :: f, c, axis, sgn, i, j, k !< Counters, face normal and side.

   if (flux_register%nfaces == 0_I4P .or. .not.allocated(flux_register%face)) return
   do f=1, flux_register%nfaces
      associate(face=>flux_register%face(f))
      if (face%coarse_realm /= self%realm_index .or. face%coarse_rank /= mpih%myrank) cycle
      call face_axis_sign(face_code=face%coarse_face, axis=axis, sgn=sgn)
      do c=1, face%nface_cells
         select case(axis)
         case(1_I4P) ; i = merge(self%ni, 1_I4P, sgn > 0) ; j = 1 + mod(c - 1, self%nj) ; k = 1 + (c - 1) / self%nj
         case(2_I4P) ; i = 1 + mod(c - 1, self%ni) ; j = merge(self%nj, 1_I4P, sgn > 0) ; k = 1 + (c - 1) / self%ni
         case default; i = 1 + mod(c - 1, self%ni) ; j = 1 + (c - 1) / self%ni ; k = merge(self%nk, 1_I4P, sgn > 0)
         endselect
         self%q(1,i,j,k,face%coarse_block) = self%q(1,i,j,k,face%coarse_block) +                              &
            sgn * self%dt / self%adam%field%dxyz(axis,face%coarse_block) * (face%F_coarse(1,c,1) - face%F_fine_sum(1,c,1))
      enddo
      endassociate
   enddo
   endsubroutine apply_reflux
```

Before calling it the forest sums the fine contributions over the ranks (`flux_register%reduce_fine_sums`, an
`MPI_ALLREDUCE`), so the rank that owns the coarse block has the complete fine flux and the correction is local.

## Oracles r2 and r3

```bash
./src/tests/ripple/check.sh --legs r2,r3
```

```text
   cells           L1         Linf  order L1 order Linf
    1792   5.7008e-04   4.7944e-03
    7168   1.2859e-04   1.1297e-03     2.148      2.085
   28672   3.1220e-05   2.4886e-04     2.042      2.183
PASS [r2]
total of u dV: 242 rows, largest drift 9.368e-17 of max|u0| times the volume
PASS [r2 conservation N=32]
total of u dV: 482 rows, largest drift 1.631e-16 of max|u0| times the volume
PASS [r2 conservation N=64]
total of u dV: 961 rows, largest drift 2.776e-16 of max|u0| times the volume
PASS [r2 conservation N=128]
28 blocks, largest |u_1 - u_2| = 0.000e+00
PASS [r3]
```

- **Order 2, not 5.** The library's 2:1 transfers (the mean restriction of the coarse ghosts and the reflux of
  point-value fluxes) are second order; across a 2:1 face the scheme is second order, as FLUME measures too.
- **Conservation to round-off.** The total of $u\,dV$ drifts by at most $2.8\cdot10^{-16}$ of $\max|u_0|$ times the
  volume: the register and the reflux make both sides of each 2:1 face exchange the same flux.
- **r3: rank invariance.** The AMR run on 1 and on 2 ranks gives bitwise identical fields (28 blocks).

::: tip What you learned
- A 2:1 face sends the forest down the staged path; the application then splits its step at the stage boundaries.
- `stage_active` is set and reset by the realm; the flux register is filled during the stages and emptied by the
  reflux after the step.
- Across 2:1 faces the scheme is second order and conservative to round-off.
:::

Previous: [2. The library's numerics](./numerics). Next: [4. The GPU twin](./gpu).
