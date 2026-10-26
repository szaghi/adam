# 2. The library's numerics

The realm of [chapter 1](./core) advances its state with the library's Runge–Kutta object; this chapter supplies what
that object needs, the residuals, and builds them from the library's WENO. Then the first oracle measures the order.

## The scheme

RIPPLE is the Shu–Osher conservative finite difference on point values at the cell centres:

$$\frac{du_{ijk}}{dt} = -\sum_d \frac{F^d_{i+1/2} - F^d_{i-1/2}}{\Delta x_d},\qquad F^d = F^{d,+} + F^{d,-},\quad
F^{d,\pm} = a_d^{\pm}\,u,$$

with $a^+ = \max(a, 0)$, $a^- = \min(a, 0)$, and each half-flux reconstructed at the face by upwind WENO from its own
side. With the library's `weno-u-5` and an SSP Runge–Kutta scheme the method is fifth order in space.

## The residuals

`compute_residuals` of `src/app/ripple/cpu/adam_ripple_cpu_object.F90`, whole:

```fortran
   subroutine compute_residuals(self, q, s, flux_register)
   !< Compute the residuals `dq` of state `q`: ghost update, WENO face fluxes of each non-null axis, flux difference. On
   !< the staged path (`flux_register` present), the seam face fluxes of stage `s` go to the register.
   class(ripple_cpu_object),    intent(inout)           :: self          !< The realm.
   real(R8P),                   intent(inout)           :: q(1:,1-self%ngc:,1-self%ngc:,1-self%ngc:,1:) !< State.
   integer(I4P),                intent(in),    optional :: s             !< Runge-Kutta stage.
   class(flux_register_object), intent(inout), optional :: flux_register !< Forest's flux register.
   integer(I4P)                                         :: b, i, j, k    !< Counters.

   call self%update_ghost(q=q)
   associate(null=>self%adam%grid%null_xyz, dxyz=>self%adam%field%dxyz)
   if (.not.null(1)) call face_fluxes(di=1, dj=0, dk=0, a=self%a(1), fl=self%flx)
   if (.not.null(2)) call face_fluxes(di=0, dj=1, dk=0, a=self%a(2), fl=self%fly)
   if (.not.null(3)) call face_fluxes(di=0, dj=0, dk=1, a=self%a(3), fl=self%flz)
   if (present(flux_register) .and. present(s)) then
      if (flux_register%nfaces > 0_I4P) call self%accumulate_seam_fluxes(s=s, flux_register=flux_register)
   endif
   do b=1, self%blocks_number
      do k=1, self%nk
         do j=1, self%nj
            do i=1, self%ni
               self%dq(1,i,j,k,b) = -(merge(0._R8P, (self%flx(1,i,j,k,b) - self%flx(1,i-1,j,k,b)) / dxyz(1,b), null(1)) + &
                                      merge(0._R8P, (self%fly(1,i,j,k,b) - self%fly(1,i,j-1,k,b)) / dxyz(2,b), null(2)) + &
                                      merge(0._R8P, (self%flz(1,i,j,k,b) - self%flz(1,i,j,k-1,b)) / dxyz(3,b), null(3)))
            enddo
         enddo
      enddo
   enddo
   endassociate
   contains
      subroutine face_fluxes(di, dj, dk, a, fl)
      !< Face fluxes along the axis of unit step (di, dj, dk): the split fluxes `a+ u` (cells `m = 1-S ... S-1`) and
      !< `a- u` (cells `m+1`) of the stencil of face `i+1/2`, in the layout of the WENO upwind primitive, reconstructed
      !< and summed.
      integer(I4P), intent(in)    :: di, dj, dk                    !< Unit step.
      real(R8P),    intent(in)    :: a                             !< Velocity component.
      real(R8P),    intent(inout) :: fl(1:,1-di:,1-dj:,1-dk:,1:)   !< Face fluxes.
      real(R8P)                   :: v(2,2*S_MAX-1)                !< Packed stencil.
      real(R8P)                   :: vr(2)                         !< Reconstructed split fluxes.
      integer(I4P)                :: S, m                          !< WENO half stencil, counter.

      S = self%weno%S
      do b=1, self%blocks_number
         do k=1-dk, self%nk
            do j=1-dj, self%nj
               do i=1-di, self%ni
                  do m=1-S, S-1
                     v(2,m+S) = max(a, 0._R8P) * q(1,i+m*di,j+m*dj,k+m*dk,b)
                     v(1,m+S) = min(a, 0._R8P) * q(1,i+(m+1)*di,j+(m+1)*dj,k+(m+1)*dk,b)
                  enddo
                  call weno_reconstruct_upwind(S=S, weno_a=self%weno%a, weno_p=self%weno%p, weno_d=self%weno%d, &
                                               weno_zeps=self%weno%zeps, weno_rmu=1._R8P, v=v, vr=vr)
                  fl(1,i,j,k,b) = vr(1) + vr(2)
               enddo
            enddo
         enddo
      enddo
      endsubroutine face_fluxes
   endsubroutine compute_residuals
```

What comes from the library:
- **The ghost cells.** `update_ghost` calls `field%update_ghost_local` (copies between blocks of the rank, and the 2:1
  transfers) and `field%update_ghost_mpi` (the exchange between ranks). Every face is periodic, so there is no
  physical boundary to apply.
- **The reconstruction.** `weno_reconstruct_upwind` is the library's upwind primitive; the WENO object holds its
  coefficients (`self%weno%a`, `p`, `d`, `zeps`) for the scheme the input names. The application only packs the
  stencil: `v(2,:)` the $a^+u$ values left of the face, `v(1,:)` the $a^-u$ values right of it.
- **The integration.** The Runge–Kutta object of chapter 1 (`initialize_stages`, `compute_stage`, `assign_stage`,
  `update_q`) combines the residuals of each stage with the scheme's SSP coefficients.

The optional `s` and `flux_register` arguments serve the staged path of [chapter 3](./amr); on the fast path they are
absent.

## The keys that select the numerics

In `src/tests/ripple/periodic.ini` (excerpt, comments added here):

```ini
[runge_kutta]
scheme = runge-kutta-ssp-54     ; SSP, 5 stages, 4th order

[weno]
scheme = weno-u-5               ; upwind WENO, S = 3

[grid]
ngc = 3                         ; at least S

[time]
CFL = 0.4
```

## Oracle r1: the order

`src/tests/ripple/check.sh` runs the oracles against the exact solution $u_0(\mathbf x - \mathbf a t)$. Leg `r1` runs
`periodic.ini` at $N = 32, 64, 128$ cells per side on 2 ranks and measures the $L_1$ order:

```bash
./src/tests/ripple/check.sh --legs r1
```

```text
   cells           L1         Linf  order L1 order Linf
    1024   1.8584e-04   5.2719e-04
    4096   5.8085e-06   1.8025e-05     5.000      4.870
   16384   1.8101e-07   5.6842e-07     5.004      4.987
PASS [r1]
```

Order 5 to three digits: at CFL 0.4 the fourth-order time integration is not visible, the spatial error dominates.

::: tip What you learned
- The application writes the residuals; the library supplies the ghost cells, the reconstruction primitive and the
  time integration.
- A WENO scheme is a stencil the application packs and a primitive the library evaluates.
- An exact solution and a ladder of grids give the order, the one check that catches a wrong stencil.
:::

Previous: [1. A realm and the core layer](./core). Next: [3. AMR and the staged layer](./amr).
