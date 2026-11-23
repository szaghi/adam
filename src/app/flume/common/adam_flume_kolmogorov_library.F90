!< ADAM, FLUME local Kolmogorov resolution of a block, the estimator of the AMR marker `mode = 7` (issue #133, phase 2).
module adam_flume_kolmogorov_library
!< ADAM, FLUME local Kolmogorov resolution of a block, the estimator of the AMR marker `mode = 7` (issue #133, phase 2).
!<
!< On a cell the ratio of the grid spacing to the Kolmogorov length,
!<```
!< r = h / eta,   eta = (nu^3 / eps)^(1/4),   eps = 2 nu S_d:S_d,   nu = mu(T) / rho,   h = (prod_d dx_d)^(1/n)
!<```
!< with `S_d = S - tr(S) I / 3` the deviatoric part of the strain rate `S = (grad u + grad u^T) / 2` (the trace over the
!< three components whatever the active axes: a 2-D run is the 2-D section of a 3-D flow), `mu(T)` the viscosity law of
!< the run (`temperature_law`, the same as the viscous fluxes') and `h` the geometric mean of the cell spacings over the
!< `n` active (not null) axes: a null axis has no resolution to measure, so its extent does not enter `h`. Hence
!< `r = h (2 S_d:S_d)^(1/4) / nu^(1/2)`, `r^2 = h^2 |dU/dy| / nu` in a pure shear `u = U(y)` (`2 S_d:S_d = U'^2`).
!< At a wall `eps+ = 1` gives `r = h+`; in the log layer `eps+ ~ 1/(kappa y+)` gives `r = h+ / (kappa y+)^(1/4)`.
!<
!< The velocity gradient is the stencil of the gradient marker (`mark_by_gradient`, issue #76): along each active axis
!< with more than one cell per block, second-order centred differences, across the block faces too (one ghost layer);
!< first-order one-sided at a face on the realm boundary (physical or seam), whose ghosts it never reads. A linear
!< velocity field is differentiated exactly everywhere. A block takes the largest `r` over its interior cells.
!<
!< Host procedures (the markers run on the host on both backends, the FNL one after copying the state back).

! FLUME modules
use :: adam_flume_dissipation_library, only : temperature_law
use :: adam_flume_parameters,          only : IA_R, IA_T, IA_U, IA_V, IA_W
! third party modules
use :: penf,                           only : I4P, R8P

implicit none
private
public :: kolmogorov_block_max
public :: kolmogorov_ratio
public :: kolmogorov_spacing

contains
   pure function kolmogorov_spacing(dxyz, active) result(h)
   !< Return the cell spacing of the estimator, the geometric mean `(prod_d dx_d)^(1/n)` over the `n` active axes (0 if
   !< none is active).
   real(R8P), intent(in) :: dxyz(3)   !< Cell spacings.
   logical,   intent(in) :: active(3) !< Active (not null) axes.
   real(R8P)             :: h         !< Spacing.
   integer(I4P)          :: n         !< Active axes.

   n = count(active)
   h = 0._R8P
   if (n > 0_I4P) h = product(dxyz, mask=active)**(1._R8P / real(n, R8P))
   endfunction kolmogorov_spacing

   pure function kolmogorov_ratio(h, nu, g) result(r)
   !< Return `r = h / eta = h (2 S_d:S_d)^(1/4) / nu^(1/2)` from the spacing `h`, the kinematic viscosity `nu > 0` and the
   !< velocity gradient `g(a,b) = d u_a / d x_b`.
   real(R8P), intent(in) :: h      !< Spacing.
   real(R8P), intent(in) :: nu     !< Kinematic viscosity.
   real(R8P), intent(in) :: g(3,3) !< Velocity gradient.
   real(R8P)             :: r      !< Spacing over the Kolmogorov length.
   real(R8P)             :: s(3,3) !< Strain rate.
   real(R8P)             :: tr     !< Its trace.
   real(R8P)             :: sd2    !< S_d:S_d.

   s = 0.5_R8P * (g + transpose(g))
   tr = s(1,1) + s(2,2) + s(3,3)
   sd2 = max(sum(s * s) - tr * tr / 3._R8P, 0._R8P) ! S_d:S_d = S:S - tr(S)^2 / 3 >= 0, clipped at round-off
   r = h * sqrt(sqrt(2._R8P * sd2)) / sqrt(nu)
   endfunction kolmogorov_ratio

   pure function kolmogorov_block_max(nv, ni, nj, nk, ngc, qa, dxyz, active, is_open, mu0, tref, suth, omega_mu) &
      result(r_max)
   !< Return the largest `r = h / eta` over the interior cells of one block from its auxiliary variables `qa` (density,
   !< velocity, temperature), one ghost layer read across its open faces.
   integer(I4P), intent(in) :: nv                    !< Auxiliary variables number.
   integer(I4P), intent(in) :: ni, nj, nk, ngc       !< Block cells, ghost cells.
   real(R8P),    intent(in) :: qa(nv,1-ngc:ni+ngc,1-ngc:nj+ngc,1-ngc:nk+ngc) !< Auxiliary variables of the block.
   real(R8P),    intent(in) :: dxyz(3)               !< Cell spacings.
   logical,      intent(in) :: active(3)             !< Active (not null) axes.
   logical,      intent(in) :: is_open(2,3)          !< The minus/plus face of each axis is internal (its ghosts valid).
   real(R8P),    intent(in) :: mu0                   !< Viscosity at the reference temperature.
   real(R8P),    intent(in) :: tref, suth, omega_mu  !< Viscosity law (`temperature_law`).
   real(R8P)                :: r_max                 !< Largest r over the interior cells.
   real(R8P)                :: g(3,3)                !< Velocity gradient of a cell, g(a,d) = d u_a / d x_d.
   real(R8P)                :: h                     !< Spacing.
   real(R8P)                :: nu                    !< Kinematic viscosity of a cell.
   logical                  :: diff(3)               !< Axes differentiated: active with more than one cell.
   integer(I4P)             :: n(3)                  !< Cells per axis.
   integer(I4P)             :: c(3)                  !< Cell index.
   integer(I4P)             :: e(3)                  !< Unit step along an axis.
   integer(I4P)             :: i, j, k, d            !< Counters.

   n = [ni, nj, nk]
   diff = active .and. n > 1_I4P
   h = kolmogorov_spacing(dxyz=dxyz, active=active)
   r_max = 0._R8P
   do k=1, nk
      do j=1, nj
         do i=1, ni
            c = [i, j, k]
            g = 0._R8P
            do d=1, 3
               if (.not.diff(d)) cycle
               e = 0_I4P ; e(d) = 1_I4P
               if (c(d) == 1_I4P .and. .not.is_open(1,d)) then
                  g(:,d) = (velocity(c + e) - velocity(c)) / dxyz(d)
               elseif (c(d) == n(d) .and. .not.is_open(2,d)) then
                  g(:,d) = (velocity(c) - velocity(c - e)) / dxyz(d)
               else
                  g(:,d) = 0.5_R8P * (velocity(c + e) - velocity(c - e)) / dxyz(d)
               endif
            enddo
            nu = mu0 * temperature_law(t=qa(IA_T,i,j,k), tref=tref, suth=suth, omega=omega_mu) / qa(IA_R,i,j,k)
            r_max = max(r_max, kolmogorov_ratio(h=h, nu=nu, g=g))
         enddo
      enddo
   enddo
   contains
      pure function velocity(ijk) result(u)
      !< Return the velocity of cell `ijk`.
      integer(I4P), intent(in) :: ijk(3) !< Cell index.
      real(R8P)                :: u(3)   !< Velocity.

      u = [qa(IA_U,ijk(1),ijk(2),ijk(3)), qa(IA_V,ijk(1),ijk(2),ijk(3)), qa(IA_W,ijk(1),ijk(2),ijk(3))]
      endfunction velocity
   endfunction kolmogorov_block_max
endmodule adam_flume_kolmogorov_library
