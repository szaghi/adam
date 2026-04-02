!< ADAM, finite difference/volume operators approximations library, FNL device backend.
!< Same as CPU backend but with variables index placed as last dimension (transposed GPU layout).
module adam_fnl_fdv_operators_library
!< ADAM, finite difference/volume operators approximations library, FNL device backend.

! ADAM modules
use :: adam_fdv_operators_library
! third party modules
use :: penf

implicit none
save
private
! public :: fdv_stencil
! finite difference
public :: compute_curl_fd_centered_dev
public :: compute_divergence_fd_centered_dev
public :: compute_gradient_fd_centered_dev
public :: compute_laplacian_fd_centered_dev
public :: compute_reconstruction_r_fd_centered_dev
! finite volume
! public :: compute_curl_fv_centered_dev
! public :: compute_divergence_fv_centered_dev
! public :: compute_gradient_fv_centered_dev
! public :: compute_laplacian_fv_centered_dev

! real(R8P) :: fdv_stencil(1-S_MAX:1+S_MAX,1-S_MAX:1+S_MAX,1-S_MAX:1+S_MAX,1:3) !< Buffer stencil for acc routine seq.

contains
   ! public methods
   ! finite difference schemes
   pure subroutine compute_curl_fd_centered_dev(s,dxyz,qsx_y,qsx_z,qsy_x,qsy_z,qsz_x,qsz_y,curl)
   !< Compute curl of q vector field with finite difference centered scheme.
   integer(I4P), intent(in)  :: s              !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)       !< Space steps [1:3].
   real(R8P),    intent(in)  :: qsx_y(1-s:)    !< Y component of vector field over the x stencil.
   real(R8P),    intent(in)  :: qsx_z(1-s:)    !< Z component of vector field over the x stencil.
   real(R8P),    intent(in)  :: qsy_x(1-s:)    !< X component of vector field over the y stencil.
   real(R8P),    intent(in)  :: qsy_z(1-s:)    !< Z component of vector field over the y stencil.
   real(R8P),    intent(in)  :: qsz_x(1-s:)    !< X component of vector field over the z stencil.
   real(R8P),    intent(in)  :: qsz_y(1-s:)    !< Y component of vector field over the z stencil.
   real(R8P),    intent(out) :: curl(1:)       !< Curl of q [1:3].
   real(R8P)                 :: dqx_dy, dqx_dz !< Derivatives of qx.
   real(R8P)                 :: dqy_dx, dqy_dz !< Derivatives of qy.
   real(R8P)                 :: dqz_dx, dqz_dy !< Derivatives of qz.
   !$acc routine seq
   !$acc routine(compute_derivative1_fd_centered)

   call compute_derivative1_fd_centered(s=s,ds=dxyz(2),q=qsy_x,dq_ds=dqx_dy)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(3),q=qsz_x,dq_ds=dqx_dz)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(1),q=qsx_y,dq_ds=dqy_dx)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(3),q=qsz_y,dq_ds=dqy_dz)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(1),q=qsx_z,dq_ds=dqz_dx)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(2),q=qsy_z,dq_ds=dqz_dy)
   curl(1) = dqz_dy - dqy_dz
   curl(2) = dqx_dz - dqz_dx
   curl(3) = dqy_dx - dqx_dy
   endsubroutine compute_curl_fd_centered_dev

   pure subroutine compute_divergence_fd_centered_dev(s,dxyz,qsx,qsy,qsz,divergence)
   !< Compute divergence of q vector field with finite difference centered scheme.
   !< FNL device backend: q has variables as last dimension [1-s:1+s,1-s:1+s,1-s:1+s,1:3].
   integer(I4P), intent(in)  :: s                   !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)            !< Space steps [1:3].
   real(R8P),    intent(in)  :: qsx(1-s:)           !< X component of vector field over the x stencil.
   real(R8P),    intent(in)  :: qsy(1-s:)           !< Y component of vector field over the y stencil.
   real(R8P),    intent(in)  :: qsz(1-s:)           !< Z component of vector field over the z stencil.
   real(R8P),    intent(out) :: divergence          !< Divergence of q.
   real(R8P)                 :: div_x, div_y, div_z !< Divergence components.
   !$acc routine seq

   call compute_derivative1_fd_centered(s=s,ds=dxyz(1),q=qsx(1-s:1+s),dq_ds=div_x)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(2),q=qsy(1-s:1+s),dq_ds=div_y)
   call compute_derivative1_fd_centered(s=s,ds=dxyz(3),q=qsz(1-s:1+s),dq_ds=div_z)
   divergence = div_x + div_y + div_z
   endsubroutine compute_divergence_fd_centered_dev

   pure subroutine compute_gradient_fd_centered_dev(s,dxyz,q,gradient)
   !< Compute gradient of q scalar field with finite difference centered scheme.
   integer(I4P), intent(in)  :: s                 !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)          !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:) !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: gradient(1:)      !< Gradient of q [1:3].
   !$acc routine seq

   call compute_derivative1_fd_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),dq_ds=gradient(1))
   call compute_derivative1_fd_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),dq_ds=gradient(2))
   call compute_derivative1_fd_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),dq_ds=gradient(3))
   endsubroutine compute_gradient_fd_centered_dev

   pure subroutine compute_laplacian_fd_centered_dev(s,dxyz,q,laplacian)
   !< Compute laplacian of q scalar field with finite difference centered scheme.
   integer(I4P), intent(in)  :: s                       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)                !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:)       !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: laplacian               !< Laplacian of q.
   real(R8P)                 :: d2q_dx2,d2q_dy2,d2q_dz2 !< Laplacian parts.
   !$acc routine seq

   call compute_derivative2_fd_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),d2q_ds2=d2q_dx2)
   call compute_derivative2_fd_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),d2q_ds2=d2q_dy2)
   call compute_derivative2_fd_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),d2q_ds2=d2q_dz2)
   laplacian = d2q_dx2 + d2q_dy2 + d2q_dz2
   endsubroutine compute_laplacian_fd_centered_dev

   pure subroutine compute_reconstruction_r_fd_centered_dev(s,q,qr)
   !< Compute reconstruction at right interface from cell center average values. Centered schemes.
   integer(I4P), intent(in)  :: s                        !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: q(1-FDV_S_MAX:FDV_S_MAX) !< Scalar field over the stencil [1-s:s].
   real(R8P),    intent(out) :: qr                       !< Reconstruction at right interface of field.
   integer(I4P)              :: m                        !< Counter.
   !$acc routine seq

   qr = 0.0_R8P
   do m=1, s
      qr = qr + FD0_CC(m,s)*(q(m) + q(1-m))
   enddo
   endsubroutine compute_reconstruction_r_fd_centered_dev

   ! finite volume schemes
   ! centered
   pure subroutine compute_curl_fv_centered_dev(s,dxyz,q,curl)
   !< Compute curl of q vector field with finite volume centered scheme.
   !< FNL device backend: q has variables as last dimension [1-s:1+s,1-s:1+s,1-s:1+s,1:3].
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)             !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:,1:) !< Vector field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s,1:3].
   real(R8P),    intent(out) :: curl(1:)             !< Curl of q [1:3].
   real(R8P)                 :: dqx_dy, dqx_dz       !< Derivatives of qx.
   real(R8P)                 :: dqy_dx, dqy_dz       !< Derivatives of qy.
   real(R8P)                 :: dqz_dx, dqz_dy       !< Derivatives of qz.
   !$acc routine seq

   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(1      ,1-s:1+s,1      ,1),dq_ds=dqx_dy)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(1      ,1      ,1-s:1+s,1),dq_ds=dqx_dz)

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1      ,1      ,2),dq_ds=dqy_dx)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(1      ,1      ,1-s:1+s,2),dq_ds=dqy_dz)

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1      ,1      ,3),dq_ds=dqz_dx)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(1      ,1-s:1+s,1      ,3),dq_ds=dqz_dy)

   curl(1) = dqz_dy - dqy_dz
   curl(2) = dqx_dz - dqz_dx
   curl(3) = dqy_dx - dqx_dy
   endsubroutine compute_curl_fv_centered_dev

   pure subroutine compute_divergence_fv_centered_dev(s,dxyz,q,divergence)
   !< Compute divergence of q vector field with finite volume centered scheme.
   !< FNL device backend: q has variables as last dimension [1-s:1+s,1-s:1+s,1-s:1+s,1:3].
   integer(I4P), intent(in)  :: s                    !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)             !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:,1:) !< Vector field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s,1:3].
   real(R8P),    intent(out) :: divergence           !< Divergence of q.
   real(R8P)                 :: div_x, div_y, div_z  !< Divergence components.
   !$acc routine seq

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1,1),dq_ds=div_x)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1,2),dq_ds=div_y)
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s,3),dq_ds=div_z)
   divergence = div_x + div_y + div_z
   endsubroutine compute_divergence_fv_centered_dev

   pure subroutine compute_gradient_fv_centered_dev(s,dxyz,q,gradient)
   !< Compute gradient of q scalar field with finite volume centered scheme.
   integer(I4P), intent(in)  :: s                 !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)          !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:) !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: gradient(1:)      !< Gradient of q [1:3].
   !$acc routine seq

   call compute_derivative1_fv_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),dq_ds=gradient(1))
   call compute_derivative1_fv_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),dq_ds=gradient(2))
   call compute_derivative1_fv_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),dq_ds=gradient(3))
   endsubroutine compute_gradient_fv_centered_dev

   pure subroutine compute_laplacian_fv_centered_dev(s,dxyz,q,laplacian)
   !< Compute laplacian of q scalar field with finite volume centered scheme.
   integer(I4P), intent(in)  :: s                       !< Stencil len, half of accuracy order.
   real(R8P),    intent(in)  :: dxyz(1:)                !< Space steps [1:3].
   real(R8P),    intent(in)  :: q(1-s:,1-s:,1-s:)       !< Scalar field over the stencil [1-s:1+s,1-s:1+s,1-s:1+s].
   real(R8P),    intent(out) :: laplacian               !< Laplacian of q.
   real(R8P)                 :: d2q_dx2,d2q_dy2,d2q_dz2 !< Laplacian parts.
   !$acc routine seq

   call compute_derivative2_fv_centered(s=s,ds=dxyz(1),q=q(1-s:1+s,1,1),d2q_ds2=d2q_dx2)
   call compute_derivative2_fv_centered(s=s,ds=dxyz(2),q=q(1,1-s:1+s,1),d2q_ds2=d2q_dy2)
   call compute_derivative2_fv_centered(s=s,ds=dxyz(3),q=q(1,1,1-s:1+s),d2q_ds2=d2q_dz2)
   laplacian = d2q_dx2 + d2q_dy2 + d2q_dz2
   endsubroutine compute_laplacian_fv_centered_dev
endmodule adam_fnl_fdv_operators_library
