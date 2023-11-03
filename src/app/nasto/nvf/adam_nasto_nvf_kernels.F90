!< ADAM, NASTO NVF application kernels.
module adam_nasto_nvf_kernels
!< ADAM, NASTO NVF application kernels.

use adam_weno_nvf_kernels
use penf, only : I4P, I8P, R8P
use cudafor

implicit none
private
public :: compute_fluxes_convective_kernel
public :: compute_fluxes_difference_cuf
public :: compute_fluxes_diffusive_cuf
public :: compute_q_aux_cuf
public :: compute_q_gradient_cuf
public :: compute_umax_cuf
public :: set_bc_q_gpu_cuf
public :: compute_rk_q_gpu_cuf

contains
   ! public procedures
   attributes(global) subroutine compute_fluxes_convective_kernel(dir,blocks_number,ni,nj,nk,ngc,nv,                &
                                                                  weno_s,weno_a_gpu,weno_p_gpu,weno_d_gpu,weno_zeps,&
                                                                  ror_number,ror_schemes_gpu,ror_threshold,         &
                                                                  ror_ivar_gpu,ror_stats_gpu,                       &
                                                                  g,q_aux_gpu,fluxes_gpu)
   !< Compute convective fluxes along x direction.
   integer(I4P), intent(in),    value  :: dir                                       !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in),    value  :: blocks_number                             !< Number of blocks.
   integer(I4P), intent(in),    value  :: ni                                        !< Grid cells number in I direction.
   integer(I4P), intent(in),    value  :: nj                                        !< Grid cells number in J direction.
   integer(I4P), intent(in),    value  :: nk                                        !< Grid cells number in K direction.
   integer(I4P), intent(in),    value  :: ngc                                       !< Ghost cells number.
   integer(I4P), intent(in),    value  :: nv                                        !< Number of conservative varibales.
   integer(I4P), intent(in),    value  :: weno_s                                    !< Weno stencils number/dimension.
   real(R8P),    intent(in),    device :: weno_a_gpu(1:,0:,1:)                      !< Optimal weights.
   real(R8P),    intent(in),    device :: weno_p_gpu(1:,0:,0:,1:)                   !< Polinomials coefficients.
   real(R8P),    intent(in),    device :: weno_d_gpu(0:,0:,0:,1:)                   !< Smoothness indicators coefficients.
   real(R8P),    intent(in),    value  :: weno_zeps                                 !< Parameter to avoid division by zero.
   integer(I4P), intent(in),    value  :: ror_number                                !< Number of ROR iterations.
   integer(I4P), intent(in),    device :: ror_schemes_gpu(1:)                       !< Scheme (S value) for each ROR step.
   real(R8P),    intent(in),    value  :: ror_threshold                             !< ROR threshold triggering.
   integer(I4P), intent(in),    device :: ror_ivar_gpu(1:)                          !< Index variables to check in ROR.
   integer(I4P), intent(inout), device :: ror_stats_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Scheme (S value) for each ROR step.
   real(R8P),    intent(in),    value  :: g                                         !< Specific heats ratio.
   real(R8P),    intent(in),    device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Auxiliary variables.
   real(R8P),    intent(inout), device :: fluxes_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Fluxes.
   integer(I4P)                        :: b, i, j, k                                !< Counter.
   integer(I4P)                        :: si(3)                                     !< Stencil increment.
   real(R8P)                           :: sir(3)                                    !< Stencil increment, real cast.

   select case(dir)
   case(1)
      b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
      j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
      if (b>blocks_number.or.j>nj) return
      si = [1,0,0]
      sir = real(si,R8P)
      do k=1-si(3), nk
      do i=1-si(1), ni
         call compute_fluxes_conv_interface_kernel(dir=dir,si=si,sir=sir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,       &
                                                   weno_s=weno_s,weno_a_gpu=weno_a_gpu,weno_p_gpu=weno_p_gpu, &
                                                   weno_d_gpu=weno_d_gpu,weno_zeps=weno_zeps,                 &
                                                   ror_number=ror_number, ror_schemes_gpu=ror_schemes_gpu,    &
                                                   ror_threshold=ror_threshold, ror_ivar_gpu=ror_ivar_gpu,    &
                                                   ror_stats_gpu=ror_stats_gpu,                               &
                                                   g=g,q_aux_gpu=q_aux_gpu,fluxes_gpu=fluxes_gpu)
      enddo
      enddo
   case(2)
      b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
      i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
      if (b>blocks_number.or.i>ni) return
      si = [0,1,0]
      sir = real(si,R8P)
      do k=1-si(3), nk
      do j=1-si(2), nj
         call compute_fluxes_conv_interface_kernel(dir=dir,si=si,sir=sir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,       &
                                                   weno_s=weno_s,weno_a_gpu=weno_a_gpu,weno_p_gpu=weno_p_gpu, &
                                                   weno_d_gpu=weno_d_gpu,weno_zeps=weno_zeps,                 &
                                                   ror_number=ror_number, ror_schemes_gpu=ror_schemes_gpu,    &
                                                   ror_threshold=ror_threshold, ror_ivar_gpu=ror_ivar_gpu,    &
                                                   ror_stats_gpu=ror_stats_gpu,                               &
                                                   g=g,q_aux_gpu=q_aux_gpu,fluxes_gpu=fluxes_gpu)
      enddo
      enddo
   case(3)
      b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
      i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
      if (b>blocks_number.or.i> ni) return
      si = [0,0,1]
      sir = real(si,R8P)
      do j=1-si(2), nj
      do k=1-si(3), nk
         call compute_fluxes_conv_interface_kernel(dir=dir,si=si,sir=sir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,       &
                                                   weno_s=weno_s,weno_a_gpu=weno_a_gpu,weno_p_gpu=weno_p_gpu, &
                                                   weno_d_gpu=weno_d_gpu,weno_zeps=weno_zeps,                 &
                                                   ror_number=ror_number, ror_schemes_gpu=ror_schemes_gpu,    &
                                                   ror_threshold=ror_threshold, ror_ivar_gpu=ror_ivar_gpu,    &
                                                   ror_stats_gpu=ror_stats_gpu,                               &
                                                   g=g,q_aux_gpu=q_aux_gpu,fluxes_gpu=fluxes_gpu)
      enddo
      enddo
   endselect
   endsubroutine compute_fluxes_convective_kernel

   subroutine compute_fluxes_difference_cuf(blocks_number, ni, nj, nk, ngc, nv, ib_eps, &
                                            dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, fl_gpu)
   !< Compute fluxes difference.
   integer(I4P), intent(in)                      :: blocks_number                       !< Number of blocks.
   integer(I4P), intent(in)                      :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)                      :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)                      :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)                      :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)                      :: nv                                  !< Number of conservative varibales.
   real(R8P),    intent(in)                      :: ib_eps                              !< Tolerance IB delta ratio.
   real(R8P),    intent(in),    device           :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)  !< Space steps.
   real(R8P),    intent(in),    device           :: flx_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< X direction fluxes.
   real(R8P),    intent(in),    device           :: fly_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Y direction fluxes.
   real(R8P),    intent(in),    device           :: flz_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Z direction fluxes.
   real(R8P),    intent(in),    device, optional :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance function.
   real(R8P),    intent(inout), device           ::  fl_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes differences.
   real(R8P)                                     :: delta_x, delta_y, delta_z           !< Space steps.
   real(R8P)                                     :: dx_locale, dy_locale, dz_locale     !< Local space steps.
   integer(I4P)                                  :: b, i, j, k, v                       !< Counter.
   integer(I4P)                                  :: iercuda                             !< Error trapping flag for CUDAFortran.
   integer(I4P)                                  :: all_solids                          !< Last phi index, all solids summary.

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$cuf kernel do(4) <<<*,*>>>
      do k=1,nk
      do j=1,nj
      do i=1,ni
      do b=1,blocks_number
         dx_locale = dx_gpu(b)
         if (phi_gpu(b,i,j,k,all_solids)<0.) then
            if (phi_gpu(b,i+1,j,k,all_solids)*phi_gpu(b,i-1,j,k,all_solids)<0) then
               if (phi_gpu(b,i+1,j,k,all_solids)>0.) then
                  delta_x= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i+1,j,k,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dx_gpu(b)
                  dx_locale = dx_gpu(b)/2 + delta_x
               else
                  delta_x= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i-1,j,k,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dx_gpu(b)
                  dx_locale = dx_gpu(b)/2 + delta_x
               endif
            endif
         endif
         dy_locale = dy_gpu(b)
         if (phi_gpu(b,i,j,k,all_solids)<0.) then
            if (phi_gpu(b,i,j+1,k,all_solids)*phi_gpu(b,i,j-1,k,all_solids)<0) then
               if (phi_gpu(b,i,j+1,k,all_solids)>0.) then
                  delta_y= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i,j+1,k,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dy_gpu(b)
                  dy_locale = dy_gpu(b)/2 + delta_y
               else
                  delta_y= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i,j-1,k,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dy_gpu(b)
                  dy_locale = dy_gpu(b)/2 + delta_y
               endif
            endif
         endif
         dz_locale = dz_gpu(b)
         if (phi_gpu(b,i,j,k,all_solids)<0.) then
            if (phi_gpu(b,i,j,k+1,all_solids)*phi_gpu(b,i,j,k-1,all_solids)<0) then
               if (phi_gpu(b,i,j,k+1,all_solids)>0.) then
                  delta_z= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i,j,k+1,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dz_gpu(b)
                  dz_locale = dz_gpu(b)/2 + delta_z
               else
                  delta_z= -phi_gpu(b,i,j,k,all_solids)/(phi_gpu(b,i,j,k-1,all_solids)-phi_gpu(b,i,j,k,all_solids)+ib_eps)*dz_gpu(b)
                  dz_locale = dz_gpu(b)/2 + delta_z
               endif
            endif
         endif
         do v=1, nv
            fl_gpu(b,i,j,k,v) = - (flx_gpu(b,i,j,k,v)-flx_gpu(b,i-1,j,k,v))/dx_locale &
                                - (fly_gpu(b,i,j,k,v)-fly_gpu(b,i,j-1,k,v))/dy_locale &
                                - (flz_gpu(b,i,j,k,v)-flz_gpu(b,i,j,k-1,v))/dz_locale
         enddo
      enddo
      enddo
      enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   else
      !$cuf kernel do(4) <<<*,*>>>
      do k=1,nk
      do j=1,nj
      do i=1,ni
      do b=1,blocks_number
         do v=1, nv
            fl_gpu(b,i,j,k,v) = - (flx_gpu(b,i,j,k,v)-flx_gpu(b,i-1,j,k,v))/dx_gpu(b) &
                                - (fly_gpu(b,i,j,k,v)-fly_gpu(b,i,j-1,k,v))/dy_gpu(b) &
                                - (flz_gpu(b,i,j,k,v)-flz_gpu(b,i,j,k-1,v))/dz_gpu(b)
         enddo
      enddo
      enddo
      enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endsubroutine compute_fluxes_difference_cuf

   subroutine compute_fluxes_diffusive_cuf(blocks_number, ni, nj, nk, ngc, nv, mu, kd, &
                                           q_aux_gpu, dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu)
   !< Compute diffusive fluxes.
   integer(I4P), intent(in)            :: blocks_number                             !< Blocks number.
   integer(I4P), intent(in)            :: ni, nj, nk                                !< Grid dimensionns.
   integer(I4P), intent(in)            :: ngc                                       !< Number of ghost cells.
   integer(I4P), intent(in)            :: nv                                        !< Number of conservative variables.
   real(R8P),    intent(in)            :: mu                                        !< Viscosity.
   real(R8P),    intent(in)            :: kd                                        !< Thermal diffusivity.
   real(R8P),    intent(in),    device :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)        !< Space steps.
   real(R8P),    intent(in),    device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:) !< Auxiliary varibales
   real(R8P),    intent(inout), device ::   flx_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:) !< Fluxes along x.
   real(R8P),    intent(inout), device ::   fly_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:) !< Fluxes along y.
   real(R8P),    intent(inout), device ::   flz_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:) !< Fluxes along z.
   real(R8P)                           :: vel_u, vel_v, vel_w                       !< (Mean) velocity.
   integer(I4P)                        :: b, i, j, k, v                             !< Counter.
   real(R8P)                           :: du_dx, dv_dx, dw_dx                       !< Velocity derivative along x.
   real(R8P)                           :: du_dy, dv_dy, dw_dy                       !< Velocity derivative along y.
   real(R8P)                           :: du_dz, dv_dz, dw_dz                       !< Velocity derivative along z.
   real(R8P)                           :: sigq, sigl                                !< Sigmas.
   real(R8P)                           :: tau_1_1, tau_2_1, tau_3_1, dT_dx          !< Stress tensor, x elements.
   real(R8P)                           :: tau_1_2, tau_2_2, tau_3_2, dT_dy          !< Stress tensor, y elements.
   real(R8P)                           :: tau_1_3, tau_2_3, tau_3_3, dT_dz          !< Stress tensor, z elements.
   integer(I4P)                        :: iercuda                                   !< CUDA error trapping flag.

   !$cuf kernel do(3) <<<*,*>>>
   do k=1,nk
      do j=1,nj
         do b=1,blocks_number
            do i=0,ni ! loop on interfaces
                du_dx = (q_aux_gpu(b,i+1,j,k,2)-q_aux_gpu(b,i,j,k,2))/dx_gpu(b)
                dv_dx = (q_aux_gpu(b,i+1,j,k,3)-q_aux_gpu(b,i,j,k,3))/dx_gpu(b)
                dw_dx = (q_aux_gpu(b,i+1,j,k,4)-q_aux_gpu(b,i,j,k,4))/dx_gpu(b)

                du_dy = (q_aux_gpu(b,i+1,j+1,k,2) - q_aux_gpu(b,i+1,j-1,k,2)+ &
                         q_aux_gpu(b,i,j+1,k,2)   - q_aux_gpu(b,i,j-1,k,2) )*0.25_R8P/dy_gpu(b)
                dv_dy = (q_aux_gpu(b,i+1,j+1,k,3) - q_aux_gpu(b,i+1,j-1,k,3)+ &
                         q_aux_gpu(b,i,j+1,k,3)   - q_aux_gpu(b,i,j-1,k,3) )*0.25_R8P/dy_gpu(b)

                du_dz = (q_aux_gpu(b,i+1,j,k+1,2) - q_aux_gpu(b,i+1,j,k-1,2)+ &
                         q_aux_gpu(b,i,j,k+1,2)   - q_aux_gpu(b,i,j,k-1,2) )*0.25_R8P/dz_gpu(b)
                dw_dz = (q_aux_gpu(b,i+1,j,k+1,4) - q_aux_gpu(b,i+1,j,k-1,4)+ &
                         q_aux_gpu(b,i,j,k+1,4)   - q_aux_gpu(b,i,j,k-1,4) )*0.25_R8P/dz_gpu(b)

                vel_u = 0.5*(q_aux_gpu(b,i,j,k,2) + q_aux_gpu(b,i+1,j,k,2))
                vel_v = 0.5*(q_aux_gpu(b,i,j,k,3) + q_aux_gpu(b,i+1,j,k,3))
                vel_w = 0.5*(q_aux_gpu(b,i,j,k,4) + q_aux_gpu(b,i+1,j,k,4))

                tau_1_1 = 2.0*mu*(du_dx-1./3.*(du_dx+dv_dy+dw_dz))
                tau_2_1 = mu*(dv_dx+du_dy)
                tau_3_1 = mu*(dw_dx+du_dz)

                dT_dx = (q_aux_gpu(b,i+1,j,k,6)-q_aux_gpu(b,i,j,k,6))/dx_gpu(b)

                sigq = kd*dT_dx
                sigl = vel_u*tau_1_1+vel_v*tau_2_1+vel_w*tau_3_1

                flx_gpu(b,i,j,k,2) = flx_gpu(b,i,j,k,2) - tau_1_1
                flx_gpu(b,i,j,k,3) = flx_gpu(b,i,j,k,3) - tau_2_1
                flx_gpu(b,i,j,k,4) = flx_gpu(b,i,j,k,4) - tau_3_1
                flx_gpu(b,i,j,k,5) = flx_gpu(b,i,j,k,5) - sigq + sigl
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   !$cuf kernel do(3) <<<*,*>>>
   do k=1,nk
      do i=1,ni
         do b=1,blocks_number
            do j=0,nj ! loop on interfaces
                du_dy = (q_aux_gpu(b,i,j+1,k,2)-q_aux_gpu(b,i,j,k,2))/dy_gpu(b)
                dv_dy = (q_aux_gpu(b,i,j+1,k,3)-q_aux_gpu(b,i,j,k,3))/dy_gpu(b)
                dw_dy = (q_aux_gpu(b,i,j+1,k,4)-q_aux_gpu(b,i,j,k,4))/dy_gpu(b)

                du_dx = (q_aux_gpu(b,i+1,j+1,k,2) - q_aux_gpu(b,i-1,j+1,k,2)+ &
                         q_aux_gpu(b,i+1,j,k,2)   - q_aux_gpu(b,i-1,j,k,2) )*0.25_R8P/dx_gpu(b)
                dv_dx = (q_aux_gpu(b,i+1,j+1,k,3) - q_aux_gpu(b,i-1,j+1,k,3)+ &
                         q_aux_gpu(b,i+1,j,k,3)   - q_aux_gpu(b,i-1,j,k,3) )*0.25_R8P/dx_gpu(b)

                dv_dz = (q_aux_gpu(b,i+1,j,k+1,3) - q_aux_gpu(b,i-1,j,k+1,3)+ &
                         q_aux_gpu(b,i+1,j,k,3)   - q_aux_gpu(b,i-1,j,k,3) )*0.25_R8P/dz_gpu(b)
                dw_dz = (q_aux_gpu(b,i+1,j,k+1,4) - q_aux_gpu(b,i-1,j,k+1,4)+ &
                         q_aux_gpu(b,i+1,j,k,4)   - q_aux_gpu(b,i-1,j,k,4) )*0.25_R8P/dz_gpu(b)

                vel_u = 0.5*(q_aux_gpu(b,i,j,k,2) + q_aux_gpu(b,i,j+1,k,2))
                vel_v = 0.5*(q_aux_gpu(b,i,j,k,3) + q_aux_gpu(b,i,j+1,k,3))
                vel_w = 0.5*(q_aux_gpu(b,i,j,k,4) + q_aux_gpu(b,i,j+1,k,4))

                tau_1_2 = mu*(du_dy+dv_dx)
                tau_2_2 = 2.0*mu*(dv_dy-1./3.*(du_dx+dv_dy+dw_dz))
                tau_3_2 = mu*(dw_dy+dv_dz)

                dT_dy = (q_aux_gpu(b,i,j+1,k,6)-q_aux_gpu(b,i,j,k,6))/dy_gpu(b)

                sigq = kd*dT_dy
                sigl = vel_u*tau_1_2+vel_v*tau_2_2+vel_w*tau_3_2

                fly_gpu(b,i,j,k,2) = fly_gpu(b,i,j,k,2) - tau_1_2
                fly_gpu(b,i,j,k,3) = fly_gpu(b,i,j,k,3) - tau_2_2
                fly_gpu(b,i,j,k,4) = fly_gpu(b,i,j,k,4) - tau_3_2
                fly_gpu(b,i,j,k,5) = fly_gpu(b,i,j,k,5) - sigq + sigl
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   !$cuf kernel do(3) <<<*,*>>>
   do j=1,nj
      do i=1,ni
         do b=1,blocks_number
            do k=0,nk ! loop on interfaces
                du_dz = (q_aux_gpu(b,i,j,k+1,2)-q_aux_gpu(b,i,j,k,2))/dz_gpu(b)
                dv_dz = (q_aux_gpu(b,i,j,k+1,3)-q_aux_gpu(b,i,j,k,3))/dz_gpu(b)
                dw_dz = (q_aux_gpu(b,i,j,k+1,4)-q_aux_gpu(b,i,j,k,4))/dz_gpu(b)

                du_dx = (q_aux_gpu(b,i+1,j,k+1,2) - q_aux_gpu(b,i-1,j,k+1,2)+ &
                         q_aux_gpu(b,i+1,j,k,2)   - q_aux_gpu(b,i-1,j,k,2) )*0.25_R8P/dx_gpu(b)
                dw_dx = (q_aux_gpu(b,i+1,j,k+1,4) - q_aux_gpu(b,i-1,j,k+1,4)+ &
                         q_aux_gpu(b,i+1,j,k,4)   - q_aux_gpu(b,i-1,j,k,4) )*0.25_R8P/dx_gpu(b)

                dv_dy = (q_aux_gpu(b,i,j+1,k+1,3) - q_aux_gpu(b,i,j-1,k+1,3)+ &
                         q_aux_gpu(b,i,j+1,k,3)   - q_aux_gpu(b,i,j-1,k,3) )*0.25_R8P/dy_gpu(b)
                dw_dy = (q_aux_gpu(b,i,j+1,k+1,4) - q_aux_gpu(b,i,j-1,k+1,4)+ &
                         q_aux_gpu(b,i,j+1,k,4)   - q_aux_gpu(b,i,j-1,k,4) )*0.25_R8P/dy_gpu(b)

                vel_u = 0.5*(q_aux_gpu(b,i,j,k,2) + q_aux_gpu(b,i,j,k+1,2))
                vel_v = 0.5*(q_aux_gpu(b,i,j,k,3) + q_aux_gpu(b,i,j,k+1,3))
                vel_w = 0.5*(q_aux_gpu(b,i,j,k,4) + q_aux_gpu(b,i,j,k+1,4))

                tau_1_3 = mu*(du_dz+dw_dx)
                tau_2_3 = mu*(dv_dz+dw_dy)
                tau_3_3 = 2.0*mu*(dw_dz-1./3.*(du_dx+dv_dy+dw_dz))

                dT_dz = (q_aux_gpu(b,i,j,k+1,6)-q_aux_gpu(b,i,j,k,6))/dz_gpu(b)

                sigq = kd*dT_dz
                sigl = vel_u*tau_1_3+vel_v*tau_2_3+vel_w*tau_3_3

                flz_gpu(b,i,j,k,2) = flz_gpu(b,i,j,k,2) - tau_1_3
                flz_gpu(b,i,j,k,3) = flz_gpu(b,i,j,k,3) - tau_2_3
                flz_gpu(b,i,j,k,4) = flz_gpu(b,i,j,k,4) - tau_3_3
                flz_gpu(b,i,j,k,5) = flz_gpu(b,i,j,k,5) - sigq + sigl
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_fluxes_diffusive_cuf

   subroutine compute_q_aux_cuf(ni, nj, nk, ngc, ns, blocks_number, R, cv, g, dha, q_gpu, q_aux_gpu)
   !< Compute auxiliary variables.
   integer(I4P), intent(in)          :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)          :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)          :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)          :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)          :: ns                                     !< Number of fluid species.
   integer(I4P), intent(in)          :: blocks_number                          !< Number of blocks.
   real(R8P),    intent(in)          :: R                                      !< Fluid constant, specific heats difference.
   real(R8P),    intent(in)          :: cv                                     !< Specific heat at constant volume.
   real(R8P),    intent(in)          :: g                                      !< Specific heats ratio.
   real(R8P),    intent(in)          :: dha                                    !< Entalpy fluid.
   real(R8P),    intent(in),  device ::     q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative variables.
   real(R8P),    intent(out), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Auxiliary variables.
   integer(I4P)                      :: b, i, j, k, s                          !< Counter.
   real(R8P)                         :: rho, uuu, vvv, www, rhe, rya, yya, tem !< State variables.
   integer(I4P)                      :: iercuda                                !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do k=1-ngc, nk+ngc
      do j=1-ngc, nj+ngc
         do i=1-ngc, ni+ngc
            do b=1, blocks_number
               rho = q_gpu(b,i,j,k,1)
               uuu = q_gpu(b,i,j,k,2)/rho
               vvv = q_gpu(b,i,j,k,3)/rho
               www = q_gpu(b,i,j,k,4)/rho
               rhe = q_gpu(b,i,j,k,5)
               if (ns==2) then
                   rya = q_gpu(b,i,j,k,ns+4)
               else
                   rya = 0._R8P
               endif
               yya = rya/rho
               tem = ((rhe-rya*dha)/rho-0.5*(uuu**2+vvv**2+www**2))/cv

               q_aux_gpu(b,i,j,k,1) = rho           ! density
               q_aux_gpu(b,i,j,k,2) = uuu           ! velocity x
               q_aux_gpu(b,i,j,k,3) = vvv           ! velocity y
               q_aux_gpu(b,i,j,k,4) = www           ! velocity z
               q_aux_gpu(b,i,j,k,5) = yya           ! mass fraction
               q_aux_gpu(b,i,j,k,6) = tem           ! temperature
               q_aux_gpu(b,i,j,k,7) = R*rho*tem     ! pressure
               q_aux_gpu(b,i,j,k,8) = rhe/rho+R*tem ! entalpy
               q_aux_gpu(b,i,j,k,9) = sqrt(g*R*tem) ! sound speed
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_q_aux_cuf

   subroutine compute_q_gradient_cuf(b, ni, nj, nk, ngc, dx, dy, dz, q_gpu, ivar, gradient)
   !< Compute gradient of q(ivar).
   integer(I4P), intent(in)         :: b                                 !< Block index.
   integer(I4P), intent(in)         :: ni                                !< Grid cells number in I direction.
   integer(I4P), intent(in)         :: nj                                !< Grid cells number in J direction.
   integer(I4P), intent(in)         :: nk                                !< Grid cells number in K direction.
   integer(I4P), intent(in)         :: ngc                               !< Ghost cells number.
   real(R8P),    intent(in)         :: dx                                !< X space step.
   real(R8P),    intent(in)         :: dy                                !< Y space step.
   real(R8P),    intent(in)         :: dz                                !< Z space step.
   real(R8P),    intent(in), device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field component to which apply gradient.
   integer(I4P), intent(in)         :: ivar                              !< Ghost cells number.
   real(R8P),    intent(out)        :: gradient                          !< Maximum gradient of q.
   real(R8P)                        :: grad                              !< Current gradient of q.
   integer(I4P)                     :: i, j, k                           !< Counter.
   integer(I4P)                     :: iercuda                           !< Error trapping flag for CUDAFortran.
   real(R8P), parameter             :: tol=1.e-12                        !< Gradient denominator tolerance.

   gradient = 0._R8P
   !$cuf kernel do(3) <<<*,*>>> reduce(max:gradient)
   do k=1, nk
      do j=1, nj
         do i=1, ni
            grad = sqrt(((q_gpu(b,i+1,j,k,ivar) - q_gpu(b,i-1,j,k,ivar))/(2*dx))**2 + &
                        ((q_gpu(b,i,j+1,k,ivar) - q_gpu(b,i,j-1,k,ivar))/(2*dy))**2 + &
                        ((q_gpu(b,i,j,k+1,ivar) - q_gpu(b,i,j,k-1,ivar))/(2*dz))**2)
            grad = grad/(abs(q_gpu(b,i,j,k,ivar))+tol)
            gradient = max(gradient, grad)
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_q_gradient_cuf

   subroutine compute_umax_cuf(ni, nj, nk, ngc, blocks_number, mu, dx_gpu, dy_gpu, dz_gpu, q_aux_gpu, umax)
   !< Compute maximum speed.
   integer(I4P), intent(in)         :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)         :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)         :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)         :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)         :: blocks_number                         !< Number of blocks.
   real(R8P),    intent(in)         :: mu                                    !< Dynamic viscosity.
   real(R8P),    intent(in), device :: dx_gpu(1:)                            !< X space step.
   real(R8P),    intent(in), device :: dy_gpu(1:)                            !< Y space step.
   real(R8P),    intent(in), device :: dz_gpu(1:)                            !< Z space step.
   real(R8P),    intent(in), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary varibales.
   real(R8P),    intent(out)        :: umax                                  !< Maximum speed.
   real(R8P)                        :: ss                                    !< Speed of sound.
   integer(I4P)                     :: i, j, k, b                            !< Counter.
   integer(I4P)                     :: iercuda                               !< Error trapping flag for CUDAFortran.
   real(R8P)                        :: dx_locale, dy_locale, dz_locale       !< Local space steps.

   umax = 0._R8P
   !$cuf kernel do(4) <<<*,*>>> reduce(max:umax)
   do k=1, nk
      do j=1, nj
         do i=1, ni
            do b=1, blocks_number
               dx_locale = dx_gpu(b)*0.5_R8P
               dy_locale = dy_gpu(b)*0.5_R8P
               dz_locale = dz_gpu(b)*0.5_R8P
               ss = q_aux_gpu(b,i,j,k,9)
               umax = max(umax, (abs(q_aux_gpu(b,i,j,k,2)) + ss)/dx_locale + 2._R8P*mu/(q_aux_gpu(b,i,j,k,1))/dx_locale**2 + &
                                (abs(q_aux_gpu(b,i,j,k,3)) + ss)/dy_locale + 2._R8P*mu/(q_aux_gpu(b,i,j,k,1))/dy_locale**2 + &
                                (abs(q_aux_gpu(b,i,j,k,4)) + ss)/dz_locale + 2._R8P*mu/(q_aux_gpu(b,i,j,k,1))/dz_locale**2)
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_umax_cuf

   subroutine set_bc_q_gpu_cuf(BC_EXTRAPOLATION, BC_INFLOW, nv, ngc, cv, R, &
                               local_map_bc_gpu, fec_1_6_array_gpu, q_bc_vars_gpu, q_gpu)
   !< Set BC over q.
   integer(I4P), intent(in)            :: BC_EXTRAPOLATION        !< Extrapolation BC parameter.
   integer(I4P), intent(in)            :: BC_INFLOW               !< Inflow BC parameter.
   integer(I4P), intent(in)            :: nv                      !< Number of variables.
   integer(I4P), intent(in)            :: ngc                     !< Ghost cells number.
   real(R8P),    intent(in)            :: cv                      !< Constant volume specific heat.
   real(R8P),    intent(in)            :: R                       !< Gas constant.
   integer(I8P), intent(in),    device :: local_map_bc_gpu(:,:,:) !< Local map for BC ghost cells.
   integer(I4P), intent(in),    device :: fec_1_6_array_gpu(:)    !< Local map for BC ghost cells.
   real(R8P),    intent(in),    device :: q_bc_vars_gpu(:,:)      !< Boundary variables.
   real(R8P),    intent(inout), device :: q_gpu(1:,    &
                                                1-ngc:,&
                                                1-ngc:,&
                                                1-ngc:,1:)        !< Conservative variables.
   integer(I4P)                        :: b                       !< Counter.
   integer(I4P)                        :: c, i, j, k, v           !< Counter.
   integer(I4P)                        :: idelta                  !< IJK delta step for extrapolation.
   integer(I4P)                        :: jdelta                  !< IJK delta step for extrapolation.
   integer(I4P)                        :: kdelta                  !< IJK delta step for extrapolation.
   integer(I4P)                        :: bc_type                 !< Boundary condition type.
   integer(I4P)                        :: crown                   !< Crown counter.
   integer(I4P)                        :: iercuda                 !< Error trapping flag for CUDAFortran.
   integer(I4P)                        :: fec                     !< Boundary fec (1 to 26).
   integer(I4P)                        :: fec_1_6                 !< Boundary fec (1 to 6).

   do crown=1, ngc
      !$cuf kernel do(1) <<<*,*>>>
      do c=1, size(local_map_bc_gpu, dim=1)
         b = local_map_bc_gpu(c, 1 ,crown)
         if (b>0) then
            i       = local_map_bc_gpu(c, 2 ,crown)
            j       = local_map_bc_gpu(c, 3 ,crown)
            k       = local_map_bc_gpu(c, 4 ,crown)
            idelta  = local_map_bc_gpu(c, 5 ,crown)
            jdelta  = local_map_bc_gpu(c, 6 ,crown)
            kdelta  = local_map_bc_gpu(c, 7 ,crown)
            bc_type = local_map_bc_gpu(c, 8 ,crown)
            fec     = local_map_bc_gpu(c, 9 ,crown)
            fec_1_6 = fec_1_6_array_gpu(fec)
            if (bc_type == BC_EXTRAPOLATION) then
               do v=1, nv
                  q_gpu(b,i,j,k,v) = q_gpu(b,i-idelta,j-jdelta,k-kdelta,v)
               enddo
            else if (bc_type == BC_INFLOW) then
                q_gpu(b,i,j,k,1) = q_bc_vars_gpu(1, fec_1_6)
                q_gpu(b,i,j,k,2) = q_bc_vars_gpu(1, fec_1_6)* q_bc_vars_gpu(2, fec_1_6)
                q_gpu(b,i,j,k,3) = q_bc_vars_gpu(1, fec_1_6)* q_bc_vars_gpu(3, fec_1_6)
                q_gpu(b,i,j,k,4) = q_bc_vars_gpu(1, fec_1_6)* q_bc_vars_gpu(4, fec_1_6)
                q_gpu(b,i,j,k,5) = q_bc_vars_gpu(1, fec_1_6)*                                   &
                                   (cv*q_bc_vars_gpu(5, fec_1_6)/(q_bc_vars_gpu(1, fec_1_6)*R)+ &
                    0.5_R8P*(q_bc_vars_gpu(2, fec_1_6)**2+q_bc_vars_gpu(3, fec_1_6)**2+q_bc_vars_gpu(4, fec_1_6)**2))
            endif
         endif
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   enddo
   endsubroutine set_bc_q_gpu_cuf

   subroutine compute_rk_q_gpu_cuf(ni,nj,nk,ngc,nv,blocks_number,dt,q_gpu,q_old_gpu,fl_gpu,phi_gpu,ark,brk,crk)
   !< Compute RK approximation over q.
   integer(I4P), intent(in)                      :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)                      :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)                      :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)                      :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)                      :: nv                                    !< Number of conservative varibales.
   integer(I4P), intent(in)                      :: blocks_number                         !< Number of blocks.
   real(R8P),    intent(in)                      :: dt                                    !< Time step.
   real(R8P),    intent(in)                      :: ark, brk, crk                         !< Time step.
   real(R8P),    intent(inout), device           ::     q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P),    intent(in),    device           ::    fl_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P),    intent(in),    device, optional ::   phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P),    intent(in),    device           :: q_old_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< RK stage.
   integer(I4P)                                  :: all_solids                            !< Last phi index, all solids summary.
   integer(I4P)                                  :: i, j, k, b, v                         !< Counter.
   integer(I4P)                                  :: iercuda                               !< Error trapping flag for CUDAFortran.

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$cuf kernel do(5) <<<*,*>>>
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     if (phi_gpu(b,i,j,k,all_solids) < 0.) then
                        q_gpu(b,i,j,k,v) = ark * q_old_gpu(b,i,j,k,v) + brk * q_gpu(b,i,j,k,v) + dt * crk * fl_gpu(b,i,j,k,v)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   else
      !$cuf kernel do(5) <<<*,*>>>
      do v=1, nv
         do k=1, nk
            do j=1, nj
               do i=1, ni
                  do b=1, blocks_number
                     q_gpu(b,i,j,k,v) = ark * q_old_gpu(b,i,j,k,v) + brk * q_gpu(b,i,j,k,v) + dt * crk * fl_gpu(b,i,j,k,v)
                  enddo
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endsubroutine compute_rk_q_gpu_cuf

   ! private procedures
   attributes(device) subroutine compute_fluxes_conv_interface_kernel(dir,si,sir,b,i,j,k,ngc,nv,                        &
                                                                      weno_s,weno_a_gpu,weno_p_gpu,weno_d_gpu,weno_zeps,&
                                                                      ror_number,ror_schemes_gpu,ror_threshold,         &
                                                                      ror_ivar_gpu,ror_stats_gpu,                       &
                                                                      g,q_aux_gpu,fluxes_gpu)
   !< Compute convective fluxes at right interface of b,i,j,k.
   integer(I4P), intent(in),    value  :: dir                                       !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)            :: si(3)                                     !< Stencil increment.
   real(R8P)   , intent(in)            :: sir(3)                                    !< Stencil increment, real cast.
   integer(I4P), intent(in),    value  :: b, i, j, k                                !< Counter.
   integer(I4P), intent(in),    value  :: ngc                                       !< Ghost cells number.
   integer(I4P), intent(in),    value  :: nv                                        !< Number of conservative varibales.
   integer(I4P), intent(in),    value  :: weno_s                                    !< Weno stencils number/dimension.
   real(R8P),    intent(in),    device :: weno_a_gpu(1:,0:,1:)                      !< Optimal weights.
   real(R8P),    intent(in),    device :: weno_p_gpu(1:,0:,0:,1:)                   !< Polinomials coefficients.
   real(R8P),    intent(in),    device :: weno_d_gpu(0:,0:,0:,1:)                   !< Smoothness indicators coefficients.
   real(R8P),    intent(in),    value  :: weno_zeps                                 !< Parameter to avoid division by zero.
   integer(I4P), intent(in),    value  :: ror_number                                !< Number of ROR iterations.
   integer(I4P), intent(in),    device :: ror_schemes_gpu(1:)                       !< Scheme (S value) for each ROR step.
   real(R8P),    intent(in),    value  :: ror_threshold                             !< ROR threshold triggering.
   integer(I4P), intent(in),    device :: ror_ivar_gpu(1:)                          !< Index variables to check in ROR.
   integer(I4P), intent(inout), device :: ror_stats_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Scheme (S value) for each ROR step.
   real(R8P),    intent(in),    value  :: g                                         !< Specific heats ratio.
   real(R8P),    intent(in),    device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Auxiliary variables.
   real(R8P),    intent(inout), device :: fluxes_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Fluxes.
   real(R8P)                           :: el(5,5), er(5,5)                          !< Left and right eigenvalues.
   real(R8P)                           :: fmpc(1:2,-2:2,1:5)                        !< Fluxes -+ decomposition in c. space.
   real(R8P)                           :: fpmr(1:2,1:5)                             !< Fluxes +- reconstructed.
   logical                             :: ror_recompute                             !< Flag to perform ROR.
   integer(I4P)                        :: r, v, vv, rv                              !< Counter.

   call compute_eigenvectors_kernel(si=si,sir=sir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,g=g,q_aux_gpu=q_aux_gpu,el=el,er=er)
   call decompose_fluxes_convective_kernel(si=si,sir=sir,el=el,weno_s=weno_s,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,g=g,q_aux_gpu=q_aux_gpu,&
                                           fmpc=fmpc)
   do v=1, nv
      call weno_reconstruct_upwind_kernel(S=weno_s,weno_a=weno_a_gpu,weno_p=weno_p_gpu,weno_d=weno_d_gpu,weno_zeps=weno_zeps,&
                                          V=fmpc(:,:,v),VR=fpmr(:,v))
   enddo
   if (ror_number>0) then
      ROR : do r=2, ror_number
         ror_recompute = .false.
         do vv=1, size(ror_ivar_gpu)
            rv = ror_ivar_gpu(vv)
            if ((abs(fpmr(1,rv)-fmpc(1,0,rv))>ror_threshold*abs(fmpc(1,0,rv))) .or. &
                (abs(fpmr(2,rv)-fmpc(2,1,rv))>ror_threshold*abs(fmpc(2,1,rv)))) ror_recompute = .true.
         enddo
         if (ror_recompute) then
            ror_stats_gpu(b,i,j,k,dir) = ror_schemes_gpu(r)
            do v=1, nv
               call weno_reconstruct_upwind_kernel(S=ror_schemes_gpu(r),weno_a=weno_a_gpu,weno_p=weno_p_gpu,weno_d=weno_d_gpu,&
                                                   weno_zeps=weno_zeps,V=fmpc(:,:,v),VR=fpmr(:,v))
            enddo
         else
            exit ROR
         endif
      enddo ROR
   endif
   ! back projection in conservative variables space
   do v=1, nv
      fluxes_gpu(b,i,j,k,v) = 0._R8P
      do vv=1,nv
         fluxes_gpu(b,i,j,k,v) = fluxes_gpu(b,i,j,k,v) + er(vv,v) * (fpmr(1,vv) + fpmr(2,vv))
      enddo
   enddo
   endsubroutine compute_fluxes_conv_interface_kernel

   attributes(device) subroutine decompose_fluxes_convective_kernel(si,sir,el,weno_s,b,i,j,k,ngc,nv,g,q_aux_gpu,fmpc)
   !< Decompose convective fluxes.
   !< Flux vector splitting by local-Lax-Friedrics (Rusanov) with projection in pseudo-characteristics psace.
   integer(I4P), intent(in)         :: si(3)                                 !< Stencil increment.
   real(R8P)   , intent(in)         :: sir(3)                                !< Stencil increment, real cast.
   real(R8P)   , intent(in)         :: el(1:,1:)                             !< Left eigeinvectors.
   integer(I4P), intent(in), value  :: weno_s                                !< Weno stencils number/dimension.
   integer(I4P), intent(in), value  :: b, i, j, k                            !< Counter.
   integer(I4P), intent(in), value  :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in), value  :: nv                                    !< Number of conservative varibales.
   real(R8P),    intent(in), value  :: g                                     !< Specific heats ratio.
   real(R8P),    intent(in), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout)      :: fmpc(1:,1-weno_s:,1:)                 !< Fluxes -+ decomposition in characteristics space.
   real(R8P)                        :: fmp(2)                                !< Fluxes -+ decomposition in each cell stencils.
   real(R8P)                        :: evmax(5)                              !< Signals speeds.
   real(R8P)                        :: q(5), f(5)                            !< Conservative variables and fluxes.
   real(R8P)                        :: gc, wc                                !< Increments for fluxes decomposition.
   integer(I4P)                     :: v, vv, s, is, js, ks                  !< Counter.

   call compute_max_eigenvalues_kernel(si=si,sir=sir,weno_s=weno_s,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,q_aux_gpu=q_aux_gpu,evmax=evmax)
   do s=1-weno_s, weno_s
      is = i + (s) * si(1) ; js = j + (s) * si(2) ; ks = k + (s) * si(3)
      q(1) =      q_aux_gpu(b,is,js,ks,1)
      q(2) = q(1)*q_aux_gpu(b,is,js,ks,2)
      q(3) = q(1)*q_aux_gpu(b,is,js,ks,3)
      q(4) = q(1)*q_aux_gpu(b,is,js,ks,4)
      q(5) = q(1)*q_aux_gpu(b,is,js,ks,8) - q_aux_gpu(b,is,js,ks,7)
      f(1) = q_aux_gpu(b,is,js,ks,1)*q_aux_gpu(b,is,js,ks,2)*sir(1) + &
             q_aux_gpu(b,is,js,ks,1)*q_aux_gpu(b,is,js,ks,3)*sir(2) + &
             q_aux_gpu(b,is,js,ks,1)*q_aux_gpu(b,is,js,ks,4)*sir(3)
      f(2) = f(1)*q_aux_gpu(b,is,js,ks,2) + q_aux_gpu(b,is,js,ks,7)*sir(1)
      f(3) = f(1)*q_aux_gpu(b,is,js,ks,3) + q_aux_gpu(b,is,js,ks,7)*sir(2)
      f(4) = f(1)*q_aux_gpu(b,is,js,ks,4) + q_aux_gpu(b,is,js,ks,7)*sir(3)
      f(5) = f(1)*q_aux_gpu(b,is,js,ks,8)
      do v=1, nv
         wc = 0._R8P
         gc = 0._R8P
         do vv=1, nv
            wc = wc + el(vv,v) * q(vv)
            gc = gc + el(vv,v) * f(vv)
         enddo
         fmp(2) = 0.5_R8P * (gc + evmax(v) * wc)
         fmp(1) = gc - fmp(2)
         if (s<weno_s)   fmpc(2,s  ,v) = fmp(2)
         if (s>1-weno_s) fmpc(1,s-1,v) = fmp(1)
      enddo
   enddo
   endsubroutine decompose_fluxes_convective_kernel

   attributes(device) subroutine compute_eigenvectors_kernel(si,sir,b,i,j,k,ngc,nv,g,q_aux_gpu,el,er)
   ! Compute eigenvectors centered in inteface i,j,k/ip,jp,kp.
   integer(I4P), intent(in)         :: si(3)                                 !< Stencil increment.
   real(R8P)   , intent(in)         :: sir(3)                                !< Stencil increment, real cast.
   integer(I4P), intent(in), value  :: b, i, j, k                            !< Counter.
   integer(I4P), intent(in), value  :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in), value  :: nv                                    !< Number of conservative varibales.
   real(R8P),    intent(in), value  :: g                                     !< Specific heats ratio.
   real(R8P),    intent(in), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout)      :: el(1:5,1:5), er(1:5,1:5)              !< Left and right eigenvectors.
   real(R8P)                           :: uu, vv, ww, h, qq, c, ci, b1, b2   !< Roe average states.
   real(R8P)                           :: uvw, uvw_r1, uvw_r2                !< Velocity rotation accordingly dir.

   call compute_roe_average_kernel(q_aux_gpu=q_aux_gpu, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i+si(1), jp=j+si(2), kp=k+si(3), &
                                   uu=uu, vv=vv, ww=ww, h=h, qq=qq, c=c, ci=ci, b1=b1, b2=b2)

   uvw    =  uu*sir(1)+vv*sir(2)+ww*sir(3)
   uvw_r1 =  uu*sir(3)+vv*sir(1)+ww*sir(2)
   uvw_r2 = -uu*sir(2)+vv*sir(3)+ww*sir(1)

   er(1,1)=1._R8P ; er(1,2)=uu-c*sir(1)   ; er(1,3)=vv-c*sir(2) ; er(1,4)=ww-c*sir(3)   ; er(1,5)=h-uvw*c
   er(2,1)=1._R8P ; er(2,2)=uu            ; er(2,3)=vv          ; er(2,4)=ww            ; er(2,5)=qq
   er(3,1)=1._R8P ; er(3,2)=uu+c*sir(1)   ; er(3,3)=vv+c*sir(2) ; er(3,4)=ww+c*sir(3)   ; er(3,5)=h+uvw*c
   er(4,1)=0._R8P ; er(4,2)=sir(2)+sir(3) ; er(4,3)=sir(1)      ; er(4,4)=0._R8P        ; er(4,5)=uvw_r1
   er(5,1)=0._R8P ; er(5,2)=0._R8P        ; er(5,3)=sir(3)      ; er(5,4)=sir(1)+sir(2) ; er(5,5)=uvw_r2

   el(1,1)= 0.5_R8P*(b1+uvw*ci)      ;el(1,2)= 1._R8P-b1;el(1,3)= 0.5_R8P*(b1-uvw*ci)      ;el(1,4)=-uvw_r1;el(1,5)=-uvw_r2
   el(2,1)=-0.5_R8P*(b2*uu+ci*sir(1));el(2,2)= b2*uu    ;el(2,3)=-0.5_R8P*(b2*uu-ci*sir(1));el(2,4)= sir(3);el(2,5)=-sir(2)
   el(3,1)=-0.5_R8P*(b2*vv+ci*sir(2));el(3,2)= b2*vv    ;el(3,3)=-0.5_R8P*(b2*vv-ci*sir(2));el(3,4)= sir(1);el(3,5)= sir(3)
   el(4,1)=-0.5_R8P*(b2*ww+ci*sir(3));el(4,2)= b2*ww    ;el(4,3)=-0.5_R8P*(b2*ww-ci*sir(3));el(4,4)= sir(2);el(4,5)= sir(1)
   el(5,1)= 0.5_R8P*b2               ;el(5,2)=-b2       ;el(5,3)= 0.5_R8P*b2               ;el(5,4)= 0._R8P;el(5,5)= 0._R8P
   endsubroutine compute_eigenvectors_kernel

   attributes(device) subroutine compute_max_eigenvalues_kernel(si,sir,weno_s,b,i,j,k,ngc,nv,q_aux_gpu,evmax)
   ! Compute maximum eigenvalues in the big stencil.
   integer(I4P), intent(in)         :: si(3)                                 !< Stencil increment.
   real(R8P)   , intent(in)         :: sir(3)                                !< Stencil increment, real cast.
   integer(I4P), intent(in), value  :: weno_s                                !< Weno stencils number/dimension.
   integer(I4P), intent(in), value  :: b, i, j, k                            !< Counter.
   integer(I4P), intent(in), value  :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in), value  :: nv                                    !< Number of conservative varibales.
   real(R8P),    intent(in), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout)      :: evmax(1:)                             !< Maximum eigenvalues in the big stencil.
   real(R8P)                        :: uu, c                                 !< Speeds.
   real(R8P)                        :: ev(5)                                 !< Signals speeds.
   integer(I4P)                     :: s, is, js, ks, v                      !< Counter.

   evmax = -1._R8P
   do s=1, 2*weno_s
      is = i + (s-weno_s) * si(1) ; js = j + (s-weno_s) * si(2) ; ks = k + (s-weno_s) * si(3)
      uu = q_aux_gpu(b,is,js,ks,1+1*si(1)+2*si(2)+3*si(3))
      c  = q_aux_gpu(b,is,js,ks,9                        )
      ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2)
      do v=1,nv
         evmax(v) = max(ev(v),evmax(v))
      enddo
   enddo
   endsubroutine compute_max_eigenvalues_kernel

   attributes(device) subroutine compute_roe_average_kernel(ngc, b, i, j, k, ip, jp, kp, g, q_aux_gpu, &
                                                            uu, vv, ww, h, qq, c, ci, b1, b2)
   !< Compute Roe averaged quantities.
   integer(I4P), intent(in)         :: ngc                                       !< Number of ghost cells.
   integer(I4P), intent(in)         :: b, i, j, k, ip, jp, kp                    !< Left/right cells indexes.
   real(R8P),    intent(in)         :: g                                         !< Specific heats ratio.
   real(R8P),    intent(in), device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:) !< Auxiliary variables.
   real(R8P),    intent(out)        :: uu, vv, ww, h, qq, c, ci, b1, b2          !< Roe state average variables.
   real(R8P)                        :: ri, up, vp, wp, hp, r, rp1, cc            !< Local varbiables.

   ! left state (node i)
   ri        =  1._R8P/q_aux_gpu(b,i,j,k,1)
   uu        =  q_aux_gpu(b,i,j,k,2)
   vv        =  q_aux_gpu(b,i,j,k,3)
   ww        =  q_aux_gpu(b,i,j,k,4)
   h         =  q_aux_gpu(b,i,j,k,8)
   ! right state (node i+1)
   up        =  q_aux_gpu(b,ip,jp,kp,2)
   vp        =  q_aux_gpu(b,ip,jp,kp,3)
   wp        =  q_aux_gpu(b,ip,jp,kp,4)
   hp        =  q_aux_gpu(b,ip,jp,kp,8)
   ! Roe average state
   r         =  sqrt(q_aux_gpu(b,ip,jp,kp,1)*ri)
   rp1       =  1._R8P/(r+1._R8P)
   uu        =  (r*up+uu)*rp1
   vv        =  (r*vp+vv)*rp1
   ww        =  (r*wp+ww)*rp1
   h         =  (r*hp+h)*rp1
   qq        =  0.5_R8P * (uu*uu+vv*vv+ww*ww)
   cc        =  (g-1._R8P) * (h - qq)
   c         =  sqrt(cc)
   ci        =  1._R8P/c
   b2        = (g-1)/cc  ! alias 1/(cp*theta)
   b1        = b2 * qq   ! alias q/(cp*theta)
   endsubroutine compute_roe_average_kernel
endmodule adam_nasto_nvf_kernels
