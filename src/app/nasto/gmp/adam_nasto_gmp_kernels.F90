!< ADAM, NASTO GMP application kernels.
module adam_nasto_gmp_kernels
!< ADAM, NASTO GMP application kernels.

use adam_weno_object
use adam_nasto_gmp_cns_kernels
use penf, only : I4P, I8P, R8P

implicit none
private
public :: compute_fluxes_convective_gmp
public :: compute_fluxes_difference_gmp
public :: compute_fluxes_diffusive_gmp
public :: compute_q_aux_gmp
public :: compute_umax_gmp
public :: set_bc_q_gpu_gmp

contains
   ! public procedures
   subroutine compute_fluxes_convective_gmp(dir,blocks_number,ni,nj,nk,ngc,nv,                &
                                            weno_s,weno_a_gpu,weno_p_gpu,weno_d_gpu,weno_zeps,&
                                            ror_number,ror_schemes_gpu,ror_threshold,         &
                                            ror_ivar_gpu,ror_stats_gpu,                       &
                                            g,q_aux_gpu,fluxes_gpu)
   !< Compute convective fluxes along x direction.
   integer(I4P), intent(in)    :: dir                                       !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: blocks_number                             !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                        !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                        !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                        !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                       !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                        !< Number of conservative varibales.
   integer(I4P), intent(in)    :: weno_s                                    !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: weno_a_gpu(1:,0:,1:)                      !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p_gpu(1:,0:,0:,1:)                   !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d_gpu(0:,0:,0:,1:)                   !< Smoothness indicators coefficients.
   real(R8P),    intent(in)    :: weno_zeps                                 !< Parameter to avoid division by zero.
   integer(I4P), intent(in)    :: ror_number                                !< Number of ROR iterations.
   integer(I4P), intent(in)    :: ror_schemes_gpu(1:)                       !< Scheme (S value) for each ROR step.
   real(R8P),    intent(in)    :: ror_threshold                             !< ROR threshold triggering.
   integer(I4P), intent(in)    :: ror_ivar_gpu(1:)                          !< Index variables to check in ROR.
   integer(I4P), intent(inout) :: ror_stats_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Scheme (S value) for each ROR step.
   real(R8P),    intent(in)    :: g                                         !< Specific heats ratio.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Auxiliary variables.
   real(R8P),    intent(inout) :: fluxes_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Fluxes.
   integer(I4P)                :: si(3)                                     !< Stencil increment.
   real(R8P)                   :: sir(3)                                    !< Stencil increment, real cast.
   integer(I4P)                :: b, i, j, k                                !< Counter.

   select case(dir)
   case(1)
      si = [1,0,0]
      sir = real(si,R8P)
   case(2)
      si = [0,1,0]
      sir = real(si,R8P)
   case(3)
      si = [0,0,1]
      sir = real(si,R8P)
   endselect

   !$omp target data map(to:si,sir)
   !$omp target teams distribute parallel do collapse(4) has_device_addr(weno_a_gpu,weno_p_gpu,weno_d_gpu,          &
   !$omp&                                                                ror_schemes_gpu,ror_ivar_gpu,ror_stats_gpu,&
   !$omp&                                                                q_aux_gpu,fluxes_gpu)
   do b=1, blocks_number
   do k=1-si(3), nk
   do j=1-si(2), nj
   do i=1-si(1), ni
      call compute_fluxes_convective_device(dir=dir,si=si,sir=sir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,       &
                                            weno_s=weno_s,weno_a_gpu=weno_a_gpu,weno_p_gpu=weno_p_gpu, &
                                            weno_d_gpu=weno_d_gpu,weno_zeps=weno_zeps,                 &
                                            ror_number=ror_number, ror_schemes_gpu=ror_schemes_gpu,    &
                                            ror_threshold=ror_threshold, ror_ivar_gpu=ror_ivar_gpu,    &
                                            ror_stats_gpu=ror_stats_gpu,                               &
                                            g=g,q_aux_gpu=q_aux_gpu,fluxes_gpu=fluxes_gpu)
   enddo
   enddo
   enddo
   enddo
   !$omp end target data
   endsubroutine compute_fluxes_convective_gmp

   subroutine compute_fluxes_difference_gmp(blocks_number, ni, nj, nk, ngc, nv, ib_eps, &
                                            dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, dq_gpu)
   !< Compute fluxes difference.
   integer(I4P), intent(in)           :: blocks_number                       !< Number of blocks.
   integer(I4P), intent(in)           :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)           :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)           :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)           :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)           :: nv                                  !< Number of conservative varibales.
   real(R8P),    intent(in)           :: ib_eps                              !< Tolerance IB delta ratio.
   real(R8P),    intent(in)           :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)  !< Space steps.
   real(R8P),    intent(in)           :: flx_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< X direction fluxes.
   real(R8P),    intent(in)           :: fly_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Y direction fluxes.
   real(R8P),    intent(in)           :: flz_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Z direction fluxes.
   real(R8P),    intent(in), optional :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< IB distance function.
   real(R8P),    intent(inout)        ::  dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes differences.
   real(R8P)                          :: delta_x, delta_y, delta_z           !< Space steps.
   real(R8P)                          :: dx_locale, dy_locale, dz_locale     !< Local space steps.
   integer(I4P)                       :: b, i, j, k, v                       !< Counter.
   integer(I4P)                       :: all_solids                          !< Last phi index, all solids summary.

   if (present(phi_gpu)) then
      all_solids = ubound(phi_gpu, dim=5)
      !$omp target teams distribute parallel do collapse(4) has_device_addr(dx_gpu,dy_gpu,dz_gpu,&
      !$omp&                                                                flx_gpu,fly_gpu,flz_gpu,phi_gpu,dq_gpu)
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
            dq_gpu(b,i,j,k,v) = - (flx_gpu(b,i,j,k,v)-flx_gpu(b,i-1,j,k,v))/dx_locale &
                                - (fly_gpu(b,i,j,k,v)-fly_gpu(b,i,j-1,k,v))/dy_locale &
                                - (flz_gpu(b,i,j,k,v)-flz_gpu(b,i,j,k-1,v))/dz_locale
         enddo
      enddo
      enddo
      enddo
      enddo
   else
      !$omp target teams distribute parallel do collapse(4) has_device_addr(dx_gpu,dy_gpu,dz_gpu,&
      !$omp&                                                                flx_gpu,fly_gpu,flz_gpu,dq_gpu)
      do k=1,nk
      do j=1,nj
      do i=1,ni
      do b=1,blocks_number
         do v=1, nv
            dq_gpu(b,i,j,k,v) = - (flx_gpu(b,i,j,k,v)-flx_gpu(b,i-1,j,k,v))/dx_gpu(b) &
                                - (fly_gpu(b,i,j,k,v)-fly_gpu(b,i,j-1,k,v))/dy_gpu(b) &
                                - (flz_gpu(b,i,j,k,v)-flz_gpu(b,i,j,k-1,v))/dz_gpu(b)
         enddo
      enddo
      enddo
      enddo
      enddo
   endif
   endsubroutine compute_fluxes_difference_gmp

   subroutine compute_fluxes_diffusive_gmp(blocks_number, ni, nj, nk, ngc, nv, mu, kd, &
                                           q_aux_gpu, dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu)
   !< Compute diffusive fluxes.
   integer(I4P), intent(in)    :: blocks_number                         !< Blocks number.
   integer(I4P), intent(in)    :: ni, nj, nk                            !< Grid dimensionns.
   integer(I4P), intent(in)    :: ngc                                   !< Number of ghost cells.
   integer(I4P), intent(in)    :: nv                                    !< Number of conservative variables.
   real(R8P),    intent(in)    :: mu                                    !< Viscosity.
   real(R8P),    intent(in)    :: kd                                    !< Thermal diffusivity.
   real(R8P),    intent(in)    :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)    !< Space steps.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary varibales
   real(R8P),    intent(inout) ::   flx_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes along x.
   real(R8P),    intent(inout) ::   fly_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes along y.
   real(R8P),    intent(inout) ::   flz_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes along z.
   real(R8P)                   :: vel_u, vel_v, vel_w                   !< (Mean) velocity.
   real(R8P)                   :: du_dx, dv_dx, dw_dx                   !< Velocity derivative along x.
   real(R8P)                   :: du_dy, dv_dy, dw_dy                   !< Velocity derivative along y.
   real(R8P)                   :: du_dz, dv_dz, dw_dz                   !< Velocity derivative along z.
   real(R8P)                   :: sigq, sigl                            !< Sigmas.
   real(R8P)                   :: tau_1_1, tau_2_1, tau_3_1, dT_dx      !< Stress tensor, x elements.
   real(R8P)                   :: tau_1_2, tau_2_2, tau_3_2, dT_dy      !< Stress tensor, y elements.
   real(R8P)                   :: tau_1_3, tau_2_3, tau_3_3, dT_dz      !< Stress tensor, z elements.
   integer(I4P)                :: b, i, j, k, v                         !< Counter.

   !$omp target teams distribute parallel do collapse(4) has_device_addr(dx_gpu,dy_gpu,dz_gpu,&
   !$omp&                                                                flx_gpu,q_aux_gpu)
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

   !$omp target teams distribute parallel do collapse(4) has_device_addr(dx_gpu,dy_gpu,dz_gpu,&
   !$omp&                                                                fly_gpu,q_aux_gpu)
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

   !$omp target teams distribute parallel do collapse(4) has_device_addr(dx_gpu,dy_gpu,dz_gpu,&
   !$omp&                                                                flz_gpu,q_aux_gpu)
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
   endsubroutine compute_fluxes_diffusive_gmp

   subroutine compute_umax_gmp(ni, nj, nk, ngc, blocks_number, mu, dx_gpu, dy_gpu, dz_gpu, q_aux_gpu, umax)
   !< Compute maximum speed.
   integer(I4P), intent(in)  :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)  :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)  :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)  :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)  :: blocks_number                         !< Number of blocks.
   real(R8P),    intent(in)  :: mu                                    !< Dynamic viscosity.
   real(R8P),    intent(in)  :: dx_gpu(1:)                            !< X space step.
   real(R8P),    intent(in)  :: dy_gpu(1:)                            !< Y space step.
   real(R8P),    intent(in)  :: dz_gpu(1:)                            !< Z space step.
   real(R8P),    intent(in)  :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary varibales.
   real(R8P),    intent(out) :: umax                                  !< Maximum speed.
   real(R8P)                 :: ss                                    !< Speed of sound.
   integer(I4P)              :: i, j, k, b                            !< Counter.
   real(R8P)                 :: dx_locale, dy_locale, dz_locale       !< Local space steps.

   ! TODO: check for phi inside, umax could be wrong otherwise
   umax = 0._R8P
   !$omp target teams distribute parallel do collapse(4) has_device_addr(dx_gpu,dy_gpu,dz_gpu,q_aux_gpu) reduction(max:umax)
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
   endsubroutine compute_umax_gmp

   subroutine set_bc_q_gpu_gmp(BC_EXTRAPOLATION, BC_INFLOW, nv, ngc, cv, R, &
                               local_map_bc_gpu, fec_1_6_array_gpu, q_bc_vars_gpu, q_gpu)
   !< Set BC over q.
   integer(I4P), intent(in)    :: BC_EXTRAPOLATION        !< Extrapolation BC parameter.
   integer(I4P), intent(in)    :: BC_INFLOW               !< Inflow BC parameter.
   integer(I4P), intent(in)    :: nv                      !< Number of variables.
   integer(I4P), intent(in)    :: ngc                     !< Ghost cells number.
   real(R8P),    intent(in)    :: cv                      !< Constant volume specific heat.
   real(R8P),    intent(in)    :: R                       !< Gas constant.
   integer(I8P), intent(in)    :: local_map_bc_gpu(:,:,:) !< Local map for BC ghost cells.
   integer(I4P), intent(in)    :: fec_1_6_array_gpu(:)    !< Local map for BC ghost cells.
   real(R8P),    intent(in)    :: q_bc_vars_gpu(:,:)      !< Boundary variables.
   real(R8P),    intent(inout) :: q_gpu(1:,    &
                                        1-ngc:,&
                                        1-ngc:,&
                                        1-ngc:,1:)        !< Conservative variables.
   integer(I4P)                :: b                       !< Counter.
   integer(I4P)                :: c, i, j, k, v           !< Counter.
   integer(I4P)                :: idelta                  !< IJK delta step for extrapolation.
   integer(I4P)                :: jdelta                  !< IJK delta step for extrapolation.
   integer(I4P)                :: kdelta                  !< IJK delta step for extrapolation.
   integer(I4P)                :: bc_type                 !< Boundary condition type.
   integer(I4P)                :: crown                   !< Crown counter.
   integer(I4P)                :: fec                     !< Boundary fec (1 to 26).
   integer(I4P)                :: fec_1_6                 !< Boundary fec (1 to 6).

   do crown=1, ngc
      !$omp target teams distribute parallel do has_device_addr(fec_1_6_array_gpu,local_map_bc_gpu,q_bc_vars_gpu,q_gpu)
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
   enddo
   endsubroutine set_bc_q_gpu_gmp

   ! private procedures
   subroutine compute_fluxes_convective_device(dir,si,sir,b,i,j,k,ngc,nv,                        &
                                               weno_s,weno_a_gpu,weno_p_gpu,weno_d_gpu,weno_zeps,&
                                               ror_number,ror_schemes_gpu,ror_threshold,         &
                                               ror_ivar_gpu,ror_stats_gpu,                       &
                                               g,q_aux_gpu,fluxes_gpu)
   !< Compute convective fluxes at right interface of b,i,j,k.
   integer(I4P), intent(in)            :: dir                                       !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)            :: si(3)                                     !< Stencil increment.
   real(R8P),    intent(in)            :: sir(3)                                    !< Stencil increment, real cast.
   integer(I4P), intent(in)            :: b, i, j, k                                !< Counter.
   integer(I4P), intent(in)            :: ngc                                       !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                        !< Number of conservative varibales.
   integer(I4P), intent(in)            :: weno_s                                    !< Weno stencils number/dimension.
   real(R8P),    intent(in)            :: weno_a_gpu(1:,0:,1:)                      !< Optimal weights.
   real(R8P),    intent(in)            :: weno_p_gpu(1:,0:,0:,1:)                   !< Polinomials coefficients.
   real(R8P),    intent(in)            :: weno_d_gpu(0:,0:,0:,1:)                   !< Smoothness indicators coefficients.
   real(R8P),    intent(in)            :: weno_zeps                                 !< Parameter to avoid division by zero.
   integer(I4P), intent(in)            :: ror_number                                !< Number of ROR iterations.
   integer(I4P), intent(in)            :: ror_schemes_gpu(1:)                       !< Scheme (S value) for each ROR step.
   real(R8P),    intent(in)            :: ror_threshold                             !< ROR threshold triggering.
   integer(I4P), intent(in)            :: ror_ivar_gpu(1:)                          !< Index variables to check in ROR.
   integer(I4P), intent(inout)         :: ror_stats_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Scheme (S value) for each ROR step.
   real(R8P),    intent(in)            :: g                                         !< Specific heats ratio.
   real(R8P),    intent(in)            :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)     !< Auxiliary variables.
   real(R8P),    intent(inout)         :: fluxes_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)    !< Fluxes.
   real(R8P)                           :: el(nv,nv), er(nv,nv)                      !< Left and right eigenvalues.
   real(R8P)                           :: fmpc(1:2,1-weno_s-1+weno_s,1:nv)          !< Fluxes -+ decomposition in c. space.
   real(R8P)                           :: fpmr(1:2,1:nv)                            !< Fluxes +- reconstructed.
   logical                             :: ror_recompute                             !< Flag to perform ROR.
   integer(I4P)                        :: r, v, vv, rv                              !< Counter.

   !$omp declare target
   call compute_eigenvectors_device(si=si,sir=sir,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,g=g,q_aux_gpu=q_aux_gpu,el=el,er=er)
   call decompose_fluxes_convective_device(si=si,sir=sir,el=el,weno_s=weno_s,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,g=g,q_aux_gpu=q_aux_gpu,&
                                           fmpc=fmpc)
   do v=1, nv
      call weno_reconstruct_upwind(S=weno_s,weno_a=weno_a_gpu,weno_p=weno_p_gpu,weno_d=weno_d_gpu,weno_zeps=weno_zeps,&
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
               call weno_reconstruct_upwind(S=ror_schemes_gpu(r),weno_a=weno_a_gpu,weno_p=weno_p_gpu,weno_d=weno_d_gpu,&
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
      ! fluxes_gpu(b,i,j,k,v) = fluxes_gpu(b,i,j,k,v) + fpmr(1,v)
      ! fluxes_gpu(b,i,j,k,v) = fluxes_gpu(b,i,j,k,v) + fpmr(2,v)
   enddo
   endsubroutine compute_fluxes_convective_device

   subroutine decompose_fluxes_convective_device(si,sir,el,weno_s,b,i,j,k,ngc,nv,g,q_aux_gpu,fmpc)
   !< Decompose convective fluxes.
   !< Flux vector splitting by local-Lax-Friedrics (Rusanov) with projection in pseudo-characteristics psace.
   integer(I4P), intent(in)    :: si(3)                                 !< Stencil increment.
   real(R8P),    intent(in)    :: sir(3)                                !< Stencil increment, real cast.
   real(R8P),    intent(in)    :: el(1:,1:)                             !< Left eigeinvectors.
   integer(I4P), intent(in)    :: weno_s                                !< Weno stencils number/dimension.
   integer(I4P), intent(in)    :: b, i, j, k                            !< Counter.
   integer(I4P), intent(in)    :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                    !< Number of conservative varibales.
   real(R8P),    intent(in)    :: g                                     !< Specific heats ratio.
   real(R8P),    intent(in)    :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fmpc(1:,1-weno_s:,1:)                 !< Fluxes -+ decomposition in characteristics space.
   real(R8P)                   :: fmp(2)                                !< Fluxes -+ decomposition in each cell stencils.
   real(R8P)                   :: evmax(nv)                             !< Signals speeds.
   real(R8P)                   :: q(nv), f(nv)                          !< Conservative variables and fluxes.
   real(R8P)                   :: gc, wc                                !< Increments for fluxes decomposition.
   integer(I4P)                :: v, vv, s, is, js, ks                  !< Counter.

   !$omp declare target
   call compute_max_eigenvalues_device(si=si,sir=sir,weno_s=weno_s,b=b,i=i,j=j,k=k,ngc=ngc,nv=nv,q_aux_gpu=q_aux_gpu,evmax=evmax)
   do s=1-weno_s, weno_s
      is = i + (s) * si(1) ; js = j + (s) * si(2) ; ks = k + (s) * si(3)
      call compute_conservatives_device(b=b,i=is,j=js,k=ks,ngc=ngc,q_aux_gpu=q_aux_gpu,q=q)
      call compute_conservative_fluxes_device(sir=sir,b=b,i=is,j=js,k=ks,ngc=ngc,q_aux_gpu=q_aux_gpu,f=f)
      do v=1, nv
         wc = 0._R8P
         gc = 0._R8P
         do vv=1, nv
            wc = wc + el(vv,v) * q(vv)
            gc = gc + el(vv,v) * f(vv)
         enddo
         fmp(2) = 0.5_R8P * (gc + evmax(v) * wc)
         fmp(1) = gc - fmp(2)
         ! fmp(2) = 0.5_R8P * (f(v) + 1.2_R8P*evmax(v) * q(v))
         ! fmp(1) = f(v) - fmp(2)
         if (s<weno_s)   fmpc(2,s  ,v) = fmp(2)
         if (s>1-weno_s) fmpc(1,s-1,v) = fmp(1)
      enddo
   enddo
   endsubroutine decompose_fluxes_convective_device
endmodule adam_nasto_gmp_kernels
