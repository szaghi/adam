!< ADAM, PRISM FNL application kernels.

#include "fundal.H"

module adam_prism_fnl_kernels
!< ADAM, PRISM FNL application kernels.

! PRSIM modules
use :: adam_prism_common_library
! ADAM modules
use :: adam_fnl_weno_kernels
! third party modules
use :: fundal
use :: penf, only : I4P, I8P, R8P

implicit none
private
public :: compute_coils_current_dev
public :: compute_div_d_b_dev
public :: compute_fluxes_convective_dev
public :: compute_fluxes_difference_dev
public :: set_bc_q_gpu_dev
public :: set_sir_dev

contains
   ! public procedures
   subroutine compute_coils_current_dev(ni,nj,nk,ngc,blocks_number,time,a_gpu,d_gpu,f_gpu,phase_gpu,coil_flag_gpu,&
                                        td,j_vec_gpu,dx_gpu,q_gpu)
   !< Compute current coils sources.
   integer(I4P), intent(in)    :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)    :: blocks_number                          !< Number of blocks.
   integer(I4P), intent(in)    :: coil_flag_gpu(1:,1-ngc:,1-ngc:,1-ngc:) !< Coils ID map.
   real(R8P),    intent(in)    :: time                                   !< Simulation time, to compute current value if AC.
   real(R8P),    intent(in)    :: a_gpu(1:)                              !< Current amplitude (A).
   real(R8P),    intent(in)    :: f_gpu(1:)                              !< Current frequency, if AC (Hz).
   real(R8P),    intent(in)    :: phase_gpu(1:)                          !< Current initial phase, if AC.
   real(R8P),    intent(in)    :: d_gpu(1:)                              !< Wire diameter.
   real(R8P),    intent(in)    :: td                                     !< Coils transitory delay.
   real(R8P),    intent(in)    :: j_vec_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Current J versors into coils.
   real(R8P),    intent(in)    :: dx_gpu                                 !< Space step in x direction.
   real(R8P),    intent(inout) :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Field variables.
   real(R8P)                   :: current_density                        !< Current density.
   real(R8P)                   :: g                                      !< Starting polynomial transitory of coils.
   real(R8P)                   :: w_, g_, f_                             !< Step function coeff to avoid if in parallel regions.
   integer(I4P)                :: coil_id                                !< Uniq coild ID.
   integer(I4P)                :: i,j,k,b                                !< Counter.

   g = 10._R8P*(time/td)**3 - 15._R8P*(time/td)**4 + 6._R8P*(time/td)**5
   !$acc parallel loop independent gang vector collapse(4)                           &
   !$acc DEVICEVAR(coil_flag_gpu,A_gpu,f_gpu,phase_gpu,d_gpu,j_vec_gpu,q_gpu,dx_gpu) &
   !$acc firstprivate(g,time,td)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      coil_id = coil_flag_gpu(b,i,j,k)
      ! use step function to avoid the following original if
      ! if (time < td) then
      !    current_density = g*A_gpu(coil_id)/((d_gpu(coil_id)-dx_gpu)**2)*cos(phase_gpu(coil_id)*PI/180.0_R8P)
      ! else
      !    current_density =   A_gpu(coil_id)/((d_gpu(coil_id)-dx_gpu)**2)*cos(2*PI*f_gpu(coil_id)*(time-td) + &
      !                      phase_gpu(coil_id)*PI/180.0_R8P)
      ! endif
      w_ = 0.5_R8P * (sign(td-time,1._R8P) + 1._R8P) ! = 1 if td>time,                                  = 0 if td<time
      g_ = w_ * g + 1._R8P - w_                      ! = g if td>time,                                  = 1 if td<time
      f_ = w_ * 2._R8P*PI*f_gpu(coil_id)*(time-td)   ! = 2._R8P*PI*f_gpu(coil_id)*(time-td) if td>time, = 0 if td<time
      current_density = g_ * A_gpu(coil_id) / ((d_gpu(coil_id)-dx_gpu)**2) * cos(f_ + phase_gpu(coil_id)*PI/180.0_R8P)
      ! the following if is not necessary because j_vec is zero everywhere except in coils
      ! if (coil_id /= 0_I4P) then
      q_gpu(b,i,j,k,VAR_JX) = current_density * j_vec_gpu(b,i,j,k,1)
      q_gpu(b,i,j,k,VAR_JY) = current_density * j_vec_gpu(b,i,j,k,2)
      q_gpu(b,i,j,k,VAR_JZ) = current_density * j_vec_gpu(b,i,j,k,3)
      ! do n=7, 9
      !    q_gpu(b,i,j,k,n) = current_density*J_vec_gpu(b,i,j,k,n-6)
      ! enddo
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_coils_current_dev

   subroutine compute_div_d_b_dev(ni, nj, nk, ngc, blocks_number, dx_gpu, q_gpu, div_gpu)
   !< Compute div(D), div(B).
   integer(I4P), intent(in)    :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in)    :: dx_gpu                              !< Space step in x direction (m)
   real(R8P),    intent(in)    ::   q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Field variables.
   real(R8P),    intent(inout) :: div_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Divergence of D, B.
   integer(I4P)                :: i,j,k,b                             !< Counter

   !$acc parallel loop independent gang vector collapse(4) DEVICEVAR(div_gpu,q_gpu,dx_gpu)
   !$omp OMPLOOP DEVICEVAR(div_gpu,q_gpu,dx_gpu)
   do b=1, blocks_number
   do k=1, nk
   do j=1, nj
   do i=1, ni
      div_gpu(b,i,j,k,1) = 0.5_R8P*((q_gpu(b,i+1,j,k,1) - q_gpu(b,i-1,j,k,1))/dx_gpu + &
                                    (q_gpu(b,i,j+1,k,2) - q_gpu(b,i,j-1,k,2))/dx_gpu + &
                                    (q_gpu(b,i,j,k+1,3) - q_gpu(b,i,j,k-1,3))/dx_gpu)
      div_gpu(b,i,j,k,2) = 0.5_R8P*((q_gpu(b,i+1,j,k,4) - q_gpu(b,i-1,j,k,4))/dx_gpu + &
                                    (q_gpu(b,i,j+1,k,5) - q_gpu(b,i,j-1,k,5))/dx_gpu + &
                                    (q_gpu(b,i,j,k+1,6) - q_gpu(b,i,j,k-1,6))/dx_gpu)

   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_div_d_b_dev

   subroutine compute_fluxes_convective_dev(dir,blocks_number,ni,nj,nk,ngc,nv_c,               &
                                            weno_s,weno_zeps,weno_a_gpu,weno_p_gpu,weno_d_gpu, &
                                            evmax,ER_GPU,EL_GPU,si_gpu,sir_gpu,q_gpu,fluxes_gpu)
   !< Compute convective fluxes along x direction.
   !< @NOTE Be carefull with `si` and `sir` variables: are not present on GPU, probably a mapping
   !< is done by compiler, check it!
   integer(I4P), intent(in)    :: dir                                    !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: blocks_number                          !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                     !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                     !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                     !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                    !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                                   !< Number of conservative varibales in q vector.
   integer(I4P), intent(in)    :: weno_s                                 !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: weno_zeps                              !< Parameter to avoid division by zero.
   real(R8P),    intent(in)    :: weno_a_gpu(1:,0:,1:)                   !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p_gpu(1:,0:,0:,1:)                !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d_gpu(0:,0:,0:,1:)                !< Smoothness indicators coefficients.
   real(R8P),    intent(in)    :: evmax                                  !< Maximum eigenvalue.
   real(R8P),    intent(in)    :: ER_GPU(:,:,:), EL_GPU(:,:,:)           !< Eigenvectors.
   integer(I4P), intent(in)    :: si_gpu(3)                              !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)    :: sir_gpu(3)                             !< Directional (1=x,2=y,3=z) increment, real cast.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Fields variables.
   real(R8P),    intent(inout) :: fluxes_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Fluxes.
   integer(I4P)                :: b, i, j, k                             !< Counter.

   select case(dir)
   case(1)
      !$acc parallel loop independent gang vector collapse(3)                                         &
      !$acc default(none)                                                                             &
      !$acc DEVICEVAR(si_gpu,sir_gpu,q_gpu,fluxes_gpu,ER_GPU,EL_GPU,weno_a_gpu,weno_p_gpu,weno_d_gpu) &
      !$acc firstprivate(dir,blocks_number,ni,nj,nk,ngc,nv_c,weno_s,weno_zeps,evmax)                  &
      !$acc private(b,i,j,k)
      do b=1, blocks_number
      do j=1, nj
      do k=1, nk
      do i=0, ni
         call compute_fluxes_convective_ri_dev(dir=dir,b=b,i=i,j=j,k=k,ngc=ngc,nv_c=nv_c,              &
                                               weno_s=weno_s, weno_zeps=weno_zeps,                     &
                                               weno_a=weno_a_gpu, weno_p=weno_p_gpu, weno_d=weno_d_gpu,&
                                               evmax=evmax,ER=ER_GPU,EL=EL_GPU,                        &
                                               si=si_gpu,sir=sir_gpu,q=q_gpu,fluxes=fluxes_gpu)
      enddo
      enddo
      enddo
      enddo
   case(2)
      !$acc parallel loop independent gang vector collapse(3)                                         &
      !$acc default(none)                                                                             &
      !$acc DEVICEVAR(si_gpu,sir_gpu,q_gpu,fluxes_gpu,ER_GPU,EL_GPU,weno_a_gpu,weno_p_gpu,weno_d_gpu) &
      !$acc firstprivate(dir,blocks_number,ni,nj,nk,ngc,nv_c,weno_s,weno_zeps,evmax)                  &
      !$acc private(b,i,j,k)
      do b=1, blocks_number
      do i=1, ni
      do k=1, nk
      do j=0, nj
         call compute_fluxes_convective_ri_dev(dir=dir,b=b,i=i,j=j,k=k,ngc=ngc,nv_c=nv_c,              &
                                               weno_s=weno_s, weno_zeps=weno_zeps,                     &
                                               weno_a=weno_a_gpu, weno_p=weno_p_gpu, weno_d=weno_d_gpu,&
                                               evmax=evmax,ER=ER_GPU,EL=EL_GPU,                        &
                                               si=si_gpu,sir=sir_gpu,q=q_gpu,fluxes=fluxes_gpu)
      enddo
      enddo
      enddo
      enddo
   case(3)
      !$acc parallel loop independent gang vector collapse(3)                                         &
      !$acc default(none)                                                                             &
      !$acc DEVICEVAR(si_gpu,sir_gpu,q_gpu,fluxes_gpu,ER_GPU,EL_GPU,weno_a_gpu,weno_p_gpu,weno_d_gpu) &
      !$acc firstprivate(dir,blocks_number,ni,nj,nk,ngc,nv_c,weno_s,weno_zeps,evmax)                  &
      !$acc private(b,i,j,k)
      do b=1, blocks_number
      do i=1, ni
      do j=1, nj
      do k=0, nk
         call compute_fluxes_convective_ri_dev(dir=dir,b=b,i=i,j=j,k=k,ngc=ngc,nv_c=nv_c,              &
                                               weno_s=weno_s, weno_zeps=weno_zeps,                     &
                                               weno_a=weno_a_gpu, weno_p=weno_p_gpu, weno_d=weno_d_gpu,&
                                               evmax=evmax,ER=ER_GPU,EL=EL_GPU,                        &
                                               si=si_gpu,sir=sir_gpu,q=q_gpu,fluxes=fluxes_gpu)
      enddo
      enddo
      enddo
      enddo
   endselect
   endsubroutine compute_fluxes_convective_dev

   subroutine compute_fluxes_difference_dev(blocks_number, ni, nj, nk, ngc, nv_c,                     &
                                            dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu, q_gpu, &
                                            eta, chi, d_divergence_cleaner, b_divergence_cleaner, dq_gpu)
   !< Compute fluxes difference.
   integer(I4P), intent(in)    :: blocks_number                       !< Number of blocks.
   integer(I4P), intent(in)    :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)    :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)    :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                                !< Number of conservative varibales in q.
   real(R8P),    intent(in)    :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)  !< Space steps.
   real(R8P),    intent(in)    :: flx_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< X direction fluxes.
   real(R8P),    intent(in)    :: fly_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Y direction fluxes.
   real(R8P),    intent(in)    :: flz_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Z direction fluxes.
   real(R8P),    intent(in)    :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Fields.
   real(R8P),    intent(in)    :: eta, chi                            !< Coefficiente modello correzione divergenza campi.
   logical,      intent(in)    :: d_divergence_cleaner                !< Flag to perform electric field divergence cleaning.
   logical,      intent(in)    :: b_divergence_cleaner                !< Flag to perform magnetic field divergence cleaning.
   real(R8P),    intent(inout) :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Fluxes differences.
   integer(I4P)                :: b, i, j, k, v                       !< Counter.

   !$acc parallel loop independent gang vector collapse(4) &
   !$acc DEVICEVAR(dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu, dq_gpu, q_gpu)
   do b=1,blocks_number
   do k=1,nk
   do j=1,nj
   do i=1,ni
      do v=1, nv_c
         dq_gpu(b,i,j,k,v) = - (flx_gpu(b,i,j,k,v)-flx_gpu(b,i-1,j,k,v))/dx_gpu(b) &
                             - (fly_gpu(b,i,j,k,v)-fly_gpu(b,i,j-1,k,v))/dy_gpu(b) &
                             - (flz_gpu(b,i,j,k,v)-flz_gpu(b,i,j,k-1,v))/dz_gpu(b)
      enddo
      ! J sources
      dq_gpu(b,i,j,k,VAR_DX) = dq_gpu(b,i,j,k,VAR_DX) - q_gpu(b,i,j,k,VAR_JX)
      dq_gpu(b,i,j,k,VAR_DY) = dq_gpu(b,i,j,k,VAR_DY) - q_gpu(b,i,j,k,VAR_JY)
      dq_gpu(b,i,j,k,VAR_DZ) = dq_gpu(b,i,j,k,VAR_DZ) - q_gpu(b,i,j,k,VAR_JZ)
      ! corrections
      ! if (d_divergence_cleaner .and. .not.b_divergence_cleaner .and. eta>0._R8P) then
      !    dq_gpu(b,i,j,k,10) = dq_gpu(b,i,j,k,10) - chi/eta*chi/eta*q_gpu(b,i,j,k,10)
      ! elseif (D_divergence_cleaner .and. B_divergence_cleaner .and. eta>0._R8P) then
      !    dq_gpu(b,i,j,k,10) = dq_gpu(b,i,j,k,10) - chi/eta*chi/eta*q_gpu(b,i,j,k,10)
      !    dq_gpu(b,i,j,k,11) = dq_gpu(b,i,j,k,11) - chi/eta*chi/eta*q_gpu(b,i,j,k,11)
      ! endif
   enddo
   enddo
   enddo
   enddo
   endsubroutine compute_fluxes_difference_dev

   subroutine set_sir_dev(si_gpu, sir_gpu)
   !< Set directional increment si and sir on device.
   integer(I4P), intent(inout) :: si_gpu(3,3)  !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(inout) :: sir_gpu(3,3) !< Directional (1=x,2=y,3=z) increment, real cast.

   !$acc serial DEVICEVAR(si_gpu, sir_gpu)
   si_gpu (1:3,1) = [1,     0,     0     ]
   sir_gpu(1:3,1) = [1._R8P,0._R8P,0._R8P]
   si_gpu (1:3,2) = [0,     1,     0     ]
   sir_gpu(1:3,2) = [0._R8P,1._R8P,0._R8P]
   si_gpu (1:3,3) = [0,     0     ,1     ]
   sir_gpu(1:3,3) = [0._R8P,0._R8P,1._R8P]
   !$acc end serial
   endsubroutine set_sir_dev

   ! private procedures
   subroutine compute_fluxes_convective_ri_dev(dir,b,i,j,k,ngc,nv_c,                  &
                                               weno_s,weno_zeps,weno_a,weno_p,weno_d, &
                                               evmax,ER,EL,si,sir,q,fluxes)
   !< Compute convective fluxes at right interface of b,i,j,k.
   integer(I4P), intent(in)    :: dir                                 !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: b, i, j, k                          !< Counter.
   integer(I4P), intent(in)    :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                                !< Number of conservative varibales in q vector.
   integer(I4P), intent(in)    :: weno_s                              !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: weno_zeps                           !< Parameter to avoid division by zero.
   real(R8P),    intent(in)    :: weno_a(1:,0:,1:)                    !< Optimal weights.
   real(R8P),    intent(in)    :: weno_p(1:,0:,0:,1:)                 !< Polinomials coefficients.
   real(R8P),    intent(in)    :: weno_d(0:,0:,0:,1:)                 !< Smoothness indicators coefficients.
   real(R8P),    intent(in)    :: evmax                               !< Maximum eigenvalue.
   real(R8P),    intent(in)    :: ER(:,:,:), EL(:,:,:)                !< Eigenvectors.
   integer(I4P), intent(in)    :: si(3)                               !< Stencil increment.
   real(R8P),    intent(in)    :: sir(3)                              !< Stencil increment, real cast.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:)       !< Fields variables.
   real(R8P),    intent(inout) :: fluxes(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Fluxes.
   real(R8P)                   :: fmpc(1:2,1-S_MAX:-1+S_MAX,1:NV_MAX) !< Fluxes -+ decomposition in c. space.
   real(R8P)                   :: fpmr(1:2,1:NV_MAX)                  !< Fluxes +- reconstructed.
   integer(I4P)                :: v, vv                               !< Counter.
   !$acc routine seq

   call decompose_fluxes_convective_dev(dir=dir, si=si, sir=sir,           &
                                        b=b,i=i,j=j,k=k,ngc=ngc,nv_c=nv_c, &
                                        weno_s=weno_s, evmax=evmax, EL=EL, &
                                        q=q, fmpc=fmpc)
   do v=1, nv_c
      call weno_reconstruct_upwind_dev(S=weno_s, weno_a=weno_a, weno_p=weno_p, weno_d=weno_d,&
                                       weno_zeps=weno_zeps, V=fmpc(:,:,v), VR=fpmr(:,v))
   enddo
   ! back projection in conservative variables space
   do v=1, nv_c
      fluxes(b,i,j,k,v) = 0._R8P
      do vv=1,nv_c
         fluxes(b,i,j,k,v) = fluxes(b,i,j,k,v) + ER(vv,v,dir) * (fpmr(1,vv) + fpmr(2,vv))
      enddo
   enddo
   endsubroutine compute_fluxes_convective_ri_dev

   subroutine decompose_fluxes_convective_dev(dir,si,sir,b,i,j,k,ngc,nv_c,weno_s,evmax,EL,q,fmpc)
   !< Decompose convective fluxes.
   integer(I4P), intent(in)    :: dir                           !< Direction, 1=X, 2=Y, 3=Z.
   integer(I4P), intent(in)    :: si(3)                         !< Stencil increment.
   real(R8P),    intent(in)    :: sir(3)                        !< Stencil increment, real cast.
   integer(I4P), intent(in)    :: b, i, j, k                    !< Counter.
   integer(I4P), intent(in)    :: ngc                           !< Ghost cells number.
   integer(I4P), intent(in)    :: nv_c                          !< Number of conservative varibales in q vector.
   integer(I4P), intent(in)    :: weno_s                        !< Weno stencils number/dimension.
   real(R8P),    intent(in)    :: evmax                         !< Maximum eigenvalue.
   real(R8P),    intent(in)    :: EL(:,:,:)                     !< Left eigenvectors.
   real(R8P),    intent(in)    :: q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: fmpc(1:,1-weno_s:,1:)         !< Fluxes -+ decomposition in characteristics space.
   real(R8P)                   :: fmp(2)                        !< Fluxes -+ decomposition in each cell stencils.
   real(R8P)                   :: gc, wc                        !< Increments for fluxes decomposition.
   integer(I4P)                :: v, vv, s, is, js, ks          !< Counter.
   real(R8P)                   :: f(NV_MAX)                     !< Conservative fluxes.
   !$acc routine seq

   do s=1-weno_s, weno_s
      is = i + (s) * si(1) ; js = j + (s) * si(2) ; ks = k + (s) * si(3)
      call compute_convective_fluxes_maxwell_dev(sir=sir,q=q(b,is,js,ks,:),f=f)
      do v=1, nv_c
         wc = 0._R8P
         gc = 0._R8P
         do vv=1, nv_c
            wc = wc + EL(vv,v,dir) * q(b,is,js,ks,vv)
            gc = gc + EL(vv,v,dir) * f(vv)
         enddo
         fmp(2) = 0.5_R8P * (gc + evmax * wc)
         fmp(1) = gc - fmp(2)
         if (s<weno_s)   fmpc(2,s  ,v) = fmp(2)
         if (s>1-weno_s) fmpc(1,s-1,v) = fmp(1)
      enddo
   enddo
   endsubroutine decompose_fluxes_convective_dev

   subroutine compute_convective_fluxes_maxwell_dev(sir, q, f)
   !< Compute convective fluxes.
   !<
   !<```
   !< Fx(1) =  0       Fy(1) = -Bz/muz  Fz(1) =  By/muy
   !< Fx(2) =  Bz/muz  Fy(2) =  0       Fz(2) = -Bx/mux
   !< Fx(3) = -By/muy  Fy(3) =  Bx/mux  Fz(3) =  0
   !< Fx(4) =  0       Fy(4) =  Dz/epsz Fz(4) = -Dy/epsy
   !< Fx(5) = -Dz/epsz Fy(5) =  0       Fz(5) =  Dx/epsx
   !< Fx(6) =  Dy/epsy Fy(6) = -Dx/epsx Fz(6) =  0
   !<```
   real(R8P), intent(in)    :: sir(3) !< Directional (1=x,2=y,3=z) increment.
   real(R8P), intent(in)    :: q(1:)  !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)  !< Conservative fluxes.
   !$acc routine seq

   f(VAR_DX) =  0.0_R8P        * sir(1) + -q(VAR_BZ)/MU0  * sir(2) +  q(VAR_BY)/MU0  * sir(3)
   f(VAR_DY) =  q(VAR_BZ)/MU0  * sir(1) +  0.0_R8P        * sir(2) + -q(VAR_BX)/MU0  * sir(3)
   f(VAR_DZ) = -q(VAR_BY)/MU0  * sir(1) +  q(VAR_BX)/MU0  * sir(2) +  0.0_R8P        * sir(3)
   f(VAR_BX) =  0.0_R8P        * sir(1) +  q(VAR_DZ)/EPS0 * sir(2) + -q(VAR_DY)/EPS0 * sir(3)
   f(VAR_BY) = -q(VAR_DZ)/EPS0 * sir(1) +  0.0_R8P        * sir(2) +  q(VAR_DX)/EPS0 * sir(3)
   f(VAR_BZ) =  q(VAR_DY)/EPS0 * sir(1) + -q(VAR_DX)/EPS0 * sir(2) +  0.0_R8P        * sir(3)
   ! if (sir(1)==1._R8P) then
   !    f(VAR_DX) =  0.0_R8P
   !    f(VAR_DY) =  q(VAR_BZ)/MU0
   !    f(VAR_DZ) = -q(VAR_BY)/MU0
   !    f(VAR_BX) =  0.0_R8P
   !    f(VAR_BY) = -q(VAR_DZ)/EPS0
   !    f(VAR_BZ) =  q(VAR_DY)/EPS0
   ! elseif(sir(2)==1._R8P) then
   !    f(VAR_DX) = -q(VAR_BZ)/MU0
   !    f(VAR_DY) =  0.0_R8P
   !    f(VAR_DZ) =  q(VAR_BX)/MU0
   !    f(VAR_BX) =  q(VAR_DZ)/EPS0
   !    f(VAR_BY) =  0.0_R8P
   !    f(VAR_BZ) = -q(VAR_DX)/EPS0
   ! elseif(sir(3)==1._R8P) then
   !    f(VAR_DX) =  q(VAR_BY)/MU0
   !    f(VAR_DY) = -q(VAR_BX)/MU0
   !    f(VAR_DZ) =  0.0_R8P
   !    f(VAR_BX) = -q(VAR_DY)/EPS0
   !    f(VAR_BY) =  q(VAR_DX)/EPS0
   !    f(VAR_BZ) =  0.0_R8P
   ! endif
   endsubroutine compute_convective_fluxes_maxwell_dev
endmodule adam_prism_fnl_kernels
