!< ADAM, NASTO NVF application kernels.
module adam_nasto_nvf_kernels
!< ADAM, NASTO NVF application kernels.

use PENF, only : I4P, I8P, R8P
use CUDAFOR

implicit none
private
! AMR procedures
public :: move_phi_cuf
! numerical procedures
public :: compute_flux_conv_x_kernel
public :: compute_flux_conv_y_kernel
public :: compute_flux_conv_z_kernel
public :: compute_flux_conv_x_central_kernel
public :: compute_flux_conv_y_central_kernel
public :: compute_flux_conv_z_central_kernel
public :: compute_fluxes_difference_cuf
public :: compute_fluxes_diffusive_cuf
public :: compute_q_aux_cuf
public :: compute_q_gradient_cuf
public :: compute_umax_cuf
public :: set_bc_q_gpu_cuf
! eikonal procedures
public :: evolve_eikonal_q_gpu_cuf
public :: invert_eikonal_q_gpu_cuf
public :: compute_eikonal_dq_gpu
! RK procedures
public :: compute_rk_linear_gpu_cuf
public :: compute_rk_q_gpu_cuf
public :: compute_rk_prhs_gpu_cuf

contains
   ! AMR procedures
   subroutine move_phi_cuf(ni, nj, nk, ngc, blocks_number, velocity, phi_gpu, dphi_gpu)
   !< Move phi and the actual ptree representation.
   integer(I4P), intent(in)            :: ni                                   !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                   !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                   !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                  !< Ghost grid number.
   integer(I4P), intent(in)            :: blocks_number                        !< Number of blocks.
   real(R8P),    intent(in)            :: velocity(3)                          !< Velocity of the movement.
   real(R8P),    intent(inout), device ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(inout), device :: dphi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function gradient.
   real(R8P)                           :: n_phi_x, n_phi_y, n_phi_z, n_phi     !< Eikonal direction.
   integer(I4P)                        :: b, i, j, k, v                        !< Counter.
   integer(I4P)                        :: iercuda                              !< Error trapping flag for CUDAFortran.

   n_phi_x = velocity(1)
   n_phi_y = velocity(2)
   n_phi_z = velocity(3)
   n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
   n_phi = 0.9_R8P / n_phi
   n_phi_x = n_phi_x * n_phi
   n_phi_y = n_phi_y * n_phi
   n_phi_z = n_phi_z * n_phi

   !$cuf kernel do(4) <<<*,*>>>
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      do v=1, 1
         dphi_gpu(b,i,j,k,v) = 0._R8P
      enddo
      if (n_phi_x > 0._R8P) then
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_x) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i-1,j,k,v))
         enddo
      else
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_x) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i+1,j,k,v))
         enddo
      endif
      if (n_phi_y > 0._R8P) then
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_y) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j-1,k,v))
         enddo
      else
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_y) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j+1,k,v))
         enddo
      endif
      if (n_phi_z > 0._R8P) then
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_z) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j,k-1,v))
         enddo
      else
         do v=1, 1
            dphi_gpu(b,i,j,k,v) = dphi_gpu(b,i,j,k,v) + abs(n_phi_z) * (phi_gpu(b,i,j,k,v) - phi_gpu(b,i,j,k+1,v))
         enddo
      endif
   enddo
   enddo
   enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()

   !$cuf kernel do(4) <<<*,*>>>
   do k=1, nk
   do j=1, nj
   do i=1, ni
   do b=1, blocks_number
      do v=1, 1
         phi_gpu(b,i,j,k,v) = phi_gpu(b,i,j,k,v) - dphi_gpu(b,i,j,k,v)
      enddo
   enddo
   enddo
   enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine move_phi_cuf

   ! numerical procedures
   attributes(global) subroutine compute_flux_conv_x_kernel(blocks_number, ni, nj, nk, ngc, nv, iweno, dha, g, R, cv,       &
                                                            ror_threshold, enable_ror_stats, cell_scheme_gpu, ror_ivar_gpu, &
                                                            ror_schemes_gpu, q_aux_gpu, ror_stats_gpu, gplus, gminus, flx_gpu)
   !< Compute convective fluxes by means of upwind WENO reconstruction, x axis direction.
   integer,      intent(in),    value  :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P),    intent(in),    value  :: dha, g, R, cv
   real(R8P),    intent(in),    value  :: ror_threshold
   logical,      intent(in),    value  :: enable_ror_stats
   integer(I4P), intent(in),    device :: cell_scheme_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), intent(in),    device :: ror_ivar_gpu(1:)
   integer(I4P), intent(in),    device :: ror_schemes_gpu(1:)
   real(R8P),    intent(in),    device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer(I4P), intent(inout), device :: ror_stats_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout), device ::  gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout), device :: gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout), device :: flx_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer                             :: b, i, j, k, l, ll, m, mm, v
   real(R8P)                           :: er(5,5), el(5,5), ev(5), evmax(5), ghat(5), gl(5), gr(5), fi(5), vi(5)
   real(R8P)                           :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                           :: gc, wc
   integer                             :: wenorec_scheme, index_var
   logical                             :: ror_to_recompute

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if (b > blocks_number .or. j > nj) return
   do k=1,nk
      do i=0,ni ! loop on faces
         ! compute Roe average
         call compute_roe_average(q_aux_gpu=q_aux_gpu, dha=dha, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i+1, jp=j, kp=k, &
                                  uu=uu, vv=vv, ww=ww, h=h, ya=ya, qq=qq, c=c, ci=ci, b1=b1, b2=b2)
         ! compute right and left eigenvectors matrices (at Roe state)
         er(1,1)=1._R8P ; er(1,2)=uu-c   ; er(1,3)=vv     ; er(1,4)=ww     ; er(1,5)=h-uu*c
         er(2,1)=1._R8P ; er(2,2)=uu     ; er(2,3)=vv     ; er(2,4)=ww     ; er(2,5)=qq
         er(3,1)=1._R8P ; er(3,2)=uu+c   ; er(3,3)=vv     ; er(3,4)=ww     ; er(3,5)=h+uu*c
         er(4,1)=0._R8P ; er(4,2)=0._R8P ; er(4,3)=1._R8P ; er(4,4)=0._R8P ; er(4,5)=vv
         er(5,1)=0._R8P ; er(5,2)=0._R8P ; er(5,3)=0._R8P ; er(5,4)=1._R8P ; er(5,5)=ww

         el(1,1)= 0.5_R8P*(b1+uu*ci) ; el(1,2)=1._R8P-b1 ; el(1,3)= 0.5_R8P*(b1-uu*ci) ; el(1,4)=-vv    ; el(1,5)=-ww
         el(2,1)=-0.5_R8P*(b2*uu+ci) ; el(2,2)=b2*uu     ; el(2,3)=-0.5_R8P*(b2*uu-ci) ; el(2,4)=0._R8P ; el(2,5)=0._R8P
         el(3,1)=-0.5_R8P*(b2*vv   ) ; el(3,2)=b2*vv     ; el(3,3)=-0.5_R8P*(b2*vv   ) ; el(3,4)=1._R8P ; el(3,5)=0._R8P
         el(4,1)=-0.5_R8P*(b2*ww   ) ; el(4,2)=b2*ww     ; el(4,3)=-0.5_R8P*(b2*ww   ) ; el(4,4)=0._R8P ; el(4,5)=1._R8P
         el(5,1)= 0.5_R8P*b2         ; el(5,2)=-b2       ; el(5,3)= 0.5_R8P*b2         ; el(5,4)=0._R8P ; el(5,5)=0._R8P

         ! Find max eigenvalues on the stencil
         do m=1,nv  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*iweno ! LLF
            ll = i + l - iweno
            uu = q_aux_gpu(b,ll,j,k,2)
            c  = q_aux_gpu(b,ll,j,k,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2)
            do m=1,nv
               evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*iweno ! loop over the stencil centered at face i
            ll = i + l - iweno
            vi(1) = q_aux_gpu(b,ll,j,k,1)
            vi(2) = vi(1)*q_aux_gpu(b,ll,j,k,2)
            vi(3) = vi(1)*q_aux_gpu(b,ll,j,k,3)
            vi(4) = vi(1)*q_aux_gpu(b,ll,j,k,4)
            vi(5) = vi(1)*(cv*q_aux_gpu(b,ll,j,k,6)+                                                      &
                    0.5_R8P*(q_aux_gpu(b,ll,j,k,2)**2+q_aux_gpu(b,ll,j,k,3)**2+q_aux_gpu(b,ll,j,k,4)**2)+ &
                    q_aux_gpu(b,ll,j,k,5)*dha)
            fi(1) = vi(2)
            fi(2) = fi(1) * q_aux_gpu(b,ll,j,k,2) + q_aux_gpu(b,ll,j,k,7)
            fi(3) = fi(1) * q_aux_gpu(b,ll,j,k,3)
            fi(4) = fi(1) * q_aux_gpu(b,ll,j,k,4)
            fi(5) = fi(1) * vi(5) / vi(1) + q_aux_gpu(b,ll,j,k,7)*q_aux_gpu(b,ll,j,k,2)
            do m=1,nv
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,nv
                  wc = wc + el(mm,m) * vi(mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,j,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,j,k,b) = gc - gplus(m,l,j,k,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         wenorec_scheme = cell_scheme_gpu(b,i,j,k,1)
         call weno_reconstruction_kernel(nvar=nv, vp=gplus(1,1,j,k,b), vm=gminus(1,1,j,k,b), vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)

         ror_x: do m = 2, size(ror_schemes_gpu)
            ror_to_recompute = .false.
            do mm = 1,size(ror_ivar_gpu)
                index_var = ror_ivar_gpu(mm)
                if((abs(gl(index_var)-gplus(index_var,iweno,j,k,b))    > ror_threshold*abs(gplus(index_var,iweno,j,k,b))   ) .or. &
                   (abs(gr(index_var)-gminus(index_var,iweno+1,j,k,b)) > ror_threshold*abs(gminus(index_var,iweno+1,j,k,b))) ) then
                   ror_to_recompute = .true.
                endif
            enddo
            if(ror_to_recompute) then
               wenorec_scheme = ror_schemes_gpu(m)
               call weno_reconstruction_kernel(nvar=nv, vp=gplus(1,1,j,k,b), vm=gminus(1,1,j,k,b), vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)
            else
               exit ror_x
            endif
         enddo ror_x

         if (enable_ror_stats) ror_stats_gpu(b,i,j,k,1) = wenorec_scheme

         ! Reassemble + and - characteristic fluxes
         do m=1,nv
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,nv
            flx_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,nv
               flx_gpu(b,i,j,k,m) = flx_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo

      enddo
   enddo
   endsubroutine compute_flux_conv_x_kernel

   attributes(global) subroutine compute_flux_conv_y_kernel(blocks_number, ni, nj, nk, ngc, nv, iweno, dha, g, R, cv,       &
                                                            ror_threshold, enable_ror_stats, cell_scheme_gpu, ror_ivar_gpu, &
                                                            ror_schemes_gpu, q_aux_gpu, ror_stats_gpu, gplus, gminus, fly_gpu)
   !< Compute convective fluxes by means of upwind WENO reconstruction, y axis direction.
   integer,      intent(in),    value  :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P),    intent(in),    value  :: dha, g, R, cv
   real(R8P),    intent(in),    value  :: ror_threshold
   logical,      intent(in),    value  :: enable_ror_stats
   integer(I4P), intent(in),    device :: cell_scheme_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), intent(in),    device :: ror_ivar_gpu(1:)
   integer(I4P), intent(in),    device :: ror_schemes_gpu(1:)
   real(R8P),    intent(in),    device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer(I4P), intent(inout), device :: ror_stats_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout), device ::  gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout), device :: gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout), device :: fly_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer                             :: b, i, j, k, l, ll, m, mm, v
   real(R8P)                           :: er(5,5), el(5,5), ev(5), evmax(5), ghat(5), gl(5), gr(5), fi(5), vi(5)
   real(R8P)                           :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                           :: gc, wc
   integer                             :: wenorec_scheme, index_var
   logical                             :: ror_to_recompute

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if (b > blocks_number .or. i > ni) return
   do k=1,nk
      do j=0,nj ! loop on faces
         ! Compute Roe average
         call compute_roe_average(q_aux_gpu=q_aux_gpu, dha=dha, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i, jp=j+1, kp=k, &
                                  uu=uu, vv=vv, ww=ww, h=h, ya=ya, qq=qq, c=c, ci=ci, b1=b1, b2=b2)
         ! Compute right and left eigenvectors matrices (at Roe state)
         er(1,1)=1._R8P ; er(1,2)=uu     ; er(1,3)=vv-c   ; er(1,4)=ww     ; er(1,5)=h-vv*c
         er(2,1)=1._R8P ; er(2,2)=uu     ; er(2,3)=vv     ; er(2,4)=ww     ; er(2,5)=qq
         er(3,1)=1._R8P ; er(3,2)=uu     ; er(3,3)=vv+c   ; er(3,4)=ww     ; er(3,5)=h+vv*c
         er(4,1)=0._R8P ; er(4,2)=1._R8P ; er(4,3)=0._R8P ; er(4,4)=0._R8P ; er(4,5)=ww
         er(5,1)=0._R8P ; er(5,2)=0._R8P ; er(5,3)=0._R8P ; er(5,4)=1._R8P ; er(5,5)=-uu

         el(1,1)= 0.5_R8P*(b1+vv*ci) ; el(1,2)=1._R8P-b1 ; el(1,3)=0.5_R8P*(b1-vv*ci)  ; el(1,4)=-ww    ; el(1,5)=uu
         el(2,1)=-0.5_R8P*(b2*uu)    ; el(2,2)=b2*uu     ; el(2,3)=-0.5_R8P*(b2*uu)    ; el(2,4)=0._R8P ; el(2,5)=-1._R8P
         el(3,1)=-0.5_R8P*(b2*vv+ci) ; el(3,2)=b2*vv     ; el(3,3)=-0.5_R8P*(b2*vv-ci) ; el(3,4)=0._R8P ; el(3,5)=0._R8P
         el(4,1)=-0.5_R8P*(b2*ww)    ; el(4,2)=b2*ww     ; el(4,3)=-0.5_R8P*(b2*ww)    ; el(4,4)=1._R8P ; el(4,5)=0._R8P
         el(5,1)= 0.5_R8P*b2         ; el(5,2)=-b2       ; el(5,3)=0.5_R8P*b2          ; el(5,4)=0._R8P ; el(5,5)=0._R8P

         ! Find max eigenvalues on the stencil
         do m=1,nv  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*iweno ! LLF
            ll = j + l - iweno
            uu = q_aux_gpu(b,i,ll,k,3)
            c  = q_aux_gpu(b,i,ll,k,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2) ;
            do m=1,nv
               evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*iweno ! loop over the stencil centered at face i
            ll = j + l - iweno
            vi(1) = q_aux_gpu(b,i,ll,k,1)
            vi(2) = vi(1)*q_aux_gpu(b,i,ll,k,2)
            vi(3) = vi(1)*q_aux_gpu(b,i,ll,k,3)
            vi(4) = vi(1)*q_aux_gpu(b,i,ll,k,4)
            vi(5) = vi(1)*(cv*q_aux_gpu(b,i,ll,k,6)+                                                  &
                0.5_R8P*(q_aux_gpu(b,i,ll,k,2)**2+q_aux_gpu(b,i,ll,k,3)**2+q_aux_gpu(b,i,ll,k,4)**2)+ &
                q_aux_gpu(b,i,ll,k,5)*dha)
            fi(1) = vi(3)
            fi(2) = fi(1) * q_aux_gpu(b,i,ll,k,2)
            fi(3) = fi(1) * q_aux_gpu(b,i,ll,k,3) + q_aux_gpu(b,i,ll,k,7)
            fi(4) = fi(1) * q_aux_gpu(b,i,ll,k,4)
            fi(5) = fi(1) * vi(5) / vi(1) + q_aux_gpu(b,i,ll,k,7)*q_aux_gpu(b,i,ll,k,3)
            do m=1,nv
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,nv
                  wc = wc + el(mm,m) * vi(mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,i,k,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,i,k,b) = gc - gplus(m,l,i,k,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         wenorec_scheme = cell_scheme_gpu(b,i,j,k,2)
         call weno_reconstruction_kernel(nvar=nv, vp=gplus(1,1,i,k,b), vm=gminus(1,1,i,k,b), vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)

         ror_y: do m = 2, size(ror_schemes_gpu)
            ror_to_recompute = .false.
            do mm = 1,size(ror_ivar_gpu)
                index_var = ror_ivar_gpu(mm)
                if((abs(gl(index_var)-gplus(index_var,iweno,i,k,b))    > ror_threshold*abs(gplus(index_var,iweno,i,k,b)   )) .or. &
                   (abs(gr(index_var)-gminus(index_var,iweno+1,i,k,b)) > ror_threshold*abs(gminus(index_var,iweno+1,i,k,b))) ) then
                   ror_to_recompute = .true.
                endif
            enddo
            if(ror_to_recompute) then
               wenorec_scheme = ror_schemes_gpu(m)
               call weno_reconstruction_kernel(nvar=nv, vp=gplus(1,1,i,k,b), vm=gminus(1,1,i,k,b), vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)
            else
               exit ror_y
            endif
         enddo ror_y

         if (enable_ror_stats) ror_stats_gpu(b,i,j,k,2) = wenorec_scheme

         ! Reassemble + and - characteristic fluxes
         do m=1,nv
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,nv
            fly_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,nv
               fly_gpu(b,i,j,k,m) = fly_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_flux_conv_y_kernel

   attributes(global) subroutine compute_flux_conv_z_kernel(blocks_number, ni, nj, nk, ngc, nv, iweno, dha, g, R, cv,       &
                                                            ror_threshold, enable_ror_stats, cell_scheme_gpu, ror_ivar_gpu, &
                                                            ror_schemes_gpu, q_aux_gpu, ror_stats_gpu, gplus, gminus, flz_gpu)
   !< Compute convective fluxes by means of upwind WENO reconstruction, z axis direction.
   integer,      intent(in),    value  :: blocks_number, ni, nj, nk, ngc, nv, iweno
   real(R8P),    intent(in),    value  :: dha, g, R, cv
   real(R8P),    intent(in),    value  :: ror_threshold
   logical,      intent(in),    value  :: enable_ror_stats
   integer(I4P), intent(in),    device :: cell_scheme_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   integer(I4P), intent(in),    device :: ror_ivar_gpu(1:)
   integer(I4P), intent(in),    device :: ror_schemes_gpu(1:)
   real(R8P),    intent(in),    device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer(I4P), intent(inout), device :: ror_stats_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)
   real(R8P),    intent(inout), device ::  gplus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout), device :: gminus(1:, 1:, 1:, 1:, 1:)
   real(R8P),    intent(inout), device :: flz_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer                             :: b, i, j, k, l, ll, m, mm, v
   real(R8P)                           :: er(5,5), el(5,5), ev(5), evmax(5), ghat(5), gl(5), gr(5), fi(5), vi(5)
   real(R8P)                           :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                           :: gc, wc
   integer                             :: wenorec_scheme, index_var
   logical                             :: ror_to_recompute

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if (b > blocks_number .or. i > ni) return
   do j=1,nj
      do k=0,nk ! loop on faces
         ! Compute Roe average
         call compute_roe_average(q_aux_gpu=q_aux_gpu, dha=dha, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i, jp=j, kp=k+1, &
                                  uu=uu, vv=vv, ww=ww, h=h, ya=ya, qq=qq, c=c, ci=ci, b1=b1, b2=b2)
         ! Compute right and left eigenvectors matrices (at Roe state)
         er(1,1)=1._R8P ; er(1,2)=uu     ; er(1,3)=vv     ; er(1,4)=ww-c   ; er(1,5)=h-ww*c
         er(2,1)=1._R8P ; er(2,2)=uu     ; er(2,3)=vv     ; er(2,4)=ww     ; er(2,5)=qq
         er(3,1)=1._R8P ; er(3,2)=uu     ; er(3,3)=vv     ; er(3,4)=ww+c   ; er(3,5)=h+ww*c
         er(4,1)=0._R8P ; er(4,2)=1._R8P ; er(4,3)=0._R8P ; er(4,4)=0._R8P ; er(4,5)=uu
         er(5,1)=0._R8P ; er(5,2)=0._R8P ; er(5,3)=1._R8P ; er(5,4)=0._R8P ; er(5,5)=vv

         el(1,1)=0.5_R8P*(b1+ww*ci)  ; el(1,2)=1._R8P-b1 ; el(1,3)=0.5_R8P*(b1-ww*ci)  ; el(1,4)=-uu    ; el(1,5)=-vv
         el(2,1)=-0.5_R8P*(b2*uu)    ; el(2,2)=b2*uu     ; el(2,3)=-0.5_R8P*(b2*uu)    ; el(2,4)=1._R8P ; el(2,5)=0._R8P
         el(3,1)=-0.5_R8P*(b2*vv)    ; el(3,2)=b2*vv     ; el(3,3)=-0.5_R8P*(b2*vv)    ; el(3,4)=0._R8P ; el(3,5)=1._R8P
         el(4,1)=-0.5_R8P*(b2*ww+ci) ; el(4,2)=b2*ww     ; el(4,3)=-0.5_R8P*(b2*ww-ci) ; el(4,4)=0._R8P ; el(4,5)=0._R8P
         el(5,1)=0.5_R8P*b2          ; el(5,2)=-b2       ; el(5,3)=0.5_R8P*b2          ; el(5,4)=0._R8P ; el(5,5)=0._R8P

         ! Find max eigenvalues on the stencil
         do m=1,nv  ! loop on characteristic fields
            evmax(m) = -1._R8P
         enddo
         do l=1,2*iweno ! LLF
            ll = k + l - iweno
            uu = q_aux_gpu(b,i,j,ll,4)
            c  = q_aux_gpu(b,i,j,ll,9)
            ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2)
            do m=1,nv
               evmax(m) = max(ev(m),evmax(m))
            enddo
         enddo

         ! Decompose fluxes as + and -
         do l=1,2*iweno ! loop over the stencil centered at face i
            ll = k + l - iweno
            vi(1) = q_aux_gpu(b,i,j,ll,1)
            vi(2) = vi(1)*q_aux_gpu(b,i,j,ll,2)
            vi(3) = vi(1)*q_aux_gpu(b,i,j,ll,3)
            vi(4) = vi(1)*q_aux_gpu(b,i,j,ll,4)
            vi(5) = vi(1)*(cv*q_aux_gpu(b,i,j,ll,6)+                                                  &
                0.5_R8P*(q_aux_gpu(b,i,j,ll,2)**2+q_aux_gpu(b,i,j,ll,3)**2+q_aux_gpu(b,i,j,ll,4)**2)+ &
                q_aux_gpu(b,i,j,ll,5)*dha)
            fi(1) = vi(4)
            fi(2) = fi(1) * q_aux_gpu(b,i,j,ll,2)
            fi(3) = fi(1) * q_aux_gpu(b,i,j,ll,3)
            fi(4) = fi(1) * q_aux_gpu(b,i,j,ll,4) + q_aux_gpu(b,i,j,ll,7)
            fi(5) = fi(1) * vi(5) / vi(1) + q_aux_gpu(b,i,j,ll,7)*q_aux_gpu(b,i,j,ll,4)
            do m=1,nv
               wc = 0._R8P
               gc = 0._R8P
               do mm=1,nv
                  wc = wc + el(mm,m) * vi(mm)
                  gc = gc + el(mm,m) * fi(mm)
               enddo
               gplus (m,l,i,j,b) = 0.5_R8P * (gc + evmax(m) * wc)
               gminus(m,l,i,j,b) = gc - gplus(m,l,i,j,b)
            enddo
         enddo

         ! Reconstruction of the + and - fluxes
         wenorec_scheme = cell_scheme_gpu(b,i,j,k,3)
         call weno_reconstruction_kernel(nvar=nv, vp=gplus(1,1,i,j,b), vm=gminus(1,1,i,j,b), vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)

         ror_z: do m = 2, size(ror_schemes_gpu)
            ror_to_recompute = .false.
            do mm = 1,size(ror_ivar_gpu)
                index_var = ror_ivar_gpu(mm)
                if((abs(gl(index_var)-gplus(index_var,iweno,i,j,b))    > ror_threshold*abs(gplus(index_var,iweno,i,j,b)   )) .or. &
                   (abs(gr(index_var)-gminus(index_var,iweno+1,i,j,b)) > ror_threshold*abs(gminus(index_var,iweno+1,i,j,b))) ) then
                   ror_to_recompute = .true.
                endif
            enddo
            if(ror_to_recompute) then
               wenorec_scheme = ror_schemes_gpu(m)
               call weno_reconstruction_kernel(nvar=nv, vp=gplus(1,1,i,j,b), vm=gminus(1,1,i,j,b), vminus=gl, vplus=gr, iweno=iweno, wenorec_ord=wenorec_scheme)
            else
               exit ror_z
            endif
         enddo ror_z

         if (enable_ror_stats) ror_stats_gpu(b,i,j,k,3) = wenorec_scheme

         ! Reassemble + and - characteristic fluxes
         do m=1,nv
            ghat(m) = gl(m) + gr(m)
         enddo

         ! Return to conservative fluxes
         do m=1,nv
            flz_gpu(b,i,j,k,m) = 0._R8P
            do mm=1,nv
               flz_gpu(b,i,j,k,m) = flz_gpu(b,i,j,k,m) + er(mm,m) * ghat(mm)
            enddo
         enddo
      enddo
   enddo
   endsubroutine compute_flux_conv_z_kernel

   attributes(global) subroutine compute_flux_conv_x_central_kernel(blocks_number, ni, nj, nk, ngc, nv, lmax, &
                                                                    fc_coeff_gpu, q_aux_gpu, dx_gpu, flx_gpu)
   !< Compute convective fluxes by means of central WENO reconstruction, x axis direction.
   integer(I4P), intent(in),    value  :: blocks_number, ni, nj, nk, ngc, nv, lmax
   real(R8P),    intent(in),    device :: fc_coeff_gpu(1:, 1:)
   real(R8P),    intent(in),    device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(in),    device :: dx_gpu(1:)
   real(R8P),    intent(inout), device :: flx_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P)                           :: rhom, uui, vvi, wwi, ppi, enti, yai
   real(R8P)                           :: uuip, vvip, wwip, ppip, entip, yaip
   real(R8P)                           :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
   real(R8P)                           :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
   integer                             :: b, i, j, k, l, v, m

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   j = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if (b > blocks_number .or. j > nj) return
   do k=1,nk
      do i=0,ni ! loop on faces
         ft1 = 0._R8P ; ft2 = 0._R8P ; ft3 = 0._R8P ; ft4 = 0._R8P ; ft5 = 0._R8P ; ft6 = 0._R8P ; ft7 = 0._R8P
         do l=1,lmax
             uvs1 = 0._R8P ; uvs2 = 0._R8P ; uvs3 = 0._R8P ; uvs4 = 0._R8P ; uvs5 = 0._R8P ; uvs6 = 0._R8P ; uvs7 = 0._R8P
             do m=0,l-1
                 rhom    = q_aux_gpu(b,i-m,j,k,1) + q_aux_gpu(b,i-m+l,j,k,1)

                 uui     = q_aux_gpu(b,i-m,j,k,2)
                 vvi     = q_aux_gpu(b,i-m,j,k,3)
                 wwi     = q_aux_gpu(b,i-m,j,k,4)
                 ppi     = q_aux_gpu(b,i-m,j,k,7)
                 enti    = q_aux_gpu(b,i-m,j,k,8)
                 yai     = q_aux_gpu(b,i-m,j,k,5)

                 uuip    = q_aux_gpu(b,i-m+l,j,k,2)
                 vvip    = q_aux_gpu(b,i-m+l,j,k,3)
                 wwip    = q_aux_gpu(b,i-m+l,j,k,4)
                 ppip    = q_aux_gpu(b,i-m+l,j,k,7)
                 entip   = q_aux_gpu(b,i-m+l,j,k,8)
                 yaip    = q_aux_gpu(b,i-m+l,j,k,5)

                 uv_part = (uui+uuip) * rhom
                 uvs1    = uvs1 + uv_part * (2._R8P)
                 uvs2    = uvs2 + uv_part * (uui+uuip)
                 uvs3    = uvs3 + uv_part * (vvi+vvip)
                 uvs4    = uvs4 + uv_part * (wwi+wwip)
                 uvs5    = uvs5 + uv_part * (enti+entip)
                 uvs6    = uvs6 + (2._R8P)*(ppi+ppip)
                 uvs7    = uvs7 + uv_part * (yai+yaip)
             enddo
             ft1 = ft1 + fc_coeff_gpu(l,lmax)*uvs1
             ft2 = ft2 + fc_coeff_gpu(l,lmax)*uvs2
             ft3 = ft3 + fc_coeff_gpu(l,lmax)*uvs3
             ft4 = ft4 + fc_coeff_gpu(l,lmax)*uvs4
             ft5 = ft5 + fc_coeff_gpu(l,lmax)*uvs5
             ft6 = ft6 + fc_coeff_gpu(l,lmax)*uvs6
             ft7 = ft7 + fc_coeff_gpu(l,lmax)*uvs7
         enddo
         flx_gpu(b,i,j,k,1) = 0.25_R8P*ft1
         flx_gpu(b,i,j,k,2) = 0.25_R8P*ft2 + 0.5_R8P*ft6
         flx_gpu(b,i,j,k,3) = 0.25_R8P*ft3
         flx_gpu(b,i,j,k,4) = 0.25_R8P*ft4
         flx_gpu(b,i,j,k,5) = 0.25_R8P*ft5
         if(nv == 6) flx_gpu(b,i,j,k,6) = 0.25_R8P*ft7
      enddo
   enddo
   endsubroutine compute_flux_conv_x_central_kernel

   attributes(global) subroutine compute_flux_conv_y_central_kernel(blocks_number, ni, nj, nk, ngc, nv, lmax, &
                                                                    fc_coeff_gpu, q_aux_gpu, dy_gpu, fly_gpu)
   !< Compute convective fluxes by means of central WENO reconstruction, y axis direction.
   integer(I4P), intent(in),    value  :: blocks_number, ni, nj, nk, ngc, nv, lmax
   real(R8P),    intent(in),    device :: fc_coeff_gpu(1:, 1:)
   real(R8P),    intent(in),    device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(in),    device :: dy_gpu(1:)
   real(R8P),    intent(inout), device :: fly_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P)                           :: rhom, uui, vvi, wwi, ppi, enti, yai
   real(R8P)                           :: uuip, vvip, wwip, ppip, entip, yaip
   real(R8P)                           :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
   real(R8P)                           :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
   integer                             :: b, i, j, k, l, v, m

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. i > ni) return
   do k=1,nk
      do j=0,nj ! loop on faces
         ft1 = 0._R8P ; ft2 = 0._R8P ; ft3 = 0._R8P ; ft4 = 0._R8P ; ft5 = 0._R8P ; ft6 = 0._R8P ; ft7 = 0._R8P
         do l=1,lmax
             uvs1 = 0._R8P ; uvs2 = 0._R8P ; uvs3 = 0._R8P ; uvs4 = 0._R8P ; uvs5 = 0._R8P ; uvs6 = 0._R8P ; uvs7 = 0._R8P
             do m=0,l-1
                 rhom    = q_aux_gpu(b,i,j-m,k,1) + q_aux_gpu(b,i,j-m+l,k,1)

                 uui     = q_aux_gpu(b,i,j-m,k,2)
                 vvi     = q_aux_gpu(b,i,j-m,k,3)
                 wwi     = q_aux_gpu(b,i,j-m,k,4)
                 ppi     = q_aux_gpu(b,i,j-m,k,7)
                 enti    = q_aux_gpu(b,i,j-m,k,8)
                 yai     = q_aux_gpu(b,i,j-m,k,5)

                 uuip    = q_aux_gpu(b,i,j-m+l,k,2)
                 vvip    = q_aux_gpu(b,i,j-m+l,k,3)
                 wwip    = q_aux_gpu(b,i,j-m+l,k,4)
                 ppip    = q_aux_gpu(b,i,j-m+l,k,7)
                 entip   = q_aux_gpu(b,i,j-m+l,k,8)
                 yaip    = q_aux_gpu(b,i,j-m+l,k,5)

                 uv_part = (vvi+vvip) * rhom
                 uvs1    = uvs1 + uv_part * (2._R8P)
                 uvs2    = uvs2 + uv_part * (uui+uuip)
                 uvs3    = uvs3 + uv_part * (vvi+vvip)
                 uvs4    = uvs4 + uv_part * (wwi+wwip)
                 uvs5    = uvs5 + uv_part * (enti+entip)
                 uvs6    = uvs6 + (2._R8P)*(ppi+ppip)
                 uvs7    = uvs7 + uv_part * (yai+yaip)
             enddo
             ft1 = ft1 + fc_coeff_gpu(l,lmax)*uvs1
             ft2 = ft2 + fc_coeff_gpu(l,lmax)*uvs2
             ft3 = ft3 + fc_coeff_gpu(l,lmax)*uvs3
             ft4 = ft4 + fc_coeff_gpu(l,lmax)*uvs4
             ft5 = ft5 + fc_coeff_gpu(l,lmax)*uvs5
             ft6 = ft6 + fc_coeff_gpu(l,lmax)*uvs6
             ft7 = ft7 + fc_coeff_gpu(l,lmax)*uvs7
         enddo
         fly_gpu(b,i,j,k,1) = 0.25_R8P*ft1
         fly_gpu(b,i,j,k,2) = 0.25_R8P*ft2
         fly_gpu(b,i,j,k,3) = 0.25_R8P*ft3 + 0.5_R8P*ft6
         fly_gpu(b,i,j,k,4) = 0.25_R8P*ft4
         fly_gpu(b,i,j,k,5) = 0.25_R8P*ft5
         if(nv == 6) fly_gpu(b,i,j,k,6) = 0.25_R8P*ft7
      enddo
   enddo
   endsubroutine compute_flux_conv_y_central_kernel

   attributes(global) subroutine compute_flux_conv_z_central_kernel(blocks_number, ni, nj, nk, ngc, nv, lmax, &
                                                                    fc_coeff_gpu, q_aux_gpu, dz_gpu, flz_gpu)
   !< Compute convective fluxes by means of central WENO reconstruction, z axis direction.
   integer(I4P), intent(in),    value  :: blocks_number, ni, nj, nk, ngc, nv, lmax
   real(R8P),    intent(in),    device :: fc_coeff_gpu(1:, 1:)
   real(R8P),    intent(in),    device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(in),    device :: dz_gpu(1:)
   real(R8P),    intent(inout), device :: flz_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P)                           :: rhom, uui, vvi, wwi, ppi, enti, yai
   real(R8P)                           :: uuip, vvip, wwip, ppip, entip, yaip
   real(R8P)                           :: ft1, ft2, ft3, ft4, ft5, ft6, ft7
   real(R8P)                           :: uvs1, uvs2, uvs3, uvs4, uvs5, uvs6, uvs7, uv_part
   integer                             :: b, i, j, k, l, v, m

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if(b > blocks_number .or. i > ni) return
   do j=1,nj
      do k=0,nk ! loop on faces
         ft1 = 0._R8P ; ft2 = 0._R8P ; ft3 = 0._R8P ; ft4 = 0._R8P ; ft5 = 0._R8P ; ft6 = 0._R8P ; ft7 = 0._R8P
         do l=1,lmax
             uvs1 = 0._R8P ; uvs2 = 0._R8P ; uvs3 = 0._R8P ; uvs4 = 0._R8P ; uvs5 = 0._R8P ; uvs6 = 0._R8P ; uvs7 = 0._R8P
             do m=0,l-1
                 rhom    = q_aux_gpu(b,i,j,k-m,1) + q_aux_gpu(b,i,j,k-m+l,1)

                 uui     = q_aux_gpu(b,i,j,k-m,2)
                 vvi     = q_aux_gpu(b,i,j,k-m,3)
                 wwi     = q_aux_gpu(b,i,j,k-m,4)
                 ppi     = q_aux_gpu(b,i,j,k-m,7)
                 enti    = q_aux_gpu(b,i,j,k-m,8)
                 yai     = q_aux_gpu(b,i,j,k-m,5)

                 uuip    = q_aux_gpu(b,i,j,k-m+l,2)
                 vvip    = q_aux_gpu(b,i,j,k-m+l,3)
                 wwip    = q_aux_gpu(b,i,j,k-m+l,4)
                 ppip    = q_aux_gpu(b,i,j,k-m+l,7)
                 entip   = q_aux_gpu(b,i,j,k-m+l,8)
                 yaip    = q_aux_gpu(b,i,j,k-m+l,5)

                 uv_part = (wwi+wwip) * rhom
                 uvs1    = uvs1 + uv_part * (2._R8P)
                 uvs2    = uvs2 + uv_part * (uui+uuip)
                 uvs3    = uvs3 + uv_part * (vvi+vvip)
                 uvs4    = uvs4 + uv_part * (wwi+wwip)
                 uvs5    = uvs5 + uv_part * (enti+entip)
                 uvs6    = uvs6 + (2._R8P)*(ppi+ppip)
                 uvs7    = uvs7 + uv_part * (yai+yaip)
             enddo
             ft1 = ft1 + fc_coeff_gpu(l,lmax)*uvs1
             ft2 = ft2 + fc_coeff_gpu(l,lmax)*uvs2
             ft3 = ft3 + fc_coeff_gpu(l,lmax)*uvs3
             ft4 = ft4 + fc_coeff_gpu(l,lmax)*uvs4
             ft5 = ft5 + fc_coeff_gpu(l,lmax)*uvs5
             ft6 = ft6 + fc_coeff_gpu(l,lmax)*uvs6
             ft7 = ft7 + fc_coeff_gpu(l,lmax)*uvs7
         enddo
         flz_gpu(b,i,j,k,1) = 0.25_R8P*ft1
         flz_gpu(b,i,j,k,2) = 0.25_R8P*ft2
         flz_gpu(b,i,j,k,3) = 0.25_R8P*ft3
         flz_gpu(b,i,j,k,4) = 0.25_R8P*ft4 + 0.5_R8P*ft6
         flz_gpu(b,i,j,k,5) = 0.25_R8P*ft5
         if(nv == 6) flz_gpu(b,i,j,k,6) = 0.25_R8P*ft7
      enddo
   enddo
   endsubroutine compute_flux_conv_z_central_kernel

   subroutine compute_fluxes_difference_cuf(blocks_number, ni, nj, nk, ngc, nv,         &
                                            fl_gpu, flx_gpu, fly_gpu, flz_gpu, phi_gpu, &
                                            dx_gpu, dy_gpu, dz_gpu, ib_eps)
   !< Compute fluxes difference.
   integer(I4P), intent(in)            :: blocks_number, ni, nj, nk, ngc, nv
   real(R8P),    intent(inout), device ::  fl_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(in)   , device :: flx_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(in)   , device :: fly_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(in)   , device :: flz_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(in)   , device :: phi_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(in)   , device :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)
   real(R8P),    intent(in)            :: ib_eps
   real(R8P)                           :: delta_x, delta_y, delta_z, dx_locale, dy_locale, dz_locale
   integer(I4P)                        :: b, i, j, k, v, iercuda

   !$cuf kernel do(4) <<<*,*>>>
   do k=1,nk
   do j=1,nj
   do i=1,ni
   do b=1,blocks_number

      dx_locale = dx_gpu(b)
      ! Update net flux (procedura alternativa all'interpolazione proposta nel paper, utilizza dx_locale).
      if(phi_gpu(b,i,j,k,1)<0.) then
          if(phi_gpu(b,i+1,j,k,1)*phi_gpu(b,i-1,j,k,1)<0) then
              if(phi_gpu(b,i+1,j,k,1)>0.) then
                  delta_x = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i+1,j,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dx_gpu(b)
                  dx_locale = dx_gpu(b)/2 + delta_x
              else !if(phi_gpu(b,i-1,j,k,1)>0) then
                  delta_x = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i-1,j,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dx_gpu(b)
                  dx_locale = dx_gpu(b)/2 + delta_x
              endif
          endif
      endif

      dy_locale = dy_gpu(b)
      if(phi_gpu(b,i,j,k,1)<0.) then
          if(phi_gpu(b,i,j+1,k,1)*phi_gpu(b,i,j-1,k,1)<0) then
              if(phi_gpu(b,i,j+1,k,1)>0.) then
                  delta_y = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j+1,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dy_gpu(b)
                  dy_locale = dy_gpu(b)/2 + delta_y
              else !if(phi_gpu(b,i-1,j,k,1)>0) then
                  delta_y = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j-1,k,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dy_gpu(b)
                  dy_locale = dy_gpu(b)/2 + delta_y
              endif
          endif
      endif

      dz_locale = dz_gpu(b)
      if(phi_gpu(b,i,j,k,1)<0.) then
          if(phi_gpu(b,i,j,k+1,1)*phi_gpu(b,i,j,k-1,1)<0) then
              if(phi_gpu(b,i,j,k+1,1)>0.) then
                  delta_z = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j,k+1,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dz_gpu(b)
                  dz_locale = dz_gpu(b)/2 + delta_z
              else !if(phi_gpu(b,i,j,k-1,1)>0) then
                  delta_z = -phi_gpu(b,i,j,k,1)/(phi_gpu(b,i,j,k-1,1)-phi_gpu(b,i,j,k,1)+ib_eps)*dz_gpu(b)
                  dz_locale = dz_gpu(b)/2 + delta_z
              endif
          endif
      endif

      do v=1,nv
         fl_gpu(b,i,j,k,v) = - (flx_gpu(b,i,j,k,v)-flx_gpu(b,i-1,j,k,v))/dx_locale &
                             - (fly_gpu(b,i,j,k,v)-fly_gpu(b,i,j-1,k,v))/dy_locale &
                             - (flz_gpu(b,i,j,k,v)-flz_gpu(b,i,j,k-1,v))/dz_locale
      enddo

   enddo
   enddo
   enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_fluxes_difference_cuf

   subroutine compute_fluxes_diffusive_cuf(blocks_number, ni, nj, nk, ngc, nv, mu, kd, &
                                           q_aux_gpu, dx_gpu, dy_gpu, dz_gpu, flx_gpu, fly_gpu, flz_gpu)
   !< Compute diffusive fluxes.
   integer(I4P), intent(in)            :: blocks_number, ni, nj, nk, ngc, nv
   real(R8P),    intent(in)            :: mu, kd
   real(R8P),    intent(in),    device :: dx_gpu(1:), dy_gpu(1:), dz_gpu(1:)
   real(R8P),    intent(in),    device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(inout), device ::   flx_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(inout), device ::   fly_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(inout), device ::   flz_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   integer(I4P)                        :: b, i, j, k, v, iercuda
   real(R8P)                           :: du_dx, dv_dx, dw_dx, du_dy, dv_dy, dw_dy, du_dz, dv_dz, dw_dz
   real(R8P)                           :: dx_locale, dy_locale, dz_locale
   real(R8P)                           :: delta_x, delta_y, delta_z
   real(R8P)                           :: sigq, sigl
   real(R8P)                           :: tau_1_1, tau_2_1, tau_3_1, dT_dx
   real(R8P)                           :: tau_1_2, tau_2_2, tau_3_2, dT_dy
   real(R8P)                           :: tau_1_3, tau_2_3, tau_3_3, dT_dz
   real(R8P)                           :: vel_u, vel_v, vel_w
   real(R8P), parameter                :: ib_eps=1.e-12_R8P

   !$cuf kernel do(3) <<<*,*>>>
   do k=1,nk
      do j=1,nj
         do b=1,blocks_number
            do i=0,ni ! loop on faces
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
            do j=0,nj ! loop on faces
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
            do k=0,nk ! loop on faces
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
   !< Gradient done by CUF threads.
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

   subroutine compute_umax_cuf(b, ni, nj, nk, ngc, ns, dx, dy, dz, mu, q_aux_gpu, umax)
   !< Compute maximum speed.
   integer(I4P), intent(in)         :: b                                     !< Block index.
   integer(I4P), intent(in)         :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)         :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)         :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)         :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)         :: ns                                    !< Number of species.
   real(R8P),    intent(in)         :: dx                                    !< X space step.
   real(R8P),    intent(in)         :: dy                                    !< Y space step.
   real(R8P),    intent(in)         :: dz                                    !< Z space step.
   real(R8P),    intent(in)         :: mu                                    !< Dynamic viscosity.
   real(R8P),    intent(in), device :: q_aux_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary varibales.
   real(R8P),    intent(out)        :: umax                                  !< Maximum speed.
   real(R8P)                        :: ss                                    !< Speed of sound.
   integer(I4P)                     :: i, j, k                               !< Counter.
   integer(I4P)                     :: iercuda                               !< Error trapping flag for CUDAFortran.
   real(R8P)                        :: dx_locale, dy_locale, dz_locale       !< Local space steps.

   umax = 0._R8P
   dx_locale = dx*0.5_R8P
   dy_locale = dy*0.5_R8P
   dz_locale = dz*0.5_R8P
   !$cuf kernel do(3) <<<*,*>>>
   do k=1, nk
      do j=1, nj
         do i=1, ni
            ss = q_aux_gpu(b,i,j,k,9)
            umax = max(umax, (abs(q_aux_gpu(b,i,j,k,2)) + ss)/dx_locale + 2._R8P*mu/(q_aux_gpu(b,i,j,k,1))/dx_locale**2 + &
                             (abs(q_aux_gpu(b,i,j,k,3)) + ss)/dy_locale + 2._R8P*mu/(q_aux_gpu(b,i,j,k,1))/dy_locale**2 + &
                             (abs(q_aux_gpu(b,i,j,k,4)) + ss)/dz_locale + 2._R8P*mu/(q_aux_gpu(b,i,j,k,1))/dz_locale**2)
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_umax_cuf

   attributes(device) subroutine compute_roe_average(q_aux_gpu, dha, g, &
                                                     ngc, b, i, j, k, ip, jp, kp, uu, vv, ww, h, ya, qq, c, ci, b1, b2)
   !< Compute Roe averaged quantities.
   real(R8P),    intent(in), device :: q_aux_gpu(1:, 1-ngc:, 1-ngc:, 1-ngc:, 1:)
   real(R8P),    intent(in)         :: dha, g
   integer(I4P), intent(in)         :: ngc, b, i, j, k, ip, jp, kp
   real(R8P),    intent(out)        :: uu, vv, ww, h, ya, qq, c, ci, b1, b2
   real(R8P)                        :: ri, up, vp, wp, hp, yap, r, rp1, cc
   ! Left state (node i)
   ri        =  1._R8P/q_aux_gpu(b,i,j,k,1)
   uu        =  q_aux_gpu(b,i,j,k,2)
   vv        =  q_aux_gpu(b,i,j,k,3)
   ww        =  q_aux_gpu(b,i,j,k,4)
   h         =  q_aux_gpu(b,i,j,k,8)
   ya        =  q_aux_gpu(b,i,j,k,5)
   ! Right state (node i+1)
   up        =  q_aux_gpu(b,ip,jp,kp,2)
   vp        =  q_aux_gpu(b,ip,jp,kp,3)
   wp        =  q_aux_gpu(b,ip,jp,kp,4)
   hp        =  q_aux_gpu(b,ip,jp,kp,8)
   yap       =  q_aux_gpu(b,ip,jp,kp,5)
   ! Average state
   r         =  sqrt(q_aux_gpu(b,ip,jp,kp,1)*ri)
   rp1       =  1._R8P/(r+1._R8P)
   uu        =  (r*up+uu)*rp1
   vv        =  (r*vp+vv)*rp1
   ww        =  (r*wp+ww)*rp1
   h         =  (r*hp+h)*rp1
   ya        =  (r*yap+ya)*rp1
   qq        =  0.5_R8P * (uu*uu+vv*vv+ww*ww)
   cc        =  (g-1._R8P) * (h - qq - ya*dha)
   !ERRATODIREIcc        =  g * (g-1._R8P) * (h - qq - ya*dha)
   c         =  sqrt(cc)
   ci        =  1._R8P/c
   b2        = (g-1)/cc  ! alias 1/(cp*theta)
   b1        = b2 * qq   ! alias q/(cp*theta)

   endsubroutine compute_roe_average

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

   ! eikonal procedures
   subroutine evolve_eikonal_q_gpu_cuf(ni, nj, nk, ngc, nv, phi_gpu, dx_gpu, dy_gpu, dz_gpu, blocks_number, dq_gpu, q_gpu)
   !< Evolve eikonal equation over q.
   integer(I4P), intent(in)            :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                  !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in),    device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(in),    device :: dx_gpu(1:)                          !< X space steps.
   real(R8P),    intent(in),    device :: dy_gpu(1:)                          !< Y space steps.
   real(R8P),    intent(in),    device :: dz_gpu(1:)                          !< Z space steps.
   real(R8P),    intent(inout), device :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !<
   real(R8P),    intent(inout), device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Conservative variables.
   integer(I4P)                        :: i, j, k, b, v                       !< Counter.
   integer(I4P)                      :: iercuda                               !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(4) <<<*,*>>>
   do k=1, nk
      do j=1, nj
         do i=1,ni
            do b=1, blocks_number
               do v=1, nv
                  if (phi_gpu(b,i,j,k,1) > 0._R8P) then
                     q_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) - dq_gpu(b,i,j,k,v)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine evolve_eikonal_q_gpu_cuf

   subroutine invert_eikonal_q_gpu_cuf(BCS_VISCOUS, BCS_EULER, ni, nj, nk, ngc, nv, blocks_number, q_gpu, q_invert_gpu, phi_gpu, bcs_type)
   !< Invert eikonal field.
   integer(I4P), intent(in)            :: BCS_VISCOUS                               !< Viscous wall BCS parameter.
   integer(I4P), intent(in)            :: BCS_EULER                                 !< Euler wall BCS parameter.
   integer(I4P), intent(in)            :: ni                                        !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                        !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                        !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                       !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                        !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                             !< Number of blocks.
   integer(I4P), intent(in)            :: bcs_type                                  !< Immersed boundary type.
   real(R8P),    intent(in),    device ::  q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)        !< Conservative field.
   real(R8P),    intent(in),    device ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)      !< Distance field.
   real(R8P),    intent(inout), device ::  q_invert_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Inverted internal field.
   integer(I4P)                        :: i, j, k, b, v, ss                         !< Counter.
   integer(I4P)                        :: iercuda                                   !< Error trapping flag for CUDAFortran.
   real(R8P)                           :: n_phi_x, n_phi_y, n_phi_z                 !< Distance function normals.
   real(R8P)                           :: n_phi_mod, un_mod                         !< Distance abs normal and normal velocity.

   if(bcs_type == BCS_VISCOUS) then
      !$cuf kernel do(4) <<<*,*>>>
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               do b=1, blocks_number
                  if(phi_gpu(b,i,j,k,1) < 0) then
                      do v=1,nv
                         q_invert_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v)
                      enddo
                  else
                      q_invert_gpu(b,i,j,k,1) =   q_gpu(b,i,j,k,1)
                      q_invert_gpu(b,i,j,k,2) = - q_gpu(b,i,j,k,2)
                      q_invert_gpu(b,i,j,k,3) = - q_gpu(b,i,j,k,3)
                      q_invert_gpu(b,i,j,k,4) = - q_gpu(b,i,j,k,4)
                      q_invert_gpu(b,i,j,k,5) =   q_gpu(b,i,j,k,5)
                      if(nv == 6) q_invert_gpu(b,i,j,k,nv) = q_gpu(b,i,j,k,nv)
                  endif
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   elseif(bcs_type == BCS_EULER) then
      !$cuf kernel do(4) <<<*,*>>>
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               do b=1, blocks_number
                  if(phi_gpu(b,i,j,k,1) < 0) then
                     do v=1,nv
                        q_invert_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v)
                     enddo
                  else
                     n_phi_x = phi_gpu(b,i+1,j,k,1)-phi_gpu(b,i-1,j,k,1)
                     n_phi_y = phi_gpu(b,i,j+1,k,1)-phi_gpu(b,i,j-1,k,1)
                     n_phi_z = phi_gpu(b,i,j,k+1,1)-phi_gpu(b,i,j,k-1,1)
                     n_phi_mod = sqrt(n_phi_x**2+n_phi_y**2+n_phi_z**2)
                     n_phi_x = n_phi_x/n_phi_mod
                     n_phi_y = n_phi_y/n_phi_mod
                     n_phi_z = n_phi_z/n_phi_mod
                     un_mod = q_gpu(b,i,j,k,2)*n_phi_x+q_gpu(b,i,j,k,3)*n_phi_y+q_gpu(b,i,j,k,4)*n_phi_z

                     q_invert_gpu(b,i,j,k,1) = q_gpu(b,i,j,k,1)
                     q_invert_gpu(b,i,j,k,2) = q_gpu(b,i,j,k,2) - 2*un_mod*n_phi_x
                     q_invert_gpu(b,i,j,k,3) = q_gpu(b,i,j,k,3) - 2*un_mod*n_phi_y
                     q_invert_gpu(b,i,j,k,4) = q_gpu(b,i,j,k,4) - 2*un_mod*n_phi_z
                     q_invert_gpu(b,i,j,k,5) = q_gpu(b,i,j,k,5)
                     if(nv == 6) q_invert_gpu(b,i,j,k,6) = q_gpu(b,i,j,k,6)
                  endif
               enddo
            enddo
         enddo
      enddo
      !@cuf iercuda=cudaDeviceSynchronize()
   endif
   endsubroutine invert_eikonal_q_gpu_cuf

   attributes(global) subroutine compute_eikonal_dq_gpu(ni, nj, nk, ngc, nv, blocks_number, &
                                                        phi_gpu, dx_gpu, dy_gpu, dz_gpu,    &
                                                        q_gpu, dq_gpu)
   !< Compute eikonal dq.
   integer(I4P), intent(in), value     :: ni                                  !< Grid cells number in I direction.
   integer(I4P), intent(in), value     :: nj                                  !< Grid cells number in J direction.
   integer(I4P), intent(in), value     :: nk                                  !< Grid cells number in K direction.
   integer(I4P), intent(in), value     :: ngc                                 !< Ghost cells number.
   integer(I4P), intent(in), value     :: nv                                  !< Number of conservative varibales.
   integer(I4P), intent(in), value     :: blocks_number                       !< Number of blocks.
   real(R8P),    intent(in),    device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Distance function.
   real(R8P),    intent(in),    device :: dx_gpu(1:)                          !< X space steps.
   real(R8P),    intent(in),    device :: dy_gpu(1:)                          !< Y space steps.
   real(R8P),    intent(in),    device :: dz_gpu(1:)                          !< Z space steps.
   real(R8P),    intent(in),    device :: q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)   !< Conservative variables.
   real(R8P),    intent(inout), device :: dq_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !<
   integer(I4P)                        :: i, j, k, b, v                       !< Counter.
   real(R8P)                           :: n_phi_x, n_phi_y, n_phi_z, n_phi    !<

   b = blockDim%x * (blockIdx%x - 1) + threadIdx%x
   i = blockDim%y * (blockIdx%y - 1) + threadIdx%y
   if (b > blocks_number .or. i > ni) return

   do k=1, nk
      do j=1, nj
         if (phi_gpu(b,i,j,k,1) > 0._R8P) then
            n_phi_x = (phi_gpu(b,i+1,j,k,1) - phi_gpu(b,i-1,j,k,1) ) ! / (2 * dx_gpu(b))
            n_phi_y = (phi_gpu(b,i,j+1,k,1) - phi_gpu(b,i,j-1,k,1) ) ! / (2 * dy_gpu(b))
            n_phi_z = (phi_gpu(b,i,j,k+1,1) - phi_gpu(b,i,j,k-1,1) ) ! / (2 * dz_gpu(b))
            n_phi = abs(n_phi_x) + abs(n_phi_y) + abs(n_phi_z) + 10e-12
            n_phi = 0.9_R8P / n_phi
            n_phi_x = n_phi_x * n_phi
            n_phi_y = n_phi_y * n_phi
            n_phi_z = n_phi_z * n_phi
            do v=1, nv
               dq_gpu(b,i,j,k,v) = 0._R8P
            enddo
            if (n_phi_x > 0._R8P) then
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_x) * (q_gpu(b,i,j,k,v) - q_gpu(b,i-1,j,k,v))
               enddo
            else
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_x) * (q_gpu(b,i,j,k,v) - q_gpu(b,i+1,j,k,v))
               enddo
            endif
            if (n_phi_y > 0._R8P) then
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_y) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j-1,k,v))
               enddo
            else
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_y) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j+1,k,v))
               enddo
            endif
            if (n_phi_z > 0._R8P) then
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_z) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j,k-1,v))
               enddo
            else
               do v=1, nv
                  dq_gpu(b,i,j,k,v) = dq_gpu(b,i,j,k,v) + abs(n_phi_z) * (q_gpu(b,i,j,k,v) - q_gpu(b,i,j,k+1,v))
               enddo
            endif
         endif
      enddo
   enddo
   endsubroutine compute_eikonal_dq_gpu

   ! RK procedures
   subroutine compute_rk_linear_gpu_cuf(ni, nj, nk, ngc, nv, blocks_number, dt, q_gpu, prhs_gpu, fl_gpu, phi_gpu, qnrk)
   !< Compute RK linear stages.
   integer(I4P), intent(in)            :: ni                                   !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                   !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                   !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                  !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                   !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                        !< Number of blocks.
   real(R8P),    intent(in)            :: dt                                   !< Time step.
   real(R8P),    intent(in)            :: qnrk                                 !< Time step.
   real(R8P),    intent(inout), device ::   q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative field.
   real(R8P),    intent(in),    device ::  fl_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative field.
   real(R8P),    intent(in),    device :: phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:)  !< Conservative field.
   real(R8P),    intent(inout), device :: prhs_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< RK stage.
   integer(I4P)                        :: i, j, k, b, v, ss                    !< Counter.
   integer(I4P)                        :: iercuda                              !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(5) <<<*,*>>>
   do v=1, nv
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  if(phi_gpu(b,i,j,k,1) < 0.) then
                     q_gpu(b,i,j,k,v) = prhs_gpu(b,i,j,k,v) + qnrk * fl_gpu(b,i,j,k,v)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_rk_linear_gpu_cuf

   subroutine compute_rk_q_gpu_cuf(ni, nj, nk, ngc, nv, blocks_number, dt, s, q_gpu, q_old_gpu, fl_gpu, phi_gpu, ark, brk, crk)
   !< Compute RK approximation over q.
   integer(I4P), intent(in)            :: ni                                    !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                    !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                    !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                   !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                    !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                         !< Number of blocks.
   real(R8P),    intent(in)            :: dt                                    !< Time step.
   real(R8P),    intent(in)            :: ark, brk, crk                         !< Time step.
   integer(I4P), intent(in)            :: s                                     !< Stage to initialize.
   real(R8P),    intent(inout), device ::     q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P),    intent(in),    device ::    fl_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P),    intent(in),    device ::   phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P),    intent(in),    device :: q_old_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< RK stage.
   integer(I4P)                        :: i, j, k, b, v, ss                     !< Counter.
   integer(I4P)                        :: iercuda                               !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(5) <<<*,*>>>
   do v=1, nv
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  if(phi_gpu(b,i,j,k,1) < 0.) then
                     q_gpu(b,i,j,k,v) = ark * q_old_gpu(b,i,j,k,v) + brk * q_gpu(b,i,j,k,v) + dt * crk * fl_gpu(b,i,j,k,v)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_rk_q_gpu_cuf

   subroutine compute_rk_prhs_gpu_cuf(ni, nj, nk, ngc, nv, blocks_number, dt, s, q_gpu, prhs_gpu, fl_gpu, phi_gpu, qnrk)
   !< Compute RK approximation with Immersed Boundary.
   integer(I4P), intent(in)            :: ni                                   !< Grid cells number in I direction.
   integer(I4P), intent(in)            :: nj                                   !< Grid cells number in J direction.
   integer(I4P), intent(in)            :: nk                                   !< Grid cells number in K direction.
   integer(I4P), intent(in)            :: ngc                                  !< Ghost cells number.
   integer(I4P), intent(in)            :: nv                                   !< Number of conservative varibales.
   integer(I4P), intent(in)            :: blocks_number                        !< Number of blocks.
   real(R8P),    intent(in)            :: dt                                   !< Time step.
   real(R8P),    intent(in)            :: qnrk                                 !< Time step.
   integer(I4P), intent(in)            :: s                                    !< Stage to initialize.
   real(R8P),    intent(in),    device ::    q_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P),    intent(in),    device ::   fl_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P),    intent(in),    device ::  phi_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative field.
   real(R8P),    intent(inout), device :: prhs_gpu(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< RK stage.
   integer(I4P)                        :: i, j, k, b, v, ss                    !< Counter.
   integer(I4P)                        :: iercuda                              !< Error trapping flag for CUDAFortran.

   !$cuf kernel do(5) <<<*,*>>>
   do v=1, nv
      do k=1, nk
         do j=1, nj
            do i=1, ni
               do b=1, blocks_number
                  if(phi_gpu(b,i,j,k,1) < 0.) then
                     prhs_gpu(b,i,j,k,v) = q_gpu(b,i,j,k,v) + qnrk * fl_gpu(b,i,j,k,v)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo
   !@cuf iercuda=cudaDeviceSynchronize()
   endsubroutine compute_rk_prhs_gpu_cuf

   ! WENO procedures
   attributes(device) subroutine weno_reconstruction_kernel(nvar, vp, vm, vminus, vplus, iweno, wenorec_ord)
   !< Compute WENO reconstruction.
   integer, intent(in)                     :: nvar, iweno, wenorec_ord
   !real(R8P), dimension(nvar,2*iweno), intent(in)  :: vm,vp
   real(R8P), dimension(1:nvar,1:*)        :: vm,vp
   real(R8P), dimension(nvar), intent(out) :: vminus,vplus
   real(R8P), dimension(-1:4)              :: dwe               ! linear weights
   real(R8P), dimension(-1:4)              :: alfp,alfm         ! alpha_l
   real(R8P), dimension(-1:4)              :: alfp_map,alfm_map ! alpha_l
   real(R8P), dimension(-1:4)              :: betap,betam       ! beta_l
   real(R8P), dimension(-1:4)              :: omp,omm           ! WENO weights
   integer                                 :: r,i,j,k,l,m
   real(R8P)                               :: c0,c1,c2,c3,c4,d0,d1,d2,d3,d4,summ,sump
   real(R8P)                               :: x,y,y2

   if (wenorec_ord==1) then ! Godunov

       i = iweno ! index of intermediate node to perform reconstruction

       vminus(1:nvar) = vp(1:nvar,i)
       vplus (1:nvar) = vm(1:nvar,i+1)

   elseif (wenorec_ord==2) then ! WENO-3

       i = iweno ! index of intermediate node to perform reconstruction

       dwe(1)   = 2._R8P/3._R8P
       dwe(0)   = 1._R8P/3._R8P

       do m=1,nvar

           betap(0)  = (vp(m,i  )-vp(m,i-1))**2
           betap(1)  = (vp(m,i+1)-vp(m,i  ))**2
           betam(0)  = (vm(m,i+2)-vm(m,i+1))**2
           betam(1)  = (vm(m,i+1)-vm(m,i  ))**2

           sump = 0._R8P
           summ = 0._R8P
           do l=0,1
               alfp(l) = dwe(l)/(0.000001_R8P+betap(l))**2
               alfm(l) = dwe(l)/(0.000001_R8P+betam(l))**2
               sump = sump + alfp(l)
               summ = summ + alfm(l)
           enddo
           do l=0,1
               omp(l) = alfp(l)/sump
               omm(l) = alfm(l)/summ
           enddo

           vminus(m) = omp(0) *(-vp(m,i-1)+3*vp(m,i  )) + omp(1) *( vp(m,i  )+ vp(m,i+1))
           vplus(m)  = omm(0) *(-vm(m,i+2)+3*vm(m,i+1)) + omm(1) *( vm(m,i  )+ vm(m,i+1))

       enddo

       do m=1,nvar
           vminus(m) = 0.5_R8P*vminus(m)
           vplus(m)  = 0.5_R8P*vplus(m)
       enddo

     elseif (wenorec_ord==3) then ! WENO-5
!
      i = iweno ! index of intermediate node to perform reconstruction
!
      dwe( 0) = 1._R8P/10._R8P
      dwe( 1) = 6._R8P/10._R8P
      dwe( 2) = 3._R8P/10._R8P
!     JS
      d0 = 13._R8P/12._R8P
      d1 = 1._R8P/4._R8P
!     Weights for polynomial reconstructions
      c0 = 1._R8P/3._R8P
      c1 = 5._R8P/6._R8P
      c2 =-1._R8P/6._R8P
      c3 =-7._R8P/6._R8P
      c4 =11._R8P/6._R8P
!
      do m=1,nvar
!
       betap(2) = d0*(     vp(m,i)-2._R8P*vp(m,i+1)+vp(m,i+2))**2+d1*(3._R8P*vp(m,i)-4._R8P*vp(m,i+1)+vp(m,i+2))**2
       betap(1) = d0*(     vp(m,i-1)-2._R8P*vp(m,i)+vp(m,i+1))**2+d1*(     vp(m,i-1)-vp(m,i+1) )**2
       betap(0) = d0*(     vp(m,i)-2._R8P*vp(m,i-1)+vp(m,i-2))**2+d1*(3._R8P*vp(m,i)-4._R8P*vp(m,i-1)+vp(m,i-2))**2
!
       betam(2) = d0*(     vm(m,i+1)-2._R8P*vm(m,i)+vm(m,i-1))**2+d1*(3._R8P*vm(m,i+1)-4._R8P*vm(m,i)+vm(m,i-1))**2
       betam(1) = d0*(     vm(m,i+2)-2._R8P*vm(m,i+1)+vm(m,i))**2+d1*(     vm(m,i+2)-vm(m,i) )**2
       betam(0) = d0*(     vm(m,i+1)-2._R8P*vm(m,i+2)+vm(m,i+3))**2+d1*(3._R8P*vm(m,i+1)-4._R8P*vm(m,i+2)+vm(m,i+3))**2
!
       sump = 0._R8P
       summ = 0._R8P
       do l=0,2
        alfp(l) = dwe(  l)/(0.000001_R8P+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_R8P+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,2
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m)   = omp(2)*(c0*vp(m,i  )+c1*vp(m,i+1)+c2*vp(m,i+2)) + &
         & omp(1)*(c2*vp(m,i-1)+c1*vp(m,i  )+c0*vp(m,i+1)) + omp(0)*(c0*vp(m,i-2)+c3*vp(m,i-1)+c4*vp(m,i  ))
       vplus(m)   = omm(2)*(c0*vm(m,i+1)+c1*vm(m,i  )+c2*vm(m,i-1)) +  &
         & omm(1)*(c2*vm(m,i+2)+c1*vm(m,i+1)+c0*vm(m,i  )) + omm(0)*(c0*vm(m,i+3)+c3*vm(m,i+2)+c4*vm(m,i+1))
!
      enddo ! end of m-loop
!
   elseif (wenorec_ord==4) then ! WENO-7
!
      i = iweno ! index of intermediate node to perform reconstruction
!
      dwe( 0) = 1._R8P/35._R8P
      dwe( 1) = 12._R8P/35._R8P
      dwe( 2) = 18._R8P/35._R8P
      dwe( 3) = 4._R8P/35._R8P
!
!     JS weights
      d1 = 1._R8P/36._R8P
      d2 = 13._R8P/12._R8P
      d3 = 781._R8P/720._R8P
!
      do m=1,nvar
!
       betap(3)= d1*(-11*vp(m,  i)+18*vp(m,i+1)- 9*vp(m,i+2)+ 2*vp(m,i+3))**2+&
       &  d2*(  2*vp(m,  i)- 5*vp(m,i+1)+ 4*vp(m,i+2)-   vp(m,i+3))**2+ &
       & d3*(   -vp(m,  i)+ 3*vp(m,i+1)- 3*vp(m,i+2)+   vp(m,i+3))**2
       betap(2)= d1*(- 2*vp(m,i-1)- 3*vp(m,i  )+ 6*vp(m,i+1)-   vp(m,i+2))**2+&
       &  d2*(    vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1)             )**2+&
       &  d3*(   -vp(m,i-1)+ 3*vp(m,i  )- 3*vp(m,i+1)+   vp(m,i+2))**2
       betap(1)= d1*(    vp(m,i-2)- 6*vp(m,i-1)+ 3*vp(m,i  )+ 2*vp(m,i+1))**2+&
       &  d2*( vp(m,i-1)- 2*vp(m,i  )+   vp(m,i+1))**2+ &
       &  d3*(   -vp(m,i-2)+ 3*vp(m,i-1)- 3*vp(m,i  )+   vp(m,i+1))**2
       betap(0)= d1*(- 2*vp(m,i-3)+ 9*vp(m,i-2)-18*vp(m,i-1)+11*vp(m,i  ))**2+&
       &  d2*(-   vp(m,i-3)+ 4*vp(m,i-2)- 5*vp(m,i-1)+ 2*vp(m,i  ))**2+&
       &  d3*(   -vp(m,i-3)+ 3*vp(m,i-2)- 3*vp(m,i-1)+   vp(m,i  ))**2
!
       betam(3)= d1*(-11*vm(m,i+1)+18*vm(m,i  )- 9*vm(m,i-1)+ 2*vm(m,i-2))**2+&
       &  d2*(  2*vm(m,i+1)- 5*vm(m,i  )+ 4*vm(m,i-1)-   vm(m,i-2))**2+&
       &  d3*(   -vm(m,i+1)+ 3*vm(m,i  )- 3*vm(m,i-1)+   vm(m,i-2))**2
       betam(2)= d1*(- 2*vm(m,i+2)- 3*vm(m,i+1)+ 6*vm(m,i  )-   vm(m,i-1))**2+&
       &  d2*(    vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  )             )**2+&
       &  d3*(   -vm(m,i+2)+ 3*vm(m,i+1)- 3*vm(m,i  )+   vm(m,i-1))**2
       betam(1)= d1*(    vm(m,i+3)- 6*vm(m,i+2)+ 3*vm(m,i+1)+ 2*vm(m,i  ))**2+&
       &  d2*(                 vm(m,i+2)- 2*vm(m,i+1)+   vm(m,i  ))**2+&
       &  d3*(   -vm(m,i+3)+ 3*vm(m,i+2)- 3*vm(m,i+1)+   vm(m,i  ))**2
       betam(0)= d1*(- 2*vm(m,i+4)+ 9*vm(m,i+3)-18*vm(m,i+2)+11*vm(m,i+1))**2+&
       &  d2*(-   vm(m,i+4)+ 4*vm(m,i+3)- 5*vm(m,i+2)+ 2*vm(m,i+1))**2+&
       &  d3*(   -vm(m,i+4)+ 3*vm(m,i+3)- 3*vm(m,i+2)+   vm(m,i+1))**2
!
       sump = 0._R8P
       summ = 0._R8P
       do l=0,3
        alfp(l) = dwe(  l)/(0.000001_R8P+betap(l))**2
        alfm(l) = dwe(  l)/(0.000001_R8P+betam(l))**2
        sump = sump + alfp(l)
        summ = summ + alfm(l)
       enddo
       do l=0,3
        omp(l) = alfp(l)/sump
        omm(l) = alfm(l)/summ
       enddo
!
       vminus(m)   = omp(3)*( 6*vp(m,i  )+26*vp(m,i+1)-10*vp(m,i+2)+ 2*vp(m,i+3))+&
        omp(2)*(-2*vp(m,i-1)+14*vp(m,i  )+14*vp(m,i+1)- 2*vp(m,i+2))+&
        omp(1)*( 2*vp(m,i-2)-10*vp(m,i-1)+26*vp(m,i  )+ 6*vp(m,i+1))+&
        omp(0)*(-6*vp(m,i-3)+26*vp(m,i-2)-46*vp(m,i-1)+50*vp(m,i  ))
       vplus(m)   =  omm(3)*( 6*vm(m,i+1)+26*vm(m,i  )-10*vm(m,i-1)+ 2*vm(m,i-2))+&
        omm(2)*(-2*vm(m,i+2)+14*vm(m,i+1)+14*vm(m,i  )- 2*vm(m,i-1))+&
        omm(1)*( 2*vm(m,i+3)-10*vm(m,i+2)+26*vm(m,i+1)+ 6*vm(m,i  ))+&
        omm(0)*(-6*vm(m,i+4)+26*vm(m,i+3)-46*vm(m,i+2)+50*vm(m,i+1))
!
      enddo ! end of m-loop
!
      vminus = vminus/24._R8P
      vplus  = vplus /24._R8P
!
   else
      write(*,*) 'Error! WENO scheme not implemented'
      stop
   endif

   endsubroutine weno_reconstruction_kernel
endmodule adam_nasto_nvf_kernels
