!< ADAM, NASTO CPU Compressible-Navier-Stokes fluidyanmics application library.
module adam_nasto_cpu_cns
!< ADAM, NASTO CPU Compressible-Navier-Stokes fluidyanmics application library.

use penf, only : I4P, R8P

implicit none
private
public :: compute_conservatives
public :: compute_conservative_fluxes
public :: compute_max_eigenvalues
public :: compute_eigenvectors
public :: compute_q_aux

contains
   ! public procedures
   pure subroutine compute_conservatives(b,i,j,k,ngc,q_aux,q)
   !< Compute convervative variables from auxiliary ones.
   integer(I4P), intent(in)    :: b, i, j, k                        !< Cell indexes.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: q(1:)                             !< Conservative varibales.

   q(1) =      q_aux(1,i,j,k,b)
   q(2) = q(1)*q_aux(2,i,j,k,b)
   q(3) = q(1)*q_aux(3,i,j,k,b)
   q(4) = q(1)*q_aux(4,i,j,k,b)
   q(5) = q(1)*q_aux(8,i,j,k,b) - q_aux(7,i,j,k,b)
   endsubroutine compute_conservatives

   pure subroutine compute_conservative_fluxes(sir,b,i,j,k,ngc,q_aux,f)
   !< Compute convervative fluxes from auxiliary variables.
   real(R8P),    intent(in)    :: sir(3)                            !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)    :: b, i, j, k                        !< Cell indexes.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: f(1:)                             !< Conservative fluxes.

   f(1) = q_aux(1,i,j,k,b)*q_aux(2,i,j,k,b)*sir(1) + &
          q_aux(1,i,j,k,b)*q_aux(3,i,j,k,b)*sir(2) + &
          q_aux(1,i,j,k,b)*q_aux(4,i,j,k,b)*sir(3)
   f(2) = f(1)*q_aux(2,i,j,k,b) + q_aux(7,i,j,k,b)*sir(1)
   f(3) = f(1)*q_aux(3,i,j,k,b) + q_aux(7,i,j,k,b)*sir(2)
   f(4) = f(1)*q_aux(4,i,j,k,b) + q_aux(7,i,j,k,b)*sir(3)
   f(5) = f(1)*q_aux(8,i,j,k,b)
   endsubroutine compute_conservative_fluxes

   pure subroutine compute_max_eigenvalues(si,sir,weno_s,b,i,j,k,ngc,nv,q_aux,evmax)
   ! Compute maximum eigenvalues in the big stencil.
   integer(I4P), intent(in)    :: si(3)                             !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)    :: sir(3)                            !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)    :: weno_s                            !< Weno stencils number/dimension.
   integer(I4P), intent(in)    :: b, i, j, k                        !< Cell indexes.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                !< Number of conservative varibales.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: evmax(1:)                         !< Maximum eigenvalues in the big stencil.
   real(R8P)                   :: uu, c                             !< Speeds.
   real(R8P)                   :: ev(nv)                            !< Signals speeds.
   integer(I4P)                :: s, is, js, ks, v                  !< Counter.

   evmax = -1._R8P
   do s=1, 2*weno_s
      is = i + (s-weno_s) * si(1) ; js = j + (s-weno_s) * si(2) ; ks = k + (s-weno_s) * si(3)
      uu = q_aux(1+1*si(1)+2*si(2)+3*si(3),is,js,ks,b)
      c  = q_aux(9                        ,is,js,ks,b)
      ev(1) = abs(uu-c) ; ev(2) = abs(uu) ; ev(3) = abs(uu+c) ; ev(4) = ev(2) ; ev(5) = ev(2)
      do v=1,nv
         evmax(v) = max(ev(v),evmax(v))
      enddo
   enddo
   endsubroutine compute_max_eigenvalues

   pure subroutine compute_eigenvectors(si,sir,b,i,j,k,ngc,nv,g,q_aux,el,er)
   ! Compute eigenvectors centered in inteface i,j,k/ip,jp,kp.
   integer(I4P), intent(in)    :: si(3)                             !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)    :: sir(3)                            !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)    :: b, i, j, k                        !< Cell indexes.
   integer(I4P), intent(in)    :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)    :: nv                                !< Number of conservative varibales.
   real(R8P),    intent(in)    :: g                                 !< Specific heats ratio.
   real(R8P),    intent(in)    :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: el(1:nv,1:nv), er(1:nv,1:nv)      !< Left and right eigenvectors.
   real(R8P)                   :: uu, vv, ww, h, qq, c, ci, b1, b2  !< Roe average states.
   real(R8P)                   :: uvw, uvw_r1, uvw_r2               !< Velocity rotation accordingly dir.

   call compute_roe_average(q_aux=q_aux, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i+si(1), jp=j+si(2), kp=k+si(3), &
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
   endsubroutine compute_eigenvectors

   subroutine compute_q_aux(ni, nj, nk, ngc, blocks_number, R, cv, g, q, q_aux)
   !< Compute auxiliary variables.
   integer(I4P), intent(in)  :: ni                                !< Grid cells number in I direction.
   integer(I4P), intent(in)  :: nj                                !< Grid cells number in J direction.
   integer(I4P), intent(in)  :: nk                                !< Grid cells number in K direction.
   integer(I4P), intent(in)  :: ngc                               !< Ghost cells number.
   integer(I4P), intent(in)  :: blocks_number                     !< Number of blocks.
   real(R8P),    intent(in)  :: R                                 !< Fluid constant, specific heats difference.
   real(R8P),    intent(in)  :: cv                                !< Specific heat at constant volume.
   real(R8P),    intent(in)  :: g                                 !< Specific heats ratio.
   real(R8P),    intent(in)  ::     q(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Conservative variables.
   real(R8P),    intent(out) :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   integer(I4P)              :: b, i, j, k, s                     !< Counter.
   real(R8P)                 :: rho, uuu, vvv, www, rhe, tem      !< State variables.

   !$omp parallel do collapse(4) default(firstprivate) shared(q,q_aux)
   do b=1, blocks_number
      do k=1-ngc, nk+ngc
         do j=1-ngc, nj+ngc
            do i=1-ngc, ni+ngc
               rho = q(1,i,j,k,b)
               uuu = q(2,i,j,k,b)/rho
               vvv = q(3,i,j,k,b)/rho
               www = q(4,i,j,k,b)/rho
               rhe = q(5,i,j,k,b)
               tem = (rhe/rho-0.5*(uuu**2+vvv**2+www**2))/cv

               q_aux(1,i,j,k,b) = rho           ! density
               q_aux(2,i,j,k,b) = uuu           ! velocity x
               q_aux(3,i,j,k,b) = vvv           ! velocity y
               q_aux(4,i,j,k,b) = www           ! velocity z
               q_aux(5,i,j,k,b) = 0._R8P        ! mass fraction
               q_aux(6,i,j,k,b) = tem           ! temperature
               q_aux(7,i,j,k,b) = R*rho*tem     ! pressure
               q_aux(8,i,j,k,b) = rhe/rho+R*tem ! entalpy
               q_aux(9,i,j,k,b) = sqrt(g*R*tem) ! sound speed
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do
   endsubroutine compute_q_aux

   ! private procedures
   pure subroutine compute_roe_average(ngc, b, i, j, k, ip, jp, kp, g, q_aux, &
                                       uu, vv, ww, h, qq, c, ci, b1, b2)
   !< Compute Roe averaged quantities.
   integer(I4P), intent(in)  :: ngc                               !< Number of ghost cells.
   integer(I4P), intent(in)  :: b, i, j, k, ip, jp, kp            !< Left/right cells indexes.
   real(R8P),    intent(in)  :: g                                 !< Specific heats ratio.
   real(R8P),    intent(in)  :: q_aux(1:,1-ngc:,1-ngc:,1-ngc:,1:) !< Auxiliary variables.
   real(R8P),    intent(out) :: uu, vv, ww, h, qq, c, ci, b1, b2  !< Roe state average variables.
   real(R8P)                 :: ri, up, vp, wp, hp, r, rp1, cc    !< Local varbiables.

   ! left state (node i)
   ri = 1._R8P/q_aux(1,i,j,k,b)
   uu = q_aux(2,i,j,k,b)
   vv = q_aux(3,i,j,k,b)
   ww = q_aux(4,i,j,k,b)
   h  = q_aux(8,i,j,k,b)
   ! right state (node i+1)
   up = q_aux(2,ip,jp,kp,b)
   vp = q_aux(3,ip,jp,kp,b)
   wp = q_aux(4,ip,jp,kp,b)
   hp = q_aux(8,ip,jp,kp,b)
   ! Roe average state
   r   = sqrt(q_aux(1,ip,jp,kp,b)*ri)
   rp1 = 1._R8P/(r+1._R8P)
   uu  = (r*up+uu)*rp1
   vv  = (r*vp+vv)*rp1
   ww  = (r*wp+ww)*rp1
   h   = (r*hp+h)*rp1
   qq  = 0.5_R8P * (uu*uu+vv*vv+ww*ww)
   cc  = (g-1._R8P) * (h - qq)
   c   = sqrt(cc)
   ci  = 1._R8P/c
   b2  = (g-1)/cc  ! alias 1/(cp*theta)
   b1  = b2 * qq   ! alias q/(cp*theta)
   endsubroutine compute_roe_average
endmodule adam_nasto_cpu_cns
