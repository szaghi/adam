!< ADAM, NASTO CPU Compressible-Navier-Stokes fluidyanmics application library.
module adam_nasto_cpu_cns
!< ADAM, NASTO CPU Compressible-Navier-Stokes fluidyanmics application library.

use penf, only : I4P, R8P

implicit none
private
public :: compute_conservatives, compute_conservatives_scalar
public :: compute_conservative_fluxes, compute_conservative_fluxes_scalar
public :: compute_max_eigenvalues
public :: compute_eigenvectors
public :: compute_q_aux
public :: compute_riemann_exact
public :: compute_riemann_exact_2
public :: compute_riemann_exact_3
public :: compute_riemann_hllc
public :: compute_riemann_llf
public :: compute_riemann_ts

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
   real(R8P), parameter        :: o_2=0.5_R8P                       !< Alias 1/2.
   real(R8P)                   :: ci1,ci2,ci3                       !< Alias sir(1)/c, sir(2)/c, sir(3)/c.

   call compute_roe_average(q_aux=q_aux, g=g, ngc=ngc, b=b, i=i, j=j, k=k, ip=i+si(1), jp=j+si(2), kp=k+si(3), &
                            uu=uu, vv=vv, ww=ww, h=h, qq=qq, c=c, ci=ci, b1=b1, b2=b2)
   ! cc  = (g-1._R8P) * (h - qq)
   ! c   = sqrt(cc)
   ! ci  = 1._R8P/c
   ! b2  = (g-1)/cc  ! alias 1/(cp*theta)
   ! b1  = b2 * qq   ! alias q/(cp*theta)

   ! uvw    =  uu*sir(1)+vv*sir(2)+ww*sir(3)
   ! uvw_r1 =  uu*sir(3)+vv*sir(1)+ww*sir(2)
   ! uvw_r2 = -uu*sir(2)+vv*sir(3)+ww*sir(1)

   ! er(1,1)=1._R8P ; er(1,2)=uu-c*sir(1)   ; er(1,3)=vv-c*sir(2) ; er(1,4)=ww-c*sir(3)   ; er(1,5)=h-uvw*c
   ! er(2,1)=1._R8P ; er(2,2)=uu            ; er(2,3)=vv          ; er(2,4)=ww            ; er(2,5)=qq
   ! er(3,1)=1._R8P ; er(3,2)=uu+c*sir(1)   ; er(3,3)=vv+c*sir(2) ; er(3,4)=ww+c*sir(3)   ; er(3,5)=h+uvw*c
   ! er(4,1)=0._R8P ; er(4,2)=sir(2)+sir(3) ; er(4,3)=sir(1)      ; er(4,4)=0._R8P        ; er(4,5)=uvw_r1
   ! er(5,1)=0._R8P ; er(5,2)=0._R8P        ; er(5,3)=sir(3)      ; er(5,4)=sir(1)+sir(2) ; er(5,5)=uvw_r2

   ! el(1,1)= 0.5_R8P*(b1+uvw*ci)      ;el(1,2)= 1._R8P-b1;el(1,3)= 0.5_R8P*(b1-uvw*ci)      ;el(1,4)=-uvw_r1;el(1,5)=-uvw_r2
   ! el(2,1)=-0.5_R8P*(b2*uu+ci*sir(1));el(2,2)= b2*uu    ;el(2,3)=-0.5_R8P*(b2*uu-ci*sir(1));el(2,4)= sir(3);el(2,5)=-sir(2)
   ! el(3,1)=-0.5_R8P*(b2*vv+ci*sir(2));el(3,2)= b2*vv    ;el(3,3)=-0.5_R8P*(b2*vv-ci*sir(2));el(3,4)= sir(1);el(3,5)= sir(3)
   ! el(4,1)=-0.5_R8P*(b2*ww+ci*sir(3));el(4,2)= b2*ww    ;el(4,3)=-0.5_R8P*(b2*ww-ci*sir(3));el(4,4)= sir(2);el(4,5)= sir(1)
   ! el(5,1)= 0.5_R8P*b2               ;el(5,2)=-b2       ;el(5,3)= 0.5_R8P*b2               ;el(5,4)= 0._R8P;el(5,5)= 0._R8P

   uvw    = uu*sir(1)+vv*sir(2)+ww*sir(3)
   uvw_r1 = (vv-uvw*sir(2))*sir(1)+(uvw*sir(1)-uu)*sir(2)+(uu-uvw*sir(1))*sir(3)
   uvw_r2 = (uvw*sir(3)-ww)*sir(1)+(ww-uvw*sir(3))*sir(2)+(uvw*sir(2)-vv)*sir(3)
   ci1    = ci*sir(1)
   ci2    = ci*sir(2)
   ci3    = ci*sir(3)

   er(1,1)= 1._R8P      ; er(1,2)= 1._R8P ; er(1,3)= 1._R8P      ; er(1,4)= 0._R8P              ; er(1,5)= 0._R8P
   er(2,1)= uu-c*sir(1) ; er(2,2)= uu     ; er(2,3)= uu+c*sir(1) ; er(2,4)= sir(2)              ; er(2,5)=-sir(3)
   er(3,1)= vv-c*sir(2) ; er(3,2)= vv     ; er(3,3)= vv+c*sir(2) ; er(3,4)=-sir(1)              ; er(3,5)= 0._R8P
   er(4,1)= ww-c*sir(3) ; er(4,2)= ww     ; er(4,3)= ww+c*sir(3) ; er(4,4)= 0._R8P              ; er(4,5)= sir(1)
   er(5,1)= h-uvw*c     ; er(5,2)= qq     ; er(5,3)= h+uvw*c     ; er(5,4)= uu*sir(2)-vv*sir(1) ; er(5,5)= ww*sir(1)-uu*sir(3)

   el(1,1)= o_2*(b1+uvw*ci);el(1,2)=-o_2*(b2*uu+ci1);el(1,3)=-o_2*(b2*vv+ci2);el(1,4)=-o_2*(b2*ww+ci3);el(1,5)= o_2*b2
   el(2,1)= 1._R8P-b1      ;el(2,2)= b2*uu          ;el(2,3)= b2*vv          ;el(2,4)= b2*ww          ;el(2,5)=-b2
   el(3,1)= o_2*(b1-uvw*ci);el(3,2)=-o_2*(b2*uu-ci1);el(3,3)=-o_2*(b2*vv-ci2);el(3,4)=-o_2*(b2*ww-ci3);el(3,5)= o_2*b2
   el(4,1)= uvw_r1         ;el(4,2)= sir(2)-sir(3)  ;el(4,3)=-sir(1)         ;el(4,4)= 0._R8P         ;el(4,5)= 0._R8P
   el(5,1)= uvw_r2         ;el(5,2)= 0._R8P         ;el(5,3)= sir(3)         ;el(5,4)= sir(1)-sir(2)  ;el(5,5)= 0._R8P
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
               q_aux(5,i,j,k,b) = g             ! specific heats ratio
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

   subroutine compute_riemann_exact(si, sir, nv, q_aux1, q_aux4, g1, g4, F, lmax, ws)
   !< Solve the Riemann problem between the state $1$ and $4$ using exact Rainkine-Hugonoit jump relations.
   integer(I4P), intent(in)            :: si(3)                   !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)                  !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: nv                      !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)              !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)              !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: g1                      !< Specific heats ratio of state 1.
   real(R8P),    intent(in)            :: g4                      !< Specific heats ratio of state 4.
   real(R8P),    intent(inout)         :: F(1:)                   !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax                    !< Maximum wave speed estimation.
   real(R8P),    intent(out), optional :: ws(8)                   !< Maximum wave speed estimation.
   real(R8P)                           :: q_auxs(1:9)             !< Intermediate state.
   real(R8P)                           :: S                       !< Velocity of the intermediate states.
   real(R8P)                           :: S1                      !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: S4                      !< Maximum wave speed of state 1 and 4.
   integer(I4P)                        :: uni, ut1, ut2           !< Index of normal and tangential velocities.
   real(R8P)                           :: pF,rF,aF                !< Primitive variables at interface.
   real(R8P)                           :: a2,a3                   !< Speed of sound of state 1, 2, 3 and 4.
   real(R8P)                           :: gm1_1,gp1_1,delta1,eta1 !< g-1, g+1, (g-1)/2, 2*g/(g-1)
   real(R8P)                           :: gm1_4,gp1_4,delta4,eta4 !< g-1, g+1, (g-1)/2, 2*g/(g-1)
   real(R8P)                           :: p2,r2                   !< Primitive variables of state 2.
   real(R8P)                           :: p3,r3                   !< Primitive variables of state 3.
   real(R8P)                           :: dp2,dp3                 !< Derivate of pessure (dp/du) of state 2 and 3.
   real(R8P)                           :: u23, p23                !< Speed and pressure of state 2 and 3.
   real(R8P)                           :: S2                      !< Upstream front of C1 wave.
   real(R8P)                           :: S3                      !< Downstream front of C2 wave.
   real(R8P)                           :: dum,alfa,beta           !< Dummies coefficients.

   uni = 1 + 1*si(1)+2*si(2)+3*si(3)
   ut1 = 1 + findloc(si, 0_I4P             , dim=1)
   ut2 = 1 + findloc(si, 0_I4P, back=.true., dim=1)
   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             v1=>q_aux1(ut1), v4=>q_aux4(ut1), &
             w1=>q_aux1(ut2), w4=>q_aux4(ut2), &
             p1=>q_aux1(7  ), p4=>q_aux4(7  ))
      gm1_1  = g1 - 1.0_R8P  ! g-1 for state 1
      gm1_4  = g4 - 1.0_R8P  ! g-1 for state 4
      gp1_1  = g1 + 1.0_R8P  ! g+1 for state 1
      gp1_4  = g4 + 1.0_R8P  ! g+1 for state 4
      delta1 = 0.5_R8P*gm1_1 ! (g-1)/2 for state 1
      delta4 = 0.5_R8P*gm1_4 ! (g-1)/2 for state 4
      eta1   = g1/delta1     ! 2*g/(g-1) for state 1
      eta4   = g4/delta4     ! 2*g/(g-1) for state 4
      if (p1<p4) then
        dum = 0.5_R8P*gm1_4/g4 ! (g-1)/(g*2)
      else
        dum = 0.5_R8P*gm1_1/g1 ! (g-1)/(g*2)
      endif
      alfa = (p1/p4)**dum
      beta = alfa*delta1/a1+delta4/a4
      ! computing first approximation of intemediate speed u23
      u23 = (alfa-1.0_R8P)/beta + &
            0.5_R8P*(u1+u4)     + &
            0.5_R8P*(u1-u4)*(alfa*delta1/a1-delta4/a4)/beta ! initializing u speed
      ! computing Newton-Rapson iterations
      Newton: do
         ! computing intermediate states using the current approximation of u23
         call compute_interstates_23u(p1=p1, u1=u1, a1=a1, g1=g1, gm1_1=gm1_1, gp1_1=gp1_1, delta1=delta1, eta1=eta1, &
                                      p4=p4, u4=u4, a4=a4, g4=g4, gm1_4=gm1_4, gp1_4=gp1_4, delta4=delta4, eta4=eta4, &
                                      u23=u23,                                                                        &
                                      r2=r2, p2=p2, a2=a2, r3=r3, p3=p3, a3=a3 ,S1=S1, S2=S2, S3=S3, S4=S4)
         ! computing dp/du for the computed intermediates states
         dp2 = -1._R8P*g1*p2/a2
         dp3 =  1._R8P*g4*p3/a3
         ! evaluating the Newton-Rapson convergence
         if (abs(1.0_R8P - (p2/p3))>=1.e-12_R8P) then ! Newton iterative step
            u23  = u23 - ((p2-p3)/(dp2-dp3))          ! updating u value
         else                                         ! Newton iterations have been converged
            exit Newton
         endif
      enddo Newton
      if (present(lmax)) lmax = max(abs(S1),abs(S4))
      S = u23
      p23 = p2
      if (present(ws)) ws = [S1,S2,S,S3,S4,p23,r2,r3]
      if     (S1> 0._R8P              ) then
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux1,f=F)
      elseif (S1<=0._R8P.and.S2>0._R8P) then
         aF = (a1 + u1*delta1)/(1._R8P + delta1)
         pF = p1*(aF/a1)**eta1
         rF = g1*pF/(aF*aF)
         q_auxs(1  ) = rF
         q_auxs(uni) = aF
         q_auxs(ut1) = v1
         q_auxs(ut2) = w1
         q_auxs(7  ) = pF
         q_auxs(8  ) = g1*pF/((g1-1._R8P)*rF)+0.5_R8P*(aF*aF+v1*v1+w1*w1)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      elseif (S2<=0._R8P.and.S >0._R8P) then
         q_auxs(1  ) = r2
         q_auxs(uni) = S
         q_auxs(ut1) = q_aux1(ut1)
         q_auxs(ut2) = q_aux1(ut2)
         q_auxs(7  ) = p23
         q_auxs(8  ) = g1*p23/((g1-1._R8P)*r2)+0.5_R8P*(S*S+v1*v1+w1*w1)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      elseif (S <=0._R8P.and.S3>0._R8P) then
         q_auxs(1  ) = r3
         q_auxs(uni) = S
         q_auxs(ut1) = q_aux4(ut1)
         q_auxs(ut2) = q_aux4(ut2)
         q_auxs(7  ) = p23
         q_auxs(8  ) = g4*p23/((g4-1._R8P)*r3)+0.5_R8P*(S*S+v4*v4+w4*w4)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      elseif (S3<=0._R8P.and.S4>0._R8P) then
         aF = (a4 - u4*delta4)/(1._R8P + delta4)
         pF = p4*(aF/a4)**eta4
         rF = g4*pF/(aF*aF)
         q_auxs(1  ) = rF
         q_auxs(uni) = -aF
         q_auxs(ut1) = v4
         q_auxs(ut2) = w4
         q_auxs(7  ) = pF
         q_auxs(8  ) = g4*pF/((g4-1._R8P)*rF)+0.5_R8P*(aF*aF+v4*v4+w4*w4)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      elseif (S4<=0._R8P              ) then
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux4,f=F)
      endif
   endassociate
   endsubroutine compute_riemann_exact

   subroutine compute_riemann_exact_2(si, sir, nv, q_aux1, q_aux4, g1, g4, F, lmax)
   !< Solve the Riemann problem between the state $1$ and $4$ using exact Rainkine-Hugonoit jump relations.
   integer(I4P), intent(in)            :: si(3)                   !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)                  !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: nv                      !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)              !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)              !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: g1                      !< Specific heats ratio of state 1.
   real(R8P),    intent(in)            :: g4                      !< Specific heats ratio of state 4.
   real(R8P),    intent(inout)         :: F(1:)                   !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax                    !< Maximum wave speed estimation.
   integer(I4P)                        :: uni, ut1, ut2           !< Index of normal and tangential velocities.
   real(R8P)                           :: q_auxs(1:9)             !< Intermediate state.
   real(R8P)                           :: S                       !< Velocity of the intermediate states.
   real(R8P)                           :: S1                      !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: S4                      !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: pF,rF,aF                !< Primitive variables at interface.
   real(R8P)                           :: a2,a3                   !< Speed of sound of state 1, 2, 3 and 4.
   real(R8P)                           :: gm1_1,gp1_1,delta1,eta1 !< g-1, g+1, (g-1)/2, 2*g/(g-1)
   real(R8P)                           :: gm1_4,gp1_4,delta4,eta4 !< g-1, g+1, (g-1)/2, 2*g/(g-1)
   real(R8P)                           :: p2,r2                   !< Primitive variables of state 2.
   real(R8P)                           :: p3,r3                   !< Primitive variables of state 3.
   real(R8P)                           :: p23_0, u23_dp           !< Derivate of pessure (dp/du) of state 2 and 3.
   real(R8P)                           :: u23, p23                !< Speed and pressure of state 2 and 3.
   real(R8P)                           :: S2                      !< Upstream front of C1 wave.
   real(R8P)                           :: S3                      !< Downstream front of C2 wave.
   real(R8P)                           :: aa1,bb1,aa4,bb4         !< Dummies coefficients.

   uni = 1 + 1*si(1)+2*si(2)+3*si(3)
   ut1 = 1 + findloc(si, 0_I4P             , dim=1)
   ut2 = 1 + findloc(si, 0_I4P, back=.true., dim=1)
   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             v1=>q_aux1(ut1), v4=>q_aux4(ut1), &
             w1=>q_aux1(ut2), w4=>q_aux4(ut2), &
             p1=>q_aux1(7  ), p4=>q_aux4(7  ))
      aa1 = 2._R8P/((g1+1._R8P)*r1) ; bb1 = p1*(g1-1._R8P)/(g1+1._R8P)
      aa4 = 2._R8P/((g4+1._R8P)*r4) ; bb4 = p4*(g4-1._R8P)/(g4+1._R8P)

      gm1_1  = g1 - 1.0_R8P  ! g-1 for state 1
      gm1_4  = g4 - 1.0_R8P  ! g-1 for state 4
      gp1_1  = g1 + 1.0_R8P  ! g+1 for state 1
      gp1_4  = g4 + 1.0_R8P  ! g+1 for state 4
      delta1 = 0.5_R8P*gm1_1 ! (g-1)/2 for state 1
      delta4 = 0.5_R8P*gm1_4 ! (g-1)/2 for state 4
      eta1   = g1/delta1     ! 2*g/(g-1) for state 1
      eta4   = g4/delta4     ! 2*g/(g-1) for state 4

      ! p23_0 = 0.5_R8P * (p1+p4)
      p23_0 = 0.5_R8P * (p1+p4) - 0.125_R8P *(u4-u1)*(r1+r4)*(a1+a4)
      p23_0 = (sqrt(aa1/(p23_0+bb1))*p1 + sqrt(aa4/(p23_0+bb4))*p4 + u1 - u4)/(sqrt(aa1/(p23_0+bb1)) + sqrt(aa4/(p23_0+bb4)))
      p23   = p23_0
      Newton: do
         u23    = integral_hugoniot(p=p23,p1=p1,u1=u1,a1=a1,g1=g1,aa1=aa1,bb1=bb1,p4=p4,u4=u4,a4=a4,g4=g4,aa4=aa4,bb4=bb4)
         u23_dp = integral_hugoniot_dp(p=p23,p1=p1,r1=r1,a1=a1,g1=g1,aa1=aa1,bb1=bb1,p4=p4,r4=r4,a4=a4,g4=g4,aa4=aa4,bb4=bb4)
         p23    = p23 - u23/u23_dp
         if (0.5_R8P*abs(p23-p23_0)/(p23+p23_0)>=1.e-8_R8P) then ! Newton iterative step
            p23_0  = p23
         else                                                    ! Newton iterations have been converged
            exit Newton
         endif
      enddo Newton
      u23 = (u1 + u4 + integral_hugoniot_4(p=p23,p4=p4,a4=a4,g4=g4,aa4=aa4,bb4=bb4) - &
                       integral_hugoniot_1(p=p23,p1=p1,a1=a1,g1=g1,aa1=aa1,bb1=bb1)) * 0.5_R8P

      call compute_interstates_23u(p1=p1, u1=u1, a1=a1, g1=g1, gm1_1=gm1_1, gp1_1=gp1_1, delta1=delta1, eta1=eta1, &
                                   p4=p4, u4=u4, a4=a4, g4=g4, gm1_4=gm1_4, gp1_4=gp1_4, delta4=delta4, eta4=eta4, &
                                   u23=u23,                                                                        &
                                   r2=r2, p2=p2, a2=a2, r3=r3, p3=p3, a3=a3 ,S1=S1, S2=S2, S3=S3, S4=S4)
      if (present(lmax)) lmax = max(abs(S1),abs(S4))
      S = u23
      select case(minloc([-S1,S1*S2,S2*S,S*S3,S3*S4,S4],dim=1))
      case(1) ! left supersonic
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux1,f=F)
      case(2) ! left transonic
         aF = (a1 + u1*delta1)/(1._R8P + delta1)
         pF = p1*(aF/a1)**eta1
         rF = g1*pF/(aF*aF)
         q_auxs(1  ) = rF
         q_auxs(uni) = aF
         q_auxs(ut1) = v1
         q_auxs(ut2) = w1
         q_auxs(7  ) = pF
         q_auxs(8  ) = g1*pF/((g1-1._R8P)*rF)+0.5_R8P*(aF*aF+v1*v1+w1*w1)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      case(3) ! left subsonic
         q_auxs(1  ) = r2
         q_auxs(uni) = S
         q_auxs(ut1) = q_aux1(ut1)
         q_auxs(ut2) = q_aux1(ut2)
         q_auxs(7  ) = p23
         q_auxs(8  ) = g1*p23/((g1-1._R8P)*r2)+0.5_R8P*(S*S+v1*v1+w1*w1)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      case(4) ! right subsonic
         q_auxs(1  ) = r3
         q_auxs(uni) = S
         q_auxs(ut1) = q_aux4(ut1)
         q_auxs(ut2) = q_aux4(ut2)
         q_auxs(7  ) = p23
         q_auxs(8  ) = g4*p23/((g4-1._R8P)*r3)+0.5_R8P*(S*S+v4*v4+w4*w4)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      case(5) ! right transonic
         aF = (a4 - u4*delta4)/(1._R8P + delta4)
         pF = p4*(aF/a4)**eta4
         rF = g4*pF/(aF*aF)
         q_auxs(1  ) = rF
         q_auxs(uni) = -aF
         q_auxs(ut1) = v4
         q_auxs(ut2) = w4
         q_auxs(7  ) = pF
         q_auxs(8  ) = g4*pF/((g4-1._R8P)*rF)+0.5_R8P*(aF*aF+v4*v4+w4*w4)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      case(6) ! right supersonic
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux4,f=F)
      endselect
   endassociate
   contains
      pure function integral_curve_1(p,p1,a1,g1) result(u)
      real(R8P), intent(in) :: p,p1,a1,g1
      real(R8P)             :: u

      u = 2._R8P*a1/(g1-1._R8P)*((p/p1)**((g1-1._R8P)/(2._R8P*g1))-1._R8P)
      endfunction integral_curve_1

      pure function hugoniot_locus_1(p,p1,aa1,bb1) result(u)
      real(R8P), intent(in) :: p,p1,aa1,bb1
      real(R8P)             :: u

      u = (p-p1)*sqrt(aa1/(p+bb1))
      endfunction hugoniot_locus_1

      pure function integral_hugoniot_1(p,p1,a1,g1,aa1,bb1) result(u)
      real(R8P), intent(in) :: p,p1,a1,g1,aa1,bb1
      real(R8P)             :: u

      if (p>p1) then
         u = hugoniot_locus_1(p=p,p1=p1,aa1=aa1,bb1=bb1)
      else
         u = integral_curve_1(p=p,p1=p1,a1=a1,g1=g1)
      endif
      endfunction integral_hugoniot_1

      pure function integral_curve_1_dp(p,p1,r1,a1,g1) result(u)
      real(R8P), intent(in) :: p,p1,r1,a1,g1
      real(R8P)             :: u

      u = (1._R8P/(r1*a1))*(p/p1)**(-(g1+1._R8P)/(2._R8P*g1))
      endfunction integral_curve_1_dp

      pure function hugoniot_locus_1_dp(p,p1,aa1,bb1) result(u)
      real(R8P), intent(in) :: p,p1,aa1,bb1
      real(R8P)             :: u

      u = sqrt(aa1/(p+bb1))*(1._R8P-0.5_R8P*(p-p1)/(p+bb1))
      endfunction hugoniot_locus_1_dp

      pure function integral_hugoniot_1_dp(p,p1,r1,a1,g1,aa1,bb1) result(u)
      real(R8P), intent(in) :: p,p1,r1,a1,g1,aa1,bb1
      real(R8P)             :: u

      if (p>p1) then
         u = hugoniot_locus_1_dp(p=p,p1=p1,aa1=aa1,bb1=bb1)
      else
         u = integral_curve_1_dp(p=p,p1=p1,r1=r1,a1=a1,g1=g1)
      endif
      endfunction integral_hugoniot_1_dp

      pure function integral_curve_4(p,p4,a4,g4) result(u)
      real(R8P), intent(in) :: p,p4,a4,g4
      real(R8P)             :: u

      u = 2._R8P*a4/(g4-1._R8P)*((p/p4)**((g4-1._R8P)/(2._R8P*g4))-1._R8P)
      endfunction integral_curve_4

      pure function hugoniot_locus_4(p,p4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,p4,aa4,bb4
      real(R8P)             :: u

      u = (p-p4)*sqrt(aa4/(p+bb4))
      endfunction hugoniot_locus_4

      pure function integral_hugoniot_4(p,p4,a4,g4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,p4,a4,g4,aa4,bb4
      real(R8P)             :: u

      if (p>p4) then
         u = hugoniot_locus_4(p=p,p4=p4,aa4=aa4,bb4=bb4)
      else
         u = integral_curve_4(p=p,p4=p4,a4=a4,g4=g4)
      endif
      endfunction integral_hugoniot_4

      pure function integral_curve_4_dp(p,p4,r4,a4,g4) result(u)
      real(R8P), intent(in) :: p,p4,r4,a4,g4
      real(R8P)             :: u

      u = (1._R8P/(r4*a4))*(p/p4)**(-(g4+1._R8P)/(2._R8P*g4))
      endfunction integral_curve_4_dp

      pure function hugoniot_locus_4_dp(p,p4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,p4,aa4,bb4
      real(R8P)             :: u

      u = sqrt(aa4/(p+bb4))*(1._R8P-0.5_R8P*(p-p4)/(p+bb4))
      endfunction hugoniot_locus_4_dp

      pure function integral_hugoniot_4_dp(p,p4,r4,a4,g4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,p4,r4,a4,g4,aa4,bb4
      real(R8P)             :: u

      if (p>p4) then
         u = hugoniot_locus_4_dp(p=p,p4=p4,aa4=aa4,bb4=bb4)
      else
         u = integral_curve_4_dp(p=p,p4=p4,r4=r4,a4=a4,g4=g4)
      endif
      endfunction integral_hugoniot_4_dp

      pure function integral_hugoniot(p,u1,p1,a1,g1,aa1,bb1,u4,p4,a4,g4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,u1,p1,a1,g1,aa1,bb1,u4,p4,a4,g4,aa4,bb4
      real(R8P)             :: u
      real(R8P)             :: u2,u3

      u2 = integral_hugoniot_1(p=p,p1=p1,a1=a1,g1=g1,aa1=aa1,bb1=bb1)
      u3 = integral_hugoniot_4(p=p,p4=p4,a4=a4,g4=g4,aa4=aa4,bb4=bb4)
      u = u4 - u1 + u2 + u3
      endfunction integral_hugoniot

      pure function integral_hugoniot_dp(p,p1,r1,a1,g1,aa1,bb1,p4,r4,a4,g4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,p1,r1,a1,g1,aa1,bb1,p4,r4,a4,g4,aa4,bb4
      real(R8P)             :: u

      u = integral_hugoniot_1_dp(p=p,p1=p1,r1=r1,a1=a1,g1=g1,aa1=aa1,bb1=bb1) + &
          integral_hugoniot_4_dp(p=p,p4=p4,r4=r4,a4=a4,g4=g4,aa4=aa4,bb4=bb4)
      endfunction integral_hugoniot_dp
   endsubroutine compute_riemann_exact_2

   subroutine compute_riemann_exact_3(si, sir, nv, q_aux1, q_aux4, g1, g4, F, lmax)
   !< Solve the Riemann problem between the state $1$ and $4$ using exact Rainkine-Hugonoit jump relations.
   integer(I4P), intent(in)            :: si(3)                   !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)                  !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: nv                      !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)              !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)              !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: g1                      !< Specific heats ratio of state 1.
   real(R8P),    intent(in)            :: g4                      !< Specific heats ratio of state 4.
   real(R8P),    intent(inout)         :: F(1:)                   !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax                    !< Maximum wave speed estimation.
   integer(I4P)                        :: uni, ut1, ut2           !< Index of normal and tangential velocities.
   real(R8P)                           :: q_auxs(1:9)             !< Intermediate state.
   real(R8P)                           :: S                       !< Velocity of the intermediate states.
   real(R8P)                           :: S1                      !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: S4                      !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: pF,rF,aF                !< Primitive variables at interface.
   real(R8P)                           :: a2,a3                   !< Speed of sound of state 1, 2, 3 and 4.
   real(R8P)                           :: gm1_1,gp1_1,delta1,eta1 !< g-1, g+1, (g-1)/2, 2*g/(g-1)
   real(R8P)                           :: gm1_4,gp1_4,delta4,eta4 !< g-1, g+1, (g-1)/2, 2*g/(g-1)
   real(R8P)                           :: p2,r2                   !< Primitive variables of state 2.
   real(R8P)                           :: p3,r3                   !< Primitive variables of state 3.
   real(R8P)                           :: p23_0, u23_dp           !< Derivate of pessure (dp/du) of state 2 and 3.
   real(R8P)                           :: u23, p23                !< Speed and pressure of state 2 and 3.
   real(R8P)                           :: S2                      !< Upstream front of C1 wave.
   real(R8P)                           :: S3                      !< Downstream front of C2 wave.
   real(R8P)                           :: aa1,bb1,aa4,bb4         !< Dummies coefficients.
   real(R8P)                           :: beta1,beta4             !< Dummies coefficients.

   uni = 1 + 1*si(1)+2*si(2)+3*si(3)
   ut1 = 1 + findloc(si, 0_I4P             , dim=1)
   ut2 = 1 + findloc(si, 0_I4P, back=.true., dim=1)
   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             v1=>q_aux1(ut1), v4=>q_aux4(ut1), &
             w1=>q_aux1(ut2), w4=>q_aux4(ut2), &
             p1=>q_aux1(7  ), p4=>q_aux4(7  ))
      aa1 = 2._R8P/((g1+1._R8P)*r1) ; bb1 = p1*(g1-1._R8P)/(g1+1._R8P)
      aa4 = 2._R8P/((g4+1._R8P)*r4) ; bb4 = p4*(g4-1._R8P)/(g4+1._R8P)

      beta1  = (g1+1._R8P)/(g1-1._R8P)
      beta4  = (g4+1._R8P)/(g4-1._R8P)

      gm1_1  = g1 - 1.0_R8P  ! g-1 for state 1
      gm1_4  = g4 - 1.0_R8P  ! g-1 for state 4
      gp1_1  = g1 + 1.0_R8P  ! g+1 for state 1
      gp1_4  = g4 + 1.0_R8P  ! g+1 for state 4
      delta1 = 0.5_R8P*gm1_1 ! (g-1)/2 for state 1
      delta4 = 0.5_R8P*gm1_4 ! (g-1)/2 for state 4
      eta1   = g1/delta1     ! 2*g/(g-1) for state 1
      eta4   = g4/delta4     ! 2*g/(g-1) for state 4

      ! p23_0 = 0.5_R8P * (p1+p4)
      p23_0 = 0.5_R8P * (p1+p4) - 0.125_R8P *(u4-u1)*(r1+r4)*(a1+a4)
      p23_0 = (sqrt(aa1/(p23_0+bb1))*p1 + sqrt(aa4/(p23_0+bb4))*p4 + u1 - u4)/(sqrt(aa1/(p23_0+bb1)) + sqrt(aa4/(p23_0+bb4)))
      p23   = p23_0
      Newton: do
         u23    = integral_hugoniot(p=p23,p1=p1,u1=u1,a1=a1,g1=g1,aa1=aa1,bb1=bb1,p4=p4,u4=u4,a4=a4,g4=g4,aa4=aa4,bb4=bb4)
         u23_dp = integral_hugoniot_dp(p=p23,p1=p1,r1=r1,a1=a1,g1=g1,aa1=aa1,bb1=bb1,p4=p4,r4=r4,a4=a4,g4=g4,aa4=aa4,bb4=bb4)
         p23    = p23 - u23/u23_dp
         if (0.5_R8P*abs(p23-p23_0)/(p23+p23_0)>=1.e-8_R8P) then ! Newton iterative step
            p23_0  = p23
         else                                                    ! Newton iterations have been converged
            exit Newton
         endif
      enddo Newton
      u23 = (u1 + u4 + integral_hugoniot_4(p=p23,p4=p4,a4=a4,g4=g4,aa4=aa4,bb4=bb4) - &
                       integral_hugoniot_1(p=p23,p1=p1,a1=a1,g1=g1,aa1=aa1,bb1=bb1)) * 0.5_R8P
      if (p23<=p1) then
         r2 = (p23/p1)**(1._R8P/g1)*r1
      else
         r2 = ((1._R8P+beta1*p23/p1)/((p23/p1)+beta1))*r1
      endif
      if (p23<=p4) then
         r3 = (p23/p4)**(1._R8P/g4)*r4
      else
         r3 = ((1._R8P+beta4*p23/p4)/((p23/p4)+beta4))*r4
      endif
      if (p23>p1) then
         S1 = (r1*u1 - r2*u23)/(r1 - r2)
         S2 = S1
      else
         S1 = u1 - a1
         S2 = u23 - sqrt(g1*p23/r2)
      endif
      if (p23>p4) then
         S4 = (r4*u4 - r3*u23)/(r4 - r3)
         S3 = S4
      else
         S3 = u23 + sqrt(g4*p23/r3)
         S4 = u4 + a4
      endif
      if (present(lmax)) lmax = max(abs(S1),abs(S4))
      S = u23
      select case(minloc([-S1,S1*S2,S2*S,S*S3,S3*S4,S4],dim=1))
      case(1) ! left supersonic
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux1,f=F)
      case(2) ! left transonic
         aF = (a1 + u1*delta1)/(1._R8P + delta1)
         pF = p1*(aF/a1)**eta1
         rF = g1*pF/(aF*aF)
         q_auxs(1  ) = rF
         q_auxs(uni) = aF
         q_auxs(ut1) = v1
         q_auxs(ut2) = w1
         q_auxs(7  ) = pF
         q_auxs(8  ) = g1*pF/((g1-1._R8P)*rF)+0.5_R8P*(aF*aF+v1*v1+w1*w1)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      case(3) ! left subsonic
         q_auxs(1  ) = r2
         q_auxs(uni) = S
         q_auxs(ut1) = q_aux1(ut1)
         q_auxs(ut2) = q_aux1(ut2)
         q_auxs(7  ) = p23
         q_auxs(8  ) = g1*p23/((g1-1._R8P)*r2)+0.5_R8P*(S*S+v1*v1+w1*w1)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      case(4) ! right subsonic
         q_auxs(1  ) = r3
         q_auxs(uni) = S
         q_auxs(ut1) = q_aux4(ut1)
         q_auxs(ut2) = q_aux4(ut2)
         q_auxs(7  ) = p23
         q_auxs(8  ) = g4*p23/((g4-1._R8P)*r3)+0.5_R8P*(S*S+v4*v4+w4*w4)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      case(5) ! right transonic
         aF = (a4 - u4*delta4)/(1._R8P + delta4)
         pF = p4*(aF/a4)**eta4
         rF = g4*pF/(aF*aF)
         q_auxs(1  ) = rF
         q_auxs(uni) = -aF
         q_auxs(ut1) = v4
         q_auxs(ut2) = w4
         q_auxs(7  ) = pF
         q_auxs(8  ) = g4*pF/((g4-1._R8P)*rF)+0.5_R8P*(aF*aF+v4*v4+w4*w4)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      case(6) ! right supersonic
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux4,f=F)
      endselect
   endassociate
   contains
      pure function integral_curve_1(p,p1,a1,g1) result(u)
      real(R8P), intent(in) :: p,p1,a1,g1
      real(R8P)             :: u

      u = 2._R8P*a1/(g1-1._R8P)*((p/p1)**((g1-1._R8P)/(2._R8P*g1))-1._R8P)
      endfunction integral_curve_1

      pure function hugoniot_locus_1(p,p1,aa1,bb1) result(u)
      real(R8P), intent(in) :: p,p1,aa1,bb1
      real(R8P)             :: u

      u = (p-p1)*sqrt(aa1/(p+bb1))
      endfunction hugoniot_locus_1

      pure function integral_hugoniot_1(p,p1,a1,g1,aa1,bb1) result(u)
      real(R8P), intent(in) :: p,p1,a1,g1,aa1,bb1
      real(R8P)             :: u

      if (p>p1) then
         u = hugoniot_locus_1(p=p,p1=p1,aa1=aa1,bb1=bb1)
      else
         u = integral_curve_1(p=p,p1=p1,a1=a1,g1=g1)
      endif
      endfunction integral_hugoniot_1

      pure function integral_curve_1_dp(p,p1,r1,a1,g1) result(u)
      real(R8P), intent(in) :: p,p1,r1,a1,g1
      real(R8P)             :: u

      u = (1._R8P/(r1*a1))*(p/p1)**(-(g1+1._R8P)/(2._R8P*g1))
      endfunction integral_curve_1_dp

      pure function hugoniot_locus_1_dp(p,p1,aa1,bb1) result(u)
      real(R8P), intent(in) :: p,p1,aa1,bb1
      real(R8P)             :: u

      u = sqrt(aa1/(p+bb1))*(1._R8P-0.5_R8P*(p-p1)/(p+bb1))
      endfunction hugoniot_locus_1_dp

      pure function integral_hugoniot_1_dp(p,p1,r1,a1,g1,aa1,bb1) result(u)
      real(R8P), intent(in) :: p,p1,r1,a1,g1,aa1,bb1
      real(R8P)             :: u

      if (p>p1) then
         u = hugoniot_locus_1_dp(p=p,p1=p1,aa1=aa1,bb1=bb1)
      else
         u = integral_curve_1_dp(p=p,p1=p1,r1=r1,a1=a1,g1=g1)
      endif
      endfunction integral_hugoniot_1_dp

      pure function integral_curve_4(p,p4,a4,g4) result(u)
      real(R8P), intent(in) :: p,p4,a4,g4
      real(R8P)             :: u

      u = 2._R8P*a4/(g4-1._R8P)*((p/p4)**((g4-1._R8P)/(2._R8P*g4))-1._R8P)
      endfunction integral_curve_4

      pure function hugoniot_locus_4(p,p4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,p4,aa4,bb4
      real(R8P)             :: u

      u = (p-p4)*sqrt(aa4/(p+bb4))
      endfunction hugoniot_locus_4

      pure function integral_hugoniot_4(p,p4,a4,g4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,p4,a4,g4,aa4,bb4
      real(R8P)             :: u

      if (p>p4) then
         u = hugoniot_locus_4(p=p,p4=p4,aa4=aa4,bb4=bb4)
      else
         u = integral_curve_4(p=p,p4=p4,a4=a4,g4=g4)
      endif
      endfunction integral_hugoniot_4

      pure function integral_curve_4_dp(p,p4,r4,a4,g4) result(u)
      real(R8P), intent(in) :: p,p4,r4,a4,g4
      real(R8P)             :: u

      u = (1._R8P/(r4*a4))*(p/p4)**(-(g4+1._R8P)/(2._R8P*g4))
      endfunction integral_curve_4_dp

      pure function hugoniot_locus_4_dp(p,p4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,p4,aa4,bb4
      real(R8P)             :: u

      u = sqrt(aa4/(p+bb4))*(1._R8P-0.5_R8P*(p-p4)/(p+bb4))
      endfunction hugoniot_locus_4_dp

      pure function integral_hugoniot_4_dp(p,p4,r4,a4,g4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,p4,r4,a4,g4,aa4,bb4
      real(R8P)             :: u

      if (p>p4) then
         u = hugoniot_locus_4_dp(p=p,p4=p4,aa4=aa4,bb4=bb4)
      else
         u = integral_curve_4_dp(p=p,p4=p4,r4=r4,a4=a4,g4=g4)
      endif
      endfunction integral_hugoniot_4_dp

      pure function integral_hugoniot(p,u1,p1,a1,g1,aa1,bb1,u4,p4,a4,g4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,u1,p1,a1,g1,aa1,bb1,u4,p4,a4,g4,aa4,bb4
      real(R8P)             :: u
      real(R8P)             :: u2,u3

      u2 = integral_hugoniot_1(p=p,p1=p1,a1=a1,g1=g1,aa1=aa1,bb1=bb1)
      u3 = integral_hugoniot_4(p=p,p4=p4,a4=a4,g4=g4,aa4=aa4,bb4=bb4)
      u = u4 - u1 + u2 + u3
      endfunction integral_hugoniot

      pure function integral_hugoniot_dp(p,p1,r1,a1,g1,aa1,bb1,p4,r4,a4,g4,aa4,bb4) result(u)
      real(R8P), intent(in) :: p,p1,r1,a1,g1,aa1,bb1,p4,r4,a4,g4,aa4,bb4
      real(R8P)             :: u

      u = integral_hugoniot_1_dp(p=p,p1=p1,r1=r1,a1=a1,g1=g1,aa1=aa1,bb1=bb1) + &
          integral_hugoniot_4_dp(p=p,p4=p4,r4=r4,a4=a4,g4=g4,aa4=aa4,bb4=bb4)
      endfunction integral_hugoniot_dp
   endsubroutine compute_riemann_exact_3

   subroutine compute_riemann_hllc(si, sir, nv, q_aux1, q_aux4, g1, g4, F, lmax)
   !< Solve the Riemann problem between the state $1$ and $4$ using the HLLC solver.
   integer(I4P), intent(in)            :: si(3)          !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)         !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: nv             !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)     !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)     !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: g1             !< Specific heats ratio of state 1.
   real(R8P),    intent(in)            :: g4             !< Specific heats ratio of state 4.
   real(R8P),    intent(inout)         :: F(1:)          !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax           !< Maximum wave speed estimation.
   real(R8P)                           :: Q1(1:nv)       !< State 1.
   real(R8P)                           :: Q4(1:nv)       !< State 4.
   real(R8P)                           :: S              !< Velocity of the intermediate states.
   real(R8P)                           :: S1             !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: S4             !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: Q1S,Q2S,Q3S    !< Intermediate state.
   real(R8P)                           :: E1,E4          !< Energy of state 1 and 4.
   integer(I4P)                        :: uni, ut1, ut2  !< Index of normal and tangential velocities.

   uni = 1 + 1*si(1)+2*si(2)+3*si(3)
   ut1 = 1 + findloc(si, 0_I4P             , dim=1)
   ut2 = 1 + findloc(si, 0_I4P, back=.true., dim=1)
   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             p1=>q_aux1(7  ), p4=>q_aux4(7  ))
   call evaluate_waves_pvrs(uni=uni, q_aux1=q_aux1, q_aux4=q_aux4, g1=g1, g4=g4, S1=S1, S=S, S4=S4, lmax=lmax)
   E1 = q_aux1(8) - q_aux1(7)/q_aux1(1)
   E4 = q_aux4(8) - q_aux4(7)/q_aux4(1)
   S = (r4*u4*(S4-u4) - r1*u1*(S1-u1) + p1 - p4)/(r4*(S4-u4) - r1*(S1-u1))
   select case(minloc([-S1,S1*S,S*S4,S4],dim=1))
   case(1)
      call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux1,f=F)
   case(2)
      call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux1,f=F)
      call compute_conservatives_scalar(q_aux=q_aux1,q=Q1)
      Q1S = r1*(S1-u1)/(S1-S)
      Q2S = Q1S*S
      Q3S = Q1S*(E1+(S-u1)*(S+p1/(r1*(S1-u1))))
      F(1  ) = F(1  ) + S1*(Q1S - Q1(1  ))
      F(uni) = F(uni) + S1*(Q2S - Q1(uni))
      F(5  ) = F(5  ) + S1*(Q3S - Q1(5  ))
      if (F(1)>0._R8P) then
         F(ut1) = F(1)*q_aux1(ut1)
         F(ut2) = F(1)*q_aux1(ut2)
         F(5  ) = F(5) + F(1) * (q_aux1(ut1)*q_aux1(ut1) + q_aux1(ut2)*q_aux1(ut2))/0.5_R8P
      else
         F(ut1) = F(1)*q_aux4(ut1)
         F(ut2) = F(1)*q_aux4(ut2)
         F(5  ) = F(5) + F(1) * (q_aux4(ut1)*q_aux4(ut1) + q_aux4(ut2)*q_aux4(ut2))/0.5_R8P
      endif
   case(3)
      call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux4,f=F)
      call compute_conservatives_scalar(q_aux=q_aux4,q=Q4)
      Q1S = r4*(S4-u4)/(S4-S)
      Q2S = Q1S*S
      Q3S = Q1S*(E4+(S-u4)*(S+p4/(r4*(S4-u4))))
      F(1  ) = F(1  ) + S4*(Q1S - Q4(1  ))
      F(uni) = F(uni) + S4*(Q2S - Q4(uni))
      F(5  ) = F(5  ) + S4*(Q3S - Q4(5  ))
      if (F(1)>0._R8P) then
         F(ut1) = F(1)*q_aux1(ut1)
         F(ut2) = F(1)*q_aux1(ut2)
         F(5  ) = F(5) + F(1) * (q_aux1(ut1)*q_aux1(ut1) + q_aux1(ut2)*q_aux1(ut2))/0.5_R8P
      else
         F(ut1) = F(1)*q_aux4(ut1)
         F(ut2) = F(1)*q_aux4(ut2)
         F(5  ) = F(5) + F(1) * (q_aux4(ut1)*q_aux4(ut1) + q_aux4(ut2)*q_aux4(ut2))/0.5_R8P
      endif
   case(4)
      call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux4,f=F)
   endselect
   endassociate
   endsubroutine compute_riemann_hllc

   subroutine compute_riemann_llf(si, sir, nv, q_aux1, q_aux4, g1, g4, F, lmax)
   !< Solve the Riemann problem between the state $1$ and $4$ using the (local) Lax Friedrichs (Rusanov) solver.
   integer(I4P), intent(in)            :: si(3)         !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)        !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: nv            !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)    !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)    !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: g1            !< Specific heats ratio of state 1.
   real(R8P),    intent(in)            :: g4            !< Specific heats ratio of state 4.
   real(R8P),    intent(inout)         :: F(1:)         !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax          !< Maximum wave speed estimation.
   real(R8P)                           :: lmax_         !< Maximum wave speed estimation, local variable.
   real(R8P)                           :: Q1(1:nv)      !< State 1.
   real(R8P)                           :: Q4(1:nv)      !< State 4.
   real(R8P)                           :: F1(1:nv)      !< State 1 fluxes.
   real(R8P)                           :: F4(1:nv)      !< State 4 fluxes.
   real(R8P)                           :: S             !< Velocity of the intermediate states.
   real(R8P)                           :: S1            !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: S4            !< Maximum wave speed of state 1 and 4.
   integer(I4P)                        :: uni, ut1, ut2 !< Index of normal and tangential velocities.

   uni = 1 + 1*si(1)+2*si(2)+3*si(3)
   ut1 = 1 + findloc(si, 0_I4P             , dim=1)
   ut2 = 1 + findloc(si, 0_I4P, back=.true., dim=1)
   call evaluate_waves_pvrs(uni=uni, q_aux1=q_aux1, q_aux4=q_aux4, g1=g1, g4=g4, S1=S1, S=S, S4=S4, lmax=lmax_)
   if (present(lmax)) lmax = lmax_
   call compute_conservatives_scalar(q_aux=q_aux1,q=Q1)
   call compute_conservatives_scalar(q_aux=q_aux4,q=Q4)
   call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux1,f=F1)
   call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux4,f=F4)
   F(1  ) = 0.5_R8P*(F1(1  ) + F4(1  ) - lmax_*(Q4(1  )-Q1(1  )))
   F(uni) = 0.5_R8P*(F1(uni) + F4(uni) - lmax_*(Q4(uni)-Q1(uni)))
   F(5  ) = 0.5_R8P*(F1(5  ) + F4(5  ) - lmax_*(Q4(5  )-Q1(5  )))
   if (F(1)>0._R8P) then
      F(ut1) = F(1)*q_aux1(ut1)
      F(ut2) = F(1)*q_aux1(ut2)
      F(5  ) = F(5) + F(1) * (q_aux1(ut1)*q_aux1(ut1) + q_aux1(ut2)*q_aux1(ut2))/0.5_R8P
   else
      F(ut1) = F(1)*q_aux4(ut1)
      F(ut2) = F(1)*q_aux4(ut2)
      F(5  ) = F(5) + F(1) * (q_aux4(ut1)*q_aux4(ut1) + q_aux4(ut2)*q_aux4(ut2))/0.5_R8P
   endif
   endsubroutine compute_riemann_llf

  subroutine compute_riemann_ts(si, sir, nv, q_aux1, q_aux4, g1, g4, F, lmax)
   !< Solve the Riemann problem between the state $1$ and $4$ using the TS (two shocks) solver.
   integer(I4P), intent(in)            :: si(3)          !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)         !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: nv             !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)     !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)     !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: g1             !< Specific heats ratio of state 1.
   real(R8P),    intent(in)            :: g4             !< Specific heats ratio of state 4.
   real(R8P),    intent(inout)         :: F(1:)          !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax           !< Maximum wave speed estimation.
   real(R8P)                           :: Q1(1:nv)       !< State 1.
   real(R8P)                           :: Q4(1:nv)       !< State 4.
   real(R8P)                           :: S              !< Velocity of the intermediate states.
   real(R8P)                           :: S1,S2          !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: S4,S3          !< Maximum wave speed of state 1 and 4.
   integer(I4P)                        :: uni, ut1, ut2  !< Index of normal and tangential velocities.
   real(R8P)                           :: q_auxs(1:9)    !< Intermediate state.
   real(R8P)                           :: r2,r3,p23      !< Velocity of the intermediate states.

   uni = 1 + 1*si(1)+2*si(2)+3*si(3)
   ut1 = 1 + findloc(si, 0_I4P             , dim=1)
   ut2 = 1 + findloc(si, 0_I4P, back=.true., dim=1)
   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             v1=>q_aux1(ut1), v4=>q_aux4(ut1), &
             w1=>q_aux1(ut2), w4=>q_aux4(ut2), &
             p1=>q_aux1(7  ), p4=>q_aux4(7  ))
      call compute_interstates_23_ts(p1=p1, r1=r1, u1=u1, a1=a1, g1=g1, &
                                     p4=p4, r4=r4, u4=u4, a4=a4, g4=g4, &
                                     p23=p23, u23=S, r2=r2, r3=r3,      &
                                     S1=S1, S2=S2, S3=S3, S4=S4)
      if (present(lmax)) lmax = max(abs(S1),abs(S4))
      select case(minloc([-S1,S1*S,S*S4,S4],dim=1))
      case(1)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux1,f=F)
      case(2)
         q_auxs(1  ) = r2
         q_auxs(uni) = S
         q_auxs(ut1) = v1
         q_auxs(ut2) = w1
         q_auxs(7  ) = p23
         q_auxs(8  ) = g1*p23/((g1-1._R8P)*r2)+0.5_R8P*(S*S+v1*v1+w1*w1)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      case(3)
         q_auxs(1  ) = r3
         q_auxs(uni) = S
         q_auxs(ut1) = v4
         q_auxs(ut2) = w4
         q_auxs(7  ) = p23
         q_auxs(8  ) = g4*p23/((g4-1._R8P)*r3)+0.5_R8P*(S*S+v4*v4+w4*w4)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_auxs,f=F)
      case(4)
         call compute_conservative_fluxes_scalar(sir=sir,q_aux=q_aux4,f=F)
      endselect
   endassociate
   endsubroutine compute_riemann_ts

   ! private procedures
   pure subroutine compute_conservatives_scalar(q_aux,q)
   !< Compute convervative variables from auxiliary ones, scalar input.
   real(R8P),    intent(in)    :: q_aux(1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: q(1:)     !< Conservative varibales.

   q(1) =      q_aux(1)
   q(2) = q(1)*q_aux(2)
   q(3) = q(1)*q_aux(3)
   q(4) = q(1)*q_aux(4)
   q(5) = q(1)*q_aux(8) - q_aux(7)
   endsubroutine compute_conservatives_scalar

   pure subroutine compute_conservative_fluxes_scalar(sir,q_aux,f)
   !< Compute convervative fluxes from auxiliary variables, scalar input.
   real(R8P), intent(in)    :: sir(3)    !< Directional (1=x,2=y,3=z) increment.
   real(R8P), intent(in)    :: q_aux(1:) !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)     !< Conservative fluxes.

   f(1) = q_aux(1)*q_aux(2)*sir(1) + &
          q_aux(1)*q_aux(3)*sir(2) + &
          q_aux(1)*q_aux(4)*sir(3)
   f(2) = f(1)*q_aux(2) + q_aux(7)*sir(1)
   f(3) = f(1)*q_aux(3) + q_aux(7)*sir(2)
   f(4) = f(1)*q_aux(4) + q_aux(7)*sir(3)
   f(5) = f(1)*q_aux(8)
   endsubroutine compute_conservative_fluxes_scalar

   elemental subroutine compute_interstates_23_ts(p1,r1,u1,a1,g1,p4,r4,u4,a4,g4,u23,p23,r2,r3,S1,S2,S3,S4)
   ! Subroutine for evaluating intermediate states using TS (Two Shocks) approximation.
   real(R8P), intent(in)  :: p1      ! Pressure of state 1.
   real(R8P), intent(in)  :: r1      ! Density of state 1.
   real(R8P), intent(in)  :: u1      ! Velocity of state 1.
   real(R8P), intent(in)  :: a1      ! Speed of sound of state 1.
   real(R8P), intent(in)  :: g1      ! Specific heats ratio of state 1.
   real(R8P), intent(in)  :: p4      ! Pressure of state 4.
   real(R8P), intent(in)  :: r4      ! Density of state 4.
   real(R8P), intent(in)  :: u4      ! Velocity of state 4.
   real(R8P), intent(in)  :: a4      ! Speed of sound of state 1.
   real(R8P), intent(in)  :: g4      ! Specific heats ratio of state 4.
   real(R8P), intent(out) :: u23     ! Velocity of intermediate states.
   real(R8P), intent(out) :: p23     ! Pressure of intermediate states.
   real(R8P), intent(out) :: r2      ! Density of state 2.
   real(R8P), intent(out) :: r3      ! Density of state 3.
   real(R8P), intent(out) :: S1,S2   ! Left signal velocities.
   real(R8P), intent(out) :: S3,S4   ! Right signal velocities.
   real(R8P)              :: gp1,gm1 ! g+1, g-1.
   real(R8P)              :: p_p     ! Pressure ratio.
   real(R8P)              :: g1p,g4p ! Dummy variables for computing TS pressure of intermediate states.
   real(R8P)              :: gm1_gp1 ! g+1, (g-1)/(g+1).

   ! evaluating intermediate states pressure and velocity
   p23 = 0.5_R8P*((p1 + p4) - 0.25_R8P*(u4 - u1)*(r1 + r4)*(a1 + a4))
   gp1 = g1 + 1._R8P ; gm1_gp1 = (g1 - 1._R8P)/gp1 ; g1p = sqrt((2._R8P/(gp1*r1))/(p23 + gm1_gp1*p1))
   gp1 = g4 + 1._R8P ; gm1_gp1 = (g4 - 1._R8P)/gp1 ; g4p = sqrt((2._R8P/(gp1*r4))/(p23 + gm1_gp1*p4))
   p23 = (g1p*p1 + g4p*p4 + u1 - u4)/(g1p + g4p)
   u23 = 0.5_R8P*(u1 + u4 + (p23 - p4)*g4p - (p23 - p1)*g1p)
   ! evaluating left state
   p_p = p23/p1
   gp1 = g1 + 1._R8P
   gm1 = g1 - 1._R8P
   r2  = r1*((gm1 + gp1*p_p)/(gp1 + gm1*p_p))
   S1  = u1 - a1*sqrt(1._R8P + 0.5_R8P*gp1/g1*(p_p-1._R8P))
   S2  = S1
   ! evaluating right state
   p_p = p23/p4
   gp1 = g4 + 1._R8P
   gm1 = g4 - 1._R8P
   r3  = r4*((gm1 + gp1*p_p)/(gp1 + gm1*p_p))
   S4  = u4 + a4*sqrt(1._R8P + 0.5_R8P*gp1/g4*(p_p-1._R8P))
   S3  = S4
   endsubroutine compute_interstates_23_ts

   elemental subroutine compute_interstates_23u(p1,u1,a1,g1,gm1_1,gp1_1,delta1,eta1, &
                                                p4,u4,a4,g4,gm1_4,gp1_4,delta4,eta4, &
                                                u23,                                 &
                                                r2,p2,a2,r3,p3,a3,S1,S2,S3,S4)
   !< Compute intermediates states knowing the value of speed (u23) of intermediates states.
   real(R8P), intent(in)  :: p1          !< Pressure of state 1.
   real(R8P), intent(in)  :: u1          !< Velocity of state 1.
   real(R8P), intent(in)  :: a1          !< Speed of sound of state 1.
   real(R8P), intent(in)  :: g1          !< Specific heats ratio of state 1.
   real(R8P), intent(in)  :: gm1_1,gp1_1 !< g-1, g+1 of state 1.
   real(R8P), intent(in)  :: delta1,eta1 !< (g-1)/2, 2*g/(g-1) of state 1.
   real(R8P), intent(in)  :: p4          !< Pressure of state 4.
   real(R8P), intent(in)  :: u4          !< Velocity of state 4.
   real(R8P), intent(in)  :: a4          !< Speed of sound of state 4.
   real(R8P), intent(in)  :: g4          !< Specific heats ratio of state 4.
   real(R8P), intent(in)  :: gm1_4,gp1_4 !< g-1, g+1 of state 4.
   real(R8P), intent(in)  :: delta4,eta4 !< (g-1)/2, 2*g/(g-1) of state 4.
   real(R8P), intent(in)  :: u23         !< Velocity of intermediate states.
   real(R8P), intent(out) :: r2,p2,a2    !< Density, pressure and speed of sound of state 2.
   real(R8P), intent(out) :: r3,p3,a3    !< Density, pressure and speed of sound of state 3.
   real(R8P), intent(out) :: S1,S2       !< Left signal velocities.
   real(R8P), intent(out) :: S3,S4       !< Right signal velocities.

   ! computing left wave
   if (abs(u23-u1)<=1.e-16_R8P) then
      call compute_rarefaction(sgn   = -1._R8P, &
                               g     = g1,      &
                               delta = delta1,  &
                               eta   = eta1,    &
                               u0    = u1,      &
                               p0    = p1,      &
                               a0    = a1,      &
                               ux    = u23,     &
                               rx    = r2,      &
                               px    = p2,      &
                               ax    = a2,      &
                               s0    = S1,      &
                               sx    = S2)
   else
      if (u23<u1) then
         call compute_shock(sgn = -1._R8P, &
                            g   = g1,      &
                            gm1 = gm1_1,   &
                            gp1 = gp1_1,   &
                            u0  = u1,      &
                            p0  = p1,      &
                            a0  = a1,      &
                            ux  = u23,     &
                            rx  = r2,      &
                            px  = p2,      &
                            ax  = a2,      &
                            ss  = S1)
         S2 = S1
      else
         call compute_rarefaction(sgn   = -1._R8P, &
                                  g     = g1,      &
                                  delta = delta1,  &
                                  eta   = eta1,    &
                                  u0    = u1,      &
                                  p0    = p1,      &
                                  a0    = a1,      &
                                  ux    = u23,     &
                                  rx    = r2,      &
                                  px    = p2,      &
                                  ax    = a2,      &
                                  s0    = S1,      &
                                  sx    = S2)
      endif
   endif
   ! computing right wave
   if (abs(u23-u4)<=1.e-16_R8P) then
      call compute_rarefaction(sgn   = 1._R8P, &
                               g     = g4,     &
                               delta = delta4, &
                               eta   = eta4,   &
                               u0    = u4,     &
                               p0    = p4,     &
                               a0    = a4,     &
                               ux    = u23,    &
                               rx    = r3,     &
                               px    = p3,     &
                               ax    = a3,     &
                               s0    = S4,     &
                               sx    = S3)
   else
      if (u23>u4) then
         call compute_shock(sgn = 1._R8P, &
                            g   = g4,     &
                            gm1 = gm1_4,  &
                            gp1 = gp1_4,  &
                            u0  = u4,     &
                            p0  = p4,     &
                            a0  = a4,     &
                            ux  = u23,    &
                            rx  = r3,     &
                            px  = p3,     &
                            ax  = a3,     &
                            ss  = S4)
         S3 = S4
      else
         call compute_rarefaction(sgn   = 1._R8P, &
                                  g     = g4,     &
                                  delta = delta4, &
                                  eta   = eta4,   &
                                  u0    = u4,     &
                                  p0    = p4,     &
                                  a0    = a4,     &
                                  ux    = u23,    &
                                  rx    = r3,     &
                                  px    = p3,     &
                                  ax    = a3,     &
                                  s0    = S4,     &
                                  sx    = S3)
      endif
   endif
   endsubroutine compute_interstates_23u

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

   elemental subroutine compute_rarefaction(sgn,g,delta,eta,u0,p0,a0,ux,rx,px,ax,s0,sx)
   !< Compute an unknown state "x" from a known state "0" when the two states are separated by a rarefaction; it is
   !< assumed that the velocity of the unknown state "ux" is known. There is an input variable that indicates if the rarefaction
   !< propagates on the "u-a" direction (left) or on the "u+a" one (right).
   real(R8P), intent(in)  :: sgn      !< Sign for distinguishing "left" (-1) from "right" (1) wave.
   real(R8P), intent(in)  :: g        !< Specific heats ratio.
   real(R8P), intent(in)  :: delta    !< (g-1)/2.
   real(R8P), intent(in)  :: eta      !< 2*g/(g-1).
   real(R8P), intent(in)  :: u0,p0,a0 !< Known state (speed, pressure and speed of sound).
   real(R8P), intent(in)  :: ux       !< Known speed of unknown state.
   real(R8P), intent(out) :: rx,px,ax !< Unknown pressure and density.
   real(R8P), intent(out) :: s0,sx    !< Wave speeds (head and back fronts).

   ax = a0 + sgn*delta*(ux - u0) ! unknown speed of sound
   px = p0*((ax/a0)**(eta))      ! unknown pressure
   rx = g*px/(ax*ax)             ! unknown density
   s0 = u0 + sgn*a0              ! left wave speed
   sx = ux + sgn*ax              ! right wave speed
   endsubroutine compute_rarefaction

   elemental subroutine compute_shock(sgn,g,gm1,gp1,u0,p0,a0,ux,rx,px,ax,ss)
   !< Compute an unknown state "x" from a known state "0" when the two states are separated by a shock; it is
   !< assumed that the velocity of the unknown state "ux" is known. There is an input variable that indicates if the shock
   !< propagates on the "u-a" direction (left) or on the "u+a" one (right).
   real(R8P), intent(in)  :: sgn      !< Sign for distinguishing "left" (-1) from "right" (1) wave.
   real(R8P), intent(in)  :: g        !< Specific heats ratio.
   real(R8P), intent(in)  :: gm1,gp1  !< Gm1 = g - 1 , gp1 = g + 1.
   real(R8P), intent(in)  :: u0,p0,a0 !< Known state (speed, pressure and speed of sound).
   real(R8P), intent(in)  :: ux       !< Unknown speed.
   real(R8P), intent(out) :: rx,px,ax !< Unknown state (density, pressure and speed of sound).
   real(R8P), intent(out) :: ss       !< Shock wave speed.
   real(R8P)              :: M0       !< Relative Mach number of known state.
   real(R8P)              :: x        !< Dummy variable.

   x   = 0.25_R8P*gp1*(ux - u0)/a0              ! dummy variable
   M0  = x + sgn*sqrt(1.0_R8P + x*x)            ! relative Mach number of known state
   x   = 1._R8P + 2._R8P*g*(M0*M0 - 1._R8P)/gp1 ! dummy variable (pressure ratio px/p0)

   ax = a0*sqrt((gp1 + gm1*x)/(gp1 + gm1/x)) ! unknown speed of sound
   px = p0*x                                 ! unknown pressure
   rx = g*px/(ax*ax)                         ! unknown density
   ss = u0 + a0*M0                           ! shock wave speed
   endsubroutine compute_shock

   pure subroutine evaluate_waves_pvrs(uni, q_aux1, q_aux4, g1, g4, S1, S, S4, lmax, p)
   !< Evaluate intermediate waves 2 and 3 from the known states 1,4 using the PVRS approximation.
   integer(I4P), intent(in)            :: uni        !< Index of normal velocity.
   real(R8P),    intent(in)            :: q_aux1(1:) !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:) !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: g1         !< Specific heats ratio of state 1.
   real(R8P),    intent(in)            :: g4         !< Specific heats ratio of state 4.
   real(R8P),    intent(out)           :: S1, S, S4  !< Waves speed estimations.
   real(R8P),    intent(out), optional :: lmax       !< Maximum wave speed estimation.
   real(R8P),    intent(out), optional :: p          !< Pressure of the intermediate states.
   real(R8P)                           :: p_         !< Pressure of the intermediate states, local variable.
   real(R8P)                           :: ram        !< Mean value of rho*a.

   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             p1=>q_aux1(7  ), p4=>q_aux4(7  ))
   ram = 0.5_R8P * (r1 + r4) * 0.5_R8P * (a1 + a4)       ! product of mean density for mean speed of sound
   S   = 0.5_R8P * (u1 + u4) - 0.5_R8P * (p4 - p1) / ram ! evaluation of the contact wave speed (velocity of intermediate states)
   p_  = 0.5_R8P * (p1 + p4) - 0.5_R8P * (u4 - u1) * ram ! evaluation of the pressure of the intermediate states
   if (present(p)) p = p_
   ! evaluation of the left wave speeds
   if (p_<=p1*(1._R8P + 1.e-12_R8P)) then ! rarefaction
     S1 = u1 - a1
   else ! shock
     S1 = u1 - a1 * sqrt(1._R8P + (g1 + 1._R8P) / (2._R8P * g1) * (p_ / p1 - 1._R8P))
   endif
   ! evaluation of the right wave speeds
   if (p_<=p4 * (1._R8P + 1.e-12_R8P)) then ! rarefaction
     S4 = u4 + a4
   else ! shock
     S4 = u4 + a4 * sqrt(1._R8P + (g4 + 1._R8P) / (2._R8P * g4) * ( p_ / p4 - 1._R8P))
   endif
   if (present(lmax)) lmax = max(abs(S1), abs(S), abs(S4))
   endassociate
   endsubroutine evaluate_waves_pvrs
endmodule adam_nasto_cpu_cns
