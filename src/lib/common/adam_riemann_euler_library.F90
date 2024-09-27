!< ADAM, Riemann Problem for Euler equations solvers and convective fluxes computations library.
module adam_riemann_euler_library
!< ADAM, Riemann Problem for Euler equations solvers and convective fluxes computations library.

use penf, only : I4P, R8P

implicit none
private
public :: compute_riemann_euler_exact
public :: compute_riemann_euler_hllc
public :: compute_riemann_euler_hllc_lm
public :: compute_riemann_euler_hllem
public :: compute_riemann_euler_llf
public :: compute_riemann_euler_ts

contains
   ! public procedures
   subroutine compute_riemann_euler_exact(si, sir, uni, ut1, ut2, nv, q_aux1, q_aux4, f, lmax)
   !< Solve the Riemann problem between the state 1 (left) and 4 (right) using exact Integral/Rainkine-Hugonoit relations.
   integer(I4P), intent(in)            :: si(3)                     !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)                    !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: uni, ut1, ut2             !< Index of normal and tangential velocities.
   integer(I4P), intent(in)            :: nv                        !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)                !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)                !< Left state, auxiliary variables.
   real(R8P),    intent(inout)         :: f(1:)                     !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax                      !< Maximum wave speed estimation.
   real(R8P)                           :: u1,r1,g1,p1,a1            !< Left state.
   real(R8P)                           :: u4,r4,g4,p4,a4            !< Right state.
   real(R8P)                           :: q_auxs(1:9)               !< Intermediate state conservative variables.
   real(R8P)                           :: S,S1,S2,S3,S4             !< Waves speeds.
   real(R8P)                           :: u23,p23,p2,r2,a2,p3,r3,a3 !< Intermediate states.
   real(R8P)                           :: pF,rF,aF                  !< Transonic state.
   real(R8P)                           :: p23_0, u23_dp             !< Netwon-Rapson variables.
   real(R8P)                           :: delta1,eta1               !< (g-1)/2, 2*g/(g-1)
   real(R8P)                           :: delta4,eta4               !< (g-1)/2, 2*g/(g-1)
   real(R8P)                           :: aa1,bb1,aa4,bb4           !< Dummies coefficients.
   real(R8P)                           :: beta1,beta4               !< Dummies coefficients.

   u1 = q_aux1(uni) ; u4 = q_aux4(uni)
   r1 = q_aux1(1  ) ; r4 = q_aux4(1  )
   g1 = q_aux1(5  ) ; g4 = q_aux4(5  )
   p1 = q_aux1(7  ) ; p4 = q_aux4(7  )
   a1 = q_aux1(9  ) ; a4 = q_aux4(9  )
   associate(v1=>q_aux1(ut1), v4=>q_aux4(ut1), &
             w1=>q_aux1(ut2), w4=>q_aux4(ut2))

      aa1 = 2._R8P/((g1+1._R8P)*r1) ; bb1 = p1*(g1-1._R8P)/(g1+1._R8P)
      aa4 = 2._R8P/((g4+1._R8P)*r4) ; bb4 = p4*(g4-1._R8P)/(g4+1._R8P)

      beta1  = (g1+1._R8P)/(g1-1._R8P)
      beta4  = (g4+1._R8P)/(g4-1._R8P)

      p23_0 = 0.5_R8P * (p1+p4) - 0.125_R8P *(u4-u1)*(r1+r4)*(a1+a4)
      p23_0 = (sqrt(aa1/(p23_0+bb1))*p1 + sqrt(aa4/(p23_0+bb4))*p4 + u1 - u4)/(sqrt(aa1/(p23_0+bb1)) + sqrt(aa4/(p23_0+bb4)))
      p23   = p23_0
      Newton: do
         u23 = integral_hugoniot(p=p23) ; u23_dp = integral_hugoniot_dp(p=p23)
         p23 = p23 - u23/u23_dp
         if (0.5_R8P*abs(p23-p23_0)/(p23+p23_0)>=1.e-12_R8P) then ! new Newton iterative step
            p23_0  = p23
         else                                                     ! Newton iterations have been converged, exit iterations
            exit Newton
         endif
      enddo Newton
      u23 = 0.5_R8P * (u1 + u4 + integral_hugoniot_4(p=p23) - integral_hugoniot_1(p=p23))
      S = u23
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
      if (p23>p1) then ! shock
         S1 = (r1*u1 - r2*u23)/(r1 - r2)
         S2 = S1
      else             ! rarefaction
         S1 = u1 - a1
         S2 = u23 - sqrt(g1*p23/r2)
      endif
      if (p23>p4) then ! shock
         S4 = (r4*u4 - r3*u23)/(r4 - r3)
         S3 = S4
      else             ! rarefaction
         S3 = u23 + sqrt(g4*p23/r3)
         S4 = u4 + a4
      endif
      if (present(lmax)) lmax = max(abs(S1),abs(S4))
      select case(minloc([-S1,S1*S2,S2*S,S*S3,S3*S4,S4],dim=1))
      case(1) ! left supersonic
         call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux1,f=F)
      case(2) ! left transonic
         delta1 = 0.5_R8P*(g1 - 1.0_R8P)
         eta1   = g1/delta1
         aF = (a1 + u1*delta1)/(1._R8P + delta1)
         pF = p1*(aF/a1)**eta1
         rF = g1*pF/(aF*aF)
         q_auxs(1  ) = rF
         q_auxs(uni) = aF
         q_auxs(ut1) = v1
         q_auxs(ut2) = w1
         q_auxs(7  ) = pF
         q_auxs(8  ) = g1*pF/((g1-1._R8P)*rF)+0.5_R8P*(aF*aF+v1*v1+w1*w1)
         call compute_convective_fluxes_euler(sir=sir,q_aux=q_auxs,f=F)
      case(3) ! left subsonic
         q_auxs(1  ) = r2
         q_auxs(uni) = S
         q_auxs(ut1) = q_aux1(ut1)
         q_auxs(ut2) = q_aux1(ut2)
         q_auxs(7  ) = p23
         q_auxs(8  ) = g1*p23/((g1-1._R8P)*r2)+0.5_R8P*(S*S+v1*v1+w1*w1)
         call compute_convective_fluxes_euler(sir=sir,q_aux=q_auxs,f=F)
      case(4) ! right subsonic
         q_auxs(1  ) = r3
         q_auxs(uni) = S
         q_auxs(ut1) = q_aux4(ut1)
         q_auxs(ut2) = q_aux4(ut2)
         q_auxs(7  ) = p23
         q_auxs(8  ) = g4*p23/((g4-1._R8P)*r3)+0.5_R8P*(S*S+v4*v4+w4*w4)
         call compute_convective_fluxes_euler(sir=sir,q_aux=q_auxs,f=F)
      case(5) ! right transonic
         delta4 = 0.5_R8P*(g4 - 1.0_R8P)
         eta4   = g4/delta4
         aF = (a4 - u4*delta4)/(1._R8P + delta4)
         pF = p4*(aF/a4)**eta4
         rF = g4*pF/(aF*aF)
         q_auxs(1  ) = rF
         q_auxs(uni) = -aF
         q_auxs(ut1) = v4
         q_auxs(ut2) = w4
         q_auxs(7  ) = pF
         q_auxs(8  ) = g4*pF/((g4-1._R8P)*rF)+0.5_R8P*(aF*aF+v4*v4+w4*w4)
         call compute_convective_fluxes_euler(sir=sir,q_aux=q_auxs,f=F)
      case(6) ! right supersonic
         call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux4,f=F)
      endselect
   endassociate
   contains
      ! definitions of u(p) and du(p)/dp for the Newton-Rapson iterative method
      pure function integral_hugoniot(p) result(u)
      !< Return u by isentropic relation or Rankine-Hugonoit jump, left/right curve-locus.
      real(R8P), intent(in) :: p !< Pressure in intermediate state.
      real(R8P)             :: u !< Velocity in intermediate state.

      u = u4 - u1 + integral_hugoniot_1(p) + integral_hugoniot_4(p)
      endfunction integral_hugoniot

      pure function integral_hugoniot_dp(p) result(du)
      !< Return du/dp by isentropic relation or Rankine-Hugonoit jump, left/right curve-locus derivative.
      real(R8P), intent(in) :: p  !< Pressure in intermediate state.
      real(R8P)             :: du !< Velocity in intermediate state.

      du = integral_hugoniot_1_dp(p) + integral_hugoniot_4_dp(p)
      endfunction integral_hugoniot_dp

      pure function integral_hugoniot_1(p) result(u)
      !< Return u by isentropic relation or Rankine-Hugonoit jump, left curve-locus.
      real(R8P), intent(in) :: p !< Pressure in right state.
      real(R8P)             :: u !< Velocity in right state.

      if (p>p1) then
         u = hugoniot_locus_1(p)
      else
         u = integral_curve_1(p)
      endif
      endfunction integral_hugoniot_1

      pure function integral_hugoniot_1_dp(p) result(du)
      !< Return du/dp by isentropic relation or Rankine-Hugonoit jump, left curve-locus derivative.
      real(R8P), intent(in) :: p  !< Pressure in right state.
      real(R8P)             :: du !< Velocity derivative in right state.

      if (p>p1) then
         du = hugoniot_locus_1_dp(p)
      else
         du = integral_curve_1_dp(p)
      endif
      endfunction integral_hugoniot_1_dp

      pure function integral_hugoniot_4(p) result(u)
      !< Return u by isentropic relation or Rankine-Hugonoit jump, right curve-locus.
      real(R8P), intent(in) :: p !< Pressure in right state.
      real(R8P)             :: u !< Velocity in right state.

      if (p>p4) then
         u = hugoniot_locus_4(p)
      else
         u = integral_curve_4(p)
      endif
      endfunction integral_hugoniot_4

      pure function integral_hugoniot_4_dp(p) result(du)
      !< Return du/dp by isentropic relation or Rankine-Hugonoit jump, right curve-locus derivative.
      real(R8P), intent(in) :: p  !< Pressure in left state.
      real(R8P)             :: du !< Velocity derivative in left state.

      if (p>p4) then
         du = hugoniot_locus_4_dp(p)
      else
         du = integral_curve_4_dp(p)
      endif
      endfunction integral_hugoniot_4_dp

      pure function integral_curve_1(p) result(u)
      !< Return u by isentropic relation, left integral curve.
      real(R8P), intent(in) :: p !< Pressure in right state.
      real(R8P)             :: u !< Velocity in right state.

      u = 2._R8P*a1/(g1-1._R8P)*((p/p1)**((g1-1._R8P)/(2._R8P*g1))-1._R8P)
      endfunction integral_curve_1

      pure function hugoniot_locus_1(p) result(u)
      !< Return u by Rankine-Hugoniot jump, left Hugonoit jump locus.
      real(R8P), intent(in) :: p !< Pressure in right state.
      real(R8P)             :: u !< Velocity in right state.

      u = (p-p1)*sqrt(aa1/(p+bb1))
      endfunction hugoniot_locus_1

      pure function integral_curve_1_dp(p) result(du)
      !< Return du/dp by isentropic relation, left integral curve derivative.
      real(R8P), intent(in) :: p  !< Pressure in right state.
      real(R8P)             :: du !< Velocity derivative in right state.

      du = (1._R8P/(r1*a1))*(p/p1)**(-(g1+1._R8P)/(2._R8P*g1))
      endfunction integral_curve_1_dp

      pure function hugoniot_locus_1_dp(p) result(du)
      !< Return du/dp by Rankine-Hugoniot jump, left Hugonoit jump locus derivative.
      real(R8P), intent(in) :: p  !< Pressure in right state.
      real(R8P)             :: du !< Velocity derivative in right state.

      du = sqrt(aa1/(p+bb1))*(1._R8P-0.5_R8P*(p-p1)/(p+bb1))
      endfunction hugoniot_locus_1_dp

      pure function integral_curve_4(p) result(u)
      !< Return u by isentropic relation, right integral curve.
      real(R8P), intent(in) :: p !< Pressure in left state.
      real(R8P)             :: u !< Velocity in left state.

      u = 2._R8P*a4/(g4-1._R8P)*((p/p4)**((g4-1._R8P)/(2._R8P*g4))-1._R8P)
      endfunction integral_curve_4

      pure function hugoniot_locus_4(p) result(u)
      !< Return u by Rankine-Hugoniot jump, right Hugonoit jump locus.
      real(R8P), intent(in) :: p !< Pressure in left state.
      real(R8P)             :: u !< Velocity in left state.

      u = (p-p4)*sqrt(aa4/(p+bb4))
      endfunction hugoniot_locus_4

      pure function integral_curve_4_dp(p) result(du)
      !< Return du/dp by isentropic relation, right integral curve derivative.
      real(R8P), intent(in) :: p  !< Pressure in left state.
      real(R8P)             :: du !< Velocity derivative in left state.

      du = (1._R8P/(r4*a4))*(p/p4)**(-(g4+1._R8P)/(2._R8P*g4))
      endfunction integral_curve_4_dp

      pure function hugoniot_locus_4_dp(p) result(du)
      !< Return du/dp by Rankine-Hugoniot jump, right Hugonoit jump locus derivative.
      real(R8P), intent(in) :: p  !< Pressure in left state.
      real(R8P)             :: du !< Velocity derivative in left state.

      du = sqrt(aa4/(p+bb4))*(1._R8P-0.5_R8P*(p-p4)/(p+bb4))
      endfunction hugoniot_locus_4_dp
   endsubroutine compute_riemann_euler_exact

   subroutine compute_riemann_euler_hllc(si, sir, uni, ut1, ut2, nv, q_aux1, q_aux4, f, lmax)
   !< Solve the Riemann problem between the state 1 (left) and 4 (right) using the HLLC solver.
   integer(I4P), intent(in)            :: si(3)          !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)         !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: uni, ut1, ut2  !< Index of normal and tangential velocities.
   integer(I4P), intent(in)            :: nv             !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)     !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)     !< Left state, auxiliary variables.
   real(R8P),    intent(inout)         :: f(1:)          !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax           !< Maximum wave speed estimation.
   real(R8P)                           :: Q1(1:nv)       !< State 1.
   real(R8P)                           :: Q4(1:nv)       !< State 4.
   real(R8P)                           :: Q1S,Q2S,Q3S    !< Intermediate state.
   real(R8P)                           :: S              !< Velocity of the intermediate states.
   real(R8P)                           :: S1             !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: S4             !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: E1,E4          !< Energy of state 1 and 4.
   real(R8P)                           :: u_roe, a_roe   !< Roe averages.
   real(R8P)                           :: sqrt_r1        !< Sqrt(r1).
   real(R8P)                           :: sqrt_r4        !< Sqrt(r4).
   real(R8P)                           :: sqrt_r14       !< Sqrt(r1) + Sqrt(r4).

   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             v1=>q_aux1(ut1), v4=>q_aux4(ut1), &
             w1=>q_aux1(ut2), w4=>q_aux4(ut2), &
             p1=>q_aux1(7  ), p4=>q_aux4(7  ))

   sqrt_r1 = sqrt(r1)
   sqrt_r4 = sqrt(r4)
   sqrt_r14 = sqrt_r1 + sqrt_r4
   u_roe = (u1*sqrt_r1 + u4*sqrt_r4)/sqrt_r14
   a_roe = sqrt((a1**2*sqrt_r1+a4**2*sqrt_r4)/sqrt_r14 + 0.5_R8P*sqrt_r1*sqrt_r4*((u4-u1)/sqrt_r14)**2)
   S1 = min(u1-a1, u_roe-a_roe)
   S4 = max(u4+a4, u_roe+a_roe)
   if (present(lmax)) lmax = max(abs(S1),abs(S4))
   S = (r4*u4*(S4-u4) - r1*u1*(S1-u1) + p1 - p4)/(r4*(S4-u4) - r1*(S1-u1))
   E1 = q_aux1(8) - q_aux1(7)/q_aux1(1)
   E4 = q_aux4(8) - q_aux4(7)/q_aux4(1)
   select case(minloc([-S1,S1*S,S*S4,S4],dim=1))
   case(1)
      call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux1,f=f)
   case(2)
      call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux1,f=f)
      call compute_conservatives_euler(q_aux=q_aux1,q=Q1)
      Q1S = r1*(S1-u1)/(S1-S)
      Q2S = Q1S*S
      Q3S = Q1S*(E1+(S-u1)*(S+p1/(r1*(S1-u1))))
      f(1  ) = f(1  ) + S1*(Q1S - Q1(1  ))
      f(uni) = f(uni) + S1*(Q2S - Q1(uni))
      f(5  ) = f(5  ) + S1*(Q3S - Q1(5  ))
      if (F(1)>0._R8P) then
         f(ut1) = f(1)*q_aux1(ut1)
         f(ut2) = f(1)*q_aux1(ut2)
         f(5  ) = f(5) + f(1) * (v1*v1 + w1*w1)/0.5_R8P
      else
         f(ut1) = f(1)*q_aux4(ut1)
         f(ut2) = f(1)*q_aux4(ut2)
         f(5  ) = f(5) + f(1) * (v4*v4 + w4*w4)/0.5_R8P
      endif
   case(3)
      call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux4,f=f)
      call compute_conservatives_euler(q_aux=q_aux4,q=Q4)
      Q1S = r4*(S4-u4)/(S4-S)
      Q2S = Q1S*S
      Q3S = Q1S*(E4+(S-u4)*(S+p4/(r4*(S4-u4))))
      f(1  ) = f(1  ) + S4*(Q1S - Q4(1  ))
      f(uni) = f(uni) + S4*(Q2S - Q4(uni))
      f(5  ) = f(5  ) + S4*(Q3S - Q4(5  ))
      if (F(1)>0._R8P) then
         f(ut1) = f(1)*q_aux1(ut1)
         f(ut2) = f(1)*q_aux1(ut2)
         f(5  ) = f(5) + f(1) * (v1*v1 + w1*w1)/0.5_R8P
      else
         f(ut1) = f(1)*q_aux4(ut1)
         f(ut2) = f(1)*q_aux4(ut2)
         f(5  ) = f(5) + f(1) * (v4*v4 + w4*w4)/0.5_R8P
      endif
   case(4)
      call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux4,f=f)
   endselect
   endassociate
   endsubroutine compute_riemann_euler_hllc

   subroutine compute_riemann_euler_hllc_lm(si, sir, uni, ut1, ut2, nv, q_aux1, q_aux4, f, lmax)
   !< Solve the Riemann problem between the state 1 (left) and 4 (right) using the HLLC-LM solver.
   !< The algorithm is based on the Low Mach (LM) modification, see "A shock-stable modification of the HLLC Riemann solver with
   !< reduced numerical dissipation", Nico Fleischmann, Stefan Adami, Nikolaus A. Adams, 2020.
   integer(I4P), intent(in)            :: si(3)                    !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)                   !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: uni, ut1, ut2            !< Index of normal and tangential velocities.
   integer(I4P), intent(in)            :: nv                       !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)               !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)               !< Left state, auxiliary variables.
   real(R8P),    intent(inout)         :: f(1:)                    !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax                     !< Maximum wave speed estimation.
   real(R8P)                           :: Q1(1:nv),F1(1:nv)        !< State/fluxes 1.
   real(R8P)                           :: Q4(1:nv),F4(1:nv)        !< State/fluxes 4.
   real(R8P)                           :: Q1S(1:3)                 !< Intermediate 1 state.
   real(R8P)                           :: Q4S(1:3)                 !< Intermediate 4 state.
   real(R8P)                           :: S                        !< Velocity of the intermediate states.
   real(R8P)                           :: S1                       !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: S4                       !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: u_roe, a_roe             !< Roe averages.
   real(R8P)                           :: sqrt_r1                  !< Sqrt(r1).
   real(R8P)                           :: sqrt_r4                  !< Sqrt(r4).
   real(R8P)                           :: sqrt_r14                 !< Sqrt(r1) + Sqrt(r4).
   real(R8P)                           :: E1,E4                    !< Energy of state 1 and 4.
   real(R8P)                           :: Mr                       !< Local Mach ratio.
   real(R8P)                           :: phi                      !< Smoothing function.
   real(R8P), parameter                :: PI_2=2._R8P*ATAN(1._R8P) !< PI/2.
   real(R8P), parameter                :: MLIMIT=0.1_R8P           !< Local Mach limit.

   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             v1=>q_aux1(ut1), v4=>q_aux4(ut1), &
             w1=>q_aux1(ut2), w4=>q_aux4(ut2), &
             p1=>q_aux1(7  ), p4=>q_aux4(7  ))
   sqrt_r1 = sqrt(r1)
   sqrt_r4 = sqrt(r4)
   sqrt_r14 = sqrt_r1 + sqrt_r4
   u_roe = (u1*sqrt_r1 + u4*sqrt_r4)/sqrt_r14
   a_roe = sqrt((a1**2*sqrt_r1+a4**2*sqrt_r4)/sqrt_r14 + 0.5_R8P*sqrt_r1*sqrt_r4*((u4-u1)/sqrt_r14)**2)
   S1 = min(u1-a1, u_roe-a_roe)
   S4 = max(u4+a4, u_roe+a_roe)
   if (present(lmax)) lmax = max(abs(S1),abs(S4))
   S = (r4*u4*(S4-u4) - r1*u1*(S1-u1) + p1 - p4)/(r4*(S4-u4) - r1*(S1-u1))
   if     (S1>0._R8P) then
      call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux1,f=f)
   elseif (S4<0._R8P) then
      call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux4,f=f)
   else
      call compute_conservatives_euler(q_aux=q_aux1,q=Q1)
      call compute_conservatives_euler(q_aux=q_aux4,q=Q4)
      call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux1,f=F1)
      call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux4,f=F4)
      E1 = q_aux1(8) - q_aux1(7)/q_aux1(1)
      E4 = q_aux4(8) - q_aux4(7)/q_aux4(1)
      Q1S(1) = r1*(S1-u1)/(S1-S)
      Q1S(2) = Q1S(1)*S
      Q1S(3) = Q1S(1)*(E1+(S-u1)*(S+p1/(r1*(S1-u1))))
      Q4S(1) = r4*(S4-u4)/(S4-S)
      Q4S(2) = Q4S(1)*S
      Q4S(3) = Q4S(1)*(E4+(S-u4)*(S+p4/(r4*(S4-u4))))
      Mr = max(abs(u1/a1),abs(u4/a4))/MLIMIT
      phi = sin(min(1._R8P,Mr)*PI_2)
      S1 = phi*S1
      S4 = phi*S4
      f(1  ) = 0.5_R8P*(F1(1  )+F4(1  )) + 0.5_R8P*(S1*(Q1S(1)-Q1(1  )) + abs(s)*(Q1S(1)-Q4S(1)) + S4*(Q4S(1)-Q4(1  )))
      f(uni) = 0.5_R8P*(F1(uni)+F4(uni)) + 0.5_R8P*(S1*(Q1S(2)-Q1(uni)) + abs(s)*(Q1S(2)-Q4S(2)) + S4*(Q4S(2)-Q4(uni)))
      f(5  ) = 0.5_R8P*(F1(5  )+F4(5  )) + 0.5_R8P*(S1*(Q1S(3)-Q1(5  )) + abs(s)*(Q1S(3)-Q4S(3)) + S4*(Q4S(3)-Q4(5  )))
      if (f(1)>0._R8P) then
         f(ut1) = f(1)*v1
         f(ut2) = f(1)*w1
         f(5  ) = f(5) + f(1)*(v1*v1+w1*w1)/0.5_R8P
      else
         f(ut1) = f(1)*v4
         f(ut2) = f(1)*w4
         f(5  ) = f(5) + f(1)*(v4*v4+w4*w4)/0.5_R8P
      endif
   endif
   endassociate
   endsubroutine compute_riemann_euler_hllc_lm

   subroutine compute_riemann_euler_hllem(si, sir, uni, ut1, ut2, nv, q_aux1, q_aux4, f, lmax)
   !< Solve the Riemann problem between the state 1 (left) and 4 (right) using the HLLEM solver.
   !< See Dumbser and Balsara "A new efficient formulation of the HLLEM Riemann solver for general conservative and non-conservative
   !< hyperbolic systems", 2016, Journal of Compuational Physics.
   integer(I4P), intent(in)            :: si(3)                   !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)                  !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: uni, ut1, ut2           !< Index of normal and tangential velocities.
   integer(I4P), intent(in)            :: nv                      !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)              !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)              !< Left state, auxiliary variables.
   real(R8P),    intent(inout)         :: f(1:)                   !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax                    !< Maximum wave speed estimation.
   real(R8P)                           :: FN(3)                   !< Normal fluxes at interface.
   real(R8P)                           :: Q1N(3),F1N(3)           !< Normal State/fluxes 1.
   real(R8P)                           :: Q4N(3),F4N(3)           !< Normal State/fluxes 4.
   real(R8P)                           :: u_roe,a_roe,g_roe,H_roe !< Roe averages.
   real(R8P)                           :: sqrt_r1                 !< Sqrt(r1).
   real(R8P)                           :: sqrt_r4                 !< Sqrt(r4).
   real(R8P)                           :: sqrt_r14                !< Sqrt(r1) + Sqrt(r4).
   real(R8P)                           :: el(3,3)                 !< Left eigenvectors matrix.
   real(R8P)                           :: er(3,3)                 !< Right eigenvectors matrix.
   real(R8P)                           :: ev(3)                   !< Eigenvalues.
   real(R8P)                           :: evp(3)                  !< Positive part of igenvalues matrix.
   real(R8P)                           :: evm(3)                  !< Negative part of igenvalues matrix.
   real(R8P)                           :: delta(3,3)              !< Delta.
   real(R8P)                           :: AD(nv)                  !< Anti-diffusive term.

   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             g1=>q_aux1(5  ), g4=>q_aux4(5  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             v1=>q_aux1(ut1), v4=>q_aux4(ut1), &
             w1=>q_aux1(ut2), w4=>q_aux4(ut2), &
             p1=>q_aux1(7  ), p4=>q_aux4(7  ), &
             H1=>q_aux1(8  ), H4=>q_aux4(8  ))
   sqrt_r1 = sqrt(r1)
   sqrt_r4 = sqrt(r4)
   sqrt_r14 = sqrt_r1 + sqrt_r4
   u_roe = (u1*sqrt_r1 + u4*sqrt_r4)/sqrt_r14
   a_roe = sqrt((a1**2*sqrt_r1+a4**2*sqrt_r4)/sqrt_r14 + 0.5_R8P*sqrt_r1*sqrt_r4*((u4-u1)/sqrt_r14)**2)
   g_roe = (g1*sqrt_r1 + g4*sqrt_r4)/sqrt_r14
   H_roe = (H1*sqrt_r1 + H4*sqrt_r4)/sqrt_r14
   call compute_eigeinvectors_left(gm1=g_roe-1._R8P, u=u_roe, a=a_roe, el=el)
   call compute_eigeinvectors_right(H=H_roe, u=u_roe, a=a_roe, er=er)
   ev(1) = minval([0._R8P, u1-a1, u_roe-a_roe])
   ev(3) = maxval([0._R8P, u4+a4, u_roe+a_roe])
   if (present(lmax)) lmax = max(abs(ev(1)),abs(ev(3)))
   ev(2) = (r4*u4*(ev(3)-u4) - r1*u1*(ev(1)-u1) + p1 - p4)/(r4*(ev(3)-u4) - r1*(ev(1)-u1))
   evp = 0.5_R8P*(ev + abs(ev)) ; evm = 0.5_R8P*(ev - abs(ev))
   delta = 0._R8P
   delta(1,1) = 1._R8P - evm(1)/(ev(1) - 1e-14_R8P) - evp(1)/(ev(3) + 1e-14_R8P)
   delta(2,2) = 1._R8P - evm(2)/(ev(1) - 1e-14_R8P) - evp(2)/(ev(3) + 1e-14_R8P)
   delta(3,3) = 1._R8P - evm(3)/(ev(1) - 1e-14_R8P) - evp(3)/(ev(3) + 1e-14_R8P)
   Q1N(1) = r1            ; Q4N(1) = r4
   Q1N(2) = r1*u1         ; Q4N(2) = r4*u4
   Q1N(3) = r1*H1 - p1    ; Q4N(3) = r4*H4 - p4
   F1N(1) = r1*u1         ; F4N(1) = r4*u4
   F1N(2) = r1*u1*u1 + p1 ; F4N(2) = r4*u4*u4 + p4
   F1N(3) = r1*u1*H1      ; F4N(3) = r4*u4*H4
   AD = ev(3)*ev(1)/(ev(3)-ev(1))*matmul(er,matmul(delta,matmul(el,Q4N-Q1N)))
   FN = ((ev(3)*F1N - ev(1)*F4N) + ev(1)*ev(3)*(Q4N-Q1N))/(ev(3)-ev(1))
   FN = FN - AD
   f(1  ) = FN(1)
   f(uni) = FN(2)
   f(5  ) = FN(3)
   if (f(1)>0._R8P) then
      f(ut1) = f(1)*v1
      f(ut2) = f(1)*w1
      f(5  ) = f(5) + f(1) * (v1*v1 + w1*w1)/0.5_R8P
   else
      f(ut1) = f(1)*v4
      f(ut2) = f(1)*w4
      f(5  ) = f(5) + f(1) * (v4*v4 + w4*w4)/0.5_R8P
   endif
   endassociate
   endsubroutine compute_riemann_euler_hllem

   subroutine compute_riemann_euler_llf(si, sir, uni, ut1, ut2, nv, q_aux1, q_aux4, f, lmax)
   !< Solve the Riemann problem between the state 1 (left) and 4 (right) using the Local-Lax-Friedrichs (LLF, Rusanov) solver.
   integer(I4P), intent(in)            :: si(3)         !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)        !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: uni, ut1, ut2 !< Index of normal and tangential velocities.
   integer(I4P), intent(in)            :: nv            !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)    !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)    !< Left state, auxiliary variables.
   real(R8P),    intent(inout)         :: f(1:)         !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax          !< Maximum wave speed estimation.
   real(R8P)                           :: lmax_         !< Maximum wave speed estimation, local variable.
   real(R8P)                           :: Q1(1:nv)      !< State 1.
   real(R8P)                           :: Q4(1:nv)      !< State 4.
   real(R8P)                           :: F1(1:nv)      !< State 1 fluxes.
   real(R8P)                           :: F4(1:nv)      !< State 4 fluxes.
   real(R8P)                           :: S1            !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: S4            !< Maximum wave speed of state 1 and 4.
   real(R8P)                           :: u_roe, a_roe  !< Roe averages.
   real(R8P)                           :: sqrt_r1       !< Sqrt(r1).
   real(R8P)                           :: sqrt_r4       !< Sqrt(r4).
   real(R8P)                           :: sqrt_r14      !< Sqrt(r1) + Sqrt(r4).

   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             v1=>q_aux1(ut1), v4=>q_aux4(ut1), &
             w1=>q_aux1(ut2), w4=>q_aux4(ut2))
   sqrt_r1 = sqrt(r1)
   sqrt_r4 = sqrt(r4)
   sqrt_r14 = sqrt_r1 + sqrt_r4
   u_roe = (u1*sqrt_r1 + u4*sqrt_r4)/sqrt_r14
   a_roe = sqrt((a1**2*sqrt_r1+a4**2*sqrt_r4)/sqrt_r14 + 0.5_R8P*sqrt_r1*sqrt_r4*((u4-u1)/sqrt_r14)**2)
   S1 = minval([0._R8P, u1-a1, u_roe-a_roe])
   S4 = maxval([0._R8P, u4+a4, u_roe+a_roe])
   lmax_ = max(abs(S1),abs(S4))
   if (present(lmax)) lmax = lmax_
   call compute_conservatives_euler(q_aux=q_aux1,q=Q1)
   call compute_conservatives_euler(q_aux=q_aux4,q=Q4)
   call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux1,f=F1)
   call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux4,f=F4)
   f(1  ) = 0.5_R8P*(F1(1  ) + F4(1  ) - lmax_*(Q4(1  )-Q1(1  )))
   f(uni) = 0.5_R8P*(F1(uni) + F4(uni) - lmax_*(Q4(uni)-Q1(uni)))
   f(5  ) = 0.5_R8P*(F1(5  ) + F4(5  ) - lmax_*(Q4(5  )-Q1(5  )))
   if (f(1)>0._R8P) then
      f(ut1) = f(1)*v1
      f(ut2) = f(1)*w1
      f(5  ) = f(5) + f(1)*(v1*v1+w1*w1)/0.5_R8P
   else
      f(ut1) = f(1)*v4
      f(ut2) = f(1)*w4
      f(5  ) = f(5) + f(1)*(v4*v4+w4*w4)/0.5_R8P
   endif
   endassociate
   endsubroutine compute_riemann_euler_llf

   subroutine compute_riemann_euler_ts(si, sir, uni, ut1, ut2, nv, q_aux1, q_aux4, f, lmax)
   !< Solve the Riemann problem between the state 1 (left) and 4 (right) using the Two-Shocks (TS) approximation.
   integer(I4P), intent(in)            :: si(3)           !< Directional (1=x,2=y,3=z) increment.
   real(R8P),    intent(in)            :: sir(3)          !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: uni, ut1, ut2   !< Index of normal and tangential velocities.
   integer(I4P), intent(in)            :: nv              !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q_aux1(1:)      !< Left state, auxiliary variables.
   real(R8P),    intent(in)            :: q_aux4(1:)      !< Left state, auxiliary variables.
   real(R8P),    intent(inout)         :: f(1:)           !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax            !< Maximum wave speed estimation.
   real(R8P)                           :: q_auxs(1:9)     !< Intermediate state.
   real(R8P)                           :: u23,p23,r2,r3   !< Intermediate states.
   real(R8P)                           :: S,S1,S2,S3,S4   !< Waves speeds.
   real(R8P)                           :: aa1,bb1,aa4,bb4 !< Dummies coefficients.
   real(R8P)                           :: g1p,g4p         !< Dummy variables for computing TS pressure of intermediate states.
   real(R8P)                           :: gp1,gm1,gm1_gp1 !< g+1, (g-1), (g-1)/(g+1).
   real(R8P)                           :: p_p             !< Pressure ratio.

   associate(a1=>q_aux1(9  ), a4=>q_aux4(9  ), &
             r1=>q_aux1(1  ), r4=>q_aux4(1  ), &
             g1=>q_aux1(5  ), g4=>q_aux4(5  ), &
             u1=>q_aux1(uni), u4=>q_aux4(uni), &
             v1=>q_aux1(ut1), v4=>q_aux4(ut1), &
             w1=>q_aux1(ut2), w4=>q_aux4(ut2), &
             p1=>q_aux1(7  ), p4=>q_aux4(7  ))
      ! PVL approximation of p23
      aa1 = 2._R8P/((g1+1._R8P)*r1) ; bb1 = p1*(g1-1._R8P)/(g1+1._R8P)
      aa4 = 2._R8P/((g4+1._R8P)*r4) ; bb4 = p4*(g4-1._R8P)/(g4+1._R8P)
      p23 = 0.5_R8P * (p1+p4) - 0.125_R8P *(u4-u1)*(r1+r4)*(a1+a4)
      p23 = (sqrt(aa1/(p23+bb1))*p1 + sqrt(aa4/(p23+bb4))*p4 + u1 - u4)/(sqrt(aa1/(p23+bb1)) + sqrt(aa4/(p23+bb4)))
      ! TS approximation of p23 u23
      gp1 = g1 + 1._R8P ; gm1_gp1 = (g1 - 1._R8P)/gp1 ; g1p = sqrt((2._R8P/(gp1*r1))/(p23 + gm1_gp1*p1))
      gp1 = g4 + 1._R8P ; gm1_gp1 = (g4 - 1._R8P)/gp1 ; g4p = sqrt((2._R8P/(gp1*r4))/(p23 + gm1_gp1*p4))
      p23 = (g1p*p1 + g4p*p4 + u1 - u4)/(g1p + g4p)
      u23 = 0.5_R8P*(u1 + u4 + (p23 - p4)*g4p - (p23 - p1)*g1p)
      ! evaluating left state
      p_p = p23/p1
      gp1 = g1 + 1._R8P
      gm1 = g1 - 1._R8P
      r2  = r1*((gm1 + gp1*p_p)/(gp1 + gm1*p_p))
      S1  = u1  - a1*sqrt(1._R8P + 0.5_R8P*gp1/g1*(p_p-1._R8P))
      S2  = S1
      ! evaluating right state
      p_p = p23/p4
      gp1 = g4 + 1._R8P
      gm1 = g4 - 1._R8P
      r3  = r4*((gm1 + gp1*p_p)/(gp1 + gm1*p_p))
      S4  = u4  + a4*sqrt(1._R8P + 0.5_R8P*gp1/g4*(p_p-1._R8P))
      S3  = S4
      if (present(lmax)) lmax = max(abs(S1),abs(S4))
      select case(minloc([-S1,S1*S,S*S4,S4],dim=1))
      case(1)
         call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux1,f=f)
      case(2)
         q_auxs(1  ) = r2
         q_auxs(uni) = S
         q_auxs(ut1) = v1
         q_auxs(ut2) = w1
         q_auxs(7  ) = p23
         q_auxs(8  ) = g1*p23/((g1-1._R8P)*r2)+0.5_R8P*(S*S+v1*v1+w1*w1)
         call compute_convective_fluxes_euler(sir=sir,q_aux=q_auxs,f=f)
      case(3)
         q_auxs(1  ) = r3
         q_auxs(uni) = S
         q_auxs(ut1) = v4
         q_auxs(ut2) = w4
         q_auxs(7  ) = p23
         q_auxs(8  ) = g4*p23/((g4-1._R8P)*r3)+0.5_R8P*(S*S+v4*v4+w4*w4)
         call compute_convective_fluxes_euler(sir=sir,q_aux=q_auxs,f=f)
      case(4)
         call compute_convective_fluxes_euler(sir=sir,q_aux=q_aux4,f=f)
      endselect
   endassociate
   endsubroutine compute_riemann_euler_ts

   ! private procedures
   pure subroutine compute_conservatives_euler(q_aux,q)
   !< Compute convervative variables from auxiliary ones.
   real(R8P),    intent(in)    :: q_aux(1:) !< Auxiliary variables.
   real(R8P),    intent(inout) :: q(1:)     !< Conservative varibales.

   q(1) =      q_aux(1)
   q(2) = q(1)*q_aux(2)
   q(3) = q(1)*q_aux(3)
   q(4) = q(1)*q_aux(4)
   q(5) = q(1)*q_aux(8) - q_aux(7)
   endsubroutine compute_conservatives_euler

   pure subroutine compute_convective_fluxes_euler(sir,q_aux,f)
   !< Compute convective fluxes for Euler equations from auxiliary variables.
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
   endsubroutine compute_convective_fluxes_euler

    pure subroutine compute_eigeinvectors_left(gm1, u, a, el)
    !< Compute left eigenvectors matrix (L) of the Jacobian fluxes matrix $A=R \Lambda L$ in primitive variables form.
    !< @note This function consider only the normal direction:
    !< $\frac{\partial P}{\partial t} + R \Lambda L \frac{\partial P}{\partial n} = 0$ where P are the primitive variables and
    !< n is the normal direction. R is the matrix of the right eigenvectors, $\Lambda$ is the diagonal matrix of the eigenvalues
    !< and L is the matrix of the left eigenvectors.
    real(R8P), intent(in)    :: gm1     !< Specific heats ration minus one, g-1.
    real(R8P), intent(in)    :: u       !< Normal velocity.
    real(R8P), intent(in)    :: a       !< Speed of sound.
    real(R8P), intent(inout) :: el(3,3) !< Left eigenvectors matrix.
    real(R8P)                :: gm1_a2  !< g-1/a^2.
    real(R8P)                :: u2      !< u^2.

    gm1_a2 = gm1/a/a
    u2 = u*u
    el(1,1) =         0.25_R8P*gm1_a2*u2 + 0.5_R8P*u/a ; el(1,2) = -gm1_a2*u/a - 0.5_R8P/a ; el(1,3) = 0.5_R8P*gm1_a2
    el(2,1) = 1._R8P - 0.5_R8P*gm1_a2*u2               ; el(2,2) =  gm1_a2*u               ; el(2,3) =        -gm1_a2
    el(3,1) =         0.25_R8P*gm1_a2*u2 + 0.5_R8P*u/a ; el(3,2) = -gm1_a2*u/a + 0.5_R8P/a ; el(3,3) = el(1,3)
    ! el = 0._R8P
    ! el(1,1) = 1._R8P
    ! el(2,2) = 1._R8P
    ! el(3,3) = 1._R8P
    endsubroutine compute_eigeinvectors_left

    pure subroutine compute_eigeinvectors_right(H, u, a, er)
    !< Compute right eigenvectors matrix (R) of the Jacobian fluxes matrix $A=R \Lambda L$ in primitive variables form.
    !< @note This function consider only the normal direction:
    !< $\frac{\partial P}{\partial t} + R \Lambda L \frac{\partial P}{\partial n} = 0$ where P are the primitive variables and
    !< n is the normal direction. R is the matrix of the right eigenvectors, $\Lambda$ is the diagonal matrix of the eigenvalues
    !< and L is the matrix of the left eigenvectors.
    real(R8P), intent(in)    :: H       !< Total specific entalpy.
    real(R8P), intent(in)    :: u       !< Normal velocity.
    real(R8P), intent(in)    :: a       !< Speed of sound.
    real(R8P), intent(inout) :: er(3,3) !< Right eigenvectors matrix.

    er(1,1) = 1._R8P  ; er(1,2) = 1._R8P      ; er(1,3) = 1._R8P
    er(2,1) = u - a   ; er(2,2) = u           ; er(2,3) = u + a
    er(3,1) = H - u*a ; er(3,2) = 0.5_R8P*u*u ; er(3,3) = H + u*a
    ! er = 0._R8P
    ! er(1,1) = 1._R8P
    ! er(2,2) = 1._R8P
    ! er(3,3) = 1._R8P
    endsubroutine compute_eigeinvectors_right
endmodule adam_riemann_euler_library
