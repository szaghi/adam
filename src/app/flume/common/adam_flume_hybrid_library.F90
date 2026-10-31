!< ADAM, FLUME pointwise central flux of the hybrid central/WENO scheme (issue #119, T2).
module adam_flume_hybrid_library
!< ADAM, FLUME pointwise central flux of the hybrid central/WENO scheme (issue #119, T2).
!<
!< The 6th-order kinetic-energy- and pressure-equilibrium-preserving face flux, the conservative form of Pirozzoli (2011) with the
!< two-point flux of Shima et al. (2021, KEEP-PE), `<a> = (a_l + a_r) / 2`:
!<```
!< F_{i+1/2} = 2 sum_{l=1..3} a_l sum_{m=0..l-1} f(i-m, i-m+l),   a = (3/4, -3/20, 1/60)
!< f(l, r): mass <rho><u_n>; momentum <rho><u_n><u_i> + <p> delta_{ni};
!<          energy <rho><u_n> (u_l . u_r) / 2 + <p><u_n> / (gamma - 1) + (u_n,l p_r + u_n,r p_l) / 2
!<```
!< A uniform velocity and pressure stay so exactly across any density jump: the form of Kennedy-Gruber-Pirozzoli `<rho><u_n><H>` does
!< not, and a 100:1 contact stopped the run at its first step. Written from the published form through the llm-wiki and memory, not from
!< the papers themselves; test_kep.F90 measures the order, the kinetic energy identity and the pressure equilibrium. The state is
!< `w(5,-2:3)`, (rho, u, v, w, p) of the 6 cells around the face, cell 0 on its left. `pure`, `!$acc routine seq`, `!$omp declare target`.

! third party modules
use :: penf, only : I4P, R8P

implicit none
private
public :: compute_ducros
public :: compute_kep_flux

contains
   pure subroutine compute_kep_flux(d, gamma, w, f)
   !< Compute the 6th-order KEP face flux of direction `d`.
   integer(I4P), intent(in)  :: d            !< Direction, 1=x, 2=y, 3=z.
   real(R8P),    intent(in)  :: gamma        !< Specific heats ratio.
   real(R8P),    intent(in)  :: w(5,-2:3)    !< Cell states.
   real(R8P),    intent(out) :: f(5)         !< Face flux.
   real(R8P)                 :: a(3)         !< Coefficients.
   real(R8P)                 :: uu           !< u_l . u_r.
   real(R8P)                 :: mf, pf       !< Two-point mass flux and pressure.
   real(R8P)                 :: rg           !< 1 / (gamma - 1).
   integer(I4P)              :: l, m, n, c, r !< Counters.
   !$acc routine seq
   !$omp declare target

   a(1) = 0.75_R8P ; a(2) = -0.15_R8P ; a(3) = 1._R8P / 60._R8P
   rg = 1._R8P / (gamma - 1._R8P)
   do n=1, 5
      f(n) = 0._R8P
   enddo
   do l=1, 3
      do m=0, l-1
         c = -m ; r = -m + l
         mf = 0.25_R8P * (w(1,c) + w(1,r)) * (w(1+d,c) + w(1+d,r))
         pf = 0.5_R8P * (w(5,c) + w(5,r))
         f(1) = f(1) + 2._R8P * a(l) * mf
         do n=1, 3
            f(1+n) = f(1+n) + 2._R8P * a(l) * mf * 0.5_R8P * (w(1+n,c) + w(1+n,r))
         enddo
         f(1+d) = f(1+d) + 2._R8P * a(l) * pf
         uu = w(2,c) * w(2,r) + w(3,c) * w(3,r) + w(4,c) * w(4,r)
         f(5) = f(5) + 2._R8P * a(l) * (0.5_R8P * mf * uu + 0.5_R8P * rg * pf * (w(1+d,c) + w(1+d,r)) + &
                                        0.5_R8P * (w(1+d,c) * w(5,r) + w(1+d,r) * w(5,c)))
      enddo
   enddo
   endsubroutine compute_kep_flux

   pure subroutine compute_ducros(eps, jump, jump_max, g, psi)
   !< Compute the Ducros et al. (1999) sensor, `theta^2 / (theta^2 + |omega|^2 + eps)`, from the velocity gradient `g(a,b) = d u_a / d x_b`;
   !< 1 where the density jumps: the Ducros ratio cannot see a contact, where the velocity is uniform and a central flux rings.
   real(R8P), intent(in)  :: eps    !< Floor of the denominator.
   real(R8P), intent(in)  :: jump   !< Largest |rho_+ - rho_-| / (rho_+ + rho_-) of the cell along its active directions.
   real(R8P), intent(in)  :: jump_max !< Its bound.
   real(R8P), intent(in)  :: g(3,3) !< Velocity gradient.
   real(R8P), intent(out) :: psi    !< Sensor.
   real(R8P)              :: th, wx, wy, wz !< Dilatation and vorticity.
   !$acc routine seq
   !$omp declare target

   th = g(1,1) + g(2,2) + g(3,3)
   wx = g(3,2) - g(2,3) ; wy = g(1,3) - g(3,1) ; wz = g(2,1) - g(1,2)
   psi = th * th / (th * th + wx * wx + wy * wy + wz * wz + eps)
   if (jump > jump_max) psi = 1._R8P
   endsubroutine compute_ducros
endmodule adam_flume_hybrid_library
