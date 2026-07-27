!< ADAM, PRISM Riemann Problem for Maxwell equations solvers and convective fluxes computations library.
module adam_prism_riemann_library
!< ADAM, PRISM Riemann Problem for Maxwell equations solvers and convective fluxes computations library.
!< @NOTE see `adam_prism_physics_object` for the field variables arrangment.

! PRISM modules
use :: adam_prism_parameters
use :: adam_prism_physics_object
! third party modules
use :: penf, only : I4P, R8P

implicit none
private
public :: compute_convective_fluxes_interface
public :: compute_riemann_maxwell_llf
public :: compute_convective_fluxes_maxwell
public :: compute_convective_fluxes_maxwell_adim
public :: compute_convective_fluxes_maxwell_div_d
public :: compute_convective_fluxes_maxwell_div_b
public :: compute_convective_fluxes_maxwell_div_d_b
public :: compute_convective_fluxes_maxwell_adim_div_d
public :: compute_convective_fluxes_maxwell_adim_div_b
public :: compute_convective_fluxes_maxwell_adim_div_d_b
public :: compute_eigenvalues_vector

interface
   subroutine compute_convective_fluxes_interface(sir,q,f,chi)
   import :: R8P
   real(R8P), intent(in)    :: sir(3) !< Directional (1=x,2=y,3=z) increment.
   real(R8P), intent(in)    :: q(1:)  !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)  !< Conservative fluxes.
   real(R8P), intent(in)    :: chi    !< Coefficiente velocità trasporto errori divergenza campi
   endsubroutine compute_convective_fluxes_interface
endinterface

contains
   subroutine compute_riemann_maxwell_llf(sir, nv, q1, q4, chi, f, lmax)
   !< Solve the Riemann problem between the state 1 (left) and 4 (right) using the Local-Lax-Friedrichs (LLF, Rusanov) solver.
   real(R8P),    intent(in)            :: sir(3)   !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: nv       !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q1(1:)   !< Left state.
   real(R8P),    intent(in)            :: q4(1:)   !< Right state.
   real(R8P),    intent(in)            :: chi      !< Coefficiente velocità trasporto errori divergenza campi (non usato)
   real(R8P),    intent(inout)         :: f(1:)    !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax     !< Maximum wave speed estimation.
   real(R8P)                           :: F1(1:nv) !< State 1 fluxes.
   real(R8P)                           :: F4(1:nv) !< State 4 fluxes.

   lmax = sqrt(1._R8P/(EPS0*MU0))
   call compute_convective_fluxes_maxwell(sir=sir,q=q1,f=F1,chi=chi)
   call compute_convective_fluxes_maxwell(sir=sir,q=q4,f=F4,chi=chi)
   f(1) = 0.5_R8P*(F1(1)+F4(1)-lmax*(q4(1)-q1(1)))
   f(2) = 0.5_R8P*(F1(2)+F4(2)-lmax*(q4(2)-q1(2)))
   f(3) = 0.5_R8P*(F1(3)+F4(3)-lmax*(q4(3)-q1(3)))
   f(4) = 0.5_R8P*(F1(4)+F4(4)-lmax*(q4(4)-q1(4)))
   f(5) = 0.5_R8P*(F1(5)+F4(5)-lmax*(q4(5)-q1(5)))
   f(6) = 0.5_R8P*(F1(6)+F4(6)-lmax*(q4(6)-q1(6)))
   endsubroutine compute_riemann_maxwell_llf

   subroutine compute_riemann_maxwell_hll(sir, nv, q1, q4, chi, f, lmax, lmin)
   !< Solve the Riemann problem between the state 1 (left) and 4 (right) using the Harten, Lax and van Leer (HLL) solver.
   real(R8P),    intent(in)            :: sir(3)   !< Directional (1=x,2=y,3=z) increment.
   integer(I4P), intent(in)            :: nv       !< Number of conservative varibales.
   real(R8P),    intent(in)            :: q1(1:)   !< Left state.
   real(R8P),    intent(in)            :: q4(1:)   !< Right state.
   real(R8P),    intent(in)            :: chi      !< Coefficiente velocità trasporto errori divergenza campi (non usato)
   real(R8P),    intent(inout)         :: f(1:)    !< Resulting fluxes.
   real(R8P),    intent(out), optional :: lmax     !< Maximum wave speed estimation.
   real(R8P),    intent(out), optional :: lmin     !< Maximum wave speed estimation.
   real(R8P)                           :: F1(1:nv) !< State 1 fluxes.
   real(R8P)                           :: F4(1:nv) !< State 4 fluxes.

   lmax = sqrt(1._R8P/(EPS0*MU0))
   lmin = -lmax
   call compute_convective_fluxes_maxwell(sir=sir,q=q1,f=F1,chi=chi)
   call compute_convective_fluxes_maxwell(sir=sir,q=q4,f=F4,chi=chi)
   f(1) = (lmax*F1(1)-lmin*F4(1)+lmax*lmin*(F4(1)-F1(1)))/(lmax-lmin)
   f(2) = (lmax*F1(2)-lmin*F4(2)+lmax*lmin*(F4(2)-F1(2)))/(lmax-lmin)
   f(3) = (lmax*F1(3)-lmin*F4(3)+lmax*lmin*(F4(3)-F1(3)))/(lmax-lmin)
   f(4) = (lmax*F1(4)-lmin*F4(4)+lmax*lmin*(F4(4)-F1(4)))/(lmax-lmin)
   f(5) = (lmax*F1(5)-lmin*F4(5)+lmax*lmin*(F4(5)-F1(5)))/(lmax-lmin)
   f(6) = (lmax*F1(6)-lmin*F4(6)+lmax*lmin*(F4(6)-F1(6)))/(lmax-lmin)
   endsubroutine compute_riemann_maxwell_hll

   subroutine compute_convective_fluxes_maxwell(sir,q,f,chi)
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
   real(R8P), intent(in)    :: sir(3) !< Direction array
   real(R8P), intent(in)    :: q(1:)  !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)  !< Conservative fluxes.
   real(R8P), intent(in)    :: chi    !< Coefficiente velocità trasporto errori divergenza campi (non usato)
   !$acc routine seq

   f(VAR_DX) =  0.0_R8P        * sir(1) + -q(VAR_BZ)/MU0  * sir(2) +  q(VAR_BY)/MU0  * sir(3)
   f(VAR_DY) =  q(VAR_BZ)/MU0  * sir(1) +  0.0_R8P        * sir(2) + -q(VAR_BX)/MU0  * sir(3)
   f(VAR_DZ) = -q(VAR_BY)/MU0  * sir(1) +  q(VAR_BX)/MU0  * sir(2) +  0.0_R8P        * sir(3)
   f(VAR_BX) =  0.0_R8P        * sir(1) +  q(VAR_DZ)/EPS0 * sir(2) + -q(VAR_DY)/EPS0 * sir(3)
   f(VAR_BY) = -q(VAR_DZ)/EPS0 * sir(1) +  0.0_R8P        * sir(2) +  q(VAR_DX)/EPS0 * sir(3)
   f(VAR_BZ) =  q(VAR_DY)/EPS0 * sir(1) + -q(VAR_DX)/EPS0 * sir(2) +  0.0_R8P        * sir(3)
   ! if (sir(1).eq.1._R8P) then !X
   !    f(VAR_DX) =  0.0_R8P
   !    f(VAR_DY) =  q(6)/MU0
   !    f(VAR_DZ) = -q(5)/MU0
   !    f(VAR_BX) =  0.0_R8P
   !    f(VAR_BY) = -q(3)/EPS0
   !    f(VAR_BZ) =  q(2)/EPS0
   ! elseif(sir(2).eq.1._R8P) then !Y
   !    f(VAR_DX) = -q(6)/MU0
   !    f(VAR_DY) =  0.0_R8P
   !    f(VAR_DZ) =  q(4)/MU0
   !    f(VAR_BX) =  q(3)/EPS0
   !    f(VAR_BY) =  0.0_R8P
   !    f(VAR_BZ) = -q(1)/EPS0
   ! elseif(sir(3).eq.1._R8P) then  !Z
   !    f(VAR_DX) =  q(5)/MU0
   !    f(VAR_DY) = -q(4)/MU0
   !    f(VAR_DZ) =  0.0_R8P
   !    f(VAR_BX) = -q(2)/EPS0
   !    f(VAR_BY) =  q(1)/EPS0
   !    f(VAR_BZ) =  0.0_R8P
   ! endif
   endsubroutine compute_convective_fluxes_maxwell

   subroutine compute_convective_fluxes_maxwell_adim(sir,q,f,chi)
   !< Compute convective fluxes for the adimensional Maxwell system.
   real(R8P), intent(in)    :: sir(3) !< Direction array
   real(R8P), intent(in)    :: q(1:)  !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)  !< Conservative fluxes.
   real(R8P), intent(in)    :: chi    !< Coefficiente velocità trasporto errori divergenza campi (non usato)
   !$acc routine seq

   f(VAR_DX) =  0.0_R8P      * sir(1) + -q(VAR_BZ) * sir(2) +  q(VAR_BY) * sir(3)
   f(VAR_DY) =  q(VAR_BZ)    * sir(1) +  0.0_R8P   * sir(2) + -q(VAR_BX) * sir(3)
   f(VAR_DZ) = -q(VAR_BY)    * sir(1) +  q(VAR_BX) * sir(2) +  0.0_R8P   * sir(3)
   f(VAR_BX) =  0.0_R8P      * sir(1) +  q(VAR_DZ) * sir(2) + -q(VAR_DY) * sir(3)
   f(VAR_BY) = -q(VAR_DZ)    * sir(1) +  0.0_R8P   * sir(2) +  q(VAR_DX) * sir(3)
   f(VAR_BZ) =  q(VAR_DY)    * sir(1) + -q(VAR_DX) * sir(2) +  0.0_R8P   * sir(3)
   endsubroutine compute_convective_fluxes_maxwell_adim

   subroutine compute_convective_fluxes_maxwell_div_d(sir,q,f,chi)
   !< Compute convective fluxes.
   !< J corrected with phi:
   !<```
   !< Fx(1)  =  phi     Fy(1)  = -Bz/muz   Fz(1)  =  By/muy
   !< Fx(2)  =  Bz/muz  Fy(2)  =  phi      Fz(2)  = -Bx/mux
   !< Fx(3)  = -By/muy  Fy(3)  =  Bx/mux   Fz(3)  =  phi
   !< Fx(4)  =  0       Fy(4)  =  Dz/epsz  Fz(4)  = -Dy/epsy
   !< Fx(5)  = -Dz/epsz Fy(5)  =  0        Fz(5)  =  Dx/epsx
   !< Fx(6)  =  Dy/epsy Fy(6)  = -Dx/epsx  Fz(6)  =  0
   !< Fx(7)  =  0       Fy(7)  =  0        Fz(7)  =  0
   !< Fx(8)  =  0       Fy(8)  =  0        Fz(8)  =  0
   !< Fx(9)  =  0       Fy(9)  =  0        Fz(9)  =  0
   !< Fx(10) =  ch^2*Dx Fy(10) =  ch^2*Dy  Fz(10) =  ch^2*Dz
   !<```
   real(R8P), intent(in)    :: sir(3) !< Direction array
   real(R8P), intent(in)    :: q(1:)  !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)  !< Conservative fluxes.
   real(R8P), intent(in)    :: chi    !< Coefficiente velocità trasporto errori divergenza campi
   real(R8P)                :: ch     !< Velocità trasporto errori divergenza campi modello iperbolico
   !$acc routine seq

   ch = chi/sqrt(EPS0*MU0)
   f(VAR_DX) =  q(7_I4P)        * sir(1) + -q(VAR_BZ)/MU0   * sir(2) +  q(VAR_BY)/MU0   * sir(3)
   f(VAR_DY) =  q(VAR_BZ)/MU0   * sir(1) +  q(7_I4P)        * sir(2) + -q(VAR_BX)/MU0   * sir(3)
   f(VAR_DZ) = -q(VAR_BY)/MU0   * sir(1) +  q(VAR_BX)/MU0   * sir(2) +  q(7_I4P)        * sir(3)
   f(VAR_BX) =  0.0_R8P         * sir(1) +  q(VAR_DZ)/EPS0  * sir(2) + -q(VAR_DY)/EPS0  * sir(3)
   f(VAR_BY) = -q(VAR_DZ)/EPS0  * sir(1) +  0.0_R8P         * sir(2) +  q(VAR_DX)/EPS0  * sir(3)
   f(VAR_BZ) =  q(VAR_DY)/EPS0  * sir(1) + -q(VAR_DX)/EPS0  * sir(2) +  0.0_R8P         * sir(3)
   f(7_I4P ) =  ch**2*q(VAR_DX) * sir(1) +  ch**2*q(VAR_DY) * sir(2) +  ch**2*q(VAR_DZ) * sir(3)
   !
   !if (sir(1).eq.1._R8P) then !X
   !   f(1)  =  q(10)
   !   f(2)  =  q(6)/MU0
   !   f(3)  = -q(5)/MU0
   !   f(4)  =  0.0_R8P
   !   f(5)  = -q(3)/EPS0
   !   f(6)  =  q(2)/EPS0
   !   f(7)  =  0._R8P
   !   f(8)  =  0._R8P
   !   f(9)  =  0._R8P
   !   f(10) =  ch*ch*q(1)
   !elseif(sir(2).eq.1._R8P) then !Y
   !   f(1)  = -q(6)/MU0
   !   f(2)  =  q(10)
   !   f(3)  =  q(4)/MU0
   !   f(4)  =  q(3)/EPS0
   !   f(5)  =  0.0_R8P
   !   f(6)  = -q(1)/EPS0
   !   f(7)  =  0._R8P
   !   f(8)  =  0._R8P
   !   f(9)  =  0._R8P
   !   f(10) =  ch*ch*q(2)
   !elseif(sir(3).eq.1._R8P) then  !Z
   !   f(1)  =  q(5)/MU0
   !   f(2)  = -q(4)/MU0
   !   f(3)  =  q(10)
   !   f(4)  = -q(2)/EPS0
   !   f(5)  =  q(1)/EPS0
   !   f(6)  =  0.0_R8P
   !   f(7)  =  0._R8P
   !   f(8)  =  0._R8P
   !   f(9)  =  0._R8P
   !   f(10) =  ch*ch*q(3)
   !endif
   endsubroutine compute_convective_fluxes_maxwell_div_d

   subroutine compute_convective_fluxes_maxwell_adim_div_d(sir,q,f,chi)
   !< Compute convective fluxes for the adimensional Maxwell system with D cleaning.
   real(R8P), intent(in)    :: sir(3) !< Direction array
   real(R8P), intent(in)    :: q(1:)  !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)  !< Conservative fluxes.
   real(R8P), intent(in)    :: chi    !< Coefficiente velocità trasporto errori divergenza campi
   !$acc routine seq

   f(VAR_DX) =  q(7_I4P)        * sir(1) + -q(VAR_BZ)     * sir(2) +  q(VAR_BY)     * sir(3)
   f(VAR_DY) =  q(VAR_BZ)       * sir(1) +  q(7_I4P)      * sir(2) + -q(VAR_BX)     * sir(3)
   f(VAR_DZ) = -q(VAR_BY)       * sir(1) +  q(VAR_BX)     * sir(2) +  q(7_I4P)      * sir(3)
   f(VAR_BX) =  0.0_R8P         * sir(1) +  q(VAR_DZ)     * sir(2) + -q(VAR_DY)     * sir(3)
   f(VAR_BY) = -q(VAR_DZ)       * sir(1) +  0.0_R8P       * sir(2) +  q(VAR_DX)     * sir(3)
   f(VAR_BZ) =  q(VAR_DY)       * sir(1) + -q(VAR_DX)     * sir(2) +  0.0_R8P       * sir(3)
   f(7_I4P ) =  chi**2*q(VAR_DX)* sir(1) +  chi**2*q(VAR_DY)* sir(2) +  chi**2*q(VAR_DZ)* sir(3)
   endsubroutine compute_convective_fluxes_maxwell_adim_div_d

   subroutine compute_convective_fluxes_maxwell_div_b(sir,q,f,chi)
   !< Compute convective fluxes.
   !< J corrected with phi:
   !<```
   !< Fx(1)  =  phi     Fy(1)  = -Bz/muz   Fz(1)  =  By/muy
   !< Fx(2)  =  Bz/muz  Fy(2)  =  phi      Fz(2)  = -Bx/mux
   !< Fx(3)  = -By/muy  Fy(3)  =  Bx/mux   Fz(3)  =  phi
   !< Fx(4)  =  0       Fy(4)  =  Dz/epsz  Fz(4)  = -Dy/epsy
   !< Fx(5)  = -Dz/epsz Fy(5)  =  0        Fz(5)  =  Dx/epsx
   !< Fx(6)  =  Dy/epsy Fy(6)  = -Dx/epsx  Fz(6)  =  0
   !< Fx(7)  =  0       Fy(7)  =  0        Fz(7)  =  0
   !< Fx(8)  =  0       Fy(8)  =  0        Fz(8)  =  0
   !< Fx(9)  =  0       Fy(9)  =  0        Fz(9)  =  0
   !< Fx(10) =  ch^2*Dx Fy(10) =  ch^2*Dy  Fz(10) =  ch^2*Dz
   !<```
   real(R8P), intent(in)    :: sir(3) !< Direction array
   real(R8P), intent(in)    :: q(1:)  !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)  !< Conservative fluxes.
   real(R8P), intent(in)    :: chi    !< Coefficiente velocità trasporto errori divergenza campi
   real(R8P)                :: ch     !< Velocità trasporto errori divergenza campi modello iperbolico
   !$acc routine seq

   ch = chi/sqrt(EPS0*MU0)
   f(VAR_DX) =  0.0_R8P         * sir(1) + -q(VAR_BZ)/MU0   * sir(2) +  q(VAR_BY)/MU0   * sir(3)
   f(VAR_DY) =  q(VAR_BZ)/MU0   * sir(1) +  0.0_R8P         * sir(2) + -q(VAR_BX)/MU0   * sir(3)
   f(VAR_DZ) = -q(VAR_BY)/MU0   * sir(1) +  q(VAR_BX)/MU0   * sir(2) +  0.0_R8P         * sir(3)
   f(VAR_BX) =  q(7_I4P)        * sir(1) +  q(VAR_DZ)/EPS0  * sir(2) + -q(VAR_DY)/EPS0  * sir(3)
   f(VAR_BY) = -q(VAR_DZ)/EPS0  * sir(1) +  q(7_I4P)        * sir(2) +  q(VAR_DX)/EPS0  * sir(3)
   f(VAR_BZ) =  q(VAR_DY)/EPS0  * sir(1) + -q(VAR_DX)/EPS0  * sir(2) +  q(7_I4P)        * sir(3)
   f(7_I4P ) =  ch**2*q(VAR_BX) * sir(1) +  ch**2*q(VAR_BY) * sir(2) +  ch**2*q(VAR_BZ) * sir(3)
   !
   !if (sir(1).eq.1._R8P) then !X
   !   f(1)  =  q(10)
   !   f(2)  =  q(6)/MU0
   !   f(3)  = -q(5)/MU0
   !   f(4)  =  0.0_R8P
   !   f(5)  = -q(3)/EPS0
   !   f(6)  =  q(2)/EPS0
   !   f(7)  =  0._R8P
   !   f(8)  =  0._R8P
   !   f(9)  =  0._R8P
   !   f(10) =  ch*ch*q(1)
   !elseif(sir(2).eq.1._R8P) then !Y
   !   f(1)  = -q(6)/MU0
   !   f(2)  =  q(10)
   !   f(3)  =  q(4)/MU0
   !   f(4)  =  q(3)/EPS0
   !   f(5)  =  0.0_R8P
   !   f(6)  = -q(1)/EPS0
   !   f(7)  =  0._R8P
   !   f(8)  =  0._R8P
   !   f(9)  =  0._R8P
   !   f(10) =  ch*ch*q(2)
   !elseif(sir(3).eq.1._R8P) then  !Z
   !   f(1)  =  q(5)/MU0
   !   f(2)  = -q(4)/MU0
   !   f(3)  =  q(10)
   !   f(4)  = -q(2)/EPS0
   !   f(5)  =  q(1)/EPS0
   !   f(6)  =  0.0_R8P
   !   f(7)  =  0._R8P
   !   f(8)  =  0._R8P
   !   f(9)  =  0._R8P
   !   f(10) =  ch*ch*q(3)
   !endif
   endsubroutine compute_convective_fluxes_maxwell_div_b

   subroutine compute_convective_fluxes_maxwell_adim_div_b(sir,q,f,chi)
   !< Compute convective fluxes for the adimensional Maxwell system with B cleaning.
   real(R8P), intent(in)    :: sir(3) !< Direction array
   real(R8P), intent(in)    :: q(1:)  !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)  !< Conservative fluxes.
   real(R8P), intent(in)    :: chi    !< Coefficiente velocità trasporto errori divergenza campi
   !$acc routine seq

   f(VAR_DX) =  0.0_R8P         * sir(1) + -q(VAR_BZ)     * sir(2) +  q(VAR_BY)     * sir(3)
   f(VAR_DY) =  q(VAR_BZ)       * sir(1) +  0.0_R8P       * sir(2) + -q(VAR_BX)     * sir(3)
   f(VAR_DZ) = -q(VAR_BY)       * sir(1) +  q(VAR_BX)     * sir(2) +  0.0_R8P       * sir(3)
   f(VAR_BX) =  q(7_I4P)        * sir(1) +  q(VAR_DZ)     * sir(2) + -q(VAR_DY)     * sir(3)
   f(VAR_BY) = -q(VAR_DZ)       * sir(1) +  q(7_I4P)      * sir(2) +  q(VAR_DX)     * sir(3)
   f(VAR_BZ) =  q(VAR_DY)       * sir(1) + -q(VAR_DX)     * sir(2) +  q(7_I4P)      * sir(3)
   f(7_I4P ) =  chi**2*q(VAR_BX)* sir(1) +  chi**2*q(VAR_BY)* sir(2) +  chi**2*q(VAR_BZ)* sir(3)
   endsubroutine compute_convective_fluxes_maxwell_adim_div_b

   subroutine compute_convective_fluxes_maxwell_div_d_b(sir,q,f,chi)
   !< Compute convective fluxes.
   !< J corrected with phi and psi:
   !<```
   !< Fx(1)  =  phi     Fy(1)  = -Bz/muz    Fz(1)  =  By/muy
   !< Fx(2)  =  Bz/muz  Fy(2)  =  phi       Fz(2)  = -Bx/mux
   !< Fx(3)  = -By/muy  Fy(3)  =  Bx/mux    Fz(3)  =  phi
   !< Fx(4)  =  csi     Fy(4)  =  Dz/epsz   Fz(4)  = -Dy/epsy
   !< Fx(5)  = -Dz/epsz Fy(5)  =  csi       Fz(5)  =  Dx/epsx
   !< Fx(6)  =  Dy/epsy Fy(6)  = -Dx/epsx   Fz(6)  =  csi
   !< Fx(7)  =  0       Fy(7)  =  0         Fz(7)  =  0
   !< Fx(8)  =  0       Fy(8)  =  0         Fz(8)  =  0
   !< Fx(9)  =  0       Fy(9)  =  0         Fz(9)  =  0
   !< Fx(10) =  ch^2*Dx Fy(10) =  ch^2*Dy   Fz(10) =  ch^2*Dz
   !< Fx(11) =  ch^2*Bx Fy(10) =  ch^2*By   Fz(10) =  ch^2*Bz
   !<```
   real(R8P), intent(in)    :: sir(3) !< Direction array
   real(R8P), intent(in)    :: q(1:)  !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)  !< Conservative fluxes.
   real(R8P), intent(in)    :: chi    !< Coefficiente velocità trasporto errori divergenza campi
   real(R8P)                :: ch     !< Velocità trasporto errori divergenza campi modello iperbolico
   !$acc routine seq

   ch = chi/sqrt(EPS0*MU0)
   f(VAR_DX) =  q(7_I4P)        * sir(1) + -q(VAR_BZ)/MU0   * sir(2) +  q(VAR_BY)/MU0   * sir(3)
   f(VAR_DY) =  q(VAR_BZ)/MU0   * sir(1) +  q(7_I4P)        * sir(2) + -q(VAR_BX)/MU0   * sir(3)
   f(VAR_DZ) = -q(VAR_BY)/MU0   * sir(1) +  q(VAR_BX)/MU0   * sir(2) +  q(7_I4P)        * sir(3)
   f(VAR_BX) =  q(8_I4P)        * sir(1) +  q(VAR_DZ)/EPS0  * sir(2) + -q(VAR_DY)/EPS0  * sir(3)
   f(VAR_BY) = -q(VAR_DZ)/EPS0  * sir(1) +  q(8_I4P)        * sir(2) +  q(VAR_DX)/EPS0  * sir(3)
   f(VAR_BZ) =  q(VAR_DY)/EPS0  * sir(1) + -q(VAR_DX)/EPS0  * sir(2) +  q(8_I4P)        * sir(3)
   f(7_I4P ) =  ch**2*q(VAR_DX) * sir(1) +  ch**2*q(VAR_DY) * sir(2) +  ch**2*q(VAR_DZ) * sir(3)
   f(8_I4P ) =  ch**2*q(VAR_BX) * sir(1) +  ch**2*q(VAR_BY) * sir(2) +  ch**2*q(VAR_BZ) * sir(3)

   !if (sir(1).eq.1._R8P) then !X
   !   f(1)  =  q(10)
   !   f(2)  =  q(6)/MU0
   !   f(3)  = -q(5)/MU0
   !   f(4)  =  q(11)
   !   f(5)  = -q(3)/EPS0
   !   f(6)  =  q(2)/EPS0
   !   f(7)  =  0._R8P
   !   f(8)  =  0._R8P
   !   f(9)  =  0._R8P
   !   f(10) =  ch*ch*q(1)
   !   f(11) =  ch*ch*q(4)
   !elseif(sir(2).eq.1._R8P) then !Y
   !   f(1)  = -q(6)/MU0
   !   f(2)  =  q(10)
   !   f(3)  =  q(4)/MU0
   !   f(4)  =  q(3)/EPS0
   !   f(5)  =  q(11)
   !   f(6)  = -q(1)/EPS0
   !   f(7)  =  0._R8P
   !   f(8)  =  0._R8P
   !   f(9)  =  0._R8P
   !   f(10) =  ch*ch*q(2)
   !   f(11) =  ch*ch*q(5)
   !elseif(sir(3).eq.1._R8P) then  !Z
   !   f(1)  =  q(5)/MU0
   !   f(2)  = -q(4)/MU0
   !   f(3)  =  q(10)
   !   f(4)  = -q(2)/EPS0
   !   f(5)  =  q(1)/EPS0
   !   f(6)  =  q(11)
   !   f(7)  =  0._R8P
   !   f(8)  =  0._R8P
   !   f(9)  =  0._R8P
   !   f(10) =  ch*ch*q(3)
   !   f(11) =  ch*ch*q(6)
   !endif
   endsubroutine compute_convective_fluxes_maxwell_div_d_b

   subroutine compute_convective_fluxes_maxwell_adim_div_d_b(sir,q,f,chi)
   !< Compute convective fluxes for the adimensional Maxwell system with D/B cleaning.
   real(R8P), intent(in)    :: sir(3) !< Direction array
   real(R8P), intent(in)    :: q(1:)  !< Auxiliary variables.
   real(R8P), intent(inout) :: f(1:)  !< Conservative fluxes.
   real(R8P), intent(in)    :: chi    !< Coefficiente velocità trasporto errori divergenza campi
   !$acc routine seq

   f(VAR_DX) =  q(7_I4P)        * sir(1) + -q(VAR_BZ)     * sir(2) +  q(VAR_BY)     * sir(3)
   f(VAR_DY) =  q(VAR_BZ)       * sir(1) +  q(7_I4P)      * sir(2) + -q(VAR_BX)     * sir(3)
   f(VAR_DZ) = -q(VAR_BY)       * sir(1) +  q(VAR_BX)     * sir(2) +  q(7_I4P)      * sir(3)
   f(VAR_BX) =  q(8_I4P)        * sir(1) +  q(VAR_DZ)     * sir(2) + -q(VAR_DY)     * sir(3)
   f(VAR_BY) = -q(VAR_DZ)       * sir(1) +  q(8_I4P)      * sir(2) +  q(VAR_DX)     * sir(3)
   f(VAR_BZ) =  q(VAR_DY)       * sir(1) + -q(VAR_DX)     * sir(2) +  q(8_I4P)      * sir(3)
   f(7_I4P ) =  chi**2*q(VAR_DX)* sir(1) +  chi**2*q(VAR_DY)* sir(2) +  chi**2*q(VAR_DZ)* sir(3)
   f(8_I4P ) =  chi**2*q(VAR_BX)* sir(1) +  chi**2*q(VAR_BY)* sir(2) +  chi**2*q(VAR_BZ)* sir(3)
   endsubroutine compute_convective_fluxes_maxwell_adim_div_d_b

   subroutine compute_eigenvalues_vector(sir, p)
   !< Compute eigenvalues vector for Maxwell equations in direction sir.
   real(R8P), intent(in)  :: sir(3) !< Direction array
   real(R8P), intent(out) :: p(:)   !< Eigenvalues vector

   p(1) = sir(2) + sir(3) !Dx
   p(2) = sir(1) + sir(3) !Dy
   p(3) = sir(1) + sir(2) !Dz
   p(4) = sir(2) + sir(3) !Bx
   p(5) = sir(1) + sir(3) !By
   p(6) = sir(1) + sir(2) !Bz
   endsubroutine compute_eigenvalues_vector
endmodule adam_prism_riemann_library
