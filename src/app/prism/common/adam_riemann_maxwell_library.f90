!< ADAM, Riemann Problem for Maxwell equations solvers and convective fluxes computations library.
module adam_riemann_maxwell_library
    !< ADAM, Riemann Problem for Maxwell equations solvers and convective fluxes computations library.
    
    use penf, only : I4P, R8P
    use adam_prism_parameters, only : MU0, EPS0 
    implicit none
    private
    public :: compute_riemann_maxwell_llf
    public :: compute_convective_fluxes_Maxwell


! conservative variables are arranged as follows (and they coincide with auxiliary ones):
! q(1): Dx
! q(2): Dy
! q(3): Dz
! q(4): Bx
! q(5): By
! q(6): Bz
! q(4): Jx
! q(5): Jy
! q(6): Jz


! fluxes:
! Fx(1) = 0        Fy(1) = -Bz/muz    Fz(1) = By/muy 
! Fx(2) = Bz/muz   Fy(2) = 0          Fz(2) = -Bx/mux
! Fx(3) = -By/muy  Fy(3) = Bx/mux     Fz(3) = 0
! Fx(4) = 0        Fy(4) = Dz/epsz    Fz(4) = -Dy/epsy
! Fx(5) = -Dz/epsz Fy(5) = 0          Fz(5) = Dx/epsx
! Fx(6) = Dy/epsy  Fy(6) = -Dx/epsx   Fz(6) = 0
! Fx(7) = 0        Fy(7) = 0          Fz(7) = 0
! Fx(8) = 0        Fy(8) = 0          Fz(8) = 0
! Fx(9) = 0        Fy(9) = 0          Fz(9) = 0

    ! ragioniamo nel vuoto e ipotizziamo di essere quindi in un mezzo isotropo e omogeno con costante dielettrica
    ! e permeabilità magnetica pari a mu0 e eps0

contains
    subroutine compute_riemann_maxwell_llf(sir, nv, q1, q4, f, lmax)
        !< Solve the Riemann problem between the state 1 (left) and 4 (right) using the Local-Lax-Friedrichs (LLF, Rusanov) solver.
        !real(R8P),    intent(in)            :: sir(3)        !< Directional (1=x,2=y,3=z) increment.
        integer(I4P), intent(in)            :: nv            !< Number of conservative varibales.
        real(R8P),    intent(in)            :: sir(3)        !< Directional (1=x,2=y,3=z) increment.
        real(R8P),    intent(in)            :: q1(1:)        !< Left state.
        real(R8P),    intent(in)            :: q4(1:)        !< Right state.
        real(R8P),    intent(inout)         :: f(1:)         !< Resulting fluxes.
        real(R8P),    intent(out), optional :: lmax          !< Maximum wave speed estimation.
        real(R8P)                           :: F1(1:nv)      !< State 1 fluxes.
        real(R8P)                           :: F4(1:nv)      !< State 4 fluxes.

        lmax = sqrt(1._R8P/(eps0*mu0))

        call compute_convective_fluxes_Maxwell(sir=sir,q=q1,f=F1)
        call compute_convective_fluxes_Maxwell(sir=sir,q=q4,f=F4)

        f(1) = 0.5_R8P*(F1(1)+F4(1)-lmax*(q4(1)-q1(1)))
        f(2) = 0.5_R8P*(F1(2)+F4(2)-lmax*(q4(2)-q1(2)))
        f(3) = 0.5_R8P*(F1(3)+F4(3)-lmax*(q4(3)-q1(3)))
        f(4) = 0.5_R8P*(F1(4)+F4(4)-lmax*(q4(4)-q1(4)))
        f(5) = 0.5_R8P*(F1(5)+F4(5)-lmax*(q4(5)-q1(5)))
        f(6) = 0.5_R8P*(F1(6)+F4(6)-lmax*(q4(6)-q1(6)))

    endsubroutine compute_riemann_maxwell_llf

    subroutine compute_riemann_maxwell_hll(sir, nv, q1, q4, f, lmax, lmin)

        !< Solve the Riemann problem between the state 1 (left) and 4 (right) using the Harten, Lax and van Leer (HLL) solver.
        !real(R8P),    intent(in)            :: sir(3)        !< Directional (1=x,2=y,3=z) increment.
        integer(I4P), intent(in)            :: nv            !< Number of conservative varibales.
        real(R8P),    intent(in)            :: sir(3)        !< Directional (1=x,2=y,3=z) increment.
        real(R8P),    intent(in)            :: q1(1:)        !< Left state.
        real(R8P),    intent(in)            :: q4(1:)        !< Right state.
        real(R8P),    intent(inout)         :: f(1:)         !< Resulting fluxes.
        real(R8P),    intent(out), optional :: lmax          !< Maximum wave speed estimation.
        real(R8P),    intent(out), optional :: lmin          !< Maximum wave speed estimation.
        real(R8P)                           :: F1(1:nv)      !< State 1 fluxes.
        real(R8P)                           :: F4(1:nv)      !< State 4 fluxes.

        lmax = sqrt(1._R8P/(EPS0*MU0))
        lmin = -lmax

        call compute_convective_fluxes_Maxwell(sir=sir,q=q1,f=F1)
        call compute_convective_fluxes_Maxwell(sir=sir,q=q4,f=F4)

        f(1) = (lmax*F1(1)-lmin*F4(1)+lmax*lmin*(F4(1)-F1(1)))/(lmax-lmin)
        f(2) = (lmax*F1(2)-lmin*F4(2)+lmax*lmin*(F4(2)-F1(2)))/(lmax-lmin)
        f(3) = (lmax*F1(3)-lmin*F4(3)+lmax*lmin*(F4(3)-F1(3)))/(lmax-lmin)
        f(4) = (lmax*F1(4)-lmin*F4(4)+lmax*lmin*(F4(4)-F1(4)))/(lmax-lmin)
        f(5) = (lmax*F1(5)-lmin*F4(5)+lmax*lmin*(F4(5)-F1(5)))/(lmax-lmin)
        f(6) = (lmax*F1(6)-lmin*F4(6)+lmax*lmin*(F4(6)-F1(6)))/(lmax-lmin)


    endsubroutine compute_riemann_maxwell_hll

    

! conservative variables are arranged as follows (and they coincide with auxiliary ones):
! q(1): Dx
! q(2): Dy
! q(3): Dz
! q(4): Bx
! q(5): By
! q(6): Bz

! fluxes:
! Fx(1) = 0        Fy(1) = -Bz/muz    Fz(1) = By/muy 
! Fx(2) = Bz/muz   Fy(2) = 0          Fz(2) = -Bx/mux
! Fx(3) = -By/muy  Fy(3) = Bx/mux     Fz(3) = 0
! Fx(4) = 0        Fy(4) = Dz/epsz    Fz(4) = -Dy/epsy
! Fx(5) = -Dz/epsz Fy(5) = 0          Fz(5) = Dx/epsx
! Fx(6) = Dy/epsy  Fy(6) = -Dx/epsx   Fz(6) = 0


    ! private procedures

    subroutine compute_convective_fluxes_Maxwell(sir,q,f)
    !< Compute convective fluxes for Euler equations from auxiliary variables.
    real(R8P),    intent(in)    :: sir(3)         !< Directional (1=x,2=y,3=z) increment.
    real(R8P),    intent(in)    :: q(1:)         !< Auxiliary variables. 
    real(R8P),    intent(inout) :: f(1:)         !< Conservative fluxes.

    !! CI METTIAMO UN IF SULLA DIREZIONE CONSIDERATA?
    !! DOVE METTIAMO L?INFO SU EPS E MU?? (LE DEVI PURE DICHIARARE)
    !select case(maxloc(si))

    if (sir(1).eq.1._R8P) then !X
    !case(1)

        f(1) =  0.0_R8P
        f(2) =  q(6)/MU0
        f(3) = -q(5)/MU0
        f(4) =  0.0_R8P
        f(5) = -q(3)/EPS0
        f(6) =  q(2)/EPS0
        f(7) = 0._R8P
        f(8) = 0._R8P
        f(9) = 0._R8P

    elseif(sir(2).eq.1._R8P) then !Y
    !case(2)

        f(1) = -q(6)/MU0
        f(2) =  0.0_R8P
        f(3) =  q(4)/MU0
        f(4) =  q(3)/EPS0
        f(5) =  0.0_R8P
        f(6) = -q(1)/EPS0
        f(7) = 0._R8P
        f(8) = 0._R8P
        f(9) = 0._R8P

    elseif(sir(3).eq.1._R8P) then  !Z
    !case(3)

        f(1) =  q(5)/MU0
        f(2) = -q(4)/MU0
        f(3) =  0.0_R8P
        f(4) = -q(2)/EPS0
        f(5) =  q(1)/EPS0  
        f(6) =  0.0_R8P 
        f(7) = 0._R8P
        f(8) = 0._R8P
        f(9) = 0._R8P

    endif
    !endselect

   endsubroutine compute_convective_fluxes_Maxwell

endmodule adam_riemann_maxwell_library