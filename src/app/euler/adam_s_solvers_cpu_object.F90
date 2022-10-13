!< ADAM, Finite Difference/Volume solvers object, CPU backend.
module adam_s_solvers_cpu_object
!< ADAM, Finite Difference/Volume solvers object, CPU backend.

use adam_eos_ic_cpu_object, only : ic_speed_of_sound, ic_total_energy, ic_total_entalpy
use adam_mpih_object,       only : mpih_object
use adam_weno_cpu_object,   only : weno_cpu_object
use finer
use penf
use mpi

implicit none
private
public :: s_solvers_cpu_object

character(len=7),  parameter :: INI_SECTION_NAME='schemes' !< INI (config) file section name containing solvers configs.
character(len=17), parameter :: FLUXES_CONVECTIVE_SCHEMES(2)=['finite-difference', &
                                                              'finite-volume    '] !< List of convective fluxes computation schemes.
character(len=12), parameter :: RIEMANN_SCHEMES(5)=['fvs-llf     ', &
                                                    'fvs-sw      ', &
                                                    'fvs-van-leer', &
                                                    'fvs-glf     ', &
                                                    'llf         '] !< List of Riemann (and FVS) solvers schemes available.

type :: s_solvers_cpu_object
   type(mpih_object)     :: mpih                     !< MPI handler.
   type(weno_cpu_object) :: weno                     !< WENO Kutta solver.
   character(999)        :: fluxes_convective_scheme !< Convective fluxes computation schemes.
   character(999)        :: riemann_solver_scheme    !< Riemann solver scheme for FV or VS solver for FD.
   ! procedures pointer
   procedure(compute_fluxes_convective_int),pass(self),pointer :: compute_fluxes_convective=>null() !< Compute convective fluxes.
   procedure(compute_fvs_int),              pass(self),pointer :: compute_fvs              =>null() !< Compute cells FVS.
   procedure(compute_cell_fvs_global_int),  nopass,    pointer :: compute_cell_fvs_global  =>null() !< Compute cell FVS, global eig.
   procedure(compute_cell_fvs_local_int),   nopass,    pointer :: compute_cell_fvs_local   =>null() !< Compute cell FVS, local  eig.
   procedure(solve_riemann_int),            nopass,    pointer :: solve_riemann            =>null() !< Solve the Riemann problem.
   contains
      ! public methods
      procedure, pass(self) :: description             !< Return pretty-printed object description.
      procedure, pass(self) :: initialize     !< Initialize solvers.
      procedure, pass(self) :: load_from_file !< Load config from file.
endtype s_solvers_cpu_object

interface
   subroutine compute_cell_fvs_local_int(ns, r, rs, u, p, g, fm, fp)
   !< Compute cell fluxes vector splitting, local eigenvalue, interface.
   import :: I4P, R8P
   integer(I4P), intent(in)  :: ns         !< Number of species.
   real(R8P),    intent(in)  :: r          !< Density.
   real(R8P),    intent(in)  :: rs(1:ns)   !< Partial densities.
   real(R8P),    intent(in)  :: u          !< Velocity.
   real(R8P),    intent(in)  :: p          !< Pressure.
   real(R8P),    intent(in)  :: g          !< Specific heats ratio.
   real(R8P),    intent(out) :: fm(1:ns+2) !< Negative (minus direction) fluxes.
   real(R8P),    intent(out) :: fp(1:ns+2) !< Positive (plus direction) fluxes.
   endsubroutine compute_cell_fvs_local_int

   subroutine compute_cell_fvs_global_int(ns, r, rs, u, p, g, lmax, fm, fp)
   !< Compute cell fluxes vector splitting, global eigenvalue, interface.
   import :: I4P, R8P
   integer(I4P), intent(in)  :: ns         !< Number of species.
   real(R8P),    intent(in)  :: r          !< Density.
   real(R8P),    intent(in)  :: rs(1:ns)   !< Partial densities.
   real(R8P),    intent(in)  :: u          !< Velocity.
   real(R8P),    intent(in)  :: p          !< Pressure.
   real(R8P),    intent(in)  :: g          !< Specific heats ratio.
   real(R8P),    intent(in)  :: lmax       !< Maximum local eigenvalue.
   real(R8P),    intent(out) :: fm(1:ns+2) !< Negative (minus direction) fluxes.
   real(R8P),    intent(out) :: fp(1:ns+2) !< Positive (plus direction) fluxes.
   endsubroutine compute_cell_fvs_global_int

   subroutine compute_fluxes_convective_int(self, gc, n, ns, np, cp0, cv0, &
                                            rhos, rho, un, ut1, ut2, g, p, &
                                            f_rho, f_rho_un, f_rho_ut1, f_rho_ut2, f_rho_e)
   !< Compute convective fluxes on a coordinate direction, interface.
   import :: s_solvers_cpu_object, I4P, R8P
   class(s_solvers_cpu_object), intent(in)    :: self                 !< Solvers.
   integer(I4P),                intent(in)    :: gc                   !< Number of ghost cells used.
   integer(I4P),                intent(in)    :: n                    !< Number of cells.
   integer(I4P),                intent(in)    :: ns                   !< Number of species.
   integer(I4P),                intent(in)    :: np                   !< Number of 1D primitive varibales.
   real(R8P),                   intent(in)    :: cp0(1:ns), cv0(1:ns) !< Specific heats of initial species.
   real(R8P),                   intent(in)    :: rhos(1:,1-gc:)       !< Partial densities       [1:ns,1-gc:n+gc].
   real(R8P),                   intent(in)    ::     rho(1-gc:)       !< Density                      [1-gc:n+gc].
   real(R8P),                   intent(in)    ::      un(1-gc:)       !< Normal velocity              [1-gc:n+gc].
   real(R8P),                   intent(in)    ::     ut1(1-gc:)       !< Tangential velocity 1        [1-gc:n+gc].
   real(R8P),                   intent(in)    ::     ut2(1-gc:)       !< Tangential velocity 2        [1-gc:n+gc].
   real(R8P),                   intent(in)    ::       g(1-gc:)       !< Specific heats ratio         [1-gc:n+gc].
   real(R8P),                   intent(in)    ::       p(1-gc:)       !< Pressure                     [1-gc:n+gc].
   real(R8P),                   intent(inout) :: f_rho(1:, 0:)        !< Flux of mass                 [0:n,1:ns].
   real(R8P),                   intent(inout) ::  f_rho_un(0:)        !< Flux normal momentums        [0:n].
   real(R8P),                   intent(inout) :: f_rho_ut1(0:)        !< Flux of tangential1 momentum [0:n].
   real(R8P),                   intent(inout) :: f_rho_ut2(0:)        !< Flux of tangential2 momentum [0:n].
   real(R8P),                   intent(inout) ::   f_rho_e(0:)        !< Flux energy                  [0:n].
   endsubroutine compute_fluxes_convective_int

   subroutine compute_fvs_int(self, ns, gc, n, r, rs, u, p, g, fm, fp)
   !< Compute fluxes vector splitting.
   import :: s_solvers_cpu_object, I4P, R8P
   class(s_solvers_cpu_object), intent(in)  :: self         !< Solvers.
   integer(I4P),                intent(in)  :: ns           !< Number of species.
   integer(I4P),                intent(in)  :: gc           !< Number of ghost cells used.
   integer(I4P),                intent(in)  :: n            !< Number of cells.
   real(R8P),                   intent(in)  :: r(1-gc:)     !< Density.
   real(R8P),                   intent(in)  :: rs(1:,1-gc:) !< Partial densities.
   real(R8P),                   intent(in)  :: u(1-gc:)     !< Velocity.
   real(R8P),                   intent(in)  :: p(1-gc:)     !< Pressure.
   real(R8P),                   intent(in)  :: g(1-gc:)     !< Specific heats ratio.
   real(R8P),                   intent(out) :: fm(1:,1-gc:) !< Negative (minus direction) fluxes.
   real(R8P),                   intent(out) :: fp(1:,1-gc:) !< Positive (plus direction) fluxes.
   endsubroutine compute_fvs_int

   subroutine solve_riemann_int(r1, u1, p1, g1, r4, u4, p4, g4, F)
   !< Solve the Riemann problem between the state $1$ and $4$, interface.
   import :: R8P
   real(R8P), intent(in)  :: r1      !< Density of state 1.
   real(R8P), intent(in)  :: u1      !< Velocity of state 1.
   real(R8P), intent(in)  :: p1      !< Pressure of state 1.
   real(R8P), intent(in)  :: g1      !< Specific heats ratio of state 1.
   real(R8P), intent(in)  :: r4      !< Density of state 4.
   real(R8P), intent(in)  :: u4      !< Velocity of state 4.
   real(R8P), intent(in)  :: p4      !< Pressure of state 4.
   real(R8P), intent(in)  :: g4      !< Specific heats ratio of state 4.
   real(R8P), intent(out) :: F(1:3)  !< Resulting fluxes.
   endsubroutine solve_riemann_int
endinterface

contains
   pure function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(s_solvers_cpu_object), intent(in) :: self             !< Solvers.
   character(len=:), allocatable           :: desc             !< Description.
   character(len=1), parameter             :: NL=new_line('a') !< New line character.

   desc =       self%mpih%myrankstr//'FD/FV schemes'//NL
   desc = desc//self%mpih%myrankstr//'  convective fluxes compuation scheme: '//trim(self%fluxes_convective_scheme)//NL
   desc = desc//self%mpih%myrankstr//'  Riemann/FVS solver:                  '//trim(self%riemann_solver_scheme   )
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize the equation.
   class(s_solvers_cpu_object), intent(inout) :: self            !< Solvers.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   call self%mpih%initialize(do_mpi_init=.false.)
   print '(A)', self%mpih%myrankstr//'s_solvers_cpu_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)

   ! set pointer procedures
   select case(trim(adjustl(self%fluxes_convective_scheme)))
   case('finite-difference')
      self%compute_fluxes_convective => compute_fluxes_convective_fd
   case('finite-volume')
      self%compute_fluxes_convective => compute_fluxes_convective_fv
   endselect

   select case(trim(adjustl(self%riemann_solver_scheme)))
   case('fvs-llf')
      self%compute_cell_fvs_local => compute_cell_fvs_llf
      self%compute_fvs            => compute_fvs_local
   case('fvs-sw')
      self%compute_cell_fvs_local => compute_cell_fvs_sw
      self%compute_fvs            => compute_fvs_local
   case('fvs-van-leer')
      self%compute_cell_fvs_local => compute_cell_fvs_van_leer
      self%compute_fvs            => compute_fvs_local
   case('fvs-glf')
      self%compute_cell_fvs_global => compute_cell_fvs_glf
      self%compute_fvs            => compute_fvs_global
   case('llf')
      self%solve_riemann => solve_riemann_llf
   endselect

   call self%weno%initialize(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', self%mpih%myrankstr//'s_solvers_cpu_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters, go_on_fail)
   !< Load config from file.
   class(s_solvers_cpu_object), intent(inout)        :: self            !< Solvers.
   type(file_ini),              intent(in)           :: file_parameters !< Simulation parameters ini file handler.
   logical,                     intent(in), optional :: go_on_fail      !< Go on if load fails.
   logical                                           :: go_on_fail_     !< Go on if load fails.
   integer(I4P)                                      :: s               !< Counter.
   integer(I4P)                                      :: error           !< Error status.
   logical                                           :: known_scheme    !< Flag to check if the scheme is known.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail
   call file_parameters%get(section_name=INI_SECTION_NAME,option_name='riemann_solver',val=self%riemann_solver_scheme,error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(riemann_solver)')
   call file_parameters%get(section_name=INI_SECTION_NAME,option_name='fluxes_convective',val=self%fluxes_convective_scheme, &
                            error=error)
   if (.not.go_on_fail_.and.error>0) call self%mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(fluxes_convective)')

   ! check if the selected schemes are known (i.e. available)
   known_scheme = .false.
   do s=1, size(FLUXES_CONVECTIVE_SCHEMES, dim=1)
      if (trim(adjustl(self%fluxes_convective_scheme))==FLUXES_CONVECTIVE_SCHEMES(s)) then
         known_scheme = .true.
         exit
      endif
   enddo
   if (.not.known_scheme) &
   call self%mpih%error_stop(msg=': convective fluxes computation scheme "'//trim(self%fluxes_convective_scheme)//'" is unknown!')

   known_scheme = .false.
   do s=1, size(RIEMANN_SCHEMES, dim=1)
      if (trim(adjustl(self%riemann_solver_scheme))==RIEMANN_SCHEMES(s)) then
         known_scheme = .true.
         exit
      endif
   enddo
   if (.not.known_scheme) &
   call self%mpih%error_stop(msg=': Riemann solver scheme "'//trim(self%riemann_solver_scheme)//'" is unknown!')
   endsubroutine load_from_file

   ! finite difference solvers
   subroutine compute_fluxes_convective_fd(self, gc, n, ns, np, cp0, cv0, &
                                           rhos, rho, un, ut1, ut2, g, p, &
                                           f_rho, f_rho_un, f_rho_ut1, f_rho_ut2, f_rho_e)
   !< Compute convective fluxes on coordinate direction by means of Flux Vector Splitting (FVS), finite difference scheme.
   class(s_solvers_cpu_object), intent(in)    :: self                 !< Solvers.
   integer(I4P),                intent(in)    :: gc                   !< Number of ghost cells used.
   integer(I4P),                intent(in)    :: n                    !< Number of cells.
   integer(I4P),                intent(in)    :: ns                   !< Number of species.
   integer(I4P),                intent(in)    :: np                   !< Number of 1D primitive varibales.
   real(R8P),                   intent(in)    :: cp0(1:ns), cv0(1:ns) !< Specific heats of initial species.
   real(R8P),                   intent(in)    :: rhos(1:,1-gc:)       !< Partial densities       [1:ns,1-gc:n+gc].
   real(R8P),                   intent(in)    ::     rho(1-gc:)       !< Density                      [1-gc:n+gc].
   real(R8P),                   intent(in)    ::      un(1-gc:)       !< Normal velocity              [1-gc:n+gc].
   real(R8P),                   intent(in)    ::     ut1(1-gc:)       !< Tangential velocity 1        [1-gc:n+gc].
   real(R8P),                   intent(in)    ::     ut2(1-gc:)       !< Tangential velocity 2        [1-gc:n+gc].
   real(R8P),                   intent(in)    ::       g(1-gc:)       !< Specific heats ratio         [1-gc:n+gc].
   real(R8P),                   intent(in)    ::       p(1-gc:)       !< Pressure                     [1-gc:n+gc].
   real(R8P),                   intent(inout) :: f_rho(1:, 0:)        !< Flux of mass                 [0:n,1:ns].
   real(R8P),                   intent(inout) ::  f_rho_un(0:)        !< Flux normal momentums        [0:n].
   real(R8P),                   intent(inout) :: f_rho_ut1(0:)        !< Flux of tangential1 momentum [0:n].
   real(R8P),                   intent(inout) :: f_rho_ut2(0:)        !< Flux of tangential2 momentum [0:n].
   real(R8P),                   intent(inout) ::   f_rho_e(0:)        !< Flux energy                  [0:n].
   real(R8P)                                  :: q (1:np,  1-gc:n+gc) !< 1D primitive variables.
   real(R8P)                                  :: fm(1:ns+2,1-gc:n+gc) !< Negative (minus direction) fluxes.
   real(R8P)                                  :: fp(1:ns+2,1-gc:n+gc) !< Positive (plus direction) fluxes.
   real(R8P)                                  :: fr(1:ns+2,0:n)       !< Reconstructed fluxes.
   real(R8P)                                  :: fluxes(3)            !< 1D fluxes.
   integer(I4P)                               :: i                    !< Counter.

   ! assembly 1D primitive variables array
   do i=1-gc, n+gc
      q(:, i) = [rhos(:,i), un(i), p(i), rho(i), g(i)]
   enddo

   ! compute fluxes vector splitting at cell center
   call self%compute_fvs(ns=ns, gc=gc, n=n, rs=rhos, r=rho, u=un, p=p, g=g, fm=fm, fp=fp)

   ! reconstruct fluxes at intefaces by WENO
   call weno_reconstruct_fvs(weno=self%weno, gc=gc, n=n, ns=ns, np=np, q=q, fp=fp, fm=fm, fr=fr)

   ! assemble 3D fluxed from 1D ones
   do i=0, n
      fluxes(1) = sum(fr(1:ns, i))
      fluxes(2) =     fr(ns+1, i)
      fluxes(3) =     fr(ns+2, i)

      f_rho(:, i) = fr(1:ns, i)
      f_rho_un(i) = fluxes(2)
      if (fluxes(1)>0._R8P) then
         f_rho_ut1(i) = fluxes(1) * ut1(i)
         f_rho_ut2(i) = fluxes(1) * ut2(i)
           f_rho_e(i) = fluxes(3) + 0.5_R8P * fluxes(1) * (ut1(i)**2 + ut2(i)**2)
      else
         f_rho_ut1(i) = fluxes(1) * ut1(i+1)
         f_rho_ut2(i) = fluxes(1) * ut2(i+1)
           f_rho_e(i) = fluxes(3) + 0.5_R8P * fluxes(1) * (ut1(i+1)**2 + ut2(i+1)**2)
      endif
   enddo
   endsubroutine compute_fluxes_convective_fd

   ! local FVS solvers
   subroutine compute_cell_fvs_llf(ns, r, rs, u, p, g, fm, fp)
   !< Compute cell fluxes vector splitting by means of (local) Lax Friedrichs (Rusanov) solver.
   integer(I4P), intent(in)  :: ns         !< Number of species.
   real(R8P),    intent(in)  :: r          !< Density.
   real(R8P),    intent(in)  :: rs(1:ns)   !< Partial densities.
   real(R8P),    intent(in)  :: u          !< Velocity.
   real(R8P),    intent(in)  :: p          !< Pressure.
   real(R8P),    intent(in)  :: g          !< Specific heats ratio.
   real(R8P),    intent(out) :: fm(1:ns+2) !< Negative (minus direction) fluxes.
   real(R8P),    intent(out) :: fp(1:ns+2) !< Positive (plus direction) fluxes.

   call compute_cell_fvs_glf(ns=ns, r=r, rs=rs, u=u, p=p, g=g, lmax=abs(u) + sqrt(g * p / r), fm=fm, fp=fp)
   endsubroutine compute_cell_fvs_llf

   subroutine compute_cell_fvs_sw(ns, r, rs, u, p, g, fm, fp)
   !< Compute cell fluxes vector splitting by means of Steger/Warming solver.
   integer(I4P), intent(in)  :: ns         !< Number of species.
   real(R8P),    intent(in)  :: r          !< Density.
   real(R8P),    intent(in)  :: rs(1:ns)   !< Partial densities.
   real(R8P),    intent(in)  :: u          !< Velocity.
   real(R8P),    intent(in)  :: p          !< Pressure.
   real(R8P),    intent(in)  :: g          !< Specific heats ratio.
   real(R8P),    intent(out) :: fm(1:ns+2) !< Negative (minus direction) fluxes.
   real(R8P),    intent(out) :: fp(1:ns+2) !< Positive (plus direction) fluxes.
   real(R8P)                 :: ss         !< Speed of sound.
   real(R8P)                 :: l(3)       !< Eigenvalues.
   real(R8P)                 :: lmp(3,2)   !< Eigenvalues decomposition.
   real(R8P)                 :: r_2g       !< `r/(2g)`.
   real(R8P)                 :: H0         !< Total entalpy.

   H0 = ic_total_entalpy(g=g, density=r, pressure=p, velocity_sq_norm=u*u)
   r_2g = 0.5_R8P * r / g
   ss = sqrt(g * p / r)
   l = [u-ss, u, u+ss]
   lmp(:,1) = 0.5_R8P * (l - abs(l))
   lmp(:,2) = 0.5_R8P * (l + abs(l))

   fm(1:ns) = r_2g * (                lmp(1,1) + 2.0_R8P * (g - 1.0_R8P)         * lmp(2,1) +                 lmp(3,1)) * rs / r
   fm(ns+1) = r_2g * ((u - ss)      * lmp(1,1) + 2.0_R8P * (g - 1.0_R8P) * u     * lmp(2,1) + (u + ss)      * lmp(3,1))
   fm(ns+2) = r_2g * ((H0 - u * ss) * lmp(1,1) +           (g - 1.0_R8P) * u * u * lmp(2,1) + (H0 + u * ss) * lmp(3,1))

   fp(1:ns) = r_2g * (                lmp(1,2) + 2.0_R8P * (g - 1.0_R8P)         * lmp(2,2) +                 lmp(3,2)) * rs / r
   fp(ns+1) = r_2g * ((u - ss)      * lmp(1,2) + 2.0_R8P * (g - 1.0_R8P) * u     * lmp(2,2) + (u + ss)      * lmp(3,2))
   fp(ns+2) = r_2g * ((H0 - u * ss) * lmp(1,2) +           (g - 1.0_R8P) * u * u * lmp(2,2) + (H0 + u * ss) * lmp(3,2))
   endsubroutine compute_cell_fvs_sw

   subroutine compute_cell_fvs_van_leer(ns, r, rs, u, p, g, fm, fp)
   !< Compute cell fluxes vector splitting by means of van Leer method.
   integer(I4P), intent(in)  :: ns         !< Number of species.
   real(R8P),    intent(in)  :: r          !< Density.
   real(R8P),    intent(in)  :: rs(1:ns)   !< Partial densities.
   real(R8P),    intent(in)  :: u          !< Velocity.
   real(R8P),    intent(in)  :: p          !< Pressure.
   real(R8P),    intent(in)  :: g          !< Specific heats ratio.
   real(R8P),    intent(out) :: fm(1:ns+2) !< Negative (minus direction) fluxes.
   real(R8P),    intent(out) :: fp(1:ns+2) !< Positive (plus direction) fluxes.
   real(R8P)                 :: frm, frp   !< Negative (minus direction) mass fluxes.
   real(R8P)                 :: a          !< Speed of sound.
   real(R8P)                 :: M          !< Mach number.
   real(R8P)                 :: gm1_2      !< `(g-1)/2`.
   real(R8P)                 :: ra_4       !< `rho*a/4`.
   real(R8P)                 :: c1         !< `2*a/g`.
   real(R8P)                 :: c2         !< `2*a**2/(g**2-1)`.

   gm1_2 = (g - 1._R8P) * 0.5_R8P
   a = sqrt(g * p / r)
   M = u / a
   ra_4 = r * a * 0.25_R8P
   c1 = 2._R8P * a / g
   c2 = 2._R8P * a**2 / (g**2 - 1._R8P)

   frm      = - ra_4 * ((1 - M)**2)                ; frp      = ra_4 * ((1 + M)**2)
   fm(1:ns) = frm * rs(1:ns)/r                     ; fp(1:ns) = frp * rs(1:ns)/r
   fm(2)    = frm * (c1 * (gm1_2 * M - 1._R8P))    ; fp(2)    = frp * (c1 * (gm1_2 * M + 1._R8P))
   fm(3)    = frm * (c2 * (gm1_2 * M - 1._R8P)**2) ; fp(3)    = frp * (c2 * (gm1_2 * M + 1._R8P)**2)
   endsubroutine compute_cell_fvs_van_leer

   subroutine compute_fvs_local(self, ns, gc, n, r, rs, u, p, g, fm, fp)
   !< Compute fluxes vector splitting, local eigenvalue estimation.
   class(s_solvers_cpu_object), intent(in)  :: self         !< Solvers.
   integer(I4P),                intent(in)  :: ns           !< Number of species.
   integer(I4P),                intent(in)  :: gc           !< Number of ghost cells used.
   integer(I4P),                intent(in)  :: n            !< Number of cells.
   real(R8P),                   intent(in)  :: r(1-gc:)     !< Density.
   real(R8P),                   intent(in)  :: rs(1:,1-gc:) !< Partial densities.
   real(R8P),                   intent(in)  :: u(1-gc:)     !< Velocity.
   real(R8P),                   intent(in)  :: p(1-gc:)     !< Pressure.
   real(R8P),                   intent(in)  :: g(1-gc:)     !< Specific heats ratio.
   real(R8P),                   intent(out) :: fm(1:,1-gc:) !< Negative (minus direction) fluxes.
   real(R8P),                   intent(out) :: fp(1:,1-gc:) !< Positive (plus direction) fluxes.
   integer(I4P)                             :: i            !< Counter.

   do i=1-gc, n+gc
      call self%compute_cell_fvs_local(ns=ns, rs=rs(:,i), r=r(i), u=u(i), p=p(i), g=g(i), fm=fm(:,i), fp=fp(:,i))
   enddo
   endsubroutine compute_fvs_local

   ! global FVS solvers
   subroutine compute_cell_fvs_glf(ns, r, rs, u, p, g, lmax, fm, fp)
   !< Compute cell fluxes vector splitting by means of global Lax Friedrichs solver.
   !<
   !< @Note The maximum eigenvalue is evalutate outside and it should have a *global* nature. If it is
   !< evaluated only into the current cell this becomes the local Lax Friedrichs (Rusanov) solver.
   integer(I4P), intent(in)  :: ns         !< Number of species.
   real(R8P),    intent(in)  :: r          !< Density.
   real(R8P),    intent(in)  :: rs(1:ns)   !< Partial densities.
   real(R8P),    intent(in)  :: u          !< Velocity.
   real(R8P),    intent(in)  :: p          !< Pressure.
   real(R8P),    intent(in)  :: g          !< Specific heats ratio.
   real(R8P),    intent(in)  :: lmax       !< Maximum global eigenvalue.
   real(R8P),    intent(out) :: fm(1:ns+2) !< Negative (minus direction) fluxes.
   real(R8P),    intent(out) :: fp(1:ns+2) !< Positive (plus direction) fluxes.
   real(R8P)                 :: Q(1:ns+2)  !< Conservative variables.
   real(R8P)                 :: F(0:ns+2)  !< Conservative fluxes.

   Q = [rs, r * u, r * ic_total_energy(gm1=g-1._R8P, density=r, pressure=p, velocity_sq_norm=u*u)]
   F(0) = r * u
   F(1:ns) = rs * u
   F(ns+1) = F(0) * u + p
   F(ns+2) = F(0) * ic_total_entalpy(g=g, density=r, pressure=p, velocity_sq_norm=u*u)
   fm = 0.5_R8P * (F(1:ns+2) - lmax * Q)
   fp = 0.5_R8P * (F(1:ns+2) + lmax * Q)
   endsubroutine compute_cell_fvs_glf

   subroutine compute_fvs_global(self, ns, gc, n, r, rs, u, p, g, fm, fp)
   !< Compute fluxes vector splitting, global eigenvalue estimation.
   class(s_solvers_cpu_object), intent(in)  :: self         !< Solvers.
   integer(I4P),                intent(in)  :: ns           !< Number of species.
   integer(I4P),                intent(in)  :: gc           !< Number of ghost cells used.
   integer(I4P),                intent(in)  :: n            !< Number of cells.
   real(R8P),                   intent(in)  :: r(1-gc:)     !< Density.
   real(R8P),                   intent(in)  :: rs(1:,1-gc:) !< Partial densities.
   real(R8P),                   intent(in)  :: u(1-gc:)     !< Velocity.
   real(R8P),                   intent(in)  :: p(1-gc:)     !< Pressure.
   real(R8P),                   intent(in)  :: g(1-gc:)     !< Specific heats ratio.
   real(R8P),                   intent(out) :: fm(1:,1-gc:) !< Negative (minus direction) fluxes.
   real(R8P),                   intent(out) :: fp(1:,1-gc:) !< Positive (plus direction) fluxes.
   real(R8P)                                :: lmax         !< Maximum local (to the stencils) eigenvalues estimation.
   integer(I4P)                             :: i            !< Counter.

   ! estimate a maxium global eigenvalue
   lmax = 0._R8P
   do i=1-gc, n+gc
      lmax = max(lmax, abs(u(i))+sqrt(g(i) * p(i) / r(i)))
   enddo
   ! compute fluxes splitting
   do i=1-gc, n+gc
      call self%compute_cell_fvs_global(ns=ns, rs=rs(:,i), r=r(i), u=u(i), p=p(i), g=g(i), lmax=lmax, fm=fm(:,i), fp=fp(:,i))
   enddo
   endsubroutine compute_fvs_global

   ! WENO FD
   subroutine weno_reconstruct_fvs(weno, gc, n, ns, np, q, fp, fm, fr)
   !< Reconstruct FVS fluxes in pseudo characteristic variables by WENO methods.
   type(weno_cpu_object), intent(in)  :: weno                     !< WENO solver.
   integer(I4P),          intent(in)  :: gc                       !< Number of ghost cells used.
   integer(I4P),          intent(in)  :: n                        !< Number of cells.
   integer(I4P),          intent(in)  :: ns                       !< Number of species.
   integer(I4P),          intent(in)  :: np                       !< Number of 1D primitive varibales.
   real(R8P),             intent(in)  :: q (1:np,  1-gc:n+gc)     !< 1D primitive variables.
   real(R8P),             intent(in)  :: fm(1:ns+2,1-gc:n+gc)     !< Negative (minus direction) fluxes.
   real(R8P),             intent(in)  :: fp(1:ns+2,1-gc:n+gc)     !< Positive (plus direction) fluxes.
   real(R8P),             intent(out) :: fr(1:ns+2,0:n)           !< Reconstructed fluxes.
   real(R8P)                          :: qa(1:np,1:2)             !< Roe averaged primitive variables.
   real(R8P)                          :: Lqm(1:ns+2,1:ns+2,1:2)   !< Left eigenvalues matrix of mean primitive var.
   real(R8P)                          :: Rqm(1:ns+2,1:ns+2,1:2)   !< Right eigenvalues matrix of mean primitive var.
   real(R8P)                          :: c(1:ns+2,1:2,1-gc:-1+gc) !< Pseudo characteristic fluxes.
   real(R8P)                          :: cr(1:ns+2,1:2)           !< Pseudo characteristic reconstructed fluxes.
   real(R8P)                          :: ffr(1:ns+2,1:2,0:n+1)    !< Reconstructed fluxes.
   integer(I4P)                       :: i, j, v, f               !< Counter.

   select case(weno%weno_s)
   case(1_I4P)
      ! first order, no reconstruction
      do i=0, n
         fr(:,i) = fp(:,i) + fm(:,i+1)
      enddo
   case(2_I4P,3_I4P,4_I4P)
      ! compute WENO reconstruction
      do i=0, n+1
         ! compute L/R eigenvectors approximation at interfaces
         do f=1, 2
            if (i==0  .and.f==1) cycle
            if (i==n+1.and.f==2) cycle
            qa(:,f) = roe_averages(ns=ns, np=np, ql=q(:,i+f-2), qr=q(:,i+f-1))
            Lqm(:,:,f) =  left_eigenvectors(ns=ns, np=np, q=qa(:,f))
            Rqm(:,:,f) = right_eigenvectors(ns=ns, np=np, q=qa(:,f))
         enddo

         ! compute pseudo characteristic fluxes
         do j=i+1-weno%weno_s, i-1+weno%weno_s
            do f=1, 2
               if (i==0  .and.f==1) cycle
               if (i==n+1.and.f==2) cycle
               if (f==1) then ! left interface, reconstruct negative fluxes
                  do v=1, ns+2
                     c(v,f,j-i) = dot_product(Lqm(v,1:ns+2,f), fm(1:ns+2,j))
                  enddo
               else           ! right interface, reconstruct positive fluxes
                  do v=1, ns+2
                     c(v,f,j-i) = dot_product(Lqm(v,1:ns+2,f), fp(1:ns+2,j))
                  enddo
               endif
            enddo
         enddo

         ! compute WENO reconstruction of pseudo charteristic fluxes
         do v=1, ns+2
            cr(v,:) = weno%reconstructed(s=weno%weno_s, v=c(v,:,1-weno%weno_s:-1+weno%weno_s))
         enddo

         ! trasform back reconstructed pseudo charteristic fluxes to primitive ones
         do f=1, 2
            if (i==0  .and.f==1) cycle
            if (i==n+1.and.f==2) cycle
            do v=1, ns+2
               ffr(v,f,i) = dot_product(Rqm(v,1:ns+2,f), cr(1:ns+2,f))
            enddo
         enddo
      enddo

      ! assemble fuxes at interfaces
      do i=0, n
         fr(:,i) = ffr(:,2,i) + ffr(:,1,i+1)
      enddo
   endselect
   endsubroutine weno_reconstruct_fvs

   ! finite volume solvers
   subroutine compute_fluxes_convective_fv(self, gc, n, ns, np, cp0, cv0, &
                                           rhos, rho, un, ut1, ut2, g, p, &
                                           f_rho, f_rho_un, f_rho_ut1, f_rho_ut2, f_rho_e)
   !< Compute convective fluxes on coordinate direction, finite volume scheme.
   class(s_solvers_cpu_object), intent(in)    :: self                 !< Solvers.
   integer(I4P),                intent(in)    :: gc                   !< Number of ghost cells used.
   integer(I4P),                intent(in)    :: n                    !< Number of cells.
   integer(I4P),                intent(in)    :: ns                   !< Number of species.
   integer(I4P),                intent(in)    :: np                   !< Number of 1D primitive varibales.
   real(R8P),                   intent(in)    :: cp0(1:ns), cv0(1:ns) !< Specific heats of initial species.
   real(R8P),                   intent(in)    :: rhos(1:,1-gc:)       !< Partial densities       [1:ns,1-gc:n+gc].
   real(R8P),                   intent(in)    ::     rho(1-gc:)       !< Density                      [1-gc:n+gc].
   real(R8P),                   intent(in)    ::      un(1-gc:)       !< Normal velocity              [1-gc:n+gc].
   real(R8P),                   intent(in)    ::     ut1(1-gc:)       !< Tangential velocity 1        [1-gc:n+gc].
   real(R8P),                   intent(in)    ::     ut2(1-gc:)       !< Tangential velocity 2        [1-gc:n+gc].
   real(R8P),                   intent(in)    ::       g(1-gc:)       !< Specific heats ratio         [1-gc:n+gc].
   real(R8P),                   intent(in)    ::       p(1-gc:)       !< Pressure                     [1-gc:n+gc].
   real(R8P),                   intent(inout) :: f_rho(1:, 0:)        !< Flux of mass                 [0:n,1:ns].
   real(R8P),                   intent(inout) ::  f_rho_un(0:)        !< Flux normal momentums        [0:n].
   real(R8P),                   intent(inout) :: f_rho_ut1(0:)        !< Flux of tangential1 momentum [0:n].
   real(R8P),                   intent(inout) :: f_rho_ut2(0:)        !< Flux of tangential2 momentum [0:n].
   real(R8P),                   intent(inout) ::   f_rho_e(0:)        !< Flux energy                  [0:n].
   real(R8P)                                  ::   fluxes(3)          !< 1D fluxes.
   real(R8P)                                  :: q (1:np, 1-gc:n+gc)  !< 1D primitive variables.
   real(R8P)                                  :: qr(1:np,1:2,0:n+1)   !< Reconstructed 1D primitive variables.
   integer(I4P)                               :: i                    !< Counter.

   ! assembly 1D primitive variables array
   do i=1-gc, n+gc
      q(:, i) = [rhos(:,i), un(i), p(i), rho(i), g(i)]
   enddo

   call weno_reconstruct_q(weno=self%weno, gc=gc, n=n, ns=ns, np=np, cp0=cp0, cv0=cv0, q=q, qr=qr)

   do i=0, n
      ! computing normal fluxes solving Riemann problem
      call self%solve_riemann(r1=qr(ns+3,2,i  ), u1=qr(ns+1,2,i  ), p1=qr(ns+2,2,i  ), g1=qr(ns+4,2,i  ), &
                              r4=qr(ns+3,1,i+1), u4=qr(ns+1,1,i+1), p4=qr(ns+2,1,i+1), g4=qr(ns+4,1,i+1), F=fluxes)
      if (fluxes(1)>0._R8P) then
           f_rho(:,i) = fluxes(1) * qr(1:ns,2,i) / qr(ns+3,2,i)
          f_rho_un(i) = fluxes(2)
         f_rho_ut1(i) = fluxes(1) * ut1(i)
         f_rho_ut2(i) = fluxes(1) * ut2(i)
           f_rho_e(i) = fluxes(3) + 0.5_R8P * fluxes(1) * (ut1(i)**2 + ut2(i)**2)
      else
           f_rho(:,i) = fluxes(1) * qr(1:ns,1,i+1) / qr(ns+3,1,i+1)
          f_rho_un(i) = fluxes(2)
         f_rho_ut1(i) = fluxes(1) * ut1(i+1)
         f_rho_ut2(i) = fluxes(1) * ut2(i+1)
           f_rho_e(i) = fluxes(3) + 0.5_R8P * fluxes(1) * (ut1(i+1)**2 + ut2(i+1)**2)
      endif
   enddo
   endsubroutine compute_fluxes_convective_fv

   ! Riemann solvers
   subroutine solve_riemann_llf(r1, u1, p1, g1, r4, u4, p4, g4, F)
   !< Solve the Riemann problem between the state $1$ and $4$ using the (local) Lax Friedrichs (Rusanov) solver.
   real(R8P), intent(in)  :: r1      !< Density of state 1.
   real(R8P), intent(in)  :: u1      !< Velocity of state 1.
   real(R8P), intent(in)  :: p1      !< Pressure of state 1.
   real(R8P), intent(in)  :: g1      !< Specific heats ratio of state 1.
   real(R8P), intent(in)  :: r4      !< Density of state 4.
   real(R8P), intent(in)  :: u4      !< Velocity of state 4.
   real(R8P), intent(in)  :: p4      !< Pressure of state 4.
   real(R8P), intent(in)  :: g4      !< Specific heats ratio of state 4.
   real(R8P), intent(out) :: F(1:3)  !< Resulting fluxes.
   real(R8P)              :: lmax    !< Maximum wave speed estimation.
   real(R8P)              :: F1(1:3) !< State 1 fluxes.
   real(R8P)              :: F4(1:3) !< State 4 fluxes.
   real(R8P)              :: u       !< Velocity of the intermediate states.
   real(R8P)              :: p       !< Pressure of the intermediate states.
   real(R8P)              :: E1, E4  !< Total energy of state 1 and 4.
   real(R8P)              :: S1, S4  !< Maximum wave speed of state 1 and 4.

   ! evaluating the intermediates states 2 and 3 from the known states U1,U4 using the PVRS approximation
   call compute_inter_states
   ! evalutaing the maximum waves speed
   lmax = max(abs(S1), abs(u), abs(S4))
   ! computing the fluxes of state 1 and 4
   F1 = fluxes(p = p1, r = r1, u = u1, g = g1)
   F4 = fluxes(p = p4, r = r4, u = u4, g = g4)
   E1 = ic_total_energy(gm1=g1-1._R8P, density=r1, pressure=p1, velocity_sq_norm=u1*u1)
   E4 = ic_total_energy(gm1=g4-1._R8P, density=r4, pressure=p4, velocity_sq_norm=u4*u4)
   ! computing the Lax-Friedrichs fluxes approximation
   F(1) = 0.5_R8P*(F1(1) + F4(1) - lmax*(r4    - r1   ))
   F(2) = 0.5_R8P*(F1(2) + F4(2) - lmax*(r4*u4 - r1*u1))
   F(3) = 0.5_R8P*(F1(3) + F4(3) - lmax*(r4*E4 - r1*E1))
   contains
      pure function fluxes(p, r, u, g) result(Fc)
      !< 1D Euler fluxes from primitive variables.
      real(R8P), intent(in) :: p       !< Pressure.
      real(R8P), intent(in) :: r       !< Density.
      real(R8P), intent(in) :: u       !< Velocity.
      real(R8P), intent(in) :: g       !< Specific heats ratio.
      real(R8P)             :: Fc(1:3) !< State fluxes.

      Fc(1) = r*u
      Fc(2) = Fc(1)*u + p
      Fc(3) = Fc(1)*ic_total_entalpy(g=g, density=r, pressure=p, velocity_sq_norm=u*u)
      ! Fc(3) = Fc(1)*H(p=p, r=r, u=u, g=g)
      endfunction fluxes

      subroutine compute_inter_states
      !< Compute inter states (23*-states) from state1 and state4.
      real(R8P)            :: a1             !< Speed of sound of state 1.
      real(R8P)            :: a4             !< Speed of sound of state 4.
      real(R8P)            :: ram            !< Mean value of rho*a.
      real(R8P), parameter :: toll=1e-10_R8P !< Tollerance.

      ! evaluation of the intermediate states pressure and velocity
      a1  = sqrt(g1 * p1 / r1)                              ! left speed of sound
      a4  = sqrt(g4 * p4 / r4)                              ! right speed of sound
      ram = 0.5_R8P * (r1 + r4) * 0.5_R8P * (a1 + a4)       ! product of mean density for mean speed of sound
      u   = 0.5_R8P * (u1 + u4) - 0.5_R8P * (p4 - p1) / ram ! evaluation of the contact wave speed (velocity of intermediate states)
      p   = 0.5_R8P * (p1 + p4) - 0.5_R8P * (u4 - u1) * ram ! evaluation of the pressure of the intermediate states
      ! evaluation of the left wave speeds
      if (p<=p1*(1._R8P + toll)) then
        ! rarefaction
        S1 = u1 - a1
      else
        ! shock
        S1 = u1 - a1 * sqrt(1._R8P + (g1 + 1._R8P) / (2._R8P * g1) * (p / p1 - 1._R8P))
      endif
      ! evaluation of the right wave speeds
      if (p<=p4 * (1._R8P + toll)) then
        ! rarefaction
        S4 = u4 + a4
      else
        ! shock
        S4 = u4 + a4 * sqrt(1._R8P + (g4 + 1._R8P) / (2._R8P * g4) * ( p / p4 - 1._R8P))
      endif
      endsubroutine compute_inter_states
   endsubroutine solve_riemann_llf

   ! WENO FD
   subroutine weno_reconstruct_q(weno, gc, n, ns, np, cp0, cv0, q, qr)
   !< Reconstruct primitive variables in pseudo characteristic variables by WENO methods.
   type(weno_cpu_object), intent(in)  :: weno                     !< WENO solver.
   integer(I4P),          intent(in)  :: gc                       !< Number of ghost cells used.
   integer(I4P),          intent(in)  :: n                        !< Number of cells.
   integer(I4P),          intent(in)  :: ns                       !< Number of species.
   integer(I4P),          intent(in)  :: np                       !< Number of 1D primitive varibales.
   real(R8P),             intent(in)  :: cp0(1:ns), cv0(1:ns)     !< Specific heats of initial species.
   real(R8P),             intent(in)  :: q (1:np, 1-gc:n+gc)      !< 1D primitive variables.
   real(R8P),             intent(out) :: qr(1:np,1:2,0:n+1)       !< Reconstructed 1D primitive variables.
   real(R8P)                          :: qm(1:np,1:2)             !< Mean primitive variables.
   real(R8P)                          :: Lqm(1:ns+2,1:ns+2,1:2)   !< Left eigenvalues matrix of mean primitive variables.
   real(R8P)                          :: Rqm(1:ns+2,1:ns+2,1:2)   !< Right eigenvalues matrix of mean primitive variables.
   real(R8P)                          :: c(1:ns+2,1:2,1-gc:-1+gc) !< Pseudo characteristic variables.
   real(R8P)                          :: cr(1:ns+2,1:2)           !< Pseudo characteristic variables reconstructed.
   integer(I4P)                       :: i, j, f, v               !< Counter.

   select case(weno%weno_s)
   case(1_I4P)
      do i=0, n+1
         qr(:,1,i) = q(:,i)
         qr(:,2,i) = qr(:,1,i)
      enddo
   case(2_I4P,3_I4P,4_I4P)
   ! compute WENO reconstruction
      do i=0, n+1
         ! compute pseudo characteristic variables
         do f=1, 2
            if (i==0  .and.f==1) cycle
            if (i==n+1.and.f==2) cycle
            qm(:,f) = 0.5_R8P * (q(:,i+f-2) + q(:,i+f-1))
         enddo
         do f=1, 2
            if (i==0  .and.f==1) cycle
            if (i==n+1.and.f==2) cycle
            Lqm(:, :, f) =  left_eigenvectors(ns=ns, np=np, q=qm(:,f))
            Rqm(:, :, f) = right_eigenvectors(ns=ns, np=np, q=qm(:,f))
         enddo
         do j=i+1-weno%weno_s, i-1+weno%weno_s
            do f=1, 2
               if (i==0  .and.f==1) cycle
               if (i==n+1.and.f==2) cycle
               do v=1, ns+2
                  c(v,f,j-i) = dot_product(Lqm(v,1:ns+2,f), q(1:ns+2,j))
               enddo
            enddo
         enddo

         ! compute WENO reconstruction of pseudo charteristic variables
         do v=1, ns+2
            cr(v,:) = weno%reconstructed(s=weno%weno_s, v=c(v,:,1-weno%weno_s:-1+weno%weno_s))
         enddo

         ! trasform back reconstructed pseudo charteristic variables to primitive ones
         do f=1, 2
            if (i==0  .and.f==1) cycle
            if (i==n+1.and.f==2) cycle
            do v=1, ns+2
               qr(v,f,i) = dot_product(Rqm(v,1:ns+2,f), cr(1:ns+2,f))
            enddo
            qr(ns+3,f,i) = sum(qr(1:ns, f, i))
            qr(ns+4,f,i) = dot_product(qr(1:ns,f,i) / qr(ns+3,f,i), cp0) / &
                           dot_product(qr(1:ns,f,i) / qr(ns+3,f,i), cv0)
         enddo
      enddo
   endselect
   endsubroutine weno_reconstruct_q

   ! EOS
   function roe_averages(ns, np, ql, qr) result(qa)
   integer(I4P), intent(in) :: ns         !< Number of species.
   integer(I4P), intent(in) :: np         !< Number of 1D primitive varibales.
   real(R8P),    intent(in) :: ql(1:np)   !< Left primitive variables.
   real(R8P),    intent(in) :: qr(1:np)   !< Right primitive variables.
   real(R8P)                :: qa(1:np)   !< Roe averaged primitive variables.
   real(R8P)                :: hl, hr, ha !< Left, rigth and average states entalpy.
   real(R8P)                :: sigma      !< Roe's sigma factor.

   hl    = ic_total_entalpy(g=ql(Ns+4), density=ql(ns+3), pressure=ql(ns+2), velocity_sq_norm=ql(ns+1)*ql(ns+1))
   hr    = ic_total_entalpy(g=qr(Ns+4), density=qr(ns+3), pressure=qr(ns+2), velocity_sq_norm=qr(ns+1)*qr(ns+1))
   sigma = sqrt(ql(ns+3)) / (sqrt(ql(ns+3)) + sqrt(qr(ns+3)))
   ha    = sigma * hl + (1._R8P - sigma) * hr

   qa(1:ns) = sqrt(ql(1:ns)*qr(1:ns))
   qa(ns+3) = sqrt(ql(ns+3)*qr(ns+3))
   qa(ns+1) = sigma * ql(ns+1) + (1._R8P - sigma) * qr(ns+1)
   qa(ns+4) = sigma * ql(ns+4) + (1._R8P - sigma) * qr(ns+4)
   qa(ns+2) = (ha - 0.5_R8P*qa(ns+1)**2)*(qa(ns+4)-1._R8P)*qa(ns+3)/qa(ns+4)
   endfunction roe_averages

   pure function left_eigenvectors(ns, np, q) result(eig)
   !< Return the left eigenvectors matrix `L` as `dF/dP = A = R ^ L` `P`` being the primitive variables and `F` the fluxes.
   integer(I4P), intent(in) :: ns                  !< Number of species.
   integer(I4P), intent(in) :: np                  !< Number of 1D primitive varibales.
   real(R8P),    intent(in) :: q(1:np)             !< Primitive variables.
   real(R8P)                :: eig(1:ns+2, 1:ns+2) !< Eigenvectors.
   real(R8P)                :: gp                  !< `g*p`.
   real(R8P)                :: gp_a                !< `g*p/a`.
   integer(I4P)             :: s                   !< Counter.

   gp = q(ns+4) * q(ns+2)
   gp_a = gp / ic_speed_of_sound(g=q(ns+4), density=q(ns+3), pressure=q(ns+2))

   eig = 0._R8P

      eig(1,    ns+1) = -gp_a      ; eig(1,    ns+2) =  1._R8P
   do s=2, ns+1
      eig(s,    s-1 ) =  gp/q(s-1) ; eig(s,    ns+2) = -1._R8P
   enddo
      eig(ns+2, ns+1) =  gp_a      ; eig(ns+2, ns+2) =  1._R8P
   endfunction left_eigenvectors

   pure function right_eigenvectors(ns, np, q) result(eig)
   !< Return the right eigenvectors matrix `R` as `dF/dP = A = R ^ L` `P`` being the primitive variables and `F` the fluxes.
   integer(I4P), intent(in) :: ns                  !< Number of species.
   integer(I4P), intent(in) :: np                  !< Number of 1D primitive varibales.
   real(R8P),    intent(in) :: q(1:np)             !< Primitive variables.
   real(R8P)                :: eig(1:ns+2, 1:ns+2) !< Eigenvectors.
   real(R8P)                :: a                   !< Speed of sound.
   real(R8P)                :: gp                  !< `g*p`.
   real(R8P)                :: gp_inv              !< `1/(g*p)`.
   integer(I4P)             :: s                   !< Counter.

   a = ic_speed_of_sound(g=q(ns+4), density=q(ns+3), pressure=q(ns+2))
   gp = q(ns+4) * q(ns+2)
   gp_inv = 1._R8P / gp
   eig = 0._R8P

   do s=1, ns
     eig(s,   1) =  0.5_R8P * q(s) * gp_inv ; eig(s,s+1) = q(s) * gp_inv ; eig(s,   ns+2) =  eig(s,   1)
   enddo
     eig(ns+1,1) = -0.5_R8P * a    * gp_inv ;                              eig(ns+1,ns+2) = -eig(ns+1,1)
     eig(ns+2,1) =  0.5_R8P                 ;                              eig(ns+2,ns+2) =  eig(ns+2,1)
   endfunction right_eigenvectors
endmodule adam_s_solvers_cpu_object
