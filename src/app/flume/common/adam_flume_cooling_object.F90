!< ADAM, FLUME optically thin radiative cooling: the `[cooling]` input (issue #93, M11).
module adam_flume_cooling_object
!< ADAM, FLUME optically thin radiative cooling: the `[cooling]` input (issue #93, M11).
!<
!< `d(rho e)/dt = -n^2 Lambda(T)`, `n = rho / particle_mass`, `T = p / (rho R)`, applied once per step on the committed
!< state by Townsend's exact integration (`adam_flume_cooling_library`; D-M11-1: operator split, first order in the
!< coupling with the dynamics, exact and stable for any `dt / t_cool`). The keys:
!<
!< | key | value |
!< |---|---|
!< | `temperatures` | the node temperatures `T_1 < ... < T_N`, `2 <= N <= COOLING_NODES_MAX` (64), all > 0 (a list) |
!< | `lambdas` | `Lambda_k > 0` at the nodes (a list, as many as the temperatures) |
!< | `particle_mass` | `m > 0`, the mass per particle of `n` (for `n_H`: `m_H / X`) |
!< | `dt_factor` | optional, default 0 (off): also `dt <= dt_factor * min t_cool`, `t_cool = rho e / (n^2 Lambda)` |
!<
!< `Lambda` is the piecewise power law through the nodes; above `T_N` its last segment continues, below `T_1` there is
!< no cooling (D-M11-2: cooling only, `T_1` the floor, reached exactly). The hydrodynamic models only (D-M11-3): with
!< `mhd-ideal` the section is fatal. Every key, the gas constant and the temperatures are in the units of the input; the
!< `[reference]` layer converts the temperatures, `Lambda` (`L0^5 u0^3 rho0`) and `particle_mass` (`rho0 L0^3`).
!< Without the section `is_active` is false and nothing changes.
!<
!< The object precomputes the segment exponents `alpha_k` and the evolution function at the nodes, `Y_k`, so a cell
!< step costs two segment searches, two powers and one logarithm. One step raises `Y` by `dt * dy_rate * rho`,
!< `dy_rate = (gamma - 1) Lambda_N / (m^2 R T_N)` (`n^2 / rho = rho / m^2`).

! ADAM singleton objects
use :: adam_mpih_global,           only : mpih
! FLUME modules
use :: adam_flume_cooling_library, only : cooling_segment_y
use :: adam_flume_parameters,      only : COOLING_NODES_MAX, is_hydrodynamic
! third party modules
use :: finer,                      only : file_ini
use :: penf,                       only : I4P, R8P, str

implicit none
private
public :: flume_cooling_object

character(len=7), parameter :: INI_SECTION_NAME="cooling" !< INI section name.

type :: flume_cooling_object
   !< FLUME radiative cooling, code units.
   logical                :: is_active=.false.     !< True with a `[cooling]` section.
   integer(I4P)           :: n_nodes=0_I4P         !< Nodes number.
   real(R8P), allocatable :: t_nodes(:)            !< Node temperatures [n_nodes].
   real(R8P), allocatable :: l_nodes(:)            !< Node values of Lambda [n_nodes].
   real(R8P), allocatable :: alpha(:)              !< Segment exponents [n_nodes], alpha(N) = alpha(N-1).
   real(R8P), allocatable :: y_nodes(:)            !< Evolution function at the nodes [n_nodes], Y_N = 0.
   real(R8P)              :: particle_mass=0._R8P  !< Mass per particle.
   real(R8P)              :: dt_factor=0._R8P      !< Time step limiter factor, 0 = off.
   real(R8P)              :: dy_rate=0._R8P        !< (gamma - 1) Lambda_N / (m^2 R T_N).
   real(R8P)              :: radiated=0._R8P       !< Energy radiated since the start of the run (all ranks).
   contains
      ! public methods
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: load_from_file !< Load `[cooling]` and precompute the tables.
endtype flume_cooling_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_cooling_object), intent(in) :: self !< Cooling.
   character(len=:), allocatable           :: desc !< Description.
   integer(I4P)                            :: k    !< Counter.

   if (.not.self%is_active) then
      desc = mpih%myrankstr//'  cooling:         none'
      return
   endif
   desc = mpih%myrankstr//'  cooling:         '//trim(str(self%n_nodes, .true.))//' nodes, particle mass '// &
          trim(str(self%particle_mass))//', dt_factor '//trim(str(self%dt_factor))
   do k=1, self%n_nodes
      desc = desc//new_line('a')//mpih%myrankstr//'    T '//trim(str(self%t_nodes(k)))//' Lambda '// &
             trim(str(self%l_nodes(k)))//' Y '//trim(str(self%y_nodes(k)))
   enddo
   endfunction description

   subroutine load_from_file(self, file_parameters, model, gamma, R)
   !< Load `[cooling]`, check it and precompute the exponents and the evolution function at the nodes.
   class(flume_cooling_object), intent(inout) :: self            !< Cooling.
   type(file_ini),              intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P),                intent(in)    :: model           !< Physical model id.
   real(R8P),                   intent(in)    :: gamma           !< Specific heats ratio.
   real(R8P),                   intent(in)    :: R               !< Gas constant.
   integer(I4P)                               :: n_lambdas       !< Number of Lambda values.
   integer(I4P)                               :: error           !< Error status.
   integer(I4P)                               :: k               !< Counter.

   self%is_active = file_parameters%has_section(section_name=INI_SECTION_NAME)
   if (.not.self%is_active) return
   if (.not.is_hydrodynamic(model)) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] needs a hydrodynamic physical_model (euler, euler-de): '// &
                               'MHD cooling is not implemented (issue #93, D-M11-3)')
   self%n_nodes = file_parameters%count_values(section_name=INI_SECTION_NAME, option_name='temperatures')
   n_lambdas = file_parameters%count_values(section_name=INI_SECTION_NAME, option_name='lambdas')
   if (self%n_nodes < 2_I4P) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(temperatures) needs at least 2 nodes, got '// &
                               trim(str(self%n_nodes, .true.)))
   if (self%n_nodes > COOLING_NODES_MAX) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(temperatures) has '//trim(str(self%n_nodes, .true.))// &
                               ' nodes, at most '//trim(str(COOLING_NODES_MAX, .true.))//' are accepted')
   if (n_lambdas /= self%n_nodes) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(lambdas) has '//trim(str(n_lambdas, .true.))// &
                               ' values, (temperatures) '//trim(str(self%n_nodes, .true.)))
   allocate(self%t_nodes(self%n_nodes), self%l_nodes(self%n_nodes), self%alpha(self%n_nodes), &
            self%y_nodes(self%n_nodes))
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='temperatures', val=self%t_nodes, error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(temperatures)')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='lambdas', val=self%l_nodes, error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(lambdas)')
   if (.not.all(self%t_nodes > 0._R8P)) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(temperatures) must all be > 0')
   if (.not.all(self%t_nodes(2:) > self%t_nodes(:self%n_nodes-1))) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(temperatures) must be strictly increasing')
   if (.not.all(self%l_nodes > 0._R8P)) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(lambdas) must all be > 0')
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='particle_mass', val=self%particle_mass, &
                            error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(particle_mass)')
   if (.not.(self%particle_mass > 0._R8P)) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(particle_mass) must be > 0, got '// &
                               trim(str(self%particle_mass)))
   self%dt_factor = 0._R8P
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='dt_factor', val=self%dt_factor, error=error)
   if (error /= 0) self%dt_factor = 0._R8P
   if (self%dt_factor < 0._R8P) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(dt_factor) must be >= 0 (0: off), got '// &
                               trim(str(self%dt_factor)))
   ! the tables: the exponents, then Y from the top node down, Y_N = 0 and Y_k = Y_k+1 - (Y(T_k+1) - Y_k on segment k)
   do k=1, self%n_nodes - 1_I4P
      self%alpha(k) = log(self%l_nodes(k+1) / self%l_nodes(k)) / log(self%t_nodes(k+1) / self%t_nodes(k))
   enddo
   self%alpha(self%n_nodes) = self%alpha(self%n_nodes - 1_I4P)
   self%y_nodes = 0._R8P
   do k=self%n_nodes - 1_I4P, 1_I4P, -1_I4P
      self%y_nodes(k) = self%y_nodes(k+1) - cooling_segment_y(n_nodes=self%n_nodes, t_nodes=self%t_nodes,          &
                                                              l_nodes=self%l_nodes, alpha=self%alpha, k=k,         &
                                                              temperature=self%t_nodes(k+1))
   enddo
   self%dy_rate = (gamma - 1._R8P) * self%l_nodes(self%n_nodes) / &
                  (self%particle_mass**2 * R * self%t_nodes(self%n_nodes))
   self%radiated = 0._R8P
   endsubroutine load_from_file
endmodule adam_flume_cooling_object
