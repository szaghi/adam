!< ADAM, FLUME sink particles: the `[sinks]` input and the replicated sink list (issue #94, M12).
module adam_flume_sinks_object
!< ADAM, FLUME sink particles: the `[sinks]` input and the replicated sink list (issue #94, M12).
!<
!< A sink is a point mass with a position, a velocity and an acceleration. The list is replicated: every rank holds every
!< sink and updates it identically from MPI-reduced sums (D-M12-1), so no sink has an owner and none migrates. Sinks
!< act through self-gravity only (D-M12-2): their mass is CIC-deposited into the density the Poisson solver sees and
!< they feel the mesh gravity at their position (the deposit and the gather live in the common object, which owns the
!< grid). The keys, lists of one value per initial sink (their number is that of `mass`; none: no initial sink):
!<
!< | key | value |
!< |---|---|
!< | `x`, `y`, `z` | positions, inside the domain |
!< | `u`, `v`, `w` | velocities (optional, default 0) |
!< | `mass` | masses, > 0 |
!< | `creation` | optional, default `.false.`: sinks form where the gas meets the Federrath et al. (2010) checks |
!< | `accretion` | optional, default `.true.`: the sinks accrete the gas above the threshold within their radius |
!< | `merging` | optional, default `.true.`: two bound, converging sinks within an accretion radius merge |
!< | `accretion_radius_cells` | optional, default 2.5: the accretion radius in cells of the sink's level, in [1, ngc] |
!<
!< The threshold density is that of a Jeans length of two accretion radii, `rho_thr = pi c_s^2 / (4 G r_acc^2)`, with
!< the cell's own `c_s^2 = gamma p / rho` (the Jeans length of the M9 marker; issue #94, D-M12-3, D-M12-4).
!<
!< With the section, `[gravity] self_gravity` is required and the model must be hydrodynamic (fatal otherwise, issue
!< #94: MHD sinks are a follow-up). The `[reference]` layer converts the positions (length), the velocities and the
!< masses (`rho0 L0^3`). Restart: `<restart_basename>.sinks`, the list in 17 significant digits (read back, the same
!< doubles). History: `<output_basename>-sinks_history.dat`, one row per sink per save (`it time id m x y z u v w`).

! ADAM singleton objects
use :: adam_mpih_global,      only : mpih
! FLUME modules
use :: adam_flume_parameters, only : is_hydrodynamic
! third party modules
use :: finer,                 only : file_ini
use :: penf,                  only : I4P, R8P, str

implicit none
private
public :: flume_sinks_object

character(len=5), parameter :: INI_SECTION_NAME="sinks" !< INI section name.

type :: flume_sinks_object
   !< FLUME sink particles, code units.
   logical                :: is_active=.false.  !< True with a `[sinks]` section.
   integer(I4P)           :: n=0_I4P            !< Sinks number.
   real(R8P), allocatable :: m(:)               !< Masses [n].
   real(R8P), allocatable :: x(:,:)             !< Positions [3,n].
   real(R8P), allocatable :: v(:,:)             !< Velocities [3,n].
   real(R8P), allocatable :: a(:,:)             !< Accelerations from the mesh gravity [3,n].
   real(R8P)              :: dt_limit=huge(1._R8P) !< Time step limit of the current list (refreshed after each update).
   logical, allocatable   :: is_new(:)          !< True for a sink formed this step (its first accretion is unchecked).
   logical                :: creation=.false.   !< Sinks form from the gas.
   logical                :: accretion=.true.   !< Sinks accrete.
   logical                :: merging=.true.     !< Sinks merge.
   real(R8P)              :: r_acc_cells=2.5_R8P !< Accretion radius in cells.
   integer(I4P)           :: formed=0_I4P       !< Sinks formed since the start of the run.
   integer(I4P)           :: merged=0_I4P       !< Mergers since the start of the run.
   integer(I4P)           :: history_unit=0_I4P !< Sinks history unit (rank 0).
   contains
      ! public methods
      procedure, pass(self) :: append         !< Append a sink.
      procedure, pass(self) :: description    !< Return pretty-printed object description.
      procedure, pass(self) :: remove         !< Remove a sink.
      procedure, pass(self) :: load_from_file !< Load `[sinks]`.
      procedure, pass(self) :: load_restart   !< Load the list from `<restart_basename>.sinks`.
      procedure, pass(self) :: open_history   !< Open the sinks history (rank 0).
      procedure, pass(self) :: save_history   !< Write one history row per sink (rank 0).
      procedure, pass(self) :: save_restart   !< Save the list to `<restart_basename>.sinks` (rank 0).
      procedure, pass(self) :: totals         !< Return the total mass, momentum and kinetic energy.
endtype flume_sinks_object

contains
   ! public methods
   subroutine append(self, m, x, v)
   !< Append a sink (every rank, identically: the list is replicated); it is new until its first accretion. The arrays
   !< grow by explicit allocation (no reallocation on assignment: under nvfortran it intermittently failed to grow an
   !< array in the sink code, issue #94).
   class(flume_sinks_object), intent(inout) :: self !< Sinks.
   real(R8P),                 intent(in)    :: m    !< Mass.
   real(R8P),                 intent(in)    :: x(3) !< Position.
   real(R8P),                 intent(in)    :: v(3) !< Velocity.
   real(R8P), allocatable                   :: m_(:), x_(:,:), v_(:,:), a_(:,:) !< Grown arrays.
   logical, allocatable                     :: is_new_(:) !< Grown flags.

   allocate(m_(self%n + 1), x_(3,self%n + 1), v_(3,self%n + 1), a_(3,self%n + 1), is_new_(self%n + 1))
   m_(1:self%n) = self%m(1:self%n) ; x_(:,1:self%n) = self%x(:,1:self%n) ; v_(:,1:self%n) = self%v(:,1:self%n)
   a_(:,1:self%n) = self%a(:,1:self%n) ; is_new_(1:self%n) = self%is_new(1:self%n)
   m_(self%n + 1) = m ; x_(:,self%n + 1) = x ; v_(:,self%n + 1) = v ; a_(:,self%n + 1) = 0._R8P
   is_new_(self%n + 1) = .true.
   call move_alloc(m_, self%m) ; call move_alloc(x_, self%x) ; call move_alloc(v_, self%v)
   call move_alloc(a_, self%a) ; call move_alloc(is_new_, self%is_new)
   self%n = self%n + 1_I4P
   endsubroutine append

   subroutine remove(self, s)
   !< Remove sink `s` (every rank, identically), by explicit allocation (see `append`).
   class(flume_sinks_object), intent(inout) :: self !< Sinks.
   integer(I4P),              intent(in)    :: s    !< Sink to remove.
   real(R8P), allocatable                   :: m_(:), x_(:,:), v_(:,:), a_(:,:) !< Shrunk arrays.
   logical, allocatable                     :: is_new_(:) !< Shrunk flags.
   integer(I4P)                             :: k, n !< Counters.

   allocate(m_(self%n - 1), x_(3,self%n - 1), v_(3,self%n - 1), a_(3,self%n - 1), is_new_(self%n - 1))
   n = 0_I4P
   do k=1, self%n
      if (k == s) cycle
      n = n + 1_I4P
      m_(n) = self%m(k) ; x_(:,n) = self%x(:,k) ; v_(:,n) = self%v(:,k) ; a_(:,n) = self%a(:,k)
      is_new_(n) = self%is_new(k)
   enddo
   call move_alloc(m_, self%m) ; call move_alloc(x_, self%x) ; call move_alloc(v_, self%v)
   call move_alloc(a_, self%a) ; call move_alloc(is_new_, self%is_new)
   self%n = self%n - 1_I4P
   endsubroutine remove

   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_sinks_object), intent(in) :: self !< Sinks.
   character(len=:), allocatable         :: desc !< Description.
   integer(I4P)                          :: s    !< Counter.

   if (.not.self%is_active) then
      desc = mpih%myrankstr//'  sinks:           none'
      return
   endif
   desc = mpih%myrankstr//'  sinks:           '//trim(str(self%n, .true.))//' initial; creation '// &
          merge('T', 'F', self%creation)//', accretion '//merge('T', 'F', self%accretion)//', merging '// &
          merge('T', 'F', self%merging)// &
          ', accretion radius '//trim(str(self%r_acc_cells))//' cells'
   do s=1, self%n
      desc = desc//new_line('a')//mpih%myrankstr//'    m '//trim(str(self%m(s)))//' x '//trim(str(self%x(:,s)))// &
             ' v '//trim(str(self%v(:,s)))
   enddo
   endfunction description

   subroutine load_from_file(self, file_parameters, model, self_gravity, domain_emin, domain_emax, ngc)
   !< Load `[sinks]` and check it.
   class(flume_sinks_object), intent(inout) :: self            !< Sinks.
   type(file_ini),            intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P),              intent(in)    :: model           !< Physical model id.
   logical,                   intent(in)    :: self_gravity    !< `[gravity] self_gravity`.
   real(R8P),                 intent(in)    :: domain_emin(3)  !< Domain minimum abscissa.
   real(R8P),                 intent(in)    :: domain_emax(3)  !< Domain maximum abscissa.
   integer(I4P),              intent(in)    :: ngc             !< Ghost cells number (bounds the accretion radius).
   character(len=1), parameter              :: XYZ(3)=['x', 'y', 'z'] !< Position keys.
   character(len=1), parameter              :: UVW(3)=['u', 'v', 'w'] !< Velocity keys.
   integer(I4P)                             :: error           !< Error status.
   integer(I4P)                             :: d, s            !< Counters.

   self%is_active = file_parameters%has_section(section_name=INI_SECTION_NAME)
   if (.not.self%is_active) return
   if (.not.is_hydrodynamic(model)) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] needs a hydrodynamic physical_model (euler, euler-de): '// &
                               'MHD sinks are not implemented (issue #94)')
   if (.not.self_gravity) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] needs [gravity].(self_gravity) = .true.: sinks act through '// &
                               'the self-gravity solve (issue #94, D-M12-2)')
   ! FiNeR's has_option searches every section: an absent key of this section counts zero values
   self%n = file_parameters%count_values(section_name=INI_SECTION_NAME, option_name='mass')
   allocate(self%m(self%n), self%x(3,self%n), self%v(3,self%n), self%a(3,self%n))
   self%x = 0._R8P ; self%v = 0._R8P ; self%a = 0._R8P
   allocate(self%is_new(self%n))
   self%is_new = .false.
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='creation', val=self%creation, error=error)
   if (error /= 0) self%creation = .false.
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='accretion', val=self%accretion, error=error)
   if (error /= 0) self%accretion = .true.
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='merging', val=self%merging, error=error)
   if (error /= 0) self%merging = .true.
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='accretion_radius_cells', val=self%r_acc_cells, &
                            error=error)
   if (error /= 0) self%r_acc_cells = 2.5_R8P
   if (.not.(self%r_acc_cells >= 1._R8P .and. self%r_acc_cells <= real(ngc, R8P))) &
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(accretion_radius_cells) must lie in [1, ngc = '// &
                               trim(str(ngc, .true.))//'], got '//trim(str(self%r_acc_cells)))
   if (self%n == 0_I4P) return
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='mass', val=self%m, error=error)
   if (error /= 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(mass)')
   if (.not.all(self%m > 0._R8P)) call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(mass) must all be > 0')
   do d=1, 3
      call load_list(key=XYZ(d), required=.true., val=self%x(d,:))
      call load_list(key=UVW(d), required=.false., val=self%v(d,:))
   enddo
   do s=1, self%n
      if (any(self%x(:,s) < domain_emin) .or. any(self%x(:,s) > domain_emax)) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] sink '//trim(str(s, .true.))//' at '// &
                                  trim(str(self%x(:,s)))//' lies outside the domain')
   enddo
   contains
      subroutine load_list(key, required, val)
      !< Load one list of `n` values; an absent optional list is zero.
      character(*), intent(in)    :: key      !< Option name.
      logical,      intent(in)    :: required !< True if the option is required.
      real(R8P),    intent(inout) :: val(:)   !< Values.
      integer(I4P)                :: count_   !< Values number.

      count_ = file_parameters%count_values(section_name=INI_SECTION_NAME, option_name=key)
      if (count_ == 0_I4P) then
         if (required) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//key//')')
         val = 0._R8P
         return
      endif
      if (count_ /= self%n) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].('//key//') has '//trim(str(count_, .true.))// &
                                  ' values, (mass) '//trim(str(self%n, .true.)))
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name=key, val=val, error=error)
      if (error /= 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].('//key//')')
      endsubroutine load_list
   endsubroutine load_from_file

   subroutine load_restart(self, basename)
   !< Load the list from `<basename>.sinks` (every rank reads it: the list is replicated).
   class(flume_sinks_object), intent(inout) :: self       !< Sinks.
   character(*),              intent(in)    :: basename   !< Restart files basename.
   integer(I4P)                             :: unit, ios, s !< File unit, I/O status, counter.
   logical                                  :: is_present !< File presence.

   inquire(file=basename//'.sinks', exist=is_present)
   if (.not.is_present) call mpih%error_stop(msg=': the restart file '//basename//'.sinks is missing')
   open(newunit=unit, file=basename//'.sinks', action='read', iostat=ios)
   if (ios == 0) read(unit, *, iostat=ios) self%n
   if (ios /= 0) call mpih%error_stop(msg=': cannot read '//basename//'.sinks')
   if (allocated(self%m)) deallocate(self%m, self%x, self%v, self%a, self%is_new)
   allocate(self%m(self%n), self%x(3,self%n), self%v(3,self%n), self%a(3,self%n), self%is_new(self%n))
   self%a = 0._R8P
   self%is_new = .false.
   do s=1, self%n
      read(unit, *, iostat=ios) self%m(s), self%x(:,s), self%v(:,s)
      if (ios /= 0) call mpih%error_stop(msg=': cannot read sink '//trim(str(s, .true.))//' of '//basename//'.sinks')
   enddo
   close(unit)
   endsubroutine load_restart

   subroutine open_history(self, output_basename, is_restart)
   !< Open the sinks history (rank 0) and write its header (append, without header, on restart).
   class(flume_sinks_object), intent(inout) :: self            !< Sinks.
   character(*),              intent(in)    :: output_basename !< Output files basename.
   logical,                   intent(in)    :: is_restart      !< Restarted run flag.

   if (.not.self%is_active .or. mpih%myrank /= 0) return
   if (is_restart) then
      open(newunit=self%history_unit, file=trim(output_basename)//'-sinks_history.dat', status='unknown', &
           position='append', action='write')
   else
      open(newunit=self%history_unit, file=trim(output_basename)//'-sinks_history.dat', status='replace', &
           action='write')
      write(self%history_unit, '(A)') 'VARIABLES="it" "time" "id" "m" "x" "y" "z" "u" "v" "w"'
   endif
   endsubroutine open_history

   subroutine save_history(self, it, time, length, velocity, mass)
   !< Write one history row per sink (rank 0), the values times the output factors.
   class(flume_sinks_object), intent(in) :: self     !< Sinks.
   integer(I4P),              intent(in) :: it       !< Time iteration.
   real(R8P),                 intent(in) :: time     !< Time (output units).
   real(R8P),                 intent(in) :: length   !< Output factor of the lengths.
   real(R8P),                 intent(in) :: velocity !< Output factor of the velocities.
   real(R8P),                 intent(in) :: mass     !< Output factor of the masses.
   integer(I4P)                          :: s        !< Counter.

   if (mpih%myrank /= 0 .or. self%history_unit == 0_I4P) return
   do s=1, self%n
      write(self%history_unit, '(I10,1X,ES24.16E3,1X,I6,*(1X,ES24.16E3))') it, time, s, self%m(s) * mass, &
         self%x(:,s) * length, self%v(:,s) * velocity
   enddo
   flush(self%history_unit)
   endsubroutine save_history

   subroutine save_restart(self, basename)
   !< Save the list to `<basename>.sinks` (rank 0), 17 significant digits: read back, the same doubles.
   class(flume_sinks_object), intent(in) :: self     !< Sinks.
   character(*),              intent(in) :: basename !< Restart files basename.
   integer(I4P)                          :: unit, s  !< File unit, counter.

   if (mpih%myrank /= 0) return
   open(newunit=unit, file=basename//'.sinks', action='write', status='replace')
   write(unit, '(I0)') self%n
   do s=1, self%n
      write(unit, '(7(ES25.17E3,1X))') self%m(s), self%x(:,s), self%v(:,s)
   enddo
   close(unit)
   endsubroutine save_restart

   pure subroutine totals(self, mass, momentum, kinetic)
   !< Return the total mass, momentum and kinetic energy of the sinks.
   class(flume_sinks_object), intent(in)  :: self        !< Sinks.
   real(R8P),                 intent(out) :: mass        !< Total mass.
   real(R8P),                 intent(out) :: momentum(3) !< Total momentum.
   real(R8P),                 intent(out) :: kinetic     !< Total kinetic energy.
   integer(I4P)                           :: s           !< Counter.

   mass = 0._R8P ; momentum = 0._R8P ; kinetic = 0._R8P
   do s=1, self%n
      mass = mass + self%m(s)
      momentum = momentum + self%m(s) * self%v(:,s)
      kinetic = kinetic + 0.5_R8P * self%m(s) * sum(self%v(:,s)**2)
   enddo
   endsubroutine totals
endmodule adam_flume_sinks_object
