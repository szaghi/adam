!< ADAM, FLUME diagnostics class definition: conservation and div(B) histories.
module adam_flume_diagnostics_object
!< ADAM, FLUME diagnostics class definition: conservation and div(B) histories.
!<
!< Writes one row per save step to `<output_basename>-conservation_history.dat`: iteration, time and the volume
!< integral of every conservative variable. MHD runs also write, on the same cadence,
!< `<output_basename>-divb_history.dat`: iteration, time, `max|div B|`, `sum |div B| dV` and the maximum over the
!< seam-local cells (issue #41, D-10, section 9). With `divb_flux` (GLM and EGLM, issue #68 A2) four columns follow:
!< the D_flux norms (`max`, `sum dV`, seam-local `max`) and `max(h |D_flux|) / max|B|`, where D_flux is the divergence the
!< cleaning acts on, that of the interface normal field `B_n = F_psi / s` of the scheme's own psi face fluxes (`s =
!< c_h^2` GLM, `c_h` EGLM), the coarse faces of the 2:1 seams taking the fine fluxes (the reflux). Rows are written as
!< they come, nothing is accumulated in memory. The values are computed (and MPI-reduced) by the backends; only rank 0
!< writes. With self-gravity (issue #81) the conservation history gains two columns: the gravitational energy
!< `W = 1/2 sum rho Phi dV` (Phi solved from the committed state) and `E + W`, the energy the scheme conserves.
!< With radiative cooling (issue #93) two more follow: the energy radiated since the start of the run and `E` plus it.
!< With sinks (issue #94) five more: their total mass, momentum and kinetic energy (`W` then includes the sinks).
!<
!< With a subgrid-scale model (#119, T1) the turbulence history gains a last column, `eps_sgs`.
!< `turbulence_history_save = N` (issue #119, T0.1; optional, 0 by default) writes `<output_basename>-turbulence_history.dat`:
!< iteration, time and the volume integrals of the kinetic energy `K`, the enstrophy `Z = rho |omega|^2 / 2`, the solenoidal
!< and dilatational dissipations `eps_s = mu |omega|^2` and `eps_d = 4/3 mu theta^2`, the viscous work `eps = tau_ij
!< du_i/dx_j` and the pressure dilatation `pd = p theta`, from 4th-order central velocity gradients (code units, not divided
!< by the volume). In a periodic box `dK/dt = pd - eps` up to the scheme's numerical dissipation. Each row refreshes the
!< ghosts and the auxiliary variables of the committed state.

! ADAM singleton objects
use :: adam_mpih_global, only : mpih
! third party modules
use :: finer,            only : file_ini
use :: penf,             only : I4P, R8P, str
use :: stringifor,       only : string

implicit none
private
public :: flume_diagnostics_object

character(len=11), parameter :: INI_SECTION_NAME="diagnostics" !< INI (config) file section name.

type :: flume_diagnostics_object
   !< FLUME diagnostics class definition: conservation history.
   integer(I4P) :: conservation_history_save=0_I4P !< Conservation history save cadence (<= 0: disabled).
   logical      :: ghost_poison=.false.            !< Poison the ghosts before each field write (verification, #65).
   logical      :: divb_flux=.false.               !< Add the D_flux columns to the div(B) history (issue #68 A2).
   integer(I4P) :: conservation_unit=0_I4P         !< Conservation history file unit.
   integer(I4P) :: divb_unit=0_I4P                 !< div(B) history file unit (MHD only).
   integer(I4P) :: turbulence_history_save=0_I4P   !< Turbulence history save cadence (<= 0: disabled, issue #119).
   integer(I4P) :: turbulence_unit=0_I4P           !< Turbulence history file unit.
   integer(I4P) :: wall_history_save=0_I4P         !< Wall-model history save cadence (<= 0: disabled, issue #119, T3.7).
   integer(I4P) :: wall_history_bins=0_I4P         !< Wall-model history: bins along its axis.
   integer(I4P) :: wall_history_axis=1_I4P         !< Wall-model history: the axis of its bins (1 x, 2 y, 3 z).
   integer(I4P) :: wall_unit(6)=0_I4P              !< Wall-model history file unit per face.
   contains
      ! public methods
      procedure, pass(self) :: close_file            !< Close the history file.
      procedure, pass(self) :: description           !< Return pretty-printed object description.
      procedure, pass(self) :: initialize            !< Initialize diagnostics.
      procedure, pass(self) :: load_from_file        !< Load config from file.
      procedure, pass(self) :: open_file             !< Open the history file and write its header.
      procedure, pass(self) :: save_conservation_row !< Write one conservation history row.
      procedure, pass(self) :: save_divb_row         !< Write one div(B) history row.
      procedure, pass(self) :: save_turbulence_row   !< Write one turbulence history row.
endtype flume_diagnostics_object

contains
   ! public methods
   subroutine close_file(self)
   !< Close the history files.
   class(flume_diagnostics_object), intent(inout) :: self    !< Diagnostics.
   logical                                        :: is_open !< Unit status.
   integer(I4P)                                   :: f       !< Face.

   if (mpih%myrank /= 0) return
   do f=1, 6
      if (self%wall_unit(f) /= 0_I4P) then
         inquire(unit=self%wall_unit(f), opened=is_open)
         if (is_open) close(self%wall_unit(f))
         self%wall_unit(f) = 0_I4P
      endif
   enddo
   if (self%turbulence_unit /= 0_I4P) then
      inquire(unit=self%turbulence_unit, opened=is_open)
      if (is_open) close(self%turbulence_unit)
      self%turbulence_unit = 0_I4P
   endif
   if (self%conservation_unit /= 0_I4P) then
      inquire(unit=self%conservation_unit, opened=is_open)
      if (is_open) close(self%conservation_unit)
      self%conservation_unit = 0_I4P
   endif
   if (self%divb_unit /= 0_I4P) then
      inquire(unit=self%divb_unit, opened=is_open)
      if (is_open) close(self%divb_unit)
      self%divb_unit = 0_I4P
   endif
   endsubroutine close_file

   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_diagnostics_object), intent(in) :: self !< Diagnostics.
   character(len=:), allocatable               :: desc !< Description.

   desc = mpih%myrankstr//'Diagnostics main data'//new_line('a')// &
          mpih%myrankstr//'  conservation_history_save: '//trim(str(self%conservation_history_save))//new_line('a')// &
          mpih%myrankstr//'  divb_flux: '//trim(str(self%divb_flux))
   endfunction description

   subroutine initialize(self, file_parameters)
   !< Initialize diagnostics.
   class(flume_diagnostics_object), intent(inout) :: self            !< Diagnostics.
   type(file_ini),                  intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   print '(A)', mpih%myrankstr//'flume_diagnostics_object%initialize start'
   call self%load_from_file(file_parameters=file_parameters)
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'flume_diagnostics_object%initialize finish'
   endsubroutine initialize

   subroutine load_from_file(self, file_parameters)
   !< Load config from file; every key is required.
   class(flume_diagnostics_object), intent(inout) :: self            !< Diagnostics.
   type(file_ini),                  intent(in)    :: file_parameters !< Simulation parameters ini file handler.
   integer(I4P)                                   :: error           !< Error status.
   character(999)                                 :: buff            !< Buffer.

   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='conservation_history_save', &
                            val=self%conservation_history_save, error=error)
   if (error > 0) call mpih%error_stop(msg=': failed to load ['//INI_SECTION_NAME//'].(conservation_history_save)')
   ! optional (default .false.): a verification instrument (issue #65): before each field write every ghost cell is set
   ! to NaN, then refilled in the order of a stage (seams, exchange, boundary conditions), so a ghost that the fill
   ! does not reach, or reads before its donor is filled, is written as NaN instead of a stale value
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='ghost_poison', val=self%ghost_poison, &
                            error=error)
   if (error > 0) self%ghost_poison = .false.
   ! optional (default .false.): the D_flux columns of the div(B) history (issue #68 A2); each history row then costs one
   ! residual evaluation of the committed state (about one RK stage)
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='divb_flux', val=self%divb_flux, error=error)
   if (error > 0) self%divb_flux = .false.
   ! optional (default 0, disabled): the turbulence history (issue #119, T0.1)
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='turbulence_history_save', &
                            val=self%turbulence_history_save, error=error)
   if (error > 0) self%turbulence_history_save = 0_I4P
   ! optional (default 0, disabled): the wall-model history (issue #119, T3.7): wall_history_bins (required with it) bins along wall_history_axis
   call file_parameters%get(section_name=INI_SECTION_NAME, option_name='wall_history_save', val=self%wall_history_save, error=error)
   if (error > 0) self%wall_history_save = 0_I4P
   if (self%wall_history_save > 0_I4P) then
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='wall_history_bins', val=self%wall_history_bins, error=error)
      if (error > 0 .or. self%wall_history_bins < 1_I4P) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(wall_history_bins) must be given and >= 1 with wall_history_save')
      call file_parameters%get(section_name=INI_SECTION_NAME, option_name='wall_history_axis', val=buff, error=error)
      if (error > 0) buff = 'x'
      select case(trim(adjustl(buff)))
      case('x')
         self%wall_history_axis = 1_I4P
      case('y')
         self%wall_history_axis = 2_I4P
      case('z')
         self%wall_history_axis = 3_I4P
      case default
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(wall_history_axis) must be x, y or z')
      endselect
   endif
   endsubroutine load_from_file

   subroutine open_file(self, output_basename, q_name, is_restart, with_divb, with_self_gravity, with_cooling, with_sinks, &
                        with_sgs)
   !< Open the history files and write their headers (append, without header, on restart).
   class(flume_diagnostics_object), intent(inout)        :: self              !< Diagnostics.
   character(*),                    intent(in)           :: output_basename   !< Output files base name.
   type(string),                    intent(in)           :: q_name(1:)        !< Conservative variables names.
   logical,                         intent(in)           :: is_restart        !< Restarted run flag.
   logical,                         intent(in)           :: with_divb         !< Open the div(B) history too (MHD).
   logical,                         intent(in), optional :: with_self_gravity !< Add the W and E + W columns.
   logical,                         intent(in), optional :: with_cooling      !< Add the radiated and E + radiated columns.
   logical,                         intent(in), optional :: with_sinks        !< Add the sinks' mass, momentum, kinetic energy.
   logical,                         intent(in), optional :: with_sgs          !< Add eps_sgs to the turbulence history (#119).
   character(:), allocatable                      :: header          !< File header.
   integer(I4P)                                   :: v               !< Counter.
   logical                                        :: with_sgs_       !< Add eps_sgs, local variable.

   if (mpih%myrank /= 0) return
   with_sgs_ = .false. ; if (present(with_sgs)) with_sgs_ = with_sgs
   if (self%turbulence_history_save > 0_I4P) then
      if (is_restart) then
         open(newunit=self%turbulence_unit, file=trim(output_basename)//'-turbulence_history.dat', &
              status='unknown', position='append', action='write')
      else
         open(newunit=self%turbulence_unit, file=trim(output_basename)//'-turbulence_history.dat', status='replace', &
              action='write')
         write(self%turbulence_unit, '(A)') 'VARIABLES="it" "time" "K" "Z" "eps_s" "eps_d" "eps" "pd"'//   &
                                            trim(merge(' "eps_sgs"', '          ', with_sgs_))
      endif
   endif
   if (self%conservation_history_save <= 0_I4P) return
   if (with_divb) then
      if (is_restart) then
         open(newunit=self%divb_unit, file=trim(output_basename)//'-divb_history.dat', &
              status='unknown', position='append', action='write')
      else
         open(newunit=self%divb_unit, file=trim(output_basename)//'-divb_history.dat', status='replace', action='write')
         header = 'VARIABLES="it" "time" "max_divb" "l1_divb" "seam_max_divb"'
         if (self%divb_flux) header = header//' "max_dflux" "l1_dflux" "seam_max_dflux" "max_h_dflux_over_bmax"'
         write(self%divb_unit, '(A)') header
      endif
   endif
   if (is_restart) then
      open(newunit=self%conservation_unit, file=trim(output_basename)//'-conservation_history.dat', &
           status='unknown', position='append', action='write')
   else
      open(newunit=self%conservation_unit, file=trim(output_basename)//'-conservation_history.dat', &
           status='replace', action='write')
      header = 'VARIABLES="it" "time"'
      do v=1, size(q_name, dim=1)
         header = header//' "int_'//q_name(v)%chars()//'"'
      enddo
      if (present(with_self_gravity)) then
         if (with_self_gravity) header = header//' "W" "E_plus_W"'
      endif
      if (present(with_cooling)) then
         if (with_cooling) header = header//' "radiated" "E_plus_radiated"' ! issue #93
      endif
      if (present(with_sinks)) then
         if (with_sinks) header = header//' "M_sinks" "Px_sinks" "Py_sinks" "Pz_sinks" "K_sinks"' ! issue #94
      endif
      write(self%conservation_unit, '(A)') header
   endif
   endsubroutine open_file

   subroutine save_conservation_row(self, it, time, integrals)
   !< Write one conservation history row (rank 0 only; integrals already MPI-reduced).
   class(flume_diagnostics_object), intent(in) :: self          !< Diagnostics.
   integer(I4P),                    intent(in) :: it            !< Time iteration.
   real(R8P),                       intent(in) :: time          !< Time.
   real(R8P),                       intent(in) :: integrals(1:) !< Volume integrals of the conservative variables.

   if (mpih%myrank /= 0 .or. self%conservation_unit == 0_I4P) return
   write(self%conservation_unit, '(I10,*(1X,ES24.16E3))') it, time, integrals
   flush(self%conservation_unit)
   endsubroutine save_conservation_row

   subroutine save_turbulence_row(self, it, time, integrals)
   !< Write one turbulence history row (rank 0 only; integrals already MPI-reduced).
   class(flume_diagnostics_object), intent(in) :: self         !< Diagnostics.
   integer(I4P),                    intent(in) :: it           !< Time iteration.
   real(R8P),                       intent(in) :: time         !< Time.
   real(R8P),                       intent(in) :: integrals(:) !< K, Z, eps_s, eps_d, eps, pd (and eps_sgs).

   if (mpih%myrank /= 0 .or. self%turbulence_unit == 0_I4P) return
   write(self%turbulence_unit, '(I10,*(1X,ES24.16E3))') it, time, integrals
   flush(self%turbulence_unit)
   endsubroutine save_turbulence_row

   subroutine save_divb_row(self, it, time, norms)
   !< Write one div(B) history row (rank 0 only; norms already MPI-reduced).
   class(flume_diagnostics_object), intent(in) :: self     !< Diagnostics.
   integer(I4P),                    intent(in) :: it       !< Time iteration.
   real(R8P),                       intent(in) :: time     !< Time.
   real(R8P),                       intent(in) :: norms(:) !< max|div B|, sum |div B| dV, seam-local max|div B|, and
                                                            !< with divb_flux the four D_flux columns.

   if (mpih%myrank /= 0 .or. self%divb_unit == 0_I4P) return
   write(self%divb_unit, '(I10,*(1X,ES24.16E3))') it, time, norms
   flush(self%divb_unit)
   endsubroutine save_divb_row
endmodule adam_flume_diagnostics_object
