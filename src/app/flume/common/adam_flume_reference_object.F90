!< ADAM, FLUME reference layer (`[reference]`): dimensional input converted to code units on the host (issue #49, N2).
module adam_flume_reference_object
!< ADAM, FLUME reference layer (`[reference]`): dimensional input converted to code units on the host (issue #49, N2).
!<
!< **Why the conversion is all there is.** Ideal Euler and MHD as FLUME writes them (`B = B_SI/sqrt(mu0)`) carry no
!< dimensionless number: with three references, density `rho0`, length `L0` and velocity `u0`, and the derived time
!< `L0/u0`, pressure `rho0 u0^2` and field `u0 sqrt(rho0)`, every term of the equations picks up the same factor. A
!< dimensional input therefore runs unchanged once each value is divided by the reference of its dimension, and the
!< kernels never see the references (issue #49, ND-1, ND-2).
!<
!< **Where.** Once, right after the input file is loaded: every dimensional option of the loaded INI object is
!< replaced by its value in code units (17 significant digits, which read back to the same double), before any
!< parser reads it. FLUME's objects and the library's (grid, AMR markers, solids, slices) then parse code units
!< without knowing about the layer. Without a `[reference]` section nothing is touched (ND-6: byte-identical runs).
!<
!< **Completeness.** With the layer active, every section and every option of the file must be classified here: an
!< unknown section or option is fatal, so no dimensional value can pass unconverted (risk R-5 of #49). Options whose
!< dimension depends on the context are resolved from it: the AMR gradient tolerance takes the dimension of the marked
!< variable per length, the linear-wave amplitude that of the eigenvector it multiplies.
!<
!< **Gas.** `cp`, `cv` are replaced by `gamma = cp / cv` (the same double the physics object computes), so that
!< `R* = 1` and the temperature unit is `u0^2 / R`; with `[physics] gamma` there is nothing to convert.
!<
!< **Restart.** Restart files stay in code units; `<restart_basename>.reference` records the references they were
!< written with (1 without the layer), and a restart under different references is refused (NV-8). Restart files
!< without the record (written before issue #49 N2b) are code units: they restart only with references 1.
!<
!< **Output.** `[reference] output_units = code` (default) writes code units; `dimensional` multiplies, at write time
!< only, the fields, the grid, the time, the slices and the histories by the reference of their dimension (ND-7:
!< normalisation and output units are separate settings). The temperature is written as `p / (rho R)` with the gas
!< constant of the input (`cp - cv` when the layer replaced them by `gamma`). `<output_basename>.units` records the
!< output units, the references and the factor of every written variable.

! ADAM singleton objects
use :: adam_mpih_global,      only : mpih
! FLUME modules
use :: adam_flume_parameters, only : strip_control
! third party modules
use :: finer,                 only : file_ini
use :: penf,                  only : I4P, R8P, str
use :: stringifor,            only : string

implicit none
private
public :: flume_reference_object

character(len=9), parameter :: INI_SECTION_NAME="reference" !< INI (config) file section name of the references.

! dimensions as the exponents of (length, velocity, density), the density one doubled (the field scales as sqrt(rho))
integer(I4P), parameter :: DIM_NONE(3)=[0_I4P, 0_I4P, 0_I4P]               !< Dimensionless.
integer(I4P), parameter :: DIM_LENGTH(3)=[1_I4P, 0_I4P, 0_I4P]             !< Length.
integer(I4P), parameter :: DIM_INV_LENGTH(3)=[-1_I4P, 0_I4P, 0_I4P]        !< Inverse length.
integer(I4P), parameter :: DIM_TIME(3)=[1_I4P, -1_I4P, 0_I4P]              !< Time.
integer(I4P), parameter :: DIM_VELOCITY(3)=[0_I4P, 1_I4P, 0_I4P]           !< Velocity.
integer(I4P), parameter :: DIM_VELOCITY2(3)=[0_I4P, 2_I4P, 0_I4P]          !< Velocity squared (specific enthalpy).
integer(I4P), parameter :: DIM_DENSITY(3)=[0_I4P, 0_I4P, 2_I4P]            !< Density.
integer(I4P), parameter :: DIM_MOMENTUM(3)=[0_I4P, 1_I4P, 2_I4P]           !< Momentum density.
integer(I4P), parameter :: DIM_PRESSURE(3)=[0_I4P, 2_I4P, 2_I4P]           !< Pressure, energy density.
integer(I4P), parameter :: DIM_FIELD(3)=[0_I4P, 1_I4P, 1_I4P]              !< Magnetic field, B_SI/sqrt(mu0).
integer(I4P), parameter :: DIM_PSI_GLM(3)=[0_I4P, 2_I4P, 1_I4P]            !< GLM psi: field times velocity.
integer(I4P), parameter :: DIM_FIELD_PER_LENGTH(3)=[-1_I4P, 1_I4P, 1_I4P]  !< div(B).
integer(I4P), parameter :: DIM_VISCOSITY(3)=[1_I4P, 1_I4P, 2_I4P]         !< Dynamic viscosity, rho0 u0 L0.
integer(I4P), parameter :: DIM_DIFFUSIVITY(3)=[1_I4P, 1_I4P, 0_I4P]       !< Diffusivity, u0 L0 (resistivity).
! option kinds
integer(I4P), parameter :: KIND_NONE=0_I4P    !< Not converted.
integer(I4P), parameter :: KIND_SCALED=1_I4P  !< Divided by the reference of its dimension.
integer(I4P), parameter :: KIND_GAS=2_I4P     !< `cp`, `cv`: replaced by `gamma`, `R* = 1`.
integer(I4P), parameter :: KIND_UNKNOWN=3_I4P !< Not classified: fatal.
integer(I4P), parameter :: KIND_TEMPERATURE=4_I4P  !< A temperature: divided by `u0^2 / R` (issue #65).
integer(I4P), parameter :: KIND_CONDUCTIVITY=5_I4P !< A thermal conductivity: divided by `rho0 u0 L0 R` (issue #65).

type :: flume_reference_object
   !< FLUME reference layer class definition.
   logical                   :: is_active=.false.    !< True if the input has a `[reference]` section.
   real(R8P)                 :: density=1._R8P       !< Reference density, rho0.
   real(R8P)                 :: length=1._R8P        !< Reference length, L0.
   real(R8P)                 :: velocity=1._R8P      !< Reference velocity, u0.
   character(:), allocatable :: velocity_preset      !< `value`, `acoustic` or `alfvenic`.
   real(R8P)                 :: ref_pressure=0._R8P  !< Reference pressure of the `acoustic` preset.
   real(R8P)                 :: ref_field=0._R8P     !< Reference field of the `alfvenic` preset.
   logical                   :: dimensional_output=.false. !< `output_units = dimensional`.
   real(R8P)                 :: gas_constant=1._R8P  !< Gas constant of the input, `cp - cv` (1 with `gamma`).
   integer(I4P)              :: converted=0_I4P      !< Options converted.
   ! context of the context-dependent options, read raw from the input
   character(:), allocatable :: physical_model       !< `[physics] physical_model`.
   character(:), allocatable :: divergence_control   !< `[mhd] divergence_control` (`none` if absent).
   character(:), allocatable :: ic_type              !< `[initial_conditions] type`.
   character(:), allocatable :: ic_wave              !< `[initial_conditions] wave` (empty if absent).
   contains
      ! public methods
      procedure, pass(self) :: description !< Return pretty-printed object description.
      procedure, pass(self) :: check_restart !< Refuse restart files written under other references.
      procedure, pass(self) :: initialize    !< Initialize the layer and convert the input.
      procedure, pass(self) :: length_output !< Return the output factor of the lengths.
      procedure, pass(self) :: save_restart  !< Record the references of the restart files.
      procedure, pass(self) :: save_units    !< Record the output units and the factor of every written variable.
      procedure, pass(self) :: scale         !< Return the reference of a dimension.
      procedure, pass(self) :: time_output   !< Return the output factor of the time.
      procedure, pass(self) :: variable_output !< Return the output factor of a written variable.
      ! private methods
      procedure, pass(self), private :: classify         !< Return the kind and dimension of an option.
      procedure, pass(self), private :: classify_marker  !< Return the dimension of an AMR gradient tolerance.
      procedure, pass(self), private :: convert_gas      !< Replace `cp`, `cv` by `gamma` (`R* = 1`).
      procedure, pass(self), private :: convert_options  !< Convert every dimensional option of the input.
      procedure, pass(self), private :: load_context     !< Load the context of the context-dependent options.
      procedure, pass(self), private :: load_from_file   !< Load `[reference]`.
endtype flume_reference_object

contains
   ! public methods
   function description(self) result(desc)
   !< Return a pretty-formatted object description.
   class(flume_reference_object), intent(in) :: self             !< Reference layer.
   character(len=:), allocatable             :: desc             !< Description.
   character(len=1), parameter               :: NL=new_line('a') !< New line character.

   desc =       mpih%myrankstr//'Reference main data'//NL
   if (.not.self%is_active) then
      desc = desc//mpih%myrankstr//'  no [reference]: the input is in code units'
      return
   endif
   desc = desc//mpih%myrankstr//'  density, length:   '//trim(str(self%density))//', '//trim(str(self%length))//NL
   desc = desc//mpih%myrankstr//'  velocity:          '//trim(str(self%velocity))//' ('//self%velocity_preset//')'//NL
   desc = desc//mpih%myrankstr//'  time, pressure:    '//trim(str(self%scale(DIM_TIME)))//', '// &
                                                       trim(str(self%scale(DIM_PRESSURE)))//NL
   desc = desc//mpih%myrankstr//'  field:             '//trim(str(self%scale(DIM_FIELD)))//NL
   if (self%dimensional_output) then
      desc = desc//mpih%myrankstr//'  output units:      dimensional'//NL
   else
      desc = desc//mpih%myrankstr//'  output units:      code'//NL
   endif
   desc = desc//mpih%myrankstr//'  options converted: '//trim(str(self%converted, .true.))
   endfunction description

   subroutine check_restart(self, basename)
   !< Refuse restart files written under other references: the record must equal the current references exactly; without
   !< a record (restart files written before issue #49 N2b, in code units) the current references must all be 1.
   class(flume_reference_object), intent(in) :: self      !< Reference layer.
   character(*),                  intent(in) :: basename  !< Restart files basename.
   character(len=*), parameter               :: ORDER='density, length, velocity' !< Record order.
   real(R8P)                                 :: saved(3)  !< References of the restart files.
   real(R8P)                                 :: current(3) !< Current references.
   character(999)                            :: key       !< Record key.
   integer(I4P)                              :: unit, ios, i !< File unit, I/O status, counter.
   logical                                   :: is_present !< Record presence.

   current = [self%density, self%length, self%velocity]
   inquire(file=basename//'.reference', exist=is_present)
   if (.not.is_present) then
      if (any(current /= 1._R8P)) &
         call mpih%error_stop(msg=': the restart files '//basename//' carry no references (code units, written before '// &
                                  'issue #49 N2b) and the run has ['//INI_SECTION_NAME//'] references not 1: refused')
      return
   endif
   open(newunit=unit, file=basename//'.reference', action='read', iostat=ios)
   do i=1, 3
      if (ios == 0) read(unit, *, iostat=ios) key, saved(i)
   enddo
   close(unit)
   if (ios /= 0) call mpih%error_stop(msg=': cannot read the references of the restart files '//basename//'.reference')
   if (any(saved /= current)) &
      call mpih%error_stop(msg=': the restart files '//basename//' were written with the references ('//ORDER//') '// &
                               trim(format_value(saved(1)))//', '//trim(format_value(saved(2)))//', '//          &
                               trim(format_value(saved(3)))//', the run has '//trim(format_value(current(1)))//   &
                               ', '//trim(format_value(current(2)))//', '//trim(format_value(current(3)))//       &
                               ': refused (restart files are in code units of their own references)')
   endsubroutine check_restart

   subroutine initialize(self, file_parameters)
   !< Initialize the layer: without `[reference]` do nothing; with it, load the references and convert every
   !< dimensional option of the loaded input to code units.
   class(flume_reference_object), intent(inout) :: self            !< Reference layer.
   type(file_ini),                intent(inout) :: file_parameters !< Simulation parameters ini file handler.

   print '(A)', mpih%myrankstr//'flume_reference_object%initialize start'
   self%is_active = file_parameters%has_section(section_name=INI_SECTION_NAME)
   if (self%is_active) then
      call self%load_context(file_parameters=file_parameters)
      call self%load_from_file(file_parameters=file_parameters)
      call self%convert_options(file_parameters=file_parameters)
   endif
   print '(A)', self%description()
   print '(A)', mpih%myrankstr//'flume_reference_object%initialize finish'
   endsubroutine initialize

   pure function length_output(self) result(s)
   !< Return the output factor of the lengths: `L0` with dimensional output units, 1 otherwise.
   class(flume_reference_object), intent(in) :: self !< Reference layer.
   real(R8P)                                 :: s    !< Output factor.

   s = 1._R8P
   if (self%dimensional_output) s = self%scale(DIM_LENGTH)
   endfunction length_output

   subroutine save_restart(self, basename)
   !< Record the references of the restart files (rank 0), 17 significant digits: read back, the same doubles.
   class(flume_reference_object), intent(in) :: self     !< Reference layer.
   character(*),                  intent(in) :: basename !< Restart files basename.
   integer(I4P)                              :: unit     !< File unit.

   if (mpih%myrank /= 0_I4P) return
   open(newunit=unit, file=basename//'.reference', action='write', status='replace')
   write(unit, '(A)') 'density '//trim(format_value(self%density))
   write(unit, '(A)') 'length '//trim(format_value(self%length))
   write(unit, '(A)') 'velocity '//trim(format_value(self%velocity))
   close(unit)
   endsubroutine save_restart

   subroutine save_units(self, basename, names)
   !< Record (rank 0) the output units in `<basename>.units`: the references, the derived ones, the gas constant and the
   !< factor by which every written variable is multiplied; the history columns follow from them (listed as comments).
   class(flume_reference_object), intent(in) :: self     !< Reference layer.
   character(*),                  intent(in) :: basename !< Output files basename.
   type(string),                  intent(in) :: names(:) !< Names of the written variables.
   integer(I4P)                              :: unit, v  !< File unit, counter.

   if (mpih%myrank /= 0_I4P) return
   open(newunit=unit, file=basename//'.units', action='write', status='replace')
   write(unit, '(A)') '# FLUME output units (issue #49, N2c): written value = code value * factor'
   if (self%dimensional_output) then
      write(unit, '(A)') 'output_units dimensional'
   else
      write(unit, '(A)') 'output_units code'
   endif
   write(unit, '(A)') 'density '//trim(format_value(self%density))
   write(unit, '(A)') 'length '//trim(format_value(self%length))
   write(unit, '(A)') 'velocity '//trim(format_value(self%velocity))
   write(unit, '(A)') 'time '//trim(format_value(self%scale(DIM_TIME)))
   write(unit, '(A)') 'pressure '//trim(format_value(self%scale(DIM_PRESSURE)))
   write(unit, '(A)') 'field '//trim(format_value(self%scale(DIM_FIELD)))
   write(unit, '(A)') 'gas_constant '//trim(format_value(self%gas_constant))
   write(unit, '(A)') 'factor_length '//trim(format_value(self%length_output()))
   write(unit, '(A)') 'factor_time '//trim(format_value(self%time_output()))
   do v=1, size(names, dim=1)
      write(unit, '(A)') 'factor '//names(v)%chars()//' '//trim(format_value(self%variable_output(names(v)%chars())))
   enddo
   write(unit, '(A)') '# grid and slice points: factor_length; time (fields, slices, histories): factor_time'
   write(unit, '(A)') '# conservation history: int_<q> = factor <q> * factor_length^3'
   write(unit, '(A)') '# div(B) history: max_divb, seam_max_divb = factor divb; l1_divb = factor bmag * factor_length^2'
   write(unit, '(A)') '# residuals history: residual of <q> = factor dq_<q>'
   close(unit)
   endsubroutine save_units

   pure function scale(self, dim) result(s)
   !< Return the reference of a dimension, `L0^a u0^b sqrt(rho0)^c`: exact when the references are powers of two
   !< (and `rho0` a power of four).
   class(flume_reference_object), intent(in) :: self   !< Reference layer.
   integer(I4P),                  intent(in) :: dim(3) !< Exponents of length, velocity, square root of density.
   real(R8P)                                 :: s      !< Reference.

   s = (self%length**dim(1) * self%velocity**dim(2)) * sqrt(self%density)**dim(3)
   endfunction scale

   pure function time_output(self) result(s)
   !< Return the output factor of the time: `L0 / u0` with dimensional output units, 1 otherwise.
   class(flume_reference_object), intent(in) :: self !< Reference layer.
   real(R8P)                                 :: s    !< Output factor.

   s = 1._R8P
   if (self%dimensional_output) s = self%scale(DIM_TIME)
   endfunction time_output

   function variable_output(self, name) result(s)
   !< Return the output factor of a written variable: the reference of its dimension with dimensional output units, 1
   !< otherwise. `dq_<name>` is the time derivative of `<name>`. An unknown name is fatal: no variable may be written in
   !< code units among dimensional ones.
   class(flume_reference_object), intent(in) :: self   !< Reference layer.
   character(*),                  intent(in) :: name   !< Variable name.
   real(R8P)                                 :: s      !< Output factor.
   character(:), allocatable                 :: base   !< Variable name without the `dq_` prefix.
   integer(I4P)                              :: dim(3) !< Variable dimension.

   s = 1._R8P
   if (.not.self%dimensional_output) return
   base = name
   if (index(name, 'dq_') == 1) base = name(4:)
   select case(base)
   case('r', 'rho')
      dim = DIM_DENSITY
   case('ru', 'rv', 'rw')
      dim = DIM_MOMENTUM
   case('rE', 'p', 'pt')
      dim = DIM_PRESSURE
   case('bx', 'by', 'bz', 'Bx', 'By', 'Bz', 'bmag')
      dim = DIM_FIELD
   case('psi')
      dim = DIM_FIELD
      if (self%divergence_control == 'glm') dim = DIM_PSI_GLM
   case('u', 'v', 'w', 'a')
      dim = DIM_VELOCITY
   case('H', 'T')
      dim = DIM_VELOCITY2
   case('beta')
      dim = DIM_NONE
   case('divb')
      dim = DIM_FIELD_PER_LENGTH
   case default
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(output_units)=dimensional: the dimension of the output '// &
                               'variable "'//name//'" is unknown')
   endselect
   if (len(base) /= len(name)) dim = dim - DIM_TIME
   s = self%scale(dim)
   if (base == 'T') s = s / self%gas_constant
   endfunction variable_output

   ! private methods
   subroutine classify(self, file_parameters, section_name, option_name, kind, dim)
   !< Return the kind and the dimension of an option. Every option FLUME reads is listed (issue #49, N2 inventory);
   !< anything else is `KIND_UNKNOWN`, fatal for the caller.
   class(flume_reference_object), intent(in)  :: self            !< Reference layer.
   type(file_ini),                intent(in)  :: file_parameters !< Simulation parameters ini file handler.
   character(*),                  intent(in)  :: section_name    !< Section name.
   character(*),                  intent(in)  :: option_name     !< Option name.
   integer(I4P),                  intent(out) :: kind            !< Option kind.
   integer(I4P),                  intent(out) :: dim(3)          !< Option dimension (if scaled).

   kind = KIND_UNKNOWN
   dim = DIM_NONE
   if (index(section_name, 'initial_conditions_region_') == 1) then
      select case(option_name)
      case('r')
         kind = KIND_SCALED ; dim = DIM_DENSITY
      case('u', 'v', 'w')
         kind = KIND_SCALED ; dim = DIM_VELOCITY
      case('p')
         kind = KIND_SCALED ; dim = DIM_PRESSURE
      case('bx', 'by', 'bz')
         kind = KIND_SCALED ; dim = DIM_FIELD
      case('emin_x', 'emin_y', 'emin_z', 'emax_x', 'emax_y', 'emax_z')
         kind = KIND_SCALED ; dim = DIM_LENGTH
      endselect
      return
   endif
   if (index(section_name, 'amr_marker_') == 1) then
      select case(option_name)
      case('mode', 'delta_type', 'geo_type', 'solid', 'target_level', 'stl_filename', 'field', 'var')
         kind = KIND_NONE
      case('delta_fine', 'delta_coarse', 'box_xmin', 'box_ymin', 'box_zmin', 'box_xmax', 'box_ymax', 'box_zmax')
         kind = KIND_SCALED ; dim = DIM_LENGTH
      case('tol')
         call self%classify_marker(file_parameters=file_parameters, section_name=section_name, kind=kind, dim=dim)
      endselect
      return
   endif
   if (index(section_name, 'solid_') == 1) then
      select case(option_name)
      case('name', 'definition', 'bc_type', 'circle_axis', 'rectangle_axis')
         kind = KIND_NONE
      case('sphere_center_x', 'sphere_center_y', 'sphere_center_z', 'sphere_radius',             &
           'circle_center_x', 'circle_center_y', 'circle_center_z', 'circle_radius',             &
           'rectangle_center_x', 'rectangle_center_y', 'rectangle_center_z', 'rectangle_edge_1', &
           'rectangle_edge_2')
         kind = KIND_SCALED ; dim = DIM_LENGTH
      endselect
      return
   endif
   if (index(section_name, 'slice_') == 1) then
      select case(option_name)
      case('itype', 'n_save', 'ni', 'nj', 'nk')
         kind = KIND_NONE
      case('emin_x', 'emin_y', 'emin_z', 'emax_x', 'emax_y', 'emax_z')
         kind = KIND_SCALED ; dim = DIM_LENGTH
      endselect
      return
   endif
   select case(section_name)
   case(INI_SECTION_NAME)
      kind = KIND_NONE
   case('numerics', 'runge_kutta', 'weno', 'linear-algebra', 'fdv', 'field', 'amr', 'solids', 'slices', 'diagnostics')
      kind = KIND_NONE ! dimensionless sections (every option a count, a flag, a name or a dimensionless number)
   case('physics')
      select case(option_name)
      case('physical_model', 'gamma', 'reynolds', 'prandtl', 'magnetic_reynolds', 'viscosity_law', 'viscosity_exponent')
         kind = KIND_NONE
      case('lundquist')
         ! S = L v_A / eta is the magnetic Reynolds number in the units of the Alfvenic preset only
         if (self%velocity_preset /= 'alfvenic') &
            call mpih%error_stop(msg=': [physics].(lundquist) needs ['//INI_SECTION_NAME//'].(velocity)=alfvenic (the '// &
                                     'Lundquist number is the magnetic Reynolds number at the Alfven speed); give '// &
                                     'magnetic_reynolds or resistivity otherwise')
         kind = KIND_NONE
      case('cp', 'cv')
         kind = KIND_GAS
      case('viscosity')
         kind = KIND_SCALED ; dim = DIM_VISCOSITY
      case('conductivity')
         kind = KIND_CONDUCTIVITY ; dim = DIM_VISCOSITY
      case('resistivity')
         kind = KIND_SCALED ; dim = DIM_DIFFUSIVITY
      case('reference_temperature')
         kind = KIND_TEMPERATURE ; dim = DIM_VELOCITY2
      endselect
   case('mhd')
      select case(option_name)
      case('divergence_control', 'glm_alpha', 'glm_ch_check', 'divb_error')
         kind = KIND_NONE
      case('glm_ch')
         kind = KIND_SCALED ; dim = DIM_VELOCITY
      case('glm_damping_length')
         kind = KIND_SCALED ; dim = DIM_LENGTH ! or the literal min-cell, left as it is
      case('rho_floor')
         kind = KIND_SCALED ; dim = DIM_DENSITY
      case('p_floor')
         kind = KIND_SCALED ; dim = DIM_PRESSURE
      case('divb_tol')
         kind = KIND_SCALED ; dim = DIM_FIELD_PER_LENGTH
      endselect
   case('IO')
      select case(option_name)
      case('output_basename', 'it_save', 'restart', 'restart_basename', 'restart_save', 'residuals_save',         &
           'divergence_history_save', 'save_memory_status', 'save_residual_fields', 'save_curl_fields',             &
           'save_divergence_fields', 'save_gradient_fields', 'save_laplacian_fields', 'seam_divB_error',            &
           'save_auxiliary_fields')
         kind = KIND_NONE
      case('seam_divB_tol')
         kind = KIND_SCALED ; dim = DIM_FIELD_PER_LENGTH
      endselect
   case('bc_x_min', 'bc_x_max', 'bc_y_min', 'bc_y_max', 'bc_z_min', 'bc_z_max')
      select case(option_name)
      case('type')
         kind = KIND_NONE
      case('r')
         kind = KIND_SCALED ; dim = DIM_DENSITY
      case('u', 'v', 'w', 'wall_u', 'wall_v', 'wall_w')
         kind = KIND_SCALED ; dim = DIM_VELOCITY
      case('wall_temperature')
         kind = KIND_TEMPERATURE ; dim = DIM_VELOCITY2
      case('p')
         kind = KIND_SCALED ; dim = DIM_PRESSURE
      case('bx', 'by', 'bz')
         kind = KIND_SCALED ; dim = DIM_FIELD
      endselect
   case('time')
      select case(option_name)
      case('it_max', 'CFL')
         kind = KIND_NONE
      case('time_max')
         kind = KIND_SCALED ; dim = DIM_TIME
      endselect
   case('grid')
      select case(option_name)
      case('ni', 'nj', 'nk', 'ngc', 'null_x', 'null_y', 'null_z')
         kind = KIND_NONE
      case('emin_x', 'emin_y', 'emin_z', 'emax_x', 'emax_y', 'emax_z')
         kind = KIND_SCALED ; dim = DIM_LENGTH
      endselect
   case('initial_conditions')
      select case(option_name)
      case('amr_iterations', 'type', 's', 'pulse_axis', 'wave', 'wave_angle', 'polarisation', 'axis', 'regions_number', &
           'normal_x', 'normal_y')
         kind = KIND_NONE
      case('x0', 'y0', 'radius', 'r0', 'r1', 'pulse_center', 'pulse_width', 'peak_x0', 'peak_y0', 'peak_radius', &
           'wavelength', 'loop_radius', 'interface_1', 'interface_2', 'period', 'interface_2_width', 'interface')
         kind = KIND_SCALED ; dim = DIM_LENGTH
      case('strength', 'kappa', 'v0')
         kind = KIND_SCALED ; dim = DIM_VELOCITY
      case('rho_in', 'rho_amplitude')
         kind = KIND_SCALED ; dim = DIM_DENSITY
      case('pulse_amplitude', 'peak_amplitude', 'b_par', 'mu', 'loop_amplitude')
         kind = KIND_SCALED ; dim = DIM_FIELD
      case('rho_wavenumber', 'gradient_x', 'gradient_y', 'gradient_z')
         kind = KIND_SCALED ; dim = DIM_INV_LENGTH
      case('wave_amplitude')
         ! cpaw: multiplies B_perp; linear wave: multiplies a Stone et al. 2008 right eigenvector, whose density
         ! component is dimensionless for the fast, slow and entropy waves, the momentum one for the Alfven wave
         kind = KIND_SCALED
         if (self%ic_type == 'mhd-cpaw') then
            dim = DIM_FIELD
         elseif (self%ic_wave == 'alfven') then
            dim = DIM_MOMENTUM
         else
            dim = DIM_DENSITY
         endif
      endselect
   endselect
   endsubroutine classify

   subroutine classify_marker(self, file_parameters, section_name, kind, dim)
   !< Return the dimension of the gradient tolerance of an AMR marker: that of the marked variable per length.
   class(flume_reference_object), intent(in)  :: self             !< Reference layer.
   type(file_ini),                intent(in)  :: file_parameters  !< Simulation parameters ini file handler.
   character(*),                  intent(in)  :: section_name     !< Marker section name.
   integer(I4P),                  intent(out) :: kind             !< Option kind.
   integer(I4P),                  intent(out) :: dim(3)           !< Option dimension.
   integer(I4P)                               :: mode, field, var !< Marker mode, marked field and variable.
   integer(I4P)                               :: error            !< Error status.

   kind = KIND_NONE
   dim = DIM_NONE
   mode = 0_I4P
   call file_parameters%get(section_name=section_name, option_name='mode', val=mode, error=error)
   if (mode /= 2_I4P) return ! the tolerance of a geometric marker is not used
   kind = KIND_SCALED
   field = 0_I4P ; var = 0_I4P
   call file_parameters%get(section_name=section_name, option_name='field', val=field, error=error)
   call file_parameters%get(section_name=section_name, option_name='var', val=var, error=error)
   select case(field)
   case(1_I4P) ! conservative variables
      select case(var)
      case(1_I4P)
         dim = DIM_DENSITY
      case(2_I4P:4_I4P)
         dim = DIM_MOMENTUM
      case(5_I4P)
         dim = DIM_PRESSURE
      case(6_I4P:8_I4P)
         dim = DIM_FIELD
      case(9_I4P)
         dim = DIM_FIELD
         if (self%divergence_control == 'glm') dim = DIM_PSI_GLM
      case default
         kind = KIND_UNKNOWN
      endselect
   case(2_I4P) ! auxiliary variables
      select case(var)
      case(1_I4P)
         dim = DIM_DENSITY
      case(2_I4P:4_I4P, 8_I4P)
         dim = DIM_VELOCITY
      case(5_I4P)
         dim = DIM_PRESSURE
      case(7_I4P)
         dim = DIM_VELOCITY2
      case(9_I4P:11_I4P)
         dim = DIM_FIELD
      case default ! 6, the temperature: its unit needs the dimensional gas constant, not converted
         kind = KIND_UNKNOWN
      endselect
   case default
      kind = KIND_UNKNOWN
   endselect
   if (kind == KIND_SCALED) dim = dim - DIM_LENGTH
   endsubroutine classify_marker

   subroutine convert_gas(self, file_parameters)
   !< Replace `cp`, `cv` by `gamma = cp / cv`: the gas constant becomes 1 (`[physics] gamma` alone means `R = 1`) and
   !< the temperature unit `u0^2 / R`. `gamma` is computed as the physics object computes it, so it is the same double
   !< and the conservative solution is unchanged bit for bit; only the temperature, `p / (rho R)`, changes unit.
   class(flume_reference_object), intent(inout) :: self            !< Reference layer.
   type(file_ini),                intent(inout) :: file_parameters !< Simulation parameters ini file handler.
   real(R8P)                                    :: cp, cv          !< Specific heats.
   integer(I4P)                                 :: error(2)        !< Error status.

   call file_parameters%get(section_name='physics', option_name='cp', val=cp, error=error(1))
   call file_parameters%get(section_name='physics', option_name='cv', val=cv, error=error(2))
   if (any(error > 0)) return ! the physics object reports the missing or inconsistent keys
   if (.not.(cp > cv .and. cv > 0._R8P)) return
   call file_parameters%del(section_name='physics', option_name='cp')
   call file_parameters%del(section_name='physics', option_name='cv')
   call file_parameters%add(section_name='physics', option_name='gamma', val=trim(format_value(cp / cv)))
   self%gas_constant = cp - cv
   self%converted = self%converted + 2_I4P
   if (mpih%myrank == 0_I4P) print '(A)', mpih%myrankstr//'  [physics] cp, cv: '//trim(format_value(cp))//', '// &
                                         trim(format_value(cv))//' -> gamma '//trim(format_value(cp / cv))//', R = 1'
   endsubroutine convert_gas

   subroutine convert_options(self, file_parameters)
   !< Convert every dimensional option of the input to code units; an unclassified section or option is fatal.
   class(flume_reference_object), intent(inout) :: self            !< Reference layer.
   type(file_ini),                intent(inout) :: file_parameters !< Simulation parameters ini file handler.
   character(len=:), allocatable                :: sections(:)     !< Sections names.
   character(len=:), allocatable                :: pair(:)         !< Option name/value pair.
   character(999), allocatable                  :: names(:)        !< Option names of a section.
   character(999)                               :: key, val_str    !< Option name and value, fixed-length buffers.
   character(999), allocatable                  :: opt_values(:)   !< Option values of a section.
   character(:), allocatable                    :: sname           !< Section name.
   integer(I4P)                                 :: kind            !< Option kind.
   integer(I4P)                                 :: dim(3)          !< Option dimension.
   real(R8P)                                    :: val             !< Option value.
   real(R8P)                                    :: converted       !< Option value in code units.
   integer(I4P)                                 :: s, o, n, ios    !< Counters, I/O status.
   logical                                      :: has_gas         !< True if `cp` or `cv` are given.

   if (self%ic_type == 'orszag-tang') &
      call mpih%error_stop(msg=': [initial_conditions].(type)=orszag-tang has its state in code units and takes no '// &
                               'parameter; it cannot be used with ['//INI_SECTION_NAME//']')
   has_gas = .false.
   call file_parameters%get_sections_list(sections)
   do s=1, size(sections, dim=1)
      key = sections(s) ! fixed-length buffer first (nvfortran, see load_from_file)
      sname = trim(key)
      if (len(sname) == 0) cycle ! global section
      ! collect the options first: the conversion rewrites them
      n = 0_I4P
      allocate(names(0), opt_values(0))
      do while (file_parameters%loop(section_name=sname, option_pairs=pair))
         key = pair(1)
         val_str = pair(2)
         names = [character(999) :: names, adjustl(strip_control(key))]
         opt_values = [character(999) :: opt_values, adjustl(strip_control(val_str))]
         n = n + 1_I4P
      enddo
      do o=1, n
         call self%classify(file_parameters=file_parameters, section_name=sname, option_name=trim(names(o)), &
                            kind=kind, dim=dim)
         select case(kind)
         case(KIND_UNKNOWN)
            call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'] is active but ['//sname//'].('//trim(names(o))// &
                                     ') is not classified: its dimension is unknown, refusing to leave it unconverted')
         case(KIND_GAS)
            has_gas = .true.
         case(KIND_SCALED, KIND_TEMPERATURE, KIND_CONDUCTIVITY)
            if (sname == 'mhd' .and. trim(names(o)) == 'glm_damping_length' .and. &
                trim(adjustl(strip_control(opt_values(o)))) == 'min-cell') cycle
            read(opt_values(o), *, iostat=ios) val
            if (ios /= 0) call mpih%error_stop(msg=': ['//sname//'].('//trim(names(o))//') "'//trim(opt_values(o))// &
                                                   '" is not a number')
            ! the temperature unit is u0^2 / R and the conductivity carries 1 / R (issue #65): the reference first
            ! (a power of two in the verifications, exact), then R, so a single rounding
            converted = val / self%scale(dim)
            if (kind == KIND_TEMPERATURE) converted = converted * self%gas_constant
            if (kind == KIND_CONDUCTIVITY) converted = converted / self%gas_constant
            call file_parameters%add(section_name=sname, option_name=trim(names(o)), &
                                     val=trim(format_value(converted)))
            self%converted = self%converted + 1_I4P
            if (mpih%myrank == 0_I4P) print '(A)', mpih%myrankstr//'  ['//sname//'] '//trim(names(o))//': '// &
                                                  trim(format_value(val))//' -> '//trim(format_value(converted))
         endselect
      enddo
      deallocate(names, opt_values)
   enddo
   if (has_gas) call self%convert_gas(file_parameters=file_parameters)
   endsubroutine convert_options

   subroutine load_context(self, file_parameters)
   !< Load, raw, the options on which the dimension of other options depends.
   class(flume_reference_object), intent(inout) :: self            !< Reference layer.
   type(file_ini),                intent(in)    :: file_parameters !< Simulation parameters ini file handler.

   self%physical_model     = raw(section_name='physics',            option_name='physical_model', default='')
   self%divergence_control = raw(section_name='mhd',                option_name='divergence_control', default='none')
   self%ic_type            = raw(section_name='initial_conditions', option_name='type', default='')
   self%ic_wave            = raw(section_name='initial_conditions', option_name='wave', default='')
   ! the gas constant of the input, needed before the options are converted (temperatures and the conductivity carry
   ! it, issue #65): cp - cv, or 1 with gamma alone (the physics object reports missing or inconsistent keys)
   self%gas_constant = 1._R8P
   block
      real(R8P)    :: cp, cv
      integer(I4P) :: error(2)
      call file_parameters%get(section_name='physics', option_name='cp', val=cp, error=error(1))
      call file_parameters%get(section_name='physics', option_name='cv', val=cv, error=error(2))
      if (all(error <= 0)) then
         if (cp > cv .and. cv > 0._R8P) self%gas_constant = cp - cv
      endif
   endblock
   contains
      function raw(section_name, option_name, default) result(val)
      !< Return an option value without control characters, or a default if the option is absent.
      character(*), intent(in)  :: section_name !< Section name.
      character(*), intent(in)  :: option_name  !< Option name.
      character(*), intent(in)  :: default      !< Default value.
      character(:), allocatable :: val          !< Option value.
      character(999)            :: buff         !< Option value buffer.
      integer(I4P)              :: error        !< Error status.

      call file_parameters%get(section_name=section_name, option_name=option_name, val=buff, error=error)
      if (error > 0) then
         val = default
      else
         val = trim(adjustl(strip_control(buff)))
      endif
      endfunction raw
   endsubroutine load_context

   subroutine load_from_file(self, file_parameters)
   !< Load `[reference]`: `density`, `length` (default 1, > 0) and `velocity` (default 1): a value > 0, or `acoustic`
   !< (`sqrt(gamma pressure / density)`, needs `pressure`) or `alfvenic` (`field / sqrt(density)`, needs `field`,
   !< MHD only); `output_units`, `code` (default) or `dimensional`. Every other option is fatal.
   class(flume_reference_object), intent(inout) :: self            !< Reference layer.
   type(file_ini),                intent(inout) :: file_parameters !< Simulation parameters ini file handler.
   character(len=:), allocatable                :: pair(:)         !< Option name/value pair.
   character(:), allocatable                    :: velocity        !< Velocity option value.
   character(:), allocatable                    :: output_units    !< Output units option value.
   character(999)                               :: key, val_str    !< Option name and value, control characters blanked.
   real(R8P)                                    :: gamma           !< Specific heats ratio (acoustic preset).
   real(R8P)                                    :: cp, cv          !< Specific heats (acoustic preset).
   integer(I4P)                                 :: error(2), ios   !< Error and I/O status.

   velocity = '1'
   output_units = 'code'
   do while (file_parameters%loop(section_name=INI_SECTION_NAME, option_pairs=pair))
      ! copy the pair into fixed-length buffers before any other use: nvfortran 26.1 (-acc -fast) passes an element of
      ! the deferred-length array to a procedure, or matches it in a select case, as garbage (issue #49, N2a)
      key = pair(1)
      val_str = pair(2)
      key = adjustl(strip_control(key))
      val_str = adjustl(strip_control(val_str))
      select case(trim(key))
      case('density')
         self%density = positive(key, val_str)
      case('length')
         self%length = positive(key, val_str)
      case('velocity')
         velocity = trim(val_str)
      case('pressure')
         self%ref_pressure = positive(key, val_str)
      case('field')
         self%ref_field = positive(key, val_str)
      case('output_units')
         output_units = trim(val_str)
      case default
         call mpih%error_stop(msg=': unknown ['//INI_SECTION_NAME//'].('//trim(key)//'); accepted: density, '// &
                                  'length, velocity, pressure, field, output_units')
      endselect
   enddo
   select case(output_units)
   case('code')
      self%dimensional_output = .false.
   case('dimensional')
      self%dimensional_output = .true.
   case default
      call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(output_units) "'//output_units//'": expected code or '// &
                               'dimensional')
   endselect
   select case(velocity)
   case('acoustic')
      self%velocity_preset = 'acoustic'
      if (self%ref_pressure <= 0._R8P) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(velocity)=acoustic needs ['//INI_SECTION_NAME// &
                                  '].(pressure) > 0')
      call file_parameters%get(section_name='physics', option_name='gamma', val=gamma, error=error(1))
      if (error(1) > 0) then
         call file_parameters%get(section_name='physics', option_name='cp', val=cp, error=error(1))
         call file_parameters%get(section_name='physics', option_name='cv', val=cv, error=error(2))
         if (any(error > 0) .or. .not.(cp > cv .and. cv > 0._R8P)) &
            call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(velocity)=acoustic needs [physics].(gamma) or '// &
                                     'cp > cv > 0')
         gamma = cp / cv
      endif
      self%velocity = sqrt(gamma * self%ref_pressure / self%density)
   case('alfvenic')
      self%velocity_preset = 'alfvenic'
      if (self%physical_model /= 'mhd-ideal') &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(velocity)=alfvenic needs [physics].(physical_model)='// &
                                  'mhd-ideal')
      if (self%ref_field <= 0._R8P) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(velocity)=alfvenic needs ['//INI_SECTION_NAME// &
                                  '].(field) > 0')
      self%velocity = self%ref_field / sqrt(self%density)
   case default
      self%velocity_preset = 'value'
      read(velocity, *, iostat=ios) self%velocity
      if (ios /= 0 .or. .not.(self%velocity > 0._R8P)) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].(velocity) "'//velocity//'": expected a value > 0, '// &
                                  'acoustic or alfvenic')
   endselect
   contains
      function positive(opt_name, opt_value) result(val)
      !< Return the value of an option, which must be a number > 0.
      character(*), intent(in) :: opt_name  !< Option name.
      character(*), intent(in) :: opt_value !< Option value.
      real(R8P)                :: val       !< Value.
      integer(I4P)             :: ios_      !< I/O status.

      read(opt_value, *, iostat=ios_) val
      if (ios_ /= 0 .or. .not.(val > 0._R8P)) &
         call mpih%error_stop(msg=': ['//INI_SECTION_NAME//'].('//trim(opt_name)//') "'//trim(opt_value)// &
                                  '": expected a value > 0')
      endfunction positive
   endsubroutine load_from_file

   ! private procedures
   pure function format_value(val) result(str_val)
   !< Return a value with 17 significant digits: read back, it is the same double.
   real(R8P), intent(in) :: val     !< Value.
   character(len=32)     :: str_val !< Formatted value.

   write(str_val, '(ES25.16E3)') val
   str_val = adjustl(str_val)
   endfunction format_value
endmodule adam_flume_reference_object
