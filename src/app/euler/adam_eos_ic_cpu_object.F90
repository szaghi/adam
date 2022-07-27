!< ADAM, **EOS** (Equation of State) of ideal compressible fluid, class definition, CPU backend.
module adam_eos_ic_cpu_object
!< ADAM, **EOS** (Equation of State) of ideal compressible fluid, class definition, CPU backend.

use finer, only : file_ini
use penf, only : I4P, R8P, str

implicit none
private
public :: eos_ic_cpu_object

character(len=3), parameter :: INI_SECTION_NAME='eos' !< INI (config) file section name containing equation of state conditions.

type :: eos_ic_cpu_object
   !< Equation of state (EOS) of ideal compressible object class.
   real(R8P) :: cp    = 0._R8P !< Specific heat at constant pressure `cp`.
   real(R8P) :: cv    = 0._R8P !< Specific heat at constant volume `cv`.
   real(R8P) :: g     = 0._R8P !< Specific heats ratio `gamma = cp / cv`.
   real(R8P) :: R     = 0._R8P !< Fluid constant `R = cp - cv`.
   real(R8P) :: gm1   = 0._R8P !< `gamma - 1`.
   real(R8P) :: gp1   = 0._R8P !< `gamma + 1`.
   real(R8P) :: delta = 0._R8P !< `(gamma - 1) / 2`.
   real(R8P) :: eta   = 0._R8P !< `2 * gamma / (gamma - 1)`.
   contains
      ! public methods
      procedure, pass(self) :: compute_derivate !< Compute derivate quantities (from `cp` and `cv`).
      procedure, pass(self) :: density          !< Return density.
      procedure, pass(self) :: description      !< Return pretty-printed object description.
      procedure, pass(self) :: destroy          !< Destroy eos.
      procedure, pass(self) :: initialize       !< initialize eos.
      procedure, pass(self) :: internal_energy  !< Return specific internal energy.
      procedure, pass(self) :: load_from_file   !< Load from finer file.
      procedure, pass(self) :: pressure         !< Return pressure.
      procedure, pass(self) :: save_into_file   !< Save into finer file.
      procedure, pass(self) :: speed_of_sound   !< Return speed of sound.
      procedure, pass(self) :: temperature      !< Return temperature.
      procedure, pass(self) :: total_entalpy    !< Return total specific entalpy.
      ! operators
      generic :: assignment(=) => eos_assign_eos !< Overload `=`.
      ! private methods
      procedure, pass(lhs), private :: eos_assign_eos !< Operator `=`.
endtype eos_ic_cpu_object

interface eos_ic_cpu_object
   !< Overload [[eos_ic_cpu_object]] name with its constructor.
   module procedure eos_ic_cpu_instance
endinterface
contains
   ! public methods
   elemental subroutine compute_derivate(self)
   !< Compute derivate quantities (from `cp` and `cv`).
   class(eos_ic_cpu_object), intent(inout) :: self !< Equation of state.

   self%g     = self%cp / self%cv
   self%R     = self%cp - self%cv
   self%gm1   = self%g - 1._R8P
   self%gp1   = self%g + 1._R8P
   self%delta = (self%g - 1._R8P) * 0.5_R8P
   self%eta   = 2._R8P * self%g / (self%g - 1._R8P)
   endsubroutine compute_derivate

   elemental function density(self, pressure, speed_of_sound) result(density_)
   !< Return density.
   class(eos_ic_cpu_object), intent(in) :: self           !< Equation of state.
   real(R8P),                intent(in) :: pressure       !< Pressure value.
   real(R8P),                intent(in) :: speed_of_sound !< Speed of sound value.
   real(R8P)                            :: density_       !< Density value.

   density_ = self%g * pressure / (speed_of_sound * speed_of_sound)
   endfunction density

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted object description.
   class(eos_ic_cpu_object), intent(in)           :: self             !< Equation of state.
   character(*),             intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                  :: desc             !< Description.
   character(len=:), allocatable                  :: prefix_          !< Prefixing string, local variable.
   character(len=1), parameter                    :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//prefix_//'cp  = '//trim(str(n=self%cp))//NL
   desc = desc//prefix_//'cv  = '//trim(str(n=self%cv))
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy eos.
   class(eos_ic_cpu_object), intent(inout) :: self  !< Equation of state.
   type(eos_ic_cpu_object)                 :: fresh !< Fresh eos.

   self = fresh
   endsubroutine destroy

   elemental subroutine initialize(self, eos, cp, cv, gam, R)
   !< Initialize equation of state.
   class(eos_ic_cpu_object), intent(inout)        :: self !< Equation of state.
   type(eos_ic_cpu_object),  intent(in), optional :: eos  !< Equation of state.
   real(R8P),                intent(in), optional :: cp   !< Specific heat at constant pressure `cp` value.
   real(R8P),                intent(in), optional :: cv   !< Specific heat at constant volume `cv` value.
   real(R8P),                intent(in), optional :: gam  !< Specific heats ratio `gamma=cp/cv` value.
   real(R8P),                intent(in), optional :: R    !< Fluid constant `R=cp-cv` value.

   if (present(eos)) then
      self = eos
   else
      self = eos_ic_cpu_object(cp=cp, cv=cv, gam=gam, R=R)
   endif
   endsubroutine initialize

   elemental function internal_energy(self, density, pressure, temperature) result(energy_)
   !< Return specific internal energy.
   class(eos_ic_cpu_object), intent(in)           :: self        !< Equation of state.
   real(R8P),                intent(in), optional :: density     !< Density value.
   real(R8P),                intent(in), optional :: pressure    !< Pressure value.
   real(R8P),                intent(in), optional :: temperature !< Temperature value.
   real(R8P)                                      :: energy_     !< Energy value.

   energy_ = 0._R8P
   if (present(density).and.present(pressure)) then
      energy_ = pressure / ((self%g - 1._R8P) * density)
   elseif (present(temperature)) then
      energy_ = self%cv * temperature
   endif
   endfunction internal_energy

   subroutine load_from_file(self, fini, go_on_fail)
   !< Load from file.
   class(eos_ic_cpu_object), intent(inout)        :: self        !< Equation of state.
   type(file_ini),          intent(in)           :: fini        !< Simulation parameters ini file handler.
   logical,                 intent(in), optional :: go_on_fail  !< Go on if load fails.
   logical                                       :: go_on_fail_ !< Go on if load fails.
   integer(I4P)                                  :: error       !< Error status.

   go_on_fail_ = .false. ; if (present(go_on_fail)) go_on_fail_ = go_on_fail

   call self%destroy
   call fini%get(section_name=INI_SECTION_NAME, option_name='cp', val=self%cp, error=error)
   if (.not.go_on_fail_.and.error>0) error stop 'error: failed to load ['//INI_SECTION_NAME//'].(cp)'
   call fini%get(section_name=INI_SECTION_NAME, option_name='cv', val=self%cv, error=error)
   if (.not.go_on_fail_.and.error>0) error stop 'error: failed to load ['//INI_SECTION_NAME//'].(cv)'
   call self%compute_derivate
   endsubroutine load_from_file

   elemental function pressure(self, density, energy, temperature) result(pressure_)
   !< Return pressure.
   class(eos_ic_cpu_object), intent(in)           :: self        !< Equation of state.
   real(R8P),                intent(in), optional :: density     !< Density value.
   real(R8P),                intent(in), optional :: energy      !< Specific internal energy value.
   real(R8P),                intent(in), optional :: temperature !< Temperature value.
   real(R8P)                                      :: pressure_   !< Pressure value.

   pressure_ = 0._R8P
   if (present(density).and.present(energy)) then
      pressure_ = density * (self%g - 1._R8P) * energy
   elseif (present(density).and.present(temperature)) then
      pressure_ = density * self%R * temperature
   endif
   endfunction pressure

   subroutine save_into_file(self, fini)
   !< Save into file.
   class(eos_ic_cpu_object), intent(inout) :: self !< Free conditions.
   type(file_ini),           intent(inout) :: fini !< Simulation parameters ini file handler.

   call fini%add(section_name=INI_SECTION_NAME, option_name='cp', val=self%cp)
   call fini%add(section_name=INI_SECTION_NAME, option_name='cv', val=self%cv)
   endsubroutine save_into_file

   elemental function speed_of_sound(self, density, pressure) result(speed_of_sound_)
   !< Return speed of sound.
   class(eos_ic_cpu_object), intent(in) :: self            !< Equation of state.
   real(R8P),                intent(in) :: density         !< Density value.
   real(R8P),                intent(in) :: pressure        !< Pressure value.
   real(R8P)                            :: speed_of_sound_ !< Speed of sound value.

   speed_of_sound_ = sqrt(self%g * pressure / density)
   endfunction speed_of_sound

   elemental function temperature(self, density, energy, pressure) result(temperature_)
   !< Return temperature.
   class(eos_ic_cpu_object), intent(in)           :: self         !< Equation of state.
   real(R8P),                intent(in), optional :: density      !< Density value.
   real(R8P),                intent(in), optional :: energy       !< Specific internal energy value.
   real(R8P),                intent(in), optional :: pressure     !< Pressure value.
   real(R8P)                                      :: temperature_ !< Temperature value.

   temperature_ = 0._R8P
   if (present(density).and.present(pressure)) then
      temperature_ = pressure / (self%R * density)
   elseif (present(energy)) then
      temperature_ = energy / self%cv
   endif
   endfunction temperature

   elemental function total_entalpy(self, density, pressure, velocity_sq_norm) result(entalpy_)
   !< Return total specific entalpy.
   class(eos_ic_cpu_object), intent(in) :: self             !< Equation of state.
   real(R8P),                intent(in) :: density          !< Density value.
   real(R8P),                intent(in) :: pressure         !< Pressure value.
   real(R8P),                intent(in) :: velocity_sq_norm !< Velocity vector square norm `||velocity||^2`.
   real(R8P)                            :: entalpy_         !< Total specific entalpy (per unit of mass).

   entalpy_ = self%g * pressure/(self%gm1 * density) + 0.5_R8P * velocity_sq_norm
   endfunction total_entalpy

   ! operators
   pure subroutine eos_assign_eos(lhs, rhs)
   !< Operator `=`.
   class(eos_ic_cpu_object), intent(inout) :: lhs !< Left hand side.
   type(eos_ic_cpu_object),  intent(in)    :: rhs !< Right hand side.

   lhs%cp    = rhs%cp
   lhs%cv    = rhs%cv
   lhs%g     = rhs%g
   lhs%R     = rhs%R
   lhs%delta = rhs%delta
   lhs%eta   = rhs%eta
   lhs%gm1   = rhs%gm1
   lhs%gp1   = rhs%gp1
   endsubroutine eos_assign_eos

   ! private non TBP
   elemental function eos_ic_cpu_instance(cp, cv, gam, R) result(instance)
   !< Return and instance of [[eos_ic_cpu_object]].
   !<
   !< @note This procedure is used for overloading [[eos_ic_cpu_object]] name.
   real(R8P), intent(in), optional :: cp       !< Specific heat at constant pressure `cp` value.
   real(R8P), intent(in), optional :: cv       !< Specific heat at constant volume `cv` value.
   real(R8P), intent(in), optional :: gam      !< Specific heats ratio `gamma=cp/cv` value.
   real(R8P), intent(in), optional :: R        !< Fluid constant `R=cp-cv` value.
   type(eos_ic_cpu_object)         :: instance !< Instance of [[eos_ic_cpu_object]].

   if (present(cp).and.present(cv)) then
      instance%cp = cp
      instance%cv = cv
   elseif (present(gam).and.present(R)) then
      instance%cv = R/(gam - 1._R8P)
      instance%cp = gam * instance%cv
   elseif (present(gam).and.present(cp)) then
      instance%cp = cp
      instance%cv = cp / gam
   elseif (present(gam).and.present(cv)) then
      instance%cp = gam * cv
      instance%cv = cv
   elseif (present(R).and.present(cp)) then
      instance%cp = cp
      instance%cv = cp - R
   elseif (present(R).and.present(cv)) then
      instance%cp = cv + R
      instance%cv = cv
   endif
   call instance%compute_derivate
   endfunction eos_ic_cpu_instance
endmodule adam_eos_ic_cpu_object
